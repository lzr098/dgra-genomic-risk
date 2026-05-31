#!/usr/bin/env python3
"""L2 unit tests for gpa_workflow.py, gpa_workflow_pm.py, gpa_workflow_runner.py"""
import pytest
from unittest.mock import MagicMock, patch
import asyncio

import gpa_workflow as wf
from gpa_workflow import WorkflowStep, FailureAction, STANDARD_WORKFLOW
import gpa_workflow_pm as pm
import gpa_workflow_runner as wr


# =============================================================================
# gpa_workflow
# =============================================================================

class TestWorkflowStep:
    def test_can_be_skipped_required_always_false(self):
        step = WorkflowStep(name="tier_classification", module="gpa_tier_classifier", required=True)
        assert step.can_be_skipped({}) == (False, "")

    def test_required_steps_never_skipped(self):
        # v0.11.0: gtex_expression, clinvar_ncbi, spliceai_prediction, gnomad_variant_frequency
        # are all required=True, so can_be_skipped always returns (False, "")
        for name in ["gtex_expression", "clinvar_ncbi", "spliceai_prediction", "gnomad_variant_frequency"]:
            step = next(s for s in STANDARD_WORKFLOW if s.name == name)
            assert step.required is True
            assert step.can_be_skipped({}) == (False, "")

    def test_skip_myvariant_all_have_af(self):
        step = next(s for s in STANDARD_WORKFLOW if s.name == "myvariant_enrichment")
        variants = [{"gnomAD_AF": 0.01}, {"gnomAD_AF": 0.02}]
        assert step.can_be_skipped({"variants": variants}) == (True, "All variants already have gnomAD AF data")
        variants = [{"gnomAD_AF": None}]
        assert step.can_be_skipped({"variants": variants}) == (False, "")

    def test_skip_vep_not_raw_vcf(self):
        step = next(s for s in STANDARD_WORKFLOW if s.name == "vep_annotation")
        assert step.can_be_skipped({"input_type": "ANNOTATED_VCF"}) == (True, "Input type is ANNOTATED_VCF, not RAW_VCF — VEP annotation not needed")
        assert step.can_be_skipped({"input_type": "RAW_VCF"}) == (False, "")

    def test_skip_transcript_not_raw_vcf(self):
        step = next(s for s in STANDARD_WORKFLOW if s.name == "transcript_selection")
        assert step.can_be_skipped({"input_type": "ANNOTATED_VCF"}) == (True, "Input type is ANNOTATED_VCF — transcript selection already done")

    def test_skip_phenotype_no_user_input(self):
        step = next(s for s in STANDARD_WORKFLOW if s.name == "phenotype_matching")
        assert step.can_be_skipped({"user_phenotypes": ""}) == (True, "No user phenotypes provided")
        assert step.can_be_skipped({"user_phenotypes": "seizures"}) == (False, "")

    def test_skip_vep_reannotation_no_discrepancy(self):
        step = next(s for s in STANDARD_WORKFLOW if s.name == "vep_reannotation")
        assert step.can_be_skipped({"discrepancy_count": 0}) == (True, "No TRANSCRIPT_DISCREPANCY variants to reannotate")
        assert step.can_be_skipped({"discrepancy_count": 5, "offline_mode": True}) == (True, "Offline mode — VEP reannotation requires API access")
        assert step.can_be_skipped({"discrepancy_count": 5, "offline_mode": False}) == (False, "")


class TestWorkflowQueries:
    def test_get_step_by_name(self):
        step = wf.get_step_by_name("tier_classification")
        assert step is not None
        assert step.module == "gpa_tier_classifier"

    def test_get_step_by_name_missing(self):
        assert wf.get_step_by_name("nonexistent") is None

    def test_get_required_steps(self):
        required = wf.get_required_steps()
        assert all(s.required for s in required)
        assert len(required) > 0

    def test_get_optional_steps(self):
        optional = wf.get_optional_steps()
        assert all(not s.required for s in optional)
        assert len(optional) > 0

    def test_get_api_steps(self):
        api_steps = wf.get_api_steps()
        assert all(s.api_name is not None for s in api_steps)


class TestValidateWorkflow:
    def test_no_errors(self):
        errors = wf.validate_workflow()
        assert errors == []

    def test_duplicate_name_detected(self):
        with patch.object(wf, "STANDARD_WORKFLOW", STANDARD_WORKFLOW + [STANDARD_WORKFLOW[0]]):
            errors = wf.validate_workflow()
            assert any("Duplicate" in e for e in errors)


# =============================================================================
# gpa_workflow_pm
# =============================================================================

class TestPhaseInfo:
    def test_progress_pct(self):
        phase = pm.PhaseInfo(name="test", start_time=0, total_items=100, processed_items=25)
        assert phase.progress_pct == 25.0

    def test_progress_pct_zero_total(self):
        phase = pm.PhaseInfo(name="test", start_time=0, total_items=0)
        assert phase.progress_pct == 0.0

    def test_eta_sec(self):
        import time
        start = time.time() - 10  # 10 seconds ago
        phase = pm.PhaseInfo(name="test", start_time=start, total_items=100, processed_items=50)
        eta = phase.eta_sec()
        assert eta is not None
        assert eta > 0  # Should be ~10s remaining

    def test_eta_sec_no_progress(self):
        phase = pm.PhaseInfo(name="test", start_time=0, total_items=100, processed_items=0)
        assert phase.eta_sec() is None


class TestWorkflowPM:
    def test_init_defaults(self):
        monitor = pm.WorkflowPM()
        assert monitor.report_interval_sec == 180
        assert monitor.milestone_intervals == [0.25, 0.50, 0.75, 1.0]

    def test_start_end_analysis(self):
        monitor = pm.WorkflowPM(report_interval_sec=3600)  # Long interval to avoid timer firing
        monitor.start_analysis()
        assert monitor._overall_start_time is not None
        monitor.end_analysis()
        assert monitor._timer is None

    def test_start_phase(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("vep_annotation", total_items=100, detail="Batch 1")
        assert "vep_annotation" in monitor.phases
        assert monitor.current_phase == "vep_annotation"

    def test_update_progress(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test", total_items=100)
        monitor.update_progress(processed=50, detail="variant 50")
        assert monitor.phases["test"].processed_items == 50

    def test_increment_progress(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test", total_items=100)
        monitor.increment_progress(amount=5)
        assert monitor.phases["test"].processed_items == 5

    def test_complete_phase(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test", total_items=100)
        monitor.complete_phase()
        assert monitor.phases["test"].status == "completed"
        assert monitor.current_phase is None

    def test_fail_phase(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test", total_items=100)
        monitor.fail_phase("something broke")
        assert monitor.phases["test"].status == "failed"

    def test_add_sub_task(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test")
        monitor.add_sub_task("sub1")
        assert len(monitor.phases["test"].sub_tasks) == 1

    def test_complete_sub_task(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test")
        monitor.add_sub_task("sub1")
        monitor.complete_sub_task("sub1")
        assert monitor.phases["test"].sub_tasks[0]["status"] == "completed"

    def test_get_current_status_idle(self):
        monitor = pm.WorkflowPM()
        status = monitor.get_current_status()
        assert status["status"] == "idle"

    def test_get_current_status_running(self):
        monitor = pm.WorkflowPM()
        monitor.start_phase("test", total_items=100, detail="running")
        monitor.update_progress(processed=50)
        status = monitor.get_current_status()
        assert status["current_phase"] == "test"
        assert status["progress_pct"] == 50.0

    def test_get_summary_report(self):
        monitor = pm.WorkflowPM()
        monitor.start_analysis()
        monitor.start_phase("test", total_items=100)
        monitor.complete_phase()
        report = monitor.get_summary_report()
        assert "Analysis Progress Summary" in report
        assert "test" in report
        assert "✓" in report  # completed icon

    def test_format_duration(self):
        monitor = pm.WorkflowPM()
        assert monitor._format_duration(30) == "30s"
        assert monitor._format_duration(120) == "2.0min"
        assert "h" in monitor._format_duration(4000)

    def test_format_eta(self):
        monitor = pm.WorkflowPM()
        assert monitor._format_eta(None) == "calculating..."
        assert monitor._format_eta(30) == "30s"
        assert "min" in monitor._format_eta(120)
        assert "h" in monitor._format_eta(4000)


class TestGPAPM:
    def test_start_step(self):
        gpm = pm.GPAPM()
        gpm.start_step("vep_annotation", total_items=100)
        assert "VEP Annotation" in gpm.phases

    def test_complete_step(self):
        gpm = pm.GPAPM()
        gpm.start_step("vep_annotation", total_items=100)
        gpm.complete_step("vep_annotation")
        assert gpm.phases["VEP Annotation"].status == "completed"

    def test_complete_step_not_found(self):
        gpm = pm.GPAPM()
        # Should not raise
        gpm.complete_step("nonexistent_step")


# =============================================================================
# gpa_workflow_runner
# =============================================================================

class TestWorkflowRunnerInit:
    def test_run_mode(self):
        runner = wr.WorkflowRunner(mode="run")
        assert runner.mode == "run"

    def test_optimize_mode(self):
        runner = wr.WorkflowRunner(mode="optimize")
        assert runner.mode == "optimize"

    def test_invalid_mode_raises(self):
        with pytest.raises(ValueError):
            wr.WorkflowRunner(mode="invalid")


class TestWorkflowRunnerExecution:
    @pytest.mark.asyncio
    async def test_execute_empty_workflow(self):
        # Use a non-empty sentinel to bypass "workflow or list(STANDARD_WORKFLOW)"
        runner = wr.WorkflowRunner(mode="run", workflow=[])
        runner.workflow = []  # Force empty after init
        report = await runner.execute({})
        assert report.steps_total == 0
        assert report.final_status == "SUCCESS"

    @pytest.mark.asyncio
    async def test_execute_required_step_success(self):
        step = WorkflowStep(name="test_step", module="test", required=True, timeout_sec=1)
        runner = wr.WorkflowRunner(mode="run", workflow=[step])
        report = await runner.execute({})
        assert report.steps_success == 1
        assert report.final_status == "SUCCESS"

    @pytest.mark.asyncio
    async def test_execute_optional_step_skipped(self):
        step = WorkflowStep(
            name="gtex_expression", module="dgra_api", required=False,
            skip_condition="No GTEx tissue", timeout_sec=1,
        )
        runner = wr.WorkflowRunner(mode="run", workflow=[step])
        report = await runner.execute({"tissue_profile": {}})
        assert report.steps_skipped == 1
        assert report.final_status == "SUCCESS"

    @pytest.mark.asyncio
    async def test_execute_timeout(self):
        async def slow_step(ctx):
            await asyncio.sleep(10)
            return {}

        step = WorkflowStep(name="slow", module="test", required=True, timeout_sec=0.01)
        runner = wr.WorkflowRunner(mode="run", workflow=[step])
        with patch.object(runner, "_execute_step", side_effect=slow_step):
            report = await runner.execute({})
        assert report.steps_timeout == 1
        # timeout with ABORT default -> FAILED final
        assert report.final_status == "FAILED"

    @pytest.mark.asyncio
    async def test_execute_step_failure_warn(self):
        step = WorkflowStep(name="test", module="test", required=False, timeout_sec=1, on_failure=FailureAction.WARN)
        runner = wr.WorkflowRunner(mode="run", workflow=[step])
        with patch.object(runner, "_execute_step", side_effect=RuntimeError("boom")):
            report = await runner.execute({})
        assert report.steps_failed == 1
        # Any failure (even non-required) makes final_status PARTIAL
        assert report.final_status == "PARTIAL"

    @pytest.mark.asyncio
    async def test_execute_step_failure_abort(self):
        step = WorkflowStep(name="test", module="test", required=True, timeout_sec=1, on_failure=FailureAction.ABORT)
        runner = wr.WorkflowRunner(mode="run", workflow=[step])
        with patch.object(runner, "_execute_step", side_effect=RuntimeError("boom")):
            report = await runner.execute({})
        assert report.steps_failed == 1
        assert report.final_status == "FAILED"

    @pytest.mark.asyncio
    async def test_execute_notifies_skip(self):
        notifications = []
        step = WorkflowStep(name="gtex_expression", module="dgra_api", required=False,
                            skip_condition="No GTEx tissue", timeout_sec=1)
        runner = wr.WorkflowRunner(mode="run", workflow=[step],
                                   notification_callback=lambda msg: notifications.append(msg))
        await runner.execute({"tissue_profile": {}})
        assert len(notifications) == 1
        assert "SKIPPED" in notifications[0]


class TestExecutionReport:
    def test_to_markdown(self):
        report = wr.ExecutionReport(
            workflow_version="0.11.0",
            mode="run",
            start_time="2024-01-01T00:00:00",
            end_time="2024-01-01T00:01:00",
            total_duration_sec=60.0,
            steps_total=2,
            steps_success=1,
            steps_skipped=1,
            steps_failed=0,
            steps_timeout=0,
            step_results=[
                wr.StepResult("step1", "SUCCESS", 0, 30, 30.0),
                wr.StepResult("step2", "SKIPPED", 30, 60, 30.0, skip_reason="no data"),
            ],
            final_status="SUCCESS",
        )
        md = report.to_markdown()
        assert "GPA Execution Report" in md
        assert "step1" in md
        assert "SKIPPED" in md
        assert "no data" in md


class TestOptimizeMode:
    @pytest.mark.asyncio
    async def test_propose_change(self):
        runner = wr.WorkflowRunner(mode="optimize")
        result = await runner.propose_workflow_change("Add new step")
        assert result is False  # pending confirmation
        assert len(runner.get_pending_changes()) == 1

    @pytest.mark.asyncio
    async def test_confirm_last_change(self):
        runner = wr.WorkflowRunner(mode="optimize")
        new_step = WorkflowStep(name="new", module="test")
        await runner.propose_workflow_change("Add step", new_steps=[new_step])
        assert runner.confirm_last_change() is True
        assert any(s.name == "new" for s in runner.workflow)

    @pytest.mark.asyncio
    async def test_reject_last_change(self):
        runner = wr.WorkflowRunner(mode="optimize")
        await runner.propose_workflow_change("Add step")
        assert runner.reject_last_change() is True
        assert len(runner.get_pending_changes()) == 0

    @pytest.mark.asyncio
    async def test_reset_workflow(self):
        runner = wr.WorkflowRunner(mode="optimize")
        new_step = WorkflowStep(name="new", module="test")
        await runner.propose_workflow_change("Add step", new_steps=[new_step])
        runner.confirm_last_change()
        runner.reset_workflow()
        assert not any(s.name == "new" for s in runner.workflow)

    def test_propose_in_run_mode_raises(self):
        runner = wr.WorkflowRunner(mode="run")
        with pytest.raises(RuntimeError):
            asyncio.run(runner.propose_workflow_change("x"))

    def test_execute_in_optimize_mode_raises(self):
        runner = wr.WorkflowRunner(mode="optimize")
        with pytest.raises(RuntimeError):
            asyncio.run(runner.execute({}))


class TestConvenienceFunctions:
    @pytest.mark.asyncio
    async def test_run_workflow(self):
        report = await wr.run_workflow({})
        assert isinstance(report, wr.ExecutionReport)

    def test_print_workflow_summary(self, capsys):
        wr.print_workflow_summary()
        captured = capsys.readouterr()
        assert "GPA Standard Workflow" in captured.out
        assert "preflight_check" in captured.out
