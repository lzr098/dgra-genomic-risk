#!/usr/bin/env python3
"""
GPA Workflow Runner (v0.11.0)
=============================
Execution engine for the GPA analysis pipeline.

Enforces workflow fidelity:
  - "run" mode: strict execution, no unauthorized skips, mandatory user notification
  - "optimize" mode: workflow modifications require explicit user confirmation

Usage (run mode):
    from gpa_workflow_runner import WorkflowRunner
    runner = WorkflowRunner(mode="run")
    result = await runner.execute(context)

Usage (optimize mode):
    runner = WorkflowRunner(mode="optimize")
    await runner.propose_workflow_change(new_steps)  # prints diff, waits for confirm
"""

import asyncio
import json
import sys
import time
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, List, Optional, Callable

from gpa_workflow import (
    STANDARD_WORKFLOW,
    WorkflowStep,
    FailureAction,
    WORKFLOW_VERSION,
)


# =============================================================================
# Execution Record
# =============================================================================

@dataclass
class StepResult:
    """Result of executing a single workflow step."""
    step_name: str
    status: str  # "SUCCESS" | "SKIPPED" | "FAILED" | "TIMEOUT"
    start_time: float
    end_time: float
    duration_sec: float
    skip_reason: str = ""
    error_message: str = ""
    error_traceback: str = ""
    produced_data: Dict[str, Any] = field(default_factory=dict)
    notification_sent: bool = False  # Whether user was notified of skip


@dataclass
class ExecutionReport:
    """Complete report of a workflow execution."""
    workflow_version: str
    mode: str
    start_time: str
    end_time: str
    total_duration_sec: float
    steps_total: int
    steps_success: int
    steps_skipped: int
    steps_failed: int
    steps_timeout: int
    step_results: List[StepResult]
    final_status: str  # "SUCCESS" | "PARTIAL" | "FAILED"
    summary: str = ""

    def to_markdown(self) -> str:
        lines = []
        lines.append(f"# GPA Execution Report (v{self.workflow_version})")
        lines.append(f"\n**Mode**: {self.mode}")
        lines.append(f"**Started**: {self.start_time}")
        lines.append(f"**Finished**: {self.end_time}")
        lines.append(f"**Duration**: {self.total_duration_sec:.1f}s")
        lines.append(f"\n## Summary")
        lines.append(f"- Total steps: {self.steps_total}")
        lines.append(f"- Success: {self.steps_success}")
        lines.append(f"- Skipped: {self.steps_skipped}")
        lines.append(f"- Failed: {self.steps_failed}")
        lines.append(f"- Timeout: {self.steps_timeout}")
        lines.append(f"- **Final Status**: {self.final_status}")

        if self.steps_skipped > 0:
            lines.append(f"\n## Skipped Steps (User Notified)")
            for r in self.step_results:
                if r.status == "SKIPPED":
                    lines.append(f"- `{r.step_name}`: {r.skip_reason}")

        if self.steps_failed > 0:
            lines.append(f"\n## Failed Steps")
            for r in self.step_results:
                if r.status == "FAILED":
                    lines.append(f"- `{r.step_name}`: {r.error_message}")

        lines.append(f"\n## Step Details")
        lines.append("| Step | Status | Duration | Note |")
        lines.append("|------|--------|----------|------|")
        for r in self.step_results:
            note = r.skip_reason if r.status == "SKIPPED" else (r.error_message[:40] if r.status == "FAILED" else "")
            lines.append(f"| {r.step_name} | {r.status} | {r.duration_sec:.1f}s | {note} |")

        return "\n".join(lines)


# =============================================================================
# Workflow Runner
# =============================================================================

class WorkflowRunner:
    """
    GPA workflow execution engine.

    Modes:
        run: Execute workflow faithfully. Optional steps may auto-skip with
             mandatory user notification. Required steps cannot be skipped.
        optimize: Propose workflow modifications. Every change requires user
                  confirmation before applying.
    """

    def __init__(
        self,
        mode: str = "run",
        workflow: Optional[List[WorkflowStep]] = None,
        progress_callback: Optional[Callable[[str, float, str], None]] = None,
        notification_callback: Optional[Callable[[str], None]] = None,
    ):
        """
        Args:
            mode: "run" or "optimize"
            workflow: Custom workflow (optimize mode only). Defaults to STANDARD_WORKFLOW.
            progress_callback: Called with (step_name, progress_0_to_1, detail)
            notification_callback: Called with notification message (for skipped steps)
        """
        if mode not in ("run", "optimize"):
            raise ValueError(f"Invalid mode: {mode}. Use 'run' or 'optimize'.")

        self.mode = mode
        self.workflow = workflow or list(STANDARD_WORKFLOW)
        self.progress_callback = progress_callback
        self.notification_callback = notification_callback or self._default_notification
        self.step_results: List[StepResult] = []
        self.context: Dict[str, Any] = {}
        self._execution_started = False
        self._optimize_changes: List[Dict[str, Any]] = []  # Track proposed changes

    def _default_notification(self, message: str) -> None:
        """Default notification: print to stdout."""
        print(f"\n[WorkflowRunner] {message}\n")

    def _notify_skip(self, step: WorkflowStep, reason: str) -> None:
        """Notify user that a step was skipped. Mandatory in run mode."""
        msg = (
            f"Step '{step.name}' SKIPPED\n"
            f"  Reason: {reason}\n"
            f"  API: {step.api_name or 'N/A'} | Module: {step.module}\n"
            f"  This is an optional step (required=False)."
        )
        self.notification_callback(msg)

    def _notify_failure(self, step: WorkflowStep, error: str, action: FailureAction) -> None:
        """Notify user that a step failed and what action was taken."""
        msg = (
            f"Step '{step.name}' FAILED\n"
            f"  Error: {error[:200]}\n"
            f"  Action taken: {action.value}\n"
            f"  Required: {step.required}"
        )
        self.notification_callback(msg)

    # ------------------------------------------------------------------
    # Run Mode
    # ------------------------------------------------------------------

    async def execute(self, context: Dict[str, Any]) -> ExecutionReport:
        """
        Execute the workflow in run mode.

        Args:
            context: Execution context containing input data, config, flags.
                     Key fields: variants, tissue_profile, offline_mode,
                     spliceai_enabled, input_type, user_phenotypes, etc.

        Returns:
            ExecutionReport with complete execution history.
        """
        if self.mode != "run":
            raise RuntimeError("execute() can only be called in 'run' mode. Use optimize_mode methods for optimization.")

        self._execution_started = True
        self.context = context
        self.step_results = []

        start_ts = time.time()
        start_dt = datetime.now().isoformat()

        # Print workflow preamble
        print(f"\n{'='*70}")
        print(f"[WorkflowRunner] Starting GPA Analysis (v{WORKFLOW_VERSION})")
        print(f"[WorkflowRunner] Mode: RUN | Steps: {len(self.workflow)}")
        print(f"[WorkflowRunner] Required: {sum(1 for s in self.workflow if s.required)} | Optional: {sum(1 for s in self.workflow if not s.required)}")
        print(f"{'='*70}\n")

        # Print optional steps that may be skipped
        optional = [s for s in self.workflow if not s.required]
        if optional:
            print("[WorkflowRunner] Optional steps (may be auto-skipped with notification):")
            for s in optional:
                print(f"  - {s.name}: {s.skip_condition}")
            print()

        aborted = False
        for i, step in enumerate(self.workflow, 1):
            if aborted:
                # Record remaining steps as ABORTED
                self.step_results.append(StepResult(
                    step_name=step.name,
                    status="ABORTED",
                    start_time=time.time(),
                    end_time=time.time(),
                    duration_sec=0.0,
                ))
                continue

            print(f"[{i}/{len(self.workflow)}] Executing: {step.name} ...", end=" ", flush=True)

            step_start = time.time()

            # Check if step should be skipped
            should_skip, skip_reason = step.can_be_skipped(self.context)
            if should_skip:
                step_end = time.time()
                result = StepResult(
                    step_name=step.name,
                    status="SKIPPED",
                    start_time=step_start,
                    end_time=step_end,
                    duration_sec=step_end - step_start,
                    skip_reason=skip_reason,
                    notification_sent=True,
                )
                self.step_results.append(result)
                self._notify_skip(step, skip_reason)
                print(f"SKIPPED ({skip_reason})")
                continue

            # Execute step with timeout
            try:
                if step.timeout_sec > 0:
                    result_data = await asyncio.wait_for(
                        self._execute_step(step),
                        timeout=step.timeout_sec,
                    )
                else:
                    result_data = await self._execute_step(step)

                step_end = time.time()
                result = StepResult(
                    step_name=step.name,
                    status="SUCCESS",
                    start_time=step_start,
                    end_time=step_end,
                    duration_sec=step_end - step_start,
                    produced_data=result_data or {},
                )
                self.step_results.append(result)
                print(f"OK ({result.duration_sec:.1f}s)")

            except asyncio.TimeoutError:
                step_end = time.time()
                result = StepResult(
                    step_name=step.name,
                    status="TIMEOUT",
                    start_time=step_start,
                    end_time=step_end,
                    duration_sec=step_end - step_start,
                    error_message=f"Timeout after {step.timeout_sec}s",
                )
                self.step_results.append(result)
                self._notify_failure(step, f"Timeout after {step.timeout_sec}s", step.on_failure)
                print(f"TIMEOUT")

                if step.on_failure == FailureAction.ABORT:
                    aborted = True

            except (OSError, RuntimeError, ValueError, ImportError, AttributeError, TypeError, KeyError) as e:  # noqa: review-gate-allow — workflow engine last-resort guard: must catch all step failures to prevent pipeline crash
                step_end = time.time()
                tb = traceback.format_exc()
                result = StepResult(
                    step_name=step.name,
                    status="FAILED",
                    start_time=step_start,
                    end_time=step_end,
                    duration_sec=step_end - step_start,
                    error_message=f"{type(e).__name__}: {e}",
                    error_traceback=tb,
                )
                self.step_results.append(result)
                self._notify_failure(step, f"{type(e).__name__}: {e}", step.on_failure)
                print(f"FAILED ({type(e).__name__})")

                if step.on_failure == FailureAction.ABORT:
                    aborted = True

            # Update progress
            if self.progress_callback:
                progress = i / len(self.workflow)
                self.progress_callback(step.name, progress, f"Step {i}/{len(self.workflow)}")

        end_ts = time.time()
        end_dt = datetime.now().isoformat()

        # Build report
        report = self._build_report(start_dt, end_dt, end_ts - start_ts)

        print(f"\n{'='*70}")
        print(f"[WorkflowRunner] Execution Complete: {report.final_status}")
        print(f"[WorkflowRunner] Success: {report.steps_success} | Skipped: {report.steps_skipped} | Failed: {report.steps_failed} | Timeout: {report.steps_timeout}")
        print(f"[WorkflowRunner] Total Duration: {report.total_duration_sec:.1f}s")
        print(f"{'='*70}\n")

        return report

    async def _execute_step(self, step: WorkflowStep) -> Dict[str, Any]:
        """
        Execute a single workflow step by importing and calling the module/function.
        This is a dispatcher — actual business logic lives in the target modules.

        For now, this is a placeholder that returns empty dict.
        In production, it would dynamically import the module and call the function.
        """
        # Placeholder: in production, dynamically import and call
        # e.g.:
        #   module = importlib.import_module(step.module)
        #   func = getattr(module, step.function)
        #   return await func(self.context)
        await asyncio.sleep(0.01)  # Simulate minimal work
        return {}

    def _build_report(self, start_dt: str, end_dt: str, duration: float) -> ExecutionReport:
        """Build the final execution report."""
        success = sum(1 for r in self.step_results if r.status == "SUCCESS")
        skipped = sum(1 for r in self.step_results if r.status == "SKIPPED")
        failed = sum(1 for r in self.step_results if r.status == "FAILED")
        timeout = sum(1 for r in self.step_results if r.status == "TIMEOUT")

        if failed > 0 or timeout > 0:
            # Check if any required step failed
            failed_required = any(
                r.status in ("FAILED", "TIMEOUT")
                for r in self.step_results
                for s in self.workflow if s.name == r.step_name and s.required
            )
            final_status = "FAILED" if failed_required else "PARTIAL"
        else:
            final_status = "SUCCESS"

        return ExecutionReport(
            workflow_version=WORKFLOW_VERSION,
            mode=self.mode,
            start_time=start_dt,
            end_time=end_dt,
            total_duration_sec=duration,
            steps_total=len(self.workflow),
            steps_success=success,
            steps_skipped=skipped,
            steps_failed=failed,
            steps_timeout=timeout,
            step_results=self.step_results,
            final_status=final_status,
        )

    # ------------------------------------------------------------------
    # Optimize Mode
    # ------------------------------------------------------------------

    async def propose_workflow_change(
        self,
        description: str,
        old_steps: Optional[List[WorkflowStep]] = None,
        new_steps: Optional[List[WorkflowStep]] = None,
    ) -> bool:
        """
        Propose a workflow modification. Prints diff and waits for user confirmation.

        Args:
            description: Human-readable description of the change
            old_steps: Steps being replaced (None = append)
            new_steps: New steps to insert/replace

        Returns:
            True if user confirmed, False if rejected.
        """
        if self.mode != "optimize":
            raise RuntimeError("propose_workflow_change() can only be called in 'optimize' mode.")

        print(f"\n{'='*70}")
        print(f"[WorkflowRunner OPTIMIZE] Proposed Workflow Change")
        print(f"{'='*70}")
        print(f"\nDescription: {description}\n")

        if old_steps and new_steps:
            print("REPLACING steps:")
            for s in old_steps:
                print(f"  - {s.name} (module={s.module}, api={s.api_name or 'N/A'})")
            print("WITH:")
            for s in new_steps:
                print(f"  + {s.name} (module={s.module}, api={s.api_name or 'N/A'}, required={s.required})")
        elif new_steps:
            print("ADDING steps:")
            for s in new_steps:
                print(f"  + {s.name} (module={s.module}, api={s.api_name or 'N/A'}, required={s.required})")
        elif old_steps:
            print("REMOVING steps:")
            for s in old_steps:
                print(f"  - {s.name} (module={s.module}, api={s.api_name or 'N/A'})")

        print(f"\n{'='*70}")

        # In optimize mode, we cannot use input() in async context easily.
        # Return False to indicate "needs user confirmation" and let the caller handle it.
        self._optimize_changes.append({
            "description": description,
            "old_steps": old_steps,
            "new_steps": new_steps,
            "status": "PENDING",
            "proposed_at": datetime.now().isoformat(),
        })

        print("[WorkflowRunner] This change is PENDING user confirmation.")
        print("[WorkflowRunner] Call confirm_last_change() to apply, or reject_last_change() to discard.\n")
        return False  # Pending confirmation

    def confirm_last_change(self) -> bool:
        """Apply the last proposed workflow change."""
        if not self._optimize_changes:
            print("[WorkflowRunner] No pending changes to confirm.")
            return False

        change = self._optimize_changes[-1]
        old_steps = change.get("old_steps", [])
        new_steps = change.get("new_steps", [])

        # Apply the change
        if old_steps and new_steps:
            # Replace
            for old in old_steps:
                for i, existing in enumerate(self.workflow):
                    if existing.name == old.name:
                        self.workflow[i] = new_steps[0]  # Simplified: 1-for-1
                        break
        elif new_steps:
            # Append
            self.workflow.extend(new_steps)
        elif old_steps:
            # Remove
            self.workflow = [s for s in self.workflow if s.name not in {o.name for o in old_steps}]

        change["status"] = "CONFIRMED"
        print(f"[WorkflowRunner] Change APPLIED: {change['description']}")
        return True

    def reject_last_change(self) -> bool:
        """Reject the last proposed workflow change."""
        if not self._optimize_changes:
            print("[WorkflowRunner] No pending changes to reject.")
            return False

        change = self._optimize_changes.pop()
        print(f"[WorkflowRunner] Change REJECTED: {change['description']}")
        return True

    def get_pending_changes(self) -> List[Dict[str, Any]]:
        """Return all pending workflow changes awaiting confirmation."""
        return self._optimize_changes

    def reset_workflow(self) -> None:
        """Reset workflow to the canonical STANDARD_WORKFLOW."""
        self.workflow = list(STANDARD_WORKFLOW)
        self._optimize_changes = []
        print("[WorkflowRunner] Workflow reset to canonical STANDARD_WORKFLOW.")


# =============================================================================
# Convenience Functions
# =============================================================================

async def run_workflow(context: Dict[str, Any], **kwargs) -> ExecutionReport:
    """
    One-shot workflow execution in run mode.

    Args:
        context: Execution context
        **kwargs: Passed to WorkflowRunner constructor

    Returns:
        ExecutionReport
    """
    runner = WorkflowRunner(mode="run", **kwargs)
    return await runner.execute(context)


def print_workflow_summary() -> None:
    """Print a human-readable summary of the standard workflow."""
    print(f"\nGPA Standard Workflow (v{WORKFLOW_VERSION})\n")
    print(f"{'Step':<30} {'Module':<25} {'API':<15} {'Required':<10} {'Timeout':<10}")
    print("-" * 95)
    for step in STANDARD_WORKFLOW:
        api = step.api_name or "-"
        req = "YES" if step.required else "NO"
        to = f"{step.timeout_sec}s" if step.timeout_sec > 0 else "-"
        print(f"{step.name:<30} {step.module:<25} {api:<15} {req:<10} {to:<10}")
    print()


if __name__ == "__main__":
    print_workflow_summary()
