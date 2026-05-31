#!/usr/bin/env python3
"""
GPA Workflow PM (Progress Monitor) (v0.11.0)
=============================================
Project Manager-style progress reporting for long-running GPA analysis tasks.

Reports progress:
  - Every 3 minutes (configurable)
  - Upon completion of each sub-task
  - When significant milestones are reached (e.g., 25%, 50%, 75%)

Usage:
    from gpa_workflow_pm import WorkflowPM
    pm = WorkflowPM(report_interval_sec=180)
    pm.start_phase("vep_annotation", total_items=39681)
    pm.update_progress(processed=1000)
    ...
    pm.complete_phase()
"""

import time
import threading
from dataclasses import dataclass, field
from datetime import datetime
from typing import Optional, Dict, Any, List, Callable


# =============================================================================
# Phase Tracking
# =============================================================================

@dataclass
class PhaseInfo:
    """Information about a running phase."""
    name: str
    start_time: float
    total_items: int = 0
    processed_items: int = 0
    current_item_detail: str = ""
    status: str = "running"  # running | paused | completed | failed
    sub_tasks: List[Dict[str, Any]] = field(default_factory=list)

    @property
    def progress_pct(self) -> float:
        if self.total_items <= 0:
            return 0.0
        return min(100.0, (self.processed_items / self.total_items) * 100)

    @property
    def elapsed_sec(self) -> float:
        return time.time() - self.start_time

    def eta_sec(self) -> Optional[int]:
        if self.total_items <= 0 or self.processed_items <= 0:
            return None
        rate = self.processed_items / self.elapsed_sec
        remaining = self.total_items - self.processed_items
        if rate > 0:
            return int(remaining / rate)
        return None


# =============================================================================
# Workflow PM
# =============================================================================

class WorkflowPM:
    """
    Progress monitor that reports GPA analysis status to the user.

    Designed to be non-intrusive: reports are printed to stdout and can be
    captured by the agent to relay to the user.
    """

    def __init__(
        self,
        report_interval_sec: int = 180,  # 3 minutes
        milestone_intervals: Optional[List[float]] = None,
        report_callback: Optional[Callable[[str], None]] = None,
    ):
        """
        Args:
            report_interval_sec: Time between automatic progress reports (seconds)
            milestone_intervals: Progress percentages at which to report [0.0-1.0]
            report_callback: Custom callback for reports (defaults to print)
        """
        self.report_interval_sec = report_interval_sec
        self.milestone_intervals = milestone_intervals or [0.25, 0.50, 0.75, 1.0]
        self.report_callback = report_callback or self._default_report

        self.phases: Dict[str, PhaseInfo] = {}
        self.current_phase: Optional[str] = None
        self._overall_start_time: Optional[float] = None
        self._last_report_time: float = 0
        self._reported_milestones: set = set()
        self._timer: Optional[threading.Timer] = None
        self._lock = threading.Lock()

    def _default_report(self, message: str) -> None:
        """Default report output: print with timestamp."""
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"\n[PM {ts}] {message}\n")

    def _format_duration(self, seconds: float) -> str:
        """Format duration as human-readable string."""
        if seconds < 60:
            return f"{seconds:.0f}s"
        elif seconds < 3600:
            return f"{seconds/60:.1f}min"
        else:
            return f"{seconds/3600:.1f}h"

    def _format_eta(self, seconds: Optional[int]) -> str:
        """Format ETA as human-readable string."""
        if seconds is None:
            return "calculating..."
        if seconds < 60:
            return f"{seconds}s"
        elif seconds < 3600:
            return f"{seconds//60}min {seconds%60}s"
        else:
            return f"{seconds//3600}h {(seconds%3600)//60}min"

    # ------------------------------------------------------------------
    # Lifecycle
    # ------------------------------------------------------------------

    def start_analysis(self) -> None:
        """Mark the start of the overall analysis."""
        self._overall_start_time = time.time()
        self._last_report_time = time.time()
        self._reported_milestones = set()
        self._start_timer()
        self.report_callback(
            "Analysis started. I will report progress every 3 minutes or when each phase completes."
        )

    def end_analysis(self, status: str = "SUCCESS") -> None:
        """Mark the end of the overall analysis."""
        self._stop_timer()
        if self._overall_start_time:
            duration = time.time() - self._overall_start_time
            self.report_callback(
                f"Analysis complete ({status}). Total duration: {self._format_duration(duration)}."
            )

    def _start_timer(self) -> None:
        """Start the periodic reporting timer."""
        self._timer = threading.Timer(self.report_interval_sec, self._timer_callback)
        self._timer.daemon = True
        self._timer.start()

    def _stop_timer(self) -> None:
        """Stop the periodic reporting timer."""
        if self._timer:
            self._timer.cancel()
            self._timer = None

    def _timer_callback(self) -> None:
        """Called by timer to trigger periodic reports."""
        with self._lock:
            now = time.time()
            if now - self._last_report_time >= self.report_interval_sec:
                self._send_periodic_report()
                self._last_report_time = now
            # Restart timer
            self._start_timer()

    # ------------------------------------------------------------------
    # Phase Management
    # ------------------------------------------------------------------

    def start_phase(self, phase_name: str, total_items: int = 0, detail: str = "") -> None:
        """
        Start tracking a new phase.

        Args:
            phase_name: Name of the phase (e.g., "vep_annotation")
            total_items: Total number of items to process (e.g., variant count)
            detail: Optional description
        """
        with self._lock:
            if phase_name in self.phases:
                # Restart existing phase
                self.phases[phase_name].start_time = time.time()
                self.phases[phase_name].processed_items = 0
                self.phases[phase_name].status = "running"
            else:
                self.phases[phase_name] = PhaseInfo(
                    name=phase_name,
                    start_time=time.time(),
                    total_items=total_items,
                    current_item_detail=detail,
                )
            self.current_phase = phase_name

        self.report_callback(
            f"Starting phase: {phase_name}"
            + (f" ({total_items} items)" if total_items > 0 else "")
            + (f" — {detail}" if detail else "")
        )

    def update_progress(
        self,
        processed: int,
        total: Optional[int] = None,
        detail: str = "",
    ) -> None:
        """
        Update progress for the current phase.

        Args:
            processed: Number of items processed so far
            total: Override total items (optional)
            detail: Current item being processed (e.g., "chr1:12345")
        """
        with self._lock:
            if not self.current_phase:
                return
            phase = self.phases[self.current_phase]
            phase.processed_items = processed
            if total is not None:
                phase.total_items = total
            if detail:
                phase.current_item_detail = detail

            # Check milestones
            if phase.total_items > 0:
                pct = phase.progress_pct / 100.0
                for m in self.milestone_intervals:
                    if pct >= m and m not in self._reported_milestones:
                        self._reported_milestones.add(m)
                        self._send_milestone_report(phase, m)
                        break

    def increment_progress(self, amount: int = 1, detail: str = "") -> None:
        """Increment processed count by amount."""
        with self._lock:
            if not self.current_phase:
                return
            phase = self.phases[self.current_phase]
            phase.processed_items += amount
            if detail:
                phase.current_item_detail = detail

            # Check milestones
            if phase.total_items > 0:
                pct = phase.progress_pct / 100.0
                for m in self.milestone_intervals:
                    if pct >= m and m not in self._reported_milestones:
                        self._reported_milestones.add(m)
                        self._send_milestone_report(phase, m)
                        break

    def complete_phase(self, detail: str = "") -> None:
        """Mark the current phase as complete and report."""
        with self._lock:
            if not self.current_phase:
                return
            phase = self.phases[self.current_phase]
            phase.status = "completed"
            phase.processed_items = phase.total_items  # Force 100%
            duration = phase.elapsed_sec

        self.report_callback(
            f"Phase complete: {phase.name} — "
            f"{phase.processed_items}/{phase.total_items} items "
            f"in {self._format_duration(duration)}"
            + (f" — {detail}" if detail else "")
        )

    def fail_phase(self, error: str) -> None:
        """Mark the current phase as failed."""
        with self._lock:
            if not self.current_phase:
                return
            phase = self.phases[self.current_phase]
            phase.status = "failed"

        self.report_callback(
            f"Phase FAILED: {phase.name} — {error}"
        )

    def add_sub_task(self, sub_task_name: str, status: str = "running", detail: str = "") -> None:
        """Add a sub-task to the current phase."""
        with self._lock:
            if not self.current_phase:
                return
            self.phases[self.current_phase].sub_tasks.append({
                "name": sub_task_name,
                "status": status,
                "detail": detail,
                "timestamp": time.time(),
            })

    def complete_sub_task(self, sub_task_name: str, detail: str = "") -> None:
        """Mark a sub-task as complete."""
        with self._lock:
            if not self.current_phase:
                return
            phase = self.phases[self.current_phase]
            for st in phase.sub_tasks:
                if st["name"] == sub_task_name and st["status"] == "running":
                    st["status"] = "completed"
                    st["detail"] = detail
                    break

        self.report_callback(
            f"Sub-task complete: {sub_task_name}"
            + (f" — {detail}" if detail else "")
        )

    # ------------------------------------------------------------------
    # Reporting
    # ------------------------------------------------------------------

    def _send_periodic_report(self) -> None:
        """Send a periodic progress report (every 3 minutes)."""
        with self._lock:
            if not self.current_phase:
                return
            phase = self.phases[self.current_phase]
            elapsed = phase.elapsed_sec
            pct = phase.progress_pct
            eta = phase.eta_sec()
            detail = phase.current_item_detail

        msg = (
            f"Progress update: {phase.name}\n"
            f"  Progress: {pct:.1f}% ({phase.processed_items}/{phase.total_items})\n"
            f"  Elapsed: {self._format_duration(elapsed)}\n"
            f"  ETA: {self._format_eta(eta)}"
        )
        if detail:
            msg += f"\n  Current: {detail}"

        # Add overall progress if available
        if self._overall_start_time:
            overall_elapsed = time.time() - self._overall_start_time
            msg += f"\n  Overall elapsed: {self._format_duration(overall_elapsed)}"

        self.report_callback(msg)

    def _send_milestone_report(self, phase: PhaseInfo, milestone: float) -> None:
        """Send a milestone report (e.g., 25%, 50%)."""
        pct = int(milestone * 100)
        eta = phase.eta_sec()

        msg = (
            f"Milestone: {phase.name} {pct}% complete\n"
            f"  {phase.processed_items}/{phase.total_items} items processed\n"
            f"  ETA: {self._format_eta(eta)}"
        )
        self.report_callback(msg)

    def get_current_status(self) -> Dict[str, Any]:
        """Get current status as a dict (for programmatic use)."""
        with self._lock:
            if not self.current_phase:
                return {"status": "idle", "phases": {}}

            phase = self.phases[self.current_phase]
            return {
                "current_phase": phase.name,
                "status": phase.status,
                "progress_pct": phase.progress_pct,
                "processed": phase.processed_items,
                "total": phase.total_items,
                "elapsed_sec": phase.elapsed_sec,
                "eta_sec": phase.eta_sec(),
                "detail": phase.current_item_detail,
                "sub_tasks": phase.sub_tasks,
                "all_phases": {
                    name: {
                        "status": p.status,
                        "progress_pct": p.progress_pct,
                        "elapsed_sec": p.elapsed_sec,
                    }
                    for name, p in self.phases.items()
                },
            }

    def get_summary_report(self) -> str:
        """Get a human-readable summary of all phases."""
        with self._lock:
            lines = []
            lines.append("Analysis Progress Summary")
            lines.append("=" * 50)

            for name, phase in self.phases.items():
                status_icon = {
                    "running": "▶",
                    "completed": "✓",
                    "failed": "✗",
                    "paused": "⏸",
                }.get(phase.status, "?")

                if phase.total_items > 0:
                    progress = f"{phase.progress_pct:.1f}% ({phase.processed_items}/{phase.total_items})"
                else:
                    progress = phase.status

                lines.append(
                    f"{status_icon} {name:<30} {progress:<20} {self._format_duration(phase.elapsed_sec)}"
                )

            if self._overall_start_time:
                total = time.time() - self._overall_start_time
                lines.append("-" * 50)
                lines.append(f"Total elapsed: {self._format_duration(total)}")

            return "\n".join(lines)


# =============================================================================
# Convenience: GPA-Specific PM Wrapper
# =============================================================================

class GPAPM(WorkflowPM):
    """
    GPA-specific PM with pre-configured phases for standard pipeline.
    """

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self._phase_map = {
            "preflight_check": ("Preflight Check", 0),
            "input_parsing": ("Input Parsing", 0),
            "vep_annotation": ("VEP Annotation", 0),
            "transcript_selection": ("Transcript Selection", 0),
            "ensembl_gene_batch": ("Ensembl Gene Query", 0),
            "uniprot_gene_batch": ("UniProt Gene Query", 0),
            "hgnc_gene_batch": ("HGNC Gene Query", 0),
            "gnomad_constraint_batch": ("gnomAD Constraint Query", 0),
            "gtex_expression": ("GTEx Expression Query", 0),
            "myvariant_enrichment": ("MyVariant.info Enrichment", 0),
            "gnomad_variant_frequency": ("gnomAD Variant Frequency", 0),
            "clinvar_ncbi": ("ClinVar NCBI Query", 0),
            "hgnc_normalization": ("HGNC Normalization", 0),
            "gene_constraint_population": ("Gene Constraint", 0),
            "nmd_prediction": ("NMD Prediction", 0),
            "transcript_correction": ("Transcript Correction", 0),
            "vep_reannotation": ("VEP Reannotation", 0),
            "tissue_relevance": ("Tissue Relevance", 0),
            "spliceai_prediction": ("SpliceAI Prediction", 0),
            "tier_classification": ("Tier Classification", 0),
            "multi_hit_detection": ("Multi-Hit Detection", 0),
            "phenotype_matching": ("Phenotype Matching", 0),
            "qc_checks": ("QC Checks", 0),
            "report_generation": ("Report Generation", 0),
        }

    def start_step(self, step_name: str, total_items: int = 0, detail: str = "") -> None:
        """Start a workflow step with friendly naming."""
        friendly_name, _ = self._phase_map.get(step_name, (step_name, 0))
        self.start_phase(friendly_name, total_items=total_items, detail=detail)

    def complete_step(self, step_name: str, detail: str = "") -> None:
        """Complete a workflow step."""
        friendly_name, _ = self._phase_map.get(step_name, (step_name, 0))
        # Find the phase by friendly name
        for phase_name, phase in self.phases.items():
            if phase_name == friendly_name:
                self.current_phase = friendly_name
                self.complete_phase(detail)
                return
        # If not found, just report
        self.report_callback(f"Step complete: {friendly_name}" + (f" — {detail}" if detail else ""))


# =============================================================================
# Demo
# =============================================================================

if __name__ == "__main__":
    import time as _time

    pm = GPAPM(report_interval_sec=5)  # 5 seconds for demo
    pm.start_analysis()

    pm.start_step("vep_annotation", total_items=100, detail="Batch 1/20")
    for i in range(101):
        pm.update_progress(processed=i, detail=f"Variant {i}/100")
        _time.sleep(0.05)
    pm.complete_step("vep_annotation")

    pm.start_step("tier_classification", total_items=50)
    for i in range(51):
        pm.update_progress(processed=i)
        _time.sleep(0.02)
    pm.complete_step("tier_classification")

    pm.end_analysis()
    print(pm.get_summary_report())
