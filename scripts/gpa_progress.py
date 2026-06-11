#!/usr/bin/env python3
"""
GPA Two-Phase Progress Tracker (v0.10.16)

Provides fine-grained, machine-readable progress logging for the GPA
two-phase pipeline so that long-running analyses can be monitored in
real time.

The tracker writes a structured JSON log file. Each line is a JSON object
with at least:
  - timestamp: ISO-8601 UTC time
  - event: one of phase_start, phase_end, step_start, step_end, progress,
           api_call, warning, error, summary
  - phase: current pipeline phase (e.g., "phase0", "phase1", "phase2")
  - step: current step within the phase (optional)
  - message: human-readable description
  - data: event-specific metrics (counts, timings, etc.)

Usage:
    tracker = GPAProgressTracker("/path/to/gpa_progress.jsonl")
    tracker.phase_start("phase1", "Fast local triage", {"total_variants": 1500000})
    tracker.step_progress("phase1", "filter", {"processed": 1000, "total": 1500000})
    tracker.phase_end("phase1", {"candidates": 500, "tier1": 50, "tier2": 450})
"""

import json
import time
import os
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional


class GPAProgressTracker:
    """Fine-grained progress tracker for GPA two-phase pipeline."""

    def __init__(self, log_path: Optional[str] = None, enabled: bool = True):
        """Initialize tracker.

        Args:
            log_path: Path to JSON Lines log file. If None, progress is only
                printed to stdout.
            enabled: If False, all methods become no-ops.
        """
        self.enabled = enabled
        self.log_path = Path(log_path) if log_path else None
        self._start_times: Dict[str, float] = {}
        self._phase_start_times: Dict[str, float] = {}
        if self.log_path:
            self.log_path.parent.mkdir(parents=True, exist_ok=True)

    def _emit(self, event: str, phase: str, step: str, message: str, data: Dict[str, Any]) -> None:
        if not self.enabled:
            return

        record = {
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "event": event,
            "phase": phase,
            "step": step,
            "message": message,
            "data": data,
        }

        # Always print key events to stdout for immediate visibility
        if event in ("phase_start", "phase_end", "step_start", "step_end", "warning", "error", "summary"):
            print(f"[GPA Progress] {phase}/{step}: {message}")
        elif event == "progress":
            # Throttle progress prints: only when pct changes by >= 10 or finished
            pct = data.get("percent")
            if pct is not None and (pct >= 100 or int(pct) % 10 == 0):
                print(f"[GPA Progress] {phase}/{step}: {message}")

        if self.log_path:
            try:
                with open(self.log_path, "a", encoding="utf-8") as f:
                    f.write(json.dumps(record, ensure_ascii=False, default=str) + "\n")
            except Exception as e:
                print(f"[GPA Progress] WARNING: failed to write progress log: {e}")

    def phase_start(self, phase: str, message: str = "", data: Optional[Dict[str, Any]] = None) -> None:
        self._phase_start_times[phase] = time.time()
        self._emit("phase_start", phase, "", message or f"{phase} started", data or {})

    def phase_end(self, phase: str, message: str = "", data: Optional[Dict[str, Any]] = None) -> None:
        elapsed = time.time() - self._phase_start_times.get(phase, time.time())
        payload = {"elapsed_seconds": round(elapsed, 2)}
        if data:
            payload.update(data)
        self._emit("phase_end", phase, "", message or f"{phase} completed", payload)

    def step_start(self, phase: str, step: str, message: str = "", data: Optional[Dict[str, Any]] = None) -> None:
        key = f"{phase}:{step}"
        self._start_times[key] = time.time()
        self._emit("step_start", phase, step, message or f"{step} started", data or {})

    def step_end(self, phase: str, step: str, message: str = "", data: Optional[Dict[str, Any]] = None) -> None:
        key = f"{phase}:{step}"
        elapsed = time.time() - self._start_times.get(key, time.time())
        payload = {"elapsed_seconds": round(elapsed, 2)}
        if data:
            payload.update(data)
        self._emit("step_end", phase, step, message or f"{step} completed", payload)

    def step_progress(self, phase: str, step: str, processed: int, total: int, message: str = "") -> None:
        pct = (processed / total * 100) if total > 0 else 0
        data = {
            "processed": processed,
            "total": total,
            "percent": round(pct, 2),
            "remaining": total - processed,
        }
        self._emit("progress", phase, step, message or f"{processed}/{total} ({pct:.1f}%)", data)

    def api_call(self, phase: str, step: str, api_name: str, status: str, detail: Optional[Dict[str, Any]] = None) -> None:
        data = {"api": api_name, "status": status}
        if detail:
            data.update(detail)
        self._emit("api_call", phase, step, f"{api_name}: {status}", data)

    def warning(self, phase: str, step: str, message: str, data: Optional[Dict[str, Any]] = None) -> None:
        self._emit("warning", phase, step, message, data or {})

    def error(self, phase: str, step: str, message: str, data: Optional[Dict[str, Any]] = None) -> None:
        self._emit("error", phase, step, message, data or {})

    def summary(self, message: str, data: Dict[str, Any]) -> None:
        self._emit("summary", "pipeline", "", message, data)

    def finish(self, message: str = "Pipeline complete", data: Optional[Dict[str, Any]] = None) -> None:
        """Emit a final summary event marking the pipeline as finished.

        This is a convenience method for callers that expect a tracker.finish()
        call at the end of a run. It simply records a summary event; the log
        file remains open-append and does not need to be closed.
        """
        self.summary(message, data or {})
