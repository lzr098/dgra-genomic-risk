#!/usr/bin/env python3
"""
GPA Workflow Engine (v0.11.1)
==============================
Minimal phase-based executor for the GPA analysis pipeline.

- Runs a list of named async phase functions in order.
- Saves a JSON checkpoint after each phase completes.
- Resumes from the latest checkpoint on restart.
- Cleans all phase checkpoints after successful completion.

This is intentionally thin: phases live in gpa_pipeline.py and are just
async def f(context: dict) -> dict. The engine only orchestrates,
serializes, and recovers.
"""

import asyncio
import json
import gzip
from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import Any, Awaitable, Callable, Dict, List, Optional, Tuple

try:
    from dgra_core import GPAConfig, Variant, Evidence
except Exception:  # pragma: no cover - engine can be imported before dgra_core in odd paths
    GPAConfig = None  # type: ignore[misc,assignment]
    Variant = None  # type: ignore[misc,assignment]
    Evidence = None  # type: ignore[misc,assignment]

try:
    from gpa_audit_trail import AuditTrail
except Exception:  # pragma: no cover
    AuditTrail = None  # type: ignore[misc,assignment]


def _variants_from_dicts(variants_data: List[Dict[str, Any]]) -> List[Any]:
    """Rehydrate Variant dataclass list from checkpoint primitives."""
    if Variant is None:
        return variants_data
    result: List[Any] = []
    for d in variants_data:
        d = dict(d)
        evidence_list = d.get("evidence_chain", [])
        if Evidence is not None and isinstance(evidence_list, list):
            d["evidence_chain"] = [Evidence(**e) for e in evidence_list]
        result.append(Variant(**d))
    return result


# Types
PhaseFn = Callable[[Dict[str, Any]], Awaitable[Dict[str, Any]]]


class WorkflowError(Exception):
    pass


class WorkflowEngine:
    """
    Execute phases sequentially with checkpoint/resume support.

    Args:
        phases: ordered list of (phase_name, phase_fn).
        checkpoint_dir: directory to store phase JSON checkpoints.
        resume: if True, load existing checkpoints instead of re-running.
        keep_checkpoints: if True, do NOT delete checkpoints on success.
    """

    def __init__(
        self,
        phases: List[Tuple[str, PhaseFn]],
        checkpoint_dir: Optional[Path] = None,
        resume: bool = False,
        keep_checkpoints: bool = False,
        audit_trail: Optional[Any] = None,
    ):
        self.phases = phases
        self.checkpoint_dir = checkpoint_dir or Path.cwd()
        self.resume = resume
        self.keep_checkpoints = keep_checkpoints
        self.audit_trail = audit_trail
        self._completed: List[str] = []

    @staticmethod
    def _checkpoint_path(checkpoint_dir: Path, phase_name: str, gz: bool = True) -> Path:
        ext = "json.gz" if gz else "json"
        return checkpoint_dir / f"gpa_checkpoint_{phase_name}.{ext}"

    def _serialize_context(self, context: Dict[str, Any]) -> Dict[str, Any]:
        """Convert context values to JSON-serializable primitives."""
        data: Dict[str, Any] = {}
        for k, v in context.items():
            if k in ("proxy_route_map", "global_config"):
                # These are reconstructed on demand; do not checkpoint.
                continue
            if GPAConfig is not None and isinstance(v, GPAConfig):
                data[k] = asdict(v)
            elif Variant is not None and isinstance(v, list) and v and isinstance(v[0], Variant):
                data[k] = [asdict(x) for x in v]
            elif isinstance(v, Path):
                data[k] = str(v)
            else:
                try:
                    json.dumps(v, ensure_ascii=False, default=str)
                    data[k] = v
                except (TypeError, ValueError):
                    # Skip unserializable runtime objects.
                    continue
        return data

    @staticmethod
    def _deserialize_context(data: Dict[str, Any]) -> Dict[str, Any]:
        """Restore GPAConfig / Variant objects from checkpoint primitives."""
        ctx = dict(data)
        if GPAConfig is not None and isinstance(ctx.get("config"), dict):
            ctx["config"] = GPAConfig(**ctx["config"])
        if Variant is not None and isinstance(ctx.get("variants"), list):
            ctx["variants"] = _variants_from_dicts(ctx["variants"])
        return ctx

    def _load_checkpoint(self, phase_name: str) -> Optional[Dict[str, Any]]:
        path = self._checkpoint_path(self.checkpoint_dir, phase_name, gz=True)
        if not path.exists():
            path = self._checkpoint_path(self.checkpoint_dir, phase_name, gz=False)
        if not path.exists():
            return None
        try:
            if path.suffix == ".gz":
                with gzip.open(path, "rt", encoding="utf-8") as f:
                    data = json.load(f)
            else:
                with open(path, "r", encoding="utf-8") as f:
                    data = json.load(f)
            ctx = self._deserialize_context(data)
            print(f"[WorkflowEngine] Resumed phase '{phase_name}' from {path}")
            return ctx
        except Exception as e:
            print(f"[WorkflowEngine] Failed to load checkpoint {path}: {e}")
            return None

    def _save_checkpoint(self, phase_name: str, context: Dict[str, Any]) -> Path:
        self.checkpoint_dir.mkdir(parents=True, exist_ok=True)
        path = self._checkpoint_path(self.checkpoint_dir, phase_name, gz=True)
        data = self._serialize_context(context)
        try:
            with gzip.open(path, "wt", encoding="utf-8") as f:
                json.dump(data, f, ensure_ascii=False, default=str)
        except Exception as e:
            print(f"[WorkflowEngine] Failed to write checkpoint {path}: {e}")
            raise
        return path

    def _cleanup_checkpoints(self) -> None:
        if self.keep_checkpoints:
            return
        removed = 0
        for phase_name, _ in self.phases:
            for gz in (True, False):
                path = self._checkpoint_path(self.checkpoint_dir, phase_name, gz=gz)
                if path.exists():
                    path.unlink(missing_ok=True)
                    removed += 1
        if removed:
            print(f"[WorkflowEngine] Cleaned {removed} checkpoint files")

    async def run(self, initial_context: Dict[str, Any]) -> Dict[str, Any]:
        context = dict(initial_context)
        context.setdefault("_workflow_meta", {})["started_at"] = datetime.now().isoformat()
        context.setdefault("_phase_status", {})
        if self.audit_trail:
            self.audit_trail.record_metric("input_variants", len(initial_context.get("variants_data", [])))
            self.audit_trail.record_metric("resume", self.resume)

        for phase_name, phase_fn in self.phases:
            status = context.get("_phase_status", {}).get(phase_name)
            if status == "completed" and self.resume:
                print(f"[WorkflowEngine] Phase '{phase_name}' already completed; skipping")
                self._completed.append(phase_name)
                if self.audit_trail:
                    self.audit_trail.record_phase_start(phase_name, {"source": "checkpoint_resume"})
                    self.audit_trail.record_phase_end(phase_name, status="resumed")
                continue

            # Try to resume from checkpoint if resuming and not already completed in context
            if self.resume:
                cached = self._load_checkpoint(phase_name)
                if cached is not None:
                    context = cached
                    context["_phase_status"][phase_name] = "completed"
                    self._completed.append(phase_name)
                    if self.audit_trail:
                        self.audit_trail.record_phase_start(phase_name, {"source": "checkpoint_resume"})
                        self.audit_trail.record_phase_end(phase_name, status="resumed")
                    continue

            print(f"[WorkflowEngine] Running phase: {phase_name}")
            if self.audit_trail:
                self.audit_trail.record_phase_start(phase_name)
            try:
                context = await phase_fn(context)
            except WorkflowError as we:
                if self.audit_trail:
                    self.audit_trail.record_phase_end(phase_name, status="failed", error=str(we))
                raise
            except Exception as e:
                # Save a failure checkpoint for inspection before re-raising
                failure_ckpt = self._checkpoint_path(self.checkpoint_dir, f"{phase_name}_FAILED", gz=True)
                try:
                    with gzip.open(failure_ckpt, "wt", encoding="utf-8") as f:
                        json.dump(context, f, ensure_ascii=False, default=str)
                except Exception:
                    pass
                if self.audit_trail:
                    self.audit_trail.record_phase_end(phase_name, status="failed", error=str(e))
                raise WorkflowError(f"Phase '{phase_name}' failed: {e}") from e

            context.setdefault("_phase_status", {})[phase_name] = "completed"
            self._completed.append(phase_name)
            self._save_checkpoint(phase_name, context)
            if self.audit_trail:
                self.audit_trail.record_phase_end(phase_name, status="success")

        context.setdefault("_workflow_meta", {})["completed_at"] = datetime.now().isoformat()
        self._cleanup_checkpoints()
        if self.audit_trail:
            written = self.audit_trail.write()
            if written:
                context.setdefault("_workflow_meta", {})["trace_path"] = str(written[0])
                context.setdefault("_workflow_meta", {})["audit_log_path"] = str(written[1])
        return context


def mark_phase_completed(context: Dict[str, Any], phase_name: str) -> None:
    """Helper for phases to self-report completion."""
    context.setdefault("_phase_status", {})[phase_name] = "completed"
