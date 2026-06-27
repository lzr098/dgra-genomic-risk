#!/usr/bin/env python3
"""
GPA Audit Trail (v0.11.2)
=========================
Minimal evidence chain for the GPA analysis pipeline.

Records:
  - Workflow phase timing and status
  - External API calls (name, URL, status, cache hit, latency)
  - Key classification decisions (tier, reason, subject variant/gene)
  - Pipeline metrics (variant counts, API coverage)

Outputs:
  - trace.json: machine-readable structured record
  - audit.log: human-readable line-oriented log

This is intentionally separate from the workflow engine so it can be
enabled only when requested and reused by other skills.
"""

import json
import time
from datetime import datetime
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple


class AuditTrail:
    """
    Lightweight audit trail for a single analysis run.

    Args:
        output_dir: directory where trace.json and audit.log will be written.
        run_id: optional run identifier; defaults to a timestamp.
        enabled: if False, all record_* calls become no-ops and write() returns None.
    """

    def __init__(
        self,
        output_dir: Optional[Path] = None,
        run_id: Optional[str] = None,
        enabled: bool = True,
    ):
        self.enabled = enabled
        self.output_dir = Path(output_dir or Path.cwd())
        self.run_id = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        self.started_at = datetime.now().isoformat()
        self.records: List[Dict[str, Any]] = []
        self.phase_stack: Dict[str, float] = {}
        self.metrics: Dict[str, Any] = {}
        self._api_call_count = 0
        self._api_cache_hits = 0
        self._api_errors = 0

    # ------------------------------------------------------------------
    # Recording
    # ------------------------------------------------------------------
    def record_phase_start(self, phase_name: str, meta: Optional[Dict[str, Any]] = None) -> None:
        if not self.enabled:
            return
        self.phase_stack[phase_name] = time.time()
        self._append(
            "phase_start",
            {
                "phase": phase_name,
                "meta": meta or {},
            },
        )

    def record_phase_end(
        self,
        phase_name: str,
        status: str = "success",
        error: Optional[str] = None,
        meta: Optional[Dict[str, Any]] = None,
    ) -> None:
        if not self.enabled:
            return
        start_ts = self.phase_stack.pop(phase_name, None)
        duration_ms = None
        if start_ts is not None:
            duration_ms = round((time.time() - start_ts) * 1000, 2)
        payload: Dict[str, Any] = {
            "phase": phase_name,
            "status": status,
            "duration_ms": duration_ms,
        }
        if error:
            payload["error"] = error
        if meta:
            payload["meta"] = meta
        self._append("phase_end", payload)

    def record_api_call(
        self,
        api_name: str,
        url: str,
        status: Optional[int] = None,
        from_cache: bool = False,
        duration_ms: Optional[float] = None,
        error: Optional[str] = None,
    ) -> None:
        if not self.enabled:
            return
        self._api_call_count += 1
        if from_cache:
            self._api_cache_hits += 1
        if error:
            self._api_errors += 1
        self._append(
            "api_call",
            {
                "api": api_name,
                "url": url,
                "status": status,
                "from_cache": from_cache,
                "duration_ms": duration_ms,
                "error": error,
            },
        )

    def record_decision(
        self,
        decision: str,
        reason: str,
        subject: Optional[str] = None,
        details: Optional[Dict[str, Any]] = None,
    ) -> None:
        if not self.enabled:
            return
        self._append(
            "decision",
            {
                "decision": decision,
                "reason": reason,
                "subject": subject,
                "details": details or {},
            },
        )

    def record_metric(self, name: str, value: Any) -> None:
        if not self.enabled:
            return
        self.metrics[name] = value
        self._append("metric", {"name": name, "value": value})

    def _append(self, event_type: str, payload: Dict[str, Any]) -> None:
        self.records.append(
            {
                "ts": datetime.now().isoformat(),
                "type": event_type,
                **payload,
            }
        )

    # ------------------------------------------------------------------
    # Output
    # ------------------------------------------------------------------
    def write(self) -> Optional[Tuple[Path, Path]]:
        """
        Write trace.json and audit.log to output_dir.

        Returns:
            (trace_path, audit_log_path) or None if disabled.
        """
        if not self.enabled:
            return None
        self.output_dir.mkdir(parents=True, exist_ok=True)
        trace_path = self.output_dir / f"gpa_trace_{self.run_id}.json"
        audit_path = self.output_dir / f"gpa_audit_{self.run_id}.log"

        summary = self._build_summary()
        trace = {
            "run_id": self.run_id,
            "started_at": self.started_at,
            "finished_at": datetime.now().isoformat(),
            "summary": summary,
            "records": self.records,
        }
        with open(trace_path, "w", encoding="utf-8") as f:
            json.dump(trace, f, ensure_ascii=False, indent=2, default=str)

        with open(audit_path, "w", encoding="utf-8") as f:
            f.write(self._format_audit_log(summary))

        return trace_path, audit_path

    def _build_summary(self) -> Dict[str, Any]:
        phase_records = [r for r in self.records if r["type"] in ("phase_start", "phase_end")]
        phases: Dict[str, Dict[str, Any]] = {}
        for r in phase_records:
            name = r.get("phase", "unknown")
            phases.setdefault(name, {"name": name})
            if r["type"] == "phase_start":
                phases[name]["started_at"] = r["ts"]
            else:
                phases[name]["ended_at"] = r["ts"]
                phases[name]["status"] = r.get("status", "unknown")
                phases[name]["duration_ms"] = r.get("duration_ms")
        return {
            "total_records": len(self.records),
            "api_calls": self._api_call_count,
            "api_cache_hits": self._api_cache_hits,
            "api_errors": self._api_errors,
            "phases": list(phases.values()),
            "metrics": self.metrics,
        }

    def _format_audit_log(self, summary: Dict[str, Any]) -> str:
        lines = [
            f"GPA Audit Log | run_id={self.run_id}",
            f"Started: {self.started_at}",
            f"Finished: {datetime.now().isoformat()}",
            f"API calls: {summary['api_calls']} | cache hits: {summary['api_cache_hits']} | errors: {summary['api_errors']}",
            "-" * 60,
        ]
        for r in self.records:
            ts = r["ts"]
            et = r["type"]
            if et == "phase_start":
                lines.append(f"[{ts}] PHASE START {r.get('phase', '')}")
            elif et == "phase_end":
                dur = r.get("duration_ms")
                dur_str = f" ({dur}ms)" if dur is not None else ""
                err = r.get("error")
                err_str = f" ERROR={err}" if err else ""
                lines.append(f"[{ts}] PHASE END {r.get('phase', '')} status={r.get('status', '')}{dur_str}{err_str}")
            elif et == "api_call":
                cache = "CACHE" if r.get("from_cache") else "LIVE"
                status = r.get("status", "-")
                dur = r.get("duration_ms")
                dur_str = f" {dur}ms" if dur is not None else ""
                err = r.get("error")
                err_str = f" ERR={err}" if err else ""
                lines.append(f"[{ts}] API {r.get('api', '')} {cache} status={status}{dur_str} {r.get('url', '')}{err_str}")
            elif et == "decision":
                subj = r.get("subject")
                subj_str = f" [{subj}]" if subj else ""
                lines.append(f"[{ts}] DECISION{subj_str}: {r.get('decision', '')} | reason={r.get('reason', '')}")
            elif et == "metric":
                lines.append(f"[{ts}] METRIC {r.get('name', '')}={r.get('value', '')}")
        return "\n".join(lines) + "\n"


# ponytail: one self-check to ensure serialization works
def _self_check() -> None:
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        trail = AuditTrail(output_dir=Path(tmp), run_id="test")
        trail.record_phase_start("preflight", {"variants_in": 10})
        trail.record_phase_end("preflight", status="success")
        trail.record_api_call("ensembl", "https://rest.ensembl.org/lookup/symbol/hs/BRCA1", status=200, duration_ms=45.0)
        trail.record_api_call("gnomad", "https://gnomad.broadinstitute.org/api", status=429, error="rate limited")
        trail.record_decision("Tier 1", "ClinVar pathogenic + low AF", subject="BRCA1 c.1234A>G")
        trail.record_metric("tier1_count", 2)
        paths = trail.write()
        assert paths is not None
        trace_path, audit_path = paths
        assert trace_path.exists()
        assert audit_path.exists()
        with open(trace_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        assert data["summary"]["api_calls"] == 2
        assert data["summary"]["api_errors"] == 1
        print("gpa_audit_trail self-check OK")


if __name__ == "__main__":
    _self_check()
