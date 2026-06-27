#!/usr/bin/env python3
"""L2 unit tests for gpa_audit_trail, WorkflowEngine + APIHub audit wiring."""
import json
import pytest
import asyncio
from pathlib import Path
from unittest.mock import MagicMock

import gpa_audit_trail as at
from gpa_workflow_engine import WorkflowEngine, WorkflowError


class TestAuditTrail:
    def test_disabled_does_nothing(self, tmp_path):
        trail = at.AuditTrail(output_dir=tmp_path, enabled=False)
        trail.record_phase_start("p")
        trail.record_api_call("x", "http://x")
        assert trail.write() is None
        assert len(list(tmp_path.iterdir())) == 0

    def test_writes_trace_and_log(self, tmp_path):
        trail = at.AuditTrail(output_dir=tmp_path, run_id="r1")
        trail.record_phase_start("preflight", {"variants_in": 3})
        trail.record_phase_end("preflight", status="success")
        trail.record_api_call("ensembl", "http://e/1", status=200, duration_ms=12.0)
        trail.record_api_call("gnomad", "http://g/1", status=429, error="rate limited")
        trail.record_decision("Tier 1", "ClinVar pathogenic", subject="BRCA1 c.1A>G")
        trail.record_metric("tier1_count", 1)
        paths = trail.write()
        assert paths is not None
        trace_path, audit_path = paths
        assert trace_path.exists()
        assert audit_path.exists()

        with open(trace_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        assert data["run_id"] == "r1"
        assert data["summary"]["api_calls"] == 2
        assert data["summary"]["api_cache_hits"] == 0
        assert data["summary"]["api_errors"] == 1
        assert data["summary"]["metrics"]["tier1_count"] == 1
        assert len(data["summary"]["phases"]) == 1

        log_text = audit_path.read_text(encoding="utf-8")
        assert "PHASE START preflight" in log_text
        assert "API ensembl LIVE status=200" in log_text
        assert "API gnomad LIVE status=429" in log_text
        assert "DECISION [BRCA1 c.1A>G]: Tier 1" in log_text

    def test_cache_hit_counts(self, tmp_path):
        trail = at.AuditTrail(output_dir=tmp_path)
        trail.record_api_call("x", "u", status=200, from_cache=True)
        trail.record_api_call("y", "u", status=200, from_cache=True)
        trail.record_api_call("z", "u", status=500, error="boom")
        summary = trail._build_summary()
        assert summary["api_calls"] == 3
        assert summary["api_cache_hits"] == 2
        assert summary["api_errors"] == 1


class TestWorkflowEngineAudit:
    @pytest.mark.asyncio
    async def test_records_phases_and_writes_output(self, tmp_path):
        trail = at.AuditTrail(output_dir=tmp_path, run_id="wf1")

        async def add_one(ctx):
            ctx["x"] = ctx.get("x", 0) + 1
            return ctx

        engine = WorkflowEngine(
            phases=[("p1", add_one)],
            checkpoint_dir=tmp_path / "ckpt",
            audit_trail=trail,
        )
        ctx = await engine.run({"variants_data": [1, 2, 3]})
        assert ctx["x"] == 1
        assert "trace_path" in ctx.get("_workflow_meta", {})
        assert "audit_log_path" in ctx.get("_workflow_meta", {})

        trace_path = Path(ctx["_workflow_meta"]["trace_path"])
        with open(trace_path, "r", encoding="utf-8") as f:
            data = json.load(f)
        assert data["summary"]["metrics"]["input_variants"] == 3
        phase_names = [p["name"] for p in data["summary"]["phases"]]
        assert "p1" in phase_names

    @pytest.mark.asyncio
    async def test_failure_recorded(self, tmp_path):
        trail = at.AuditTrail(output_dir=tmp_path, run_id="wf_fail")

        async def boom(ctx):
            raise ValueError("nope")

        engine = WorkflowEngine(
            phases=[("p1", boom)],
            checkpoint_dir=tmp_path / "ckpt",
            audit_trail=trail,
        )
        with pytest.raises(WorkflowError):
            await engine.run({})

        # Even on failure the trail may not have been written by the engine,
        # but the failed phase should be recorded in memory.
        phases = [r for r in trail.records if r["type"] == "phase_end"]
        assert any(r.get("status") == "failed" and "nope" in (r.get("error") or "") for r in phases)


class TestAPIHubAudit:
    @pytest.mark.asyncio
    async def test_records_live_and_cached_calls(self, tmp_path):
        from api_hub import APIHub
        from dgra_config import DGRAGlobalConfig

        trail = at.AuditTrail(output_dir=tmp_path, run_id="hub1")
        cache = MagicMock()
        cache.get = MagicMock(return_value={"data": {"ok": True}, "http_status": 200, "confidence": "high"})
        cache.set = MagicMock()

        config = DGRAGlobalConfig()
        config.offline_mode = False
        hub = APIHub(config, cache=cache, detect_proxy=False, audit_trail=trail)
        await hub.setup()
        try:
            result = await hub.request("ensembl", "/lookup/symbol/human/BRCA1")
        finally:
            await hub.close()

        assert result["from_cache"] is True
        assert result["data"]["ok"] is True
        api_records = [r for r in trail.records if r["type"] == "api_call"]
        assert len(api_records) == 1
        assert api_records[0]["from_cache"] is True
        assert api_records[0]["status"] == 200

    @pytest.mark.asyncio
    async def test_records_offline_call(self, tmp_path):
        from api_hub import APIHub
        from dgra_config import DGRAGlobalConfig

        trail = at.AuditTrail(output_dir=tmp_path, run_id="hub_off")
        config = DGRAGlobalConfig()
        config.offline_mode = True
        hub = APIHub(config, cache=None, detect_proxy=False, audit_trail=trail)
        await hub.setup()
        try:
            result = await hub.request("ensembl", "/lookup/symbol/human/BRCA1")
        finally:
            await hub.close()

        assert "error" in result
        api_records = [r for r in trail.records if r["type"] == "api_call"]
        assert len(api_records) == 1
        assert api_records[0]["error"] == "offline mode"
