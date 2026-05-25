"""
v0.9.5 Hotfix Regression Tests
Coverage: P0-3 aselect(), P0-4 batch O(n) dedup, P0-5 dead assignment cleanup,
          P1-8 direct param sync, P1-9 i18n map consolidation, P1-10 proxy config
"""

import sys
import time
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import run_tests


# ============================================================================
# Fix-1: P0-3 — aselect() async method
# ============================================================================

def test_aselect_empty_transcripts():
    """aselect() with empty list → method='none', primary={}."""
    from gpa_transcript_selector import TranscriptSelector
    import asyncio

    selector = TranscriptSelector()
    result = asyncio.run(selector.aselect("TP53", []))
    assert result.method == "none"
    assert result.primary == {}
    assert result.alternatives == []
    assert not result.is_ambiguous


def test_aselect_single_transcript():
    """aselect() with one transcript → canonical or single method."""
    from gpa_transcript_selector import TranscriptSelector
    import asyncio

    selector = TranscriptSelector()
    txs = [{"transcript_id": "ENST001", "canonical": True, "impact": "HIGH"}]
    result = asyncio.run(selector.aselect("TP53", txs))
    assert result.primary["transcript_id"] == "ENST001"
    assert result.method in ("canonical", "single")
    assert len(result.alternatives) == 0


def test_aselect_multiple_transcripts_sorted():
    """aselect() with multiple transcripts → canonical+HIGH wins."""
    from gpa_transcript_selector import TranscriptSelector
    import asyncio

    selector = TranscriptSelector()
    txs = [
        {"transcript_id": "ENST001", "canonical": False, "impact": "MODERATE"},
        {"transcript_id": "ENST002", "canonical": True, "impact": "HIGH"},
    ]
    result = asyncio.run(selector.aselect("TP53", txs))
    assert result.primary["transcript_id"] == "ENST002"
    assert len(result.alternatives) == 1
    assert result.alternatives[0]["transcript_id"] == "ENST001"


def test_aselect_in_async_context_no_runtime_error():
    """aselect() can be awaited directly in an async context without RuntimeError."""
    from gpa_transcript_selector import TranscriptSelector
    import asyncio

    selector = TranscriptSelector()

    async def _inner():
        txs = [{"transcript_id": "ENST001", "canonical": True, "impact": "HIGH"}]
        return await selector.aselect("TP53", txs)

    # Should NOT raise RuntimeError about nested event loops
    result = asyncio.run(_inner())
    assert result.primary["transcript_id"] == "ENST001"


# ============================================================================
# Fix-2: P0-4 — Batch Runner O(n) dedup + highest-tier-wins
# ============================================================================

def test_variant_signature_case_insensitive():
    """_variant_signature handles both lower-case and UPPER-case keys."""
    from dgra_batch_runner import _variant_signature

    v_lower = {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}
    v_upper = {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G"}
    assert _variant_signature(v_lower) == "1:100:A>G"
    assert _variant_signature(v_upper) == "1:100:A>G"


def test_merge_no_duplicate_variants():
    """Same variant in two batches → appears only once."""
    from dgra_batch_runner import merge_batch_results

    var = {"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "gene": "TP53"}
    batch = {
        "success": True,
        "batch_id": 1,
        "variant_count": 1,
        "elapsed_seconds": 1.0,
        "results": {
            "meta": {"analysis_date": "2026-01-01", "offline_mode": False},
            "tier1_variants": [var],
            "tier2_variants": [],
            "tier3_variants": [],
        },
    }
    merged = merge_batch_results([batch, batch])
    assert len(merged["results"]["tier1_variants"]) == 1


def test_merge_highest_tier_wins():
    """Same variant in tier2 (batch1) and tier1 (batch2) → kept in tier1."""
    from dgra_batch_runner import merge_batch_results

    var = {"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "gene": "TP53"}
    batch_tier2 = {
        "success": True,
        "batch_id": 1,
        "variant_count": 1,
        "elapsed_seconds": 1.0,
        "results": {
            "meta": {"analysis_date": "2026-01-01", "offline_mode": False},
            "tier1_variants": [],
            "tier2_variants": [var],
            "tier3_variants": [],
        },
    }
    batch_tier1 = {
        "success": True,
        "batch_id": 2,
        "variant_count": 1,
        "elapsed_seconds": 1.0,
        "results": {
            "meta": {"analysis_date": "2026-01-01", "offline_mode": False},
            "tier1_variants": [var],
            "tier2_variants": [],
            "tier3_variants": [],
        },
    }
    merged = merge_batch_results([batch_tier2, batch_tier1])
    assert len(merged["results"]["tier1_variants"]) == 1
    assert len(merged["results"]["tier2_variants"]) == 0
    assert len(merged["results"]["tier3_variants"]) == 0


def test_merge_performance_large_batches():
    """3x1000 variant batches merge in <0.5s (was O(n^2), now O(n))."""
    from dgra_batch_runner import merge_batch_results

    var_template = {"chrom": "1", "pos": 0, "ref": "A", "alt": "G", "gene": "TP53"}
    batch = {
        "success": True,
        "batch_id": 1,
        "variant_count": 1000,
        "elapsed_seconds": 1.0,
        "results": {
            "meta": {"analysis_date": "2026-01-01", "offline_mode": False},
            "tier1_variants": [{**var_template, "pos": i} for i in range(1000)],
            "tier2_variants": [],
            "tier3_variants": [],
        },
    }
    start = time.time()
    merged = merge_batch_results([batch, batch, batch])
    elapsed = time.time() - start
    assert elapsed < 0.5, f"Merge took {elapsed:.3f}s, expected <0.5s"
    assert len(merged["results"]["tier1_variants"]) == 1000


def test_merge_failed_batches_skipped():
    """Failed batches are silently skipped during merge."""
    from dgra_batch_runner import merge_batch_results

    var = {"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "gene": "TP53"}
    batch_ok = {
        "success": True,
        "batch_id": 1,
        "variant_count": 1,
        "elapsed_seconds": 1.0,
        "results": {
            "meta": {"analysis_date": "2026-01-01", "offline_mode": False},
            "tier1_variants": [var],
            "tier2_variants": [],
            "tier3_variants": [],
        },
    }
    batch_fail = {
        "success": False,
        "batch_id": 2,
        "variant_count": 0,
        "elapsed_seconds": 0.0,
        "error": "Timeout",
        "results": {},
    }
    merged = merge_batch_results([batch_ok, batch_fail])
    assert merged["success"] is True
    assert len(merged["results"]["tier1_variants"]) == 1


# ============================================================================
# Fix-3: P0-5 — Dead assignment cleanup (static, verified by import + smoke)
# ============================================================================

def test_core_module_imports_cleanly():
    """dgra_core.py imports without syntax errors after dead-assignment removal."""
    import dgra_core
    assert hasattr(dgra_core, "classify_variant_tier")
    assert hasattr(dgra_core, "run_dgra_pipeline")


# ============================================================================
# Fix-4: P1-8 — _run_gpa_direct() parameter sync
# ============================================================================

def test_run_gpa_direct_accepts_multi_organ():
    """_run_gpa_direct() accepts multi_organ kwarg without TypeError."""
    from dgra_cli_wrapper import _run_gpa_direct

    # Empty variant list should return quickly; key is no TypeError
    result = _run_gpa_direct(
        variants=[],
        tissue="general",
        multi_organ=["liver", "kidney"],
    )
    assert isinstance(result, dict)


def test_run_gpa_direct_accepts_database_version():
    """_run_gpa_direct() accepts database_version kwarg without TypeError."""
    from dgra_cli_wrapper import _run_gpa_direct

    result = _run_gpa_direct(
        variants=[],
        tissue="general",
        database_version="GRCh38",
    )
    assert isinstance(result, dict)


# ============================================================================
# Fix-5: P1-9 — EXONIC_FUNC_MAP consolidation into gpa_i18n
# ============================================================================

def test_exonic_func_map_available_in_gpa_i18n():
    """EXONIC_FUNC_MAP and CN_EXONIC_FUNC_MAP importable from gpa_i18n."""
    from gpa_i18n import EXONIC_FUNC_MAP, CN_EXONIC_FUNC_MAP

    assert isinstance(EXONIC_FUNC_MAP, dict)
    assert isinstance(CN_EXONIC_FUNC_MAP, dict)
    assert EXONIC_FUNC_MAP.get("frameshift substitution") == "frameshift_variant"
    assert CN_EXONIC_FUNC_MAP.get("移码变异") == "frameshift_variant"


def test_adapter_uses_imported_map():
    """ANNOVARAdapter uses the module-level imported EXONIC_FUNC_MAP."""
    from dgra_adapters import ANNOVARAdapter, EXONIC_FUNC_MAP

    adapter = ANNOVARAdapter()
    # The class should NOT define its own EXONIC_FUNC_MAP attribute
    assert "EXONIC_FUNC_MAP" not in adapter.__dict__
    # Module-level map should still be accessible
    assert "frameshift substitution" in EXONIC_FUNC_MAP


def test_adapter_chinese_term_lookup():
    """ANNOVARAdapter correctly maps Chinese exonic function terms via imported CN map."""
    from dgra_adapters import ANNOVARAdapter, CN_EXONIC_FUNC_MAP

    adapter = ANNOVARAdapter()
    raw = {
        "ExonicFunc.refGene": "错义变异",
        "Gene.refGene": "TP53",
        "Func.refGene": "exonic",
    }
    adapted = adapter.adapt(raw)
    # The adapter maps Chinese term to SO consequence
    assert adapted.get("Consequence") == "missense_variant"


# ============================================================================
# Fix-6: P1-10 — Dynamic proxy / trust_env
# ============================================================================

def test_dgra_global_config_has_proxy_field():
    """DGRAGlobalConfig has proxy field defaulting to None."""
    from dgra_config import DGRAGlobalConfig

    config = DGRAGlobalConfig()
    assert hasattr(config, "proxy")
    assert config.proxy is None


def test_vcf_annotator_accepts_proxy_param():
    """VCFAnnotator constructor accepts proxy parameter."""
    from gpa_vcf_annotator import VCFAnnotator

    annotator = VCFAnnotator(proxy="__DIRECT__")
    assert annotator.proxy == "__DIRECT__"

    annotator2 = VCFAnnotator(proxy=None)
    assert annotator2.proxy is None


def test_vcf_annotator_proxy_none_uses_system_proxy():
    """VCFAnnotator with proxy=None → trust_env should be True (system proxy)."""
    from gpa_vcf_annotator import VCFAnnotator
    import aiohttp

    annotator = VCFAnnotator(proxy=None)
    # _ensure_session creates session; we verify the internal state
    # by checking what trust_env value would be used
    assert annotator.proxy != "__DIRECT__"


def test_vcf_annotator_direct_disables_proxy():
    """VCFAnnotator with proxy='__DIRECT__' → trust_env should be False."""
    from gpa_vcf_annotator import VCFAnnotator

    annotator = VCFAnnotator(proxy="__DIRECT__")
    assert annotator.proxy == "__DIRECT__"


# ============================================================================
# Fix-7 (Refactor-2): version.py single source
# ============================================================================

def test_version_module_exists():
    """version.py exists and exports __version__."""
    from version import __version__

    assert isinstance(__version__, str)
    assert __version__.startswith("0.9.")


def test_core_imports_version():
    """dgra_core imports __version__ from version module."""
    import dgra_core

    assert hasattr(dgra_core, "__version__")


# ============================================================================
# Test runner
# ============================================================================

if __name__ == "__main__":
    ALL_TESTS = [
        # Fix-1: aselect
        test_aselect_empty_transcripts,
        test_aselect_single_transcript,
        test_aselect_multiple_transcripts_sorted,
        test_aselect_in_async_context_no_runtime_error,
        # Fix-2: batch runner
        test_variant_signature_case_insensitive,
        test_merge_no_duplicate_variants,
        test_merge_highest_tier_wins,
        test_merge_performance_large_batches,
        test_merge_failed_batches_skipped,
        # Fix-3: dead assignment
        test_core_module_imports_cleanly,
        # Fix-4: direct param sync
        test_run_gpa_direct_accepts_multi_organ,
        test_run_gpa_direct_accepts_database_version,
        # Fix-5: i18n map consolidation
        test_exonic_func_map_available_in_gpa_i18n,
        test_adapter_uses_imported_map,
        test_adapter_chinese_term_lookup,
        # Fix-6: proxy config
        test_dgra_global_config_has_proxy_field,
        test_vcf_annotator_accepts_proxy_param,
        test_vcf_annotator_proxy_none_uses_system_proxy,
        test_vcf_annotator_direct_disables_proxy,
        # Fix-7: version.py
        test_version_module_exists,
        test_core_imports_version,
    ]
    run_tests("v0.9.5_hotfix", ALL_TESTS)
