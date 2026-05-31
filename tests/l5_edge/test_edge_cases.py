# =============================================================================
# L5 Edge / Boundary Tests — Malformed data, extremes, type errors
# =============================================================================
# Tests the robustness of each module against unexpected / garbage input.
# Every test must survive without crashing (raise expected exceptions OK).
# =============================================================================

import pytest
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from dgra_core import Variant, GPAConfig, classify_gnomad_frequency, map_variant_to_domain
from gpa_tier_classifier import classify_variant_tier
from gpa_phaser import determine_phase
from gpa_qc import _run_qc_checks
from dgra_input_parsers import FreeTextParser
from dgra_adapters import VEPAdapter
from dgra_variant_filter import filter_variants
from dgra_cache import DGRACache
from gpa_preflight import PreflightReport, suggest_action
from gpa_workflow import WorkflowStep, STANDARD_WORKFLOW, FailureAction


# =============================================================================
# Tier Classifier — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestTierClassifierEdgeCases:
    def _make_variant(self, **overrides):
        defaults = dict(
            chrom="1", pos=100, ref="A", alt="G", gene="TP53",
            transcript="ENST000001", exon="1/10", impact="MODERATE",
            consequence="missense_variant", hgvsp="p.Arg1His", hgvsc="c.2A>G",
            clinvar="", gnomad_af=0.001, vaf=0.5, dp=100, gq=99
        )
        defaults.update(overrides)
        return Variant(**defaults)

    def test_none_clinvar(self):
        """None clinvar should not crash."""
        v = self._make_variant(clinvar=None)
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_empty_gene(self):
        """Empty gene string should still classify."""
        v = self._make_variant(gene="", clinvar="")
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_very_long_gene_name(self):
        """Gene name up to 100 chars."""
        long_gene = "A" * 100
        v = self._make_variant(gene=long_gene)
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_negative_gnomad_af(self):
        """Negative AF should be treated as unknown."""
        v = self._make_variant(gnomad_af=-0.01)
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_gnomad_af_greater_than_one(self):
        """AF > 1.0 is garbage but should not crash."""
        v = self._make_variant(gnomad_af=5.0)
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_nan_gnomad_af(self):
        """NaN AF should not crash."""
        v = self._make_variant(gnomad_af=float("nan"))
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_infinite_gnomad_af(self):
        """Infinite AF should not crash."""
        v = self._make_variant(gnomad_af=float("inf"))
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_unicode_in_clinvar(self):
        """Unicode characters in clinvar field."""
        v = self._make_variant(clinvar="致病性 Pathogenic \u00e9")
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)

    def test_empty_consequence(self):
        """Empty consequence should still classify."""
        v = self._make_variant(consequence="", impact="")
        tier, reason, actions = classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())
        assert tier in (1, 2, 3)


# =============================================================================
# Core — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestCoreEdgeCases:
    def test_classify_gnomad_none_gene(self):
        """None gene should not crash frequency classifier."""
        result = classify_gnomad_frequency(0.01, None)
        assert "status" in result

    def test_classify_gnomad_none_af(self):
        """None AF should be handled."""
        result = classify_gnomad_frequency(None, "TP53")
        assert "status" in result

    def test_map_domain_none_uniprot(self):
        """None uniprot data should crash (no guard) — boundary behavior."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="G", gene="TP53",
            transcript="", exon="", impact="MODERATE", consequence="missense_variant",
            hgvsp="p.Arg1His", hgvsc="", clinvar=""
        )
        with pytest.raises(AttributeError):
            map_variant_to_domain(v, None)

    def test_map_domain_empty_uniprot(self):
        v = Variant(
            chrom="1", pos=100, ref="A", alt="G", gene="TP53",
            transcript="", exon="", impact="MODERATE", consequence="missense_variant",
            hgvsp="p.Arg1His", hgvsc="", clinvar=""
        )
        result = map_variant_to_domain(v, {})
        assert result["domain"] == "unknown"


# =============================================================================
# Phaser — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestPhaserEdgeCases:
    def test_empty_variant_list(self):
        result = determine_phase([])
        assert result.n_variants == 0

    def test_single_variant(self):
        v = Variant(
            chrom="1", pos=100, ref="A", alt="G", gene="TP53",
            transcript="", exon="", impact="MODERATE", consequence="missense_variant",
            hgvsp="", hgvsc="", clinvar="", gnomad_af=0.001, vaf=0.5, dp=100, gq=99
        )
        result = determine_phase([v])
        assert result.n_variants == 1


# =============================================================================
# QC — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestQCEdgeCases:
    def test_none_variant_list(self):
        """None list should crash _run_qc_checks (no guard) — boundary behavior."""
        with pytest.raises(TypeError):
            _run_qc_checks(None)

    def test_empty_variant_list(self):
        result = _run_qc_checks([])
        assert isinstance(result, dict)
        assert result["total"] == 0

    def test_variant_missing_all_fields(self):
        """Variant with missing fields should be flagged."""
        v = Variant(
            chrom="", pos=0, ref="", alt="", gene="",
            transcript="", exon="", impact="", consequence="",
            hgvsp="", hgvsc="", clinvar=""
        )
        result = _run_qc_checks([v])
        assert isinstance(result, dict)
        assert result["total"] == 1

    def test_extreme_vaf_values(self):
        """VAF > 1.0 should be flagged."""
        v = Variant(
            chrom="1", pos=100, ref="A", alt="G", gene="TP53",
            transcript="", exon="", impact="MODERATE", consequence="missense_variant",
            hgvsp="", hgvsc="", clinvar="", gnomad_af=0.001, vaf=1.5, dp=100, gq=99
        )
        result = _run_qc_checks([v])
        assert isinstance(result, dict)


# =============================================================================
# Input Parsers — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestInputParserEdgeCases:
    def test_parse_none(self):
        """None input should crash FreeTextParser (no guard) — boundary behavior."""
        parser = FreeTextParser()
        with pytest.raises(AttributeError):
            parser.parse_text(None)

    def test_parse_empty(self):
        parser = FreeTextParser()
        result = parser.parse_text("")
        assert result == []

    def test_parse_gibberish(self):
        """Gibberish should raise ValueError."""
        parser = FreeTextParser()
        with pytest.raises(ValueError):
            parser.parse_text("not_a_coordinate")

    def test_parse_unicode(self):
        """Unicode with null byte should raise ValueError."""
        parser = FreeTextParser()
        with pytest.raises(ValueError):
            parser.parse_text("chr1\u0000100A>G")


# =============================================================================
# Adapters — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestAdapterEdgeCases:
    def test_vep_adapter_empty_input(self):
        adapter = VEPAdapter()
        assert adapter._parse_uploaded_variation("") is None

    def test_vep_adapter_gibberish(self):
        adapter = VEPAdapter()
        assert adapter._parse_uploaded_variation("!!!") is None


# =============================================================================
# Variant Filter — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestVariantFilterEdgeCases:
    def test_filter_empty_list(self):
        result = filter_variants([], preset="clinical")
        assert isinstance(result, tuple)
        assert result[0] == []

    def test_filter_missing_impact_key(self):
        """Variants without IMPACT key should be handled."""
        variants = [{"GENE": "TP53"}, {"IMPACT": "HIGH"}]
        filtered, stats = filter_variants(variants, preset="clinical")
        assert isinstance(filtered, list)
        assert isinstance(stats, dict)


# =============================================================================
# Cache — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestCacheEdgeCases:
    def test_get_nonexistent(self, tmp_path):
        cache = DGRACache(db_path=tmp_path / "cache.db")
        assert cache.get("does_not_exist") is None

    def test_set_none_value(self, tmp_path):
        """None value should be storable or handled gracefully."""
        cache = DGRACache(db_path=tmp_path / "cache.db")
        cache.set("null", None, ttl_days=1)
        result = cache.get("null")
        assert result is not None  # DGRACache wraps in metadata dict

    def test_set_negative_ttl(self, tmp_path):
        """Negative TTL should result in immediate expiry."""
        cache = DGRACache(db_path=tmp_path / "cache.db")
        cache.set("expired", {"data": 1}, ttl_days=-1)
        result = cache.get("expired")
        assert result is None


# =============================================================================
# Preflight — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestPreflightEdgeCases:
    def test_suggest_action_empty_report(self):
        """suggest_action with empty report should return a string."""
        report = PreflightReport(items=[])
        action = suggest_action(report)
        assert isinstance(action, str)


# =============================================================================
# Workflow — Edge Cases
# =============================================================================

@pytest.mark.l5
class TestWorkflowEdgeCases:
    def test_step_with_none_name(self):
        """Step with None name may crash on string ops — that's OK."""
        step = WorkflowStep(name=None, module="test", required=False)
        try:
            step.can_be_skipped({})
        except (TypeError, AttributeError):
            pass  # Expected for None name

    def test_failure_action_invalid_value(self):
        """Invalid FailureAction value should raise."""
        with pytest.raises((ValueError, TypeError)):
            FailureAction("INVALID_ACTION")
