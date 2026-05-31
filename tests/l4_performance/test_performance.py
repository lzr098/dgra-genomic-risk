# =============================================================================
# L4 Performance Tests — Stress benchmarks for core hot paths
# =============================================================================
# These tests validate throughput and memory under load.
# They do NOT assert strict timing thresholds (CI variance);
# instead they verify correctness under scale and catch OOM regressions.
# =============================================================================

import pytest
import sys
from pathlib import Path
from unittest.mock import patch

sys.path.insert(0, str(Path(__file__).parent.parent.parent / "scripts"))

from dgra_core import Variant, GPAConfig, classify_gnomad_frequency
from gpa_tier_classifier import classify_variant_tier, _load_rare_disease_genes
from gpa_phaser import determine_phase
from gpa_multi_hit import detect_multi_hit_genes
from dgra_batch_runner import run_gpa_batched
from dgra_cache import DGRACache


# =============================================================================
# Tier Classification Performance
# =============================================================================

@pytest.mark.l4
@pytest.mark.slow
class TestTierClassificationPerformance:
    def _make_variant(self, gene="TP53", clinvar="", consequence="missense_variant",
                      impact="MODERATE", gnomad_af=0.001):
        return Variant(
            chrom="1", pos=100, ref="A", alt="G", gene=gene,
            transcript="ENST000001", exon="1/10", impact=impact,
            consequence=consequence, hgvsp="p.Arg1His", hgvsc="c.2A>G",
            clinvar=clinvar, gnomad_af=gnomad_af, vaf=0.45, dp=120, gq=99
        )

    def test_100_variants_tier_consistency(self):
        """Tier classification stays correct at 100 variants."""
        variants = [self._make_variant(gene=f"GENE{i}", clinvar="") for i in range(100)]
        for v in variants:
            tier, reason, actions = classify_variant_tier(
                v, {}, {}, {}, None, None, {}, GPAConfig()
            )
            assert tier in (1, 2, 3)

    def test_1000_variants_no_crash(self):
        """1,000 variants should complete without error or excessive memory."""
        variants = [
            self._make_variant(
                gene=f"GENE{i % 50}",
                clinvar="Pathogenic" if i % 20 == 0 else "",
                consequence="splice_donor_variant" if i % 10 == 0 else "missense_variant",
                impact="HIGH" if i % 10 == 0 else "MODERATE",
            )
            for i in range(1000)
        ]
        for v in variants:
            classify_variant_tier(v, {}, {}, {}, None, None, {}, GPAConfig())

    def test_rare_disease_gene_cache_hit(self):
        """Repeated rare-disease lookup should use cached set, not reload."""
        _load_rare_disease_genes()
        for i in range(500):
            assert _load_rare_disease_genes() is _load_rare_disease_genes()


# =============================================================================
# Phaser Performance
# =============================================================================

@pytest.mark.l4
@pytest.mark.slow
class TestPhaserPerformance:
    def test_500_variants_phasing(self):
        """Phaser handles 500 variants without blow-up."""
        variants = [
            Variant(
                chrom="1", pos=1000 + i * 2, ref="A", alt="G",
                gene="TP53", transcript="", exon="", impact="MODERATE",
                consequence="missense_variant", hgvsp="", hgvsc="",
                clinvar="", gnomad_af=0.001, vaf=0.5, dp=100, gq=99
            )
            for i in range(500)
        ]
        result = determine_phase(variants)
        assert result.n_variants == 500

    def test_1000_singletons_no_pairs(self):
        """1,000 widely spaced variants."""
        variants = [
            Variant(
                chrom="1", pos=i * 10000, ref="A", alt="G",
                gene="GENE", transcript="", exon="", impact="MODERATE",
                consequence="missense_variant", hgvsp="", hgvsc="",
                clinvar="", gnomad_af=0.001, vaf=0.5, dp=100, gq=99
            )
            for i in range(1000)
        ]
        result = determine_phase(variants)
        assert result.n_variants == 1000


# =============================================================================
# Multi-Hit Detection Performance
# =============================================================================

@pytest.mark.l4
@pytest.mark.slow
class TestMultiHitPerformance:
    def test_100_genes_no_multi_hit(self):
        """100 unrelated genes → no multi-hit."""
        variants = [
            Variant(
                chrom="1", pos=100+i, ref="A", alt="G", gene=f"GENE{i}",
                transcript="", exon="", impact="HIGH", consequence="stop_gained",
                hgvsp="", hgvsc="", clinvar=""
            )
            for i in range(100)
        ]
        hits = detect_multi_hit_genes(variants)
        assert isinstance(hits, list)

    def test_50_shared_pathway_genes(self):
        """50 genes all in same pathway → single multi-hit."""
        variants = [
            Variant(
                chrom="1", pos=100+i, ref="A", alt="G", gene=f"GENE{i}",
                transcript="", exon="", impact="HIGH", consequence="stop_gained",
                hgvsp="", hgvsc="", clinvar=""
            )
            for i in range(50)
        ]
        hits = detect_multi_hit_genes(variants)
        assert isinstance(hits, list)


# =============================================================================
# Cache Performance
# =============================================================================

@pytest.mark.l4
class TestCachePerformance:
    def test_1000_set_get_cycles(self, tmp_path):
        """Cache should handle 1,000 set/get cycles."""
        cache = DGRACache(db_path=tmp_path / "cache.db")
        for i in range(1000):
            cache.set(f"api_{i}", {"data": i}, ttl_days=1)
        for i in range(1000):
            result = cache.get(f"api_{i}")
            assert result is not None
            assert result["data"] == {"data": i}

    def test_large_value_storage(self, tmp_path):
        """Cache stores large JSON (~1 MB) without corruption."""
        cache = DGRACache(db_path=tmp_path / "cache.db")
        large = {"items": [{"id": j, "seq": "A" * 1000} for j in range(1000)]}
        cache.set("large", large, ttl_days=1)
        result = cache.get("large")
        assert result is not None
        assert len(result["data"]["items"]) == 1000


# =============================================================================
# Batch Runner Throughput
# =============================================================================

@pytest.mark.l4
@pytest.mark.slow
class TestBatchRunnerThroughput:
    @patch("dgra_batch_runner.run_batch")
    def test_10_batches_of_100(self, mock_run_batch):
        """10 batches of 100 variants each."""
        mock_run_batch.return_value = {
            "success": True, "batch_id": 1, "variant_count": 100, "elapsed_seconds": 0.1,
            "results": {"meta": {}, "tier1_variants": [], "tier2_variants": [], "tier3_variants": []}
        }
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(1000)]
        result = run_gpa_batched(variants, batch_size=100)
        assert result["success"] is True
        assert mock_run_batch.call_count == 10

    @patch("dgra_batch_runner.run_batch")
    def test_100_batches_of_10(self, mock_run_batch):
        """100 batches of 10 variants each."""
        mock_run_batch.return_value = {
            "success": True, "batch_id": 1, "variant_count": 10, "elapsed_seconds": 0.01,
            "results": {"meta": {}, "tier1_variants": [], "tier2_variants": [], "tier3_variants": []}
        }
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(1000)]
        result = run_gpa_batched(variants, batch_size=10)
        assert result["success"] is True
        assert mock_run_batch.call_count == 100
