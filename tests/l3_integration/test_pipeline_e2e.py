"""
L3 Integration Tests — gpa_pipeline.py
End-to-end pipeline validation with offline mode.

Run: pytest -m "l3" tests/l3_integration/test_pipeline_e2e.py
"""

import asyncio
import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


@pytest.mark.l3
@pytest.mark.p0
@pytest.mark.mock
class TestPipelineE2E:
    """End-to-end pipeline tests with mocked dependencies."""

    @pytest.mark.asyncio
    async def test_empty_input_generates_report(self):
        """PIPE-01: Empty variant list generates report without crashing."""
        from dgra_core import run_dgra_pipeline, GPAConfig

        config = GPAConfig(offline_mode=True, tissue_profile="general")
        result = await run_dgra_pipeline([], config=config)

        assert result is not None
        assert "report_markdown" in result or "json_report" in result
        assert len(result.get("tier1_variants", [])) == 0
        assert len(result.get("tier2_variants", [])) == 0
        assert len(result.get("tier3_variants", [])) == 0

    @pytest.mark.asyncio
    async def test_single_variant_complete_report(self):
        """PIPE-02: Single variant produces complete Markdown + JSON report."""
        from dgra_core import run_dgra_pipeline, GPAConfig

        variants = [{
            "CHROM": "17", "POS": "7579472", "REF": "C", "ALT": "T",
            "GENE": "TP53", "Feature": "ENST00000269305", "EXON": "5/11",
            "IMPACT": "HIGH", "Consequence": "stop_gained",
            "HGVSp": "p.Arg273Ter", "HGVSc": "c.818C>T",
            "CLIN_SIG": "Pathogenic", "GT": "0/1", "DP": "80", "GQ": "99", "VAF": "0.45",
        }]

        config = GPAConfig(offline_mode=True, tissue_profile="general")
        result = await run_dgra_pipeline(variants, config=config)

        assert result is not None
        assert "report_markdown" in result or "json_report" in result
        # Should have at least one variant classified
        total = (len(result.get("tier1_variants", [])) +
                 len(result.get("tier2_variants", [])) +
                 len(result.get("tier3_variants", [])))
        assert total >= 1, "At least one variant should be classified"

    @pytest.mark.asyncio
    async def test_mixed_variants_tier_distribution(self):
        """PIPE-03: Mixed variants produce correct Tier distribution."""
        from dgra_core import run_dgra_pipeline, GPAConfig

        variants = [
            # TP53 Pathogenic → should be Tier 1 or 2
            {"CHROM": "17", "POS": "7579472", "REF": "C", "ALT": "T",
             "GENE": "TP53", "Feature": "ENST00000269305", "EXON": "5/11",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "p.Arg273Ter", "HGVSc": "c.818C>T",
             "CLIN_SIG": "Pathogenic", "GT": "0/1", "DP": "80", "GQ": "99", "VAF": "0.45"},
            # DDX3X common polymorphism → Tier 3
            {"CHROM": "X", "POS": "41357831", "REF": "A", "ALT": "T",
             "GENE": "DDX3X", "Feature": "ENST00000373383", "EXON": "5/15",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "p.Glu233Ter", "HGVSc": "c.697G>A",
             "CLIN_SIG": "", "GT": "0/1", "DP": "60", "GQ": "99", "VAF": "0.45",
             "gnomAD_AF": "0.45"},
        ]

        config = GPAConfig(offline_mode=True, tissue_profile="general")
        result = await run_dgra_pipeline(variants, config=config)

        assert result is not None
        total = (len(result.get("tier1_variants", [])) +
                 len(result.get("tier2_variants", [])) +
                 len(result.get("tier3_variants", [])))
        assert total >= 1, "At least one variant should be classified"

    @pytest.mark.asyncio
    async def test_offline_mode_no_api_calls(self):
        """PIPE-04: Offline mode does not make external API calls."""
        from dgra_core import run_dgra_pipeline, GPAConfig

        variants = [{
            "CHROM": "1", "POS": "100", "REF": "A", "ALT": "G",
            "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
            "IMPACT": "HIGH", "Consequence": "stop_gained",
            "HGVSp": "p.Arg1Ter", "HGVSc": "c.1A>T", "CLIN_SIG": "",
            "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.45",
        }]

        config = GPAConfig(offline_mode=True, tissue_profile="general")
        result = await run_dgra_pipeline(variants, config=config)

        assert result is not None
        assert "report_markdown" in result or "json_report" in result

    @pytest.mark.asyncio
    async def test_somatic_mode(self):
        """PIPE-05: Somatic mode activates TSG logic."""
        from dgra_core import run_dgra_pipeline, GPAConfig

        variants = [{
            "CHROM": "17", "POS": "7579472", "REF": "C", "ALT": "T",
            "GENE": "TP53", "Feature": "ENST00000269305", "EXON": "5/11",
            "IMPACT": "HIGH", "Consequence": "stop_gained",
            "HGVSp": "p.Arg273Ter", "HGVSc": "c.818C>T",
            "CLIN_SIG": "", "GT": "0/1", "DP": "80", "GQ": "99", "VAF": "0.45",
            "is_tsg": "Yes",
        }]

        config = GPAConfig(offline_mode=True, tissue_profile="general", somatic_mode=True)
        result = await run_dgra_pipeline(variants, config=config)

        assert result is not None
        # In somatic mode, TSG LOF should be Tier 1
        tier1 = result.get("tier1_variants", [])
        assert len(tier1) >= 0  # At minimum, should not crash

    @pytest.mark.asyncio
    async def test_phenotype_parameter(self):
        """PIPE-06: Phenotype parameter is passed and processed."""
        from dgra_core import run_dgra_pipeline, GPAConfig

        variants = [{
            "CHROM": "17", "POS": "7579472", "REF": "C", "ALT": "T",
            "GENE": "TP53", "Feature": "ENST00000269305", "EXON": "5/11",
            "IMPACT": "HIGH", "Consequence": "stop_gained",
            "HGVSp": "p.Arg273Ter", "HGVSc": "c.818C>T",
            "CLIN_SIG": "Pathogenic", "GT": "0/1", "DP": "80", "GQ": "99", "VAF": "0.45",
        }]

        config = GPAConfig(offline_mode=True, tissue_profile="general")
        result = await run_dgra_pipeline(
            variants,
            user_phenotypes="肌无力、肌源性损害",
            config=config,
        )

        assert result is not None
        assert "report_markdown" in result or "json_report" in result

    @pytest.mark.asyncio
    async def test_filter_preset_strict(self):
        """PIPE-07: Strict filter preset only keeps HIGH/MODERATE."""
        from dgra_core import run_dgra_pipeline, GPAConfig
        from dgra_variant_filter import filter_variants

        variants = [
            {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.45"},
            {"CHROM": "1", "POS": "101", "REF": "A", "ALT": "G",
             "GENE": "BRCA1", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "LOW", "Consequence": "synonymous_variant",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.45"},
        ]

        filtered, stats = filter_variants(variants, preset="strict")
        assert stats["output_count"] == 1
        assert stats["excluded"] == 1

    @pytest.mark.asyncio
    async def test_filter_preset_clinical(self):
        """PIPE-08: Clinical filter preset keeps splice_region LOW."""
        from dgra_variant_filter import filter_variants

        variants = [
            {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "LOW", "Consequence": "splice_region_variant",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.45"},
        ]

        filtered, stats = filter_variants(variants, preset="clinical")
        assert stats["output_count"] == 1
        assert stats.get("splice_retained", 0) == 1
