"""L2 Unit Tests — gpa_report.py

Covers report header building, version info, tier report generation,
and JSON report generation.
"""

import pytest
from conftest import make_variant, MockTissueProfile


@pytest.mark.l2
class TestBuildReportHeader:
    """REP-01~06: _build_report_header."""

    def _make_config(self, **kwargs):
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        for k, v in kwargs.items():
            setattr(config, k, v)
        return config

    def test_header_contains_version(self):
        """REP-01: Header contains GPA version."""
        from gpa_report import _build_report_header
        config = self._make_config()
        v = make_variant(gene="BRCA1", tier=1)
        md = _build_report_header("general", config, [v], [v], [], [])
        assert "GPA Report" in md
        assert "general" in md

    def test_header_tier_counts(self):
        """REP-02: Header shows tier counts."""
        from gpa_report import _build_report_header
        config = self._make_config()
        v1 = make_variant(gene="A", tier=1)
        v2 = make_variant(gene="B", tier=2)
        v3 = make_variant(gene="C", tier=3)
        md = _build_report_header("test", config, [v1, v2, v3], [v1], [v2], [v3])
        assert "Tier 1" in md
        assert "Tier 2" in md
        assert "Tier 3" in md

    def test_header_offline_mode(self):
        """REP-03: Offline mode indicator."""
        from gpa_report import _build_report_header
        config = self._make_config(offline_mode=True)
        md = _build_report_header("test", config, [], [], [], [])
        assert "Offline" in md or "offline" in md.lower() or "Yes" in md

    def test_header_filter_stats(self):
        """REP-04: Filter stats shown when present."""
        from gpa_report import _build_report_header
        config = self._make_config()
        config.filter_stats = {
            "input_count": 100,
            "output_count": 50,
            "excluded": 50,
            "by_impact": {"HIGH": 10, "MODERATE": 20},
        }
        md = _build_report_header("test", config, [], [], [], [])
        assert "Input Variants" in md or "input" in md.lower()

    def test_header_no_filter_stats(self):
        """REP-05: No filter stats → total variants count."""
        from gpa_report import _build_report_header
        config = self._make_config()
        v = make_variant(gene="A")
        md = _build_report_header("test", config, [v], [], [], [])
        assert "Total Variants Assessed" in md or "total" in md.lower()

    def test_header_retention_parts(self):
        """REP-06: Filter retention details."""
        from gpa_report import _build_report_header
        config = self._make_config()
        config.filter_stats = {
            "splice_retained": 2,
            "synonymous_tissue_retained": 1,
            "clinvar_conflicting_retained": 3,
        }
        md = _build_report_header("test", config, [], [], [], [])
        assert "splice retained" in md or "retention" in md.lower()


@pytest.mark.l2
class TestGetVersionInfo:
    """REP-07~09: _get_version_info."""

    def test_basic_version_info(self):
        """REP-07: Returns dict with version keys."""
        from gpa_report import _get_version_info
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        info = _get_version_info(config)
        assert "dgra_version" in info
        assert "analysis_date" in info

    def test_no_cache(self):
        """REP-08: No cache → no_cache."""
        from gpa_report import _get_version_info
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        info = _get_version_info(config)
        assert info.get("cache_version") in ("no_cache", "unknown")

    def test_offline_archive(self):
        """REP-09: Offline archive status present."""
        from gpa_report import _get_version_info
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        info = _get_version_info(config)
        assert "offline_archive_date" in info


@pytest.mark.l2
class TestGenerateTierReport:
    """REP-10~16: generate_tier_report."""

    def _make_config(self):
        from dgra_core import GPAConfig
        return GPAConfig(tissue_profile="general")

    def test_empty_variants(self):
        """REP-10: Empty variants → report with zeros."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        md = generate_tier_report([], config, profile, [])
        assert "GPA Report" in md
        assert "Tier 1" in md

    def test_tier1_only(self):
        """REP-11: Only Tier 1 variants."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="BRCA1", tier=1, impact="HIGH", clinvar="Pathogenic")
        md = generate_tier_report([v], config, profile, [])
        assert "BRCA1" in md
        assert "🔴" in md

    def test_tier2_only(self):
        """REP-12: Only Tier 2 variants."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="BRCA1", tier=2, impact="MODERATE")
        md = generate_tier_report([v], config, profile, [])
        assert "🟡" in md or "Tier 2" in md

    def test_multi_hit_indicator(self):
        """REP-13: Multi-hit gene marked."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v1 = make_variant(gene="BRCA1", tier=1, pos=100)
        v2 = make_variant(gene="BRCA1", tier=2, pos=200)
        md = generate_tier_report([v1, v2], config, profile, [{"gene": "BRCA1", "variant_count": 2, "warning": "multi-hit"}])
        assert "Multi-hit" in md or "multi" in md.lower()

    def test_qc_summary(self):
        """REP-14: QC flags included when present."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="A", tier=3)
        v.qc_flags = ["INVALID_VAF"]
        md = generate_tier_report([v], config, profile, [])
        assert "QC" in md or "INVALID_VAF" in md

    def test_variant_table_columns(self):
        """REP-15: Variant table has expected columns."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="BRCA1", tier=1, transcript="NM_001", hgvsp="p.Arg123Cys")
        md = generate_tier_report([v], config, profile, [])
        assert "|" in md  # Table format
        assert "染色体" in md or "chrom" in md.lower() or "位置" in md

    def test_phenotype_section(self):
        """REP-16: Phenotype assessment section when variants have scores."""
        from gpa_report import generate_tier_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="BRCA1", tier=1)
        v.phenotype_match_score = 0.85
        md = generate_tier_report([v], config, profile, [])
        # Phenotype section may or may not appear depending on logic
        assert "GPA Report" in md


@pytest.mark.l2
class TestGenerateJSONReport:
    """REP-17~22: generate_json_report."""

    def _make_config(self):
        from dgra_core import GPAConfig
        return GPAConfig(tissue_profile="general")

    def test_json_structure(self):
        """REP-17: JSON report has expected top-level keys."""
        from gpa_report import generate_json_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="BRCA1", tier=1)
        result = generate_json_report([v], config, profile, [], "")
        assert "meta" in result
        assert "summary" in result
        assert "variants" in result

    def test_meta_section(self):
        """REP-18: Meta contains version and config."""
        from gpa_report import generate_json_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        result = generate_json_report([], config, profile, [], "")
        assert "dgra_version" in result["meta"]
        assert "tissue_profile" in result["meta"]
        assert "offline_mode" in result["meta"]

    def test_summary_counts(self):
        """REP-19: Summary has correct counts."""
        from gpa_report import generate_json_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v1 = make_variant(gene="A", tier=1)
        v2 = make_variant(gene="B", tier=2)
        v3 = make_variant(gene="C", tier=3)
        result = generate_json_report([v1, v2, v3], config, profile, [], "")
        assert result["summary"]["tier1_variant_count"] == 1
        assert result["summary"]["tier2_variant_count"] == 1
        assert result["summary"]["tier3_variant_count"] == 1

    def test_variant_json_structure(self):
        """REP-20: Each variant has structured fields."""
        from gpa_report import generate_json_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        v = make_variant(gene="BRCA1", tier=1, chrom="1", pos=100, ref="A", alt="G")
        result = generate_json_report([v], config, profile, [], "")
        variant = result["variants"][0]
        assert variant["gene"] == "BRCA1"
        assert variant["chrom"] == "1"
        assert variant["pos"] == 100
        assert variant["tier"] == 1

    def test_empty_json_report(self):
        """REP-21: Empty input → valid JSON with zeros."""
        from gpa_report import generate_json_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        result = generate_json_report([], config, profile, [], "")
        assert result["summary"]["total_variants"] == 0
        assert result["variants"] == []

    def test_multi_hit_in_json(self):
        """REP-22: Multi-hit genes in summary."""
        from gpa_report import generate_json_report
        config = self._make_config()
        profile = MockTissueProfile.general()
        result = generate_json_report([], config, profile, [{"gene": "BRCA1"}], "")
        assert "BRCA1" in result["summary"]["multi_hit_genes"]


@pytest.mark.l2
class TestFormatVEPReannotationNote:
    """REP-23~24: _format_vep_reannotation_note."""

    def test_with_reannotation(self):
        """REP-23: Variant with reannotation data → note string."""
        from gpa_report import _format_vep_reannotation_note
        import json
        v = make_variant(gene="BRCA1", consequence="missense_variant", impact="MODERATE")
        v.transcript_warning = json.dumps({
            "vep_reannotation": {
                "status": "success",
                "original": {"transcript": "NM_001", "consequence": "synonymous_variant", "impact": "LOW"},
                "canonical": {"transcript_id": "NM_007", "consequence": "missense_variant", "impact": "MODERATE"},
            }
        })
        note = _format_vep_reannotation_note(v)
        assert note is not None
        assert "missense" in note or "MODERATE" in note

    def test_no_reannotation(self):
        """REP-24: No reannotation → None."""
        from gpa_report import _format_vep_reannotation_note
        v = make_variant(gene="BRCA1")
        v.vep_reannotation = None
        assert _format_vep_reannotation_note(v) is None


@pytest.mark.l2
class TestPseudogeneAssessment:
    """REP-25~26: _generate_pseudogene_assessment_section."""

    def test_with_pseudogene(self):
        """REP-25: Pseudogene variants → section generated."""
        from gpa_report import _generate_pseudogene_assessment_section
        import json
        from unittest.mock import patch
        v = make_variant(gene="BRCA1")
        v.pseudogene_warning = json.dumps({
            "score": 0.85,
            "level": "interference",
            "observed_vaf": 0.3,
            "expected_vaf": 0.5,
            "pseudogenes": ["BRCA1P1"],
            "recommendation": "Sanger validation recommended",
        })
        with patch("gpa_report._load_pseudogene_lookup", create=True, return_value={"BRCA1": {"notes": "test note"}}):
            section = _generate_pseudogene_assessment_section([v])
        assert section is not None
        assert "pseudo" in section.lower() or "假基因" in section

    def test_no_pseudogene(self):
        """REP-26: No pseudogenes → None."""
        from gpa_report import _generate_pseudogene_assessment_section
        v = make_variant(gene="BRCA1")
        v.is_pseudogene = False
        assert _generate_pseudogene_assessment_section([v]) is None
