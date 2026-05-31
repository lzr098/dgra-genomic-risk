"""
L2 Unit Tests — gpa_multi_hit.py
Compound heterozygosity detection and pairwise phase analysis.

Run: pytest -m "l2 and multihit" tests/l2_unit/test_multi_hit.py
"""

import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import make_variant


@pytest.mark.l2
@pytest.mark.multihit
@pytest.mark.p0
class TestDetectMultiHit:
    """Test detect_multi_hit_genes: compound heterozygosity detection."""

    def test_single_variant_no_multihit(self):
        """MH-01: Single variant in gene → no multi-hit."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        result = detect_multi_hit_genes([v1])
        assert result == []

    def test_two_pathogenic_same_gene(self):
        """MH-02: Two HIGH-impact variants in same gene → multi-hit detected."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=500)
        result = detect_multi_hit_genes([v1, v2])
        assert len(result) == 1
        assert result[0]["gene"] == "CFTR"
        assert result[0]["pathogenic_count"] == 2
        assert result[0]["warning"] == "MULTI_HIT_GENE"

    def test_two_variants_different_genes(self):
        """MH-03: Two variants in different genes → no multi-hit."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="BRCA1", gt="0/1", impact="HIGH", pos=500)
        result = detect_multi_hit_genes([v1, v2])
        assert result == []

    def test_benign_excluded(self):
        """MH-04: Benign variants are excluded from pathogenic count."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(
            gene="CFTR", gt="0/1", impact="LOW",
            clinvar="Benign", consequence="synonymous_variant", pos=500
        )
        result = detect_multi_hit_genes([v1, v2])
        # Only v1 has pathogenic evidence → no multi-hit
        assert result == []

    def test_three_pathogenic_same_gene(self):
        """MH-05: Three pathogenic variants → multi-hit with pairwise analysis."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=200)
        v3 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=800)
        result = detect_multi_hit_genes([v1, v2, v3])
        assert len(result) == 1
        assert result[0]["pathogenic_count"] == 3
        # Pairwise analysis: (v1,v2)=100bp, (v1,v3)=700bp, (v2,v3)=600bp
        pairwise = result[0].get("pairwise_phase_analysis", [])
        assert len(pairwise) == 3  # All pairs <= 1000bp

    def test_pairwise_excludes_distant(self):
        """MH-06: Pairwise analysis excludes pairs > 1000bp."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=200)
        v3 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=1500)
        result = detect_multi_hit_genes([v1, v2, v3])
        assert len(result) == 1
        pairwise = result[0].get("pairwise_phase_analysis", [])
        # Only (v1,v2)=100bp <= 1000; (v1,v3)=1400 and (v2,v3)=1300 excluded
        assert len(pairwise) == 1
        assert pairwise[0]["distance_bp"] == 100

    def test_phase_result_structure(self):
        """MH-07: Phase result has all required fields."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="1|0", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0|1", impact="HIGH", pos=200)
        result = detect_multi_hit_genes([v1, v2])
        assert len(result) == 1
        phase = result[0]["phase_result"]
        assert "status" in phase
        assert "confidence" in phase
        assert "method" in phase
        assert "evidence" in phase
        assert "max_gap_bp" in phase
        assert "min_gap_bp" in phase
        assert "n_variants" in phase

    def test_trans_elevates_risk(self):
        """MH-08: Trans phase indicates compound heterozygosity risk."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="1|0", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0|1", impact="HIGH", pos=200)
        result = detect_multi_hit_genes([v1, v2])
        assert result[0]["phase_result"]["status"] == "trans"
        assert "compound" in result[0]["phases"]["trans"].lower()

    def test_clinical_significance_populated(self):
        """MH-09: phase_clinical_significance is populated."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="1|0", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="1|0", impact="HIGH", pos=200)
        result = detect_multi_hit_genes([v1, v2])
        assert result[0]["phase_clinical_significance"] != "未知"
        assert "单倍型" in result[0]["phase_clinical_significance"]

    def test_action_and_required_evidence(self):
        """MH-10: Action and required_evidence fields are present."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=200)
        result = detect_multi_hit_genes([v1, v2])
        assert "Priority P1" in result[0]["action"]
        assert len(result[0]["required_evidence"]) >= 2
        assert any("trio" in e.lower() for e in result[0]["required_evidence"])

    def test_variant_details_populated(self):
        """MH-11: pathogenic_variants details include hgvsp, impact, clinvar."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(
            gene="CFTR", gt="0/1", impact="HIGH", pos=100,
            hgvsp="p.Phe508del", hgvsc="c.1521_1523del"
        )
        v2 = make_variant(
            gene="CFTR", gt="0/1", impact="HIGH", pos=500,
            hgvsp="p.Gly551Asp", hgvsc="c.1652G>A"
        )
        result = detect_multi_hit_genes([v1, v2])
        details = result[0]["pathogenic_variants"]
        assert len(details) == 2
        assert details[0]["hgvsp"] == "p.Phe508del"
        assert details[0]["impact"] == "HIGH"

    def test_total_vs_pathogenic_count(self):
        """MH-12: variant_count = total, pathogenic_count = filtered."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=200)
        v3 = make_variant(
            gene="CFTR", gt="0/1", impact="LOW",
            clinvar="Benign", consequence="synonymous_variant", pos=500
        )
        result = detect_multi_hit_genes([v1, v2, v3])
        assert result[0]["variant_count"] == 3  # total
        assert result[0]["pathogenic_count"] == 2  # only HIGH-impact ones

    def test_empty_input(self):
        """MH-13: Empty variant list → empty result."""
        from gpa_multi_hit import detect_multi_hit_genes
        result = detect_multi_hit_genes([])
        assert result == []

    def test_gtex_data_passed_through(self):
        """MH-14: gtex_data parameter is accepted (does not crash)."""
        from gpa_multi_hit import detect_multi_hit_genes
        v1 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", impact="HIGH", pos=200)
        gtex = {"CFTR": {"median_tpm": 5.0}}
        result = detect_multi_hit_genes([v1, v2], gtex_data=gtex)
        assert len(result) == 1
