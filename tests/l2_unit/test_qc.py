"""
L2 Unit Tests — gpa_qc.py
Quality control: VAF range, depth, repeat region, gene symbol, VAF-GT consistency.

Run: pytest -m "l2 and qc" tests/l2_unit/test_qc.py
"""

import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import make_variant


@pytest.mark.l2
@pytest.mark.qc
@pytest.mark.p0
class TestQCChecks:
    """Test _run_qc_checks: quality control flagging."""

    def test_valid_variant_no_flags(self):
        """QC-01: Normal variant → no flags."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=0.45, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert result["total"] == 1
        assert result["flagged"] == 0
        assert v.qc_flags == []

    def test_invalid_vaf_negative(self):
        """QC-02: VAF < 0 → INVALID_VAF."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=-0.1, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_VAF" in v.qc_flags
        assert result["flagged"] == 1

    def test_invalid_vaf_over_one(self):
        """QC-03: VAF > 1 → INVALID_VAF."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=1.5, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_VAF" in v.qc_flags

    def test_low_depth(self):
        """QC-04: DP < 10 → LOW_DEPTH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=0.45, dp=5, pos=100)
        result = _run_qc_checks([v])
        assert "LOW_DEPTH" in v.qc_flags

    def test_sufficient_depth_no_flag(self):
        """QC-05: DP >= 10 → no LOW_DEPTH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=0.45, dp=10, pos=100)
        result = _run_qc_checks([v])
        assert "LOW_DEPTH" not in v.qc_flags

    def test_vaf_gt_mismatch_heterozygous(self):
        """QC-06: GT=0/1, VAF < 0.20 → VAF_GT_MISMATCH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=0.10, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" in v.qc_flags

    def test_vaf_gt_mismatch_heterozygous_high(self):
        """QC-07: GT=0/1, VAF > 0.80 → VAF_GT_MISMATCH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=0.95, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" in v.qc_flags

    def test_vaf_gt_match_heterozygous(self):
        """QC-08: GT=0/1, VAF 0.20-0.80 → no mismatch."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=0.45, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" not in v.qc_flags

    def test_vaf_gt_mismatch_homozygous(self):
        """QC-09: GT=1/1, VAF < 0.70 → VAF_GT_MISMATCH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="1/1", vaf=0.50, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" in v.qc_flags

    def test_vaf_gt_match_homozygous(self):
        """QC-10: GT=1/1, VAF >= 0.70 → no mismatch."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="1/1", vaf=0.85, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" not in v.qc_flags

    def test_vaf_gt_mismatch_ref(self):
        """QC-11: GT=0/0, VAF > 0.10 → VAF_GT_MISMATCH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/0", vaf=0.30, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" in v.qc_flags

    def test_vaf_gt_match_ref(self):
        """QC-12: GT=0/0, VAF <= 0.10 → no mismatch."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/0", vaf=0.05, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" not in v.qc_flags

    def test_invalid_gene_symbol_starts_with_digit(self):
        """QC-13: Gene starts with digit → INVALID_GENE_SYMBOL."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="1ABC", gt="0/1", dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_GENE_SYMBOL" in v.qc_flags

    def test_invalid_gene_symbol_too_long(self):
        """QC-14: Gene > 50 chars → INVALID_GENE_SYMBOL."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="A" * 51, gt="0/1", dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_GENE_SYMBOL" in v.qc_flags

    def test_invalid_gene_symbol_illegal_chars(self):
        """QC-15: Gene with illegal chars → INVALID_GENE_SYMBOL."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="ABC@DEF", gt="0/1", dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_GENE_SYMBOL" in v.qc_flags

    def test_valid_gene_symbol(self):
        """QC-16: Valid gene symbol → no flag."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_GENE_SYMBOL" not in v.qc_flags

    def test_valid_gene_with_hyphen(self):
        """QC-17: Gene with hyphen (e.g., HLA-A) → valid."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="HLA-A", gt="0/1", dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_GENE_SYMBOL" not in v.qc_flags

    def test_multiple_flags_same_variant(self):
        """QC-18: Variant can have multiple flags."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="1ABC", gt="0/1", vaf=-0.5, dp=5, pos=100)
        result = _run_qc_checks([v])
        assert "INVALID_VAF" in v.qc_flags
        assert "LOW_DEPTH" in v.qc_flags
        assert "INVALID_GENE_SYMBOL" in v.qc_flags
        assert result["flagged"] == 1

    def test_summary_counts(self):
        """QC-19: QC summary counts are correct."""
        from gpa_qc import _run_qc_checks
        v1 = make_variant(gene="BRCA1", gt="0/1", vaf=0.45, dp=50, pos=100)  # clean
        v2 = make_variant(gene="BRCA1", gt="./.", vaf=-0.1, dp=50, pos=200)  # INVALID_VAF (missing GT skips VAF check)
        v3 = make_variant(gene="BRCA1", gt="1/1", vaf=0.50, dp=5, pos=300)   # LOW_DEPTH + mismatch
        result = _run_qc_checks([v1, v2, v3])
        assert result["total"] == 3
        assert result["flagged"] == 2
        assert result["by_flag"]["INVALID_VAF"] == 1
        assert result["by_flag"]["LOW_DEPTH"] == 1
        assert result["by_flag"]["VAF_GT_MISMATCH"] == 1
        assert len(result["flagged_variants"]) == 2

    def test_missing_gt_skips_vaf_check(self):
        """QC-20: Missing GT → no VAF_GT_MISMATCH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt=".", vaf=0.45, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" not in v.qc_flags

    def test_none_vaf_skips_vaf_check(self):
        """QC-21: None VAF → no VAF_GT_MISMATCH."""
        from gpa_qc import _run_qc_checks
        v = make_variant(gene="BRCA1", gt="0/1", vaf=None, dp=50, pos=100)
        result = _run_qc_checks([v])
        assert "VAF_GT_MISMATCH" not in v.qc_flags
        assert "INVALID_VAF" not in v.qc_flags

    def test_empty_variant_list(self):
        """QC-22: Empty list → summary with total=0."""
        from gpa_qc import _run_qc_checks
        result = _run_qc_checks([])
        assert result["total"] == 0
        assert result["flagged"] == 0
