"""
L2 Unit Tests — gpa_phaser.py
Phase analysis: cis/trans/ambiguous detection from VCF GT fields and distances.

Run: pytest -m "l2 and phaser" tests/l2_unit/test_phaser.py
"""

import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import make_variant


# =============================================================================
# GT field parser
# =============================================================================

@pytest.mark.l2
@pytest.mark.phaser
@pytest.mark.p0
class TestParseGTField:
    """Test _parse_gt_field: VCF GT string parsing."""

    def test_phased_homozygous(self):
        """PHASER-01: 1|1 → phased, both alleles ALT."""
        from gpa_phaser import _parse_gt_field
        result = _parse_gt_field("1|1")
        assert result["is_phased"] is True
        assert result["allele_0"] == 1
        assert result["allele_1"] == 1

    def test_phased_heterozygous_ref_alt(self):
        """PHASER-02: 0|1 → phased, hap0=REF, hap1=ALT."""
        from gpa_phaser import _parse_gt_field
        result = _parse_gt_field("0|1")
        assert result["is_phased"] is True
        assert result["allele_0"] == 0
        assert result["allele_1"] == 1

    def test_phased_heterozygous_alt_ref(self):
        """PHASER-03: 1|0 → phased, hap0=ALT, hap1=REF."""
        from gpa_phaser import _parse_gt_field
        result = _parse_gt_field("1|0")
        assert result["is_phased"] is True
        assert result["allele_0"] == 1
        assert result["allele_1"] == 0

    def test_unphased_heterozygous(self):
        """PHASER-04: 0/1 → unphased."""
        from gpa_phaser import _parse_gt_field
        result = _parse_gt_field("0/1")
        assert result["is_phased"] is False
        assert result["allele_0"] == 0
        assert result["allele_1"] == 1

    def test_homozygous_ref(self):
        """PHASER-05: 0/0 → unphased, both REF."""
        from gpa_phaser import _parse_gt_field
        result = _parse_gt_field("0/0")
        assert result["is_phased"] is False
        assert result["allele_0"] == 0
        assert result["allele_1"] == 0

    def test_missing_gt(self):
        """PHASER-06: ./. → is_phased=False, alleles=-1."""
        from gpa_phaser import _parse_gt_field
        for missing in [".", "./.", ".|.", "nan", None, ""]:
            result = _parse_gt_field(missing)
            assert result["is_phased"] is False
            assert result["allele_0"] == -1
            assert result["allele_1"] == -1

    def test_haploid(self):
        """PHASER-07: Single allele (haploid)."""
        from gpa_phaser import _parse_gt_field
        result = _parse_gt_field("1")
        assert result["is_phased"] is False
        assert result["allele_0"] == 1
        assert result["allele_1"] == 1


# =============================================================================
# Level 1: GATK phased GT phase determination
# =============================================================================

@pytest.mark.l2
@pytest.mark.phaser
@pytest.mark.p0
class TestLevel1GATKPhase:
    """Test _level1_gatk_phase: phased GT-based phase calls."""

    def test_cis_both_homozygous(self):
        """PHASER-08: Two variants both 1|1 → cis_both."""
        from gpa_phaser import _level1_gatk_phase
        v1 = make_variant(gene="CFTR", gt="1|1", pos=100)
        v2 = make_variant(gene="CFTR", gt="1|1", pos=200)
        result = _level1_gatk_phase([v1, v2])
        assert result is not None
        assert result.phase_status == "cis_both"
        assert result.confidence == "high"
        assert result.method == "gatk_phased_gt"

    def test_cis_hap0_alt(self):
        """PHASER-09: Hap0 all ALT, Hap1 all REF → cis."""
        from gpa_phaser import _level1_gatk_phase
        v1 = make_variant(gene="CFTR", gt="1|0", pos=100)
        v2 = make_variant(gene="CFTR", gt="1|0", pos=200)
        result = _level1_gatk_phase([v1, v2])
        assert result is not None
        assert result.phase_status == "cis"
        assert "Hap0" in result.evidence

    def test_cis_hap1_alt(self):
        """PHASER-10: Hap0 all REF, Hap1 all ALT → cis."""
        from gpa_phaser import _level1_gatk_phase
        v1 = make_variant(gene="CFTR", gt="0|1", pos=100)
        v2 = make_variant(gene="CFTR", gt="0|1", pos=200)
        result = _level1_gatk_phase([v1, v2])
        assert result is not None
        assert result.phase_status == "cis"
        assert "Hap1" in result.evidence

    def test_trans_hap0(self):
        """PHASER-11: Hap0 has both REF and ALT → trans."""
        from gpa_phaser import _level1_gatk_phase
        v1 = make_variant(gene="CFTR", gt="1|0", pos=100)
        v2 = make_variant(gene="CFTR", gt="0|1", pos=200)
        result = _level1_gatk_phase([v1, v2])
        assert result is not None
        assert result.phase_status == "trans"
        assert result.confidence == "high"

    def test_unphased_returns_none(self):
        """PHASER-12: Unphased GT (0/1) → None (fallback to Level 2)."""
        from gpa_phaser import _level1_gatk_phase
        v1 = make_variant(gene="CFTR", gt="0/1", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", pos=200)
        result = _level1_gatk_phase([v1, v2])
        assert result is None

    def test_mixed_phased_unphased_returns_none(self):
        """PHASER-13: Mixed phased/unphased → None."""
        from gpa_phaser import _level1_gatk_phase
        v1 = make_variant(gene="CFTR", gt="1|0", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", pos=200)
        result = _level1_gatk_phase([v1, v2])
        assert result is None


# =============================================================================
# Level 2: Distance assessment
# =============================================================================

@pytest.mark.l2
@pytest.mark.phaser
@pytest.mark.p0
class TestLevel2Distance:
    """Test _level2_distance_assessment: gap-based feasibility."""

    def test_short_reads_overlap(self):
        """PHASER-14: Gap < 50bp → short_reads_overlap, feasible, high confidence."""
        from gpa_phaser import _level2_distance_assessment
        v1 = make_variant(gene="CFTR", pos=100)
        v2 = make_variant(gene="CFTR", pos=130)
        result = _level2_distance_assessment([v1, v2])
        assert result["feasible"] is True
        assert result["confidence"] == "high"
        assert result["method"] == "short_reads_overlap"

    def test_paired_end_or_overlap(self):
        """PHASER-15: Gap 50-150bp → short_reads_overlap_or_paired_end."""
        from gpa_phaser import _level2_distance_assessment
        v1 = make_variant(gene="CFTR", pos=100)
        v2 = make_variant(gene="CFTR", pos=220)
        result = _level2_distance_assessment([v1, v2])
        assert result["feasible"] is True
        assert result["confidence"] == "high"
        assert result["method"] == "short_reads_overlap_or_paired_end"

    def test_paired_end_only(self):
        """PHASER-16: Gap 150-500bp → paired_end_only, medium confidence."""
        from gpa_phaser import _level2_distance_assessment
        v1 = make_variant(gene="CFTR", pos=100)
        v2 = make_variant(gene="CFTR", pos=400)
        result = _level2_distance_assessment([v1, v2])
        assert result["feasible"] is True
        assert result["confidence"] == "medium"
        assert result["method"] == "paired_end_only"

    def test_infeasible_long_distance(self):
        """PHASER-17: Gap > 500bp → infeasible_short_reads."""
        from gpa_phaser import _level2_distance_assessment
        v1 = make_variant(gene="CFTR", pos=100)
        v2 = make_variant(gene="CFTR", pos=1000)
        result = _level2_distance_assessment([v1, v2])
        assert result["feasible"] is False
        assert result["confidence"] == "none"
        assert result["method"] == "infeasible_short_reads"

    def test_single_variant_no_gap(self):
        """PHASER-18: Single variant → gaps empty, feasible=True (default)."""
        from gpa_phaser import _level2_distance_assessment
        v1 = make_variant(gene="CFTR", pos=100)
        result = _level2_distance_assessment([v1])
        assert result["feasible"] is True


# =============================================================================
# Main function: determine_phase
# =============================================================================

@pytest.mark.l2
@pytest.mark.phaser
@pytest.mark.p0
class TestDeterminePhase:
    """Test determine_phase: full phase determination pipeline."""

    def test_phased_trans(self):
        """PHASER-19: Phased trans variants → trans, high confidence."""
        from gpa_phaser import determine_phase
        v1 = make_variant(gene="CFTR", gt="1|0", pos=100)
        v2 = make_variant(gene="CFTR", gt="0|1", pos=200)
        result = determine_phase([v1, v2])
        assert result.phase_status == "trans"
        assert result.confidence == "high"
        assert result.method == "gatk_phased_gt"
        assert result.max_gap_bp == 100
        assert result.n_variants == 2

    def test_unphased_short_gap_ambiguous(self):
        """PHASER-20: Unphased, gap < 50bp → ambiguous (cannot determine cis/trans)."""
        from gpa_phaser import determine_phase
        v1 = make_variant(gene="CFTR", gt="0/1", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", pos=120)
        result = determine_phase([v1, v2])
        assert result.phase_status == "ambiguous"
        assert result.confidence == "high"
        assert result.method == "short_reads_overlap"

    def test_unphased_medium_gap_ambiguous(self):
        """PHASER-21: Unphased, gap 50-150bp → ambiguous, medium confidence."""
        from gpa_phaser import determine_phase
        v1 = make_variant(gene="CFTR", gt="0/1", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", pos=220)
        result = determine_phase([v1, v2])
        assert result.phase_status == "ambiguous"
        assert result.confidence == "medium"

    def test_unphased_long_gap_unphased(self):
        """PHASER-22: Unphased, gap > 500bp → unphased (infeasible)."""
        from gpa_phaser import determine_phase
        v1 = make_variant(gene="CFTR", gt="0/1", pos=100)
        v2 = make_variant(gene="CFTR", gt="0/1", pos=1000)
        result = determine_phase([v1, v2])
        assert result.phase_status == "unphased"
        assert result.confidence == "none"
        assert "trio" in result.evidence.lower() or "long" in result.evidence.lower()

    def test_three_variants_cis_both(self):
        """PHASER-23: Three variants all 1|1 → cis_both."""
        from gpa_phaser import determine_phase
        v1 = make_variant(gene="CFTR", gt="1|1", pos=100)
        v2 = make_variant(gene="CFTR", gt="1|1", pos=200)
        v3 = make_variant(gene="CFTR", gt="1|1", pos=300)
        result = determine_phase([v1, v2, v3])
        assert result.phase_status == "cis_both"
        assert result.n_variants == 3
        assert result.max_gap_bp == 100  # max(200-100, 300-200) = 100

    def test_gap_fields_populated(self):
        """PHASER-24: max_gap_bp and min_gap_bp are correctly calculated."""
        from gpa_phaser import determine_phase
        v1 = make_variant(gene="CFTR", gt="1|0", pos=100)
        v2 = make_variant(gene="CFTR", gt="1|0", pos=500)
        v3 = make_variant(gene="CFTR", gt="1|0", pos=550)
        result = determine_phase([v1, v2, v3])
        assert result.max_gap_bp == 400  # max(500-100, 550-500) = 400
        assert result.min_gap_bp == 50   # 550 - 500
