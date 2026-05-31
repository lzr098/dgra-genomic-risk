"""
L2 Unit Tests — gpa_tier_classifier.py
Core tier classification rules with mocked dependencies.

Run: pytest -m "l2 and tier" tests/l2_unit/test_tier_classifier.py
"""

import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import (
    MockGnomAD,
    MockTissueAssessment,
    MockTissueProfile,
    make_variant,
)


# =============================================================================
# Tier 1 — Must-intervene pathogenic variants
# =============================================================================

@pytest.mark.l2
@pytest.mark.tier
@pytest.mark.p0
@pytest.mark.hematopoietic
class TestTier1:
    """Tier 1 classification: must-intervene variants."""

    def test_clinvar_pathogenic_high_primary_hom_tier1(self):
        """TIER-01: ClinVar Pathogenic + HIGH + primary tissue + homozygous → Tier 1."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="RUNX1", gt="1/1", impact="HIGH",
            clinvar="Pathogenic",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 1, f"Expected Tier 1, got {tier}. Reason: {reason}"

    def test_clinvar_pathogenic_high_primary_het_tier1(self):
        """TIER-02: ClinVar Pathogenic + HIGH + primary tissue + heterozygous → Tier 1.

        v0.5.2 FIX: Heterozygous pathogenic truncating variants in tissue-relevant
        genes are Tier 1 regardless of zygosity (heterozygous pathogenic = actionable).
        """
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="RUNX1", gt="0/1", impact="HIGH",
            clinvar="Pathogenic",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 1, f"Expected Tier 1 for heterozygous pathogenic, got {tier}. Reason: {reason}"

    def test_homozygous_lof_primary_rare_tier1(self):
        """TIER-08: Homozygous LoF + primary tissue + rare AF → Tier 1."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="CFTR", gt="1/1", impact="HIGH",
            gnomad_af=0.0001, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.success_rare("1", 100, "A", "G", af=0.0001)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 1, f"Expected Tier 1 for homozygous LoF, got {tier}. Reason: {reason}"

    def test_not_captured_homozygous_tier1(self):
        """TIER-11: NOT_CAPTURED + homozygous → Tier 1."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="RUNX1", gt="1/1", impact="HIGH",
            gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 1, f"Expected Tier 1 for NOT_CAPTURED hom, got {tier}. Reason: {reason}"

    def test_somatic_tsg_lof_tier1(self):
        """TIER-26: Somatic TSG LoF → Tier 1."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="TP53", gt="0/1", impact="HIGH",
            gnomad_af=0.0001, gnomad_status="SUCCESS",
        )
        # Mark as TSG for somatic mode
        v.is_tsg = "Yes"
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general", somatic_mode=True)
        gnomad_info = MockGnomAD.success_rare("17", 7579472, "C", "T", af=0.0001)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 1, f"Expected Tier 1 for somatic TSG LoF, got {tier}. Reason: {reason}"


# =============================================================================
# Tier 2 — Monitor / carrier variants
# =============================================================================

@pytest.mark.l2
@pytest.mark.tier
@pytest.mark.p0
@pytest.mark.hematopoietic
class TestTier2:
    """Tier 2 classification: monitor or carrier variants."""

    def test_heterozygous_lof_primary_tier2(self):
        """TIER-09: Heterozygous LoF + primary tissue → Tier 2."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="RUNX1", gt="0/1", impact="HIGH",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 2, f"Expected Tier 2 for heterozygous LoF, got {tier}. Reason: {reason}"

    def test_gnomad_api_failed_downgrade_tier2(self):
        """TIER-10: gnomAD API_FAILED → downgrade Tier 1 candidate to Tier 2."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="DDX3X", gt="1/1", impact="HIGH",
            gnomad_status="API_FAILED", gnomad_error_msg="GraphQL 400",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.api_failed("X", 41357831, "A", "T")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 2, f"Expected Tier 2 for API_FAILED, got {tier}. Reason: {reason}"
        assert "API_FAILED" in reason or any("Downgraded" in str(a) for a in actions)

    def test_likely_pathogenic_fa_pathway_tier1(self):
        """TIER-03: ClinVar Likely_pathogenic + HIGH + FA pathway (BRCA1) → Tier 1.

        _clinvar_pathogenic() treats Likely_pathogenic as True.
        FA pathway genes (BRCA1, BRCA2, PALB2, etc.) auto-Tier 1 for pathogenic
        variants due to marrow failure risk — priority over standard tier rules.
        """
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="BRCA1", gt="0/1", impact="HIGH",
            clinvar="Likely_pathogenic",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 1, f"Expected Tier 1 for FA pathway likely_pathogenic, got {tier}. Reason: {reason}"

    def test_pharmacogenomics_drug_metabolism_tier2(self):
        """TIER-27: Drug metabolism gene polymorphism → Tier 2.

        Uses EAS AF=0.005 (< 1% threshold) to bypass frequency guard and reach
        the drug_metabolism special-list handler (line 830-837).
        """
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="CYP2D6", gt="0/1", impact="MODERATE",
            consequence="missense_variant",
            clinvar="Uncertain_significance",
            gnomad_af=0.005, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.secondary()
        profile = MockTissueProfile.general()
        profile["special_gene_lists"]["drug_metabolism"] = {"CYP2D6", "TPMT", "DPYD"}
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.success_rare("1", 100, "A", "G", af=0.005)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 2, f"Expected Tier 2 for drug metabolism, got {tier}. Reason: {reason}"


# =============================================================================
# Tier 3 — Benign / common / no relevance
# =============================================================================

@pytest.mark.l2
@pytest.mark.tier
@pytest.mark.p0
class TestTier3:
    """Tier 3 classification: benign, common, or no tissue relevance."""

    def test_eas_af_over_50_percent_tier3(self):
        """TIER-06: EAS AF > 50% → Tier 3."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="OR2B11", gt="0/1", impact="HIGH",
            gnomad_af=0.52, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.none()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.success_common("1", 100, "A", "G", af=0.52)
        gnomad_info["af_populations"]["EAS"] = {"af": 0.55, "ac": 100, "an": 222}

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 3, f"Expected Tier 3 for EAS AF>50%, got {tier}. Reason: {reason}"

    def test_global_af_over_80_percent_tier3(self):
        """TIER-07: Global AF > 80% → Tier 3."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="MAD2L2", gt="0/1", impact="HIGH",
            gnomad_af=0.85, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.none()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.success_common("1", 100, "A", "G", af=0.85)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 3, f"Expected Tier 3 for global AF>80%, got {tier}. Reason: {reason}"

    def test_clinvar_benign_tier3(self):
        """TIER-04: ClinVar Benign → Tier 3."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="CFTR", gt="0/1", impact="HIGH",
            clinvar="Benign",
            gnomad_af=0.02, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.success_common("1", 100, "A", "G", af=0.02)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 3, f"Expected Tier 3 for Benign, got {tier}. Reason: {reason}"

    def test_no_tissue_relevance_fast_track_tier3(self):
        """TIER-28: No tissue relevance → fast-track Tier 3."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="BRCA1", gt="0/1", impact="HIGH",
            gnomad_af=0.45, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.none()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.success_common("1", 100, "A", "G", af=0.45)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 3, f"Expected Tier 3 for no relevance, got {tier}. Reason: {reason}"

    def test_somatic_vaf_over_50_tier3(self):
        """TIER-25: Somatic VAF > 0.5 → Tier 3 (germline contamination)."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="TP53", gt="0/1", impact="HIGH",
            vaf=0.98,
            gnomad_af=0.0001, gnomad_status="SUCCESS",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general", somatic_mode=True)
        gnomad_info = MockGnomAD.success_rare("17", 7579472, "C", "T", af=0.0001)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier == 3, f"Expected Tier 3 for somatic VAF>0.5, got {tier}. Reason: {reason}"


# =============================================================================
# Special cases
# =============================================================================

@pytest.mark.l2
@pytest.mark.tier
@pytest.mark.p0
class TestSpecialCases:
    """Special classification rules: conflicting, NMD, SpliceAI, etc."""

    def test_clinvar_conflicting_not_upgraded(self):
        """TIER-05: ClinVar Conflicting → not upgraded, flagged."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="BRCA1", gt="0/1", impact="HIGH",
            clinvar="Conflicting interpretations of pathogenicity",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert tier != 1, f"Conflicting should NOT be Tier 1, got {tier}"
        assert "CLINVAR_CONFLICTING" in v.qc_flags or "CONFLICTING" in str(v.qc_flags)

    def test_review_status_practice_guideline_high_weight(self):
        """TIER-14: practice_guideline → high confidence weight."""
        from dgra_core import classify_variant_tier, GPAConfig, Evidence

        v = make_variant(
            gene="BRCA1", gt="0/1", impact="HIGH",
            clinvar="Pathogenic",
            clinvar_review_status="practice_guideline",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        # practice_guideline should result in high-confidence Tier 1 or 2
        assert tier in (1, 2), f"Expected Tier 1/2 for practice_guideline, got {tier}"
        # Check evidence weight
        clinvar_evidence = [e for e in v.evidence_chain if e.source == "ClinVar"]
        if clinvar_evidence:
            assert clinvar_evidence[0].weight >= 0.9, "practice_guideline should have high weight"

    def test_review_status_single_submitter_tier1_low_confidence(self):
        """TIER-15: single_submitter → Tier 1 but with reduced confidence.

        Review status affects evidence weight (0.40 for single_submitter) and
        tier_confidence, but does NOT prevent Tier 1 for Pathogenic + HIGH +
        tissue-relevant variants.
        """
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="RUNX1", gt="0/1", impact="HIGH",
            clinvar="Pathogenic",
            clinvar_review_status="criteria_provided,_single_submitter",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        # single_submitter Pathogenic + HIGH + tissue-relevant = still Tier 1
        assert tier == 1, f"Expected Tier 1 for single_submitter Pathogenic, got {tier}. Reason: {reason}"
        # But confidence should be lower
        assert v.tier_confidence in ("LOW", "MEDIUM"), f"single_submitter should have LOW/MEDIUM confidence, got {v.tier_confidence}"

    def test_phenotype_match_fa_pathway_priority_tier1(self):
        """TIER-12: FA pathway genes override phenotype match downgrade.

        BRCA1 is in the fa_dna_repair special list; pathogenic variants in FA
        pathway genes are Tier 1 regardless of phenotype match score.
        Phenotype match downgrade (score < 0.6) applies only to the standard
        Pathogenic + HIGH + tissue-relevant path (line 619-633), which is
        evaluated AFTER the FA pathway check (line 563-570).
        """
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="BRCA1", gt="0/1", impact="HIGH",
            clinvar="Pathogenic",
            phenotype_match_score=0.3,
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        # FA pathway priority: still Tier 1 despite low phenotype match
        assert tier == 1, f"Expected Tier 1 (FA pathway priority), got {tier}. Reason: {reason}"

    def test_unknown_impact_conservative_high(self):
        """TIER-13: UNKNOWN impact → conservatively treated as HIGH."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="TEST1", gt="1/1", impact="UNKNOWN",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        # UNKNOWN treated as HIGH + homozygous + primary → should be Tier 1 or 2
        assert tier in (1, 2), f"UNKNOWN impact should be conservatively HIGH, got tier {tier}"

    def test_evidence_chain_populated(self):
        """TIER-29: evidence_chain is populated after classification."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="RUNX1", gt="0/1", impact="HIGH",
            clinvar="Pathogenic",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.hematopoietic()
        config = GPAConfig(tissue_profile="hematopoietic")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        assert len(v.evidence_chain) > 0, "evidence_chain should not be empty after classification"

    def test_qc_flags_set(self):
        """TIER-30: qc_flags are correctly set for special cases."""
        from dgra_core import classify_variant_tier, GPAConfig

        v = make_variant(
            gene="BRCA1", gt="0/1", impact="HIGH",
            clinvar="Conflicting interpretations of pathogenicity",
            gnomad_af=None, gnomad_status="NOT_CAPTURED",
        )
        tissue = MockTissueAssessment.primary()
        profile = MockTissueProfile.general()
        config = GPAConfig(tissue_profile="general")
        gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

        classify_variant_tier(v, {}, tissue, gnomad_info, {}, None, profile, config)
        assert len(v.qc_flags) > 0, "qc_flags should be set for conflicting ClinVar"
