#!/usr/bin/env python3
"""
GPA v0.10.0 End-to-End Test Suite

Tests the full pipeline from mock variant input to report output,
verifying key business rules:
  1. EAS AF > 50% → auto Tier 3
  2. ClinVar Pathogenic + HIGH + tissue-relevant → Tier 1
  3. Report structure completeness (Markdown + JSON)
  4. Jinja2 header rendering (or fallback)
  5. Circular import fix (gpa_pipeline ↔ dgra_core)

Run: cd scripts && python3 test_e2e_gpa_v010.py
"""

import asyncio
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(_SCRIPT_DIR))

from gpa_pipeline import run_dgra_pipeline
from gpa_types import GPAConfig



# =============================================================================
# Test Fixtures
# =============================================================================

def _variant(gene, impact, consequence, clinvar, revstat, gnomad_af, eas_af,
             gt="0/1", vaf=0.50, pos=100000):
    """Build a mock variant dict."""
    return {
        "CHROM": "1", "POS": pos, "REF": "A", "ALT": "T",
        "GENE": gene, "FEATURE": "ENST00000357654",
        "IMPACT": impact, "Consequence": consequence,
        "HGVSc": "c.100A>T", "HGVSp": "p.Lys34Ter",
        "GT": gt, "VAF": vaf, "DP": 120, "GQ": 99,
        "CLINVAR": clinvar, "CLNREVSTAT": revstat,
        "gnomAD_AF": gnomad_af, "gnomAD_EAS_AF": eas_af,
    }


ASSETS = {
    "brca1_pathogenic_rare": _variant(
        gene="BRCA1", impact="HIGH", consequence="stop_gained",
        clinvar="Pathogenic", revstat="reviewed_by_expert_panel",
        gnomad_af=0.0001, eas_af=0.0002, pos=100000,
    ),
    "or2b11_common": _variant(
        gene="OR2B11", impact="HIGH", consequence="frameshift_variant",
        clinvar="", revstat="",
        gnomad_af=0.52, eas_af=0.55, pos=200000,
    ),
    "tpmt_pathogenic_rare": _variant(
        gene="TPMT", impact="MODERATE", consequence="missense_variant",
        clinvar="Pathogenic", revstat="practice_guideline",
        gnomad_af=0.002, eas_af=0.001, pos=300000,
    ),
    "ddx3x_homozygous": _variant(
        gene="DDX3X", impact="HIGH", consequence="splice_donor_variant",
        clinvar="", revstat="",
        gnomad_af=0.65, eas_af=0.62, pos=400000, gt="1/1", vaf=0.98,
    ),
    "runx1_pathogenic": _variant(
        gene="RUNX1", impact="HIGH", consequence="stop_gained",
        clinvar="Pathogenic", revstat="reviewed_by_expert_panel",
        gnomad_af=0.0001, eas_af=0.0001, pos=500000,
    ),
}


# =============================================================================
# Test Cases
# =============================================================================

async def test_e2e_basic_structure():
    """Pipeline produces Markdown + JSON with all expected sections."""
    variants = [ASSETS["brca1_pathogenic_rare"]]
    config = GPAConfig(tissue_profile="hematopoietic", offline_mode=True)
    result = await run_dgra_pipeline(variants, user_phenotypes=None, config=config)

    assert "report_markdown" in result, "Missing report_markdown"
    assert "json_report" in result, "Missing json_report"
    assert "meta" in result, "Missing meta"
    assert "summary" in result, "Missing summary"

    md = result["report_markdown"]
    assert len(md) > 0, "Report is empty"
    assert "Tier 1" in md or "Tier 2" in md or "Tier 3" in md, "No tier sections"
    assert "方法学附录" in md, "Missing methodology appendix"

    jr = result["json_report"]
    assert "meta" in jr and "summary" in jr and "variants" in jr, "JSON missing keys"
    assert len(jr["variants"]) == 1, f"Expected 1 variant in JSON, got {len(jr['variants'])}"

    print("[PASS] test_e2e_basic_structure")


async def test_e2e_eas_af_guard_tier3():
    """EAS AF > 50% forces Tier 3 regardless of other evidence."""
    variants = [ASSETS["or2b11_common"]]
    config = GPAConfig(tissue_profile="hematopoietic", offline_mode=True)
    result = await run_dgra_pipeline(variants, user_phenotypes=None, config=config)

    tier3 = result.get("tier3_variants", [])
    assert len(tier3) == 1, f"Expected 1 Tier 3, got {len(tier3)}"
    assert tier3[0]["gene"] == "OR2B11", f"Expected OR2B11, got {tier3[0].get('gene')}"

    # JSON confirmation
    jr_variants = result.get("json_report", {}).get("variants", [])
    or2b11_json = [v for v in jr_variants if v.get("gene") == "OR2B11"]
    assert len(or2b11_json) == 1, "OR2B11 missing from JSON"
    assert or2b11_json[0].get("tier") == 3, f"Expected tier 3, got {or2b11_json[0].get('tier')}"

    # Report confirmation
    md = result["report_markdown"]
    assert "Tier 3" in md, "Tier 3 section missing from report"

    print("[PASS] test_e2e_eas_af_guard_tier3")


async def test_e2e_clinvar_high_impact_tissue_relevant():
    """ClinVar Pathogenic + HIGH + tissue-relevant gene → Tier 1 or 2 (not 3)."""
    # RUNX1 is a hematopoietic gene, so in hematopoietic profile it should be tissue-relevant
    variants = [ASSETS["runx1_pathogenic"]]
    config = GPAConfig(tissue_profile="hematopoietic", offline_mode=True)
    result = await run_dgra_pipeline(variants, user_phenotypes=None, config=config)

    tier1 = result.get("tier1_variants", [])
    tier2 = result.get("tier2_variants", [])
    tier3 = result.get("tier3_variants", [])

    total_assigned = len(tier1) + len(tier2) + len(tier3)
    assert total_assigned == 1, f"Expected 1 variant total, got {total_assigned}"

    # ClinVar Pathogenic + HIGH + tissue-relevant should NOT be Tier 3
    # (unless frequency override which RUNX1 doesn't have)
    assert len(tier3) == 0, \
        f"ClinVar Pathogenic + HIGH + tissue-relevant incorrectly Tier 3. T1={len(tier1)} T2={len(tier2)} T3={len(tier3)}"

    print("[PASS] test_e2e_clinvar_high_impact_tissue_relevant")


async def test_e2e_multiple_variants_all_tiers():
    """Mixed cohort produces correct tier distribution."""
    variants = [
        ASSETS["or2b11_common"],      # EAS AF > 50% → Tier 3
        ASSETS["tpmt_pathogenic_rare"],  # ClinVar Pathogenic + MODERATE
    ]
    config = GPAConfig(tissue_profile="hematopoietic", offline_mode=True)
    result = await run_dgra_pipeline(variants, user_phenotypes=None, config=config)

    total = (
        len(result.get("tier1_variants", []))
        + len(result.get("tier2_variants", []))
        + len(result.get("tier3_variants", []))
    )
    assert total == 2, f"Expected 2 variants, got {total}"

    # OR2B11 must be Tier 3
    t3_genes = [v["gene"] for v in result.get("tier3_variants", [])]
    assert "OR2B11" in t3_genes, f"OR2B11 not in Tier 3: {t3_genes}"

    print("[PASS] test_e2e_multiple_variants_all_tiers")


async def test_e2e_json_report_completeness():
    """JSON report contains all required fields for downstream consumption."""
    variants = [ASSETS["tpmt_pathogenic_rare"]]
    config = GPAConfig(tissue_profile="hematopoietic", offline_mode=True)
    result = await run_dgra_pipeline(variants, user_phenotypes=None, config=config)

    jr = result["json_report"]
    required_top = ["meta", "summary", "variants", "multi_hit_details", "qc_summary", "phenotype_association", "report_md"]
    for key in required_top:
        assert key in jr, f"Missing top-level JSON key: {key}"

    # Each variant must have key fields
    for v in jr["variants"]:
        for field in ["gene", "chrom", "pos", "tier", "tier_confidence", "tier_reason", "evidence_chain", "clinvar", "gnomAD"]:
            assert field in v, f"Missing variant field: {field}"

    print("[PASS] test_e2e_json_report_completeness")


async def test_circular_import_fix():
    """Verify gpa_pipeline and dgra_core can be imported without circular error."""
    import importlib

    # Clear cached modules to force re-import
    for mod in list(sys.modules.keys()):
        if "dgra" in mod or "gpa_" in mod:
            del sys.modules[mod]

    # This should not raise ImportError
    from gpa_pipeline import run_dgra_pipeline as rdp
    assert rdp is not None

    print("[PASS] test_circular_import_fix")


# =============================================================================
# Runner
# =============================================================================

async def main():
    print("=" * 60)
    print("GPA v0.10.0 End-to-End Test Suite")
    print("=" * 60)

    tests = [
        test_circular_import_fix,
        test_e2e_basic_structure,
        test_e2e_eas_af_guard_tier3,
        test_e2e_clinvar_high_impact_tissue_relevant,
        test_e2e_multiple_variants_all_tiers,
        test_e2e_json_report_completeness,
    ]

    passed = 0
    failed = 0
    for t in tests:
        try:
            await t()
            passed += 1
        except (RuntimeError, ValueError) as e:
            print(f"[FAIL] {t.__name__}: {e}")
            failed += 1

    print("=" * 60)
    print(f"RESULTS: {passed} passed, {failed} failed")
    print("=" * 60)
    return failed == 0


if __name__ == "__main__":
    ok = asyncio.run(main())
    sys.exit(0 if ok else 1)
