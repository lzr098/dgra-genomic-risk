#!/usr/bin/env python3
"""
HLA Region Reliability Warning Tests — v0.11.4

Tests the P0 HLA region marking strategy:
1. Coordinate-based HLA region detection
2. qc_flag injection in classify_variant_tier
3. Tier 1/2 action injection in pipeline
4. Report section generation
"""

import sys
sys.path.insert(0, "/Users/zhaorongli/.workbuddy/skills/dgra-genomic-risk/scripts")

from dgra_core import _is_hla_region, Variant
from gpa_tier_classifier import classify_variant_tier
from gpa_report import _generate_hla_region_section


def test_is_hla_region_grch38():
    """HLA region detection for GRCh38."""
    # Inside HLA (classic positions)
    assert _is_hla_region("chr6", 30_000_000, "GRCh38") is True
    assert _is_hla_region("6", 31_000_000, "GRCh38") is True
    assert _is_hla_region("chr06", 32_000_000, "GRCh38") is True

    # Boundaries
    assert _is_hla_region("chr6", 29_690_683, "GRCh38") is True   # start
    assert _is_hla_region("chr6", 33_079_845, "GRCh38") is True   # end

    # Outside HLA
    assert _is_hla_region("chr6", 29_000_000, "GRCh38") is False  # before
    assert _is_hla_region("chr6", 34_000_000, "GRCh38") is False  # after
    assert _is_hla_region("chr1", 30_000_000, "GRCh38") is False  # wrong chrom
    assert _is_hla_region("chrX", 30_000_000, "GRCh38") is False  # wrong chrom


def test_is_hla_region_grch37():
    """HLA region detection for GRCh37."""
    assert _is_hla_region("chr6", 30_000_000, "GRCh37") is True
    assert _is_hla_region("6", 30_000_000, "GRCh37") is True

    # GRCh37 boundaries
    assert _is_hla_region("chr6", 29_941_030, "GRCh37") is True
    assert _is_hla_region("chr6", 33_330_192, "GRCh37") is True

    # Outside
    assert _is_hla_region("chr6", 29_000_000, "GRCh37") is False


def test_is_hla_region_default_build():
    """Default genome build is GRCh38."""
    assert _is_hla_region("chr6", 31_000_000) is True
    assert _is_hla_region("chr6", 29_000_000) is False


def test_tier_classifier_adds_hla_flag():
    """classify_variant_tier adds HLA_REGION_SHORT_READ_UNRELIABLE flag."""
    # HLA-region variant
    v_hla = Variant(
        chrom="chr6", pos=31_000_000, ref="A", alt="G",
        gene="HLA-A", transcript="ENST00000376809",
        exon="E3/8", impact="MODERATE", consequence="missense_variant",
        hgvsp="p.Ala123Val", hgvsc="c.368C>T",
        clinvar="", gnomad_af=0.001,
    )
    # Non-HLA variant
    v_non = Variant(
        chrom="chr17", pos=43_045_752, ref="T", alt="G",
        gene="BRCA1", transcript="ENST00000357654",
        exon="E5/22", impact="HIGH", consequence="frameshift_variant",
        hgvsp="p.Ser1253ArgfsTer5", hgvsc="c.3756_3759del",
        clinvar="Pathogenic", gnomad_af=0.0001,
    )

    tier_hla, _, _ = classify_variant_tier(
        v_hla, {}, {"relevance": "none"}, {"af": 0.001}, None, None, {}
    )
    tier_non, _, _ = classify_variant_tier(
        v_non, {}, {"relevance": "primary"}, {"af": 0.0001}, None, None, {}
    )

    assert "HLA_REGION_SHORT_READ_UNRELIABLE" in v_hla.qc_flags
    assert "HLA_REGION_SHORT_READ_UNRELIABLE" not in v_non.qc_flags


def test_tier_classifier_hla_flag_does_not_change_tier():
    """HLA flag should NOT alter the tier result."""
    v = Variant(
        chrom="chr6", pos=31_000_000, ref="A", alt="G",
        gene="HLA-A", transcript="ENST00000376809",
        exon="E3/8", impact="HIGH", consequence="frameshift_variant",
        hgvsp="p.Ala123ValfsTer5", hgvsc="c.368del",
        clinvar="Pathogenic", gnomad_af=0.0001,
    )
    tier, reason, actions = classify_variant_tier(
        v, {}, {"relevance": "primary"}, {"af": 0.0001}, None, None, {}
    )
    # HIGH + ClinVar pathogenic + primary relevance → Tier 1
    assert tier == 1
    assert "HLA_REGION_SHORT_READ_UNRELIABLE" in v.qc_flags


def test_report_section_with_hla_tier12():
    """Report section highlights Tier 1/2 HLA variants."""
    variants = [
        Variant(
            chrom="chr6", pos=31_000_000, ref="A", alt="G",
            gene="HLA-A", transcript="ENST00000376809",
            exon="E3/8", impact="HIGH", consequence="frameshift_variant",
            hgvsp="p.Ala123ValfsTer5", hgvsc="c.368del",
            clinvar="Pathogenic", gnomad_af=0.0001,
            tier=1, qc_flags=["HLA_REGION_SHORT_READ_UNRELIABLE"],
        ),
        Variant(
            chrom="chr6", pos=31_500_000, ref="C", alt="T",
            gene="HLA-B", transcript="ENST00000376809",
            exon="E2/7", impact="MODERATE", consequence="missense_variant",
            hgvsp="p.Thr456Met", hgvsc="c.1367C>T",
            clinvar="", gnomad_af=0.5,
            tier=3, qc_flags=["HLA_REGION_SHORT_READ_UNRELIABLE", "BENIGN_POLYMORPHISM_FREQUENCY"],
        ),
    ]

    section = _generate_hla_region_section(variants)
    assert section is not None
    assert "HLA 区域变异可靠性提示" in section
    assert "HLA-A" in section
    # HLA-B is Tier 3 — not in the highlight table, only in count
    assert "Tier 1" in section
    assert "Tier 3" not in section  # Tier 3 not highlighted
    assert "建议验证方法" in section
    assert "HLA-HD" in section


def test_report_section_no_hla_variants():
    """Report section returns None when no HLA variants."""
    variants = [
        Variant(
            chrom="chr17", pos=43_045_752, ref="T", alt="G",
            gene="BRCA1", transcript="ENST00000357654",
            exon="E5/22", impact="HIGH", consequence="frameshift_variant",
            hgvsp="p.Ser1253ArgfsTer5", hgvsc="c.3756_3759del",
            clinvar="Pathogenic", gnomad_af=0.0001,
            tier=1, qc_flags=[],
        ),
    ]
    assert _generate_hla_region_section(variants) is None


def test_report_section_empty_list():
    """Report section returns None for empty variant list."""
    assert _generate_hla_region_section([]) is None


def test_hla_chr_prefix_variants():
    """HLA detection works with various chromosome formats."""
    assert _is_hla_region("chr6", 31_000_000) is True
    assert _is_hla_region("6", 31_000_000) is True
    assert _is_hla_region("chr06", 31_000_000) is True
    assert _is_hla_region("06", 31_000_000) is True


if __name__ == "__main__":
    print("=" * 60)
    print("HLA Region Reliability Warning Tests — v0.11.4")
    print("=" * 60)

    tests = [
        test_is_hla_region_grch38,
        test_is_hla_region_grch37,
        test_is_hla_region_default_build,
        test_tier_classifier_adds_hla_flag,
        test_tier_classifier_hla_flag_does_not_change_tier,
        test_report_section_with_hla_tier12,
        test_report_section_no_hla_variants,
        test_report_section_empty_list,
        test_hla_chr_prefix_variants,
    ]

    passed = 0
    failed = 0
    for t in tests:
        try:
            t()
            print(f"  ✅ {t.__name__}")
            passed += 1
        except AssertionError as e:
            print(f"  ❌ {t.__name__}: {e}")
            failed += 1
        except Exception as e:
            print(f"  ❌ {t.__name__}: {type(e).__name__}: {e}")
            failed += 1

    print(f"\nHLA: {passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
