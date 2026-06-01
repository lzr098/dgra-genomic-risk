#!/usr/bin/env python3
"""
Biotype-aware Transcript Selection Tests — v0.11.5

Tests three fixes from the diagnostic report:
1. VEP API response parses biotype field
2. Local VEP CSQ parsing respects BIOTYPE for transcript selection
3. TranscriptSelector scores penalize NMD/pseudogene and boost canonical
"""

import sys
sys.path.insert(0, "/Users/zhaorongli/.workbuddy/skills/dgra-genomic-risk/scripts")

from gpa_transcript_selector import TranscriptSelector
from dgra_input_parsers import VCFParser


def test_vep_api_biotype_extracted():
    """Simulated VEP API response includes biotype in parsed output."""
    from gpa_vcf_annotator import VCFAnnotator

    # Mock VEP response with biotype
    mock_response = [
        {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST00000372759",
                    "gene_symbol": "ZMPSTE24",
                    "biotype": "protein_coding",
                    "impact": "MODIFIER",
                    "canonical": 1,
                },
                {
                    "transcript_id": "ENST00000674703",
                    "gene_symbol": "ZMPSTE24",
                    "biotype": "nonsense_mediated_decay",
                    "impact": "HIGH",
                    "canonical": 0,
                },
            ]
        }
    ]
    mock_batch = [{"chrom": "1", "pos": 40262857, "ref": "G", "alt": "A"}]

    results = VCFAnnotator._parse_vep_response(mock_response, mock_batch)
    txs = results[0]["transcript_consequences"]
    assert txs[0]["biotype"] == "protein_coding"
    assert txs[1]["biotype"] == "nonsense_mediated_decay"


def test_selector_prefers_protein_coding_over_nmd():
    """ZMPSTE24 scenario: canonical protein_coding vs NMD HIGH impact."""
    selector = TranscriptSelector()

    transcripts = [
        {
            "transcript_id": "ENST00000372759",
            "gene_symbol": "ZMPSTE24",
            "biotype": "protein_coding",
            "impact": "MODIFIER",
            "canonical": 1,
            "mane_select": 0,
            "protein_domains": [],
        },
        {
            "transcript_id": "ENST00000674703",
            "gene_symbol": "ZMPSTE24",
            "biotype": "nonsense_mediated_decay",
            "impact": "HIGH",
            "canonical": 0,
            "mane_select": 0,
            "protein_domains": [],
        },
    ]

    result = selector.select("ZMPSTE24", transcripts)
    assert result.primary["transcript_id"] == "ENST00000372759"
    assert result.primary["biotype"] == "protein_coding"
    assert not result.is_ambiguous


def test_selector_nmd_pseudogene_penalty():
    """NMD and pseudogene biotypes receive large negative penalties."""
    selector = TranscriptSelector()

    score_pc, _ = selector._score_transcript(
        {"biotype": "protein_coding", "canonical": 0, "impact": "MODERATE", "protein_domains": []},
        "TEST"
    )
    score_nmd, _ = selector._score_transcript(
        {"biotype": "nonsense_mediated_decay", "canonical": 0, "impact": "HIGH", "protein_domains": []},
        "TEST"
    )
    score_pseudo, _ = selector._score_transcript(
        {"biotype": "processed_pseudogene", "canonical": 1, "impact": "HIGH", "protein_domains": []},
        "TEST"
    )

    # protein_coding gets +5 bonus
    assert score_pc >= 3  # +5 protein_coding + 3 MODERATE = 8
    # NMD gets -20 penalty, HIGH impact only +5
    assert score_nmd < score_pc
    # Even canonical (+15) can't save a processed_pseudogene (-25)
    assert score_pseudo < 0


def test_selector_canonical_beats_high_impact():
    """canonical (+15) + protein_coding (+5) = 20 should beat NMD HIGH (+5 -20 = -15)."""
    selector = TranscriptSelector()

    score_canonical, _ = selector._score_transcript(
        {"biotype": "protein_coding", "canonical": 1, "impact": "MODIFIER", "protein_domains": []},
        "TEST"
    )
    score_nmd_high, _ = selector._score_transcript(
        {"biotype": "nonsense_mediated_decay", "canonical": 0, "impact": "HIGH", "protein_domains": []},
        "TEST"
    )

    assert score_canonical > score_nmd_high
    assert score_canonical >= 15  # 15 canonical + 5 protein_coding + 1 LOW/MODIFIER


def test_local_vcf_prefers_protein_coding_canonical():
    """Local VCF parser _pick_csq prefers CANONICAL + protein_coding."""
    parser = VCFParser()

    csq_map = {
        "Allele": 0, "Consequence": 1, "IMPACT": 2, "SYMBOL": 3,
        "Feature": 4, "Feature_type": 5, "EXON": 6, "INTRON": 7,
        "HGVSc": 8, "HGVSp": 9, "CANONICAL": 10, "MANE_SELECT": 11,
        "BIOTYPE": 12,
    }

    # Entry 0: canonical=YES, biotype=processed_pseudogene
    # Entry 1: canonical=NO,  biotype=protein_coding
    # Entry 2: canonical=YES, biotype=protein_coding
    entries = [
        ["A", "missense_variant", "MODERATE", "ZMPSTE24", "ENST00000447743", "Transcript", "3/8", "", "c.123A>G", "p.Arg41His", "YES", "", "processed_pseudogene"],
        ["A", "synonymous_variant", "LOW", "ZMPSTE24", "ENST00000567508", "Transcript", "2/8", "", "c.100C>T", "p.Pro34Pro", "", "", "protein_coding"],
        ["A", "missense_variant", "MODERATE", "ZMPSTE24", "ENST00000372759", "Transcript", "3/8", "", "c.123A>G", "p.Arg41His", "YES", "", "protein_coding"],
    ]

    chosen = parser._pick_csq(entries, csq_map)
    assert chosen[4] == "ENST00000372759"  # canonical + protein_coding


def test_local_vcf_fallback_to_any_protein_coding():
    """If no canonical+protein_coding, fallback to any protein_coding."""
    parser = VCFParser()

    csq_map = {
        "Allele": 0, "Consequence": 1, "IMPACT": 2, "SYMBOL": 3,
        "Feature": 4, "Feature_type": 5, "EXON": 6, "INTRON": 7,
        "HGVSc": 8, "HGVSp": 9, "CANONICAL": 10, "MANE_SELECT": 11,
        "BIOTYPE": 12,
    }

    entries = [
        ["A", "missense_variant", "MODERATE", "GENE", "ENST001", "Transcript", "1/5", "", "c.1A>G", "p.Met1Val", "YES", "", "processed_pseudogene"],
        ["A", "synonymous_variant", "LOW", "GENE", "ENST002", "Transcript", "2/5", "", "c.100C>T", "p.Pro34Pro", "", "", "protein_coding"],
    ]

    chosen = parser._pick_csq(entries, csq_map)
    assert chosen[4] == "ENST002"  # protein_coding, not canonical pseudogene


def test_local_vcf_no_biotype_backward_compat():
    """If BIOTYPE field is missing from CSQ, fallback to original logic."""
    parser = VCFParser()

    csq_map_no_biotype = {
        "Allele": 0, "Consequence": 1, "IMPACT": 2, "SYMBOL": 3,
        "Feature": 4, "Feature_type": 5, "CANONICAL": 6, "MANE_SELECT": 7,
    }

    entries = [
        ["A", "missense_variant", "MODERATE", "GENE", "ENST001", "Transcript", "YES", "NM_001"],
        ["A", "synonymous_variant", "LOW", "GENE", "ENST002", "Transcript", "", ""],
    ]

    chosen = parser._pick_csq(entries, csq_map_no_biotype)
    assert chosen[4] == "ENST001"  # canonical fallback works


def test_selector_lncrna_penalty():
    """lncRNA biotype also receives penalty."""
    selector = TranscriptSelector()

    score_lnc, reasons = selector._score_transcript(
        {"biotype": "lncRNA", "canonical": 1, "impact": "HIGH", "protein_domains": []},
        "TEST"
    )
    # canonical (+15) + lncRNA (-15) + HIGH (+5) = 5
    # But we want it to be lower than a canonical protein_coding MODERATE
    score_pc, _ = selector._score_transcript(
        {"biotype": "protein_coding", "canonical": 1, "impact": "MODERATE", "protein_domains": []},
        "TEST"
    )
    assert score_lnc < score_pc


if __name__ == "__main__":
    print("=" * 60)
    print("Biotype-aware Transcript Selection Tests — v0.11.5")
    print("=" * 60)

    tests = [
        test_vep_api_biotype_extracted,
        test_selector_prefers_protein_coding_over_nmd,
        test_selector_nmd_pseudogene_penalty,
        test_selector_canonical_beats_high_impact,
        test_local_vcf_prefers_protein_coding_canonical,
        test_local_vcf_fallback_to_any_protein_coding,
        test_local_vcf_no_biotype_backward_compat,
        test_selector_lncrna_penalty,
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

    print(f"\nBiotype: {passed} passed, {failed} failed")
    sys.exit(0 if failed == 0 else 1)
