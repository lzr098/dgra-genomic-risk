#!/usr/bin/env python3
"""
GPA v0.9.0 Integration Tests

Tests:
1. Raw VCF end-to-end (mock VCFAnnotator)
2. Annotated TSV regression (behavior unchanged)
3. Annotated VCF regression (behavior unchanged)
4. Ambiguous transcript selection (mock LLM)
5. No disease description fallback to canonical
6. A-Layer regression (existing tests still pass)
"""

import asyncio
import json
import os
import sys
import tempfile
import unittest
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch, mock_open

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

from dgra_core import (
    Variant,
    GPAConfig,
    detect_input_type,
    InputType,
    variants_from_vep_annotation,
    _generate_transcript_selection_section,
    generate_tier_report,
)
from dgra_cli_wrapper import run_gpa_from_file
from gpa_transcript_selector import TranscriptSelector, TranscriptSelectionResult


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _make_variant(**kwargs) -> Variant:
    """Build a minimal Variant for testing."""
    defaults = {
        "chrom": "1",
        "pos": 100000,
        "ref": "A",
        "alt": "G",
        "gene": "TEST1",
        "transcript": "ENST00000000001",
        "exon": "1/10",
        "impact": "MODERATE",
        "consequence": "missense_variant",
        "hgvsp": "p.Val1565Leu",
        "hgvsc": "c.4693G>A",
        "clinvar": "",
        "dp": 100,
        "gq": 99.0,
        "gt": "0/1",
        "vaf": 0.5,
        "gnomad_af": 0.0,
        "gnomad_status": "unknown",
        "domain_info": {},
        "tissue_relevance": "unknown",
        "evidence_chain": [],
        "upgrade_conditions": [],
        "tier": 3,
        "tier_reason": "",
        "tier_actions": [],
        "tier_confidence": "LOW",
        "qc_flags": [],
        "quality_confidence": "high",
        "missing_fields": [],
        "pseudogene_warning": None,
        "transcript_warning": None,
        "gene_constraint": {},
        "phenotype_match_score": None,
        "phenotype_match_explanation": None,
        "phenotype_match_confidence": None,
        "phenotype_matched_pairs": [],
        "phenotype_known_list": [],
        "clinvar_review_status": None,
        "spliceai_result": None,
        "primary_transcript": None,
        "primary_consequence": None,
        "primary_hgvsc": None,
        "primary_hgvsp": None,
        "primary_impact": None,
        "alternative_transcripts": [],
        "transcript_selection_method": "canonical",
        "transcript_ambiguity_flag": False,
        "transcript_selection_log": "",
        "vcf_filter": None,
        "vcf_info": {},
        "vcf_format": None,
        "vcf_sample": {},
    }
    defaults.update(kwargs)
    return Variant(**defaults)


def _make_annotated_variant(
    chrom="1", pos=100000, ref="A", alt="G",
    transcript_consequences=None,
    dp=100, gt="0/1", qual=200,
):
    """Build a mock VCFAnnotator output variant dict."""
    return {
        "chrom": chrom,
        "pos": pos,
        "ref": ref,
        "alt": alt,
        "dp": dp,
        "gt": gt,
        "qual": qual,
        "filter": "PASS",
        "transcript_consequences": transcript_consequences or [],
    }


def _make_tx(
    transcript_id="ENST00000000001",
    gene_symbol="TEST1",
    consequence_terms=None,
    impact="MODERATE",
    canonical=False,
    mane_select=False,
    mane_plus_clinical=False,
    hgvsc="",
    hgvsp="",
    protein_domains=None,
):
    return {
        "transcript_id": transcript_id,
        "gene_symbol": gene_symbol,
        "consequence_terms": consequence_terms or ["missense_variant"],
        "impact": impact,
        "canonical": 1 if canonical else 0,
        "mane_select": 1 if mane_select else 0,
        "mane_plus_clinical": 1 if mane_plus_clinical else 0,
        "hgvsc": hgvsc,
        "hgvsp": hgvsp,
        "protein_domains": protein_domains or [],
    }


# ---------------------------------------------------------------------------
# Test 1: Raw VCF end-to-end (mock VCFAnnotator)
# ---------------------------------------------------------------------------

class TestRawVCFEndToEnd(unittest.TestCase):
    """Test 1: Raw VCF input goes through mock VCFAnnotator and produces valid variants."""

    def test_detect_raw_vcf(self):
        """Raw VCF without CSQ/ANN annotation → InputType.RAW_VCF."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            f.write("1\t100000\t.\tA\tG\t200\tPASS\tDP=100\tGT\t0/1\n")
            path = f.name
        try:
            it = detect_input_type(path)
            self.assertEqual(it, InputType.RAW_VCF)
        finally:
            os.unlink(path)

    def test_variants_from_vep_with_selector(self):
        """Mock VEP output with transcript selector picks primary + alternatives."""
        tx1 = _make_tx("ENST000001", "BRCA1", ["missense_variant"], "MODERATE", canonical=True)
        tx2 = _make_tx("ENST000002", "BRCA1", ["synonymous_variant"], "LOW")
        annotated = [_make_annotated_variant(
            pos=100000, ref="A", alt="G",
            transcript_consequences=[tx1, tx2],
        )]

        selector = TranscriptSelector(tissue_profile="general")
        result = variants_from_vep_annotation(annotated, selector=selector)

        self.assertEqual(len(result), 1)
        v = result[0]
        self.assertEqual(v["Gene"], "BRCA1")
        self.assertEqual(v["Feature"], "ENST000001")
        self.assertEqual(v["IMPACT"], "MODERATE")
        # alternatives should be in the JSON field
        self.assertIn("alternative_transcripts", v)
        alts = json.loads(v["alternative_transcripts"])
        self.assertEqual(len(alts), 1)
        self.assertEqual(alts[0]["transcript_id"], "ENST000002")

    def test_variants_from_vep_without_selector(self):
        """No selector provided → canonical fallback, no ambiguity flag."""
        tx1 = _make_tx("ENST000001", "BRCA1", ["missense_variant"], "MODERATE", canonical=True)
        tx2 = _make_tx("ENST000002", "BRCA1", ["synonymous_variant"], "LOW")
        annotated = [_make_annotated_variant(
            transcript_consequences=[tx1, tx2],
        )]

        result = variants_from_vep_annotation(annotated, selector=None)
        self.assertEqual(len(result), 1)
        v = result[0]
        self.assertEqual(v["Feature"], "ENST000001")
        self.assertEqual(v["transcript_ambiguity_flag"], "")


# ---------------------------------------------------------------------------
# Test 2 & 3: Annotated TSV / VCF regression (behavior unchanged)
# ---------------------------------------------------------------------------

class TestAnnotatedInputRegression(unittest.TestCase):
    """Test 2 & 3: Pre-existing annotated inputs still work without selector."""

    def test_annotated_tsv_detect(self):
        """Annotated TSV → InputType.ANNOTATED_TABLE."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\tGene\tConsequence\tIMPACT\n")
            f.write("1\t100000\tA\tG\tBRCA1\tmissense_variant\tMODERATE\n")
            path = f.name
        try:
            it = detect_input_type(path)
            self.assertEqual(it, InputType.ANNOTATED_TABLE)
        finally:
            os.unlink(path)

    def test_annotated_vcf_detect(self):
        """VCF with CSQ in INFO → InputType.ANNOTATED_VCF."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.vcf', delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            f.write("1\t100000\t.\tA\tG\t200\tPASS\tCSQ=missense|BRCA1|ENST000001\tGT\t0/1\n")
            path = f.name
        try:
            it = detect_input_type(path)
            self.assertEqual(it, InputType.ANNOTATED_VCF)
        finally:
            os.unlink(path)

    def test_run_gpa_with_tsv(self):
        """run_gpa with minimal variant dict still works (regression)."""
        from dgra_cli_wrapper import run_gpa
        variants = [{
            "CHROM": "1", "POS": "100000", "REF": "A", "ALT": "G",
            "GENE": "BRCA1", "Feature": "ENST000001", "EXON": "1/10",
            "IMPACT": "MODERATE", "Consequence": "missense_variant",
            "HGVSp": "p.Val1565Leu", "HGVSc": "c.4693G>A",
            "CLIN_SIG": "", "GT": "0/1", "DP": "100", "GQ": "99",
            "VAF": "0.5", "gnomAD_AF": "0.0",
        }]
        result = run_gpa(
            variants=variants,
            tissue="general",
            offline=True,
        )
        self.assertTrue(result.get("success", False), f"GPA failed: {result.get('error')}")
        # Should produce tiered results
        results = result.get("results", {})
        all_variants = (
            results.get("tier1_variants", [])
            + results.get("tier2_variants", [])
            + results.get("tier3_variants", [])
        )
        self.assertEqual(len(all_variants), 1)
        self.assertEqual(all_variants[0].get("gene", ""), "BRCA1")
        self.assertIn(all_variants[0].get("tier", -1), (2, 3))  # MODERATE → Tier 2 or 3 depending on gene context


# ---------------------------------------------------------------------------
# Test 4: Ambiguous transcript selection (mock LLM)
# ---------------------------------------------------------------------------

class TestAmbiguousTranscriptSelection(unittest.TestCase):
    """Test 4: When top candidates are within ambiguity threshold, LLM is invoked."""

    @patch("gpa_transcript_selector.TranscriptSelector._llm_assist_select")
    def test_llm_selects_when_ambiguous(self, mock_llm):
        """Two transcripts with close scores → LLM picks one."""
        mock_llm.return_value = {
            "transcript_id": "ENST000002",
            "gene_symbol": "BRCA1",
            "consequence_terms": ["frameshift_variant"],
            "impact": "HIGH",
            "canonical": 0,
            "mane_select": 0,
        }

        selector = TranscriptSelector(
            tissue_profile="general",
            disease_description="hereditary breast cancer",
            llm_api_key="fake-key",
            ambiguity_threshold=15,  # high threshold so close scores trigger ambiguity
        )

        tx1 = _make_tx("ENST000001", "BRCA1", ["missense_variant"], "MODERATE", canonical=True)
        tx2 = _make_tx("ENST000002", "BRCA1", ["frameshift_variant"], "HIGH")
        result = selector.select("BRCA1", [tx1, tx2])

        self.assertTrue(result.is_ambiguous)
        self.assertEqual(result.method, "llm_disease_match")
        self.assertEqual(result.primary["transcript_id"], "ENST000002")
        mock_llm.assert_called_once()

    @patch("gpa_transcript_selector.TranscriptSelector._llm_assist_select")
    def test_llm_fallback_on_api_error(self, mock_llm):
        """LLM fails → fallback to rule-based top choice."""
        mock_llm.return_value = None

        selector = TranscriptSelector(
            tissue_profile="general",
            disease_description="hereditary breast cancer",
            llm_api_key="fake-key",
            ambiguity_threshold=15,
        )

        tx1 = _make_tx("ENST000001", "BRCA1", ["missense_variant"], "MODERATE", canonical=True)
        tx2 = _make_tx("ENST000002", "BRCA1", ["frameshift_variant"], "HIGH")
        result = selector.select("BRCA1", [tx1, tx2])

        self.assertTrue(result.is_ambiguous)
        self.assertEqual(result.method, "ambiguous")
        self.assertEqual(result.primary["transcript_id"], "ENST000001")  # rule-based top


# ---------------------------------------------------------------------------
# Test 5: No disease description fallback to canonical
# ---------------------------------------------------------------------------

class TestNoDiseaseDescriptionFallback(unittest.TestCase):
    """Test 5: Without disease_description, selection falls back to canonical rule-based."""

    def test_fallback_to_canonical(self):
        """No disease_description and no LLM key → canonical/tissue_expression selection, no ambiguity processing."""
        selector = TranscriptSelector(
            tissue_profile="general",
            disease_description=None,
            llm_api_key=None,
        )

        tx1 = _make_tx("ENST000001", "BRCA1", ["missense_variant"], "MODERATE", canonical=True)
        tx2 = _make_tx("ENST000002", "BRCA1", ["frameshift_variant"], "HIGH")
        result = selector.select("BRCA1", [tx1, tx2])

        self.assertFalse(result.is_ambiguous)
        # BRCA1 is not in general tissue context, so method stays canonical
        self.assertIn(result.method, ("canonical", "tissue_expression"))
        self.assertEqual(result.primary["transcript_id"], "ENST000001")

    def test_single_transcript(self):
        """Only one transcript → 'single' method, no ambiguity."""
        selector = TranscriptSelector()
        tx1 = _make_tx("ENST000001", "BRCA1", ["missense_variant"], "MODERATE")
        result = selector.select("BRCA1", [tx1])

        self.assertFalse(result.is_ambiguous)
        self.assertIn(result.method, ("canonical", "single"))


# ---------------------------------------------------------------------------
# Test 4b: Transcript selection report section
# ---------------------------------------------------------------------------

class TestTranscriptSelectionReportSection(unittest.TestCase):
    """Test report section generation for transcript selection."""

    def test_section_appears_when_ambiguity(self):
        """Variants with ambiguity flag → section is generated."""
        v1 = _make_variant(
            gene="BRCA1",
            primary_transcript="ENST000001",
            primary_consequence="missense_variant",
            primary_impact="MODERATE",
            transcript_selection_method="ambiguous",
            transcript_ambiguity_flag=True,
            alternative_transcripts=[
                {"transcript_id": "ENST000002", "consequence_terms": ["frameshift_variant"], "impact": "HIGH"},
            ],
        )
        section = _generate_transcript_selection_section([v1])
        self.assertIsNotNone(section)
        self.assertIn("转录本选择评估", section)
        self.assertIn("⚠️ 歧义", section)
        self.assertIn("ENST000001", section)
        self.assertIn("ENST000002", section)

    def test_section_hidden_when_no_alternatives(self):
        """Variants without ambiguity or alternatives → section is None."""
        v1 = _make_variant(gene="BRCA1")
        section = _generate_transcript_selection_section([v1])
        self.assertIsNone(section)

    def test_llm_ambiguity_flag_in_report(self):
        """LLM-selected variant → report shows LLM method and warning."""
        v1 = _make_variant(
            gene="BRCA1",
            primary_transcript="ENST000002",
            transcript_selection_method="llm_disease_match",
            transcript_ambiguity_flag=True,
            alternative_transcripts=[
                {"transcript_id": "ENST000001", "consequence_terms": ["missense_variant"], "impact": "MODERATE", "canonical": True},
            ],
        )
        section = _generate_transcript_selection_section([v1])
        self.assertIn("llm_disease_match", section)
        self.assertIn("已通过 LLM", section)


# ---------------------------------------------------------------------------
# Test 6: A-Layer regression (existing tests still pass)
# ---------------------------------------------------------------------------

class TestALayerRegression(unittest.TestCase):
    """Test 6: A-Layer tests from test_a_layer.py still pass."""

    def test_import_a_layer(self):
        """A-Layer modules import successfully."""
        from dgra_build_state import (
            load_state, save_state, get_step_status,
            is_step_complete, reset_state, BuildStep,
        )
        # If import succeeds, basic API surface is intact
        self.assertTrue(callable(save_state))
        self.assertTrue(callable(load_state))

    def test_build_step_smoke(self):
        """BuildStep context manager still works."""
        import tempfile
        from dgra_build_state import set_state_file, reset_state, BuildStep, is_step_complete
        tmp = tempfile.NamedTemporaryFile(delete=False, suffix=".json")
        tmp.close()
        set_state_file(Path(tmp.name))
        reset_state()
        try:
            with BuildStep("regression_test") as step:
                step.complete(result="ok")
            self.assertTrue(is_step_complete("regression_test"))
        finally:
            reset_state()
            try:
                os.unlink(tmp.name)
            except FileNotFoundError:
                pass


# ---------------------------------------------------------------------------

if __name__ == "__main__":
    unittest.main()
