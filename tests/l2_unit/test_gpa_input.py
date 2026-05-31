#!/usr/bin/env python3
"""L2 unit tests for gpa_input.py"""
import pytest
from pathlib import Path

import gpa_input as gi


class TestDetectInputType:
    def test_vcf_extension(self, tmp_path):
        f = tmp_path / "test.vcf"
        f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t100\n")
        assert gi.detect_input_type(str(f)) == gi.InputType.RAW_VCF

    def test_vcf_gz_extension(self, tmp_path):
        import gzip
        f = tmp_path / "test.vcf.gz"
        with gzip.open(f, "wt") as fh:
            fh.write("##fileformat=VCFv4.2\n#CHROM\tPOS\n1\t100\n")
        assert gi.detect_input_type(str(f)) == gi.InputType.RAW_VCF

    def test_annotated_vcf_with_csq(self, tmp_path):
        f = tmp_path / "test.vcf"
        f.write_text('##fileformat=VCFv4.2\n##INFO=<ID=CSQ,Number=.,Type=String>\n#CHROM\tPOS\n1\t100\n')
        assert gi.detect_input_type(str(f)) == gi.InputType.ANNOTATED_VCF

    def test_tsv_extension(self, tmp_path):
        f = tmp_path / "test.tsv"
        f.write_text("CHROM\tPOS\tGene\n")
        assert gi.detect_input_type(str(f)) == gi.InputType.ANNOTATED_TABLE

    def test_csv_extension(self, tmp_path):
        f = tmp_path / "test.csv"
        f.write_text("CHROM,POS,Gene\n")
        assert gi.detect_input_type(str(f)) == gi.InputType.ANNOTATED_TABLE

    def test_txt_extension(self, tmp_path):
        f = tmp_path / "test.txt"
        f.write_text("some text")
        assert gi.detect_input_type(str(f)) == gi.InputType.FREE_TEXT

    def test_content_detection_vcf(self, tmp_path):
        f = tmp_path / "unknown"
        f.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\n")
        assert gi.detect_input_type(str(f)) == gi.InputType.RAW_VCF

    def test_content_detection_table(self, tmp_path):
        f = tmp_path / "unknown"
        f.write_text("CHROM\tPOS\tGene\tConsequence\n")
        assert gi.detect_input_type(str(f)) == gi.InputType.ANNOTATED_TABLE

    def test_unknown_extension_no_content(self, tmp_path):
        f = tmp_path / "unknown"
        f.write_text("random stuff")
        assert gi.detect_input_type(str(f)) == gi.InputType.FREE_TEXT

    def test_unreadable_returns_unknown(self, tmp_path):
        f = tmp_path / "unknown"
        # Don't write anything - empty file
        f.write_text("")
        # Actually it will try to read and get empty -> free_text
        assert gi.detect_input_type(str(f)) == gi.InputType.FREE_TEXT


class TestVariantsFromVepAnnotation:
    def test_no_transcript_consequences(self):
        annotated = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "dp": 30, "gt": "0/1"}]
        result = gi.variants_from_vep_annotation(annotated)
        assert len(result) == 1
        assert result[0]["Gene"] == ""
        assert result[0]["Consequence"] == ""

    def test_single_transcript(self):
        annotated = [{
            "chrom": "1", "pos": 100, "ref": "A", "alt": "G", "dp": 30, "gt": "0/1",
            "transcript_consequences": [
                {"gene_symbol": "TP53", "transcript_id": "ENST001",
                 "consequence_terms": ["missense_variant"], "impact": "MODERATE",
                 "hgvsc": "c.100A>G", "hgvsp": "p.Ala100Val"}
            ]
        }]
        result = gi.variants_from_vep_annotation(annotated)
        assert len(result) == 1
        assert result[0]["Gene"] == "TP53"
        assert result[0]["Feature"] == "ENST001"
        assert result[0]["Consequence"] == "missense_variant"
        assert result[0]["IMPACT"] == "MODERATE"

    def test_multiple_genes(self):
        annotated = [{
            "chrom": "1", "pos": 100, "ref": "A", "alt": "G", "dp": 30, "gt": "0/1",
            "transcript_consequences": [
                {"gene_symbol": "TP53", "transcript_id": "ENST001",
                 "consequence_terms": ["missense_variant"], "impact": "MODERATE"},
                {"gene_symbol": "BRCA1", "transcript_id": "ENST002",
                 "consequence_terms": ["synonymous_variant"], "impact": "LOW"},
            ]
        }]
        result = gi.variants_from_vep_annotation(annotated)
        assert len(result) == 2
        genes = {r["Gene"] for r in result}
        assert genes == {"TP53", "BRCA1"}

    def test_with_selector(self):
        class FakeSelector:
            method = "disease_aware"
            def select(self, gene, txs):
                class Result:
                    primary = txs[0]
                    alternatives = txs[1:] if len(txs) > 1 else []
                return Result()

        annotated = [{
            "chrom": "1", "pos": 100, "ref": "A", "alt": "G", "dp": 30, "gt": "0/1",
            "transcript_consequences": [
                {"gene_symbol": "TP53", "transcript_id": "ENST001",
                 "consequence_terms": ["missense_variant"], "impact": "MODERATE",
                 "hgvsc": "c.100A>G", "hgvsp": "p.Ala100Val"}
            ]
        }]
        result = gi.variants_from_vep_annotation(annotated, selector=FakeSelector())
        assert result[0]["transcript_selection_method"] == "disease_aware"

    def test_canonical_fallback(self):
        annotated = [{
            "chrom": "1", "pos": 100, "ref": "A", "alt": "G", "dp": 30, "gt": "0/1",
            "transcript_consequences": [
                {"gene_symbol": "TP53", "transcript_id": "ENST001",
                 "consequence_terms": ["missense_variant"], "impact": "MODERATE",
                 "canonical": True},
                {"gene_symbol": "TP53", "transcript_id": "ENST002",
                 "consequence_terms": ["synonymous_variant"], "impact": "LOW"},
            ]
        }]
        result = gi.variants_from_vep_annotation(annotated)
        assert len(result) == 1  # grouped by gene
        assert result[0]["Feature"] == "ENST001"  # canonical picked
        assert result[0]["transcript_selection_method"] == "canonical"

    def test_alternative_transcripts_json(self):
        annotated = [{
            "chrom": "1", "pos": 100, "ref": "A", "alt": "G", "dp": 30, "gt": "0/1",
            "transcript_consequences": [
                {"gene_symbol": "TP53", "transcript_id": "ENST001",
                 "consequence_terms": ["missense_variant"], "impact": "MODERATE"},
                {"gene_symbol": "TP53", "transcript_id": "ENST002",
                 "consequence_terms": ["synonymous_variant"], "impact": "LOW"},
            ]
        }]
        result = gi.variants_from_vep_annotation(annotated)
        import json
        alt = json.loads(result[0]["alternative_transcripts"])
        assert len(alt) == 1
        assert alt[0]["transcript_id"] == "ENST002"
