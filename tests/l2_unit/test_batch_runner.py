#!/usr/bin/env python3
"""L2 unit tests for dgra_batch_runner.py"""
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

import dgra_batch_runner as br


class TestWriteTsv:
    def test_basic(self, tmp_path):
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TP53"}]
        out = tmp_path / "out.tsv"
        br._write_tsv(variants, out)
        text = out.read_text()
        assert "CHROM\tPOS\tREF\tALT" in text
        assert "1\t100\tA\tG" in text

    def test_missing_optional_defaults(self, tmp_path):
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TP53"}]
        out = tmp_path / "out.tsv"
        br._write_tsv(variants, out)
        text = out.read_text()
        lines = text.strip().split("\n")
        assert len(lines) == 2  # header + 1 row

    def test_critical_fields_empty(self, tmp_path):
        """Critical fields missing should be written as empty strings."""
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TP53"}]
        out = tmp_path / "out.tsv"
        br._write_tsv(variants, out)
        text = out.read_text()
        # IMPACT, Consequence, CLIN_SIG, VAF, DP, GQ, gnomAD_AF should be empty
        assert "\t\t" in text


class TestVariantSignature:
    def test_upper_keys(self):
        v = {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G"}
        assert br._variant_signature(v) == "1:100:A>G"

    def test_lower_keys(self):
        v = {"chrom": "X", "pos": "200", "ref": "C", "alt": "T"}
        assert br._variant_signature(v) == "X:200:C>T"

    def test_missing_defaults(self):
        v = {}
        assert br._variant_signature(v) == "::>"


class TestMergeBatchResults:
    def test_empty(self):
        merged = br.merge_batch_results([])
        assert merged["success"] is True
        assert merged["results"]["summary"]["tier1_variant_count"] == 0

    def test_single_batch(self):
        batch = {
            "success": True,
            "batch_id": 0,
            "variant_count": 1,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01", "offline_mode": False},
                "tier1_variants": [{"chrom": "1", "pos": "100", "ref": "A", "alt": "G", "gene": "TP53"}],
                "tier2_variants": [],
                "tier3_variants": [],
            },
        }
        merged = br.merge_batch_results([batch])
        assert merged["results"]["summary"]["tier1_variant_count"] == 1
        assert merged["results"]["summary"]["tier1_gene_count"] == 1

    def test_tier_priority_upgrade(self):
        """Same variant in tier2 and tier1 -> keep tier1."""
        v = {"chrom": "1", "pos": "100", "ref": "A", "alt": "G", "gene": "TP53"}
        batch1 = {
            "success": True,
            "batch_id": 0,
            "variant_count": 1,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01"},
                "tier1_variants": [],
                "tier2_variants": [v],
                "tier3_variants": [],
            },
        }
        batch2 = {
            "success": True,
            "batch_id": 1,
            "variant_count": 1,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01"},
                "tier1_variants": [v],
                "tier2_variants": [],
                "tier3_variants": [],
            },
        }
        merged = br.merge_batch_results([batch1, batch2])
        assert merged["results"]["summary"]["tier1_variant_count"] == 1
        assert merged["results"]["summary"]["tier2_variant_count"] == 0

    def test_multi_hit_detection(self):
        v1 = {"chrom": "1", "pos": "100", "ref": "A", "alt": "G", "gene": "TP53"}
        v2 = {"chrom": "1", "pos": "200", "ref": "C", "alt": "T", "gene": "TP53"}
        batch = {
            "success": True,
            "batch_id": 0,
            "variant_count": 2,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01"},
                "tier1_variants": [v1, v2],
                "tier2_variants": [],
                "tier3_variants": [],
            },
        }
        merged = br.merge_batch_results([batch])
        assert "TP53" in merged["results"]["summary"]["multi_hit_genes"]

    def test_failed_batch_ignored(self):
        batch = {
            "success": False,
            "batch_id": 0,
            "error": "timeout",
            "variant_count": 1,
        }
        merged = br.merge_batch_results([batch])
        assert merged["results"]["summary"]["tier1_variant_count"] == 0

    def test_report_md_generated(self):
        batch = {
            "success": True,
            "batch_id": 0,
            "variant_count": 1,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01"},
                "tier1_variants": [{"chrom": "1", "pos": "100", "ref": "A", "alt": "G", "gene": "TP53"}],
                "tier2_variants": [],
                "tier3_variants": [],
            },
        }
        merged = br.merge_batch_results([batch])
        assert "GPA Batch Analysis Report" in merged["report_md"]
        assert "Tier 1" in merged["report_md"]


class TestRunGpaBatched:
    def test_empty_variants(self):
        result = br.run_gpa_batched([])
        assert result["success"] is False
        assert "empty" in result["error"].lower()

    def test_invalid_tissue(self):
        result = br.run_gpa_batched([{"CHROM": "1", "POS": "100"}], tissue="invalid")
        assert result["success"] is False
        assert "Invalid tissue" in result["error"]

    @patch("dgra_batch_runner.run_batch")
    def test_single_batch_under_threshold(self, mock_run_batch):
        mock_run_batch.return_value = {
            "success": True,
            "batch_id": 0,
            "variant_count": 1,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01"},
                "tier1_variants": [],
                "tier2_variants": [],
                "tier3_variants": [],
            },
        }
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(10)]
        result = br.run_gpa_batched(variants, batch_size=500)
        assert result["success"] is True
        mock_run_batch.assert_called_once()

    @patch("dgra_batch_runner.run_batch")
    def test_multiple_batches(self, mock_run_batch):
        mock_run_batch.return_value = {
            "success": True,
            "batch_id": 1,
            "variant_count": 5,
            "elapsed_seconds": 1.0,
            "results": {
                "meta": {"analysis_date": "2024-01-01"},
                "tier1_variants": [{"chrom": "1", "pos": "1", "ref": "A", "alt": "G", "gene": "TP53"}],
                "tier2_variants": [],
                "tier3_variants": [],
            },
        }
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(10)]
        result = br.run_gpa_batched(variants, batch_size=5)
        assert result["success"] is True
        assert mock_run_batch.call_count == 2

    @patch("dgra_batch_runner.run_batch")
    def test_batch_retry_then_success(self, mock_run_batch):
        # 2 batches, batch1 retries once (2 calls), batch2 succeeds (1 call)
        mock_run_batch.side_effect = [
            {"success": False, "batch_id": 1, "error": "timeout"},
            {"success": True, "batch_id": 1, "variant_count": 5, "elapsed_seconds": 1.0,
             "results": {"meta": {}, "tier1_variants": [], "tier2_variants": [], "tier3_variants": []}},
            {"success": True, "batch_id": 2, "variant_count": 5, "elapsed_seconds": 1.0,
             "results": {"meta": {}, "tier1_variants": [], "tier2_variants": [], "tier3_variants": []}},
        ]
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(10)]
        result = br.run_gpa_batched(variants, batch_size=5, max_retries=1)
        assert result["success"] is True
        assert mock_run_batch.call_count == 3

    @patch("dgra_batch_runner.run_batch")
    def test_all_batches_fail(self, mock_run_batch):
        mock_run_batch.return_value = {"success": False, "batch_id": 1, "error": "timeout"}
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(10)]
        result = br.run_gpa_batched(variants, batch_size=5, max_retries=0)
        assert result["success"] is False
        assert "All" in result["error"] or "failed" in result["error"].lower()

    @patch("dgra_batch_runner.run_batch")
    def test_partial_failure_with_warning(self, mock_run_batch):
        calls = [
            {"success": True, "batch_id": 1, "variant_count": 5, "elapsed_seconds": 1.0,
             "results": {"meta": {}, "tier1_variants": [{"chrom": "1", "pos": "1", "ref": "A", "alt": "G", "gene": "TP53"}],
                         "tier2_variants": [], "tier3_variants": []}},
            {"success": False, "batch_id": 2, "error": "timeout"},
        ]
        mock_run_batch.side_effect = calls
        variants = [{"CHROM": "1", "POS": str(i), "REF": "A", "ALT": "G"} for i in range(10)]
        result = br.run_gpa_batched(variants, batch_size=5, max_retries=0)
        assert result["success"] is True
        assert "warning" in result
        assert "failed" in result["warning"].lower()
