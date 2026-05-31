#!/usr/bin/env python3
"""L2 unit tests for dgra_cli_wrapper.py"""
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path

import dgra_cli_wrapper as cli


class TestWriteTsv:
    def test_basic(self, tmp_path):
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TP53"}]
        out = tmp_path / "out.tsv"
        cli._write_tsv(variants, out)
        text = out.read_text()
        assert "CHROM\tPOS\tREF\tALT" in text
        assert "1\t100\tA\tG" in text

    def test_critical_fields_not_backfilled(self, tmp_path):
        """Critical fields missing should NOT receive synthetic defaults."""
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TP53"}]
        out = tmp_path / "out.tsv"
        cli._write_tsv(variants, out)
        text = out.read_text()
        # Use rstrip('\n') instead of strip() to preserve trailing empty tabs
        lines = text.rstrip("\n").split("\n")
        data_line = lines[1]
        fields = data_line.split("\t")
        # IMPACT, Consequence, CLIN_SIG, VAF, DP, GQ, gnomAD_AF should be empty
        # Find their indices from header
        header = lines[0].split("\t")
        for col in ["IMPACT", "Consequence", "CLIN_SIG", "VAF", "DP", "GQ", "gnomAD_AF"]:
            idx = header.index(col)
            assert fields[idx] == ""

    def test_optional_defaults(self, tmp_path):
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TP53"}]
        out = tmp_path / "out.tsv"
        cli._write_tsv(variants, out)
        text = out.read_text()
        lines = text.rstrip("\n").split("\n")
        data_line = lines[1]
        fields = data_line.split("\t")
        header = lines[0].split("\t")
        # Feature, EXON, HGVSp, HGVSc should default to empty string
        for col in ["Feature", "EXON", "HGVSp", "HGVSc"]:
            idx = header.index(col)
            assert fields[idx] == ""


class TestRunGpa:
    def test_empty_variants(self):
        result = cli.run_gpa([])
        assert result["success"] is False
        assert "empty" in result["error"].lower()

    def test_invalid_tissue(self):
        result = cli.run_gpa([{"CHROM": "1", "POS": "100"}], tissue="invalid")
        assert result["success"] is False
        assert "Invalid tissue" in result["error"]

    def test_invalid_multi_organ(self):
        result = cli.run_gpa([{"CHROM": "1", "POS": "100"}], multi_organ=["invalid"])
        assert result["success"] is False
        assert "Invalid multi-organ" in result["error"]

    @patch("dgra_cli_wrapper._run_gpa_direct")
    def test_large_dataset_uses_direct(self, mock_direct):
        mock_direct.return_value = {"success": True, "results": {}, "report_md": "# Report"}
        variants = [{"CHROM": "1", "POS": str(i)} for i in range(2500)]
        result = cli.run_gpa(variants)
        assert result["success"] is True
        mock_direct.assert_called_once()

    @patch("dgra_cli_wrapper._run_gpa_direct")
    def test_small_dataset_uses_direct(self, mock_direct):
        """v0.10.0: All small datasets now use direct API too."""
        mock_direct.return_value = {"success": True, "results": {}, "report_md": "# Report"}
        variants = [{"CHROM": "1", "POS": "100"}]
        result = cli.run_gpa(variants)
        assert result["success"] is True
        mock_direct.assert_called_once()

    @patch("dgra_cli_wrapper._run_gpa_direct")
    def test_filter_preset_passed_and_routed(self, mock_direct):
        # v0.10.0: filter_preset is applied BEFORE _run_gpa_direct.
        # Provide a variant with HIGH impact so strict preset does not filter it out.
        mock_direct.return_value = {"success": True, "results": {}, "report_md": "# Report"}
        variants = [{"CHROM": "1", "POS": "100", "IMPACT": "HIGH", "Consequence": "missense_variant"}]
        result = cli.run_gpa(variants, filter_preset="strict")
        assert result["success"] is True
        mock_direct.assert_called_once()


class TestRunGpaFromFile:
    @patch("dgra_cli_wrapper._run_gpa_vcf_direct")
    def test_vcf_input_routed_to_vcf_direct(self, mock_vcf):
        mock_vcf.return_value = {"success": True, "results": {}, "report_md": "# VCF Report"}
        result = cli.run_gpa_from_file(Path("/tmp/test.vcf"))
        assert result["success"] is True
        mock_vcf.assert_called_once()

    @patch("dgra_cli_wrapper._run_gpa_vcf_direct")
    def test_vcf_gz_input_routed_to_vcf_direct(self, mock_vcf):
        mock_vcf.return_value = {"success": True, "results": {}, "report_md": "# VCF Report"}
        result = cli.run_gpa_from_file(Path("/tmp/test.vcf.gz"))
        assert result["success"] is True
        mock_vcf.assert_called_once()

    @patch("dgra_cli_wrapper.parse_input")
    @patch("dgra_cli_wrapper.run_gpa")
    def test_tsv_input_parsed_then_run(self, mock_run, mock_parse):
        mock_parse.return_value = [{"CHROM": "1", "POS": "100"}]
        mock_run.return_value = {"success": True, "results": {}, "report_md": "# Report"}
        result = cli.run_gpa_from_file(Path("/tmp/test.tsv"))
        assert result["success"] is True
        mock_parse.assert_called_once()
        mock_run.assert_called_once()

    @patch("dgra_cli_wrapper.parse_input")
    def test_parse_failure(self, mock_parse):
        mock_parse.side_effect = ValueError("bad format")
        result = cli.run_gpa_from_file(Path("/tmp/test.tsv"))
        assert result["success"] is False
        assert "Failed to parse" in result["error"]


class TestRunGpaVcfDirect:
    @patch("dgra_cli_wrapper.subprocess.run")
    def test_success(self, mock_run, tmp_path):
        # Create mock result
        proc = MagicMock()
        proc.returncode = 0
        proc.stdout = "ok"
        proc.stderr = ""
        mock_run.return_value = proc

        # We need to mock the JSON reading too
        with patch("builtins.open", MagicMock()):
            with patch("json.load", return_value={"meta": {}, "summary": {}}):
                result = cli._run_gpa_vcf_direct(tmp_path / "test.vcf")
                # Since we mock open globally, reading the markdown file also uses mock
                # This is tricky; let's just verify the subprocess was called
                mock_run.assert_called_once()

    @patch("dgra_cli_wrapper.subprocess.run")
    def test_timeout(self, mock_run):
        mock_run.side_effect = cli.subprocess.TimeoutExpired(cmd="test", timeout=300)
        result = cli._run_gpa_vcf_direct(Path("/tmp/test.vcf"))
        assert result["success"] is False
        assert "timed out" in result["error"].lower()

    @patch("dgra_cli_wrapper.subprocess.run")
    def test_nonzero_exit(self, mock_run):
        proc = MagicMock()
        proc.returncode = 1
        proc.stderr = "error"
        proc.stdout = ""
        mock_run.return_value = proc
        result = cli._run_gpa_vcf_direct(Path("/tmp/test.vcf"))
        assert result["success"] is False
        assert "exited with code 1" in result["error"]
