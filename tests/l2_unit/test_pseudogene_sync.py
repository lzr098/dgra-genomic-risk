#!/usr/bin/env python3
"""L2 unit tests for dgra_pseudogene_sync.py"""
import pytest
from unittest.mock import patch, MagicMock
from pathlib import Path
import json

import dgra_pseudogene_sync as pgs


class TestParseGtfAttributes:
    def test_basic(self):
        attr = 'gene_id "ENSG000001"; gene_name "TP53"; gene_type "protein_coding";'
        result = pgs._parse_gtf_attributes(attr)
        assert result["gene_id"] == "ENSG000001"
        assert result["gene_name"] == "TP53"
        assert result["gene_type"] == "protein_coding"

    def test_empty(self):
        assert pgs._parse_gtf_attributes("") == {}

    def test_no_space(self):
        assert pgs._parse_gtf_attributes("gene_id") == {}


class TestInferParentGene:
    def test_pattern_p_digits(self):
        assert pgs._infer_parent_gene("GUSBP1", "processed_pseudogene") == "GUSB"
        assert pgs._infer_parent_gene("CICP27", "processed_pseudogene") == "CIC"

    def test_pattern_p_only(self):
        assert pgs._infer_parent_gene("GUSBP", "processed_pseudogene") == "GUSB"

    def test_no_match(self):
        assert pgs._infer_parent_gene("SURF6P1", "processed_pseudogene") == "SURF6"

    def test_unprocessed_returns_none(self):
        assert pgs._infer_parent_gene("SOME", "unprocessed_pseudogene") is None


class TestLoadGencodePseudogenes:
    def test_missing_file_returns_empty(self, tmp_path):
        result = pgs.load_gencode_pseudogenes(tmp_path)
        assert result == {}
        # Reset global cache
        pgs._GENCODE_PSEUDOGENE_DB = None

    def test_loads_valid_file(self, tmp_path):
        data = {
            "pseudogenes": [{"gene_name": "PG1"}],
            "parent_pseudogene_pairs": {"GENE1": ["PG1"]},
        }
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        result = pgs.load_gencode_pseudogenes(tmp_path)
        assert len(result["pseudogenes"]) == 1
        # Reset global cache
        pgs._GENCODE_PSEUDOGENE_DB = None


class TestGetPseudogenesForGene:
    def test_found(self, tmp_path):
        data = {"parent_pseudogene_pairs": {"VWF": ["VWFP1", "VWFP2"]}}
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        result = pgs.get_pseudogenes_for_gene("VWF", tmp_path)
        assert result == ["VWFP1", "VWFP2"]
        pgs._GENCODE_PSEUDOGENE_DB = None

    def test_not_found(self, tmp_path):
        data = {"parent_pseudogene_pairs": {}}
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        result = pgs.get_pseudogenes_for_gene("UNKNOWN", tmp_path)
        assert result == []
        pgs._GENCODE_PSEUDOGENE_DB = None

    def test_deduplication(self, tmp_path):
        data = {"parent_pseudogene_pairs": {"VWF": ["VWFP1", "VWFP1"]}}
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        result = pgs.get_pseudogenes_for_gene("VWF", tmp_path)
        assert result == ["VWFP1"]
        pgs._GENCODE_PSEUDOGENE_DB = None


class TestIsGencodePseudogene:
    def test_true(self, tmp_path):
        data = {"pseudogenes": [{"gene_name": "PG1"}]}
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        assert pgs.is_gencode_pseudogene("PG1", tmp_path) is True
        pgs._GENCODE_PSEUDOGENE_DB = None

    def test_false(self, tmp_path):
        data = {"pseudogenes": [{"gene_name": "PG1"}]}
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        assert pgs.is_gencode_pseudogene("TP53", tmp_path) is False
        pgs._GENCODE_PSEUDOGENE_DB = None

    def test_empty_db(self, tmp_path):
        data = {"pseudogenes": []}
        (tmp_path / "gencode_pseudogenes.json").write_text(json.dumps(data))
        assert pgs.is_gencode_pseudogene("PG1", tmp_path) is False
        pgs._GENCODE_PSEUDOGENE_DB = None


class TestSyncGencodePseudogenes:
    @pytest.mark.asyncio
    async def test_ttl_not_expired(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        output = refs / "gencode_pseudogenes.json"
        output.write_text(json.dumps({"version": "48"}))
        result = await pgs.sync_gencode_pseudogenes(refs)
        assert result == output

    @pytest.mark.asyncio
    async def test_force_overrides_ttl(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        output = refs / "gencode_pseudogenes.json"
        output.write_text(json.dumps({"version": "48"}))
        # force=True should trigger download even if TTL not expired
        # But we can't actually download in tests, so we patch _download
        with patch("dgra_pseudogene_sync._download_gtf_streaming") as mock_dl:
            mock_dl.side_effect = RuntimeError("no network")
            # Should fallback to existing file
            result = await pgs.sync_gencode_pseudogenes(refs, force=True)
            assert result == output

    @pytest.mark.asyncio
    async def test_build_state_complete_skips(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        output = refs / "gencode_pseudogenes.json"
        output.write_text(json.dumps({"version": "48"}))
        with patch("dgra_build_state.is_step_complete", return_value=True):
            result = await pgs.sync_gencode_pseudogenes(refs)
            assert result == output

    @pytest.mark.asyncio
    async def test_download_failure_fallback(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        output = refs / "gencode_pseudogenes.json"
        output.write_text(json.dumps({"version": "48"}))
        with patch("dgra_pseudogene_sync._download_gtf_streaming") as mock_dl:
            mock_dl.side_effect = RuntimeError("network down")
            result = await pgs.sync_gencode_pseudogenes(refs, force=True)
            assert result == output

    @pytest.mark.asyncio
    async def test_download_failure_no_fallback_raises(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        with patch("dgra_pseudogene_sync._download_gtf_streaming") as mock_dl:
            mock_dl.side_effect = RuntimeError("network down")
            with pytest.raises(RuntimeError):
                await pgs.sync_gencode_pseudogenes(refs, force=True)


class TestSyncWrapper:
    @patch("dgra_pseudogene_sync.asyncio.run")
    def test_sync_wrapper(self, mock_run, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        mock_run.return_value = refs / "gencode_pseudogenes.json"
        result = pgs.sync_gencode_pseudogenes_sync(refs)
        assert result == refs / "gencode_pseudogenes.json"
