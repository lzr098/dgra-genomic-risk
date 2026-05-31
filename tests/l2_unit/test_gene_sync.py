#!/usr/bin/env python3
"""L2 unit tests for dgra_gene_sync.py"""
import pytest
from unittest.mock import patch, AsyncMock, MagicMock
from pathlib import Path
import json

import dgra_gene_sync as gs


class TestGeneListSynchronizerInit:
    def test_default_paths(self, tmp_path):
        refs = tmp_path / "references"
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        assert syncer.references_dir == refs
        assert syncer.offline_mode is False
        assert syncer.ttl_days == 7

    def test_custom_params(self, tmp_path):
        refs = tmp_path / "references"
        syncer = gs.GeneListSynchronizer(
            references_dir=refs,
            offline_mode=True,
            sync_enabled=False,
            ttl_days=14,
        )
        assert syncer.offline_mode is True
        assert syncer.sync_enabled is False
        assert syncer.ttl_days == 14


class TestLoadStaticLists:
    def test_missing_file_returns_empty(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        assert syncer._load_static_lists("general") == {}

    def test_loads_profile_gene_lists(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        data = {
            "profiles": {
                "general": {
                    "special_gene_lists": {
                        "coagulation": ["F8", "F9"]
                    }
                }
            }
        }
        (refs / "tissue_context.json").write_text(json.dumps(data))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        result = syncer._load_static_lists("general")
        assert result["coagulation"] == ["F8", "F9"]

    def test_unknown_profile_returns_empty(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        data = {"profiles": {"general": {"special_gene_lists": {}}}}
        (refs / "tissue_context.json").write_text(json.dumps(data))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        assert syncer._load_static_lists("neurological") == {}


class TestCacheOperations:
    def test_save_and_load_cached(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        syncer._save_cached_lists("general", {"coagulation": ["F8"]})
        cached = syncer._load_cached_lists("general")
        assert cached["coagulation"] == ["F8"]

    def test_cache_expired_when_missing(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        assert syncer._is_cache_expired("general") is True

    def test_cache_not_expired(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        syncer._save_cached_lists("general", {"coagulation": ["F8"]})
        assert syncer._is_cache_expired("general") is False


class TestApplyUserOverrides:
    def test_add_genes(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        user_cfg = {"add": {"coagulation": ["NEW1"]}}
        (refs / "user_gene_lists.json").write_text(json.dumps(user_cfg))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        merged = {"coagulation": ["F8"]}
        result = syncer._apply_user_overrides(merged)
        assert "NEW1" in result["coagulation"]
        assert "F8" in result["coagulation"]

    def test_remove_genes(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        user_cfg = {"remove": {"coagulation": ["F8"]}}
        (refs / "user_gene_lists.json").write_text(json.dumps(user_cfg))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        merged = {"coagulation": ["F8", "F9"]}
        result = syncer._apply_user_overrides(merged)
        assert "F8" not in result["coagulation"]
        assert "F9" in result["coagulation"]

    def test_custom_lists(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        user_cfg = {"custom_lists": {"my_list": {"genes": ["GENE1", "GENE2"]}}}
        (refs / "user_gene_lists.json").write_text(json.dumps(user_cfg))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        merged = {}
        result = syncer._apply_user_overrides(merged)
        assert result["my_list"] == ["GENE1", "GENE2"]

    def test_skip_underscore_keys(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        user_cfg = {"add": {"_metadata": ["X"]}, "_version": "1.0"}
        (refs / "user_gene_lists.json").write_text(json.dumps(user_cfg))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        merged = {}
        result = syncer._apply_user_overrides(merged)
        assert "_metadata" not in result


class TestValidateSymbols:
    def test_valid_symbols_preserved(self, tmp_path):
        syncer = gs.GeneListSynchronizer(references_dir=tmp_path)
        result = syncer._validate_symbols({"list": ["BRCA1", "TP53", "F8", "A1BG-AS1"]})
        assert result["list"] == ["BRCA1", "TP53", "F8", "A1BG-AS1"]

    def test_invalid_symbols_dropped(self, tmp_path):
        syncer = gs.GeneListSynchronizer(references_dir=tmp_path)
        result = syncer._validate_symbols({"list": ["BRCA1", "", "123", "A" * 30]})
        assert result["list"] == ["BRCA1"]


class TestFilterSyncedForProfile:
    def test_no_config_returns_all(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        synced = {"coagulation": ["F8"]}
        assert syncer._filter_synced_for_profile(synced, "general") == synced

    def test_filter_by_profile(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        config = {"mapping_rules": {"coagulation": ["hematopoietic"], "cardio": ["cardiovascular"]}}
        (refs / "gene_list_sources.json").write_text(json.dumps(config))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        synced = {"coagulation": ["F8"], "cardio": ["MYH7"]}
        result = syncer._filter_synced_for_profile(synced, "hematopoietic")
        assert "coagulation" in result
        assert "cardio" not in result


class TestExtractGenesFromOrphanet:
    def test_list_of_dicts(self):
        data = [{"geneSymbol": "BRCA1"}, {"symbol": "TP53"}]
        assert sorted(gs.GeneListSynchronizer._extract_genes_from_orphanet(data)) == ["BRCA1", "TP53"]

    def test_nested_gene_ref(self):
        data = [{"gene": {"symbol": "F8"}}]
        assert gs.GeneListSynchronizer._extract_genes_from_orphanet(data) == ["F8"]

    def test_single_dict(self):
        data = {"geneSymbol": "BRCA1"}
        assert gs.GeneListSynchronizer._extract_genes_from_orphanet(data) == ["BRCA1"]

    def test_nested_list(self):
        data = {"genes": [{"symbol": "F8"}, {"symbol": "F9"}]}
        assert sorted(gs.GeneListSynchronizer._extract_genes_from_orphanet(data)) == ["F8", "F9"]

    def test_empty_returns_empty(self):
        assert gs.GeneListSynchronizer._extract_genes_from_orphanet([]) == []
        assert gs.GeneListSynchronizer._extract_genes_from_orphanet({}) == []


class TestGetMergedGeneLists:
    @pytest.mark.asyncio
    async def test_offline_mode_returns_core_plus_static(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        data = {
            "profiles": {
                "general": {
                    "special_gene_lists": {
                        "custom": ["GENE1"]
                    }
                }
            }
        }
        (refs / "tissue_context.json").write_text(json.dumps(data))
        syncer = gs.GeneListSynchronizer(references_dir=refs, offline_mode=True)
        result = await syncer.get_merged_gene_lists("general")
        assert "coagulation" in result  # from CORE
        assert "custom" in result       # from static

    @pytest.mark.asyncio
    async def test_cache_hit(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs, offline_mode=True)
        syncer._save_cached_lists("general", {"cached_list": ["GENE1"]})
        result = await syncer.get_merged_gene_lists("general")
        assert "cached_list" in result

    @pytest.mark.asyncio
    async def test_force_sync_bypasses_cache(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs, offline_mode=True)
        syncer._save_cached_lists("general", {"cached_list": ["GENE1"]})
        result = await syncer.get_merged_gene_lists("general", force_sync=True)
        # Should rebuild from core + static, cache list may still be there if no static
        # but at minimum core lists should be present
        assert "coagulation" in result


class TestSyncAllSources:
    @pytest.mark.asyncio
    async def test_no_config_returns_empty(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        result = await syncer._sync_all_sources()
        assert result == {}

    @pytest.mark.asyncio
    async def test_orphanet_disabled(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        config = {"sources": {"orphanet": {"enabled": False}}}
        (refs / "gene_list_sources.json").write_text(json.dumps(config))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        result = await syncer._sync_all_sources()
        assert result == {}

    @pytest.mark.asyncio
    async def test_omim_no_api_key(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        config = {"sources": {"omim": {"enabled": True, "api_key": ""}}}
        (refs / "gene_list_sources.json").write_text(json.dumps(config))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        result = await syncer._sync_all_sources()
        assert result == {}

    @pytest.mark.asyncio
    async def test_omim_with_api_key_logs_skip(self, tmp_path):
        refs = tmp_path / "references"
        refs.mkdir()
        config = {"sources": {"omim": {"enabled": True, "api_key": "test_key"}}}
        (refs / "gene_list_sources.json").write_text(json.dumps(config))
        syncer = gs.GeneListSynchronizer(references_dir=refs)
        result = await syncer._sync_all_sources()
        # OMIM returns empty with a log note
        assert result == {}


class TestConvenienceWrapper:
    @patch("dgra_gene_sync.asyncio.run")
    def test_sync_wrapper(self, mock_run, tmp_path):
        mock_run.return_value = {"coagulation": ["F8"]}
        refs = tmp_path / "references"
        refs.mkdir()
        result = gs.get_merged_gene_lists_sync("general", references_dir=refs, offline_mode=True)
        assert result == {"coagulation": ["F8"]}
