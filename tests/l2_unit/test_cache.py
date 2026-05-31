"""L2 unit tests for dgra_cache.py — DGRACache SQLite operations."""
import json
import time
import pytest
from pathlib import Path
from dgra_cache import DGRACache
from unittest.mock import patch


@pytest.fixture
def temp_cache(tmp_path):
    db = tmp_path / "test_cache.db"
    return DGRACache(db_path=db, default_ttl_days=1)


class TestCacheSetGet:
    def test_set_and_get(self, temp_cache):
        temp_cache.set("ensembl", {"gene": "TP53"}, gene="TP53")
        result = temp_cache.get("ensembl", gene="TP53")
        assert result is not None
        assert result["data"]["gene"] == "TP53"
        assert result["from_cache"] is True

    def test_get_miss(self, temp_cache):
        result = temp_cache.get("uniprot", gene="NOT_FOUND")
        assert result is None

    def test_ttl_expiration(self, temp_cache):
        # Use negative TTL by mocking time so entry is already expired
        base_time = 1000000.0
        with patch("dgra_cache.time.time", return_value=base_time):
            temp_cache.set("ensembl", {"x": 1}, ttl_days=-1, gene="A")
        # Now advance time
        with patch("dgra_cache.time.time", return_value=base_time + 10):
            result = temp_cache.get("ensembl", gene="A")
        assert result is None

    def test_overwrite_existing(self, temp_cache):
        temp_cache.set("ensembl", {"v": 1}, gene="G")
        temp_cache.set("ensembl", {"v": 2}, gene="G")
        result = temp_cache.get("ensembl", gene="G")
        assert result["data"]["v"] == 2

    def test_short_ttl(self, temp_cache):
        base_time = 1000000.0
        with patch("dgra_cache.time.time", return_value=base_time):
            temp_cache.set_short_ttl("ensembl", {"x": 1}, ttl_minutes=-1, gene="B")
        with patch("dgra_cache.time.time", return_value=base_time + 10):
            assert temp_cache.get("ensembl", gene="B") is None


class TestCacheStats:
    def test_hit_recorded(self, temp_cache):
        temp_cache.set("ensembl", {"gene": "TP53"}, gene="TP53")
        temp_cache.get("ensembl", gene="TP53")
        stats = temp_cache.get_stats()
        assert "ensembl" in stats
        assert stats["ensembl"]["hits"] >= 1

    def test_miss_recorded(self, temp_cache):
        temp_cache.get("ensembl", gene="MISS")
        stats = temp_cache.get_stats()
        assert stats["ensembl"]["misses"] >= 1

    def test_hit_rate(self, temp_cache):
        temp_cache.set("uniprot", {"x": 1}, protein="P53")
        temp_cache.get("uniprot", protein="P53")
        stats = temp_cache.get_stats()
        assert stats["uniprot"]["hit_rate"] == 1.0


class TestCacheClear:
    def test_clear_expired(self, temp_cache):
        base_time = 1000000.0
        with patch("dgra_cache.time.time", return_value=base_time):
            temp_cache.set("a", {"x": 1}, ttl_days=-1, key="k1")
        with patch("dgra_cache.time.time", return_value=base_time + 10):
            removed = temp_cache.clear_expired()
        assert removed >= 1

    def test_clear_all(self, temp_cache):
        temp_cache.set("a", {"x": 1}, key="k1")
        temp_cache.clear_all()
        assert temp_cache.get("a", key="k1") is None
        # cache_stats table is NOT cleared by clear_all()
        assert "a" in temp_cache.get_stats()


class TestCacheInvalidate:
    def test_invalidate_specific(self, temp_cache):
        temp_cache.set("a", {"x": 1}, key="k1")
        assert temp_cache.invalidate("a", key="k1") is True
        assert temp_cache.get("a", key="k1") is None

    def test_invalidate_missing(self, temp_cache):
        assert temp_cache.invalidate("a", key="missing") is False

    def test_invalidate_pattern(self, temp_cache):
        temp_cache.set("a", {"x": 1}, key="k1")
        temp_cache.set("a", {"x": 2}, key="k2")
        removed = temp_cache.invalidate_pattern("a", pattern="%")
        assert removed == 2


class TestCacheWarm:
    def test_warm_cache(self, temp_cache):
        entries = [
            {"params": {"gene": "TP53"}, "data": {"domains": []}},
            {"params": {"gene": "BRCA1"}, "data": {"domains": []}},
        ]
        count = temp_cache.warm_cache("uniprot", entries, ttl_days=1)
        assert count == 2
        assert temp_cache.get("uniprot", gene="TP53") is not None


class TestCacheExportImport:
    def test_dump_and_load_json(self, temp_cache, tmp_path):
        # Use long TTL so load_from_json int() doesn't truncate to 0
        temp_cache.set("ensembl", {"gene": "TP53"}, ttl_days=30, gene="TP53")
        dump_path = tmp_path / "dump.json"
        temp_cache.dump_to_json(dump_path)
        assert dump_path.exists()

        # Fresh cache
        db2 = tmp_path / "cache2.db"
        cache2 = DGRACache(db_path=db2)
        loaded = cache2.load_from_json(dump_path)
        assert loaded >= 1
        assert cache2.get("ensembl", gene="TP53") is not None

    def test_get_all_for_api(self, temp_cache):
        temp_cache.set("ensembl", {"g": "A"}, gene="A")
        temp_cache.set("ensembl", {"g": "B"}, gene="B")
        all_entries = temp_cache.get_all_for_api("ensembl")
        assert len(all_entries) == 2

    def test_corrupted_entry_deleted(self, temp_cache):
        # Manually insert corrupted JSON
        cache_key = temp_cache._make_key("bad", key="k")
        with temp_cache._connection() as conn:
            conn.execute(
                "INSERT INTO api_cache (cache_key, api_name, query_params, response_data, created_at, expires_at) VALUES (?, ?, ?, ?, ?, ?)",
                (cache_key, "bad", "{}", "not-json", time.time(), time.time() + 1000),
            )
            conn.commit()
        assert temp_cache.get("bad", key="k") is None
