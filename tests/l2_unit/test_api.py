"""
L2 Unit Tests — dgra_api.py (comprehensive coverage)
API client: proxy detection, cache integration, retry logic, all query methods.

Run: pytest -m "l2 and api" tests/l2_unit/test_api.py
"""

import sys
import tempfile
import asyncio
import json
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch, PropertyMock

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


@pytest.fixture
def temp_cache():
    """Provide a temporary DGRACache instance."""
    from dgra_cache import DGRACache
    with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as f:
        db_path = Path(f.name)
    cache = DGRACache(db_path)
    yield cache
    db_path.unlink(missing_ok=True)


@pytest.fixture
def client(temp_cache):
    """Provide a DGRAAPIClient with default config and temp cache."""
    from dgra_api import DGRAAPIClient, DGRAGlobalConfig
    config = DGRAGlobalConfig()
    return DGRAAPIClient(config, temp_cache)


# =============================================================================
# Helpers
# =============================================================================

def _mock_aiohttp_response(status=200, json_data=None, text_data=None, headers=None):
    """Build a mock aiohttp response that supports async context manager."""
    mock_resp = MagicMock()
    mock_resp.status = status
    mock_resp.headers = headers or {}
    if json_data is not None:
        mock_resp.json = AsyncMock(return_value=json_data)
    if text_data is not None:
        mock_resp.text = AsyncMock(return_value=text_data)
    mock_resp.__aenter__ = AsyncMock(return_value=mock_resp)
    mock_resp.__aexit__ = AsyncMock(return_value=False)
    return mock_resp


def _mock_aiohttp_session(mock_resp):
    """Build a mock aiohttp ClientSession that returns mock_resp from request/get/post."""
    mock_session = MagicMock()
    mock_session.__aenter__ = AsyncMock(return_value=mock_session)
    mock_session.__aexit__ = AsyncMock(return_value=False)
    # request/get/post must return mock_resp directly (not a coroutine)
    # because mock_resp already has __aenter__/__aexit__ for async-with.
    mock_session.request = MagicMock(return_value=mock_resp)
    mock_session.get = MagicMock(return_value=mock_resp)
    mock_session.post = MagicMock(return_value=mock_resp)
    mock_session.close = AsyncMock(return_value=None)
    return mock_session


# =============================================================================
# Initialization
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestDGRAAPIClientInit:
    """Test DGRAAPIClient initialization."""

    def test_init_basic(self, temp_cache):
        """API-01: Basic init stores config and cache."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        config = DGRAGlobalConfig()
        client = DGRAAPIClient(config, temp_cache)
        assert client.config is config
        assert client.cache is temp_cache

    def test_init_with_proxy_route_map(self, temp_cache):
        """API-02: proxy_route_map is stored."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        config = DGRAGlobalConfig()
        route_map = {"ncbi": "http://127.0.0.1:7897"}
        client = DGRAAPIClient(config, temp_cache, proxy_route_map=route_map)
        assert client._proxy_route_map is route_map

    def test_init_fallback_config_route_map(self, temp_cache):
        """API-02b: Falls back to config._proxy_route_map when none provided."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        config = DGRAGlobalConfig()
        config._proxy_route_map = {"ncbi": "http://127.0.0.1:7897"}
        client = DGRAAPIClient(config, temp_cache)
        assert client._proxy_route_map is config._proxy_route_map


# =============================================================================
# Proxy Detection
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestProxyDetection:
    """Test proxy auto-detection logic."""

    async def test_direct_connection_works(self, temp_cache):
        """API-03: Direct connection OK → no proxy."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        client = DGRAAPIClient(DGRAGlobalConfig(), temp_cache)

        with patch.object(client, "_probe_endpoint", return_value=True):
            proxy = await client._detect_proxy()
            assert proxy is None

    async def test_all_proxies_fail(self, temp_cache):
        """API-05: All routes fail → returns None (best effort)."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        client = DGRAAPIClient(DGRAGlobalConfig(), temp_cache)

        with patch.object(client, "_probe_endpoint", return_value=False):
            proxy = await client._detect_proxy()
            assert proxy is None

    async def test_first_working_proxy_returned(self, temp_cache):
        """API-05b: Returns first working proxy from COMMON_PROXIES list."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        client = DGRAAPIClient(DGRAGlobalConfig(), temp_cache)

        async def _fake_probe(proxy, timeout):
            return proxy == "http://127.0.0.1:7890"

        with patch.object(client, "_probe_endpoint", side_effect=_fake_probe):
            proxy = await client._detect_proxy()
            assert proxy == "http://127.0.0.1:7890"

    async def test_probe_endpoint_success(self):
        """API-03b: _probe_endpoint returns True on valid NCBI response."""
        from dgra_api import DGRAAPIClient
        mock_resp = _mock_aiohttp_response(
            status=200,
            json_data={"esearchresult": {"count": "5"}}
        )
        mock_session = MagicMock()
        mock_session.get = MagicMock(return_value=mock_resp)
        mock_session.close = AsyncMock()

        with patch("aiohttp.ClientSession", return_value=mock_session):
            ok = await DGRAAPIClient._probe_endpoint(None, timeout=3.0)
            assert ok is True

    async def test_probe_endpoint_failure(self):
        """API-04: _probe_endpoint returns False on error/empty response."""
        from dgra_api import DGRAAPIClient
        mock_resp = _mock_aiohttp_response(
            status=200,
            json_data={"esearchresult": {"count": "0"}}
        )
        mock_session = MagicMock()
        mock_session.get = MagicMock(return_value=mock_resp)
        mock_session.close = AsyncMock()

        with patch("aiohttp.ClientSession", return_value=mock_session):
            ok = await DGRAAPIClient._probe_endpoint(None, timeout=3.0)
            assert ok is False

    async def test_probe_endpoint_exception(self):
        """API-04b: _probe_endpoint returns False on exception."""
        from dgra_api import DGRAAPIClient
        mock_session = MagicMock()
        mock_session.get = MagicMock(side_effect=Exception("Connection refused"))
        mock_session.close = AsyncMock()

        with patch("aiohttp.ClientSession", return_value=mock_session):
            ok = await DGRAAPIClient._probe_endpoint("http://127.0.0.1:7897", timeout=3.0)
            assert ok is False


# =============================================================================
# Rate Limiting
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestRateLimit:
    """Test _rate_limit token bucket behavior."""

    async def test_rate_limit_sleeps_when_needed(self, client):
        """API-10: Rate limit sleeps when requests are too close."""
        with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
            # First call should not sleep
            await client._rate_limit("ensembl")
            assert mock_sleep.call_count == 0

            # Immediate second call should sleep
            await client._rate_limit("ensembl")
            assert mock_sleep.call_count == 1

    async def test_rate_limit_different_apis_independent(self, client):
        """API-10b: Rate limits are per-API independent."""
        with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
            await client._rate_limit("ensembl")
            await client._rate_limit("uniprot")
            # Different APIs, so no sleep needed
            assert mock_sleep.call_count == 0


# =============================================================================
# Request With Retry (Core)
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestRequestWithRetry:
    """Test _request_with_retry: cache, retry, offline, errors."""

    async def test_cache_hit_returns_cached_data(self, client, temp_cache):
        """API-11: Cache hit returns data immediately without HTTP call."""
        temp_cache.set(
            api_name="ensembl",
            response_data={"gene": "BRCA1"},
            http_status=200,
            confidence="high",
            url="https://rest.ensembl.org/lookup/symbol/homo_sapiens/BRCA1",
            expand="1",
        )
        client._session = AsyncMock()

        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/lookup/symbol/homo_sapiens/BRCA1",
            params={"expand": "1"},
        )
        assert result["from_cache"] is True
        assert result["data"]["gene"] == "BRCA1"
        assert result["confidence"] == "high"
        # No HTTP call made
        client._session.request.assert_not_called()

    async def test_offline_mode_no_cache_returns_none(self, client):
        """API-12: Offline mode with cache miss returns None-like result."""
        client.config.offline_mode = True
        client._session = AsyncMock()

        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/lookup/symbol/homo_sapiens/BRCA1",
        )
        assert result["data"] is None
        assert result["from_cache"] is False
        assert "Offline mode" in result.get("error", "")
        client._session.request.assert_not_called()

    async def test_success_200_caches_result(self, client):
        """API-13: HTTP 200 parses JSON, caches, and returns data."""
        mock_resp = _mock_aiohttp_response(
            status=200, json_data={"id": "ENSG000001"}
        )
        mock_resp.__aenter__ = AsyncMock(return_value=mock_resp)
        mock_resp.__aexit__ = AsyncMock(return_value=False)
        client._session = _mock_aiohttp_session(mock_resp)

        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/lookup/symbol/homo_sapiens/BRCA1",
        )
        assert result["http_status"] == 200
        assert result["data"]["id"] == "ENSG000001"
        assert result["from_cache"] is False
        assert result["confidence"] == "medium"

    async def test_404_caches_negative_result(self, client):
        """API-14: HTTP 404 caches negative result with short TTL."""
        mock_resp = _mock_aiohttp_response(status=404)
        mock_resp.__aenter__ = AsyncMock(return_value=mock_resp)
        mock_resp.__aexit__ = AsyncMock(return_value=False)
        client._session = _mock_aiohttp_session(mock_resp)

        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/lookup/symbol/homo_sapiens/FAKEGENE",
        )
        assert result["http_status"] == 404
        assert result["error"] == "Not found"

    async def test_429_retries_then_success(self, client):
        """API-15: HTTP 429 triggers retry with backoff, then success."""
        resp_429 = _mock_aiohttp_response(status=429)
        resp_200 = _mock_aiohttp_response(status=200, json_data={"ok": True})

        mock_session = _mock_aiohttp_session(resp_429)
        mock_session.request = MagicMock(side_effect=[resp_429, resp_200])
        client._session = mock_session

        with patch("asyncio.sleep", new_callable=AsyncMock):
            result = await client._request_with_retry(
                api_name="ensembl",
                endpoint="/test",
            )
        assert result["http_status"] == 200
        assert result["data"]["ok"] is True
        assert client._session.request.call_count == 2

    async def test_502_retries_then_gives_up(self, client):
        """API-16: HTTP 502 retries up to max_retries then returns error."""
        resp_502 = _mock_aiohttp_response(status=502)

        mock_session = _mock_aiohttp_session(resp_502)
        client._session = mock_session

        with patch("asyncio.sleep", new_callable=AsyncMock):
            result = await client._request_with_retry(
                api_name="ensembl",
                endpoint="/test",
            )
        assert result["http_status"] == 502
        assert "Server temporarily unavailable (502)" in result.get("error", "")

    async def test_timeout_retries_then_gives_up(self, client):
        """API-17: Timeout triggers retry then returns error."""
        mock_session = _mock_aiohttp_session(_mock_aiohttp_response(status=200))
        mock_session.request = MagicMock(side_effect=asyncio.TimeoutError)
        client._session = mock_session

        with patch("asyncio.sleep", new_callable=AsyncMock):
            result = await client._request_with_retry(
                api_name="ensembl",
                endpoint="/test",
            )
        assert result["data"] is None
        assert "Timeout" in result.get("error", "")

    async def test_client_error_retries_then_gives_up(self, client):
        """API-18: aiohttp.ClientError triggers retry then returns error."""
        import aiohttp
        mock_session = _mock_aiohttp_session(_mock_aiohttp_response(status=200))
        mock_session.request = MagicMock(side_effect=aiohttp.ClientError("Connection reset"))
        client._session = mock_session

        with patch("asyncio.sleep", new_callable=AsyncMock):
            result = await client._request_with_retry(
                api_name="ensembl",
                endpoint="/test",
            )
        assert result["data"] is None
        assert "Connection error" in result.get("error", "")

    async def test_non_json_200_returns_error(self, client):
        """API-19: HTTP 200 with non-JSON response returns parse error."""
        mock_resp = _mock_aiohttp_response(status=200, text_data="not json")
        mock_resp.json = AsyncMock(side_effect=Exception("bad json"))
        client._session = _mock_aiohttp_session(mock_resp)

        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/test",
        )
        assert result["http_status"] == 200
        assert result["data"] is None
        assert "not valid JSON" in result.get("error", "")

    async def test_400_client_error_no_retry(self, client):
        """API-20: HTTP 400 (client error) does not retry."""
        mock_resp = _mock_aiohttp_response(status=400, text_data="Bad Request")
        client._session = _mock_aiohttp_session(mock_resp)

        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/test",
        )
        assert result["http_status"] == 400
        assert client._session.request.call_count == 1

    async def test_proxy_route_map_used(self, client):
        """API-21: _proxy_route_map.get_proxy is called per API."""
        mock_proxy = MagicMock()
        mock_proxy.get_proxy = MagicMock(return_value="http://proxy:8080")
        client._proxy_route_map = mock_proxy

        mock_resp = _mock_aiohttp_response(status=200, json_data={})
        client._session = _mock_aiohttp_session(mock_resp)

        await client._request_with_retry(
            api_name="ensembl",
            endpoint="/test",
        )
        mock_proxy.get_proxy.assert_called_once_with("ensembl")

    async def test_429_with_retry_after_header(self, client):
        """API-22: HTTP 429 with Retry-After header uses that value."""
        resp_429 = _mock_aiohttp_response(status=429, headers={"Retry-After": "2"})
        resp_200 = _mock_aiohttp_response(status=200, json_data={"ok": True})

        mock_session = _mock_aiohttp_session(resp_429)
        mock_session.request = MagicMock(side_effect=[resp_429, resp_200])
        client._session = mock_session

        with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
            await client._request_with_retry(
                api_name="ensembl",
                endpoint="/test",
            )
        # Should sleep for 2 seconds (Retry-After) among the calls
        calls = [c for c in mock_sleep.call_args_list if c.args == (2,)]
        assert len(calls) >= 1


# =============================================================================
# Ensembl REST API
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestEnsemblQueries:
    """Test Ensembl gene, transcript, and VEP queries."""

    async def test_query_ensembl_gene_success(self, client):
        """API-30: query_ensembl_gene parses canonical transcript."""
        raw = {
            "Transcript": [
                {"id": "ENST000001", "is_canonical": 1},
                {"id": "ENST000002", "is_canonical": 0},
            ],
            "biotype": "protein_coding",
            "description": "BRCA1 DNA repair associated",
            "seq_region_name": "17",
            "start": 43044295,
            "end": 43125364,
            "strand": -1,
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ensembl_gene("BRCA1")
            assert result["canonical_transcript"] == "ENST000001"
            assert result["biotype"] == "protein_coding"
            assert result["source"] == "ensembl"

    async def test_query_ensembl_gene_fallback_first_transcript(self, client):
        """API-31: Falls back to first transcript when none is canonical."""
        raw = {
            "Transcript": [
                {"id": "ENST000003", "is_canonical": 0},
            ],
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ensembl_gene("BRCA1")
            assert result["canonical_transcript"] == "ENST000003"

    async def test_query_ensembl_gene_failure(self, client):
        """API-32: Failed request returns low-confidence result."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 500, "from_cache": False, "confidence": "low",
            "error": "Server error"
        }):
            result = await client.query_ensembl_gene("BRCA1")
            assert result["source"] == "failed"
            assert result["confidence"] == "low"

    async def test_query_ensembl_transcript_info_success(self, client):
        """API-33: query_ensembl_transcript_info parses CDS/exons."""
        raw = {
            "id": "ENST000001",
            "display_name": "BRCA1-001",
            "biotype": "protein_coding",
            "CDS": [{}, {}, {}],
            "Exon": [{}, {}],
            "Translation": {"id": "ENSP000001"},
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ensembl_transcript_info("ENST000001")
            assert result["transcript_id"] == "ENST000001"
            assert result["cds_length"] == 3
            assert result["exon_count"] == 2
            assert result["translation_id"] == "ENSP000001"

    async def test_query_ensembl_transcript_info_failure(self, client):
        """API-34: Failed transcript query returns error result."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 404, "from_cache": False, "confidence": "low", "error": "Not found"
        }):
            result = await client.query_ensembl_transcript_info("ENST_FAKE")
            assert result["transcript_id"] == "ENST_FAKE"
            assert result["source"] == "failed"

    async def test_query_ensembl_vep_region_success(self, client):
        """API-35: query_ensembl_vep_region parses VEP response."""
        vep_data = [{
            "input": "1 100 A G . . .",
            "transcript_consequences": [
                {
                    "transcript_id": "ENST000001",
                    "gene_symbol": "TP53",
                    "consequence_terms": ["missense_variant"],
                    "impact": "MODERATE",
                    "hgvsc": "c.818C>T",
                    "hgvsp": "p.Arg273Trp",
                    "canonical": 1,
                    "protein_domains": [
                        {"name": "P53", "db": "Pfam", "start": 100, "end": 300},
                    ],
                }
            ]
        }]
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": vep_data, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ensembl_vep_region("1", 100, "A", "G")
            assert result["transcript_id"] == "ENST000001"
            assert result["gene_symbol"] == "TP53"
            assert result["impact"] == "MODERATE"
            assert len(result["protein_domains"]) == 1

    async def test_query_ensembl_vep_region_invalid_response(self, client):
        """API-36: VEP non-list response returns error."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": {"error": "bad"}, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ensembl_vep_region("1", 100, "A", "G")
            assert result["source"] == "failed"
            assert "Invalid VEP response format" in result.get("error", "")

    async def test_query_ensembl_vep_region_empty_result(self, client):
        """API-37: VEP empty list returns error result."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": [], "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ensembl_vep_region("1", 100, "A", "G")
            assert result["source"] == "failed"


# =============================================================================
# VEP Parsing
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestVEPParsing:
    """Test _parse_vep_batch_response logic."""

    def test_parse_vep_canonical_priority(self, client):
        """API-40: Priority 1 — canonical=1 transcript."""
        data = [{
            "transcript_consequences": [
                {"transcript_id": "T1", "canonical": 0, "consequence_terms": ["synonymous_variant"]},
                {"transcript_id": "T2", "canonical": 1, "consequence_terms": ["missense_variant"]},
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["transcript_id"] == "T2"

    def test_parse_vep_mane_select_priority(self, client):
        """API-41: Priority 2 — MANE Select transcript."""
        data = [{
            "transcript_consequences": [
                {"transcript_id": "T1", "canonical": 0, "mane_select": 0, "consequence_terms": ["synonymous_variant"]},
                {"transcript_id": "T2", "canonical": 0, "mane_select": 1, "consequence_terms": ["missense_variant"]},
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["transcript_id"] == "T2"

    def test_parse_vep_protein_coding_priority(self, client):
        """API-42: Priority 3 — protein_coding biotype."""
        data = [{
            "transcript_consequences": [
                {"transcript_id": "T1", "biotype": "lncRNA", "consequence_terms": ["non_coding_transcript_variant"]},
                {"transcript_id": "T2", "biotype": "protein_coding", "consequence_terms": ["missense_variant"]},
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["transcript_id"] == "T2"

    def test_parse_vep_first_fallback(self, client):
        """API-43: Priority 4 — first available transcript."""
        data = [{
            "transcript_consequences": [
                {"transcript_id": "T1", "consequence_terms": ["upstream_gene_variant"]},
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["transcript_id"] == "T1"

    def test_parse_vep_no_consequences(self, client):
        """API-44: No transcript consequences returns low-confidence entry."""
        data = [{"transcript_consequences": []}]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert "error" in results[0]

    def test_parse_vep_protein_domains_dict(self, client):
        """API-45: Protein domains as dicts are parsed."""
        data = [{
            "transcript_consequences": [
                {
                    "transcript_id": "T1",
                    "consequence_terms": ["missense_variant"],
                    "protein_domains": [
                        {"name": "DomainA", "db": "Pfam", "start": 10, "end": 100},
                    ],
                }
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["protein_domains"][0]["name"] == "DomainA"

    def test_parse_vep_protein_domains_string(self, client):
        """API-46: Protein domains as "Db:ID:Name" strings are parsed."""
        data = [{
            "transcript_consequences": [
                {
                    "transcript_id": "T1",
                    "consequence_terms": ["missense_variant"],
                    "protein_domains": ["Pfam:IPR001:DomainA"],
                }
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["protein_domains"][0]["name"] == "DomainA"
        assert results[0]["protein_domains"][0]["db"] == "Pfam"

    def test_parse_vep_protein_domains_string_two_parts(self, client):
        """API-47: Two-part domain strings are parsed."""
        data = [{
            "transcript_consequences": [
                {
                    "transcript_id": "T1",
                    "consequence_terms": ["missense_variant"],
                    "protein_domains": ["Pfam:DomainA"],
                }
            ]
        }]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert results[0]["protein_domains"][0]["name"] == "DomainA"

    def test_parse_vep_invalid_variant_result(self, client):
        """API-48: Non-dict variant result returns error entry."""
        data = ["invalid"]
        results = client._parse_vep_batch_response(data, [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
        assert "error" in results[0]


# =============================================================================
# Batch VEP
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestBatchVEP:
    """Test batch_query_vep_region."""

    async def test_batch_query_empty_list(self, client):
        """API-50: Empty variant list returns empty dict."""
        result = await client.batch_query_vep_region([])
        assert result == {}

    async def test_batch_query_single_variant(self, client):
        """API-51: Single variant batch query."""
        vep_data = [{
            "transcript_consequences": [
                {
                    "transcript_id": "ENST000001",
                    "gene_symbol": "TP53",
                    "consequence_terms": ["missense_variant"],
                    "impact": "MODERATE",
                    "canonical": 1,
                }
            ]
        }]
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": vep_data, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.batch_query_vep_region([{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}])
            assert "1:100_A>G" in result
            assert result["1:100_A>G"]["transcript_id"] == "ENST000001"

    async def test_batch_query_chunking(self, client):
        """API-52: Variants are chunked correctly (max 50 per chunk)."""
        call_count = 0
        async def _fake_request(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            body = kwargs.get("json_body", [])
            n = len(body) if isinstance(body, list) else 1
            return {
                "data": [
                    {
                        "transcript_consequences": [
                            {"transcript_id": "ENST000001", "canonical": 1, "consequence_terms": ["synonymous_variant"]}
                        ]
                    }
                    for _ in range(n)
                ],
                "http_status": 200, "from_cache": False, "confidence": "medium"
            }

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            variants = [{"chrom": "1", "pos": i, "ref": "A", "alt": "G"} for i in range(55)]
            result = await client.batch_query_vep_region(variants)
            assert call_count == 2  # 50 + 5
            assert len(result) == 55

    async def test_batch_query_chunk_exception_handled(self, client):
        """API-53: Chunk exceptions are handled gracefully."""
        async def _fake_request(*args, **kwargs):
            raise Exception("Simulated failure")

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
            result = await client.batch_query_vep_region(variants)
            assert result == {}


# =============================================================================
# UniProt
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestUniProt:
    """Test query_uniprot_by_gene."""

    async def test_query_uniprot_success(self, client):
        """API-60: query_uniprot_by_gene returns parsed protein data."""
        search_data = {
            "results": [
                {
                    "primaryAccession": "P04637",
                    "entryType": "UniProtKB reviewed (Swiss-Prot)",
                    "sequence": {"length": 393},
                }
            ]
        }
        entry_data = {
            "primaryAccession": "P04637",
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Cellular tumor antigen p53"}}
            },
            "sequence": {"length": 393},
            "features": [
                {"type": "DOMAIN", "description": "Transactivation", "location": {"start": {"value": 1}, "end": {"value": 42}}}
            ],
            "uniProtKBCrossReferences": [
                {"database": "GO", "id": "GO:0006915"},
                {"database": "InterPro", "id": "IPR002117"},
            ],
        }

        async def _fake_request(api_name, endpoint, **kwargs):
            if "search" in endpoint:
                return {"data": search_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": entry_data, "http_status": 200, "from_cache": False, "confidence": "medium"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            result = await client.query_uniprot_by_gene("TP53")
            assert result["uniprot_id"] == "P04637"
            assert result["protein_name"] == "Cellular tumor antigen p53"
            assert result["sequence_length"] == 393
            assert len(result["domains"]) == 1
            assert "GO:0006915" in result["go_terms"]
            assert "IPR002117" in result["interpro_ids"]

    async def test_query_uniprot_no_results(self, client):
        """API-61: No search results returns appropriate error."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": {"results": []}, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_uniprot_by_gene("FAKEGENE")
            assert result["uniprot_id"] is None
            assert "No UniProt entry found" in result.get("error", "")

    async def test_query_uniprot_search_fails(self, client):
        """API-62: Search failure returns low-confidence result."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "Server error"
        }):
            result = await client.query_uniprot_by_gene("TP53")
            assert result["source"] == "failed"

    async def test_query_uniprot_entry_fails(self, client):
        """API-63: Entry fetch failure returns partial result with uniprot_id."""
        search_data = {
            "results": [
                {"primaryAccession": "P04637", "entryType": "UniProtKB reviewed (Swiss-Prot)", "sequence": {"length": 393}}
            ]
        }
        async def _fake_request(api_name, endpoint, **kwargs):
            if "search" in endpoint:
                return {"data": search_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "fail"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            result = await client.query_uniprot_by_gene("TP53")
            assert result["uniprot_id"] == "P04637"
            assert result["source"] == "failed"

    async def test_query_uniprot_prefers_reviewed(self, client):
        """API-64: Prefers reviewed (Swiss-Prot) over unreviewed."""
        search_data = {
            "results": [
                {"primaryAccession": "Q9XXX1", "entryType": "UniProtKB unreviewed (TrEMBL)", "sequence": {"length": 200}},
                {"primaryAccession": "P04637", "entryType": "UniProtKB reviewed (Swiss-Prot)", "sequence": {"length": 100}},
            ]
        }
        entry_data = {
            "primaryAccession": "P04637",
            "proteinDescription": {"recommendedName": {"fullName": {"value": "p53"}}},
            "sequence": {"length": 100},
            "features": [],
            "uniProtKBCrossReferences": [],
        }
        async def _fake_request(api_name, endpoint, **kwargs):
            if "search" in endpoint:
                return {"data": search_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": entry_data, "http_status": 200, "from_cache": False, "confidence": "medium"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            result = await client.query_uniprot_by_gene("TP53")
            assert result["uniprot_id"] == "P04637"


# =============================================================================
# GTEx
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestGTEx:
    """Test GTEx expression queries."""

    async def test_query_gtex_expression_success(self, client, tmp_path):
        """API-70: query_gtex_expression resolves gencodeId and returns TPM."""
        # Avoid writing to real offline_data dir
        client._load_gencode_cache = lambda: {}
        client._save_gencode_cache = lambda c: None

        gene_data = {"data": [{"geneSymbol": "VWF", "gencodeId": "ENSG00000110799.13"}]}
        expr_data = {"data": [{"tissueSiteDetailId": "Whole_Blood", "median": 268.7}]}

        call_log = []
        async def _fake_request(api_name, endpoint, **kwargs):
            call_log.append(endpoint)
            if "reference/gene" in endpoint:
                return {"data": gene_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": expr_data, "http_status": 200, "from_cache": False, "confidence": "medium"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            result = await client.query_gtex_expression("VWF", "Whole_Blood")
            assert result["median_tpm"] == 268.7
            assert result["tissue"] == "Whole_Blood"
            assert result["gencode_id"] == "ENSG00000110799.13"

    async def test_query_gtex_expression_no_gencode(self, client):
        """API-71: Cannot resolve gencodeId returns error."""
        client._load_gencode_cache = lambda: {}
        client._save_gencode_cache = lambda c: None

        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": {"data": []}, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gtex_expression("FAKEGENE", "Whole_Blood")
            assert result["median_tpm"] is None
            assert "Could not resolve gencodeId" in result.get("error", "")

    async def test_query_gtex_expression_multi_success(self, client):
        """API-72: query_gtex_expression_multi queries multiple tissues in one call."""
        client._load_gencode_cache = lambda: {"VWF": "ENSG00000110799.13"}
        client._save_gencode_cache = lambda c: None

        expr_data = {
            "data": [
                {"tissueSiteDetailId": "Whole_Blood", "median": 268.7},
                {"tissueSiteDetailId": "Heart_Left_Ventricle", "median": 50.2},
            ]
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": expr_data, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gtex_expression_multi("VWF", ["Whole Blood", "Heart - Left Ventricle"])
            assert len(result) == 2
            assert result[0]["median_tpm"] == 268.7
            assert result[1]["median_tpm"] == 50.2

    async def test_query_gtex_expression_multi_empty_tissues(self, client):
        """API-73: Empty tissue list returns empty list."""
        result = await client.query_gtex_expression_multi("VWF", [])
        assert result == []

    async def test_query_gtex_expression_uses_cache(self, client):
        """API-74: GTEx uses cached gencodeId if available."""
        client._load_gencode_cache = lambda: {"VWF": "ENSG00000110799.13"}
        client._save_gencode_cache = lambda c: None

        expr_data = {"data": [{"tissueSiteDetailId": "Whole_Blood", "median": 268.7}]}
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": expr_data, "http_status": 200, "from_cache": False, "confidence": "medium"
        }) as mock_req:
            result = await client.query_gtex_expression("VWF", "Whole_Blood")
            assert result["median_tpm"] == 268.7
            # Should NOT call reference/gene since gencodeId is cached
            calls = [c for c in mock_req.call_args_list if "reference/gene" in str(c)]
            assert len(calls) == 0


# =============================================================================
# gnomAD
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestGnomAD:
    """Test gnomAD variant and gene constraint queries."""

    async def test_query_gnomad_variant_success(self, client):
        """API-80: query_gnomad_variant parses AF and populations."""
        raw = {
            "data": {
                "variant": {
                    "variantId": "1-12345-A-G",
                    "exome": {
                        "an": 100000, "ac": 100, "homozygote_count": 2,
                        "populations": [
                            {"id": "EAS", "ac": 10, "an": 15000, "homozygote_count": 0},
                        ]
                    },
                    "genome": {
                        "an": 50000, "ac": 50, "homozygote_count": 1,
                        "populations": [
                            {"id": "EAS", "ac": 5, "an": 7000, "homozygote_count": 0},
                        ]
                    },
                }
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gnomad_variant("1", 12345, "A", "G")
            assert result["status"] == "SUCCESS"
            assert result["af"] is not None
            assert "EAS" in result["af_populations"]
            assert result["hom_count"] == 3

    async def test_query_gnomad_variant_not_captured(self, client):
        """API-81: Variant not in gnomAD returns NOT_CAPTURED."""
        raw = {"data": {"variant": None}, "errors": []}
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gnomad_variant("1", 99999, "A", "G")
            assert result["status"] == "NOT_CAPTURED"

    async def test_query_gnomad_variant_query_error(self, client):
        """API-82: GraphQL errors return QUERY_ERROR status."""
        raw = {
            "data": {"variant": None},
            "errors": [{"message": "Field 'foo' not found"}]
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gnomad_variant("1", 99999, "A", "G")
            assert result["status"] == "QUERY_ERROR"
            assert "graphql_errors" in result

    async def test_query_gnomad_variant_api_failed(self, client):
        """API-83: HTTP failure returns API_FAILED."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "Server error"
        }):
            result = await client.query_gnomad_variant("1", 12345, "A", "G")
            assert result["status"] == "API_FAILED"

    async def test_query_gnomad_variant_strips_chr_prefix(self, client):
        """API-84: chr prefix is stripped from variant ID."""
        raw = {
            "data": {
                "variant": {
                    "variantId": "1-12345-A-G",
                    "exome": {"an": 100, "ac": 1, "homozygote_count": 0, "populations": []},
                    "genome": None,
                }
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gnomad_variant("chr1", 12345, "A", "G")
            assert result["variant_id"] == "1-12345-A-G"

    async def test_query_gnomad_gene_constraint_success(self, client):
        """API-85: query_gnomad_gene_constraint parses pLI/LOEUF."""
        raw = {
            "data": {
                "gene": {
                    "gnomad_constraint": {
                        "pLI": 1.0,
                        "oe_lof": 0.08,
                        "oe_lof_upper": 0.12,
                        "oe_lof_lower": 0.05,
                        "lof_z": 5.2,
                        "mis_z": 3.1,
                    }
                }
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gnomad_gene_constraint("BRCA1")
            assert result["status"] == "SUCCESS"
            assert result["pLI"] == 1.0
            assert result["loeuf"] == 0.12
            assert result["lof_z"] == 5.2

    async def test_query_gnomad_gene_constraint_no_data(self, client):
        """API-86: No constraint data returns NO_CONSTRAINT_DATA."""
        raw = {"data": {"gene": {"gnomad_constraint": None}}}
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_gnomad_gene_constraint("FAKEGENE")
            assert result["status"] == "NO_CONSTRAINT_DATA"

    async def test_query_gnomad_gene_constraint_api_failed(self, client):
        """API-87: HTTP failure returns API_FAILED."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "fail"
        }):
            result = await client.query_gnomad_gene_constraint("BRCA1")
            assert result["status"] == "API_FAILED"


# =============================================================================
# ClinVar Parsing
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestClinVarParsing:
    """Test static ClinVar parsers."""

    def test_parse_clinvar_efetch_json_basic(self, client):
        """API-90: _parse_clinvar_efetch_json extracts significance and traits."""
        raw = {
            "ClinVarSet": [{
                "ReferenceClinVarAssertion": {
                    "ClinicalSignificance": {
                        "Description": "Pathogenic",
                        "ReviewStatus": "reviewed by expert panel"
                    },
                    "TraitSet": {
                        "Trait": [{
                            "Name": {"ElementValue": {"$": "Breast-ovarian cancer"}}
                        }]
                    }
                },
                "ClinVarAssertion": [
                    {"ClinicalSignificance": {"Description": "Pathogenic"}},
                ]
            }]
        }
        result = client._parse_clinvar_efetch_json(raw)
        assert result["clinical_significance"] == "Pathogenic"
        assert result["review_status"] == "reviewed by expert panel"
        assert "Breast-ovarian cancer" in result["trait_names"]
        assert result["submitter_count"] == 1
        assert result["conflicting"] is False

    def test_parse_clinvar_efetch_json_conflicting(self, client):
        """API-91: Conflicting submissions are detected."""
        raw = {
            "ClinVarSet": [{
                "ReferenceClinVarAssertion": {
                    "ClinicalSignificance": {"Description": "Pathogenic", "ReviewStatus": "criteria provided, multiple submitters, no conflicts"}
                },
                "ClinVarAssertion": [
                    {"ClinicalSignificance": {"Description": "Pathogenic"}},
                    {"ClinicalSignificance": {"Description": "Benign"}},
                ]
            }]
        }
        result = client._parse_clinvar_efetch_json(raw)
        assert result["conflicting"] is True

    def test_parse_clinvar_efetch_json_no_clinvarset(self, client):
        """API-92: Missing ClinVarSet returns mostly None."""
        result = client._parse_clinvar_efetch_json({})
        assert result["clinical_significance"] is None
        assert result["submitter_count"] == 0

    def test_parse_clinvar_esummary_json_basic(self, client):
        """API-93: _parse_clinvar_esummary_json extracts germline classification."""
        raw = {
            "result": {
                "12345": {
                    "germline_classification": {
                        "description": "Likely pathogenic",
                        "review_status": "criteria provided, single submitter",
                        "trait_set": [{"trait_name": "Cardiomyopathy"}]
                    },
                    "supporting_submissions": {"scv": ["SCV001", "SCV002"]}
                }
            }
        }
        result = client._parse_clinvar_esummary_json(raw, "12345")
        assert result["clinical_significance"] == "Likely pathogenic"
        assert result["review_status"] == "criteria provided, single submitter"
        assert "Cardiomyopathy" in result["trait_names"]
        assert result["submitter_count"] == 2
        assert result["conflicting"] is False

    def test_parse_clinvar_esummary_json_missing_doc(self, client):
        """API-94: Missing doc returns empty defaults."""
        result = client._parse_clinvar_esummary_json({"result": {}}, "99999")
        assert result["clinical_significance"] is None
        assert result["submitter_count"] == 0


# =============================================================================
# NCBI ClinVar Query
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestNCBIClinVar:
    """Test query_ncbi_clinvar."""

    async def test_query_ncbi_clinvar_success(self, client):
        """API-100: query_ncbi_clinvar with position finds ClinVar record."""
        esearch_data = {"esearchresult": {"idlist": ["12345"]}}
        esummary_data = {
            "result": {
                "12345": {
                    "germline_classification": {
                        "description": "Pathogenic",
                        "review_status": "practice_guideline",
                        "trait_set": [{"trait_name": "Disease"}]
                    },
                    "supporting_submissions": {"scv": ["SCV001"]}
                }
            }
        }
        call_log = []
        async def _fake_request(api_name, endpoint, **kwargs):
            call_log.append(endpoint)
            if "esearch" in endpoint:
                return {"data": esearch_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": esummary_data, "http_status": 200, "from_cache": False, "confidence": "medium"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            result = await client.query_ncbi_clinvar("BRCA1", chrom="17", pos=43044295)
            assert result["clinvar_id"] == "12345"
            assert result["clinical_significance"] == "Pathogenic"
            assert result["source"] == "clinvar"

    async def test_query_ncbi_clinvar_no_results(self, client):
        """API-101: No ClinVar records found returns not_found note."""
        esearch_data = {"esearchresult": {"idlist": []}}
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": esearch_data, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_ncbi_clinvar("FAKEGENE")
            assert result["clinvar_id"] is None
            assert "No ClinVar records found" in result.get("note", "")

    async def test_query_ncbi_clinvar_esearch_fails(self, client):
        """API-102: ESearch failure returns failed status."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "fail"
        }):
            result = await client.query_ncbi_clinvar("BRCA1")
            assert result["source"] == "failed"

    async def test_query_ncbi_clinvar_esummary_fails(self, client):
        """API-103: ESummary failure returns partial result."""
        esearch_data = {"esearchresult": {"idlist": ["12345"]}}
        async def _fake_request(api_name, endpoint, **kwargs):
            if "esearch" in endpoint:
                return {"data": esearch_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "fail"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            result = await client.query_ncbi_clinvar("BRCA1", chrom="17", pos=43044295)
            assert result["clinvar_id"] == "12345"
            assert result["source"] == "failed"

    async def test_query_ncbi_clinvar_with_hgvs(self, client):
        """API-104: HGVS search term is used when provided."""
        esearch_data = {"esearchresult": {"idlist": ["12345"]}}
        esummary_data = {
            "result": {
                "12345": {
                    "germline_classification": {"description": "Pathogenic", "review_status": "practice_guideline", "trait_set": []},
                    "supporting_submissions": {"scv": []}
                }
            }
        }
        captured = []
        async def _fake_request(api_name, endpoint, params=None, **kwargs):
            if params:
                captured.append(params.get("term", ""))
            if "esearch" in endpoint:
                return {"data": esearch_data, "http_status": 200, "from_cache": False, "confidence": "medium"}
            return {"data": esummary_data, "http_status": 200, "from_cache": False, "confidence": "medium"}

        with patch.object(client, "_request_with_retry", side_effect=_fake_request):
            await client.query_ncbi_clinvar("BRCA1", hgvs="NM_007294.4:c.68_69del")
            # Should include the hgvs in search terms
            assert any("NM_007294" in t for t in captured)


# =============================================================================
# HGNC
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestHGNC:
    """Test HGNC symbol normalization."""

    async def test_query_hgnc_symbol_approved(self, client):
        """API-110: Approved symbol returns approved status."""
        raw = {
            "response": {
                "docs": [
                    {"hgnc_id": "HGNC:1100", "symbol": "BRCA1", "status": "Approved", "locus_type": "gene with protein product"}
                ]
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_hgnc_symbol("BRCA1")
            assert result["approved_symbol"] == "BRCA1"
            assert result["status"] == "approved"
            assert result["hgnc_id"] == "HGNC:1100"

    async def test_query_hgnc_symbol_not_found(self, client):
        """API-111: No docs returns not_found status."""
        raw = {"response": {"docs": []}}
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_hgnc_symbol("FAKEGENE")
            assert result["status"] == "not_found"

    async def test_query_hgnc_symbol_previous(self, client):
        """API-112: Input matches previous symbol → status previous."""
        raw = {
            "response": {
                "docs": [
                    {"hgnc_id": "HGNC:1100", "symbol": "BRCA1", "status": "Approved",
                     "prev_symbol": ["BRCC1"], "alias_symbol": []}
                ]
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_hgnc_symbol("BRCC1")
            assert result["status"] == "previous"

    async def test_query_hgnc_symbol_alias(self, client):
        """API-113: Input matches alias symbol → status alias."""
        raw = {
            "response": {
                "docs": [
                    {"hgnc_id": "HGNC:1100", "symbol": "BRCA1", "status": "Approved",
                     "prev_symbol": [], "alias_symbol": ["RNF53"]}
                ]
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_hgnc_symbol("RNF53")
            assert result["status"] == "alias"

    async def test_query_hgnc_symbol_withdrawn(self, client):
        """API-114: Withdrawn symbol returns withdrawn status."""
        raw = {
            "response": {
                "docs": [
                    {"hgnc_id": "HGNC:99999", "symbol": "WITHDRAWN1", "status": "Entry Withdrawn",
                     "prev_symbol": [], "alias_symbol": []}
                ]
            }
        }
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": raw, "http_status": 200, "from_cache": False, "confidence": "medium"
        }):
            result = await client.query_hgnc_symbol("WITHDRAWN1")
            assert result["status"] == "withdrawn"

    async def test_query_hgnc_symbol_api_failed(self, client):
        """API-115: API failure returns query_failed."""
        with patch.object(client, "_request_with_retry", new_callable=AsyncMock, return_value={
            "data": None, "http_status": 500, "from_cache": False, "confidence": "low", "error": "fail"
        }):
            result = await client.query_hgnc_symbol("BRCA1")
            assert result["status"] == "query_failed"
            assert result["source"] == "failed"

    async def test_batch_normalize_gene_symbols(self, client):
        """API-116: batch_normalize_gene_symbols delegates to batch_query_genes."""
        with patch.object(client, "batch_query_genes", new_callable=AsyncMock, return_value={
            "BRCA1": {"approved_symbol": "BRCA1", "status": "approved"}
        }) as mock_batch:
            result = await client.batch_normalize_gene_symbols(["BRCA1"])
            mock_batch.assert_awaited_once_with(["BRCA1"], query_type="hgnc")
            assert result["BRCA1"]["status"] == "approved"


# =============================================================================
# Batch Query Genes
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestBatchQueryGenes:
    """Test batch_query_genes."""

    async def test_batch_query_uniprot(self, client):
        """API-120: batch_query_genes with uniprot type."""
        with patch.object(client, "query_uniprot_by_gene", new_callable=AsyncMock, return_value={
            "uniprot_id": "P04637", "source": "uniprot"
        }) as mock_q:
            result = await client.batch_query_genes(["TP53", "BRCA1"], query_type="uniprot")
            assert len(result) == 2
            assert mock_q.call_count == 2

    async def test_batch_query_ensembl(self, client):
        """API-121: batch_query_genes with ensembl type."""
        with patch.object(client, "query_ensembl_gene", new_callable=AsyncMock, return_value={
            "canonical_transcript": "ENST000001", "source": "ensembl"
        }) as mock_q:
            result = await client.batch_query_genes(["TP53"], query_type="ensembl")
            assert result["TP53"]["source"] == "ensembl"

    async def test_batch_query_gtex(self, client):
        """API-122: batch_query_genes with gtex type."""
        with patch.object(client, "query_gtex_expression", new_callable=AsyncMock, return_value={
            "median_tpm": 100.0, "source": "gtex"
        }) as mock_q:
            result = await client.batch_query_genes(["VWF"], query_type="gtex", tissue="Whole Blood")
            mock_q.assert_awaited_once_with("VWF", "Whole Blood")

    async def test_batch_query_gnomad_constraint(self, client):
        """API-123: batch_query_genes with gnomad_constraint type."""
        with patch.object(client, "query_gnomad_gene_constraint", new_callable=AsyncMock, return_value={
            "pLI": 1.0, "source": "gnomad"
        }) as mock_q:
            result = await client.batch_query_genes(["BRCA1"], query_type="gnomad_constraint")
            assert result["BRCA1"]["pLI"] == 1.0

    async def test_batch_query_unknown_type(self, client):
        """API-124: Unknown query_type returns failed result per gene."""
        result = await client.batch_query_genes(["BRCA1"], query_type="unknown")
        assert result["BRCA1"]["source"] == "failed"
        assert "Unknown query_type" in result["BRCA1"]["error"]

    async def test_batch_query_exception_handled(self, client):
        """API-125: Individual gene exceptions are captured gracefully."""
        async def _fake(gene):
            if gene == "BRCA1":
                raise Exception("Simulated failure")
            return {"source": "ok"}

        with patch.object(client, "query_uniprot_by_gene", side_effect=_fake):
            result = await client.batch_query_genes(["BRCA1", "TP53"], query_type="uniprot")
            assert result["BRCA1"]["source"] == "failed"
            assert "Simulated failure" in result["BRCA1"]["error"]
            assert result["TP53"]["source"] == "ok"


# =============================================================================
# API Health
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestAPIHealth:
    """Test probe_api_health and format_health_report."""

    async def test_probe_api_health_all_ok(self, client):
        """API-130: All APIs OK when responses are 200 and fast."""
        mock_resp = _mock_aiohttp_response(status=200, json_data={})
        mock_resp.__aenter__ = AsyncMock(return_value=mock_resp)
        mock_resp.__aexit__ = AsyncMock(return_value=False)
        client._session = _mock_aiohttp_session(mock_resp)

        health = await client.probe_api_health()
        for api_name in ["ensembl", "uniprot", "gtex", "gnomad", "ncbi_eutils", "hgnc"]:
            assert health[api_name]["status"] == "OK", f"{api_name} should be OK"

    async def test_probe_api_health_slow(self, client):
        """API-131: Slow API returns SLOW status."""
        mock_resp = _mock_aiohttp_response(status=200, json_data={})
        client._session = _mock_aiohttp_session(mock_resp)

        # Make every request take a very long time
        def _time_gen():
            t = 0
            while True:
                yield t
                t += 10000

        with patch("time.time", side_effect=_time_gen()):
            health = await client.probe_api_health()
            for api_name in health:
                if api_name in client.config.apis:
                    assert health[api_name]["status"] == "SLOW"

    async def test_probe_api_health_rate_limited(self, client):
        """API-132: Rate limited API returns ERROR."""
        mock_resp = _mock_aiohttp_response(status=429)
        client._session = _mock_aiohttp_session(mock_resp)

        health = await client.probe_api_health()
        for api_name in health:
            if api_name in client.config.apis:
                assert health[api_name]["status"] == "ERROR"
                assert health[api_name]["http_status"] == 429

    async def test_probe_api_health_timeout(self, client):
        """API-133: Timeout returns ERROR with timeout message."""
        mock_session = _mock_aiohttp_session(_mock_aiohttp_response(status=200))
        mock_session.request = MagicMock(side_effect=asyncio.TimeoutError)
        client._session = mock_session

        health = await client.probe_api_health()
        for api_name in health:
            if api_name in client.config.apis:
                assert health[api_name]["status"] == "ERROR"
                assert "Timeout" in health[api_name]["error_msg"]

    async def test_probe_api_health_client_error(self, client):
        """API-134: Client error returns ERROR."""
        import aiohttp
        mock_session = _mock_aiohttp_session(_mock_aiohttp_response(status=200))
        mock_session.request = MagicMock(side_effect=aiohttp.ClientError("DNS fail"))
        client._session = mock_session

        health = await client.probe_api_health()
        for api_name in health:
            if api_name in client.config.apis:
                assert health[api_name]["status"] == "ERROR"
                assert "Connection error" in health[api_name]["error_msg"]

    async def test_probe_api_health_unexpected_error(self, client):
        """API-135: Unexpected exception returns ERROR."""
        mock_session = _mock_aiohttp_session(_mock_aiohttp_response(status=200))
        mock_session.request = MagicMock(side_effect=RuntimeError("Boom"))
        client._session = mock_session

        health = await client.probe_api_health()
        for api_name in health:
            if api_name in client.config.apis:
                assert health[api_name]["status"] == "ERROR"
                assert "Unexpected" in health[api_name]["error_msg"]

    async def test_probe_api_health_skips_unconfigured(self, client):
        """API-136: APIs not in config are skipped."""
        # Remove one API from config
        del client.config.apis["gtex"]
        mock_resp = _mock_aiohttp_response(status=200, json_data={})
        mock_resp.__aenter__ = AsyncMock(return_value=mock_resp)
        mock_resp.__aexit__ = AsyncMock(return_value=False)
        client._session = _mock_aiohttp_session(mock_resp)

        health = await client.probe_api_health()
        assert health["gtex"]["status"] == "SKIPPED"

    def test_format_health_report(self, client):
        """API-137: format_health_report produces readable output."""
        health = {
            "ensembl": {"status": "OK", "latency_ms": 120.5, "details": "Fast"},
            "gnomad": {"status": "SLOW", "latency_ms": 5000, "details": "Slow"},
            "gtex": {"status": "ERROR", "latency_ms": 0, "error_msg": "Timeout", "details": "Down"},
            "uniprot": {"status": "SKIPPED", "reason": "Not configured"},
        }
        report = client.format_health_report(health)
        assert "GPA API Health Check" in report
        assert "ensembl" in report
        assert "gnomad" in report
        assert "gtex" in report
        assert "Summary:" in report


# =============================================================================
# Cache Integration (existing)
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestCacheIntegration:
    """Test cache lookup and storage."""

    def test_cache_set_and_get(self, temp_cache):
        """API-06: Cache set/get roundtrip."""
        temp_cache.set("gnomad", {"af": 0.01}, chrom="1", pos=100)
        result = temp_cache.get("gnomad", chrom="1", pos=100)
        assert result is not None
        assert result["data"] == {"af": 0.01}
        assert result["from_cache"] is True

    def test_cache_miss(self, temp_cache):
        """API-07: Cache miss → None."""
        result = temp_cache.get("nonexistent", chrom="1", pos=100)
        assert result is None


# =============================================================================
# Offline Mode (existing)
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestOfflineMode:
    """Test offline mode behavior."""

    def test_offline_mode_flag(self, temp_cache):
        """API-08: Offline mode flag is accessible."""
        from dgra_api import DGRAAPIClient, DGRAGlobalConfig
        config = DGRAGlobalConfig()
        config.offline_mode = True
        client = DGRAAPIClient(config, temp_cache)
        assert config.offline_mode is True

    async def test_offline_mode_request_with_retry(self, client):
        """API-09: Offline mode returns no-cache error on miss."""
        client.config.offline_mode = True
        client._session = AsyncMock()
        result = await client._request_with_retry(
            api_name="ensembl",
            endpoint="/lookup/symbol/homo_sapiens/BRCA1",
        )
        assert result["data"] is None
        assert "Offline mode" in result.get("error", "")
        client._session.request.assert_not_called()


# =============================================================================
# Async Context Manager
# =============================================================================

@pytest.mark.l2
@pytest.mark.api
@pytest.mark.mock
class TestAsyncContextManager:
    """Test __aenter__ / __aexit__."""

    async def test_aenter_detects_proxy(self, client):
        """API-140: __aenter__ calls _detect_proxy and creates session."""
        with patch.object(client, "_detect_proxy", new_callable=AsyncMock, return_value="http://127.0.0.1:7897"):
            with patch("aiohttp.ClientSession") as mock_session_cls:
                mock_session = AsyncMock()
                mock_session.close = AsyncMock()
                mock_session_cls.return_value = mock_session
                async with client as c:
                    assert c is client
                    assert c._proxy_url == "http://127.0.0.1:7897"
                mock_session.close.assert_awaited_once()

    async def test_aexit_closes_session(self, client):
        """API-141: __aexit__ closes the session."""
        mock_session = AsyncMock()
        mock_session.close = AsyncMock()
        client._session = mock_session
        await client.__aexit__(None, None, None)
        mock_session.close.assert_awaited_once()
        assert client._session is None
