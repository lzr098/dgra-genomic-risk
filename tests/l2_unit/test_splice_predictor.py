"""
L2 Unit Tests — dgra_splice_predictor.py
SpliceAI integration and threshold classification tests.

Run: pytest -m "l2 and spliceai" tests/l2_unit/test_splice_predictor.py
"""

import asyncio
import sys
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from dgra_splice_predictor import (
    SpliceAIPredictor,
    SpliceAIResult,
    should_query_spliceai,
    reset_spliceai_cache,
    SPLICEAI_THRESHOLDS,
    _cache_key,
    _get_predictor,
    query_spliceai_batch,
)


# =============================================================================
# Helpers
# =============================================================================

def _mock_aiohttp_response(status=200, json_data=None):
    """Build a mock aiohttp response."""
    resp = AsyncMock()
    resp.status = status
    resp.headers = {}
    resp.json = AsyncMock(return_value=json_data or {})
    return resp


def _mock_aiohttp_session(response):
    """Build a mock aiohttp ClientSession context manager that returns the given response.

    The source code uses:
        async with aiohttp.ClientSession(...) as session:
            async with session.get(...) as resp:
    We need to mock both levels of context managers.
    """
    get_cm = MagicMock()
    get_cm.__aenter__ = AsyncMock(return_value=response)
    get_cm.__aexit__ = AsyncMock(return_value=False)

    session = MagicMock()
    session.get = MagicMock(return_value=get_cm)

    client_cm = MagicMock()
    client_cm.__aenter__ = AsyncMock(return_value=session)
    client_cm.__aexit__ = AsyncMock(return_value=False)
    return client_cm


@pytest.mark.l2
@pytest.mark.spliceai
class TestShouldQuerySpliceAI:
    """Tests for should_query_spliceai() function."""

    @pytest.mark.parametrize("consequence", [
        "splice_acceptor_variant",
        "splice_donor_variant",
        "splice_region_variant",
        "splice_polypyrimidine_tract_variant",
        "splice_donor_5th_base_variant",
    ])
    def test_splice_consequences_should_query(self, consequence: str):
        """SPL-01: Splice-related consequences should trigger SpliceAI query."""
        assert should_query_spliceai(consequence) is True

    @pytest.mark.parametrize("consequence", [
        "missense_variant",
        "stop_gained",
        "frameshift_variant",
        "synonymous_variant",
        "inframe_deletion",
        "inframe_insertion",
    ])
    def test_non_splice_consequences_should_not_query(self, consequence: str):
        """SPL-02: Non-splice consequences should NOT trigger SpliceAI query."""
        assert should_query_spliceai(consequence) is False

    def test_should_query_chinese_terms(self):
        """SPL-08: Chinese consequence terms should trigger query."""
        assert should_query_spliceai("剪接受体位点变异") is True
        assert should_query_spliceai("剪接供体位点变异") is True
        assert should_query_spliceai("剪接区域变异") is True

    def test_should_query_synonymous_near_boundary(self):
        """SPL-09: Synonymous variant near exon boundary triggers query."""
        assert should_query_spliceai("synonymous_variant", is_near_exon_boundary=True) is True
        assert should_query_spliceai("synonymous_variant", is_near_exon_boundary=False) is False

    def test_should_query_string_normalization(self):
        """SPL-10: String input is normalized via gpa_i18n."""
        # When normalize_consequence returns empty list
        with patch("gpa_i18n.normalize_consequence", return_value=[]):
            assert should_query_spliceai("unknown_term") is False
        # When normalize_consequence returns a splice term
        with patch("gpa_i18n.normalize_consequence", return_value=["splice_donor_variant"]):
            assert should_query_spliceai("some_term") is True


@pytest.mark.l2
@pytest.mark.spliceai
class TestSpliceAIThresholds:
    """Tests for SpliceAI threshold classification."""

    def test_canonical_none(self):
        """SPL-03: canonical delta=0 → none."""
        assert SpliceAIPredictor.determine_impact(0.0, "canonical") == "none"
        assert SpliceAIPredictor.determine_impact(0.05, "canonical") == "none"

    def test_canonical_weak(self):
        """canonical delta=0.15 → weak (>= 0.1)."""
        assert SpliceAIPredictor.determine_impact(0.15, "canonical") == "weak"

    def test_canonical_moderate(self):
        """canonical delta=0.25 → moderate (>= 0.2)."""
        assert SpliceAIPredictor.determine_impact(0.25, "canonical") == "moderate"

    def test_canonical_strong(self):
        """SPL-04: canonical delta=0.55 → strong (>= 0.5)."""
        assert SpliceAIPredictor.determine_impact(0.55, "canonical") == "strong"

    def test_splice_region_strong(self):
        """SPL-05: splice_region delta=0.25 → strong (>= 0.2)."""
        assert SpliceAIPredictor.determine_impact(0.25, "splice_region") == "strong"

    def test_splice_region_none(self):
        """splice_region delta=0.03 → none (< 0.05)."""
        assert SpliceAIPredictor.determine_impact(0.03, "splice_region") == "none"

    def test_unknown_threshold_type_fallback(self):
        """SPL-06: Unknown threshold_type falls back to splice_region thresholds."""
        assert SpliceAIPredictor.determine_impact(0.25, "unknown_type") == "strong"
        assert SpliceAIPredictor.determine_impact(0.15, "unknown_type") == "moderate"
        assert SpliceAIPredictor.determine_impact(0.03, "unknown_type") == "none"


@pytest.mark.l2
@pytest.mark.spliceai
class TestThresholdConstants:
    """Tests for threshold constant values."""

    def test_canonical_thresholds(self):
        """Canonical splice thresholds are correct."""
        canonical = SPLICEAI_THRESHOLDS["canonical"]
        assert canonical["strong"] == 0.50
        assert canonical["moderate"] == 0.20
        assert canonical["weak"] == 0.10

    def test_splice_region_thresholds(self):
        """Splice region thresholds are correct."""
        region = SPLICEAI_THRESHOLDS["splice_region"]
        assert region["strong"] == 0.20
        assert region["moderate"] == 0.10
        assert region["weak"] == 0.05


@pytest.mark.l2
@pytest.mark.spliceai
class TestCanonicalVsRegion:
    """Tests for canonical vs splice_region classification."""

    def test_is_canonical_splice_acceptor(self):
        """splice_acceptor_variant is canonical."""
        assert SpliceAIPredictor.is_canonical_splice(["splice_acceptor_variant"]) is True

    def test_is_canonical_splice_donor(self):
        """splice_donor_variant is canonical."""
        assert SpliceAIPredictor.is_canonical_splice(["splice_donor_variant"]) is True

    def test_is_not_canonical_splice_region(self):
        """splice_region_variant is NOT canonical but still queried."""
        assert SpliceAIPredictor.is_canonical_splice(["splice_region_variant"]) is False
        assert should_query_spliceai("splice_region_variant") is True

    def test_is_canonical_chinese_term(self):
        """Chinese canonical terms are recognized."""
        assert SpliceAIPredictor.is_canonical_splice(["剪接受体位点变异"]) is True
        assert SpliceAIPredictor.is_canonical_splice(["剪接供体位点变异"]) is True


@pytest.mark.l2
@pytest.mark.spliceai
class TestSpliceAIResult:
    """Tests for SpliceAIResult dataclass."""

    def test_result_creation(self):
        """SpliceAIResult can be created with all fields."""
        result = SpliceAIResult(
            delta_score=0.55,
            delta_acceptor_gain=0.0,
            delta_acceptor_loss=0.55,
            delta_donor_gain=0.0,
            delta_donor_loss=0.0,
            predicted_impact="strong",
            threshold_type="canonical",
            source="spliceai_lookup",
            raw_response={"delta_scores": {"AG": 0.0, "AL": 0.55}},
        )
        assert result.delta_score == 0.55
        assert result.predicted_impact == "strong"
        assert result.threshold_type == "canonical"
        assert result.source == "spliceai_lookup"

    def test_result_defaults(self):
        """SpliceAIResult defaults are correct."""
        result = SpliceAIResult()
        assert result.delta_score == 0.0
        assert result.predicted_impact == "none"
        assert result.threshold_type == "canonical"
        assert result.source == "unknown"


@pytest.mark.l2
@pytest.mark.spliceai
class TestCache:
    """Tests for SpliceAI cache behavior."""

    def test_reset_cache(self):
        """SPL-07: Cache can be reset."""
        reset_spliceai_cache()

    def test_cache_hit_after_insert(self):
        """Cached result is returned on subsequent queries."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor(max_concurrency=5, timeout=30)

        cache_key = predictor._cache_key("1", 100000, "A", "G")
        result = SpliceAIResult(
            delta_score=0.55,
            delta_acceptor_gain=0.0,
            delta_acceptor_loss=0.55,
            delta_donor_gain=0.0,
            delta_donor_loss=0.0,
            predicted_impact="strong",
            threshold_type="canonical",
            source="spliceai_lookup",
            raw_response={},
        )
        predictor._cache[cache_key] = result

        cached = predictor._cache.get(cache_key)
        assert cached is not None
        assert cached.delta_score == 0.55


# =============================================================================
# NEW: SpliceAIPredictor.query
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestSpliceAIPredictorQuery:
    """SPL-11~14: SpliceAIPredictor.query cache and network behavior."""

    @pytest.mark.asyncio
    async def test_query_cache_hit(self):
        """SPL-11: Cache hit returns cached result without network call."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        cache_key = predictor._cache_key("1", 100, "A", "G")
        cached = SpliceAIResult(delta_score=0.5, source="spliceai_lookup")
        predictor._cache[cache_key] = cached

        with patch("aiohttp.ClientSession") as mock_session:
            result = await predictor.query("1", 100, "A", "G")
        assert result == cached
        mock_session.assert_not_called()

    @pytest.mark.asyncio
    async def test_query_cache_miss_broad_success(self):
        """SPL-12: Cache miss → Broad Institute API success."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=200, json_data={
            "scores": [{"DS_AG": "0.1", "DS_AL": "0.2", "DS_DG": "0.0", "DS_DL": "0.0"}],
            "source": "lookup",
        })
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor.query("1", 100, "A", "G")
        assert result.delta_score == 0.2
        assert result.source == "spliceai_lookup"
        assert result.delta_acceptor_loss == 0.2

    @pytest.mark.asyncio
    async def test_query_chr_prefix_stripped(self):
        """SPL-13: chr prefix is stripped from chromosome."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=200, json_data={
            "scores": [{"DS_AG": "0.0", "DS_AL": "0.0", "DS_DG": "0.0", "DS_DL": "0.0"}],
        })

        get_cm = MagicMock()
        get_cm.__aenter__ = AsyncMock(return_value=resp)
        get_cm.__aexit__ = AsyncMock(return_value=False)
        inner_session = MagicMock()
        inner_session.get = MagicMock(return_value=get_cm)
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            result = await predictor.query("chr1", 100, "A", "G")
        assert result.delta_score == 0.0
        # Verify the URL params included variant without chr prefix
        call_args = inner_session.get.call_args
        assert "1-100-A-G" in call_args[1]["params"]["variant"]

    @pytest.mark.asyncio
    async def test_query_grch37_url(self):
        """SPL-14: GRCh37 uses the correct base URL."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=200, json_data={"scores": []})

        get_cm = MagicMock()
        get_cm.__aenter__ = AsyncMock(return_value=resp)
        get_cm.__aexit__ = AsyncMock(return_value=False)
        inner_session = MagicMock()
        inner_session.get = MagicMock(return_value=get_cm)
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            await predictor.query("1", 100, "A", "G", genome="GRCh37")
        call_args = inner_session.get.call_args
        assert call_args[1]["params"]["hg"] == "37"


# =============================================================================
# NEW: _query_with_retry
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestQueryWithRetry:
    """SPL-15~23: _query_with_retry retry and fallback logic."""

    @pytest.mark.asyncio
    async def test_broad_404_not_in_db(self):
        """SPL-15: Broad 404 → not_in_db result."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=404, json_data={})
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor.query("1", 100, "A", "G")
        assert result.source == "not_in_db"

    @pytest.mark.asyncio
    async def test_broad_429_retry_then_success(self):
        """SPL-16: 429 retries then succeeds."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp_429 = _mock_aiohttp_response(status=429, json_data={})
        resp_429.headers = {"Retry-After": "1"}
        resp_ok = _mock_aiohttp_response(status=200, json_data={"scores": [{"DS_AG": "0.5"}]})

        def _make_get_cm(resp):
            cm = MagicMock()
            cm.__aenter__ = AsyncMock(return_value=resp)
            cm.__aexit__ = AsyncMock(return_value=False)
            return cm

        inner_session = MagicMock()
        inner_session.get = MagicMock(side_effect=[
            _make_get_cm(resp_429),
            _make_get_cm(resp_ok),
        ])
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await predictor.query("1", 100, "A", "G")
        assert result.delta_score == 0.5
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_broad_timeout_then_success(self):
        """SPL-17: Timeout retries then succeeds."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp_ok = _mock_aiohttp_response(status=200, json_data={"scores": [{"DS_AG": "0.3"}]})

        def _make_get_cm(resp):
            cm = MagicMock()
            cm.__aenter__ = AsyncMock(return_value=resp)
            cm.__aexit__ = AsyncMock(return_value=False)
            return cm

        inner_session = MagicMock()
        inner_session.get = MagicMock(side_effect=[
            asyncio.TimeoutError(),
            _make_get_cm(resp_ok),
        ])
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await predictor.query("1", 100, "A", "G")
        assert result.delta_score == 0.3
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_broad_client_error_then_success(self):
        """SPL-18: ClientError retries then succeeds."""
        import aiohttp
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp_ok = _mock_aiohttp_response(status=200, json_data={"scores": [{"DS_AG": "0.4"}]})

        def _make_get_cm(resp):
            cm = MagicMock()
            cm.__aenter__ = AsyncMock(return_value=resp)
            cm.__aexit__ = AsyncMock(return_value=False)
            return cm

        inner_session = MagicMock()
        inner_session.get = MagicMock(side_effect=[
            aiohttp.ClientError("connection reset"),
            _make_get_cm(resp_ok),
        ])
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await predictor.query("1", 100, "A", "G")
        assert result.delta_score == 0.4
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_broad_all_fail_vep_fallback(self):
        """SPL-19: All Broad retries fail → VEP REST fallback."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor(vep_enabled=True)
        resp_500 = _mock_aiohttp_response(status=500, json_data={})
        broad_session = _mock_aiohttp_session(resp_500)
        vep_resp = _mock_aiohttp_response(status=200, json_data=[{
            "transcript_consequences": [
                {"transcript_id": "ENST001", "spliceai": {"DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.8, "DS_DL": 0.0}}
            ]
        }])
        vep_session = _mock_aiohttp_session(vep_resp)

        # Broad loop calls ClientSession 4 times (once per attempt), then VEP calls it once.
        call_count = 0
        def session_factory(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            if call_count <= 4:
                return broad_session
            return vep_session

        with patch("aiohttp.ClientSession", side_effect=session_factory):
            with patch("asyncio.sleep", new_callable=AsyncMock):
                result = await predictor.query("1", 100, "A", "G")
        assert result.source == "vep_rest"
        assert result.delta_score == 0.8

    @pytest.mark.asyncio
    async def test_broad_all_fail_vep_disabled(self):
        """SPL-20: All retries fail and vep_enabled=False → api_error."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor(vep_enabled=False)
        resp_500 = _mock_aiohttp_response(status=500, json_data={})
        session = _mock_aiohttp_session(resp_500)

        with patch("aiohttp.ClientSession", return_value=session):
            with patch("asyncio.sleep", new_callable=AsyncMock):
                result = await predictor.query("1", 100, "A", "G")
        assert result.source == "api_error"

    @pytest.mark.asyncio
    async def test_broad_empty_scores(self):
        """SPL-21: Broad returns empty scores list → delta_score 0."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=200, json_data={"scores": [], "source": "lookup"})
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor.query("1", 100, "A", "G")
        assert result.delta_score == 0.0
        assert result.source == "spliceai_lookup"

    @pytest.mark.asyncio
    async def test_broad_multiple_scores_max_delta(self):
        """SPL-22: Multiple score entries → max delta is selected."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=200, json_data={
            "scores": [
                {"DS_AG": "0.1", "DS_AL": "0.0", "DS_DG": "0.0", "DS_DL": "0.0"},
                {"DS_AG": "0.0", "DS_AL": "0.5", "DS_DG": "0.0", "DS_DL": "0.0"},
            ],
        })
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor.query("1", 100, "A", "G")
        assert result.delta_score == 0.5
        assert result.delta_acceptor_loss == 0.5

    @pytest.mark.asyncio
    async def test_broad_429_all_retries_exhausted(self):
        """SPL-23: 429 on all retries → fallback or api_error."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor(vep_enabled=False)
        resp_429 = _mock_aiohttp_response(status=429, json_data={})
        resp_429.headers = {"Retry-After": "1"}
        session = _mock_aiohttp_session(resp_429)

        with patch("aiohttp.ClientSession", return_value=session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await predictor.query("1", 100, "A", "G")
        assert result.source == "api_error"
        assert mock_sleep.await_count >= 3


# =============================================================================
# NEW: _query_vep_rest
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestQueryVEPRest:
    """SPL-24~28: _query_vep_rest fallback method."""

    @pytest.mark.asyncio
    async def test_vep_rest_success(self):
        """SPL-24: VEP REST fallback succeeds."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=200, json_data=[{
            "transcript_consequences": [
                {"transcript_id": "ENST001", "spliceai": {"DS_AG": 0.0, "DS_AL": 0.3, "DS_DG": 0.0, "DS_DL": 0.0}}
            ]
        }])
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor._query_vep_rest("1", 100, "A", "G")
        assert result.source == "vep_rest"
        assert result.delta_score == 0.3

    @pytest.mark.asyncio
    async def test_vep_rest_400_ref_mismatch(self):
        """SPL-25: VEP REST 400 → api_error with REF mismatch."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=400, json_data={})
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor._query_vep_rest("1", 100, "A", "G")
        assert result.source == "api_error"
        assert result.raw_response.get("reason") == "REF mismatch"

    @pytest.mark.asyncio
    async def test_vep_rest_other_http_error(self):
        """SPL-26: VEP REST 500 → api_error with status."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        resp = _mock_aiohttp_response(status=500, json_data={})
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            result = await predictor._query_vep_rest("1", 100, "A", "G")
        assert result.source == "api_error"
        assert result.raw_response.get("status") == 500

    @pytest.mark.asyncio
    async def test_vep_rest_timeout(self):
        """SPL-27: VEP REST timeout → api_error."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        inner_session = MagicMock()
        inner_session.get = MagicMock(side_effect=asyncio.TimeoutError())
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            result = await predictor._query_vep_rest("1", 100, "A", "G")
        assert result.source == "api_error"

    @pytest.mark.asyncio
    async def test_vep_rest_client_error(self):
        """SPL-28: VEP REST ClientError → api_error."""
        import aiohttp
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        inner_session = MagicMock()
        inner_session.get = MagicMock(side_effect=aiohttp.ClientError("dns failed"))
        client_cm = MagicMock()
        client_cm.__aenter__ = AsyncMock(return_value=inner_session)
        client_cm.__aexit__ = AsyncMock(return_value=False)

        with patch("aiohttp.ClientSession", return_value=client_cm):
            result = await predictor._query_vep_rest("1", 100, "A", "G")
        assert result.source == "api_error"


# =============================================================================
# NEW: _parse_response (Broad API)
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestParseResponse:
    """SPL-29~33: _parse_response Broad Institute API format."""

    def test_parse_new_api_format(self):
        """SPL-29: New API format with scores list."""
        predictor = SpliceAIPredictor()
        data = {
            "scores": [
                {"DS_AG": "0.1", "DS_AL": "0.2", "DS_DG": "0.3", "DS_DL": "0.4"},
            ],
            "source": "lookup",
        }
        result = predictor._parse_response(data, "spliceai_lookup")
        assert result.delta_score == 0.4
        assert result.delta_acceptor_gain == 0.1
        assert result.delta_acceptor_loss == 0.2
        assert result.delta_donor_gain == 0.3
        assert result.delta_donor_loss == 0.4
        assert result.source == "spliceai_lookup"

    def test_parse_empty_scores(self):
        """SPL-30: Empty scores → delta_score 0."""
        predictor = SpliceAIPredictor()
        data = {"scores": [], "source": "lookup"}
        result = predictor._parse_response(data, "spliceai_lookup")
        assert result.delta_score == 0.0
        assert result.source == "spliceai_lookup"

    def test_parse_none_scores(self):
        """SPL-31: Missing scores key → delta_score 0."""
        predictor = SpliceAIPredictor()
        data = {"source": "lookup"}
        result = predictor._parse_response(data, "spliceai_lookup")
        assert result.delta_score == 0.0

    def test_parse_exception_handling(self):
        """SPL-32: Exception during parsing → api_error."""
        predictor = SpliceAIPredictor()
        data = {"scores": "not_a_list"}
        result = predictor._parse_response(data, "spliceai_lookup")
        assert result.source == "api_error"

    def test_parse_string_scores_with_none(self):
        """SPL-33: Score values that are None or empty strings."""
        predictor = SpliceAIPredictor()
        data = {
            "scores": [
                {"DS_AG": None, "DS_AL": "", "DS_DG": "0.5", "DS_DL": None},
            ],
        }
        result = predictor._parse_response(data, "spliceai_lookup")
        assert result.delta_score == 0.5
        assert result.delta_donor_gain == 0.5


# =============================================================================
# NEW: _parse_vep_response
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestParseVEPResponse:
    """SPL-34~39: _parse_vep_response VEP REST format."""

    def test_parse_vep_with_data(self):
        """SPL-34: VEP response with spliceai data."""
        predictor = SpliceAIPredictor()
        data = {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST001",
                    "spliceai": {
                        "DS_AG": 0.1, "DS_AL": 0.2, "DS_DG": 0.0, "DS_DL": 0.0,
                        "SYMBOL": "TP53",
                    },
                },
                {
                    "transcript_id": "ENST002",
                    "spliceai": {
                        "DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.5, "DS_DL": 0.0,
                        "SYMBOL": "TP53",
                    },
                },
            ]
        }
        result = predictor._parse_vep_response(data)
        assert result.delta_score == 0.5
        assert result.delta_donor_gain == 0.5
        assert result.source == "vep_rest"
        assert result.raw_response["symbol"] == "TP53"
        assert result.raw_response["transcript_id"] == "ENST002"

    def test_parse_vep_empty_transcripts(self):
        """SPL-35: VEP response with no transcript consequences."""
        predictor = SpliceAIPredictor()
        data = {"transcript_consequences": []}
        result = predictor._parse_vep_response(data)
        assert result.source == "not_in_db"
        assert result.raw_response == data

    def test_parse_vep_no_spliceai_in_transcripts(self):
        """SPL-36: Transcripts exist but none have spliceai → delta 0."""
        predictor = SpliceAIPredictor()
        data = {
            "transcript_consequences": [
                {"transcript_id": "ENST001"},
                {"transcript_id": "ENST002"},
            ]
        }
        result = predictor._parse_vep_response(data)
        assert result.delta_score == 0.0
        assert result.source == "vep_rest"

    def test_parse_vep_all_ds_zero(self):
        """SPL-37: All transcripts have DS=0 → delta 0 with source vep_rest."""
        predictor = SpliceAIPredictor()
        data = {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST001",
                    "spliceai": {"DS_AG": 0.0, "DS_AL": 0.0, "DS_DG": 0.0, "DS_DL": 0.0},
                },
            ]
        }
        result = predictor._parse_vep_response(data)
        assert result.delta_score == 0.0
        assert result.source == "vep_rest"

    def test_parse_vep_exception_handling(self):
        """SPL-38: Exception during VEP parsing → api_error."""
        predictor = SpliceAIPredictor()
        data = {"transcript_consequences": "not_a_list"}
        result = predictor._parse_vep_response(data)
        assert result.source == "api_error"

    def test_parse_vep_float_conversion_with_strings(self):
        """SPL-39: String scores in VEP response are converted to float."""
        predictor = SpliceAIPredictor()
        data = {
            "transcript_consequences": [
                {
                    "transcript_id": "ENST001",
                    "spliceai": {"DS_AG": "0.25", "DS_AL": "0.0", "DS_DG": "0.0", "DS_DL": "0.0"},
                },
            ]
        }
        result = predictor._parse_vep_response(data)
        assert result.delta_score == 0.25


# =============================================================================
# NEW: batch_query
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestBatchQuery:
    """SPL-40~43: batch_query behavior."""

    @pytest.mark.asyncio
    async def test_batch_query_all_cache_hits(self):
        """SPL-40: All variants in cache → no network calls."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        v = MagicMock()
        v.chrom = "1"
        v.pos = 100
        v.ref = "A"
        v.alt = "G"
        predictor._cache["1:100:A:G"] = SpliceAIResult(delta_score=0.5)

        with patch.object(predictor, "query") as mock_query:
            results = await predictor.batch_query([v])
        mock_query.assert_not_called()
        assert results["1:100:A:G"].delta_score == 0.5

    @pytest.mark.asyncio
    async def test_batch_query_mixed_cache(self):
        """SPL-41: Mix of cache hits and misses."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        v1 = MagicMock(chrom="1", pos=100, ref="A", alt="G")
        v2 = MagicMock(chrom="1", pos=101, ref="C", alt="T")
        predictor._cache["1:100:A:G"] = SpliceAIResult(delta_score=0.5)

        with patch.object(predictor, "query", return_value=SpliceAIResult(delta_score=0.3)) as mock_query:
            results = await predictor.batch_query([v1, v2])
        assert mock_query.await_count == 1
        assert results["1:100:A:G"].delta_score == 0.5
        assert results["1:101:C:T"].delta_score == 0.3

    @pytest.mark.asyncio
    async def test_batch_query_exception_in_query(self):
        """SPL-42: Exception in single query → api_error for that variant."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        v1 = MagicMock(chrom="1", pos=100, ref="A", alt="G")
        v2 = MagicMock(chrom="1", pos=101, ref="C", alt="T")

        async def side_effect(*args, **kwargs):
            if args[1] == 100:
                raise RuntimeError("boom")
            return SpliceAIResult(delta_score=0.3)

        with patch.object(predictor, "query", side_effect=side_effect):
            results = await predictor.batch_query([v1, v2])
        assert results["1:100:A:G"].source == "api_error"
        assert results["1:101:C:T"].delta_score == 0.3

    @pytest.mark.asyncio
    async def test_batch_query_empty_input(self):
        """SPL-43: Empty input returns current cache dict."""
        reset_spliceai_cache()
        predictor = SpliceAIPredictor()
        results = await predictor.batch_query([])
        assert results == {}


# =============================================================================
# NEW: Module-level helpers
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestModuleLevelHelpers:
    """SPL-44~48: Module-level functions and global predictor."""

    def test_cache_key_strips_chr(self):
        """SPL-44: _cache_key strips chr prefix."""
        assert _cache_key("chr1", 100, "A", "G") == "1:100:A:G"
        assert _cache_key("X", 100, "A", "G") == "X:100:A:G"

    def test_get_predictor_singleton(self):
        """SPL-45: _get_predictor returns singleton."""
        reset_spliceai_cache()
        p1 = _get_predictor()
        p2 = _get_predictor()
        assert p1 is p2

    def test_reset_cache_clears_singleton(self):
        """SPL-46: reset_spliceai_cache clears cache and resets singleton."""
        reset_spliceai_cache()
        p1 = _get_predictor()
        p1._cache["key"] = SpliceAIResult()
        reset_spliceai_cache()
        p2 = _get_predictor()
        assert p1 is not p2
        assert p2._cache == {}

    def test_predictor_should_query_list(self):
        """SPL-47: SpliceAIPredictor.should_query with list input."""
        assert SpliceAIPredictor.should_query(["missense_variant"]) is False
        assert SpliceAIPredictor.should_query(["splice_acceptor_variant"]) is True
        assert SpliceAIPredictor.should_query(["synonymous_variant"], is_near_exon_boundary=True) is True
        assert SpliceAIPredictor.should_query(["synonymous_variant"], is_near_exon_boundary=False) is False

    def test_predictor_should_query_chinese(self):
        """SPL-48: SpliceAIPredictor.should_query with Chinese terms."""
        assert SpliceAIPredictor.should_query(["剪接区域变异"]) is True
        assert SpliceAIPredictor.should_query(["剪接供体第5位碱基变异"]) is True


# =============================================================================
# NEW: query_spliceai_batch pipeline function
# =============================================================================

@pytest.mark.l2
@pytest.mark.spliceai
class TestQuerySpliceaiBatch:
    """SPL-49~52: query_spliceai_batch pipeline-level function."""

    @pytest.mark.asyncio
    async def test_query_spliceai_batch_success(self):
        """SPL-49: Batch query succeeds for all variants."""
        reset_spliceai_cache()
        v1 = MagicMock(chrom="1", pos=100, ref="A", alt="G")
        v2 = MagicMock(chrom="1", pos=101, ref="C", alt="T")
        semaphore = asyncio.Semaphore(2)

        resp = _mock_aiohttp_response(status=200, json_data={
            "scores": [{"DS_AG": "0.1", "DS_AL": "0.0", "DS_DG": "0.0", "DS_DL": "0.0"}],
        })
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            results = await query_spliceai_batch([v1, v2], semaphore)
        assert len(results) == 2
        assert "1:100:A:G" in results
        assert "1:101:C:T" in results

    @pytest.mark.asyncio
    async def test_query_spliceai_batch_with_cache(self):
        """SPL-50: Cached variants are skipped in batch query."""
        reset_spliceai_cache()
        predictor = _get_predictor()
        predictor._cache["1:100:A:G"] = SpliceAIResult(delta_score=0.9)
        v1 = MagicMock(chrom="1", pos=100, ref="A", alt="G")
        v2 = MagicMock(chrom="1", pos=101, ref="C", alt="T")
        semaphore = asyncio.Semaphore(2)

        resp = _mock_aiohttp_response(status=200, json_data={
            "scores": [{"DS_AG": "0.1", "DS_AL": "0.0", "DS_DG": "0.0", "DS_DL": "0.0"}],
        })
        session = _mock_aiohttp_session(resp)

        with patch("aiohttp.ClientSession", return_value=session):
            results = await query_spliceai_batch([v1, v2], semaphore)
        assert results["1:100:A:G"].delta_score == 0.9
        assert results["1:101:C:T"].delta_score == 0.1

    @pytest.mark.asyncio
    async def test_query_spliceai_batch_exception(self):
        """SPL-51: Exception in batch query → api_error for that variant."""
        reset_spliceai_cache()
        v1 = MagicMock(chrom="1", pos=100, ref="A", alt="G")
        semaphore = asyncio.Semaphore(2)

        session = AsyncMock()
        session.get = MagicMock(side_effect=RuntimeError("network down"))

        with patch("aiohttp.ClientSession", return_value=session):
            results = await query_spliceai_batch([v1], semaphore)
        assert results["1:100:A:G"].source == "api_error"

    @pytest.mark.asyncio
    async def test_query_spliceai_batch_empty_input(self):
        """SPL-52: Empty input returns empty dict."""
        semaphore = asyncio.Semaphore(2)
        results = await query_spliceai_batch([], semaphore)
        assert results == {}
