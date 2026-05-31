"""
L2 Unit Tests — gpa_phenotype_match.py
Phenotype association: local DB lookup, fallback keyword match, LLM mock.

Run: pytest -m "l2 and phenotype" tests/l2_unit/test_phenotype_match.py
"""

import sys
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


@pytest.mark.l2
@pytest.mark.phenotype
@pytest.mark.p0
class TestPhenotypeMatcherInit:
    """Test PhenotypeMatcher initialization and local DB loading."""

    def test_init_default(self):
        """PM-01: Default init loads local DB if present."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        assert pm.model == "gpt-4o-mini"
        assert pm._local_db is not None  # dict, may be empty

    def test_init_with_api_key(self):
        """PM-02: api_key parameter is stored."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key="sk-test123")
        assert pm.api_key == "sk-test123"

    def test_init_no_local_db(self):
        """PM-03: Missing local DB → empty dict, no crash."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(refs_dir=Path("/nonexistent"))
        assert pm._local_db == {}


@pytest.mark.l2
@pytest.mark.phenotype
@pytest.mark.p0
class TestFallbackKeywordMatch:
    """Test _fallback_keyword_match: no-LLM degraded matching."""

    def test_exact_keyword_match(self):
        """PM-04: Exact keyword overlap → score > 0."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        result = pm._fallback_keyword_match(
            user_phenotypes="肌无力",
            known_phenotypes=["远端肌无力", "肌源性损害"]
        )
        assert result["score"] > 0
        assert len(result["matched_pairs"]) > 0
        assert result["confidence"] == "low"

    def test_no_overlap(self):
        """PM-05: No keyword overlap → score = 0."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        result = pm._fallback_keyword_match(
            user_phenotypes="视力丧失",
            known_phenotypes=["远端肌无力", "肌源性损害"]
        )
        assert result["score"] == 0.0
        assert result["matched_pairs"] == []

    def test_empty_known_phenotypes(self):
        """PM-06: Empty known phenotypes → score = 0."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        result = pm._fallback_keyword_match(
            user_phenotypes="肌无力",
            known_phenotypes=[]
        )
        assert result["score"] == 0.0
        assert "No known phenotypes" in result["explanation"]

    def test_multiple_delimiters(self):
        """PM-07: Handles comma, semicolon, Chinese comma delimiters."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        result = pm._fallback_keyword_match(
            user_phenotypes="肌无力,肌萎缩;肌痛、疲劳",
            known_phenotypes=["肌无力", "疲劳感"]
        )
        assert result["score"] > 0
        # Should match "肌无力" and "疲劳"
        assert len(result["matched_pairs"]) >= 1

    def test_case_insensitive(self):
        """PM-08: Matching is case-insensitive."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        result = pm._fallback_keyword_match(
            user_phenotypes="MUSCLE WEAKNESS",
            known_phenotypes=["muscle weakness", "myopathy"]
        )
        assert result["score"] > 0

    def test_warning_in_result(self):
        """PM-09: Fallback result includes warning about degraded accuracy."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        result = pm._fallback_keyword_match(
            user_phenotypes="肌无力",
            known_phenotypes=["远端肌无力"]
        )
        assert "No LLM API key" in result["reasoning"]


@pytest.mark.l2
@pytest.mark.phenotype
@pytest.mark.p0
class TestBuildMatchPrompt:
    """Test _build_match_prompt: prompt construction."""

    def test_prompt_contains_gene(self):
        """PM-10: Prompt contains gene symbol."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        prompt = pm._build_match_prompt("CAPN3", "肌无力", ["远端肌无力"])
        assert "CAPN3" in prompt
        assert "肌无力" in prompt
        assert "远端肌无力" in prompt

    def test_prompt_requests_json(self):
        """PM-11: Prompt requests JSON output with expected schema."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher()
        prompt = pm._build_match_prompt("CAPN3", "肌无力", ["远端肌无力"])
        assert "json_object" in prompt or "JSON" in prompt
        assert "score" in prompt
        assert "matched_pairs" in prompt


@pytest.mark.l2
@pytest.mark.phenotype
@pytest.mark.asyncio
class TestAsyncMatch:
    """Test match() async flow with mocked LLM."""

    async def test_no_api_key_uses_fallback(self):
        """PM-12: No API key → fallback keyword match, includes warning."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key=None)
        # Seed local db with a known entry
        pm._local_db = {"CAPN3": {"phenotypes": [{"name": "远端肌无力"}]}}
        result = await pm.match("CAPN3", "肌无力")
        assert "warning" in result
        assert "fallback" in result["warning"].lower()
        assert result["gene"] == "CAPN3"
        assert result["user_phenotypes"] == "肌无力"

    async def test_local_db_cache(self):
        """PM-13: Second call for same gene uses cache."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key=None)
        pm._local_db = {"CAPN3": {"phenotypes": [{"name": "远端肌无力"}]}}
        # First call
        r1 = await pm.match("CAPN3", "肌无力")
        # Second call
        r2 = await pm.match("CAPN3", "肌萎缩")
        assert "CAPN3" in pm.gene_phenotype_cache
        assert r2["known_phenotypes"] == ["远端肌无力"]

    async def test_unknown_gene_empty_phenotypes(self):
        """PM-14: Gene not in DB → empty known_phenotypes, score=0."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key=None)
        pm._local_db = {}
        result = await pm.match("UNKNOWN_GENE", "some phenotype")
        assert result["known_phenotypes"] == []
        assert result["score"] == 0.0

    async def test_llm_mock_success(self):
        """PM-15: Mocked LLM API returns parsed JSON result."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key="sk-test")
        pm._local_db = {"CAPN3": {"phenotypes": [{"name": "远端肌无力"}]}}

        # Patch _llm_semantic_match directly to avoid aiohttp mocking complexity
        async def mock_llm_match(gene, user_phenotypes, known_phenotypes):
            return {
                "score": 0.85,
                "matched_pairs": [["肌无力", "远端肌无力"]],
                "explanation": "semantic match",
                "confidence": "high",
                "reasoning": "LLM judged high similarity",
            }

        pm._llm_semantic_match = mock_llm_match
        result = await pm.match("CAPN3", "肌无力")

        assert result["score"] == 0.85
        assert result["confidence"] == "high"
        assert result["gene"] == "CAPN3"

    async def test_llm_api_error(self):
        """PM-16: LLM API error → graceful fallback with score=0."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key="sk-test")
        pm._local_db = {"CAPN3": {"phenotypes": [{"name": "远端肌无力"}]}}

        async def mock_llm_error(gene, user_phenotypes, known_phenotypes):
            return {
                "score": 0.0,
                "matched_pairs": [],
                "explanation": "LLM API error (HTTP 429): Rate limited",
                "confidence": "low",
                "reasoning": "API request failed.",
            }

        pm._llm_semantic_match = mock_llm_error
        result = await pm.match("CAPN3", "肌无力")

        assert result["score"] == 0.0
        assert "429" in result["explanation"]
        assert result["confidence"] == "low"
        assert result["confidence"] == "low"

    async def test_match_batch(self):
        """PM-17: Batch match returns results for all genes."""
        from gpa_phenotype_match import PhenotypeMatcher
        pm = PhenotypeMatcher(llm_api_key=None)
        pm._local_db = {
            "CAPN3": {"phenotypes": [{"name": "远端肌无力"}]},
            "DYSF": {"phenotypes": [{"name": "肢带型肌营养不良"}]},
        }
        results = await pm.match_batch(["CAPN3", "DYSF"], "肌无力")
        assert len(results) == 2
        assert results[0]["gene"] == "CAPN3"
        assert results[1]["gene"] == "DYSF"
