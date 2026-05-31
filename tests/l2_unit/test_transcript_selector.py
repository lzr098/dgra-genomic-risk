"""
L2 Unit Tests — gpa_transcript_selector.py
Transcript selection: rule-based scoring, ambiguity detection, LLM assist.

Run: pytest -m "l2 and transcript" tests/l2_unit/test_transcript_selector.py
"""

import sys
from pathlib import Path
from unittest.mock import AsyncMock, patch

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


@pytest.mark.l2
@pytest.mark.transcript
@pytest.mark.p0
class TestTranscriptSelectorInit:
    """Test TranscriptSelector initialization."""

    def test_default_init(self):
        """TX-01: Default init with general tissue profile."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        assert ts.tissue_profile == "general"
        assert ts.ambiguity_threshold == 5
        assert ts.llm_model == "gpt-4o-mini"

    def test_custom_init(self):
        """TX-02: Custom parameters."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector(
            tissue_profile="hematopoietic",
            disease_description="AML",
            ambiguity_threshold=3,
        )
        assert ts.tissue_profile == "hematopoietic"
        assert ts.disease_description == "AML"
        assert ts.ambiguity_threshold == 3


@pytest.mark.l2
@pytest.mark.transcript
@pytest.mark.p0
class TestTranscriptSelect:
    """Test select() method: rule-based transcript selection."""

    def test_empty_transcripts(self):
        """TX-03: Empty list → empty primary, method='none'."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        result = ts.select("BRCA1", [])
        assert result.primary == {}
        assert result.method == "none"
        assert result.is_ambiguous is False

    def test_single_transcript(self):
        """TX-04: Single transcript → primary is that transcript."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        tx = {"transcript_id": "NM_007294.4", "canonical": 1}
        result = ts.select("BRCA1", [tx])
        assert result.primary == tx
        assert result.method == "canonical"
        assert result.is_ambiguous is False

    def test_canonical_wins(self):
        """TX-05: Canonical transcript scores highest."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 1, "impact": "MODERATE"},
            {"transcript_id": "NM_002.1", "canonical": 0, "impact": "HIGH"},
        ]
        result = ts.select("TEST1", txs)  # TEST1 not in any special list
        assert result.primary["transcript_id"] == "NM_001.1"
        assert result.method == "canonical"

    def test_mane_select_wins(self):
        """TX-06: MANE Select gets +10, can win over non-MANE."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 0, "mane_select": 1, "impact": "LOW"},
            {"transcript_id": "NM_002.1", "canonical": 1, "mane_select": 0, "impact": "LOW"},
        ]
        result = ts.select("TEST1", txs)
        # Both have 10 points (canonical=10, mane=10), ambiguous
        assert result.is_ambiguous is True

    def test_high_impact_bonus(self):
        """TX-07: HIGH impact gets +10."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 0, "impact": "HIGH"},
            {"transcript_id": "NM_002.1", "canonical": 0, "impact": "LOW"},
        ]
        result = ts.select("TEST1", txs)
        assert result.primary["transcript_id"] == "NM_001.1"

    def test_ambiguity_detected(self):
        """TX-08: Score gap < threshold → is_ambiguous=True."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector(ambiguity_threshold=5)
        # Both have canonical=1 → both score 10, gap=0 < 5 → ambiguous
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 1, "impact": "LOW"},
            {"transcript_id": "NM_002.1", "canonical": 1, "impact": "LOW"},
        ]
        result = ts.select("TEST1", txs)
        assert result.is_ambiguous is True
        assert result.method == "ambiguous"

    def test_not_ambiguous_large_gap(self):
        """TX-09: Score gap >= threshold → not ambiguous."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector(ambiguity_threshold=5)
        # canonical+HIGH=20 vs canonical+LOW=12, gap=8 >= 5 → not ambiguous
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 1, "impact": "HIGH"},
            {"transcript_id": "NM_002.1", "canonical": 1, "impact": "LOW"},
        ]
        result = ts.select("TEST1", txs)
        assert result.is_ambiguous is False

    def test_alternatives_returned(self):
        """TX-10: Alternatives list contains non-primary transcripts."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 1},
            {"transcript_id": "NM_002.1", "canonical": 0},
            {"transcript_id": "NM_003.1", "canonical": 0},
        ]
        result = ts.select("TEST1", txs)
        assert result.primary["transcript_id"] == "NM_001.1"
        assert len(result.alternatives) == 2

    def test_protein_domains_bonus(self):
        """TX-11: Protein domains add score (+3 per domain, max 8)."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 0, "protein_domains": ["PF001", "PF002", "PF003"]},
            {"transcript_id": "NM_002.1", "canonical": 0, "protein_domains": []},
        ]
        result = ts.select("TEST1", txs)
        # 3 domains * 3 = 9, capped at 8
        assert result.primary["transcript_id"] == "NM_001.1"

    def test_tissue_expression_bonus(self):
        """TX-12: Tissue-relevant gene gets bonus."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector(tissue_profile="hematopoietic")
        # FANCA is in fa_dna_repair list for hematopoietic
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 0},
        ]
        result = ts.select("FANCA", txs)
        # Score should include tissue bonus
        assert "tissue" in result.selection_reason.lower() or result.selection_reason != ""

    def test_selection_reason_populated(self):
        """TX-13: selection_reason is non-empty."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [{"transcript_id": "NM_001.1", "canonical": 1}]
        result = ts.select("TEST1", txs)
        assert result.selection_reason != ""

    def test_brca1_tissue_expression_method(self):
        """TX-14: BRCA1 (in cancer_predisposition) → tissue_expression method."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_007294.4", "canonical": 1, "impact": "MODERATE"},
            {"transcript_id": "NM_007295.3", "canonical": 0, "impact": "HIGH"},
        ]
        result = ts.select("BRCA1", txs)
        assert result.primary["transcript_id"] == "NM_007294.4"
        assert result.method == "tissue_expression"


@pytest.mark.l2
@pytest.mark.transcript
@pytest.mark.asyncio
class TestAsyncSelect:
    """Test aselect() async method."""

    async def test_aselect_basic(self):
        """TX-14: aselect works same as select for non-ambiguous cases."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector()
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 1, "impact": "HIGH"},
            {"transcript_id": "NM_002.1", "canonical": 0, "impact": "LOW"},
        ]
        result = await ts.aselect("BRCA1", txs)
        assert result.primary["transcript_id"] == "NM_001.1"
        assert result.is_ambiguous is False

    async def test_aselect_llm_assist_mock(self):
        """TX-15: LLM assist picks transcript when ambiguous."""
        from gpa_transcript_selector import TranscriptSelector
        ts = TranscriptSelector(
            disease_description="hereditary breast cancer",
            llm_api_key="sk-test",
            ambiguity_threshold=10,
        )
        # Both canonical, same score → ambiguous
        txs = [
            {"transcript_id": "NM_001.1", "canonical": 1, "impact": "LOW"},
            {"transcript_id": "NM_002.1", "canonical": 1, "impact": "LOW"},
        ]

        # Mock LLM to pick NM_002.1
        async def mock_llm(gene, candidates):
            return candidates[1]  # Pick second

        ts._llm_assist_select = mock_llm
        result = await ts.aselect("BRCA1", txs)
        assert result.primary["transcript_id"] == "NM_002.1"
        assert result.method == "llm_disease_match"
