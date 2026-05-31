"""
L2 Unit Tests — dgra_clinvar.py
ClinVar query logic: consequence filtering, XML parsing, position matching.

Run: pytest -m "l2 and clinvar" tests/l2_unit/test_clinvar.py
"""

import sys
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


@pytest.mark.l2
@pytest.mark.clinvar
@pytest.mark.p0
class TestVariantNeedsClinvar:
    """Test variant_needs_clinvar: consequence-based filtering."""

    def test_missense_needs_clinvar(self):
        """CV-01: missense_variant → True."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("missense_variant") is True

    def test_splice_donor_needs_clinvar(self):
        """CV-02: splice_donor_variant → True."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("splice_donor_variant") is True

    def test_synonymous_skips_clinvar(self):
        """CV-03: synonymous_variant → False."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("synonymous_variant") is False

    def test_intron_skips_clinvar(self):
        """CV-04: intron_variant → False."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("intron_variant") is False

    def test_compound_relevant(self):
        """CV-05: Compound with relevant term → True."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("missense_variant,splice_region_variant") is True

    def test_compound_irrelevant(self):
        """CV-06: Compound with only irrelevant terms → False."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("synonymous_variant,intron_variant") is False

    def test_unknown_consequence_defaults_true(self):
        """CV-07: Unknown consequence → True (conservative)."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("some_new_consequence") is True

    def test_empty_consequence(self):
        """CV-08: Empty consequence → False."""
        from dgra_clinvar import variant_needs_clinvar
        assert variant_needs_clinvar("") is False
        assert variant_needs_clinvar("_UNKNOWN_") is False


@pytest.mark.l2
@pytest.mark.clinvar
@pytest.mark.p0
class TestParseClinvarXML:
    """Test _parse_clinvar_xml: XML parsing."""

    def test_parse_valid_xml(self):
        """CV-09: Valid VCV XML → parsed dict."""
        from dgra_clinvar import _parse_clinvar_xml
        xml = """<?xml version="1.0"?>
        <ClinVarResult>
            <VariationArchive Accession="VCV000012345" VariationName="NM_001.1(GENE):c.100A>G" VariationType="SNV"/>
            <ClassifiedRecord>
                <ReviewStatus>criteria provided, single submitter</ReviewStatus>
                <ClinicalAssertionList>
                    <ClinicalAssertion>
                        <Interpretation>
                            <Description DateLastEvaluated="2024-01-01">Pathogenic</Description>
                        </Interpretation>
                    </ClinicalAssertion>
                </ClinicalAssertionList>
            </ClassifiedRecord>
            <SequenceLocation Assembly="GRCh38" Chr="1" start="100" stop="100"
                referenceAlleleVCF="A" alternateAlleleVCF="G" positionVCF="100"/>
        </ClinVarResult>"""
        result = _parse_clinvar_xml(xml)
        assert result is not None
        assert result["accession"] == "VCV000012345"
        assert result["clinical_significance"] == "Pathogenic"
        assert result["review_status"] == "criteria provided, single submitter"
        assert len(result["positions"]) == 1
        assert result["positions"][0]["chr"] == "1"

    def test_parse_invalid_xml(self):
        """CV-10: Invalid XML → None."""
        from dgra_clinvar import _parse_clinvar_xml
        result = _parse_clinvar_xml("not xml at all")
        assert result is None

    def test_parse_no_clinical_significance(self):
        """CV-11: XML without clinical significance → None for that field."""
        from dgra_clinvar import _parse_clinvar_xml
        xml = """<?xml version="1.0"?>
        <ClinVarResult>
            <VariationArchive Accession="VCV000012345"/>
        </ClinVarResult>"""
        result = _parse_clinvar_xml(xml)
        assert result is not None
        assert result["clinical_significance"] is None
        assert result["accession"] == "VCV000012345"

    def test_parse_grch38_filter(self):
        """CV-12: Only GRCh38 positions are extracted."""
        from dgra_clinvar import _parse_clinvar_xml
        xml = """<?xml version="1.0"?>
        <ClinVarResult>
            <SequenceLocation Assembly="GRCh37" Chr="1" start="100" stop="100"/>
            <SequenceLocation Assembly="GRCh38" Chr="1" start="110" stop="110"/>
        </ClinVarResult>"""
        result = _parse_clinvar_xml(xml)
        assert len(result["positions"]) == 1
        assert result["positions"][0]["start"] == "110"


@pytest.mark.l2
@pytest.mark.clinvar
@pytest.mark.p0
class TestCheckPositionMatch:
    """Test _check_position_match: position matching logic."""

    def test_exact_match(self):
        """CV-13: Exact position match → True, exact."""
        from dgra_clinvar import _check_position_match
        positions = [{"chr": "1", "start": "100", "stop": "100", "position_vcf": "100"}]
        is_match, quality = _check_position_match("1", 100, positions)
        assert is_match is True
        assert quality == "exact"

    def test_overlap_match(self):
        """CV-14: Position within span → True, overlap."""
        from dgra_clinvar import _check_position_match
        positions = [{"chr": "1", "start": "90", "stop": "110", "position_vcf": None}]
        is_match, quality = _check_position_match("1", 100, positions)
        assert is_match is True
        assert quality == "overlap"

    def test_nearby_match(self):
        """CV-15: Position within 5bp → True, nearby."""
        from dgra_clinvar import _check_position_match
        # start/stop far from query_pos, only position_vcf is nearby
        positions = [{"chr": "1", "start": "1000", "stop": "1000", "position_vcf": "103"}]
        is_match, quality = _check_position_match("1", 100, positions)
        assert is_match is True
        assert quality == "nearby"

    def test_no_match(self):
        """CV-16: Position far away → False, no_match."""
        from dgra_clinvar import _check_position_match
        positions = [{"chr": "1", "start": "1000", "stop": "1000", "position_vcf": "1000"}]
        is_match, quality = _check_position_match("1", 100, positions)
        assert is_match is False
        assert quality == "no_match"

    def test_chromosome_mismatch(self):
        """CV-17: Different chromosome → no_match."""
        from dgra_clinvar import _check_position_match
        positions = [{"chr": "2", "start": "100", "stop": "100", "position_vcf": "100"}]
        is_match, quality = _check_position_match("1", 100, positions)
        assert is_match is False

    def test_chr_prefix_stripped(self):
        """CV-18: chr prefix is stripped for comparison."""
        from dgra_clinvar import _check_position_match
        positions = [{"chr": "1", "start": "100", "stop": "100", "position_vcf": "100"}]
        is_match, quality = _check_position_match("chr1", 100, positions)
        assert is_match is True
        assert quality == "exact"

    def test_empty_positions(self):
        """CV-19: Empty positions list → no_match."""
        from dgra_clinvar import _check_position_match
        is_match, quality = _check_position_match("1", 100, [])
        assert is_match is False
        assert quality == "no_match"


@pytest.mark.l2
@pytest.mark.clinvar
@pytest.mark.asyncio
class TestQueryClinvarVariant:
    """Test query_clinvar_variant: full flow with mocked HTTP."""

    async def test_skips_irrelevant_consequence(self):
        """CV-20: Irrelevant consequence → skipped without API call."""
        from dgra_clinvar import query_clinvar_variant
        mock_session = MagicMock()
        result = await query_clinvar_variant("1", 100, "A", "G", "BRCA1", "synonymous_variant", mock_session)
        assert result.status == "skipped"
        assert "not relevant" in result.error.lower()
        mock_session.get.assert_not_called()

    async def test_not_found_empty_ids(self):
        """CV-21: ESearch returns empty → not_found."""
        from dgra_clinvar import query_clinvar_variant

        async def mock_esearch(*args, **kwargs):
            return [], None

        with patch("dgra_clinvar._ncbi_esearch_clinvar", side_effect=mock_esearch):
            mock_session = MagicMock()
            result = await query_clinvar_variant("1", 100, "A", "G", "BRCA1", "missense_variant", mock_session)
        assert result.status == "not_found"

    async def test_esearch_error(self):
        """CV-22: ESearch HTTP error → error status."""
        from dgra_clinvar import query_clinvar_variant

        async def mock_esearch(*args, **kwargs):
            return [], "HTTP 503: Service Unavailable"

        with patch("dgra_clinvar._ncbi_esearch_clinvar", side_effect=mock_esearch):
            mock_session = MagicMock()
            result = await query_clinvar_variant("1", 100, "A", "G", "BRCA1", "missense_variant", mock_session)
        assert result.status == "error"
        assert "ESearch failed" in result.error

    async def test_success_exact_match(self):
        """CV-23: Full flow with exact position match → success."""
        from dgra_clinvar import query_clinvar_variant

        async def mock_esearch(*args, **kwargs):
            return ["12345"], None

        async def mock_efetch(*args, **kwargs):
            return {
                "accession": "VCV000012345",
                "variant_name": "NM_001.1(BRCA1):c.100A>G",
                "variant_type": "SNV",
                "clinical_significance": "Pathogenic",
                "review_status": "criteria provided, single submitter",
                "positions": [
                    {"chr": "1", "start": "100", "stop": "100", "position_vcf": "100",
                     "ref": "A", "alt": "G"}
                ],
            }

        with patch("dgra_clinvar._ncbi_esearch_clinvar", side_effect=mock_esearch), \
             patch("dgra_clinvar._ncbi_efetch_clinvar", side_effect=mock_efetch):
            mock_session = MagicMock()
            result = await query_clinvar_variant("1", 100, "A", "G", "BRCA1", "missense_variant", mock_session)

        assert result.status == "success"
        assert result.position_match is True
        assert result.match_quality == "exact"
        assert result.clinvar_significance == "Pathogenic"
        assert result.clinvar_accession == "VCV000012345"

    async def test_batch_query(self):
        """CV-24: Batch query returns results for all variants."""
        from dgra_clinvar import query_clinvar_batch, ClinVarResult

        # Mock query_clinvar_variant
        async def mock_query(chrom, pos, ref, alt, gene, consequence, session):
            return ClinVarResult(
                variant_id=f"{chrom}:{pos}_{ref}>{alt}",
                gene=gene,
                status="success",
                clinvar_significance="Pathogenic",
            )

        mock_session = MagicMock()

        with patch("dgra_clinvar.query_clinvar_variant", side_effect=mock_query):
            variants = [
                ("1", 100, "A", "G", "BRCA1", "missense_variant"),
                ("17", 7579472, "C", "T", "TP53", "stop_gained"),
            ]
            results = await query_clinvar_batch(variants, mock_session)

        assert len(results) == 2
        assert results["1:100_A>G"].status == "success"
        assert results["17:7579472_C>T"].status == "success"
