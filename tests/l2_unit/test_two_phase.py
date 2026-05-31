"""L2 Unit Tests — gpa_two_phase.py

Covers Phase 1 parsing, pathogenic pre-filter, fast tier classification,
candidate priority scoring, variant need checks, Phase 2 enrichment,
and the full two-phase pipeline entry points.
"""

import pytest
from unittest.mock import AsyncMock, MagicMock, patch
from conftest import make_variant


@pytest.mark.l2
class TestParseVariantsPhase1:
    """TP-01~07: _parse_variants_phase1."""

    def test_basic_parse(self):
        """TP-01: Basic variant dict → Variant object."""
        from gpa_two_phase import _parse_variants_phase1
        data = [{
            "CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
            "GENE": "BRCA1", "Feature": "NM_001.1",
            "IMPACT": "HIGH", "Consequence": "stop_gained",
            "CLIN_SIG": "Pathogenic", "DP": 30, "GQ": 99,
            "GT": "0/1", "VAF": "0.45", "gnomAD_AF": "0.001",
        }]
        variants = _parse_variants_phase1(data)
        assert len(variants) == 1
        v = variants[0]
        assert v.chrom == "1"
        assert v.pos == 100
        assert v.gene == "BRCA1"
        assert v.impact == "HIGH"
        assert v.consequence == "stop_gained"
        assert v.clinvar == "Pathogenic"
        assert v.dp == 30
        assert v.gq == 99.0
        assert v.gt == "0/1"
        assert v.vaf == 0.45
        assert v.gnomad_af == 0.001

    def test_chinese_impact_mapping(self):
        """TP-02: Chinese impact terms mapped."""
        from gpa_two_phase import _parse_variants_phase1
        data = [{"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
                 "GENE": "T", "IMPACT": "高", "Consequence": "错义变异"}]
        v = _parse_variants_phase1(data)[0]
        assert v.impact == "HIGH"
        assert v.consequence == "missense_variant"

    def test_missing_fields(self):
        """TP-03: Missing fields → defaults."""
        from gpa_two_phase import _parse_variants_phase1, _UNKNOWN
        data = [{"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G", "GENE": "T"}]
        v = _parse_variants_phase1(data)[0]
        assert v.impact == _UNKNOWN
        assert v.consequence == _UNKNOWN
        assert v.clinvar == _UNKNOWN
        assert v.dp == 0
        assert v.gq == 0.0
        assert v.vaf is None
        assert v.gnomad_af is None

    def test_invalid_dp_gq(self):
        """TP-04: Invalid DP/GQ → 0."""
        from gpa_two_phase import _parse_variants_phase1
        data = [{"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G", "GENE": "T",
                 "DP": "abc", "GQ": "xyz", "VAF": "bad"}]
        v = _parse_variants_phase1(data)[0]
        assert v.dp == 0
        assert v.gq == 0.0
        assert v.vaf is None

    def test_dot_vaf(self):
        """TP-05: VAF='.' → None."""
        from gpa_two_phase import _parse_variants_phase1
        data = [{"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G", "GENE": "T", "VAF": "."}]
        v = _parse_variants_phase1(data)[0]
        assert v.vaf is None

    def test_gnomad_na(self):
        """TP-06: gnomAD_AF='N/A' → None."""
        from gpa_two_phase import _parse_variants_phase1
        data = [{"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G", "GENE": "T", "gnomAD_AF": "N/A"}]
        v = _parse_variants_phase1(data)[0]
        assert v.gnomad_af is None

    def test_empty_list(self):
        """TP-07: Empty list → empty result."""
        from gpa_two_phase import _parse_variants_phase1
        assert _parse_variants_phase1([]) == []


@pytest.mark.l2
class TestIsPotentiallyPathogenic:
    """TP-08~14: _is_potentially_pathogenic."""

    def test_high_impact_true(self):
        """TP-08: HIGH impact → True."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="RUNX1", impact="HIGH")
        assert _is_potentially_pathogenic(v) is True

    def test_high_common_homozygous_false(self):
        """TP-09: HIGH + common homozygous (AF>0.95, 1/1) → False."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="HIGH", gnomad_af=0.99, gt="1/1")
        assert _is_potentially_pathogenic(v) is False

    def test_moderate_rare_missense_true(self):
        """TP-10: MODERATE + rare missense + unknown AF → True (conservative include)."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="MODERATE", consequence="missense_variant", gnomad_af=None)
        assert _is_potentially_pathogenic(v) is True

    def test_moderate_common_missense_false(self):
        """TP-11: MODERATE + common missense (AF>1%) → False."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="MODERATE", consequence="missense_variant", gnomad_af=0.05)
        assert _is_potentially_pathogenic(v) is False

    def test_neuromuscular_exempt(self):
        """TP-12: Neuromuscular gene exempt from common filter."""
        from gpa_two_phase import _is_potentially_pathogenic, _NEUROMUSCULAR_DISEASE_GENES
        # Pick a known neuromuscular gene if available
        gene = next(iter(_NEUROMUSCULAR_DISEASE_GENES), "DMD")
        v = make_variant(gene=gene, impact="MODERATE", consequence="missense_variant", gnomad_af=0.05)
        assert _is_potentially_pathogenic(v) is True

    def test_splice_relevant_true(self):
        """TP-13: Splice-relevant consequence → True."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="LOW", consequence="splice_region_variant")
        assert _is_potentially_pathogenic(v) is True

    def test_low_modifier_false(self):
        """TP-14: LOW/MODIFIER non-splice → False."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="LOW", consequence="synonymous_variant")
        assert _is_potentially_pathogenic(v) is False


@pytest.mark.l2
class TestFastTierClassification:
    """TP-15~20: _fast_tier_classification."""

    def test_prefiltered_tier3(self):
        """TP-15: Not potentially pathogenic → Tier 3."""
        from gpa_two_phase import _fast_tier_classification
        v = make_variant(gene="TEST", impact="LOW", consequence="synonymous_variant")
        tier, reason = _fast_tier_classification(v)
        assert tier == 3
        assert "Pre-filtered" in reason

    def test_high_impact_tier1(self):
        """TP-16: HIGH impact → Tier 1 preliminary."""
        from gpa_two_phase import _fast_tier_classification
        v = make_variant(gene="RUNX1", impact="HIGH")
        tier, reason = _fast_tier_classification(v)
        assert tier == 1
        assert "HIGH impact" in reason
        assert "PRELIMINARY" in reason

    def test_pathogenic_clinvar_tier1(self):
        """TP-17: ClinVar Pathogenic → Tier 1."""
        from gpa_two_phase import _fast_tier_classification
        v = make_variant(gene="BRCA1", impact="MODERATE", clinvar="Pathogenic")
        tier, reason = _fast_tier_classification(v)
        assert tier == 1
        assert "ClinVar Pathogenic" in reason

    def test_likely_pathogenic_tier1(self):
        """TP-18: ClinVar Likely_pathogenic → Tier 1."""
        from gpa_two_phase import _fast_tier_classification
        v = make_variant(gene="BRCA1", impact="MODERATE", clinvar="Likely_pathogenic")
        tier, reason = _fast_tier_classification(v)
        assert tier == 1

    def test_splice_relevant_tier1(self):
        """TP-19: Splice-relevant consequence → Tier 1."""
        from gpa_two_phase import _fast_tier_classification
        v = make_variant(gene="TEST", impact="LOW", consequence="splice_region_variant")
        tier, reason = _fast_tier_classification(v)
        assert tier == 1
        assert "Splice-relevant" in reason

    def test_candidate_tier2(self):
        """TP-20: Candidate but not Tier 1 criteria → Tier 2."""
        from gpa_two_phase import _fast_tier_classification
        v = make_variant(gene="TEST", impact="MODERATE", consequence="missense_variant")
        tier, reason = _fast_tier_classification(v)
        assert tier == 2
        assert "PRELIMINARY" in reason


@pytest.mark.l2
class TestCandidatePriorityScore:
    """TP-21~26: _candidate_priority_score."""

    def test_tier_priority(self):
        """TP-21: Tier 1 < Tier 2 in tuple (lower = higher priority)."""
        from gpa_two_phase import _candidate_priority_score
        v1 = make_variant(gene="A", impact="HIGH")
        v1.tier = 1
        v2 = make_variant(gene="B", impact="HIGH")
        v2.tier = 2
        assert _candidate_priority_score(v1) < _candidate_priority_score(v2)

    def test_impact_priority(self):
        """TP-22: HIGH < MODERATE < LOW."""
        from gpa_two_phase import _candidate_priority_score
        v1 = make_variant(gene="A", impact="HIGH")
        v2 = make_variant(gene="B", impact="MODERATE")
        v3 = make_variant(gene="C", impact="LOW")
        assert _candidate_priority_score(v1) < _candidate_priority_score(v2)
        assert _candidate_priority_score(v2) < _candidate_priority_score(v3)

    def test_clinvar_priority(self):
        """TP-23: Pathogenic < non-pathogenic."""
        from gpa_two_phase import _candidate_priority_score
        v1 = make_variant(gene="A", impact="MODERATE", clinvar="Pathogenic")
        v2 = make_variant(gene="B", impact="MODERATE", clinvar="Uncertain_significance")
        assert _candidate_priority_score(v1) < _candidate_priority_score(v2)

    def test_af_priority(self):
        """TP-24: Unknown < rare < uncommon < common."""
        from gpa_two_phase import _candidate_priority_score
        v1 = make_variant(gene="A", impact="MODERATE")
        v1.gnomad_af = None
        v2 = make_variant(gene="B", impact="MODERATE")
        v2.gnomad_af = 0.0005
        v3 = make_variant(gene="C", impact="MODERATE")
        v3.gnomad_af = 0.005
        v4 = make_variant(gene="D", impact="MODERATE")
        v4.gnomad_af = 0.05
        assert _candidate_priority_score(v1) < _candidate_priority_score(v2)
        assert _candidate_priority_score(v2) < _candidate_priority_score(v3)
        assert _candidate_priority_score(v3) < _candidate_priority_score(v4)

    def test_chinese_pathogenic(self):
        """TP-25: Chinese pathogenic term recognized."""
        from gpa_two_phase import _candidate_priority_score
        v = make_variant(gene="A", impact="MODERATE", clinvar="致病")
        score = _candidate_priority_score(v)
        assert score[2] == 0  # clinvar_score = 0 for pathogenic

    def test_score_tuple_length(self):
        """TP-26: Score is a 4-tuple."""
        from gpa_two_phase import _candidate_priority_score
        v = make_variant(gene="A", impact="MODERATE")
        score = _candidate_priority_score(v)
        assert len(score) == 4


@pytest.mark.l2
class TestVariantNeedsChecks:
    """TP-27~30: _variant_needs_clinvar / _variant_needs_spliceai."""

    def test_needs_clinvar_pathogenic(self):
        """TP-27: Pathogenic variant does NOT need ClinVar (already known)."""
        from gpa_two_phase import _variant_needs_clinvar
        v = make_variant(gene="BRCA1", clinvar="Pathogenic")
        assert _variant_needs_clinvar(v) is False

    def test_needs_clinvar_unknown(self):
        """TP-28: Missense consequence → needs ClinVar lookup (consequence-based)."""
        from gpa_two_phase import _variant_needs_clinvar
        v = make_variant(gene="BRCA1", impact="MODERATE", consequence="missense_variant")
        assert _variant_needs_clinvar(v) is True

    def test_needs_spliceai_splice(self):
        """TP-29: Splice consequence → needs SpliceAI."""
        from gpa_two_phase import _variant_needs_spliceai
        v = make_variant(gene="BRCA1", consequence="splice_donor_variant")
        assert _variant_needs_spliceai(v) is True

    def test_needs_spliceai_missense(self):
        """TP-30: Missense → does NOT need SpliceAI."""
        from gpa_two_phase import _variant_needs_spliceai
        v = make_variant(gene="BRCA1", consequence="missense_variant")
        assert _variant_needs_spliceai(v) is False


# =============================================================================
# Expanded coverage — Phase 2 enrichment + pipeline entry points
# =============================================================================

@pytest.mark.l2
@pytest.mark.mock
class TestIsPotentiallyPathogenicExtended:
    """TP-31~38: Additional edge cases for _is_potentially_pathogenic."""

    def test_moderate_clinvar_pathogenic_true(self):
        """TP-31: MODERATE + ClinVar pathogenic → True."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="MODERATE", clinvar="Pathogenic")
        assert _is_potentially_pathogenic(v) is True

    def test_moderate_rare_af_true(self):
        """TP-32: MODERATE + AF < 0.001 → True."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="MODERATE", consequence="missense_variant", gnomad_af=0.0001)
        assert _is_potentially_pathogenic(v) is True

    def test_moderate_neuromuscular_common_above_5pct_false(self):
        """TP-33: Neuromuscular gene + AF > 5% → False."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="DMD", impact="MODERATE", consequence="missense_variant", gnomad_af=0.1)
        assert _is_potentially_pathogenic(v) is False

    def test_moderate_neuromuscular_2pct_true(self):
        """TP-34: Neuromuscular gene + AF 2% (below 5%) → True."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="DMD", impact="MODERATE", consequence="missense_variant", gnomad_af=0.02)
        assert _is_potentially_pathogenic(v) is True

    def test_low_clinvar_pathogenic_true(self):
        """TP-35: LOW impact + ClinVar pathogenic → True."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="LOW", consequence="synonymous_variant", clinvar="Pathogenic")
        assert _is_potentially_pathogenic(v) is True

    def test_high_af_99_hom_false(self):
        """TP-36: HIGH + AF 0.99 + 1/1 → False."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="HIGH", gnomad_af=0.99, gt="1/1")
        assert _is_potentially_pathogenic(v) is False

    def test_high_af_99_het_true(self):
        """TP-37: HIGH + AF 0.99 + 0/1 → True (not homozygous)."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="HIGH", gnomad_af=0.99, gt="0/1")
        assert _is_potentially_pathogenic(v) is True

    def test_modifier_empty_clinvar_false(self):
        """TP-38: MODIFIER + empty ClinVar → False."""
        from gpa_two_phase import _is_potentially_pathogenic
        v = make_variant(gene="TEST", impact="MODIFIER", consequence="upstream_gene_variant", clinvar="")
        assert _is_potentially_pathogenic(v) is False


@pytest.mark.l2
@pytest.mark.mock
class TestEnrichCandidateGenes:
    """TP-39~42: _enrich_candidate_genes async enrichment."""

    @pytest.mark.asyncio
    async def test_empty_candidates(self):
        """TP-39: Empty candidate_indices → returns empty dicts."""
        from gpa_two_phase import _enrich_candidate_genes
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        result = await _enrich_candidate_genes([], [], config, {})
        assert result == ({}, {}, {})

    @pytest.mark.asyncio
    async def test_offline_mode_loads_archive(self):
        """TP-40: Offline mode loads archived data per gene."""
        from gpa_two_phase import _enrich_candidate_genes
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        config.offline_mode = True
        v = make_variant(gene="TP53")
        archive = {
            "ensembl": {"gene": "TP53"},
            "uniprot": {"entry": "P04637"},
            "gnomad_constraint": {"status": "CAPTURED", "pLI": 1.0},
        }
        with patch("gpa_two_phase._load_offline_archive", return_value=archive):
            ensembl, uniprot, constraint = await _enrich_candidate_genes([v], [0], config, {})
        assert ensembl.get("TP53") == {"gene": "TP53"}
        assert constraint.get("TP53", {}).get("pLI") == 1.0

    @pytest.mark.asyncio
    async def test_online_mode_api_success(self):
        """TP-41: Online mode queries APIs and maps results by gene."""
        from gpa_two_phase import _enrich_candidate_genes
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        config.offline_mode = False
        v = make_variant(gene="TP53")

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.batch_query_genes = AsyncMock(return_value={"TP53": {"status": "SUCCESS"}})

        with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
            ensembl, uniprot, constraint = await _enrich_candidate_genes([v], [0], config, {})
        assert ensembl.get("TP53") == {"status": "SUCCESS"}
        assert uniprot.get("TP53") == {"status": "SUCCESS"}

    @pytest.mark.asyncio
    async def test_online_mode_api_partial(self):
        """TP-42: API returns only some genes → missing genes get empty dict."""
        from gpa_two_phase import _enrich_candidate_genes
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        config.offline_mode = False
        v1 = make_variant(gene="TP53")
        v2 = make_variant(gene="BRCA1")

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.batch_query_genes = AsyncMock(return_value={"TP53": {"status": "SUCCESS"}})

        with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
            ensembl, uniprot, constraint = await _enrich_candidate_genes([v1, v2], [0, 1], config, {})
        assert "TP53" in ensembl
        assert "BRCA1" in ensembl
        assert ensembl["BRCA1"] == {}


@pytest.mark.l2
@pytest.mark.mock
class TestEnrichVariantFrequencies:
    """TP-43~48: _enrich_variant_frequencies async gnomAD + MyVariant."""

    @pytest.mark.asyncio
    async def test_all_candidates_have_af(self):
        """TP-43: All candidates already have gnomAD AF → skips queries."""
        from gpa_two_phase import _enrich_variant_frequencies
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", gnomad_af=0.001)
        # Should complete without any API calls
        await _enrich_variant_frequencies([v], [0], config)
        assert v.gnomad_af == 0.001

    @pytest.mark.asyncio
    async def test_myvariant_batch_success(self):
        """TP-44: MyVariant.info batch fills AF for candidates."""
        from gpa_two_phase import _enrich_variant_frequencies
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", gnomad_af=None)

        mock_mv_stats = {"gnomad_filled": 1, "clinvar_filled": 0, "cadd_filled": 0, "not_found": 0, "errors": 0}

        with patch("dgra_myvariant.query_myvariant_batch", new_callable=AsyncMock, return_value={}):
            with patch("dgra_myvariant.apply_myvariant_results", return_value=mock_mv_stats):
                await _enrich_variant_frequencies([v], [0], config)

    @pytest.mark.asyncio
    async def test_myvariant_failure_gnomad_fallback(self):
        """TP-45: MyVariant fails → gnomAD GraphQL fallback fills AF."""
        from gpa_two_phase import _enrich_variant_frequencies
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", gnomad_af=None)

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.query_gnomad_variant = AsyncMock(return_value={
            "af": 0.0001, "af_populations": {}, "status": "SUCCESS", "source": "gnomad"
        })

        with patch("dgra_myvariant.query_myvariant_batch", side_effect=Exception("network error")):
            with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
                await _enrich_variant_frequencies([v], [0], config)
        assert v.gnomad_af == 0.0001
        assert v.gnomad_status == "SUCCESS"

    @pytest.mark.asyncio
    async def test_gnomad_not_captured(self):
        """TP-46: gnomAD returns NOT_CAPTURED → variant keeps None AF."""
        from gpa_two_phase import _enrich_variant_frequencies
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", gnomad_af=None)

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.query_gnomad_variant = AsyncMock(return_value={
            "af": None, "status": "NOT_CAPTURED", "source": "gnomad"
        })

        with patch("dgra_myvariant.query_myvariant_batch", side_effect=Exception("network error")):
            with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
                await _enrich_variant_frequencies([v], [0], config)
        assert v.gnomad_af is None

    @pytest.mark.asyncio
    async def test_gnomad_api_failed(self):
        """TP-47: gnomAD API_FAILED → error handled gracefully."""
        from gpa_two_phase import _enrich_variant_frequencies
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", gnomad_af=None)

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.query_gnomad_variant = AsyncMock(return_value={
            "status": "API_FAILED", "error": "timeout", "source": "failed"
        })

        with patch("dgra_myvariant.query_myvariant_batch", side_effect=Exception("network error")):
            with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
                await _enrich_variant_frequencies([v], [0], config)
        # Should not crash; AF remains None
        assert v.gnomad_af is None

    @pytest.mark.asyncio
    async def test_clinvar_ncbi_query(self):
        """TP-48: NCBI ClinVar direct query for variants still unknown after MyVariant."""
        from gpa_two_phase import _enrich_variant_frequencies
        from dgra_core import GPAConfig, _UNKNOWN
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", gnomad_af=None, clinvar=_UNKNOWN, consequence="missense_variant")

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.query_gnomad_variant = AsyncMock(return_value={
            "af": None, "status": "NOT_CAPTURED", "source": "gnomad"
        })
        mock_client.query_ncbi_clinvar = AsyncMock(return_value={
            "clinical_significance": "Likely_pathogenic",
            "review_status": "criteria provided, single submitter",
        })

        mock_mv_stats = {"gnomad_filled": 0, "clinvar_filled": 0, "cadd_filled": 0, "not_found": 0, "errors": 0}

        with patch("dgra_myvariant.query_myvariant_batch", new_callable=AsyncMock, return_value={}):
            with patch("dgra_myvariant.apply_myvariant_results", return_value=mock_mv_stats):
                with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
                    await _enrich_variant_frequencies([v], [0], config)
        # ClinVar should be updated by NCBI query
        assert v.clinvar == "Likely_pathogenic"


@pytest.mark.l2
@pytest.mark.mock
class TestEnrichSpliceAI:
    """TP-49~52: _enrich_spliceai async enrichment."""

    @pytest.mark.asyncio
    async def test_no_splice_candidates(self):
        """TP-49: No splice-relevant candidates → skips query."""
        from gpa_two_phase import _enrich_spliceai
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", consequence="missense_variant")
        await _enrich_spliceai([v], [0], config)
        assert not hasattr(v, "spliceai_result") or v.spliceai_result is None

    @pytest.mark.asyncio
    async def test_splice_candidates_queried(self):
        """TP-50: Splice candidates trigger batch query and attach results."""
        from gpa_two_phase import _enrich_spliceai
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", consequence="splice_donor_variant")

        with patch("dgra_splice_predictor.query_spliceai_batch", new_callable=AsyncMock, return_value={
            "1:100000:A:G": {"delta_score": 0.85, "predicted_impact": "HIGH"}
        }):
            await _enrich_spliceai([v], [0], config)
        assert v.spliceai_result["delta_score"] == 0.85

    @pytest.mark.asyncio
    async def test_splice_candidate_not_in_db(self):
        """TP-51: Splice candidate not in SpliceAI DB → marked not_in_db."""
        from gpa_two_phase import _enrich_spliceai
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", consequence="splice_donor_variant")

        with patch("dgra_splice_predictor.query_spliceai_batch", new_callable=AsyncMock, return_value={}):
            await _enrich_spliceai([v], [0], config)
        assert v.spliceai_result == {"source": "not_in_db", "delta_score": None, "predicted_impact": None}

    @pytest.mark.asyncio
    async def test_empty_candidate_indices(self):
        """TP-52: Empty candidate_indices → returns immediately."""
        from gpa_two_phase import _enrich_spliceai
        from dgra_core import GPAConfig
        config = GPAConfig(tissue_profile="general")
        v = make_variant(gene="TP53", consequence="splice_donor_variant")
        await _enrich_spliceai([v], [], config)
        assert not hasattr(v, "spliceai_result") or v.spliceai_result is None


@pytest.mark.l2
@pytest.mark.mock
class TestEnrichPhenotype:
    """TP-53~57: _enrich_phenotype async phenotype matching."""

    @pytest.mark.asyncio
    async def test_no_user_phenotypes(self):
        """TP-53: No user_phenotypes → returns immediately."""
        from gpa_two_phase import _enrich_phenotype
        v = make_variant(gene="TP53")
        await _enrich_phenotype([v], [0], None)
        assert v.phenotype_match_score is None

    @pytest.mark.asyncio
    async def test_empty_candidates(self):
        """TP-54: Empty candidate_indices → returns immediately."""
        from gpa_two_phase import _enrich_phenotype
        v = make_variant(gene="TP53")
        await _enrich_phenotype([v], [], "肌无力")
        assert v.phenotype_match_score is None

    @pytest.mark.asyncio
    async def test_with_local_data(self):
        """TP-55: Genes with local phenotype data get matched."""
        from gpa_two_phase import _enrich_phenotype
        v = make_variant(gene="DMD")

        mock_matcher = MagicMock()
        mock_matcher._local_db = {"DMD": {"phenotypes": [{"name": "Duchenne muscular dystrophy"}]}}
        mock_matcher.match_batch = AsyncMock(return_value=[{
            "score": 0.95, "explanation": "High match", "confidence": "high",
            "matched_pairs": [], "known_phenotypes": [],
        }])

        with patch("gpa_phenotype_match.PhenotypeMatcher", return_value=mock_matcher):
            await _enrich_phenotype([v], [0], "肌无力")
        assert v.phenotype_match_score == 0.95

    @pytest.mark.asyncio
    async def test_without_local_data(self):
        """TP-56: Genes without local data get zero score."""
        from gpa_two_phase import _enrich_phenotype
        v = make_variant(gene="FAKEGENE")

        mock_matcher = MagicMock()
        mock_matcher._local_db = {"OTHER": {"phenotypes": [{"name": "something"}]}}
        mock_matcher.match_batch = AsyncMock(return_value=[])

        with patch("gpa_phenotype_match.PhenotypeMatcher", return_value=mock_matcher):
            await _enrich_phenotype([v], [0], "肌无力")
        assert v.phenotype_match_score == 0.0
        assert "No known phenotypes" in v.phenotype_match_explanation

    @pytest.mark.asyncio
    async def test_matcher_exception_fallback(self):
        """TP-57: PhenotypeMatcher exception → handled gracefully."""
        from gpa_two_phase import _enrich_phenotype
        v = make_variant(gene="DMD")

        with patch("gpa_phenotype_match.PhenotypeMatcher", side_effect=Exception("DB error")):
            await _enrich_phenotype([v], [0], "肌无力")
        # Should not crash; phenotype fields may remain unset


@pytest.mark.l2
@pytest.mark.mock
class TestRunTwoPhasePipeline:
    """TP-58~66: run_two_phase_pipeline and _run_two_phase_pipeline_impl."""

    @pytest.mark.asyncio
    async def test_empty_input(self):
        """TP-58: Empty variants_data → no candidates, quick return."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general"}
        config.to_global.return_value = MagicMock()

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_report.generate_tier_report", return_value="# Report"):
                with patch("gpa_report.generate_json_report", return_value={}):
                    with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                        with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 0, "by_flag": {}}):
                            result = await run_two_phase_pipeline([], config=config)
        assert "summary" in result

    @pytest.mark.asyncio
    async def test_all_tier3_no_candidates(self):
        """TP-59: All variants pre-filtered to Tier 3 → Phase 2 skipped."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general"}
        config.to_global.return_value = MagicMock()
        config.target_population = None

        variants_data = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
             "GENE": "TEST", "IMPACT": "LOW", "Consequence": "synonymous_variant",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": "0.5"},
        ]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_report.generate_tier_report", return_value="# Report"):
                with patch("gpa_report.generate_json_report", return_value={}):
                    with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                        with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 1, "by_flag": {}}):
                            result = await run_two_phase_pipeline(variants_data, config=config)
        assert result["summary"]["tier3_variant_count"] == 1

    @pytest.mark.asyncio
    async def test_normal_flow_with_candidates(self):
        """TP-60: Normal flow with Tier 1/2 candidates and mocked Phase 2."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general", "tier_rules": {}, "special_gene_lists": {}, "tissue_genes": set()}
        config.to_global.return_value = MagicMock()
        config.target_population = None

        variants_data = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
             "GENE": "TP53", "IMPACT": "HIGH", "Consequence": "stop_gained",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": ""},
        ]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_tier_classifier.classify_variant_tier", return_value=(1, "test", ["action"])):
                with patch("gpa_report.generate_tier_report", return_value="# Report"):
                    with patch("gpa_report.generate_json_report", return_value={}):
                        with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                            with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 1, "by_flag": {}}):
                                result = await run_two_phase_pipeline(variants_data, config=config)
        assert result["summary"]["tier1_variant_count"] == 1
        assert result["meta"]["total_variants"] == 1

    @pytest.mark.asyncio
    async def test_exception_wrapper_returns_error_dict(self):
        """TP-61: Exception in pipeline → wrapper returns error dict, does not raise."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general"}
        config.to_global.return_value = MagicMock()

        variants_data = [{"CHROM": "1", "POS": 100}]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            # Force _parse_variants_phase1 to raise by passing bad data
            result = await run_two_phase_pipeline(variants_data, config=config)
        assert "error" in result
        assert result["summary"]["tier1_variant_count"] == 0

    @pytest.mark.asyncio
    async def test_preflight_offline_switch(self):
        """TP-62: Preflight suggests offline → config switches to offline_mode."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = False
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = False
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general"}
        config.to_global.return_value = MagicMock()

        variants_data = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
             "GENE": "TEST", "IMPACT": "LOW", "Consequence": "synonymous_variant",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": "0.5"},
        ]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_preflight.suggest_action", return_value="offline"):
                result = await run_two_phase_pipeline(variants_data, config=config)
        assert config.offline_mode is True

    @pytest.mark.asyncio
    async def test_max_candidates_warning(self):
        """TP-63: Candidate count > max_candidates prints warning but proceeds."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general", "tier_rules": {}, "special_gene_lists": {}, "tissue_genes": set()}
        config.to_global.return_value = MagicMock()
        config.target_population = None

        # 200 HIGH impact variants → exceeds default max_candidates=150
        variants_data = [
            {"CHROM": "1", "POS": i, "REF": "A", "ALT": "G",
             "GENE": f"GENE{i}", "IMPACT": "HIGH", "Consequence": "stop_gained",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": ""}
            for i in range(200)
        ]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_tier_classifier.classify_variant_tier", return_value=(1, "test", ["action"])):
                with patch("gpa_report.generate_tier_report", return_value="# Report"):
                    with patch("gpa_report.generate_json_report", return_value={}):
                        with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                            with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 200, "by_flag": {}}):
                                result = await run_two_phase_pipeline(variants_data, config=config, max_candidates=150)
        assert result["summary"]["total_variants"] == 200
        assert result["summary"]["tier1_variant_count"] == 200

    @pytest.mark.asyncio
    async def test_phase2_enrichment_online(self):
        """TP-64: Online Phase 2 calls API enrichment for candidates."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = False
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general", "tier_rules": {}, "special_gene_lists": {}, "tissue_genes": set()}
        config.to_global.return_value = MagicMock()
        config.target_population = None
        config.spliceai_concurrency = 5

        variants_data = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
             "GENE": "TP53", "IMPACT": "HIGH", "Consequence": "stop_gained",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": ""},
        ]

        mock_client = AsyncMock()
        mock_client.__aenter__ = AsyncMock(return_value=mock_client)
        mock_client.__aexit__ = AsyncMock(return_value=False)
        mock_client.batch_query_genes = AsyncMock(return_value={"TP53": {"status": "SUCCESS"}})

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_two_phase.DGRAAPIClient", return_value=mock_client):
                with patch("gpa_tier_classifier.classify_variant_tier", return_value=(1, "test", ["action"])):
                    with patch("gpa_report.generate_tier_report", return_value="# Report"):
                        with patch("gpa_report.generate_json_report", return_value={}):
                            with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                                with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 1, "by_flag": {}}):
                                    result = await run_two_phase_pipeline(variants_data, config=config)
        assert result["summary"]["tier1_variant_count"] == 1

    @pytest.mark.asyncio
    async def test_large_input_all_tier3(self):
        """TP-65: Large input (>5000) with all Tier 3 → efficient handling."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general"}
        config.to_global.return_value = MagicMock()

        # 1000 common LOW variants
        variants_data = [
            {"CHROM": "1", "POS": i, "REF": "A", "ALT": "G",
             "GENE": "GENE", "IMPACT": "LOW", "Consequence": "synonymous_variant",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": "0.5"}
            for i in range(1000)
        ]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_report.generate_tier_report", return_value="# Report"):
                with patch("gpa_report.generate_json_report", return_value={}):
                    with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                        with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 1000, "by_flag": {}}):
                            result = await run_two_phase_pipeline(variants_data, config=config)
        assert result["summary"]["tier3_variant_count"] == 1000

    @pytest.mark.asyncio
    async def test_user_phenotypes_passed(self):
        """TP-66: user_phenotypes passed through to phenotype enrichment."""
        from gpa_two_phase import run_two_phase_pipeline
        report = MagicMock()
        report.is_ready.return_value = True
        route_map = MagicMock()
        route_map.to_markdown.return_value = ""

        config = MagicMock()
        config.offline_mode = True
        config.tissue_profile = "general"
        config.get_tissue_profile.return_value = {"display_name": "general", "tier_rules": {}, "special_gene_lists": {}, "tissue_genes": set()}
        config.to_global.return_value = MagicMock()
        config.target_population = None

        variants_data = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G",
             "GENE": "TP53", "IMPACT": "HIGH", "Consequence": "stop_gained",
             "CLIN_SIG": "", "DP": 30, "GQ": 99, "GT": "0/1", "VAF": "0.45", "gnomAD_AF": ""},
        ]

        with patch("gpa_preflight.run_preflight_check", return_value=(report, route_map)):
            with patch("gpa_tier_classifier.classify_variant_tier", return_value=(1, "test", ["action"])):
                with patch("gpa_report.generate_tier_report", return_value="# Report"):
                    with patch("gpa_report.generate_json_report", return_value={}):
                        with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                            with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 1, "by_flag": {}}):
                                with patch("gpa_two_phase._enrich_phenotype", new_callable=AsyncMock) as mock_pheno:
                                    result = await run_two_phase_pipeline(variants_data, config=config, user_phenotypes="肌无力")
        mock_pheno.assert_awaited_once()


@pytest.mark.l2
@pytest.mark.mock
class TestBuildOutput:
    """TP-67~70: _build_output report generation."""

    def test_build_output_basic(self):
        """TP-67: _build_output produces correct summary counts."""
        from gpa_two_phase import _build_output
        v1 = make_variant(gene="TP53", impact="HIGH")
        v1.tier = 1
        v2 = make_variant(gene="BRCA1", impact="MODERATE")
        v2.tier = 2
        v3 = make_variant(gene="TEST", impact="LOW")
        v3.tier = 3

        with patch("gpa_report.generate_tier_report", return_value="# Report"):
            with patch("gpa_report.generate_json_report", return_value={}):
                with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                    with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 3, "by_flag": {}}):
                        result = _build_output([v1, v2, v3], [0, 1], {}, MagicMock(), "general", "test")

        assert result["summary"]["tier1_variant_count"] == 1
        assert result["summary"]["tier2_variant_count"] == 1
        assert result["summary"]["tier3_variant_count"] == 1
        assert result["summary"]["tier1_gene_count"] == 1
        assert result["meta"]["profile_display_name"] == "general"

    def test_build_output_no_enrichment(self):
        """TP-68: No enriched indices → all variants marked preliminary."""
        from gpa_two_phase import _build_output
        v = make_variant(gene="TEST", impact="LOW")
        v.tier = 3
        v.tier_reason = "low impact"

        with patch("gpa_report.generate_tier_report", return_value="# Report"):
            with patch("gpa_report.generate_json_report", return_value={}):
                with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                    with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 1, "by_flag": {}}):
                        result = _build_output([v], [], {}, MagicMock(), "general", "test")

        assert result["summary"]["tier3_variant_count"] == 1
        assert "[PRELIMINARY]" in result["tier3_variants"][0]["tier_reason"]

    def test_build_output_multi_hit(self):
        """TP-69: Multi-hit genes included in output."""
        from gpa_two_phase import _build_output
        v1 = make_variant(gene="TP53", impact="HIGH")
        v1.tier = 1
        v2 = make_variant(gene="TP53", impact="MODERATE", pos=100001)
        v2.tier = 1

        with patch("gpa_report.generate_tier_report", return_value="# Report"):
            with patch("gpa_report.generate_json_report", return_value={}):
                with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[{"gene": "TP53", "count": 2}]):
                    with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 2, "by_flag": {}}):
                        result = _build_output([v1, v2], [0, 1], {}, MagicMock(), "general", "test")

        assert "TP53" in result["summary"]["multi_hit_genes"]

    def test_build_output_empty(self):
        """TP-70: Empty variant list → zero counts."""
        from gpa_two_phase import _build_output
        with patch("gpa_report.generate_tier_report", return_value="# Report"):
            with patch("gpa_report.generate_json_report", return_value={}):
                with patch("gpa_multi_hit.detect_multi_hit_genes", return_value=[]):
                    with patch("gpa_qc._run_qc_checks", return_value={"flagged": 0, "total": 0, "by_flag": {}}):
                        result = _build_output([], [], {}, MagicMock(), "general", "test")
        assert result["summary"]["total_variants"] == 0


@pytest.mark.l2
@pytest.mark.mock
class TestVariantNeedsChecksExtended:
    """TP-71~76: Extended consequence-based need checks."""

    def test_needs_clinvar_frameshift(self):
        """TP-71: Frameshift → needs ClinVar."""
        from gpa_two_phase import _variant_needs_clinvar
        v = make_variant(consequence="frameshift_variant")
        assert _variant_needs_clinvar(v) is True

    def test_needs_clinvar_intron(self):
        """TP-72: Intron variant → does NOT need ClinVar."""
        from gpa_two_phase import _variant_needs_clinvar
        v = make_variant(consequence="intron_variant")
        assert _variant_needs_clinvar(v) is False

    def test_needs_clinvar_comma_separated(self):
        """TP-73: Comma-separated consequences with relevant term → True."""
        from gpa_two_phase import _variant_needs_clinvar
        v = make_variant(consequence="intron_variant,missense_variant")
        assert _variant_needs_clinvar(v) is True

    def test_needs_spliceai_intron_near_exon(self):
        """TP-74: intron_variant with SPLICE_DIST <= 50 → True."""
        from gpa_two_phase import _variant_needs_spliceai
        v = make_variant(consequence="intron_variant")
        v.vcf_info = {"SPLICE_DIST": 30}
        assert _variant_needs_spliceai(v) is True

    def test_needs_spliceai_intron_far(self):
        """TP-75: intron_variant with SPLICE_DIST > 50 → False."""
        from gpa_two_phase import _variant_needs_spliceai
        v = make_variant(consequence="intron_variant")
        v.vcf_info = {"SPLICE_DIST": 100}
        assert _variant_needs_spliceai(v) is False

    def test_needs_spliceai_combined_consequences(self):
        """TP-76: Combined consequences with splice term → True."""
        from gpa_two_phase import _variant_needs_spliceai
        v = make_variant(consequence="missense_variant,splice_region_variant")
        assert _variant_needs_spliceai(v) is True
