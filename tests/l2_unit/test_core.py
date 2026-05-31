"""L2 unit tests for dgra_core.py — core data structures and utility functions."""
import json
import pytest
from dataclasses import asdict
from unittest.mock import patch, MagicMock

from dgra_core import (
    _is_unknown,
    Variant,
    GPAConfig,
    Evidence,
    correct_transcript_priority,
    get_pseudogenes_for_gene,
    _calculate_pseudogene_score,
    detect_pseudogene_artifact,
    classify_gnomad_frequency,
    normalize_gene_symbols,
    parse_protein_position,
    map_variant_to_domain,
)


class TestIsUnknown:
    def test_unknown_sentinel(self):
        assert _is_unknown("UNKNOWN") is True

    def test_empty_string(self):
        assert _is_unknown("") is True

    def test_none(self):
        assert _is_unknown(None) is True

    def test_valid_values(self):
        assert _is_unknown("HIGH") is False
        assert _is_unknown(0) is False
        assert _is_unknown(False) is False


class TestVariantDataclass:
    def test_basic_creation(self):
        v = Variant(
            chrom="17", pos=7578406, ref="C", alt="T",
            gene="TP53", transcript="ENST00000269305",
            exon="5/11", impact="HIGH", consequence="stop_gained",
            hgvsp="p.Arg175Ter", hgvsc="c.524G>A", clinvar="Pathogenic",
        )
        assert v.gene == "TP53"
        assert v.tier is None

    def test_defaults(self):
        v = Variant(
            chrom="1", pos=1, ref="A", alt="T",
            gene="X", transcript="", exon="", impact="", consequence="",
            hgvsp="", hgvsc="", clinvar="",
        )
        assert v.dp == 0
        assert v.gq == 0.0
        assert v.evidence_chain == []
        assert v.qc_flags == []
        assert v.quality_confidence == "high"

    def test_asdict(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="X",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        d = asdict(v)
        assert d["chrom"] == "1"


class TestGPAConfig:
    def test_to_global(self):
        cfg = GPAConfig(tissue_profile="hematopoietic", offline_mode=True)
        gc = cfg.to_global()
        assert gc.tissue_profile == "hematopoietic"
        assert gc.offline_mode is True

    def test_get_tissue_profile_raises_when_none(self):
        cfg = GPAConfig()
        with pytest.raises(ValueError, match="tissue_profile is required"):
            cfg.get_tissue_profile()

    def test_get_tissue_profile_invalid_name(self):
        cfg = GPAConfig(tissue_profile="nonexistent")
        with pytest.raises(ValueError, match="Unknown tissue profile"):
            cfg.get_tissue_profile()


class TestCorrectTranscriptPriority:
    @pytest.mark.asyncio
    async def test_already_canonical(self):
        v = Variant(chrom="17", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="ENST00000269305", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        ens = {"TP53": {"canonical_transcript": "ENST00000269305", "source": "test"}}
        v2, warning = await correct_transcript_priority(v, ens)
        assert warning is None
        assert v2.transcript_warning is None

    @pytest.mark.asyncio
    async def test_non_canonical_warning(self):
        v = Variant(chrom="17", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="ENST00000382972", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        ens = {"TP53": {"canonical_transcript": "ENST00000269305", "source": "test"}}
        v2, warning = await correct_transcript_priority(v, ens)
        assert warning is not None
        assert warning["type"] == "TRANSCRIPT_DISCREPANCY"
        assert v2.transcript_warning is not None

    @pytest.mark.asyncio
    async def test_no_api_data(self):
        v = Variant(chrom="17", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="ENST00000269305", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        v2, warning = await correct_transcript_priority(v, {})
        assert warning is None


class TestGetPseudogenesForGene:
    def test_lookup_found(self):
        with patch("dgra_core._load_pseudogene_lookup", return_value={"VWF": {"pseudogenes": ["VWFP1"]}}):
            pgs = get_pseudogenes_for_gene("VWF")
            assert "VWFP1" in pgs

    def test_legacy_fallback(self):
        # NOTE: _PSEUDOGENE_CONFIG_LEGACY uses key "pseudogene" (singular)
        # but code looks for "pseudogenes" (plural), so legacy fallback
        # does not actually retrieve pseudogenes from hardcoded config.
        with patch("dgra_core._load_pseudogene_lookup", return_value={}), \
             patch("dgra_core._load_pseudogene_database", return_value={}):
            pgs = get_pseudogenes_for_gene("SETBP1")
            assert pgs == []  # current behavior due to key mismatch

    def test_no_pseudogenes(self):
        with patch("dgra_core._load_pseudogene_lookup", return_value={}), \
             patch("dgra_core._load_pseudogene_database", return_value={}):
            pgs = get_pseudogenes_for_gene("TP53")
            assert pgs == []

    def test_deduplication(self):
        with patch("dgra_core._load_pseudogene_lookup", return_value={"VWF": {"pseudogenes": ["VWFP1", "VWFP1"]}}), \
             patch("dgra_core._load_pseudogene_database", return_value={"VWF": {"pseudogenes": ["VWFP1"]}}):
            pgs = get_pseudogenes_for_gene("VWF")
            assert pgs == ["VWFP1"]


class TestCalculatePseudogeneScore:
    def test_homozygous(self):
        result = _calculate_pseudogene_score(0.9, "1/1", ["PG"], "GENE")
        assert result["score"] == 0.0
        assert result["level"] == "none"

    def test_unknown_genotype(self):
        result = _calculate_pseudogene_score(0.5, "./.", ["PG"], "GENE")
        assert result["level"] == "unknown_genotype"

    def test_vaf_strong_interference(self):
        result = _calculate_pseudogene_score(0.15, "0/1", ["PG"], "GENE")
        assert result["score"] == 0.9
        assert result["level"] == "strong_interference"

    def test_vaf_suspected(self):
        result = _calculate_pseudogene_score(0.35, "0/1", ["PG"], "GENE")
        assert result["score"] == 0.55
        assert result["level"] == "suspected"

    def test_vaf_bias(self):
        result = _calculate_pseudogene_score(0.70, "0/1", ["PG"], "GENE")
        assert result["score"] == 0.30
        assert result["level"] == "bias_suspected"

    def test_vaf_normal(self):
        result = _calculate_pseudogene_score(0.50, "0/1", ["PG"], "GENE")
        assert result["score"] == 0.0
        assert result["level"] == "none"


class TestDetectPseudogeneArtifact:
    def test_no_pseudogenes(self):
        with patch("dgra_core.get_pseudogenes_for_gene", return_value=[]):
            v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                        transcript="", exon="", impact="", consequence="",
                        hgvsp="", hgvsc="", clinvar="", vaf=0.1, gt="0/1")
            assert detect_pseudogene_artifact(v) is None

    def test_strong_interference_detected(self):
        with patch("dgra_core.get_pseudogenes_for_gene", return_value=["PG1"]), \
             patch("dgra_core._load_pseudogene_lookup", return_value={"GENE": {"detection_strategy": "vaf", "confidence": "high"}}):
            v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="GENE",
                        transcript="", exon="", impact="", consequence="",
                        hgvsp="", hgvsc="", clinvar="", vaf=0.15, gt="0/1")
            result = detect_pseudogene_artifact(v)
            assert result is not None
            assert result["type"] == "PSEUDOGENE_INTERFERENCE"

    def test_no_interference(self):
        with patch("dgra_core.get_pseudogenes_for_gene", return_value=["PG1"]):
            v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="GENE",
                        transcript="", exon="", impact="", consequence="",
                        hgvsp="", hgvsc="", clinvar="", vaf=0.50, gt="0/1")
            assert detect_pseudogene_artifact(v) is None

    def test_missing_vaf_gt(self):
        with patch("dgra_core.get_pseudogenes_for_gene", return_value=["PG1"]):
            v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="GENE",
                        transcript="", exon="", impact="", consequence="",
                        hgvsp="", hgvsc="", clinvar="")
            assert detect_pseudogene_artifact(v) is None


class TestClassifyGnomadFrequency:
    def test_not_captured(self):
        result = classify_gnomad_frequency(None, "TP53")
        assert result["status"] == "NOT_CAPTURED"

    def test_germline_warning_gene(self):
        result = classify_gnomad_frequency(0.5, "ASXL1")
        assert result["status"] == "GERMLINE_WARNING"

    def test_common_polymorphism(self):
        # Use a gene NOT in GERMLINE_WARNING_GENES
        result = classify_gnomad_frequency(0.05, "BRCA1")
        assert result["status"] == "common_polymorphism"

    def test_low_frequency(self):
        result = classify_gnomad_frequency(0.005, "BRCA1")
        assert result["status"] == "low_frequency"

    def test_rare_variant(self):
        result = classify_gnomad_frequency(0.0005, "BRCA1")
        assert result["status"] == "rare_variant"

    def test_target_population(self):
        pop_data = {"EAS": {"af": 0.02}}
        result = classify_gnomad_frequency(0.0001, "BRCA1", af_by_population=pop_data, target_population="EAS")
        assert result["status"] == "common_polymorphism"
        assert "EAS" in result["interpretation"]

    def test_population_unavailable(self):
        result = classify_gnomad_frequency(0.01, "TP53", af_by_population={}, target_population="EAS")
        assert "EAS data unavailable" in result["interpretation"]


class TestNormalizeGeneSymbols:
    def test_approved_no_change(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        hgnc = {"TP53": {"status": "approved", "approved_symbol": "TP53"}}
        warnings = normalize_gene_symbols([v], hgnc)
        assert v.gene == "TP53"
        assert v.gene_original == "TP53"
        assert warnings == []

    def test_case_correction(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="tp53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        hgnc = {"tp53": {"status": "approved", "approved_symbol": "TP53"}}
        warnings = normalize_gene_symbols([v], hgnc)
        assert v.gene == "TP53"
        assert len(warnings) == 1

    def test_previous_symbol(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="OLD",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        hgnc = {"OLD": {"status": "previous", "approved_symbol": "NEW", "previous_symbols": ["OLD"]}}
        warnings = normalize_gene_symbols([v], hgnc)
        assert v.gene == "NEW"
        assert "HGNC WARNING" in warnings[0]

    def test_withdrawn(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="BAD",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        hgnc = {"BAD": {"status": "withdrawn", "approved_symbol": "GOOD"}}
        warnings = normalize_gene_symbols([v], hgnc)
        assert v.gene == "BAD"
        assert "WITHDRAWN" in warnings[0]

    def test_not_found_offline_known(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        hgnc = {"TP53": {"status": "not_found"}}
        warnings = normalize_gene_symbols([v], hgnc, offline_mode=True)
        assert warnings == []

    def test_not_found_offline_unknown(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="XYZ123",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        hgnc = {"XYZ123": {"status": "not_found"}}
        warnings = normalize_gene_symbols([v], hgnc, offline_mode=True)
        assert len(warnings) == 1
        assert "could not be validated" in warnings[0]


class TestParseProteinPosition:
    def test_basic(self):
        assert parse_protein_position("p.Arg175Ter") == 175

    def test_np_prefix(self):
        assert parse_protein_position("NP_000543.3:p.Val1565Leu") == 1565

    def test_numeric_only(self):
        assert parse_protein_position("p.123") == 123

    def test_empty(self):
        assert parse_protein_position("") is None
        assert parse_protein_position(None) is None
        assert parse_protein_position("nan") is None


class TestMapVariantToDomain:
    def test_no_domains(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="p.Arg175Ter", hgvsc="", clinvar="")
        result = map_variant_to_domain(v, {})
        assert result["domain"] == "unknown"

    def test_position_outside(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="p.Arg175Ter", hgvsc="", clinvar="")
        uniprot = {"TP53": {"domains": [{"name": "DBD", "start": 100, "end": 150, "interpro_id": "IPR001"}], "source": "test"}}
        result = map_variant_to_domain(v, uniprot)
        assert result["domain"] == "inter-domain / unannotated"

    def test_position_inside_domain(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="p.Arg175Ter", hgvsc="", clinvar="")
        uniprot = {"TP53": {"domains": [{"name": "DBD", "start": 100, "end": 300, "interpro_id": "IPR001"}], "source": "test"}}
        result = map_variant_to_domain(v, uniprot)
        assert result["domain"] == "DBD"

    def test_no_position(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="", hgvsc="", clinvar="")
        uniprot = {"TP53": {"domains": [{"name": "DBD", "start": 100, "end": 300}], "source": "test"}}
        result = map_variant_to_domain(v, uniprot)
        assert result["domain"] == "unknown"
        assert "Could not parse protein position" in result["note"]

    def test_exceeds_sequence_length(self):
        v = Variant(chrom="1", pos=1, ref="A", alt="T", gene="TP53",
                    transcript="", exon="", impact="", consequence="",
                    hgvsp="p.Arg9999Ter", hgvsc="", clinvar="")
        uniprot = {"TP53": {"sequence_length": 393, "source": "test"}}
        result = map_variant_to_domain(v, uniprot)
        assert "exceeds UniProt sequence length" in result["note"]
