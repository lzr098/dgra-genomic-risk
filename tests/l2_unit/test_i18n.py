"""L2 unit tests for gpa_i18n.py — consequence normalization, header translation."""
import pytest
from gpa_i18n import (
    normalize_consequence,
    infer_impact_from_consequence,
    translate_chinese_headers,
    is_chinese_header,
    normalize_clinvar,
    _translate_single_cn_header,
    EXONIC_FUNC_MAP,
    CN_EXONIC_FUNC_MAP,
)


class TestNormalizeConsequence:
    def test_empty_and_unknown(self):
        assert normalize_consequence("") == []
        assert normalize_consequence("UNKNOWN") == []
        assert normalize_consequence("N/A") == []
        assert normalize_consequence(".") == []
        assert normalize_consequence(None) == []

    def test_english_terms_pass_through(self):
        assert normalize_consequence("missense_variant") == ["missense_variant"]
        assert normalize_consequence("missense_variant,splice_region_variant") == [
            "missense_variant",
            "splice_region_variant",
        ]

    def test_chinese_to_english(self):
        assert normalize_consequence("错义变异") == ["missense_variant"]
        assert normalize_consequence("移码变异") == ["frameshift_variant"]

    def test_ampersand_delimiter(self):
        assert normalize_consequence("missense_variant&splice_region_variant") == [
            "missense_variant",
            "splice_region_variant",
        ]

    def test_mixed_chinese_english(self):
        result = normalize_consequence("错义变异,stop_gained")
        assert "missense_variant" in result
        assert "stop_gained" in result

    def test_whitespace_trimmed(self):
        assert normalize_consequence("  missense_variant  ") == ["missense_variant"]


class TestInferImpact:
    def test_high_impact(self):
        assert infer_impact_from_consequence("stop_gained") == "HIGH"
        assert infer_impact_from_consequence("frameshift_variant") == "HIGH"
        assert infer_impact_from_consequence("splice_acceptor_variant") == "HIGH"

    def test_moderate_impact(self):
        assert infer_impact_from_consequence("missense_variant") == "MODERATE"
        assert infer_impact_from_consequence("inframe_deletion") == "MODERATE"

    def test_low_impact(self):
        assert infer_impact_from_consequence("synonymous_variant") == "LOW"
        assert infer_impact_from_consequence("intron_variant") == "LOW"

    def test_modifier(self):
        assert infer_impact_from_consequence("regulatory_region_variant") == "MODIFIER"
        assert infer_impact_from_consequence("intergenic_variant") == "MODIFIER"

    def test_chinese_input(self):
        assert infer_impact_from_consequence("错义变异") == "MODERATE"
        assert infer_impact_from_consequence("移码变异") == "HIGH"

    def test_empty_input(self):
        assert infer_impact_from_consequence("") == ""
        assert infer_impact_from_consequence("UNKNOWN") == ""

    def test_combined_terms(self):
        # If any term is HIGH, overall is HIGH
        assert infer_impact_from_consequence("missense_variant,stop_gained") == "HIGH"


class TestTranslateHeaders:
    def test_exact_match(self):
        assert translate_chinese_headers(["基因"]) == ["Gene"]
        assert translate_chinese_headers(["变异后果"]) == ["Consequence"]

    def test_english_preserved(self):
        assert translate_chinese_headers(["HGVSp"]) == ["HGVSp"]

    def test_unknown_preserved(self):
        assert translate_chinese_headers(["FooBar"]) == ["FooBar"]

    def test_is_chinese_header_true(self):
        assert is_chinese_header(["基因", "位置"]) is True

    def test_is_chinese_header_false(self):
        assert is_chinese_header(["Gene", "POS"]) is False


class TestTranslateSingleHeader:
    def test_sample_column(self):
        assert _translate_single_cn_header("样本: P008") == "GT"

    def test_exact_mapping(self):
        assert _translate_single_cn_header("染色体") == "CHROM"
        assert _translate_single_cn_header("基因符号") == "GENE"

    def test_prefix_match(self):
        assert _translate_single_cn_header("质量值(QUAL)") == "_qual"

    def test_unmapped_pass_through(self):
        assert _translate_single_cn_header("UnknownColumn") == "UnknownColumn"


class TestNormalizeClinVar:
    def test_basic(self):
        assert normalize_clinvar("Pathogenic&Likely_pathogenic") == "Pathogenic/Likely pathogenic"

    def test_empty(self):
        assert normalize_clinvar("") == ""
        assert normalize_clinvar(None) == ""

    def test_underscore_replacement(self):
        assert normalize_clinvar("Likely_pathogenic") == "Likely pathogenic"


class TestExonicFuncMaps:
    def test_english_map(self):
        assert EXONIC_FUNC_MAP["nonsynonymous snv"] == "missense_variant"
        assert EXONIC_FUNC_MAP["stopgain"] == "stop_gained"

    def test_chinese_map(self):
        assert CN_EXONIC_FUNC_MAP["错义变异"] == "missense_variant"
        assert CN_EXONIC_FUNC_MAP["移码变异"] == "frameshift_variant"
