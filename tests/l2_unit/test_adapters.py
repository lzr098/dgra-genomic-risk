"""L2 unit tests for dgra_adapters.py — VEP, ANNOVAR, SnpEff adapters."""
import pytest
from dgra_adapters import (
    VEPAdapter,
    ANNOVARAdapter,
    SnpEffAdapter,
    auto_detect_adapter,
    adapt_rows,
    DGRA_COLS,
)


class TestVEPAdapter:
    def test_direct_pass_through(self):
        adapter = VEPAdapter()
        row = {
            "CHROM": "17",
            "POS": "7578406",
            "REF": "C",
            "ALT": "T",
            "SYMBOL": "TP53",
            "Feature": "ENST00000269305",
            "EXON": "5/11",
            "IMPACT": "HIGH",
            "Consequence": "stop_gained",
            "HGVSp": "p.Arg175Ter",
            "HGVSc": "c.524G>A",
            "CLIN_SIG": "Pathogenic",
            "GT": "0/1",
            "DP": "30",
            "GQ": "99",
            "VAF": "0.45",
            "gnomAD_AF": "0.00001",
        }
        out = adapter.adapt(row)
        assert out["CHROM"] == "17"
        assert out["GENE"] == "TP53"
        assert out["IMPACT"] == "HIGH"
        assert out["HGVSp"] == "p.Arg175Ter"

    def test_clin_sig_normalization(self):
        adapter = VEPAdapter()
        row = {"CLIN_SIG": "Pathogenic&Likely_pathogenic", "Consequence": "missense_variant"}
        out = adapter.adapt(row)
        assert "/" in out["CLIN_SIG"]
        assert "&" not in out["CLIN_SIG"]

    def test_hgvsp_prefix_added(self):
        adapter = VEPAdapter()
        row = {"HGVSp": "Arg175Ter"}
        out = adapter.adapt(row)
        assert out["HGVSp"].startswith("p.")

    def test_parse_uploaded_variation(self):
        adapter = VEPAdapter()
        assert adapter._parse_uploaded_variation("1_12345_A/G") == ("1", "12345")
        assert adapter._parse_uploaded_variation("1-12345-A/G") == ("1", "12345")
        assert adapter._parse_uploaded_variation("rs123") is None

    def test_missing_columns_filled(self):
        adapter = VEPAdapter()
        out = adapter.adapt({"SYMBOL": "BRCA1"})
        for col in DGRA_COLS:
            assert col in out

    def test_chinese_impact_translation(self):
        adapter = VEPAdapter()
        out = adapter.adapt({"IMPACT": "高"})
        assert out["IMPACT"] == "HIGH"
        out2 = adapter.adapt({"IMPACT": "中等"})
        assert out2["IMPACT"] == "MODERATE"

    def test_chinese_consequence_normalization(self):
        adapter = VEPAdapter()
        out = adapter.adapt({"Consequence": "错义变异"})
        assert "missense_variant" in out["Consequence"]

    def test_sample_column_parsing(self):
        adapter = VEPAdapter()
        out = adapter.adapt({"GT": "GT:0/1 DP:48 GQ:99"})
        assert out["GT"] == "0/1"
        assert out["DP"] == "48"
        assert out["GQ"] == "99"

    def test_supports_headers(self):
        assert VEPAdapter.supports_headers(["Consequence", "IMPACT", "HGVSp"])
        assert not VEPAdapter.supports_headers(["Chr", "Start"])


class TestANNOVARAdapter:
    def test_basic_mapping(self):
        adapter = ANNOVARAdapter()
        row = {
            "Chr": "17",
            "Start": "7578406",
            "Ref": "C",
            "Alt": "T",
            "Gene.refGene": "TP53",
        }
        out = adapter.adapt(row)
        assert out["CHROM"] == "17"
        assert out["GENE"] == "TP53"

    def test_aachange_parsing(self):
        adapter = ANNOVARAdapter()
        row = {"AAChange.refGene": "NM_001:exon5:c.123A>G:p.Arg41Cys"}
        out = adapter.adapt(row)
        assert out["HGVSc"] == "c.123A>G"
        assert out["HGVSp"] == "p.Arg41Cys"

    def test_exonic_func_mapping(self):
        adapter = ANNOVARAdapter()
        out = adapter.adapt({"ExonicFunc.refGene": "nonsynonymous snv"})
        assert out["Consequence"] == "missense_variant"

    def test_func_fallback(self):
        adapter = ANNOVARAdapter()
        out = adapter.adapt({"Func.refGene": "splicing"})
        assert out["Consequence"] == "splice_region_variant"

    def test_impact_inference(self):
        adapter = ANNOVARAdapter()
        out = adapter.adapt({"ExonicFunc.refGene": "frameshift deletion"})
        assert out["IMPACT"] == "HIGH"

    def test_supports_headers(self):
        assert ANNOVARAdapter.supports_headers(["Chr", "Start", "Gene.refGene"])
        assert not ANNOVARAdapter.supports_headers(["Consequence", "IMPACT"])


class TestSnpEffAdapter:
    def test_structured_columns(self):
        adapter = SnpEffAdapter()
        row = {
            "ANN[0].GENE": "TP53",
            "ANN[0].IMPACT": "HIGH",
            "ANN[0].EFFECT": "stop_gained",
            "ANN[0].HGVS_C": "c.524G>A",
            "ANN[0].HGVS_P": "p.Arg175Ter",
        }
        out = adapter.adapt(row)
        assert out["GENE"] == "TP53"
        assert out["IMPACT"] == "HIGH"

    def test_raw_ann_parsing_format_b(self):
        adapter = SnpEffAdapter()
        row = {
            "ANN": "missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|5/11|c.524G>A|p.Arg175His",
        }
        out = adapter.adapt(row)
        assert out["GENE"] == "TP53"
        assert out["Consequence"] == "missense_variant"

    def test_raw_ann_parsing_format_a(self):
        adapter = SnpEffAdapter()
        row = {
            "ANN": "A|missense_variant|MODERATE|TP53|ENSG00000141510|transcript|ENST00000269305|protein_coding|5/11|c.524G>A|p.Arg175His",
        }
        out = adapter.adapt(row)
        assert out["GENE"] == "TP53"
        assert out["IMPACT"] == "MODERATE"

    def test_supports_headers(self):
        assert SnpEffAdapter.supports_headers(["ANN[0].EFFECT", "ANN[0].GENE"])
        assert not SnpEffAdapter.supports_headers(["Chr", "Start"])

    def test_missing_columns_filled(self):
        adapter = SnpEffAdapter()
        out = adapter.adapt({"ANN[0].GENE": "BRCA1"})
        for col in DGRA_COLS:
            assert col in out


class TestAutoDetect:
    def test_detect_vep(self):
        adapter = auto_detect_adapter(["Consequence", "IMPACT", "HGVSp", "SYMBOL"])
        assert isinstance(adapter, VEPAdapter)

    def test_detect_annovar(self):
        adapter = auto_detect_adapter(["Chr", "Start", "Gene.refGene", "AAChange.refGene"])
        assert isinstance(adapter, ANNOVARAdapter)

    def test_detect_snpeff(self):
        adapter = auto_detect_adapter(["ANN[0].EFFECT", "ANN[0].GENE"])
        assert isinstance(adapter, SnpEffAdapter)

    def test_default_to_vep(self):
        adapter = auto_detect_adapter(["Foo", "Bar"])
        assert isinstance(adapter, VEPAdapter)

    def test_chinese_headers(self):
        adapter = auto_detect_adapter(["染色体", "位置", "基因符号"])
        # After translation these map to CHROM, POS, GENE — not enough for any adapter,
        # so it falls back to VEP
        assert isinstance(adapter, VEPAdapter)


class TestAdaptRows:
    def test_empty_rows(self):
        assert adapt_rows([]) == []

    def test_auto_detect_and_adapt(self):
        rows = [{"CHROM": "1", "POS": "123", "Consequence": "missense_variant"}]
        result = adapt_rows(rows)
        assert result[0]["CHROM"] == "1"
        assert result[0]["Consequence"] == "missense_variant"

    def test_explicit_adapter(self):
        rows = [{"Chr": "1", "Start": "123", "Gene.refGene": "TP53"}]
        result = adapt_rows(rows, adapter=ANNOVARAdapter())
        assert result[0]["GENE"] == "TP53"
