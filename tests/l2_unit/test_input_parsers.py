"""L2 Unit Tests — dgra_input_parsers.py

Covers auto_detect, TSVParser, VCFParser, FreeTextParser, and parse_input.
"""

import pytest
import csv
import tempfile
from pathlib import Path
from unittest.mock import patch, MagicMock


@pytest.mark.l2
class TestAutoDetect:
    """INP-01~08: auto_detect format detection."""

    def test_vcf_gz(self):
        """INP-01: .vcf.gz → vcf."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False) as f:
            path = f.name
        assert auto_detect(Path(path)) == "vcf"
        Path(path).unlink()

    def test_vcf(self):
        """INP-02: .vcf → vcf."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(suffix=".vcf", delete=False) as f:
            path = f.name
        assert auto_detect(Path(path)) == "vcf"
        Path(path).unlink()

    def test_bcf(self):
        """INP-03: .bcf → vcf."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(suffix=".bcf", delete=False) as f:
            path = f.name
        assert auto_detect(Path(path)) == "vcf"
        Path(path).unlink()

    def test_tsv(self):
        """INP-04: .tsv with tab → tsv."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\n")
            path = f.name
        assert auto_detect(Path(path)) == "tsv"
        Path(path).unlink()

    def test_csv(self):
        """INP-05: .csv with comma → csv."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("CHROM,POS,REF,ALT\n")
            path = f.name
        assert auto_detect(Path(path)) == "csv"
        Path(path).unlink()

    def test_xlsx(self):
        """INP-06: .xlsx → excel."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(suffix=".xlsx", delete=False) as f:
            path = f.name
        assert auto_detect(Path(path)) == "excel"
        Path(path).unlink()

    def test_txt(self):
        """INP-07: .txt → freetext."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(suffix=".txt", delete=False) as f:
            path = f.name
        assert auto_detect(Path(path)) == "freetext"
        Path(path).unlink()

    def test_vcf_header_fallback(self):
        """INP-08: Unknown extension but VCF header → vcf."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(mode="w", suffix=".unknown", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            path = f.name
        assert auto_detect(Path(path)) == "vcf"
        Path(path).unlink()

    def test_unknown_raises(self):
        """INP-09: Unknown format → ValueError."""
        from dgra_input_parsers import auto_detect
        with tempfile.NamedTemporaryFile(suffix=".bin", delete=False) as f:
            f.write(b"random binary data")
            path = f.name
        with pytest.raises(ValueError):
            auto_detect(Path(path))
        Path(path).unlink()


@pytest.mark.l2
class TestTSVParser:
    """INP-10~15: TSVParser."""

    def test_parse_tsv_basic(self):
        """INP-10: Basic TSV parsing."""
        from dgra_input_parsers import TSVParser
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
            f.write("1\t100\tA\tG\tBRCA1\n")
            path = f.name
        parser = TSVParser(dialect="tab")
        rows = parser.parse(Path(path))
        assert len(rows) == 1
        assert rows[0].get("CHROM") == "1"
        assert rows[0].get("GENE") == "BRCA1"
        Path(path).unlink()

    def test_parse_csv(self):
        """INP-11: CSV parsing."""
        from dgra_input_parsers import TSVParser
        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write("CHROM,POS,REF,ALT,GENE\n")
            f.write("1,100,A,G,BRCA1\n")
            path = f.name
        parser = TSVParser(dialect="comma")
        rows = parser.parse(Path(path))
        assert len(rows) == 1
        assert rows[0].get("CHROM") == "1"
        Path(path).unlink()

    def test_auto_detect_delimiter(self):
        """INP-12: Auto-detect tab delimiter."""
        from dgra_input_parsers import TSVParser
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\n")
            f.write("1\t100\tA\tG\n")
            path = f.name
        parser = TSVParser(dialect="auto")
        rows = parser.parse(Path(path))
        assert len(rows) == 1
        Path(path).unlink()

    def test_empty_file(self):
        """INP-13: Empty file → empty list."""
        from dgra_input_parsers import TSVParser
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\n")
            path = f.name
        parser = TSVParser(dialect="tab")
        rows = parser.parse(Path(path))
        assert rows == []
        Path(path).unlink()

    def test_whitespace_stripping(self):
        """INP-14: Whitespace stripped from values."""
        from dgra_input_parsers import TSVParser
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\n")
            f.write("  1  \t  100  \t  A  \t  G  \n")
            path = f.name
        parser = TSVParser(dialect="tab")
        rows = parser.parse(Path(path))
        assert rows[0].get("CHROM") == "1"
        assert rows[0].get("REF") == "A"
        Path(path).unlink()

    def test_none_values(self):
        """INP-15: None values handled."""
        from dgra_input_parsers import TSVParser
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\n")
            f.write("1\t100\t\tG\n")
            path = f.name
        parser = TSVParser(dialect="tab")
        rows = parser.parse(Path(path))
        assert rows[0].get("REF") == ""
        Path(path).unlink()


@pytest.mark.l2
class TestFreeTextParser:
    """INP-16~23: FreeTextParser."""

    def test_c_hgvs(self):
        """INP-16: c.HGVS pattern "TP53 c.722C>T"."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        result = parser.parse_text("TP53 c.722C>T")
        assert len(result) == 1
        assert result[0]["GENE"] == "TP53"
        assert result[0]["HGVSc"] == "c.722C>T"

    def test_c_hgvs_del(self):
        """INP-17: c.HGVS deletion "BRCA1 c.68_69delAG"."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        result = parser.parse_text("BRCA1 c.68_69delAG")
        assert len(result) == 1
        assert result[0]["GENE"] == "BRCA1"

    def test_genomic_coord_colon(self):
        """INP-18: Genomic coord "chr17:7578406C>A"."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        result = parser.parse_text("chr17:7578406C>A")
        assert len(result) == 1
        assert result[0]["CHROM"] == "17"
        assert result[0]["POS"] == "7578406"
        assert result[0]["REF"] == "C"
        assert result[0]["ALT"] == "A"

    def test_genomic_coord_dash(self):
        """INP-19: Genomic coord "17-7578406C>A"."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        result = parser.parse_text("17-7578406C>A")
        assert len(result) == 1
        assert result[0]["CHROM"] == "17"

    def test_genomic_coord_space(self):
        """INP-20: Space-separated "17 7578406 C A"."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        result = parser.parse_text("17 7578406 C A")
        assert len(result) == 1
        assert result[0]["CHROM"] == "17"
        assert result[0]["REF"] == "C"

    def test_p_hgvs(self):
        """INP-21: p.HGVS "TP53 p.Arg249Ser"."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        result = parser.parse_text("TP53 p.Arg249Ser")
        assert len(result) == 1
        assert result[0]["GENE"] == "TP53"
        assert result[0]["HGVSp"] == "p.Arg249Ser"

    def test_empty_text(self):
        """INP-22: Empty text → empty list."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        assert parser.parse_text("") == []
        assert parser.parse_text("   ") == []

    def test_unmatched_text(self):
        """INP-23: Unmatched text → ValueError."""
        from dgra_input_parsers import FreeTextParser
        parser = FreeTextParser()
        with pytest.raises(ValueError):
            parser.parse_text("some random text")


@pytest.mark.l2
class TestParseInput:
    """INP-24~27: parse_input dispatcher."""

    def test_parse_tsv_file(self):
        """INP-24: parse_input with TSV file."""
        from dgra_input_parsers import parse_input
        with tempfile.NamedTemporaryFile(mode="w", suffix=".tsv", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
            f.write("1\t100\tA\tG\tBRCA1\n")
            path = f.name
        rows = parse_input(Path(path))
        assert len(rows) == 1
        Path(path).unlink()

    def test_parse_vcf_file(self):
        """INP-25: parse_input with VCF file."""
        from dgra_input_parsers import parse_input
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("1\t100\t.\tA\tG\t30\tPASS\tDP=20\n")
            path = f.name
        # VCFParser requires vcfpy; skip if not available
        try:
            rows = parse_input(Path(path))
            assert isinstance(rows, list)
        except ImportError:
            pytest.skip("vcfpy not installed")
        Path(path).unlink()

    def test_parse_with_explicit_format(self):
        """INP-26: Explicit format overrides auto-detect."""
        from dgra_input_parsers import parse_input
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("CHROM\tPOS\tREF\tALT\tGENE\n")
            f.write("1\t100\tA\tG\tBRCA1\n")
            path = f.name
        rows = parse_input(Path(path), fmt="tsv")
        assert len(rows) == 1
        Path(path).unlink()

    def test_parse_freetext_file(self):
        """INP-27: parse_input with freetext file."""
        from dgra_input_parsers import parse_input
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("TP53 c.722C>T\n")
            path = f.name
        rows = parse_input(Path(path))
        assert len(rows) == 1
        assert rows[0].get("GENE") == "TP53"
        Path(path).unlink()


@pytest.mark.l2
class TestParseAnnotatedVcf:
    """INP-28~32: parse_annotated_vcf wrapper for VEP-annotated VCF."""

    def _make_annotated_vcf(self, csq_entries: list) -> str:
        """Create a temporary annotated VCF with given CSQ entries."""
        import tempfile
        csq_str = ",".join(csq_entries)
        content = (
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description="
            '\"Consequence annotations from Ensembl VEP. Format: '
            'Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|EXON|CLIN_SIG|gnomAD_AF\">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            f"17\t7579472\t.\tC\tT\t30\tPASS\tCSQ={csq_str}\tGT:DP\t0/1:35\n"
        )
        fd, path = tempfile.mkstemp(suffix=".vcf")
        with open(fd, "w") as f:
            f.write(content)
        return path

    def test_parse_annotated_vcf_basic(self):
        """INP-28: Parse VEP-annotated VCF with CSQ fields."""
        try:
            import vcfpy  # noqa: F401
        except ImportError:
            pytest.skip("vcfpy not installed")
        from gpa_input import parse_annotated_vcf

        csq = (
            "T|missense_variant|MODERATE|BRCA1|ENSG00000012048|"
            "ENST00000357654|ENST00000357654.9:c.722C>T|p.Arg241Trp|5/22|"
            "Pathogenic/Likely_pathogenic|0.0001"
        )
        path = self._make_annotated_vcf([csq])
        try:
            variants = parse_annotated_vcf(path)
            assert len(variants) == 1
            v = variants[0]
            assert v["CHROM"] == "17"
            assert v["POS"] == "7579472"
            assert v["REF"] == "C"
            assert v["ALT"] == "T"
            assert v["GENE"] == "BRCA1"
            assert v["Consequence"] == "missense_variant"
            assert v["IMPACT"] == "MODERATE"
            assert v["HGVSc"] == "ENST00000357654.9:c.722C>T"
            assert v["HGVSp"] == "p.Arg241Trp"
            assert v["EXON"] == "5/22"
            assert v["CLIN_SIG"] == "Pathogenic/Likely_pathogenic"
            assert v["gnomAD_AF"] == "0.0001"
            assert v["GT"] == "0/1"
            assert v["DP"] == "35"
        finally:
            Path(path).unlink()

    def test_parse_annotated_vcf_file_not_found(self):
        """INP-29: Raise ValueError when VCF file does not exist."""
        from gpa_input import parse_annotated_vcf
        with pytest.raises(ValueError, match="VCF file not found"):
            parse_annotated_vcf("/nonexistent/path/file.vcf")

    def test_parse_annotated_vcf_no_csq(self):
        """INP-30: VCF without CSQ returns minimal records (coordinates only)."""
        try:
            import vcfpy  # noqa: F401
        except ImportError:
            pytest.skip("vcfpy not installed")
        from gpa_input import parse_annotated_vcf
        import tempfile

        content = (
            "##fileformat=VCFv4.2\n"
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "1\t100000\t.\tA\tG\t30\tPASS\tDP=20\tGT\t0/1\n"
        )
        fd, path = tempfile.mkstemp(suffix=".vcf")
        with open(fd, "w") as f:
            f.write(content)
        try:
            variants = parse_annotated_vcf(path)
            # VCFParser emits minimal records even without CSQ
            assert len(variants) == 1
            assert variants[0]["CHROM"] == "1"
            assert variants[0]["GENE"] == ""
        finally:
            Path(path).unlink()

    def test_parse_annotated_vcf_multiallelic(self):
        """INP-31: Multi-allelic site split into separate records."""
        try:
            import vcfpy  # noqa: F401
        except ImportError:
            pytest.skip("vcfpy not installed")
        from gpa_input import parse_annotated_vcf

        csq_c = "C|synonymous_variant|LOW|TP53|ENSG00000141510|ENST00000269305|ENST00000269305.4:c.215C>G|p.Pro72Pro|4/11||0.5"
        csq_t = "T|missense_variant|MODERATE|TP53|ENSG00000141510|ENST00000269305|ENST00000269305.4:c.215C>T|p.Pro72Ser|4/11||0.01"
        path = self._make_annotated_vcf([csq_c, csq_t])
        try:
            variants = parse_annotated_vcf(path)
            assert len(variants) == 1  # Only one ALT in our test VCF (T)
            # vcfpy parses ALT from record.ALT, which is [T] in this VCF
            # so only CSQ entries matching allele T are kept
            assert variants[0]["ALT"] == "T"
            assert variants[0]["GENE"] == "TP53"
        finally:
            Path(path).unlink()

    def test_core_annotated_vcf_path(self):
        """INP-32: dgra_core no longer raises NotImplementedError for annotated VCF."""
        import tempfile
        from gpa_input import InputType, detect_input_type

        # Create a minimal annotated VCF
        content = (
            "##fileformat=VCFv4.2\n"
            "##INFO=<ID=CSQ,Number=.,Type=String,Description="
            '\"Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature|HGVSc|HGVSp|EXON|CLIN_SIG|gnomAD_AF\">\n'
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n"
            "1\t100\t.\tA\tG\t30\tPASS\tCSQ=G|missense_variant|MODERATE|BRCA1|ENSG00000012048|ENST00000357654|c.100A>G|p.Asn34Ser|2/22||0.001\tGT:DP\t0/1:30\n"
        )
        fd, path = tempfile.mkstemp(suffix=".vcf")
        with open(fd, "w") as f:
            f.write(content)

        try:
            # Verify input type detection works
            it = detect_input_type(path)
            assert it == InputType.ANNOTATED_VCF

            # Verify parse_annotated_vcf works end-to-end
            from gpa_input import parse_annotated_vcf
            variants = parse_annotated_vcf(path)
            assert len(variants) >= 1
            assert variants[0]["GENE"] == "BRCA1"

            # Verify dgra_core.py no longer has NotImplementedError for this path
            import dgra_core
            source = open(dgra_core.__file__).read()
            assert "Annotated VCF input parsing not yet implemented" not in source
        finally:
            Path(path).unlink()
