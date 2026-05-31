"""L2 Unit Tests — gpa_vcf_annotator.py

Covers VCFAnnotator initialization, VCF parsing, genome detection,
pre-filtering, annotator resolution, VEP API annotation, local VEP
annotation, shard management, checkpoint handling, and static helpers.
"""

import asyncio
import gzip
import json
import os
import tempfile
from pathlib import Path
from unittest.mock import AsyncMock, MagicMock, patch, mock_open

import pytest


# =============================================================================
# Helpers
# =============================================================================

def _make_vcf_content(lines, header_extra=None):
    """Create a temporary VCF file with given variant lines."""
    with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
        f.write("##fileformat=VCFv4.2\n")
        if header_extra:
            for h in header_extra:
                f.write(h + "\n")
        f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
        for line in lines:
            f.write(line + "\n")
        return f.name


@pytest.mark.l2
class TestVCFAnnotatorInit:
    """VAN-01~04: VCFAnnotator initialization."""

    def test_default_init(self):
        """VAN-01: Default parameters."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        assert va.annotator == "auto"
        assert va.genome == "auto"
        assert va.batch_size == 200
        assert va.max_concurrency == 5
        assert va.timeout == 30
        assert va.vep_params == {}
        assert va.interactive is True

    def test_custom_init(self):
        """VAN-02: Custom parameters."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(
            annotator="vep_api",
            genome="GRCh37",
            batch_size=500,
            max_concurrency=10,
            timeout=60,
            vep_params={"SIFT": "1"},
            interactive=False,
        )
        assert va.annotator == "vep_api"
        assert va.genome == "GRCh37"
        assert va.batch_size == 500
        assert va.max_concurrency == 10
        assert va.timeout == 60
        assert va.vep_params == {"SIFT": "1"}
        assert va.interactive is False

    def test_checkpoint_path(self):
        """VAN-03: Checkpoint path stored."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(checkpoint_path="/tmp/check.json")
        assert va.checkpoint_path == "/tmp/check.json"

    def test_shard_dir(self):
        """VAN-04: Shard dir stored."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(shard_dir="/tmp/shards", resume=True)
        assert va.shard_dir == "/tmp/shards"
        assert va.resume is True


@pytest.mark.l2
class TestGenomeDetection:
    """VAN-05~07: Genome build detection from VCF header."""

    def _make_vcf(self, header_lines):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            for h in header_lines:
                f.write(h + "\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("1\t100\t.\tA\tG\t30\tPASS\tDP=20\n")
            return f.name

    def test_detect_grch38(self):
        """VAN-05: GRCh38 in header."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf(["##reference=GRCh38"])
        va = VCFAnnotator()
        assert va._detect_genome(Path(path)) == "GRCh38"
        Path(path).unlink()

    def test_detect_grch37(self):
        """VAN-06: GRCh37 in header."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf(["##reference=GRCh37"])
        va = VCFAnnotator()
        assert va._detect_genome(Path(path)) == "GRCh37"
        Path(path).unlink()

    def test_detect_fallback(self):
        """VAN-07: No genome info → GRCh38 fallback."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf([])
        va = VCFAnnotator()
        assert va._detect_genome(Path(path)) == "GRCh38"
        Path(path).unlink()


@pytest.mark.l2
class TestVCFOpener:
    """VAN-08~09: VCF file opener (plain vs gzipped)."""

    def test_plain_opener(self):
        """VAN-08: Plain VCF → open."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            path = f.name
        va = VCFAnnotator()
        opener = va._vcf_opener(Path(path))
        assert opener is open
        Path(path).unlink()

    def test_gzipped_opener(self):
        """VAN-09: Gzipped VCF → gzip.open."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".vcf.gz", delete=False) as f:
            f.write(gzip.compress(b"##fileformat=VCFv4.2\n"))
            path = f.name
        va = VCFAnnotator()
        opener = va._vcf_opener(Path(path))
        assert opener is gzip.open
        Path(path).unlink()


@pytest.mark.l2
class TestVCFParsing:
    """VAN-10~17: VCF parsing and static helpers."""

    def _make_vcf_content(self, lines):
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE\n")
            for line in lines:
                f.write(line + "\n")
            return f.name

    def test_parse_basic(self):
        """VAN-10: Basic VCF parsing."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf_content(["1\t100\t.\tA\tG\t30\tPASS\tDP=20\tGT\t0/1"])
        va = VCFAnnotator()
        variants = va._parse_vcf(Path(path))
        assert len(variants) == 1
        assert variants[0]["chrom"] == "1"
        assert variants[0]["pos"] == 100
        assert variants[0]["ref"] == "A"
        assert variants[0]["alt"] == "G"
        assert variants[0]["qual"] == 30.0
        assert variants[0]["filter"] == "PASS"
        assert variants[0]["dp"] == 20
        assert variants[0]["gt"] == "0/1"
        Path(path).unlink()

    def test_parse_multi_alt(self):
        """VAN-11: Multiple ALT alleles → multiple variants."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf_content(["1\t100\t.\tA\tG,T\t30\tPASS\tDP=20\tGT\t0/1"])
        va = VCFAnnotator()
        variants = va._parse_vcf(Path(path))
        assert len(variants) == 2
        assert variants[0]["alt"] == "G"
        assert variants[1]["alt"] == "T"
        Path(path).unlink()

    def test_parse_dot_alt_skipped(self):
        """VAN-12: ALT='.' → skipped."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf_content(["1\t100\t.\tA\t.\t30\tPASS\tDP=20\tGT\t0/0"])
        va = VCFAnnotator()
        variants = va._parse_vcf(Path(path))
        assert len(variants) == 0
        Path(path).unlink()

    def test_parse_missing_dp(self):
        """VAN-13: No DP in INFO → dp=None."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf_content(["1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"])
        va = VCFAnnotator()
        variants = va._parse_vcf(Path(path))
        assert variants[0]["dp"] is None
        Path(path).unlink()

    def test_parse_no_sample(self):
        """VAN-14: No FORMAT/SAMPLE → gt='./.'."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("1\t100\t.\tA\tG\t30\tPASS\tDP=20\n")
            path = f.name
        va = VCFAnnotator()
        variants = va._parse_vcf(Path(path))
        assert variants[0]["gt"] == "./."
        Path(path).unlink()

    def test_parse_short_line(self):
        """VAN-14b: Line with fewer than 8 columns is skipped."""
        from gpa_vcf_annotator import VCFAnnotator
        path = self._make_vcf_content(["1\t100\t.\tA\tG\t30\tPASS"])
        va = VCFAnnotator()
        variants = va._parse_vcf(Path(path))
        assert len(variants) == 0
        Path(path).unlink()

    def test_normalize_chrom(self):
        """VAN-15: chr prefix stripped, M→MT."""
        from gpa_vcf_annotator import VCFAnnotator
        assert VCFAnnotator._normalize_chrom("chr1") == "1"
        assert VCFAnnotator._normalize_chrom("chrX") == "X"
        assert VCFAnnotator._normalize_chrom("M") == "MT"
        assert VCFAnnotator._normalize_chrom("chrM") == "MT"

    def test_extract_dp(self):
        """VAN-16: DP extraction from INFO."""
        from gpa_vcf_annotator import VCFAnnotator
        assert VCFAnnotator._extract_dp("DP=30;AF=0.5") == 30
        assert VCFAnnotator._extract_dp("AF=0.5") is None
        assert VCFAnnotator._extract_dp("DP=.") is None

    def test_extract_gt(self):
        """VAN-17: GT extraction from FORMAT/SAMPLE."""
        from gpa_vcf_annotator import VCFAnnotator
        assert VCFAnnotator._extract_gt("GT:DP", "0/1:30") == "0/1"
        assert VCFAnnotator._extract_gt("DP", "30") == "./."
        assert VCFAnnotator._extract_gt("GT:DP", "0/1") == "0/1"


@pytest.mark.l2
class TestPreFilter:
    """VAN-18~21: Pre-filtering logic."""

    def test_qual_filter(self):
        """VAN-18: QUAL < 20 → excluded."""
        from gpa_vcf_annotator import VCFAnnotator
        variants = [
            {"qual": 15, "dp": 20},
            {"qual": 25, "dp": 20},
        ]
        result = VCFAnnotator._prefilter(variants)
        assert len(result) == 1
        assert result[0]["qual"] == 25

    def test_dp_filter(self):
        """VAN-19: DP < 10 → excluded."""
        from gpa_vcf_annotator import VCFAnnotator
        variants = [
            {"qual": 30, "dp": 5},
            {"qual": 30, "dp": 15},
        ]
        result = VCFAnnotator._prefilter(variants)
        assert len(result) == 1
        assert result[0]["dp"] == 15

    def test_none_dp_passes(self):
        """VAN-20: DP=None → not filtered (only checked if present)."""
        from gpa_vcf_annotator import VCFAnnotator
        variants = [{"qual": 30, "dp": None}]
        result = VCFAnnotator._prefilter(variants)
        assert len(result) == 1

    def test_both_filters(self):
        """VAN-21: Both QUAL and DP filters applied."""
        from gpa_vcf_annotator import VCFAnnotator
        variants = [
            {"qual": 15, "dp": 5},
            {"qual": 15, "dp": 20},
            {"qual": 30, "dp": 5},
            {"qual": 30, "dp": 20},
        ]
        result = VCFAnnotator._prefilter(variants)
        assert len(result) == 1
        assert result[0]["qual"] == 30 and result[0]["dp"] == 20


@pytest.mark.l2
class TestAnnotatorResolution:
    """VAN-22~25: Annotator resolution."""

    def test_explicit_vep_api(self):
        """VAN-22: Explicit vep_api → no auto-resolve."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api")
        assert va._resolve_annotator(100) == "vep_api"

    def test_large_dataset_forces_api(self):
        """VAN-23: >5000 variants → vep_api."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="auto")
        assert va._resolve_annotator(6000) == "vep_api"

    def test_small_dataset_local_available(self):
        """VAN-24: <5000 variants + local VEP available → vep_local."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="auto")
        with patch.object(va, "_vep_local_available", return_value=True):
            assert va._resolve_annotator(100) == "vep_local"

    def test_small_dataset_local_unavailable(self):
        """VAN-25: <5000 variants + no local VEP → vep_api."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="auto")
        with patch.object(va, "_vep_local_available", return_value=False):
            assert va._resolve_annotator(100) == "vep_api"


@pytest.mark.l2
class TestIsAnnotatedVCF:
    """VAN-26~27: Check if VCF already contains annotation."""

    def test_annotated_vcf(self):
        """VAN-26: CSQ in header → True."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            path = f.name
        va = VCFAnnotator()
        assert va.is_annotated_vcf(path) is True
        Path(path).unlink()

    def test_unannotated_vcf(self):
        """VAN-27: No CSQ → False."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            path = f.name
        va = VCFAnnotator()
        assert va.is_annotated_vcf(path) is False
        Path(path).unlink()

    def test_is_annotated_exception_returns_false(self):
        """VAN-27b: Exception during file read → False."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        def broken_opener(path, mode):
            raise OSError("boom")
        with patch.object(va, "_vcf_opener", return_value=broken_opener):
            assert va.is_annotated_vcf("/nonexistent.vcf") is False


@pytest.mark.l2
class TestCheckpoint:
    """VAN-28~30: Checkpoint loading and saving."""

    @pytest.mark.asyncio
    async def test_load_from_checkpoint(self):
        """VAN-28: Existing checkpoint → load directly, skip VEP."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            json.dump([{"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "gene": "BRCA1", "extra_data": "x" * 200}], f)
            checkpoint = f.name
        va = VCFAnnotator(checkpoint_path=checkpoint)
        result = await va.annotate("/nonexistent.vcf")
        assert len(result) == 1
        assert result[0]["chrom"] == "1"
        Path(checkpoint).unlink()

    @pytest.mark.asyncio
    async def test_corrupt_checkpoint_rerun(self):
        """VAN-29: Corrupt checkpoint → rerun VEP (raises FileNotFound for VCF)."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False) as f:
            f.write("not json")
            checkpoint = f.name
        va = VCFAnnotator(checkpoint_path=checkpoint)
        with pytest.raises(FileNotFoundError):
            await va.annotate("/nonexistent.vcf")
        Path(checkpoint).unlink()

    @pytest.mark.asyncio
    async def test_save_checkpoint_on_success(self):
        """VAN-30: Successful annotation → checkpoint saved."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = Path(tmpdir) / "check.json"
            vcf_path = Path(tmpdir) / "test.vcf"
            vcf_path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t100\t.\tA\tG\t30\tPASS\tDP=20\n"
            )
            va = VCFAnnotator(checkpoint_path=str(checkpoint))
            with patch.object(va, "_annotate_internal", return_value=[{"chrom": "1", "pos": 100}]):
                result = await va.annotate(str(vcf_path))
            assert checkpoint.exists()
            data = json.loads(checkpoint.read_text())
            assert len(data) == 1

    @pytest.mark.asyncio
    async def test_checkpoint_save_failure_nonfatal(self):
        """VAN-30b: Checkpoint save failure is non-fatal."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = Path(tmpdir) / "check.json"
            vcf_path = Path(tmpdir) / "test.vcf"
            vcf_path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t100\t.\tA\tG\t30\tPASS\tDP=20\n"
            )
            va = VCFAnnotator(checkpoint_path=str(checkpoint))
            with patch.object(va, "_annotate_internal", return_value=[{"chrom": "1", "pos": 100}]):
                with patch("gpa_vcf_annotator.open", side_effect=OSError("disk full")):
                    result = await va.annotate(str(vcf_path))
            assert result == [{"chrom": "1", "pos": 100}]


@pytest.mark.l2
class TestShardManagement:
    """VAN-31~35: Shard-based incremental annotation."""

    def test_shard_path(self):
        """VAN-31: Shard path formatting."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(shard_dir="/tmp/shards")
        assert va._shard_path(0) == Path("/tmp/shards/shard_00000.json")
        assert va._shard_path(123) == Path("/tmp/shards/shard_00123.json")

    def test_index_path(self):
        """VAN-32: Index path."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(shard_dir="/tmp/shards")
        assert va._index_path() == Path("/tmp/shards/missing_by_shard.json")

    def test_no_shard_dir_returns_empty(self):
        """VAN-33: No shard_dir → empty completed set."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        assert va._get_completed_shards() == set()

    def test_save_shard_atomic(self):
        """VAN-34: Atomic shard write (tmp → rename)."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            va = VCFAnnotator(shard_dir=tmpdir)
            va._save_shard_atomic(0, [{"data": 1}])
            assert va._shard_path(0).exists()
            data = json.loads(va._shard_path(0).read_text())
            assert data == [{"data": 1}]

    def test_mark_shard_complete(self):
        """VAN-35: Mark shard complete updates index."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            va = VCFAnnotator(shard_dir=tmpdir)
            va._mark_shard_complete(0)
            index = json.loads(va._index_path().read_text())
            assert 0 in index.get("completed_shards", [])

    def test_init_shard_dir(self):
        """VAN-35b: Initialize shard directory with index."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            va = VCFAnnotator(shard_dir=tmpdir)
            va._init_shard_dir(2500)
            index_path = va._index_path()
            assert index_path.exists()
            data = json.loads(index_path.read_text())
            assert data["total_shards"] == 3
            assert data["shard_size"] == 1000
            assert data["completed_shards"] == []

    def test_init_shard_dir_no_resume(self):
        """VAN-35c: Re-initialize when resume=False overwrites index."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            va = VCFAnnotator(shard_dir=tmpdir, resume=False)
            va._init_shard_dir(100)
            first = json.loads(va._index_path().read_text())["created_at"]
            import time
            time.sleep(0.01)
            va._init_shard_dir(100)
            second = json.loads(va._index_path().read_text())["created_at"]
            assert second != first

    def test_get_completed_shards_with_index(self):
        """VAN-35d: Read completed shards from existing index."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            va = VCFAnnotator(shard_dir=tmpdir)
            va._index_path().write_text(json.dumps({"total_shards": 5, "completed_shards": [0, 2]}))
            assert va._get_completed_shards() == {0, 2}

    def test_get_completed_shards_corrupt_index(self):
        """VAN-35e: Corrupt index returns empty set."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            va = VCFAnnotator(shard_dir=tmpdir)
            va._index_path().write_text("not json")
            assert va._get_completed_shards() == set()


@pytest.mark.l2
class TestVEPBatchFailureError:
    """VAN-36~37: VEPBatchFailureError exception."""

    def test_exception_message(self):
        """VAN-36: Default message includes batch count."""
        from gpa_vcf_annotator import VEPBatchFailureError
        err = VEPBatchFailureError([{"batch": 1}, {"batch": 2}])
        assert "2 VEP batch(es) failed" in str(err)
        assert len(err.failed_batches) == 2

    def test_custom_message(self):
        """VAN-37: Custom message override."""
        from gpa_vcf_annotator import VEPBatchFailureError
        err = VEPBatchFailureError([], message="Custom error")
        assert str(err) == "Custom error"


# =============================================================================
# NEW: _annotate_internal flow
# =============================================================================

@pytest.mark.l2
class TestAnnotateInternal:
    """VAN-38~42: _annotate_internal full flow."""

    @pytest.mark.asyncio
    async def test_annotate_internal_file_not_found(self):
        """VAN-38: File not found raises FileNotFoundError."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        with pytest.raises(FileNotFoundError):
            await va._annotate_internal("/nonexistent/path.vcf")

    @pytest.mark.asyncio
    async def test_annotate_internal_empty_after_filter(self):
        """VAN-39: All variants filtered out → empty list."""
        from gpa_vcf_annotator import VCFAnnotator
        path = _make_vcf_content(["1\t100\t.\tA\tG\t10\tPASS\tDP=5"])
        va = VCFAnnotator()
        result = await va._annotate_internal(path)
        assert result == []
        Path(path).unlink()

    @pytest.mark.asyncio
    async def test_annotate_internal_vep_api_flow(self):
        """VAN-40: Full flow through VEP API annotator."""
        from gpa_vcf_annotator import VCFAnnotator
        path = _make_vcf_content(["1\t100\t.\tA\tG\t30\tPASS\tDP=20"])
        va = VCFAnnotator(annotator="vep_api")
        mock_annotated = [{"transcript_consequences": [{"gene": "TP53"}], "vep_summary": {}}]
        with patch.object(va, "_annotate_vep_api", return_value=mock_annotated):
            result = await va._annotate_internal(path)
        assert len(result) == 1
        assert result[0]["transcript_consequences"] == [{"gene": "TP53"}]
        assert result[0]["chrom"] == "1"
        Path(path).unlink()

    @pytest.mark.asyncio
    async def test_annotate_internal_vep_local_flow(self):
        """VAN-41: Full flow through local VEP annotator."""
        from gpa_vcf_annotator import VCFAnnotator
        path = _make_vcf_content(["1\t100\t.\tA\tG\t30\tPASS\tDP=20"])
        va = VCFAnnotator(annotator="vep_local")
        mock_annotated = [{"transcript_consequences": [{"gene": "BRCA1"}], "vep_summary": {}}]
        with patch.object(va, "_annotate_vep_local", return_value=mock_annotated):
            result = await va._annotate_internal(path)
        assert len(result) == 1
        assert result[0]["transcript_consequences"] == [{"gene": "BRCA1"}]
        Path(path).unlink()

    @pytest.mark.asyncio
    async def test_annotate_internal_unsupported_annotator(self):
        """VAN-42: Unsupported annotator raises NotImplementedError."""
        from gpa_vcf_annotator import VCFAnnotator
        path = _make_vcf_content(["1\t100\t.\tA\tG\t30\tPASS\tDP=20"])
        va = VCFAnnotator(annotator="annovar")
        with pytest.raises(NotImplementedError):
            await va._annotate_internal(path)
        Path(path).unlink()


# =============================================================================
# NEW: VEP API annotation
# =============================================================================

def _mock_aiohttp_session(status=200, json_data=None, headers=None):
    """Build a mock aiohttp ClientSession that returns given status/json."""
    mock_resp = AsyncMock()
    mock_resp.status = status
    mock_resp.headers = headers or {}
    mock_resp.json = AsyncMock(return_value=json_data or [])

    mock_session = AsyncMock()
    mock_session.post = MagicMock(return_value=AsyncMock(__aenter__=AsyncMock(return_value=mock_resp), __aexit__=AsyncMock(return_value=False)))
    mock_session.close = AsyncMock()
    return mock_session


@pytest.mark.l2
class TestAnnotateVEPAPI:
    """VAN-43~56: _annotate_vep_api async methods."""

    @pytest.mark.asyncio
    async def test_annotate_vep_api_success(self):
        """VAN-43: Successful VEP API annotation."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        vep_response = [
            {
                "most_severe_consequence": "missense_variant",
                "variant_class": "SNV",
                "transcript_consequences": [
                    {
                        "transcript_id": "ENST000001",
                        "gene_symbol": "TP53",
                        "gene_id": "ENSG000001",
                        "consequence_terms": ["missense_variant"],
                        "impact": "MODERATE",
                        "hgvsc": "c.818C>T",
                        "hgvsp": "p.Arg273Ter",
                        "canonical": 1,
                        "mane_select": 1,
                        "exon": "5/11",
                        "intron": "",
                        "domains": ["Pfam:PF123"],
                    }
                ],
            }
        ]
        mock_session = _mock_aiohttp_session(status=200, json_data=vep_response)
        with patch("aiohttp.ClientSession", return_value=mock_session):
            result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        assert result[0]["transcript_consequences"][0]["gene_symbol"] == "TP53"
        assert result[0]["vep_summary"]["most_severe_consequence"] == "missense_variant"

    @pytest.mark.asyncio
    async def test_annotate_vep_api_grch37_refseq(self):
        """VAN-44: GRCh37 adds refseq=1 parameter."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        vep_response = [{"transcript_consequences": [], "most_severe_consequence": ""}]
        mock_session = _mock_aiohttp_session(status=200, json_data=vep_response)
        with patch("aiohttp.ClientSession", return_value=mock_session):
            await va._annotate_vep_api(variants, "GRCh37")
        call_args = mock_session.post.call_args
        assert call_args[1]["params"]["refseq"] == "1"

    @pytest.mark.asyncio
    async def test_annotate_vep_api_429_retry_then_success(self):
        """VAN-45: 429 rate limit retries then succeeds."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_429 = AsyncMock()
        resp_429.status = 429
        resp_429.headers = {"Retry-After": "1"}
        resp_ok = AsyncMock()
        resp_ok.status = 200
        resp_ok.headers = {}
        resp_ok.json = AsyncMock(return_value=[{"transcript_consequences": []}])

        mock_session = AsyncMock()
        mock_session.post = MagicMock(side_effect=[
            AsyncMock(__aenter__=AsyncMock(return_value=resp_429), __aexit__=AsyncMock(return_value=False)),
            AsyncMock(__aenter__=AsyncMock(return_value=resp_ok), __aexit__=AsyncMock(return_value=False)),
        ])

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_annotate_vep_api_timeout_then_success(self):
        """VAN-46: Timeout retries then succeeds."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_ok = AsyncMock()
        resp_ok.status = 200
        resp_ok.headers = {}
        resp_ok.json = AsyncMock(return_value=[{"transcript_consequences": []}])

        mock_session = AsyncMock()
        post_mock = MagicMock(side_effect=[
            asyncio.TimeoutError(),
            AsyncMock(__aenter__=AsyncMock(return_value=resp_ok), __aexit__=AsyncMock(return_value=False)),
        ])
        mock_session.post = post_mock

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_annotate_vep_api_client_error_then_success(self):
        """VAN-47: ClientError retries then succeeds."""
        import aiohttp
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_ok = AsyncMock()
        resp_ok.status = 200
        resp_ok.headers = {}
        resp_ok.json = AsyncMock(return_value=[{"transcript_consequences": []}])

        mock_session = AsyncMock()
        mock_session.post = MagicMock(side_effect=[
            aiohttp.ClientError("connection reset"),
            AsyncMock(__aenter__=AsyncMock(return_value=resp_ok), __aexit__=AsyncMock(return_value=False)),
        ])

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_annotate_vep_api_os_error_then_success(self):
        """VAN-48: OSError retries then succeeds."""
        import aiohttp
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_ok = AsyncMock()
        resp_ok.status = 200
        resp_ok.headers = {}
        resp_ok.json = AsyncMock(return_value=[{"transcript_consequences": []}])

        mock_session = AsyncMock()
        mock_session.post = MagicMock(side_effect=[
            OSError("network unreachable"),
            AsyncMock(__aenter__=AsyncMock(return_value=resp_ok), __aexit__=AsyncMock(return_value=False)),
        ])

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_annotate_vep_api_unexpected_error_then_success(self):
        """VAN-49: Unexpected exception retries then succeeds."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=10)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_ok = AsyncMock()
        resp_ok.status = 200
        resp_ok.headers = {}
        resp_ok.json = AsyncMock(return_value=[{"transcript_consequences": []}])

        mock_session = AsyncMock()
        mock_session.post = MagicMock(side_effect=[
            ValueError("unexpected"),
            AsyncMock(__aenter__=AsyncMock(return_value=resp_ok), __aexit__=AsyncMock(return_value=False)),
        ])

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        mock_sleep.assert_awaited()

    @pytest.mark.asyncio
    async def test_annotate_vep_api_progress_callback(self):
        """VAN-50: Progress callback is invoked during annotation."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=1)
        variants = [
            {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
            {"chrom": "1", "pos": 101, "ref": "C", "alt": "T"},
        ]
        vep_response = [
            {"transcript_consequences": [], "most_severe_consequence": ""},
            {"transcript_consequences": [], "most_severe_consequence": ""},
        ]
        mock_session = _mock_aiohttp_session(status=200, json_data=vep_response)
        progress_calls = []

        def progress(done, total):
            progress_calls.append((done, total))

        with patch("aiohttp.ClientSession", return_value=mock_session):
            result = await va._annotate_vep_api(variants, "GRCh38", progress_callback=progress)
        assert len(result) == 2
        assert len(progress_calls) >= 1
        assert progress_calls[-1][1] == 2

    @pytest.mark.asyncio
    async def test_annotate_vep_api_failed_batches_skip(self):
        """VAN-51: Failed batches with interactive=False → auto-skip."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=1, interactive=False)
        variants = [
            {"chrom": "1", "pos": 100, "ref": "A", "alt": "G"},
            {"chrom": "1", "pos": 101, "ref": "C", "alt": "T"},
        ]

        resp_500 = AsyncMock()
        resp_500.status = 500
        resp_500.headers = {}

        mock_session = AsyncMock()
        mock_session.post = MagicMock(return_value=AsyncMock(
            __aenter__=AsyncMock(return_value=resp_500),
            __aexit__=AsyncMock(return_value=False),
        ))

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock):
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 2
        assert result[0]["vep_summary"].get("error", "").startswith("VEP API failed")

    @pytest.mark.asyncio
    async def test_annotate_vep_api_failed_batches_abort(self):
        """VAN-52: Failed batches with interactive abort choice."""
        from gpa_vcf_annotator import VCFAnnotator, VEPBatchFailureError
        va = VCFAnnotator(annotator="vep_api", batch_size=1, interactive=True)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_500 = AsyncMock()
        resp_500.status = 500
        resp_500.headers = {}

        mock_session = AsyncMock()
        mock_session.post = MagicMock(return_value=AsyncMock(
            __aenter__=AsyncMock(return_value=resp_500),
            __aexit__=AsyncMock(return_value=False),
        ))

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock):
                with patch("asyncio.to_thread", return_value="a") as mock_input:
                    with pytest.raises(VEPBatchFailureError):
                        await va._annotate_vep_api(variants, "GRCh38")

    @pytest.mark.asyncio
    async def test_annotate_vep_api_proxy_route_map(self):
        """VAN-53: Per-API proxy routing via proxy_route_map."""
        from gpa_vcf_annotator import VCFAnnotator
        route_map = MagicMock()
        route_map.get_proxy = MagicMock(return_value="http://proxy:8080")
        va = VCFAnnotator(annotator="vep_api", proxy_route_map=route_map)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        mock_session = _mock_aiohttp_session(status=200, json_data=[{"transcript_consequences": []}])
        with patch("aiohttp.ClientSession", return_value=mock_session):
            await va._annotate_vep_api(variants, "GRCh38")
        route_map.get_proxy.assert_called_once_with("ensembl")

    @pytest.mark.asyncio
    async def test_annotate_vep_api_direct_proxy(self):
        """VAN-54: __DIRECT__ proxy disables system proxy."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", proxy="__DIRECT__")
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        mock_session = _mock_aiohttp_session(status=200, json_data=[{"transcript_consequences": []}])
        with patch("aiohttp.ClientSession", return_value=mock_session) as mock_cls:
            await va._annotate_vep_api(variants, "GRCh38")
        assert mock_cls.call_args[1]["trust_env"] is False

    @pytest.mark.asyncio
    async def test_annotate_vep_api_reuse_session(self):
        """VAN-55: Existing session is reused."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api")
        mock_session = _mock_aiohttp_session(status=200, json_data=[{"transcript_consequences": []}])
        va._session = mock_session
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        with patch("aiohttp.ClientSession") as mock_cls:
            await va._annotate_vep_api(variants, "GRCh38")
        mock_cls.assert_not_called()

    @pytest.mark.asyncio
    async def test_annotate_vep_api_500_all_retries_fail(self):
        """VAN-56: HTTP 500 after all retries → failed batch."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api", batch_size=1, interactive=False)
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]

        resp_500 = AsyncMock()
        resp_500.status = 500
        resp_500.headers = {}

        mock_session = AsyncMock()
        mock_session.post = MagicMock(return_value=AsyncMock(
            __aenter__=AsyncMock(return_value=resp_500),
            __aexit__=AsyncMock(return_value=False),
        ))

        with patch("aiohttp.ClientSession", return_value=mock_session):
            with patch("asyncio.sleep", new_callable=AsyncMock) as mock_sleep:
                result = await va._annotate_vep_api(variants, "GRCh38")
        assert len(result) == 1
        assert "VEP API failed" in result[0]["vep_summary"]["error"]
        assert mock_sleep.await_count >= 3


@pytest.mark.l2
class TestQuerySingleBatch:
    """VAN-57~58: _query_single_batch retry logic."""

    @pytest.mark.asyncio
    async def test_query_single_batch_success(self):
        """VAN-57: Retry batch succeeds."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api")
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        annotated = [None]
        semaphore = asyncio.Semaphore(1)
        mock_session = _mock_aiohttp_session(status=200, json_data=[{"transcript_consequences": []}])
        va._session = mock_session
        success = await va._query_single_batch(0, variants, "GRCh38", semaphore, None, annotated)
        assert success is True
        assert annotated[0] is not None

    @pytest.mark.asyncio
    async def test_query_single_batch_failure(self):
        """VAN-58: Retry batch fails with non-200 status."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(annotator="vep_api")
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G"}]
        annotated = [None]
        semaphore = asyncio.Semaphore(1)
        mock_session = _mock_aiohttp_session(status=500, json_data=[])
        va._session = mock_session
        success = await va._query_single_batch(0, variants, "GRCh38", semaphore, None, annotated)
        assert success is False
        assert annotated[0] is None


@pytest.mark.l2
class TestParseVEPResponse:
    """VAN-59~62: _parse_vep_response static method."""

    def test_parse_normal(self):
        """VAN-59: Normal VEP response parsing."""
        from gpa_vcf_annotator import VCFAnnotator
        data = [
            {
                "most_severe_consequence": "missense_variant",
                "variant_class": "SNV",
                "colocated_variants": [{"id": "rs123"}],
                "transcript_consequences": [
                    {
                        "transcript_id": "ENST001",
                        "gene_symbol": "TP53",
                        "gene_id": "ENSG001",
                        "consequence_terms": ["missense_variant"],
                        "impact": "MODERATE",
                        "hgvsc": "c.1A>G",
                        "hgvsp": "p.M1V",
                        "canonical": 1,
                        "mane_select": 1,
                        "mane_plus_clinical": 0,
                        "exon": "2/10",
                        "intron": "",
                        "domains": [],
                    }
                ],
            }
        ]
        batch = [{"chrom": "1", "pos": 100}]
        result = VCFAnnotator._parse_vep_response(data, batch)
        assert len(result) == 1
        assert result[0]["vep_summary"]["most_severe_consequence"] == "missense_variant"
        assert "colocated_variants" in result[0]["vep_summary"]
        tx = result[0]["transcript_consequences"][0]
        assert tx["transcript_id"] == "ENST001"
        assert tx["protein_domains"] == []

    def test_parse_invalid_entry(self):
        """VAN-60: Non-dict entry → error annotation."""
        from gpa_vcf_annotator import VCFAnnotator
        data = ["not_a_dict"]
        batch = [{"chrom": "1", "pos": 100}]
        result = VCFAnnotator._parse_vep_response(data, batch)
        assert result[0]["vep_summary"]["error"] == "Invalid VEP response format"

    def test_parse_shorter_response(self):
        """VAN-61: VEP returns fewer results than batch → padded with errors."""
        from gpa_vcf_annotator import VCFAnnotator
        data = [{"transcript_consequences": []}]
        batch = [{"chrom": "1", "pos": 100}, {"chrom": "1", "pos": 101}]
        result = VCFAnnotator._parse_vep_response(data, batch)
        assert len(result) == 2
        assert result[1]["vep_summary"]["error"] == "VEP response shorter than input batch"

    def test_parse_empty_transcript_consequences(self):
        """VAN-62: No transcript consequences → empty list."""
        from gpa_vcf_annotator import VCFAnnotator
        data = [{"most_severe_consequence": "intergenic_variant", "variant_class": "SNV"}]
        batch = [{"chrom": "1", "pos": 100}]
        result = VCFAnnotator._parse_vep_response(data, batch)
        assert result[0]["transcript_consequences"] == []


@pytest.mark.l2
class TestAnnotateVEPLocal:
    """VAN-63~66: _annotate_vep_local async wrapper."""

    @pytest.mark.asyncio
    async def test_annotate_vep_local_success(self):
        """VAN-63: Local VEP subprocess succeeds."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(vep_cache="/tmp/vep_cache")
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "qual": 30}]

        mock_proc = AsyncMock()
        mock_proc.returncode = 0
        mock_proc.communicate = AsyncMock(return_value=(b"", b""))

        with patch("asyncio.create_subprocess_exec", return_value=mock_proc):
            with patch.object(va, "_parse_vep_local_output", return_value=[{"transcript_consequences": []}]):
                result = await va._annotate_vep_local(variants, "GRCh38")
        assert len(result) == 1

    @pytest.mark.asyncio
    async def test_annotate_vep_local_failure_exit_code(self):
        """VAN-64: Local VEP non-zero exit code → error annotations."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "qual": 30}]

        mock_proc = AsyncMock()
        mock_proc.returncode = 1
        mock_proc.communicate = AsyncMock(return_value=(b"", b"VEP ERROR"))

        with patch("asyncio.create_subprocess_exec", return_value=mock_proc):
            result = await va._annotate_vep_local(variants, "GRCh38")
        assert len(result) == 1
        assert "Local VEP exit 1" in result[0]["vep_summary"]["error"]

    @pytest.mark.asyncio
    async def test_annotate_vep_local_execution_error(self):
        """VAN-65: Local VEP subprocess raises exception."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "qual": 30}]

        with patch("asyncio.create_subprocess_exec", side_effect=OSError("vep not found")):
            result = await va._annotate_vep_local(variants, "GRCh38")
        assert len(result) == 1
        assert "vep not found" in result[0]["vep_summary"]["error"]

    @pytest.mark.asyncio
    async def test_annotate_vep_local_grch37(self):
        """VAN-66: GRCh37 uses correct assembly parameter."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator(vep_cache="/tmp/vep_cache")
        variants = [{"chrom": "1", "pos": 100, "ref": "A", "alt": "G", "qual": 30}]

        mock_proc = AsyncMock()
        mock_proc.returncode = 0
        mock_proc.communicate = AsyncMock(return_value=(b"", b""))

        with patch("asyncio.create_subprocess_exec", return_value=mock_proc) as mock_exec:
            with patch.object(va, "_parse_vep_local_output", return_value=[{"transcript_consequences": []}]):
                await va._annotate_vep_local(variants, "GRCh37")
        cmd_args = mock_exec.call_args[0]
        assert "GRCh37" in cmd_args


@pytest.mark.l2
class TestParseVEPLocalOutput:
    """VAN-67~69: _parse_vep_local_output CSQ parsing."""

    def test_parse_local_with_csq(self):
        """VAN-67: Output VCF with CSQ field → parsed transcripts.

        The code maps CSQ fields by hardcoded positions:
        0=allele, 1=consequence, 2=impact, 3=gene_symbol, 4=gene_id,
        5=unused, 6=transcript_id, 7=unused, 8=unused, 9=hgvsc, 10=hgvsp
        """
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence">\n')
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            # Position 6 must be transcript_id, position 3 must be gene_symbol
            f.write("1\t100\t.\tA\tG\t30\tPASS\tCSQ=G|missense_variant|MODERATE|TP53|ENSG000001||ENST000001|||c.1A>G|p.M1V\n")
            path = f.name
        result = va._parse_vep_local_output(path, [{"chrom": "1", "pos": 100}])
        assert len(result) == 1
        assert result[0]["transcript_consequences"][0]["gene_symbol"] == "TP53"
        assert result[0]["transcript_consequences"][0]["transcript_id"] == "ENST000001"
        assert result[0]["transcript_consequences"][0]["hgvsc"] == "c.1A>G"
        assert result[0]["transcript_consequences"][0]["hgvsp"] == "p.M1V"
        Path(path).unlink()

    def test_parse_local_without_csq(self):
        """VAN-68: Output VCF without CSQ → empty transcripts."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        with tempfile.NamedTemporaryFile(mode="w", suffix=".vcf", delete=False) as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            f.write("1\t100\t.\tA\tG\t30\tPASS\tDP=20\n")
            path = f.name
        result = va._parse_vep_local_output(path, [{"chrom": "1", "pos": 100}])
        assert len(result) == 1
        assert result[0]["transcript_consequences"] == []
        Path(path).unlink()

    def test_parse_local_gzipped(self):
        """VAN-69: Gzipped output VCF is parsed correctly."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        with tempfile.NamedTemporaryFile(mode="wb", suffix=".vcf.gz", delete=False) as f:
            content = (
                "##fileformat=VCFv4.2\n"
                + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                + "1\t100\t.\tA\tG\t30\tPASS\tCSQ=G|synonymous_variant|LOW|BRCA1|ENSG000002|ENST000002||||c.1A>G\n"
            )
            f.write(gzip.compress(content.encode()))
            path = f.name
        result = va._parse_vep_local_output(path, [{"chrom": "1", "pos": 100}])
        assert len(result) == 1
        assert result[0]["transcript_consequences"][0]["gene_symbol"] == "BRCA1"
        Path(path).unlink()


@pytest.mark.l2
class TestClose:
    """VAN-70: Session cleanup."""

    @pytest.mark.asyncio
    async def test_close_session(self):
        """VAN-70: close() clears the session."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        mock_session = AsyncMock()
        va._session = mock_session
        await va.close()
        mock_session.close.assert_awaited_once()
        assert va._session is None

    @pytest.mark.asyncio
    async def test_close_no_session(self):
        """VAN-70b: close() with no session is a no-op."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        await va.close()
        assert va._session is None


@pytest.mark.l2
class TestAnnotateEntryPoint:
    """VAN-71~73: annotate() entry point behavior."""

    @pytest.mark.asyncio
    async def test_annotate_calls_internal(self):
        """VAN-71: annotate() delegates to _annotate_internal."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        with patch.object(va, "_annotate_internal", return_value=[{"chrom": "1"}]) as mock_internal:
            with patch.object(va, "close", new_callable=AsyncMock):
                result = await va.annotate("/tmp/test.vcf")
        mock_internal.assert_awaited_once()
        assert result == [{"chrom": "1"}]

    @pytest.mark.asyncio
    async def test_annotate_cleans_up_on_error(self):
        """VAN-72: Exception in annotate() triggers close()."""
        from gpa_vcf_annotator import VCFAnnotator
        va = VCFAnnotator()
        with patch.object(va, "_annotate_internal", side_effect=RuntimeError("boom")):
            with patch.object(va, "close", new_callable=AsyncMock) as mock_close:
                with pytest.raises(RuntimeError):
                    await va.annotate("/tmp/test.vcf")
        mock_close.assert_awaited_once()

    @pytest.mark.asyncio
    async def test_annotate_checkpoint_empty_result(self):
        """VAN-73: Empty result does not write checkpoint."""
        from gpa_vcf_annotator import VCFAnnotator
        with tempfile.TemporaryDirectory() as tmpdir:
            checkpoint = Path(tmpdir) / "check.json"
            vcf_path = Path(tmpdir) / "test.vcf"
            vcf_path.write_text(
                "##fileformat=VCFv4.2\n"
                "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
                "1\t100\t.\tA\tG\t30\tPASS\tDP=20\n"
            )
            va = VCFAnnotator(checkpoint_path=str(checkpoint))
            with patch.object(va, "_annotate_internal", return_value=[]):
                result = await va.annotate(str(vcf_path))
            assert result == []
            assert not checkpoint.exists()
