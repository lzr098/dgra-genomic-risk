"""
GPA Test Fixtures — pytest configuration + mock utilities + recording-replay infrastructure.

This file serves dual purposes:
1. pytest fixtures for the new pytest-based test suite
2. Backward-compatible mock utilities for the legacy assert-based tests

Usage:
    pytest tests/                          # Run all pytest tests
    python tests/l2_unit/test_xxx.py       # Run legacy standalone tests
"""

from __future__ import annotations

import asyncio
import hashlib
import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Union
from unittest.mock import AsyncMock, MagicMock, patch

import pytest

# =============================================================================
# Path setup (always run, needed by both pytest and standalone)
# =============================================================================

SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

TESTS_DIR = Path(__file__).parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))


# =============================================================================
# Recording-Replay Infrastructure
# =============================================================================

class MissingRecordingError(Exception):
    """Raised when a recording is required but not found in strict mode."""


class APIRecorder:
    """
    Record and replay API responses for stable, fast tests.

    Modes:
        - "playback" (default): Load from recording files. Fail if missing in strict mode.
        - "record": Call real API, save response to file.
        - "refresh": Re-record all, overwriting existing files.
        - "strict": Like playback, but raise MissingRecordingError if file missing.

    Example:
        recorder = APIRecorder(recording_dir=Path("tests/recording"), mode="playback")
        with recorder.record("gnomad", "1-100000-A-G") as response:
            if response is None:
                # Make real API call
                real_response = await fetch_from_gnomad(...)
                recorder.save("gnomad", "1-100000-A-G", real_response)
            else:
                # Use recorded response
                assert response["status"] == "SUCCESS"
    """

    def __init__(self, recording_dir: Path, mode: str = "playback"):
        self.recording_dir = recording_dir.resolve()
        self.mode = mode
        self._index: Dict[str, str] = {}  # key -> filepath
        self._load_index()

    def _load_index(self):
        """Load the recording index file if it exists."""
        index_path = self.recording_dir / ".index.json"
        if index_path.exists():
            self._index = json.loads(index_path.read_text(encoding="utf-8"))

    def save_index(self):
        """Save the recording index to disk."""
        index_path = self.recording_dir / ".index.json"
        index_path.write_text(json.dumps(self._index, indent=2, ensure_ascii=False), encoding="utf-8")

    def _recording_path(self, api_name: str, key: str) -> Path:
        """Compute the recording file path for a given API and key."""
        safe_key = hashlib.md5(key.encode()).hexdigest()[:12]
        subdir = self.recording_dir / api_name
        subdir.mkdir(parents=True, exist_ok=True)
        return subdir / f"{safe_key}.json"

    def _key(self, api_name: str, variant_signature: str, endpoint: str = "") -> str:
        """Generate a unique cache key."""
        return f"{api_name}:{variant_signature}:{endpoint}"

    def get(self, api_name: str, variant_signature: str, endpoint: str = "") -> Optional[Dict[str, Any]]:
        """
        Retrieve a recorded response.
        Returns None if not found (caller should make real API call in record mode).
        Raises MissingRecordingError in strict mode if not found.
        """
        key = self._key(api_name, variant_signature, endpoint)
        path = self._recording_path(api_name, key)

        if path.exists():
            return json.loads(path.read_text(encoding="utf-8"))

        if self.mode == "strict":
            raise MissingRecordingError(
                f"Missing recording for {api_name}/{variant_signature}. "
                f"Run with --record-mode=record to create it."
            )

        return None

    def save(self, api_name: str, variant_signature: str, response: Dict[str, Any], endpoint: str = ""):
        """Save a response to the recording directory."""
        key = self._key(api_name, variant_signature, endpoint)
        path = self._recording_path(api_name, key)

        # Wrap with metadata
        wrapped = {
            "meta": {
                "api_name": api_name,
                "variant_id": variant_signature,
                "endpoint": endpoint,
                "recorded_at": _iso_now(),
                "mode": self.mode,
            },
            "response": response,
        }

        path.write_text(json.dumps(wrapped, indent=2, ensure_ascii=False), encoding="utf-8")
        self._index[key] = str(path.relative_to(self.recording_dir))

    def record(self, api_name: str, variant_signature: str, endpoint: str = ""):
        """
        Context manager for record/replay pattern.

        Usage:
            with recorder.record("gnomad", "1-100000-A-G") as resp:
                if resp is None:
                    resp = await real_api_call()
                    recorder.save("gnomad", "1-100000-A-G", resp)
                process(resp)
        """
        return _RecordContext(self, api_name, variant_signature, endpoint)


class _RecordContext:
    """Context manager helper for APIRecorder.record()."""

    def __init__(self, recorder: APIRecorder, api_name: str, variant_signature: str, endpoint: str):
        self.recorder = recorder
        self.api_name = api_name
        self.variant_signature = variant_signature
        self.endpoint = endpoint
        self._response: Optional[Dict] = None
        self._from_recording = False

    def __enter__(self) -> Optional[Dict]:
        if self.recorder.mode in ("playback", "strict"):
            self._response = self.recorder.get(self.api_name, self.variant_signature, self.endpoint)
            self._from_recording = self._response is not None
            return self._response
        # record / refresh mode: always return None to trigger real call
        return None

    def __exit__(self, exc_type, exc_val, exc_tb):
        # If real call was made and succeeded, save it (in record/refresh mode)
        if self.recorder.mode in ("record", "refresh") and not exc_type:
            if self._response is not None and not self._from_recording:
                self.recorder.save(self.api_name, self.variant_signature, self._response, self.endpoint)
        return False


def _iso_now() -> str:
    """Return current ISO 8601 timestamp."""
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()


# =============================================================================
# pytest CLI Option: --record-mode
# =============================================================================

def pytest_addoption(parser):
    parser.addoption(
        "--record-mode",
        action="store",
        default="playback",
        choices=["playback", "record", "refresh", "strict"],
        help="API recording mode: playback (default), record, refresh, strict",
    )


# =============================================================================
# pytest Fixtures
# =============================================================================

@pytest.fixture(scope="session")
def record_mode(request) -> str:
    """Return the --record-mode CLI option value."""
    return request.config.getoption("--record-mode")


@pytest.fixture(scope="session")
def recording_dir() -> Path:
    """Return the path to the recording directory."""
    return TESTS_DIR / "recording"


@pytest.fixture
def api_recorder(record_mode: str, recording_dir: Path) -> APIRecorder:
    """Provide an APIRecorder instance for test functions."""
    recorder = APIRecorder(recording_dir=recording_dir, mode=record_mode)
    yield recorder
    recorder.save_index()


@pytest.fixture
def mock_gnomad_common() -> Dict[str, Any]:
    """Return a mocked gnomAD common variant response."""
    return MockGnomAD.success_common("X", 41357831, "A", "T", af=0.45)


@pytest.fixture
def mock_gnomad_rare() -> Dict[str, Any]:
    """Return a mocked gnomAD rare variant response."""
    return MockGnomAD.success_rare("17", 7579472, "C", "T", af=0.00001)


@pytest.fixture
def mock_gnomad_failed() -> Dict[str, Any]:
    """Return a mocked gnomAD API_FAILED response."""
    return MockGnomAD.api_failed("1", 100, "A", "G")


@pytest.fixture
def mock_tissue_profile_hematopoietic() -> Dict[str, Any]:
    """Return a mocked hematopoietic tissue profile."""
    return MockTissueProfile.hematopoietic()


@pytest.fixture
def mock_tissue_profile_general() -> Dict[str, Any]:
    """Return a mocked general tissue profile."""
    return MockTissueProfile.general()


@pytest.fixture
def mock_tissue_primary() -> Dict[str, Any]:
    """Return a mocked primary tissue assessment."""
    return MockTissueAssessment.primary()


@pytest.fixture
def mock_tissue_none() -> Dict[str, Any]:
    """Return a mocked 'none' tissue assessment (fast-track to Tier 3)."""
    return MockTissueAssessment.none()


@pytest.fixture
def mock_variant_factory() -> Callable:
    """Return the make_variant factory function."""
    return make_variant


@pytest.fixture
def mock_ensembl_vep() -> MockEnsemblVEP:
    """Return the MockEnsemblVEP class."""
    return MockEnsemblVEP()


# =============================================================================
# Async helpers
# =============================================================================

@pytest.fixture
def event_loop():
    """Provide a fresh event loop for each test."""
    loop = asyncio.new_event_loop()
    yield loop
    loop.close()


# =============================================================================
# Legacy Mock Utilities (backward compatible with standalone tests)
# =============================================================================

class MockGnomAD:
    """Mock gnomAD API responses."""

    @staticmethod
    def success_common(chrom: str, pos: int, ref: str, alt: str, af: float = 0.45):
        """Return a SUCCESS response for a common variant (DDX3X-like)."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": af,
            "af_popmax": af,
            "af_populations": {
                "EAS": {"af": af, "ac": 100, "an": 222},
                "NFE": {"af": 0.35, "ac": 80, "an": 228},
            },
            "af_exome": af,
            "af_genome": None,
            "an_exome": 222,
            "an_genome": None,
            "hom_count": 12,
            "status": "SUCCESS",
            "source": "gnomad",
            "confidence": "medium",
            "raw": {},
        }

    @staticmethod
    def success_rare(chrom: str, pos: int, ref: str, alt: str, af: float = 0.0001):
        """Return a SUCCESS response for a rare variant."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": af,
            "af_popmax": af,
            "af_populations": {
                "EAS": {"af": af, "ac": 1, "an": 10000},
            },
            "af_exome": af,
            "af_genome": None,
            "an_exome": 10000,
            "an_genome": None,
            "hom_count": 0,
            "status": "SUCCESS",
            "source": "gnomad",
            "confidence": "medium",
            "raw": {},
        }

    @staticmethod
    def api_failed(chrom: str, pos: int, ref: str, alt: str):
        """Return an API_FAILED response."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": None,
            "af_popmax": None,
            "af_populations": {},
            "status": "API_FAILED",
            "source": "failed",
            "confidence": "low",
            "error": "GraphQL 400: unknown error",
        }

    @staticmethod
    def not_captured(chrom: str, pos: int, ref: str, alt: str):
        """Return a NOT_CAPTURED response."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": None,
            "af_popmax": None,
            "af_populations": {},
            "status": "NOT_CAPTURED",
            "source": "gnomad",
            "confidence": "medium",
            "note": "Variant not in gnomAD dataset",
            "raw": {},
        }


class MockEnsemblVEP:
    """Mock Ensembl VEP responses."""

    @staticmethod
    def vep_result(variant_id: str, gene: str, consequence: str, impact: str,
                   hgvsp: str = "", hgvsc: str = "", exon: str = ""):
        """Return a minimal VEP result dict."""
        return {
            "input": variant_id,
            "transcript_consequences": [{
                "gene_symbol": gene,
                "consequence_terms": [consequence],
                "impact": impact,
                "hgvsp": hgvsp,
                "hgvsc": hgvsc,
                "exon": exon,
            }]
        }

    @staticmethod
    def batch_results(variants: List[Dict]) -> List[Dict]:
        """Return mocked VEP results for a batch of variants."""
        return [MockEnsemblVEP.vep_result(**v) for v in variants]


class MockTissueProfile:
    """Mock tissue profile for testing tier classification."""

    @staticmethod
    def hematopoietic():
        return {
            "display_name": "hematopoietic",
            "tier_rules": {},
            "special_gene_lists": {
                "coagulation": {"VWF", "F8", "F9"},
                "fa_dna_repair": {"BRCA1", "BRCA2", "FANCA"},
                "drug_metabolism": {"CYP2D6", "TPMT", "DPYD"},
            },
            "tissue_genes": {
                "RUNX1", "CEBPA", "GATA2", "ASXL1", "BCOR", "BCORL1", "PHF6",
                "FLT3", "NPM1", "IDH1", "IDH2", "DNMT3A", "TET2",
                "TP53", "KIT", "NRAS", "KRAS", "PTPN11",
            },
        }

    @staticmethod
    def general():
        return {
            "display_name": "general",
            "tier_rules": {},
            "special_gene_lists": {
                "coagulation": {"VWF", "F8", "F9"},
                "fa_dna_repair": {"BRCA1", "BRCA2"},
            },
            "tissue_genes": set(),
        }

    @staticmethod
    def cardiovascular():
        return {
            "display_name": "cardiovascular",
            "tier_rules": {},
            "special_gene_lists": {},
            "tissue_genes": {"VWF", "F8", "F9", "MYH7", "MYBPC3", "TTN"},
        }

    @staticmethod
    def neurological():
        return {
            "display_name": "neurological",
            "tier_rules": {},
            "special_gene_lists": {},
            "tissue_genes": {"DMD", "SMN1", "HTT", "C9orf72"},
        }

    @staticmethod
    def hepatic():
        return {
            "display_name": "hepatic",
            "tier_rules": {},
            "special_gene_lists": {},
            "tissue_genes": {"CFTR", "HFE", "ALDOB"},
        }

    @staticmethod
    def renal():
        return {
            "display_name": "renal",
            "tier_rules": {},
            "special_gene_lists": {},
            "tissue_genes": {"PKD1", "PKD2", "COL4A5"},
        }


class MockTissueAssessment:
    """Build tissue_assessment dicts for classify_variant_tier."""

    @staticmethod
    def primary(gtex_tpm: float = 50.0):
        return {
            "relevance": "primary",
            "gtex_tpm": gtex_tpm,
            "fast_track": False,
            "tier_suggestion": None,
            "reason": "Primary tissue gene",
        }

    @staticmethod
    def secondary(gtex_tpm: float = 10.0):
        return {
            "relevance": "secondary",
            "gtex_tpm": gtex_tpm,
            "fast_track": False,
            "tier_suggestion": None,
            "reason": "Secondary tissue gene",
        }

    @staticmethod
    def none(gtex_tpm: float = 0.0):
        return {
            "relevance": "none",
            "gtex_tpm": gtex_tpm,
            "fast_track": True,
            "tier_suggestion": 3,
            "reason": "No tissue relevance",
        }


def make_variant(**kwargs) -> Any:
    """Factory: create a dgra_core.Variant with sensible defaults."""
    from dgra_core import Variant
    defaults = {
        "chrom": "1",
        "pos": 100000,
        "ref": "A",
        "alt": "G",
        "gene": "TP53",
        "transcript": "ENST00000269305",
        "exon": "5/11",
        "impact": "HIGH",
        "consequence": "stop_gained",
        "hgvsp": "p.Arg273Ter",
        "hgvsc": "c.818C>T",
        "clinvar": "",
        "gnomad_af": None,
        "dp": 50,
        "gq": 99.0,
        "gt": "0/1",
        "vaf": 0.45,
    }
    defaults.update(kwargs)
    return Variant(**defaults)


def run_tests(module_name: str, test_funcs: List):
    """Simple test runner — no pytest dependency (legacy support)."""
    passed = 0
    failed = 0
    errors = []
    for fn in test_funcs:
        name = getattr(fn, "__name__", str(fn))
        try:
            import asyncio
            if asyncio.iscoroutinefunction(fn):
                asyncio.run(fn())
            else:
                fn()
            passed += 1
            print(f"  ✅ {name}")
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"  ❌ {name}: {e}")
    print(f"\n{module_name}: {passed} passed, {failed} failed")
    if errors:
        for name, err in errors:
            print(f"    {name}: {err}")
    return passed, failed
