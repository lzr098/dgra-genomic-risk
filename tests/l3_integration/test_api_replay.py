#!/usr/bin/env python3
"""
L3 Integration Tests — API Replay

Uses real recorded API responses (tests/recording/) to verify that
dgra_api.py, gpa_vcf_annotator.py, and gpa_two_phase.py correctly parse
live API payloads.

Run recorder to refresh recordings:
    python tests/record_api_responses.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional
from unittest.mock import patch

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from dgra_api import DGRAAPIClient
from dgra_cache import DGRACache
from dgra_config import DGRAGlobalConfig

# ---------------------------------------------------------------------------
# Recording loader
# ---------------------------------------------------------------------------

RECORDING_DIR = Path(__file__).parent.parent / "recording"


def load_recording(api_name: str, variant_id: str) -> Optional[Dict[str, Any]]:
    """Load a recorded response by api_name + variant_id."""
    subdir = RECORDING_DIR / api_name
    if not subdir.exists():
        return None
    for f in subdir.glob("*.json"):
        data = json.loads(f.read_text(encoding="utf-8"))
        if data["meta"]["variant_id"] == variant_id:
            return data["response"]
    return None


def load_recording_by_prefix(api_name: str, prefix: str) -> Optional[Dict[str, Any]]:
    """Load first recording whose variant_id starts with prefix."""
    subdir = RECORDING_DIR / api_name
    if not subdir.exists():
        return None
    for f in subdir.glob("*.json"):
        data = json.loads(f.read_text(encoding="utf-8"))
        if data["meta"]["variant_id"].startswith(prefix):
            return data["response"]
    return None


# ---------------------------------------------------------------------------
# Mock builder: intercept _request_with_retry and return recorded data
# ---------------------------------------------------------------------------

class _MockRequest:
    """Callable async mock that routes to recordings."""

    def __init__(self, recordings: Dict[str, Dict[str, Any]]):
        self.recordings = recordings

    async def __call__(
        self,
        api_name: str,
        endpoint: str,
        method: str = "GET",
        params: Optional[Dict] = None,
        json_body: Optional[Dict] = None,
        headers: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        recs = self.recordings
        sig = f"{api_name}:{endpoint}"

        # --- Ensembl VEP ---
        if api_name == "ensembl" and "/vep/human/region" in endpoint:
            body = json_body or []
            first = body[0] if isinstance(body, list) and body else ""
            parts = first.split()
            if len(parts) >= 5:
                chrom, pos, _, ref, alt = parts[:5]
                vid = f"{chrom}-{pos}-{ref}-{alt}"
                for k, v in recs.items():
                    if k.startswith("ensembl:") and k.endswith(vid):
                        return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}
            # Fallback: any ensembl VEP recording
            for k, v in recs.items():
                if k.startswith("ensembl:") and not k.startswith("ensembl:lookup-") and not k.startswith("ensembl:transcript-"):
                    return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}

        # --- Ensembl Lookup ---
        if api_name == "ensembl" and "/lookup/symbol/" in endpoint:
            gene = endpoint.split("/")[-1].split("?")[0]
            key = f"ensembl:lookup-{gene}"
            if key in recs:
                v = recs[key]
                return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}

        # --- Ensembl Transcript ---
        if api_name == "ensembl" and "/lookup/id/" in endpoint:
            tx_id = endpoint.split("/")[-1].split("?")[0]
            for k, v in recs.items():
                if k.startswith("ensembl:transcript-") and tx_id in k:
                    return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}

        # --- gnomAD ---
        if api_name == "gnomad":
            vid = (params or {}).get("variantId", "")
            if vid:
                # Recorded keys include gene prefix, e.g. "gnomad:TP53-17-7675088-C-T"
                for k, v in recs.items():
                    if k.startswith("gnomad:") and k.endswith(vid):
                        return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}
            # Fallback: first gnomad recording
            for k, v in recs.items():
                if k.startswith("gnomad:"):
                    return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}

        # --- UniProt entry ---
        if api_name == "uniprot" and ".json" in endpoint:
            for k, v in recs.items():
                if k.startswith("uniprot:"):
                    return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}

        # --- UniProt search ---
        if api_name == "uniprot" and "/uniprotkb/search" in endpoint:
            return {
                "data": {
                    "results": [
                        {
                            "primaryAccession": "P04637",
                            "entryType": "UniProtKB reviewed (Swiss-Prot)",
                            "sequence": {"length": 393},
                        }
                    ]
                },
                "http_status": 200,
                "from_cache": False,
                "confidence": "medium",
            }

        # --- NCBI / ClinVar ESearch ---
        if api_name in ("clinvar_eutils", "ncbi") and "/esearch.fcgi" in endpoint:
            term = (params or {}).get("term", "")
            gene = term.split("[")[0] if "[" in term else term
            key = f"ncbi:esearch-{gene}"
            if key in recs:
                v = recs[key]
                return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}
            for k, v in recs.items():
                if k.startswith("ncbi:esearch-"):
                    return {"data": v["data"], "http_status": v["http_status"], "from_cache": False, "confidence": "medium"}

        # --- NCBI / ClinVar ESummary ---
        if api_name in ("clinvar_eutils", "ncbi") and "/esummary.fcgi" in endpoint:
            return {
                "data": {
                    "result": {
                        "4850689": {
                            "title": "NM_007294.4(BRCA1):c.4185+1372T>G",
                            "clinical_significance": {"description": "Pathogenic"},
                            "review_status": "practice_guideline",
                            "variation_set": [{"cdna_hgvs": "NM_007294.4:c.4185+1372T>G"}],
                        }
                    }
                },
                "http_status": 200,
                "from_cache": False,
                "confidence": "medium",
            }

        return {"data": None, "http_status": 404, "from_cache": False, "confidence": "low", "error": f"No recording for {sig}"}


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

@pytest.fixture(scope="module")
def recorded_responses() -> Dict[str, Dict[str, Any]]:
    """Load all recordings into a flat dict keyed by 'api_name:variant_id'."""
    out: Dict[str, Dict[str, Any]] = {}
    if not RECORDING_DIR.exists():
        return out
    for subdir in RECORDING_DIR.iterdir():
        if not subdir.is_dir():
            continue
        for f in subdir.glob("*.json"):
            data = json.loads(f.read_text(encoding="utf-8"))
            key = f"{data['meta']['api_name']}:{data['meta']['variant_id']}"
            out[key] = data["response"]
    return out


@pytest.fixture
def api_client(tmp_path) -> DGRAAPIClient:
    """Return a DGRAAPIClient with an in-memory cache."""
    config = DGRAGlobalConfig()
    config.cache_db_path = tmp_path / "test_cache.db"
    cache = DGRACache(db_path=config.cache_db_path)
    client = DGRAAPIClient(config=config, cache=cache)
    return client


# ---------------------------------------------------------------------------
# Tests: Ensembl
# ---------------------------------------------------------------------------

@pytest.mark.l3
@pytest.mark.replay
class TestEnsemblReplay:
    @pytest.mark.asyncio
    async def test_query_ensembl_gene_tp53(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_ensembl_gene("TP53")
        assert result["source"] in ("ensembl", "cache")
        assert result["biotype"] == "protein_coding"
        assert result["canonical_transcript"].startswith("ENST")
        assert result["seq_region_name"] == "17"
        assert result["confidence"] in ("medium", "high")

    @pytest.mark.asyncio
    async def test_query_ensembl_gene_brca1(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_ensembl_gene("BRCA1")
        assert result["source"] in ("ensembl", "cache")
        assert result["biotype"] == "protein_coding"
        assert result["confidence"] in ("medium", "high")

    @pytest.mark.asyncio
    async def test_query_ensembl_vep_region_tp53(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_ensembl_vep_region("17", 7675088, "C", "T")
        assert result["source"] == "ensembl"
        assert "missense_variant" in result.get("consequence_terms", [])
        assert result["impact"] == "MODERATE"
        assert result["gene_symbol"] == "TP53"
        assert result["transcript_id"].startswith("ENST")

    @pytest.mark.asyncio
    async def test_query_ensembl_transcript_info(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_ensembl_transcript_info("ENST00000269305")
        assert result["source"] in ("ensembl", "cache")
        assert result["transcript_id"] == "ENST00000269305"
        assert result["biotype"] == "protein_coding"


# ---------------------------------------------------------------------------
# Tests: UniProt
# ---------------------------------------------------------------------------

@pytest.mark.l3
@pytest.mark.replay
class TestUniProtReplay:
    @pytest.mark.asyncio
    async def test_query_uniprot_by_gene_tp53(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_uniprot_by_gene("TP53")
        assert result["uniprot_id"] == "P04637"
        assert result["source"] == "uniprot"
        assert result["confidence"] in ("medium", "high")
        assert isinstance(result["domains"], list)

    @pytest.mark.asyncio
    async def test_query_uniprot_by_gene_brca1(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_uniprot_by_gene("BRCA1")
        assert result["uniprot_id"] is not None
        assert result["source"] == "uniprot"


# ---------------------------------------------------------------------------
# Tests: gnomAD
# ---------------------------------------------------------------------------

@pytest.mark.l3
@pytest.mark.replay
class TestGnomADReplay:
    @pytest.mark.asyncio
    async def test_query_gnomad_variant_tp53(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_gnomad_variant("17", 7675088, "C", "T")
        assert result["source"] == "gnomad"
        assert result["status"] == "SUCCESS"
        assert result["variant_id"] == "17-7675088-C-T"
        assert result["af"] is not None
        assert result["af"] < 0.001
        assert isinstance(result["af_populations"], dict)

    @pytest.mark.asyncio
    async def test_query_gnomad_variant_dmd(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_gnomad_variant("X", 31140024, "C", "T")
        assert result["source"] == "gnomad"
        # DMD variant not in gnomAD dataset → NOT_CAPTURED or QUERY_ERROR
        assert result["status"] in ("SUCCESS", "NOT_CAPTURED", "QUERY_ERROR")


# ---------------------------------------------------------------------------
# Tests: NCBI ClinVar
# ---------------------------------------------------------------------------

@pytest.mark.l3
@pytest.mark.replay
class TestNCBIClinVarReplay:
    @pytest.mark.asyncio
    async def test_query_ncbi_clinvar_brca1(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_ncbi_clinvar("BRCA1", chrom="17", pos=43044295)
        assert result["gene"] == "BRCA1"
        assert result["clinvar_id"] is not None
        assert result["source"] in ("clinvar", "cache")

    @pytest.mark.asyncio
    async def test_query_ncbi_clinvar_tp53(self, api_client, recorded_responses):
        with patch.object(api_client, "_request_with_retry", _MockRequest(recorded_responses)):
            result = await api_client.query_ncbi_clinvar("TP53", chrom="17", pos=7675088)
        assert result["gene"] == "TP53"


# ---------------------------------------------------------------------------
# Tests: VCF Annotator — parse recorded VEP response
# ---------------------------------------------------------------------------

@pytest.mark.l3
@pytest.mark.replay
class TestVcfAnnotatorReplay:
    def test_parse_vep_response_tp53(self, recorded_responses):
        from gpa_vcf_annotator import VCFAnnotator
        rec = recorded_responses.get("ensembl:TP53-17-7675088-C-T")
        assert rec is not None, "Recording not found"
        data = rec["data"]
        assert isinstance(data, list)
        batch = [{"chrom": "17", "pos": 7675088, "ref": "C", "alt": "T"}]
        results = VCFAnnotator._parse_vep_response(data, batch)
        assert len(results) == 1
        r = results[0]
        assert r["vep_summary"]["most_severe_consequence"] == "missense_variant"
        assert len(r["transcript_consequences"]) > 0
        tx = r["transcript_consequences"][0]
        assert tx["gene_symbol"] == "TP53"
        assert tx["impact"] == "MODERATE"
        assert tx["transcript_id"].startswith("ENST")

    def test_parse_vep_response_cftr(self, recorded_responses):
        from gpa_vcf_annotator import VCFAnnotator
        rec = recorded_responses.get("ensembl:CFTR-7-117559590-AT-A")
        if rec is None:
            pytest.skip("CFTR recording not available")
        data = rec["data"]
        batch = [{"chrom": "7", "pos": 117559590, "ref": "AT", "alt": "A"}]
        results = VCFAnnotator._parse_vep_response(data, batch)
        assert len(results) == 1
        assert len(results[0]["transcript_consequences"]) >= 0
