#!/usr/bin/env python3
"""
Record real API responses for GPA test suite replay.

Usage:
    cd /Users/zhaorongli/.workbuddy/skills/dgra-genomic-risk
    python tests/record_api_responses.py

Requires: aiohttp (pip install aiohttp)
Network:  Uses proxy 127.0.0.1:7897 if available, else direct.
"""

import asyncio
import hashlib
import json
import sys
from pathlib import Path
from typing import Any, Dict, Optional

import aiohttp

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

RECORDING_DIR = Path(__file__).parent / "recording"
PROXY = "http://127.0.0.1:7897"

# Variants to record (GRCh38)
VARIANTS = [
    {"id": "TP53-17-7675088-C-T",   "chrom": "17", "pos": 7675088,  "ref": "C", "alt": "T", "gene": "TP53"},
    {"id": "BRCA1-17-43044295-T-G", "chrom": "17", "pos": 43044295, "ref": "T", "alt": "G", "gene": "BRCA1"},
    {"id": "DMD-X-31140024-C-T",    "chrom": "X",  "pos": 31140024, "ref": "C", "alt": "T", "gene": "DMD"},
    {"id": "CFTR-7-117559590-AT-A", "chrom": "7",  "pos": 117559590,"ref": "AT","alt": "A","gene": "CFTR"},
]

GENES = ["TP53", "BRCA1", "DMD"]

UNIPROT_IDS = {
    "TP53":  "P04637",
    "BRCA1": "P38398",
    "DMD":   "P11532",
}

# ---------------------------------------------------------------------------
# Recorder
# ---------------------------------------------------------------------------

class ResponseRecorder:
    def __init__(self, recording_dir: Path):
        self.recording_dir = recording_dir
        self.index: Dict[str, str] = {}
        self._load_index()

    def _load_index(self):
        index_path = self.recording_dir / ".index.json"
        if index_path.exists():
            self.index = json.loads(index_path.read_text(encoding="utf-8"))

    def save_index(self):
        index_path = self.recording_dir / ".index.json"
        index_path.write_text(json.dumps(self.index, indent=2, ensure_ascii=False), encoding="utf-8")

    def _recording_path(self, api_name: str, key: str) -> Path:
        safe_key = hashlib.md5(key.encode()).hexdigest()[:12]
        subdir = self.recording_dir / api_name
        subdir.mkdir(parents=True, exist_ok=True)
        return subdir / f"{safe_key}.json"

    def save(self, api_name: str, key: str, response: Dict[str, Any]):
        path = self._recording_path(api_name, key)
        wrapped = {
            "meta": {
                "api_name": api_name,
                "variant_id": key,
                "recorded_at": _iso_now(),
                "mode": "record",
            },
            "response": response,
        }
        path.write_text(json.dumps(wrapped, indent=2, ensure_ascii=False), encoding="utf-8")
        rel = str(path.relative_to(self.recording_dir))
        self.index[key] = rel
        print(f"  [SAVED] {api_name}/{key} -> {rel}")


def _iso_now() -> str:
    from datetime import datetime, timezone
    return datetime.now(timezone.utc).isoformat()


# ---------------------------------------------------------------------------
# API callers
# ---------------------------------------------------------------------------

async def _session(proxy: Optional[str] = None):
    return aiohttp.ClientSession(
        trust_env=False,
        timeout=aiohttp.ClientTimeout(total=30),
    )


async def record_ensembl_vep(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[Ensembl VEP]")
    async with await _session(proxy) as session:
        for v in VARIANTS:
            key = v["id"]
            body = {
                "variants": [f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ."],
                "canonical": "1",
                "mane_select": "1",
                "domains": "1",
                "protein": "1",
                "hgvs": "1",
                "numbers": "1",
            }
            try:
                async with session.post(
                    "https://rest.ensembl.org/vep/human/region",
                    headers={"Content-Type": "application/json"},
                    json=body,
                    proxy=proxy,
                ) as resp:
                    data = await resp.json()
                    recorder.save("ensembl", key, {
                        "http_status": resp.status,
                        "data": data,
                    })
            except Exception as e:
                print(f"  [FAIL] ensembl/{key}: {e}")


async def record_gnomad(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[gnomAD]")
    async with await _session(proxy) as session:
        for v in VARIANTS:
            key = v["id"]
            query = (
                f'{{ variant(variantId: "{v["chrom"]}-{v["pos"]}-{v["ref"]}-{v["alt"]}"'
                f', dataset: gnomad_r4) '
                f'{{ variantId exome {{ af ac an }} genome {{ af ac an }} }} }}'
            )
            try:
                async with session.post(
                    "https://gnomad.broadinstitute.org/api/",
                    headers={"Content-Type": "application/json"},
                    json={"query": query},
                    proxy=proxy,
                ) as resp:
                    data = await resp.json()
                    recorder.save("gnomad", key, {
                        "http_status": resp.status,
                        "data": data,
                    })
            except Exception as e:
                print(f"  [FAIL] gnomad/{key}: {e}")


async def record_uniprot(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[UniProt]")
    async with await _session(proxy) as session:
        for gene, acc in UNIPROT_IDS.items():
            key = gene
            try:
                async with session.get(
                    f"https://rest.uniprot.org/uniprotkb/{acc}.json",
                    proxy=proxy,
                ) as resp:
                    data = await resp.json()
                    recorder.save("uniprot", key, {
                        "http_status": resp.status,
                        "data": data,
                    })
            except Exception as e:
                print(f"  [FAIL] uniprot/{key}: {e}")


async def record_ncbi_esearch(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[NCBI ESearch]")
    async with await _session(proxy) as session:
        for gene in GENES:
            key = gene
            url = (
                "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
                f"?db=clinvar&term={gene}%5BGene%5D&retmode=json&retmax=5"
            )
            try:
                async with session.get(url, proxy=proxy) as resp:
                    data = await resp.json()
                    recorder.save("ncbi", f"esearch-{key}", {
                        "http_status": resp.status,
                        "data": data,
                    })
            except Exception as e:
                print(f"  [FAIL] ncbi/esearch-{key}: {e}")


async def record_ncbi_efetch(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[NCBI EFetch]")
    await asyncio.sleep(2)  # Rate-limit buffer between esearch and efetch
    async with await _session(proxy) as session:
        # Fetch first ClinVar record for BRCA1
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            "?db=clinvar&id=4850689&rettype=vcv&is_variationid"
        )
        try:
            async with session.get(url, proxy=proxy) as resp:
                text = await resp.text()
                recorder.save("ncbi", "efetch-BRCA1-4850689", {
                    "http_status": resp.status,
                    "data": text,
                })
        except Exception as e:
            print(f"  [FAIL] ncbi/efetch-BRCA1-4850689: {e}")


async def record_ensembl_lookup(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[Ensembl Gene Lookup]")
    async with await _session(proxy) as session:
        for gene in GENES:
            key = gene
            url = (
                f"https://rest.ensembl.org/lookup/symbol/homo_sapiens/{gene}"
                "?expand=1&content-type=application/json"
            )
            try:
                async with session.get(url, proxy=proxy) as resp:
                    data = await resp.json()
                    recorder.save("ensembl", f"lookup-{key}", {
                        "http_status": resp.status,
                        "data": data,
                    })
            except Exception as e:
                print(f"  [FAIL] ensembl/lookup-{key}: {e}")


async def record_ensembl_transcript(recorder: ResponseRecorder, proxy: Optional[str]):
    print("\n[Ensembl Transcript Lookup]")
    transcripts = {
        "TP53":  "ENST00000269305",
        "BRCA1": "ENST00000357654",
    }
    async with await _session(proxy) as session:
        for gene, tx in transcripts.items():
            key = f"{gene}-{tx}"
            url = (
                f"https://rest.ensembl.org/lookup/id/{tx}"
                "?expand=1&content-type=application/json"
            )
            try:
                async with session.get(url, proxy=proxy) as resp:
                    data = await resp.json()
                    recorder.save("ensembl", f"transcript-{key}", {
                        "http_status": resp.status,
                        "data": data,
                    })
            except Exception as e:
                print(f"  [FAIL] ensembl/transcript-{key}: {e}")


async def record_gtex(recorder: ResponseRecorder, proxy: Optional[str]):
    """GTEx v2 API — try multiple endpoints."""
    print("\n[GTEx]")
    async with await _session(proxy) as session:
        for gene in GENES:
            key = gene
            urls = [
                f"https://gtexportal.org/api/v2/reference/gene?geneId={gene}&format=json",
                f"https://gtexportal.org/rest/v1/reference/gene?geneId={gene}&format=json",
                f"https://gtexportal.org/rest/v1/reference/geneSearch?geneId={gene}",
            ]
            for url in urls:
                try:
                    async with session.get(url, proxy=proxy, timeout=aiohttp.ClientTimeout(total=10)) as resp:
                        data = await resp.json()
                        if data:
                            recorder.save("gtex", key, {
                                "http_status": resp.status,
                                "data": data,
                            })
                            break
                except Exception:
                    continue
            else:
                print(f"  [SKIP] gtex/{key}: all endpoints empty or failed")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

async def main():
    print("=" * 60)
    print("GPA API Response Recorder")
    print("=" * 60)

    recorder = ResponseRecorder(RECORDING_DIR)

    # Probe proxy
    proxy = None
    try:
        async with await _session() as session:
            async with session.get(
                "https://rest.ensembl.org/info/ping?content-type=application/json",
                proxy=PROXY,
                timeout=aiohttp.ClientTimeout(total=5),
            ) as resp:
                if resp.status == 200:
                    proxy = PROXY
                    print(f"[Proxy] Using {PROXY}")
    except Exception:
        print("[Proxy] Direct connection (no proxy)")

    await record_ensembl_vep(recorder, proxy)
    await record_gnomad(recorder, proxy)
    await record_uniprot(recorder, proxy)
    await record_ncbi_esearch(recorder, proxy)
    await record_ncbi_efetch(recorder, proxy)
    await record_ensembl_lookup(recorder, proxy)
    await record_ensembl_transcript(recorder, proxy)
    await record_gtex(recorder, proxy)

    recorder.save_index()
    print("\n[Done] Index saved to .index.json")
    print(f"[Done] Recordings in {RECORDING_DIR}")


if __name__ == "__main__":
    asyncio.run(main())
