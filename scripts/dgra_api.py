#!/usr/bin/env python3
"""
DGRA API Query Layer
Phase 1 - v0.4 Architecture

Async API wrappers for Ensembl, UniProt, GTEx, gnomAD, and NCBI E-utilities.
Each wrapper implements: query -> cache check -> API call -> cache store -> return.
All functions are async and use aiohttp for concurrent requests.

Skeleton implementation: defines interfaces, implements basic HTTP logic,
leaves full response parsing for Phase 2.
"""

import asyncio
import aiohttp
import hashlib
import json
import time
from typing import Optional, Dict, Any, List
from pathlib import Path

from dgra_config import DGRAGlobalConfig, APIConfig
from dgra_cache import DGRACache


class DGRAAPIError(Exception):
    """Base exception for API errors."""
    def __init__(self, api_name: str, message: str, status: Optional[int] = None, 
                 response: Optional[str] = None):
        self.api_name = api_name
        self.status = status
        self.response = response
        super().__init__(f"[{api_name}] {message}")


class DGRAAPIClient:
    """
    Unified async API client for all DGRA external data sources.

    Handles:
    - Cache lookup before API call
    - Rate limiting (per-API token bucket)
    - Retry with exponential backoff
    - Response caching on success
    - Offline mode fallback (skip API, return cached or None)
    - Dynamic proxy auto-detection (v0.10.5)
    """

    # Common proxy endpoints to probe (in priority order)
    _COMMON_PROXIES: List[str] = [
        "http://127.0.0.1:7897",
        "http://127.0.0.1:7890",
        "http://127.0.0.1:7891",
        "http://127.0.0.1:1080",
        "http://127.0.0.1:10808",
        "http://127.0.0.1:10809",
        "http://127.0.0.1:8000",
        "http://127.0.0.1:8080",
        "http://127.0.0.1:3128",
    ]

    def __init__(self, config: DGRAGlobalConfig, cache: DGRACache, proxy_route_map: Optional[Any] = None):
        self.config = config
        self.cache = cache
        self._proxy_route_map = proxy_route_map
        # v0.10.12: fallback to config-attached route map (set by gpa_two_phase pipeline)
        if self._proxy_route_map is None and hasattr(config, '_proxy_route_map'):
            self._proxy_route_map = config._proxy_route_map
        self._session: Optional[aiohttp.ClientSession] = None
        self._last_request_time: Dict[str, float] = {}  # api_name -> timestamp
        self._proxy_url: Optional[str] = None  # detected working proxy or None for direct

    @staticmethod
    async def _probe_endpoint(proxy: Optional[str], timeout: float = 3.0) -> bool:
        """Quickly probe NCBI esearch to test if a route works."""
        test_url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
            "?db=clinvar&term=BRCA1%5BGene%5D&retmode=json&retmax=1"
        )
        try:
            session = aiohttp.ClientSession(
                trust_env=False,
                timeout=aiohttp.ClientTimeout(total=timeout),
            )
            try:
                async with session.get(test_url, proxy=proxy) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        count = data.get("esearchresult", {}).get("count")
                        return bool(count and int(count) > 0)
            finally:
                await session.close()
        except Exception:
            return False
        return False

    async def _detect_proxy(self) -> Optional[str]:
        """
        Auto-detect the best network route for scientific APIs.

        Returns:
            Proxy URL (e.g. "http://127.0.0.1:7897") if a proxy works,
            None if direct connection works best.
        """
        # 1. Try direct connection first (fastest if it works)
        print("[ProxyProbe] Testing direct connection to NCBI...")
        if await self._probe_endpoint(proxy=None, timeout=3.0):
            print("[ProxyProbe] Direct connection OK — no proxy needed")
            return None
        print("[ProxyProbe] Direct connection FAILED")

        # 2. Try each common proxy in parallel (limited concurrency)
        semaphore = asyncio.Semaphore(3)
        results: Dict[str, bool] = {}

        async def _test_one(proxy: str) -> None:
            async with semaphore:
                ok = await self._probe_endpoint(proxy=proxy, timeout=5.0)
                results[proxy] = ok
                if ok:
                    print(f"[ProxyProbe] Proxy {proxy} OK")
                else:
                    print(f"[ProxyProbe] Proxy {proxy} FAILED")

        tasks = [asyncio.create_task(_test_one(p)) for p in self._COMMON_PROXIES]
        await asyncio.gather(*tasks)

        # Return first working proxy
        for proxy in self._COMMON_PROXIES:
            if results.get(proxy):
                return proxy

        print("[ProxyProbe] No working proxy found — will use direct (expect failures)")
        return None

    async def __aenter__(self):
        # v0.10.5: Dynamic proxy auto-detection
        # Probe once at session start, then reuse the route for all APIs.
        self._proxy_url = await self._detect_proxy()
        self._session = aiohttp.ClientSession(
            connector=aiohttp.TCPConnector(limit=50, limit_per_host=20),
            timeout=aiohttp.ClientTimeout(total=120),
            trust_env=False,  # Never trust env; use our detected proxy explicitly
        )
        return self
    
    async def __aexit__(self, exc_type, exc_val, exc_tb):
        if self._session:
            await self._session.close()
            self._session = None
    
    async def _rate_limit(self, api_name: str):
        """Enforce per-API rate limit using token bucket logic."""
        cfg = self.config.apis[api_name]
        min_interval = 1.0 / cfg.rate_limit_per_sec
        now = time.time()
        last = self._last_request_time.get(api_name, 0)
        elapsed = now - last
        if elapsed < min_interval:
            await asyncio.sleep(min_interval - elapsed)
        self._last_request_time[api_name] = time.time()
    
    async def _request_with_retry(
        self, 
        api_name: str,
        endpoint: str,
        method: str = "GET",
        params: Optional[Dict] = None,
        json_body: Optional[Dict] = None,
        headers: Optional[Dict] = None,
    ) -> Dict[str, Any]:
        """
        Execute HTTP request with cache check, rate limiting, and retry.
        
        Returns dict with keys:
        - data: parsed JSON response
        - http_status: HTTP status code
        - from_cache: bool
        - confidence: 'high' if cache hit, 'medium' if API success, 'low' if partial
        """
        cfg = self.config.apis[api_name]
        url = f"{cfg.base_url}/{endpoint.lstrip('/')}"
        
        # Phase 1: Check cache (skip if offline mode - we already checked before calling)
        # v0.12.2 FIX: include json_body in cache key for POST/GraphQL requests.
        # Without this, all POSTs to the same endpoint share one cache key,
        # causing cache collisions (e.g. all gnomAD GraphQL queries overwrite each other).
        cache_key_params = {"url": url, **(params or {})}
        if json_body is not None:
            # Use a stable, compact representation of the JSON body
            body_key = json.dumps(json_body, sort_keys=True, separators=(',', ':'))
            cache_key_params["_body"] = body_key
        cached = self.cache.get(api_name, **cache_key_params)
        
        if cached:
            return {
                "data": cached["data"],
                "http_status": cached["http_status"],
                "from_cache": True,
                "confidence": cached["confidence"],
            }
        
        # Offline mode: no cache hit = return None
        if self.config.offline_mode:
            return {
                "data": None,
                "http_status": None,
                "from_cache": False,
                "confidence": "low",
                "error": "Offline mode: no cached data available",
            }
        
        # Phase 2: API call with retry
        last_error = None
        for attempt in range(cfg.max_retries):
            try:
                await self._rate_limit(api_name)

                # v0.10.12: per-API proxy routing
                proxy = self._proxy_url
                if self._proxy_route_map is not None:
                    proxy = self._proxy_route_map.get_proxy(api_name)

                async with self._session.request(
                    method=method,
                    url=url,
                    params=params,
                    json=json_body,
                    headers=headers,
                    proxy=proxy,
                    timeout=aiohttp.ClientTimeout(total=cfg.timeout),
                ) as response:
                    http_status = response.status
                    
                    if http_status == 200:
                        try:
                            data = await response.json()
                        except Exception:
                            # Try to parse text as JSON fallback
                            text = await response.text()
                            try:
                                data = json.loads(text)
                            except Exception:
                                # Not valid JSON — don't cache, return as error
                                return {
                                    "data": None,
                                    "http_status": http_status,
                                    "from_cache": False,
                                    "confidence": "low",
                                    "error": f"HTTP 200 but response is not valid JSON (len={len(text)})",
                                }
                        
                        # Cache successful response
                        self.cache.set(
                            api_name=api_name,
                            response_data=data,
                            http_status=http_status,
                            confidence="medium",
                            **cache_key_params
                        )
                        
                        return {
                            "data": data,
                            "http_status": http_status,
                            "from_cache": False,
                            "confidence": "medium",
                        }
                    
                    elif http_status == 404:
                        # Not found - cache the negative result with shorter TTL
                        self.cache.set(
                            api_name=api_name,
                            response_data={"error": "not_found", "status": 404},
                            http_status=404,
                            confidence="medium",
                            ttl_days=7,  # Shorter TTL for negatives
                            **cache_key_params
                        )
                        return {
                            "data": None,
                            "http_status": 404,
                            "from_cache": False,
                            "confidence": "medium",
                            "error": "Not found",
                        }
                    
                    elif http_status == 429:
                        # Rate limited - read Retry-After header and wait
                        retry_after = response.headers.get("Retry-After")
                        if retry_after:
                            wait = int(retry_after)
                        else:
                            wait = cfg.retry_delay * (2 ** attempt)
                        last_error = DGRAAPIError(api_name, f"Rate limited (429), retry after {wait}s", http_status)
                        print(f"[DGRA API] {api_name}: Rate limited (429), waiting {wait}s before retry {attempt+1}/{cfg.max_retries}")
                        await asyncio.sleep(wait)
                        continue
                    
                    elif http_status in (502, 503, 504):
                        # Gateway/server temporarily unavailable - retry
                        last_error = DGRAAPIError(api_name, f"Server temporarily unavailable ({http_status})", http_status)
                        wait = cfg.retry_delay * (2 ** attempt)
                        print(f"[DGRA API] {api_name}: HTTP {http_status}, retrying in {wait}s (attempt {attempt+1}/{cfg.max_retries})")
                        await asyncio.sleep(wait)
                        continue
                    
                    elif http_status >= 500:
                        # Other server errors - retry with backoff
                        last_error = DGRAAPIError(api_name, f"Server error {http_status}", http_status)
                        wait = cfg.retry_delay * (2 ** attempt)
                        print(f"[DGRA API] {api_name}: HTTP {http_status}, retrying in {wait}s (attempt {attempt+1}/{cfg.max_retries})")
                        await asyncio.sleep(wait)
                        continue
                    
                    else:
                        # Client error or other - don't retry
                        text = await response.text()
                        return {
                            "data": None,
                            "http_status": http_status,
                            "from_cache": False,
                            "confidence": "low",
                            "error": f"HTTP {http_status}: {text[:200]}",
                        }
            
            except asyncio.TimeoutError:
                last_error = DGRAAPIError(api_name, f"Timeout after {cfg.timeout}s")
                wait = cfg.retry_delay * (2 ** attempt)
                print(f"[DGRA API] {api_name}: Timeout, retrying in {wait}s (attempt {attempt+1}/{cfg.max_retries})")
                await asyncio.sleep(wait)
                continue
            
            except aiohttp.ClientError as e:
                last_error = DGRAAPIError(api_name, f"Connection error: {e}")
                wait = cfg.retry_delay * (2 ** attempt)
                print(f"[DGRA API] {api_name}: Connection error ({e}), retrying in {wait}s (attempt {attempt+1}/{cfg.max_retries})")
                await asyncio.sleep(wait)
                continue
        
        # All retries exhausted
        return {
            "data": None,
            "http_status": last_error.status if last_error else None,
            "from_cache": False,
            "confidence": "low",
            "error": str(last_error) if last_error else "All retries failed",
        }
    
    # =====================================================================
    # Ensembl REST API
    # =====================================================================
    
    async def query_ensembl_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Query Ensembl for gene canonical transcript, biotype, and basic info.
        
        Endpoint: GET /lookup/symbol/homo_sapiens/{gene_symbol}?expand=1
        
        Returns:
        {
            "canonical_transcript": "ENST...",
            "biotype": "protein_coding",
            "description": "...",
            "seq_region_name": "chr...",
            "start": 12345,
            "end": 67890,
            "strand": 1,
            "source": "ensembl|cache",
            "confidence": "high|medium|low",
        }
        """
        result = await self._request_with_retry(
            api_name="ensembl",
            endpoint=f"/lookup/symbol/homo_sapiens/{gene_symbol}",
            params={"expand": "1", "content-type": "application/json"},
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            # Extract canonical transcript (usually the first or longest)
            transcripts = data.get("Transcript", [])
            canonical = None
            for tx in transcripts:
                if tx.get("is_canonical", 0) == 1:
                    canonical = tx["id"]
                    break
            if not canonical and transcripts:
                canonical = transcripts[0]["id"]  # Fallback to first
            
            return {
                "canonical_transcript": canonical,
                "biotype": data.get("biotype"),
                "description": data.get("description"),
                "seq_region_name": data.get("seq_region_name"),
                "start": data.get("start"),
                "end": data.get("end"),
                "strand": data.get("strand"),
                "source": "cache" if result["from_cache"] else "ensembl",
                "confidence": result["confidence"],
                "raw": data,  # Keep full response for Phase 2 parsing
            }
        
        return {
            "canonical_transcript": None,
            "biotype": None,
            "description": None,
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
            "raw": None,
        }
    
    async def query_ensembl_transcript_info(self, transcript_id: str) -> Dict[str, Any]:
        """
        Query Ensembl for transcript details (CDS, exons, translation).
        
        Endpoint: GET /lookup/id/{transcript_id}?expand=1
        """
        result = await self._request_with_retry(
            api_name="ensembl",
            endpoint=f"/lookup/id/{transcript_id}",
            params={"expand": "1"},
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            return {
                "transcript_id": data.get("id"),
                "display_name": data.get("display_name"),
                "biotype": data.get("biotype"),
                "cds_length": len(data.get("CDS", [])),
                "exon_count": len(data.get("Exon", [])),
                "translation_id": data.get("Translation", {}).get("id") if data.get("Translation") else None,
                "source": "cache" if result["from_cache"] else "ensembl",
                "confidence": result["confidence"],
                "raw": data,
            }
        
        return {
            "transcript_id": transcript_id,
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
            "raw": None,
        }
    

    async def query_ensembl_vep_region(self, chrom: str, pos: int, ref: str, alt: str) -> Dict[str, Any]:
        """
        Query Ensembl VEP for canonical transcript annotation of a specific variant.

        Endpoint: POST /vep/human/region
        Body: ["{chrom} {pos} . {ref} {alt} . . ."]
        Params: canonical=1&domains=1&protein=1&hgvs=1&mane_select=1

        Parses VEP JSON response to extract canonical transcript's:
        - consequence_terms
        - impact
        - hgvsc
        - hgvsp
        - transcript_id
        - protein_domains

        Returns structured dict with these fields.
        """
        variant_string = f"{chrom} {pos} . {ref} {alt} . . ."
        # Include body hash in params for cache key uniqueness (server ignores unknown params)
        body_hash = hashlib.md5(json.dumps([variant_string], sort_keys=True).encode()).hexdigest()[:16]
        params = {
            "canonical": "1",
            "domains": "1",
            "protein": "1",
            "hgvs": "1",
            "mane_select": "1",
            "_body_hash": body_hash,
        }
        result = await self._request_with_retry(
            api_name="ensembl",
            endpoint="/vep/human/region",
            method="POST",
            json_body=[variant_string],
            params=params,
        )

        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            if not data or not isinstance(data, list):
                return {
                    "consequence_terms": [],
                    "impact": None,
                    "hgvsc": None,
                    "hgvsp": None,
                    "transcript_id": None,
                    "protein_domains": [],
                    "source": "failed",
                    "confidence": "low",
                    "error": "Invalid VEP response format (expected list)",
                }

            parsed = self._parse_vep_batch_response(data, [{"chrom": chrom, "pos": pos, "ref": ref, "alt": alt}])
            if parsed:
                return parsed[0]

        return {
            "consequence_terms": [],
            "impact": None,
            "hgvsc": None,
            "hgvsp": None,
            "transcript_id": None,
            "protein_domains": [],
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
        }

    def _parse_vep_batch_response(self, data: List[Dict], variants: List[Dict]) -> List[Dict]:
        """Parse VEP batch response, matching results to input variants.

        For each variant result, finds the canonical transcript consequence
        (prioritizing: canonical=1 > MANE Select > protein_coding > first).
        """
        results = []
        for idx, variant_result in enumerate(data):
            if not isinstance(variant_result, dict):
                results.append({
                    "error": "Invalid VEP response format",
                    "source": "failed",
                    "confidence": "low",
                })
                continue

            tx_consequences = variant_result.get("transcript_consequences", [])
            canonical_tx = None
            # Priority 1: canonical=1
            for tx in tx_consequences:
                if tx.get("canonical") == 1:
                    canonical_tx = tx
                    break
            # Priority 2: MANE Select
            if not canonical_tx:
                for tx in tx_consequences:
                    if tx.get("mane_select"):
                        canonical_tx = tx
                        break
            # Priority 3: protein_coding biotype
            if not canonical_tx:
                for tx in tx_consequences:
                    if tx.get("biotype") == "protein_coding":
                        canonical_tx = tx
                        break
            # Priority 4: first available
            if not canonical_tx and tx_consequences:
                canonical_tx = tx_consequences[0]

            if canonical_tx:
                # Parse protein_domains
                # VEP returns protein_domains as list of strings: "Db:ID:Name" or dicts
                domains = []
                raw_domains = canonical_tx.get("protein_domains", [])
                for d in raw_domains:
                    if isinstance(d, dict):
                        domains.append({
                            "name": d.get("name", d.get("description", "unnamed")),
                            "start": d.get("start", d.get("beg")),
                            "end": d.get("end"),
                            "db": d.get("db"),
                        })
                    elif isinstance(d, str):
                        parts = d.split(":")
                        if len(parts) >= 3:
                            domains.append({
                                "name": parts[2],
                                "db": parts[0],
                                "interpro_id": parts[1] if len(parts) > 1 else None,
                            })
                        elif len(parts) == 2:
                            domains.append({
                                "name": parts[1],
                                "db": parts[0],
                            })

                results.append({
                    "consequence_terms": canonical_tx.get("consequence_terms", []),
                    "impact": canonical_tx.get("impact"),
                    "hgvsc": canonical_tx.get("hgvsc"),
                    "hgvsp": canonical_tx.get("hgvsp"),
                    "transcript_id": canonical_tx.get("transcript_id"),
                    "gene_symbol": canonical_tx.get("gene_symbol"),
                    "protein_domains": domains,
                    "source": "ensembl",
                    "confidence": "medium",
                    "raw": variant_result,
                })
            else:
                results.append({
                    "error": "No transcript consequences found",
                    "source": "ensembl",
                    "confidence": "low",
                })

        return results

    async def batch_query_vep_region(self, variants: List[Dict[str, Any]]) -> Dict[str, Dict[str, Any]]:
        """
        Batch query Ensembl VEP for multiple variants with controlled concurrency.

        VEP API accepts multiple variants per POST request.
        Strategy: chunk variants (max 50 per request) + semaphore (max 5 concurrent).

        Args:
            variants: List of dicts with keys chrom, pos, ref, alt, and optionally key.

        Returns:
            Dict mapping variant key -> VEP result dict.
        """
        if not variants:
            return {}

        CHUNK_SIZE = 50
        MAX_CONCURRENT = 5
        semaphore = asyncio.Semaphore(MAX_CONCURRENT)

        async def _query_chunk(chunk: List[Dict]) -> List[Dict]:
            async with semaphore:
                body = [f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ." for v in chunk]
                # Body hash for cache key uniqueness
                body_hash = hashlib.md5(json.dumps(body, sort_keys=True).encode()).hexdigest()[:16]
                params = {
                    "canonical": "1",
                    "domains": "1",
                    "protein": "1",
                    "hgvs": "1",
                    "mane_select": "1",
                    "_body_hash": body_hash,
                }
                result = await self._request_with_retry(
                    api_name="ensembl",
                    endpoint="/vep/human/region",
                    method="POST",
                    json_body=body,
                    params=params,
                )
                if result["data"] and result["http_status"] == 200:
                    return self._parse_vep_batch_response(result["data"], chunk)
                else:
                    return [{"error": result.get("error"), "source": "failed", "confidence": "low"} for _ in chunk]

        all_results = {}
        total = len(variants)
        chunks = [variants[i:i + CHUNK_SIZE] for i in range(0, total, CHUNK_SIZE)]
        
        # v0.9.3: True concurrency — Semaphore(5) limits to 5 chunks at a time
        tasks = [_query_chunk(chunk) for chunk in chunks]
        all_chunk_results = await asyncio.gather(*tasks, return_exceptions=True)
        
        for chunk, chunk_results in zip(chunks, all_chunk_results):
            if isinstance(chunk_results, Exception):
                # Log error and continue — individual chunk failures don't stop the batch
                continue
            for v, r in zip(chunk, chunk_results):
                key = v.get("key", f"{v['chrom']}:{v['pos']}_{v['ref']}>{v['alt']}")
                all_results[key] = r

        return all_results

    # =====================================================================
    # UniProt REST API
    # =====================================================================
    
    async def query_uniprot_by_gene(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Query UniProt for protein entry by gene symbol.
        
        Step 1: Search gene -> uniprot ID mapping
        Step 2: Query /uniprotkb/{id}.json for full entry
        
        Returns:
        {
            "uniprot_id": "P12345",
            "protein_name": "...",
            "sequence_length": 1250,
            "domains": [
                {"name": "Motor domain", "start": 1, "end": 780, "type": "DOMAIN"},
                ...
            ],
            "go_terms": ["GO:0005524", ...],
            "source": "uniprot|cache",
            "confidence": "high|medium|low",
        }
        """
        # Step 1: Search — fetch up to 5 results, prefer reviewed/canonical with longest sequence
        search_result = await self._request_with_retry(
            api_name="uniprot",
            endpoint="/uniprotkb/search",
            params={
                "query": f"gene:{gene_symbol} AND organism_id:9606",
                "format": "json",
                "size": 5,
            },
        )
        
        if not (search_result["data"] and search_result["http_status"] == 200):
            return {
                "uniprot_id": None,
                "protein_name": None,
                "domains": [],
                "go_terms": [],
                "interpro_ids": [],
                "source": "failed",
                "confidence": "low",
                "error": search_result.get("error"),
            }
        
        # Extract uniprot ID from search results — prefer reviewed + longest sequence
        search_data = search_result["data"]
        results = search_data.get("results", [])
        if not results:
            return {
                "uniprot_id": None,
                "protein_name": None,
                "domains": [],
                "go_terms": [],
                "interpro_ids": [],
                "source": "uniprot",
                "confidence": "medium",
                "error": "No UniProt entry found for gene",
            }
        
        # Pick best entry: reviewed (Swiss-Prot) preferred, then longest sequence
        def _entry_score(entry):
            is_reviewed = 1 if entry.get("entryType") == "UniProtKB reviewed (Swiss-Prot)" else 0
            seq_len = entry.get("sequence", {}).get("length", 0) or 0
            return (is_reviewed, seq_len)
        
        sorted_results = sorted(results, key=_entry_score, reverse=True)
        best_entry = sorted_results[0]
        uniprot_id = best_entry.get("primaryAccession")
        
        # Step 2: Fetch full entry
        entry_result = await self._request_with_retry(
            api_name="uniprot",
            endpoint=f"/uniprotkb/{uniprot_id}.json",
        )
        
        if entry_result["data"] and entry_result["http_status"] == 200:
            data = entry_result["data"]
            
            # Extract domains from features
            domains = []
            for feature in data.get("features", []):
                if feature.get("type", "").lower() in ("domain", "region", "repeat", "zn_fing", "dna_bind"):
                    loc = feature.get("location", {})
                    start = loc.get("start", {}).get("value")
                    end = loc.get("end", {}).get("value")
                    if start and end:
                        domains.append({
                            "name": feature.get("description", "unnamed"),
                            "start": start,
                            "end": end,
                            "type": feature.get("type"),
                        })
            
            # Extract GO terms
            go_terms = []
            for ref in data.get("uniProtKBCrossReferences", []):
                if ref.get("database") == "GO":
                    go_id = ref.get("id")
                    if go_id:
                        go_terms.append(go_id)
            
            # Extract InterPro IDs
            interpro_ids = []
            for ref in data.get("uniProtKBCrossReferences", []):
                if ref.get("database") == "InterPro":
                    ip_id = ref.get("id")
                    if ip_id:
                        interpro_ids.append(ip_id)
            
            # Protein name
            protein_desc = data.get("proteinDescription", {})
            rec_name = protein_desc.get("recommendedName", {}).get("fullName", {}).get("value", "")
            
            seq_length = None
            seq_info = data.get("sequence", {})
            if seq_info:
                seq_length = seq_info.get("length")
            
            return {
                "uniprot_id": uniprot_id,
                "protein_name": rec_name,
                "sequence_length": seq_length,
                "domains": domains,
                "go_terms": go_terms,
                "interpro_ids": interpro_ids,
                "source": "cache" if entry_result["from_cache"] else "uniprot",
                "confidence": entry_result["confidence"],
                "raw": data,
            }
        
        return {
            "uniprot_id": uniprot_id,
            "protein_name": None,
            "domains": [],
            "go_terms": [],
            "interpro_ids": [],
            "source": "failed",
            "confidence": "low",
            "error": entry_result.get("error"),
        }
    
    # =====================================================================
    # GTEx Portal API
    # =====================================================================
    
    async def query_gtex_expression(self, gene_id: str, tissue: str) -> Dict[str, Any]:
        """
        Query GTEx v2 API for median gene expression in a specific tissue.
        
        GTEx v2 requires versioned gencodeIds (e.g. ENSG00000110799.13).
        Two-step process:
          1. Resolve gene symbol -> versioned gencodeId (cached)
          2. Query medianGeneExpression endpoint
        
        Returns:
        {
            "gene": "VWF",
            "tissue": "Whole_Blood",
            "median_tpm": 268.7,
            "unit": "TPM",
            "source": "gtex|cache",
            "confidence": "medium",
        }
        """
        # --- Step 1: Resolve gene symbol -> versioned gencodeId ---
        gencode_map = self._load_gencode_cache()
        gencode_id = gencode_map.get(gene_id)
        
        if not gencode_id:
            # Query GTEx get_genes endpoint
            search_result = await self._request_with_retry(
                api_name="gtex",
                endpoint="/reference/gene",
                params={
                    "geneId": gene_id,
                    "page": 0,
                    "itemsPerPage": 10,
                },
            )
            
            if search_result["data"] and search_result["http_status"] == 200:
                data = search_result["data"]
                items = data.get("data", [])
                if items:
                    # Prefer exact symbol match
                    for item in items:
                        if item.get("geneSymbol") == gene_id:
                            gencode_id = item.get("gencodeId")
                            break
                    if not gencode_id:
                        gencode_id = items[0].get("gencodeId")
                
                if gencode_id:
                    gencode_map[gene_id] = gencode_id
                    self._save_gencode_cache(gencode_map)
        
        if not gencode_id:
            return {
                "gene": gene_id,
                "tissue": tissue,
                "median_tpm": None,
                "unit": "TPM",
                "source": "failed",
                "confidence": "low",
                "error": f"Could not resolve gencodeId for {gene_id}",
            }
        
        # --- Step 2: Query medianGeneExpression ---
        # GTEx tissue IDs use underscores: "Whole_Blood", "Heart_Left_Ventricle"
        gtex_tissue = tissue.replace(" - ", "_").replace(" ", "_")
        
        result = await self._request_with_retry(
            api_name="gtex",
            endpoint="/expression/medianGeneExpression",
            params={
                "gencodeId": gencode_id,
                "tissueSiteDetailIds": gtex_tissue,
                "datasetId": "gtex_v8",
            },
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            items = data.get("data", [])
            
            # Find the tissue-specific record
            median_val = None
            for item in items:
                if item.get("tissueSiteDetailId") == gtex_tissue:
                    median_val = item.get("median")
                    break
            
            # Fallback: if tissue not found but data exists, return first (for debugging)
            if median_val is None and items:
                median_val = items[0].get("median")
            
            return {
                "gene": gene_id,
                "tissue": tissue,
                "median_tpm": median_val,
                "unit": "TPM",
                "gencode_id": gencode_id,
                "source": "cache" if result["from_cache"] else "gtex",
                "confidence": result["confidence"],
                "raw": data,
            }
        
        return {
            "gene": gene_id,
            "tissue": tissue,
            "median_tpm": None,
            "unit": "TPM",
            "gencode_id": gencode_id,
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
        }
    
    def _load_gencode_cache(self) -> Dict[str, str]:
        """Load gene symbol -> gencodeId mapping cache."""
        _script_dir = Path(__file__).resolve().parent
        cache_path = _script_dir / ".." / "references" / "offline_data" / "gtex_gencode_map.json"
        if cache_path.exists():
            try:
                with open(cache_path, "r") as f:
                    return json.load(f)
            except Exception:
                pass
        return {}
    
    def _save_gencode_cache(self, cache: Dict[str, str]) -> None:
        """Save gene symbol -> gencodeId mapping cache."""
        _script_dir = Path(__file__).resolve().parent
        cache_path = _script_dir / ".." / "references" / "offline_data" / "gtex_gencode_map.json"
        cache_path.parent.mkdir(parents=True, exist_ok=True)
        with open(cache_path, "w") as f:
            json.dump(cache, f, indent=2, sort_keys=True)
    
    async def query_gtex_expression_multi(self, gene_id: str, tissues: List[str]) -> List[Dict[str, Any]]:
        """
        Query GTEx v2 API for median gene expression across multiple tissues.
        
        Uses asyncio.gather to concurrently query all tissues.
        Returns a list of individual tissue results.
        
        Args:
            gene_id: Gene symbol (e.g., "VWF")
            tissues: List of GTEx tissue names (e.g., ["Bone Marrow", "Whole Blood"])
            
        Returns:
            List of result dicts, one per tissue (same format as query_gtex_expression)
        """
        if not tissues:
            return []
        
        # v0.10.5 FIX: Single API call with comma-separated tissueSiteDetailIds
        # instead of N concurrent calls. GTEx v2 API supports multiple tissues.
        gencode_map = self._load_gencode_cache()
        gencode_id = gencode_map.get(gene_id)

        if not gencode_id:
            search_result = await self._request_with_retry(
                api_name="gtex",
                endpoint="/reference/gene",
                params={"geneId": gene_id, "page": 0, "itemsPerPage": 10},
            )
            if search_result["data"] and search_result["http_status"] == 200:
                data = search_result["data"]
                items = data.get("data", [])
                if items:
                    for item in items:
                        if item.get("geneSymbol") == gene_id:
                            gencode_id = item.get("gencodeId")
                            break
                    if not gencode_id:
                        gencode_id = items[0].get("gencodeId")
                if gencode_id:
                    gencode_map[gene_id] = gencode_id
                    self._save_gencode_cache(gencode_map)

        if not gencode_id:
            return [
                {
                    "gene": gene_id, "tissue": t, "median_tpm": None,
                    "unit": "TPM", "source": "failed", "confidence": "low",
                    "error": f"Could not resolve gencodeId for {gene_id}",
                }
                for t in tissues
            ]

        gtex_tissues = [t.replace(" - ", "_").replace(" ", "_") for t in tissues]
        tissue_param = ",".join(gtex_tissues)

        result = await self._request_with_retry(
            api_name="gtex",
            endpoint="/expression/medianGeneExpression",
            params={
                "gencodeId": gencode_id,
                "tissueSiteDetailIds": tissue_param,
                "datasetId": "gtex_v8",
            },
        )

        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            items = data.get("data", [])
            tissue_median: Dict[str, float] = {}
            for item in items:
                tid = item.get("tissueSiteDetailId")
                if tid in gtex_tissues:
                    tissue_median[tid] = item.get("median")

            processed = []
            for tissue in tissues:
                gtex_t = tissue.replace(" - ", "_").replace(" ", "_")
                median_val = tissue_median.get(gtex_t)
                processed.append({
                    "gene": gene_id, "tissue": tissue, "median_tpm": median_val,
                    "unit": "TPM", "gencode_id": gencode_id,
                    "source": "cache" if result["from_cache"] else "gtex",
                    "confidence": result["confidence"],
                })
            return processed

        return [
            {
                "gene": gene_id, "tissue": t, "median_tpm": None,
                "unit": "TPM", "gencode_id": gencode_id,
                "source": "failed", "confidence": "low",
                "error": result.get("error"),
            }
            for t in tissues
        ]
    
    # =====================================================================
    # gnomAD GraphQL API
    # =====================================================================
    
    async def query_gnomad_variant(self, chrom: str, pos: int, ref: str, alt: str,
                                    dataset: str = "gnomad_r4",
                                    populations: Optional[List[str]] = None) -> Dict[str, Any]:
        """
        Query gnomAD for variant allele frequency and constraint metrics.
        v0.5 P1-1: Population subgroup frequencies (EAS, AMR, AFR, NFE, SAS, etc.)

        Uses GraphQL API. Returns structured frequency data including per-population AFs.

        NOTE (v0.10.13): Callers should filter variants BEFORE calling this method.
        Only query gnomAD for Tier 1/2 candidate variants — skip Tier 3 variants
        (common SNPs, low impact, high AF) to reduce unnecessary API calls.
        The two-phase pipeline (_enrich_variant_frequencies) handles this filtering.

        Args:
            populations: List of population codes to query. Default: ["EAS", "AMR", "AFR", "NFE", "SAS", "ASJ", "FIN"]
        
        Returns:
        {
            "variant_id": "1-12345-A-G",
            "af": 0.00123,           # overall combined AF
            "af_popmax": 0.00234,   # max across populations
            "af_populations": {      # per-population AFs
                "EAS": {"af": 0.0001, "ac": 1, "an": 15278},
                "AMR": {"af": 0.0023, "ac": 3, "an": 13000},
                ...
            },
            "an": 152000,
            "hom_count": 2,
            "gene_constraint": {"lof_z": 2.5, "pLI": 0.99},
            "source": "gnomad|cache",
            "confidence": "medium",
        }
        """
        if populations is None:
            populations = ["EAS", "AMR", "AFR", "NFE", "SAS", "ASJ", "FIN", "MID", "OTH"]
        
        # Build populations query fragment
        # v0.9.3: gnomAD removed population-specific fields (EAS, AMR, etc.) from VariantPopulation type.
        # Now populations is an array of {id, ac, an, homozygote_count} objects.
        pop_query = """                id
                    ac
                    an
                    homozygote_count"""
        
        query = f"""
        query VariantQuery($variantId: String!, $datasetId: DatasetId!) {{
            variant(variantId: $variantId, dataset: $datasetId) {{
                variantId
                exome {{
                    an
                    ac
                    homozygote_count
                    populations {{
{pop_query}
                    }}
                }}
                genome {{
                    an
                    ac
                    homozygote_count
                    populations {{
{pop_query}
                    }}
                }}
            }}
        }}
        """
        # v0.9.2: Strip chr prefix for gnomAD variant ID format
        chrom_std = chrom.replace("chr", "").replace("CHR", "") if chrom.upper().startswith("CHR") else chrom
        variant_id = f"{chrom_std}-{pos}-{ref}-{alt}"
        
        result = await self._request_with_retry(
            api_name="gnomad",
            endpoint="/",
            method="POST",
            json_body={
                "query": query,
                "variables": {
                    "variantId": variant_id,
                    "datasetId": dataset,
                },
            },
            # v0.10.1: Include variantId in cache key so each variant gets its own cache entry.
            # Without this, all gnomAD GraphQL queries share the same cache key (url only)
            # because params=None and json_body is not used for cache key generation.
            params={"variantId": variant_id, "datasetId": dataset},
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            variant_data = data.get("data", {}).get("variant", {})

            # === P0 FIX 2026-05-24: Check for GraphQL errors before declaring NOT_CAPTURED ===
            # Distinguish "variant truly not in gnomAD" from "query failed (e.g., removed field)"
            gql_errors = data.get("errors", [])
            gql_error_msgs = []
            if gql_errors:
                gql_error_msgs = [e.get("message", str(e)) for e in gql_errors]
                # Cache this failed result with very short TTL to avoid persistent false negatives
                self.cache.set(
                    api_name="gnomad",
                    response_data={"error": "GraphQL error", "errors": gql_error_msgs},
                    http_status=200,
                    confidence="low",
                    ttl_days=0.0035,  # ~5 minutes for query errors
                    url=f"{self.config.apis['gnomad'].base_url}/",
                    **{"queryVariant": variant_id, "dataset": dataset},
                )

            if not variant_data:
                if gql_error_msgs:
                    return {
                        "variant_id": variant_id,
                        "af": None,
                        "af_popmax": None,
                        "af_populations": {},
                        "status": "QUERY_ERROR",
                        "source": "gnomad",
                        "confidence": "low",
                        "note": f"gnomAD GraphQL query returned errors: {'; '.join(gql_error_msgs[:3])}. "
                                "Data may exist but query failed — fallback recommended.",
                        "raw": data,
                        "graphql_errors": gql_error_msgs,
                    }
                return {
                    "variant_id": variant_id,
                    "af": None,
                    "af_popmax": None,
                    "af_populations": {},
                    "status": "NOT_CAPTURED",
                    "source": "gnomad",
                    "confidence": "medium",
                    "note": "Variant not in gnomAD dataset (no GraphQL errors — confirmed absent)",
                    "raw": data,
                }
            
            # Combine exome + genome AF (v0.9.3: gnomAD removed exome.af/genome.af fields)
            exome = variant_data.get("exome", {}) or {}
            genome = variant_data.get("genome", {}) or {}
            
            ex_ac = exome.get("ac", 0) or 0
            ex_an = exome.get("an", 0) or 0
            gen_ac = genome.get("ac", 0) or 0
            gen_an = genome.get("an", 0) or 0
            
            exome_af = ex_ac / ex_an if ex_an > 0 else None
            genome_af = gen_ac / gen_an if gen_an > 0 else None
            
            # Use whichever is available, prefer combined
            if exome_af is not None and genome_af is not None:
                combined_af = (exome_af * ex_an + genome_af * gen_an) / (ex_an + gen_an)
            elif exome_af is not None:
                combined_af = exome_af
            elif genome_af is not None:
                combined_af = genome_af
            else:
                combined_af = None
            
            # v0.5 P1-1: Aggregate per-population frequencies across exome + genome
            # v0.9.3: gnomAD populations now returns array of {id, ac, an, homozygote_count}
            af_populations = {}
            ex_pops_raw = exome.get("populations", []) or []
            gen_pops_raw = genome.get("populations", []) or []
            
            # Convert list format to dict format keyed by population id
            ex_pops = {p.get("id"): p for p in ex_pops_raw if p.get("id")}
            gen_pops = {p.get("id"): p for p in gen_pops_raw if p.get("id")}
            
            all_pops = set(ex_pops.keys()) | set(gen_pops.keys())
            popmax_af = 0.0
            
            for pop in all_pops:
                ex_pop = ex_pops.get(pop, {}) or {}
                gen_pop = gen_pops.get(pop, {}) or {}
                
                ex_ac = ex_pop.get("ac", 0) or 0
                ex_an = ex_pop.get("an", 0) or 0
                ex_hom = ex_pop.get("homozygote_count", 0) or 0
                
                gen_ac = gen_pop.get("ac", 0) or 0
                gen_an = gen_pop.get("an", 0) or 0
                gen_hom = gen_pop.get("homozygote_count", 0) or 0
                
                total_ac = ex_ac + gen_ac
                total_an = ex_an + gen_an
                total_hom = ex_hom + gen_hom
                
                pop_af = total_ac / total_an if total_an > 0 else None
                
                if pop_af is not None:
                    af_populations[pop] = {
                        "af": pop_af,
                        "ac": total_ac,
                        "an": total_an,
                        "homozygote_count": total_hom,
                    }
                    if pop_af > popmax_af:
                        popmax_af = pop_af
            
            return {
                "variant_id": variant_id,
                "af": combined_af,
                "af_popmax": popmax_af if af_populations else None,
                "af_populations": af_populations,
                "af_exome": exome_af,
                "af_genome": genome_af,
                "an_exome": exome.get("an"),
                "an_genome": genome.get("an"),
                "hom_count": (exome.get("homozygote_count") or 0) + (genome.get("homozygote_count") or 0),
                "status": "SUCCESS",
                "source": "cache" if result["from_cache"] else "gnomad",
                "confidence": result["confidence"],
                "raw": data,
            }
        
        return {
            "variant_id": variant_id,
            "af": None,
            "af_popmax": None,
            "af_populations": {},
            "status": "API_FAILED",
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
        }
    
    # =====================================================================
    # gnomAD Gene Constraint (v0.5 P1-4: pLI/LOEUF)
    # =====================================================================
    
    async def query_gnomad_gene_constraint(self, gene_symbol: str) -> Dict[str, Any]:
        """
        Query gnomAD gene constraint metrics (pLI, LOEUF, lof_z, mis_z).
        
        Uses GraphQL API. Returns structured constraint data.
        
        Returns:
        {
            "gene": "BRCA1",
            "pLI": 1.0,
            "lof_z": 5.2,
            "mis_z": 3.1,
            "loeuf": 0.12,       # Loss-of-function observed/expected upper bound fraction
            "oe_lof": 0.08,       # Observed/expected LOF
            "source": "gnomad|cache",
            "confidence": "medium",
        }
        """
        query = """
        query GeneConstraint($geneSymbol: String!, $datasetId: DatasetId!) {
            gene(gene_symbol: $geneSymbol, reference_genome: GRCh38) {
                gnomad_constraint {
                    pLI
                    oe_lof
                    oe_lof_upper  # LOEUF
                    oe_lof_lower
                    lof_z
                    mis_z
                }
            }
        }
        """
        
        result = await self._request_with_retry(
            api_name="gnomad",
            endpoint="/",
            method="POST",
            json_body={
                "query": query,
                "variables": {
                    "geneSymbol": gene_symbol,
                    "datasetId": "gnomad_r4",
                },
            },
            params={"geneSymbol": gene_symbol, "datasetId": "gnomad_r4"},
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            gene_data = data.get("data", {}).get("gene", {})
            constraint = gene_data.get("gnomad_constraint", {}) if gene_data else {}
            
            if not constraint:
                return {
                    "gene": gene_symbol,
                    "pLI": None,
                    "lof_z": None,
                    "mis_z": None,
                    "loeuf": None,
                    "oe_lof": None,
                    "status": "NO_CONSTRAINT_DATA",
                    "source": "gnomad",
                    "confidence": "medium",
                    "note": f"No constraint data for {gene_symbol}",
                }
            
            return {
                "gene": gene_symbol,
                "pLI": constraint.get("pLI"),
                "lof_z": constraint.get("lof_z"),
                "mis_z": constraint.get("mis_z"),
                "loeuf": constraint.get("oe_lof_upper"),
                "oe_lof": constraint.get("oe_lof"),
                "oe_lof_lower": constraint.get("oe_lof_lower"),
                "status": "SUCCESS",
                "source": "cache" if result["from_cache"] else "gnomad",
                "confidence": result["confidence"],
            }
        
        return {
            "gene": gene_symbol,
            "pLI": None,
            "lof_z": None,
            "mis_z": None,
            "loeuf": None,
            "oe_lof": None,
            "status": "API_FAILED",
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
        }
    
    # =====================================================================
    # NCBI E-utilities (ClinVar, Gene)
    # =====================================================================

    @staticmethod
    def _parse_clinvar_efetch_json(raw: Dict[str, Any]) -> Dict[str, Any]:
        """
        Parse ClinVar efetch JSON (retmode=json) into structured data.

        v0.10.5: Robust parser handling both PascalCase and camelCase keys,
        conflict detection across submitters, and review-status mapping.

        Returns:
            {
                "clinical_significance": "Pathogenic" | "Likely pathogenic" | ...,
                "review_status": "practice_guideline" | "reviewed by expert panel" | ...,
                "conflicting": bool,
                "trait_names": List[str],
                "submitter_count": int,
            }
        """
        # Helper: case-insensitive deep get
        def _get(data: Any, *keys: str) -> Any:
            if not isinstance(data, dict):
                return None
            for k in keys:
                if k in data:
                    return data[k]
                # Try case variants
                for variant in (k.lower(), k.upper(), k.title(), k):
                    if variant in data:
                        return data[variant]
            return None

        # Navigate to ClinVarSet (may be list or dict)
        clinvar_set = _get(raw, "ClinVarSet", "clinvarSet")
        if isinstance(clinvar_set, list) and clinvar_set:
            entry = clinvar_set[0]
        elif isinstance(clinvar_set, dict):
            entry = clinvar_set
        else:
            entry = raw  # Fallback: raw may already be the inner dict

        # --- ReferenceClinVarAssertion (the authoritative assertion) ---
        ref_assertion = _get(entry, "ReferenceClinVarAssertion", "referenceClinVarAssertion")

        significance = None
        review_status = None
        trait_names = []

        if ref_assertion:
            cs = _get(ref_assertion, "ClinicalSignificance", "clinicalSignificance")
            if cs:
                significance = _get(cs, "Description", "description")
                review_status = _get(cs, "ReviewStatus", "reviewStatus")

            # Traits
            trait_set = _get(ref_assertion, "TraitSet", "traitSet")
            if trait_set:
                traits = _get(trait_set, "Trait", "trait")
                if isinstance(traits, dict):
                    traits = [traits]
                if isinstance(traits, list):
                    for t in traits:
                        name_data = _get(t, "Name", "name")
                        if name_data:
                            val = _get(name_data, "ElementValue", "elementValue")
                            if isinstance(val, dict):
                                name = val.get("$", val.get("value", val.get("Value")))
                            elif isinstance(val, str):
                                name = val
                            else:
                                name = None
                            if name:
                                trait_names.append(name)

        # --- ClinVarAssertion (submitter assertions) — check for conflicts ---
        assertions = _get(entry, "ClinVarAssertion", "clinVarAssertion")
        submitter_sigs = []
        if isinstance(assertions, dict):
            assertions = [assertions]
        if isinstance(assertions, list):
            for a in assertions:
                a_cs = _get(a, "ClinicalSignificance", "clinicalSignificance")
                if a_cs:
                    desc = _get(a_cs, "Description", "description")
                    if desc:
                        submitter_sigs.append(desc)

        # Detect conflict: mixed pathogenic/benign directions
        patho_like = {"pathogenic", "likely pathogenic", "pathogenic/likely pathogenic"}
        benign_like = {"benign", "likely benign", "benign/likely benign"}
        has_patho = any(s.lower() in patho_like for s in submitter_sigs)
        has_benign = any(s.lower() in benign_like for s in submitter_sigs)
        conflicting = has_patho and has_benign

        # Normalize review status to lowercase standard form
        if review_status:
            review_status = review_status.lower()

        return {
            "clinical_significance": significance,
            "review_status": review_status,
            "conflicting": conflicting,
            "trait_names": trait_names,
            "submitter_count": len(submitter_sigs),
        }

    @staticmethod
    def _parse_clinvar_esummary_json(raw: Dict[str, Any], cv_id: str) -> Dict[str, Any]:
        """
        Parse ClinVar esummary JSON (retmode=json) into structured data.

        v0.10.6 FIX: NCBI ClinVar efetch does NOT support retmode=json — it
        always returns XML. ESummary is the only E-utilities endpoint that
        returns structured JSON for ClinVar records.

        Returns:
            {
                "clinical_significance": "Pathogenic" | "Likely pathogenic" | ...,
                "review_status": "criteria provided, multiple submitters, no conflicts" | ...,
                "conflicting": bool,
                "trait_names": List[str],
                "submitter_count": int,
            }
        """
        result = raw.get("result", {})
        doc = result.get(str(cv_id), {})

        if not doc:
            return {
                "clinical_significance": None,
                "review_status": None,
                "conflicting": False,
                "trait_names": [],
                "submitter_count": 0,
            }

        # Germline classification (the aggregated clinical significance)
        gc = doc.get("germline_classification", {})
        significance = gc.get("description") if isinstance(gc, dict) else None
        review_status = gc.get("review_status") if isinstance(gc, dict) else None

        # Traits
        trait_names = []
        trait_set = gc.get("trait_set", []) if isinstance(gc, dict) else []
        if isinstance(trait_set, list):
            for t in trait_set:
                name = t.get("trait_name") if isinstance(t, dict) else None
                if name:
                    trait_names.append(name)

        # Submitter count from supporting_submissions
        submissions = doc.get("supporting_submissions", {})
        scv_list = submissions.get("scv", []) if isinstance(submissions, dict) else []
        submitter_count = len(scv_list)

        # Conflict detection: esummary only provides the aggregated classification,
        # so we cannot detect submitter conflicts from esummary alone.
        # For conflict detection, full efetch XML parsing would be needed.
        conflicting = False

        return {
            "clinical_significance": significance,
            "review_status": review_status,
            "conflicting": conflicting,
            "trait_names": trait_names,
            "submitter_count": submitter_count,
        }

    async def query_ncbi_clinvar(self, gene: str, hgvs: Optional[str] = None,
                                  chrom: Optional[str] = None, pos: Optional[int] = None) -> Dict[str, Any]:
        """
        Query ClinVar via NCBI E-utilities for clinical significance.
        
        v0.10.6 FIX: 
        - ESearch: [pos] → [chrpos] (zero-padded for proxy stability)
        - EFetch → ESummary (efetch does not support retmode=json for ClinVar)
        
        Previously gene-level search (`{gene}[Gene] AND ClinVar[Title]`) could only find
        gene records, not specific variant records. Now supports:
        1. Position-based search: `{chrom}[chr] AND {pos}[chrpos]` — finds exact variant
        2. Gene + position: `{gene}[Gene] AND {chrom}[chr] AND {pos}[chrpos]`
        3. Fallback to gene-level: `{gene}[Gene] AND ClinVar[Title]`
        
        Uses esearch to find ClinVar records, then esummary for JSON details.
        
        Args:
            gene: Gene symbol
            hgvs: HGVS notation (e.g., "NM_007294.4:c.68_69del")
            chrom: Chromosome (e.g., "1", "chr1", "X")
            pos: Genomic position (1-based)
        
        Returns:
        {
            "gene": "VWF",
            "clinvar_id": "RCV000012345.6",
            "clinical_significance": "Pathogenic",
            "review_status": "practice_guideline",
            "source": "clinvar|cache",
            "confidence": "medium",
        }
        """
        # Step 1: esearch with variant-level search strategies
        # Strategy priority: position-based > gene+position > gene+HGVS > gene-only
        search_terms = []
        
        if chrom and pos:
            # Strip "chr" prefix for NCBI queries
            chrom_std = chrom.replace("chr", "").replace("CHR", "")
            # v0.10.6 FIX: [pos] is not a valid ClinVar ESearch field — use [chrpos].
            # Also add zero-padded variant because raw position may fail in proxy mode.
            pos_padded = str(pos).zfill(9)
            search_terms.append(f"{chrom_std}[chr] AND {pos}[chrpos]")
            search_terms.append(f"{gene}[Gene] AND {chrom_std}[chr] AND {pos}[chrpos]")
            search_terms.append(f"{chrom_std}[chr] AND {pos_padded}[chrpos]")
            search_terms.append(f"{gene}[Gene] AND {chrom_std}[chr] AND {pos_padded}[chrpos]")
        
        if hgvs:
            search_terms.append(f"{gene}[Gene] AND {hgvs}")
        
        # Gene-level fallback (original behavior)
        search_terms.append(f"{gene}[Gene] AND ClinVar[Title]")
        
        search_data = None
        idlist = []
        used_search = ""
        
        for search_term in search_terms:
            search_result = await self._request_with_retry(
                api_name="clinvar_eutils",
                endpoint="/esearch.fcgi",
                params={
                    "db": "clinvar",
                    "term": search_term,
                    "retmode": "json",
                    "retmax": 5,
                },
            )
            
            if not (search_result["data"] and search_result["http_status"] == 200):
                continue
            
            search_data = search_result["data"]
            idlist = search_data.get("esearchresult", {}).get("idlist", [])
            used_search = search_term
            
            if idlist:
                # Found results with this search strategy — stop here
                break
        
        if not search_data or not (search_result.get("data") and search_result.get("http_status") == 200):
            return {
                "gene": gene,
                "clinvar_id": None,
                "clinical_significance": None,
                "source": "failed",
                "confidence": "low",
                "error": search_result.get("error") if search_result else "All ClinVar search strategies failed",
            }
        
        if not idlist:
            return {
                "gene": gene,
                "clinvar_id": None,
                "clinical_significance": None,
                "source": "clinvar",
                "confidence": "medium",
                "note": f"No ClinVar records found (searched: {used_search})",
            }
        
        # Step 2: esummary (JSON) — efetch only returns XML for ClinVar
        # v0.10.6 FIX: NCBI ClinVar efetch does not support retmode=json.
        # Use esummary.fcgi with retmode=json instead.
        clinvar_id = idlist[0]
        summary_result = await self._request_with_retry(
            api_name="clinvar_eutils",
            endpoint="/esummary.fcgi",
            params={
                "db": "clinvar",
                "id": clinvar_id,
                "retmode": "json",
            },
        )
        
        if summary_result["data"] and summary_result["http_status"] == 200:
            parsed = self._parse_clinvar_esummary_json(summary_result["data"], clinvar_id)
            return {
                "gene": gene,
                "clinvar_id": clinvar_id,
                "clinical_significance": parsed.get("clinical_significance"),
                "review_status": parsed.get("review_status"),
                "conflicting": parsed.get("conflicting", False),
                "trait_names": parsed.get("trait_names", []),
                "submitter_count": parsed.get("submitter_count", 0),
                "source": "cache" if summary_result["from_cache"] else "clinvar",
                "confidence": summary_result["confidence"],
                "raw": summary_result["data"],
            }
        
        return {
            "gene": gene,
            "clinvar_id": clinvar_id,
            "clinical_significance": None,
            "source": "failed",
            "confidence": "low",
            "error": summary_result.get("error"),
        }
    
    # =====================================================================
    # HGNC Symbol Normalization (v0.5 P1-2)
    # =====================================================================
    
    async def query_hgnc_symbol(self, symbol: str) -> Dict[str, Any]:
        """
        Query HGNC REST API to validate and normalize a gene symbol.
        
        Uses /search/symbol:{symbol} endpoint. Returns approval status,
        previous/alias symbols, and HGNC ID.
        
        Returns:
        {
            "input": "BRCA1",
            "approved_symbol": "BRCA1",
            "hgnc_id": "HGNC:1100",
            "status": "approved",       # approved | previous | alias | withdrawn | not_found
            "previous_symbols": [],
            "alias_symbols": [],
            "locus_type": "gene with protein product",
            "source": "hgnc|cache",
            "confidence": "medium",
        }
        """
        result = await self._request_with_retry(
            api_name="hgnc",
            endpoint=f"/search/symbol:{symbol}",
            headers={"Accept": "application/json"},
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            docs = data.get("response", {}).get("docs", [])
            
            if not docs:
                return {
                    "input": symbol,
                    "approved_symbol": symbol,
                    "hgnc_id": None,
                    "status": "not_found",
                    "previous_symbols": [],
                    "alias_symbols": [],
                    "locus_type": None,
                    "source": "hgnc",
                    "confidence": "medium",
                }
            
            doc = docs[0]  # Best match
            hgnc_id = doc.get("hgnc_id")
            approved = doc.get("symbol")
            status = doc.get("status", "unknown")
            
            # HGNC returns status as strings like "Approved", "Entry Withdrawn"
            status_norm = status.lower().replace(" ", "_")
            if "approved" in status_norm and "withdrawn" not in status_norm:
                status_code = "approved"
            elif "withdrawn" in status_norm:
                status_code = "withdrawn"
            elif status_norm == "previous":
                status_code = "previous"
            elif status_norm == "alias":
                status_code = "alias"
            else:
                status_code = "unknown"
            
            # If input doesn't match approved symbol, check if it's a previous/alias
            if symbol.upper() != approved.upper():
                # Check previous symbols
                prev_syms = doc.get("prev_symbol", [])
                if symbol.upper() in [s.upper() for s in prev_syms]:
                    status_code = "previous"
                
                # Check alias symbols
                alias_syms = doc.get("alias_symbol", [])
                if symbol.upper() in [s.upper() for s in alias_syms]:
                    status_code = "alias"
            
            return {
                "input": symbol,
                "approved_symbol": approved,
                "hgnc_id": hgnc_id,
                "status": status_code,
                "previous_symbols": doc.get("prev_symbol", []),
                "alias_symbols": doc.get("alias_symbol", []),
                "locus_type": doc.get("locus_type"),
                "source": "cache" if result["from_cache"] else "hgnc",
                "confidence": result["confidence"],
            }
        
        return {
            "input": symbol,
            "approved_symbol": symbol,
            "hgnc_id": None,
            "status": "query_failed",
            "previous_symbols": [],
            "alias_symbols": [],
            "locus_type": None,
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
        }
    
    async def batch_normalize_gene_symbols(
        self,
        symbols: List[str],
    ) -> Dict[str, Dict[str, Any]]:
        """
        Batch normalize gene symbols via HGNC API.
        
        Returns dict mapping original symbol -> normalization result.
        Uses batch_query_genes with hgnc query_type internally.
        """
        return await self.batch_query_genes(symbols, query_type="hgnc")
    
    # =====================================================================
    # Batch Query Support
    # =====================================================================
    
    async def batch_query_genes(
        self,
        gene_symbols: List[str],
        query_type: str = "uniprot",
        **kwargs
    ) -> Dict[str, Dict[str, Any]]:
        """
        Execute batch queries with controlled concurrency.
        
        Strategy: semaphore (max 3 concurrent) + chunked batches (30 per batch)
        to avoid overwhelming public APIs while maintaining throughput.
        """
        CHUNK_SIZE = 30
        MAX_CONCURRENT = 3  # v0.9.2: reduced from 20 → 3 to respect rate limits (especially gnomAD)
        semaphore = asyncio.Semaphore(MAX_CONCURRENT)
        
        async def _query_one(gene: str) -> Dict[str, Any]:
            async with semaphore:
                if query_type == "uniprot":
                    return await self.query_uniprot_by_gene(gene)
                elif query_type == "ensembl":
                    return await self.query_ensembl_gene(gene)
                elif query_type == "gtex":
                    tissue = kwargs.get("tissue", "Whole Blood")
                    return await self.query_gtex_expression(gene, tissue)
                elif query_type == "gtex_multi":
                    tissues = kwargs.get("tissues", ["Whole Blood"])
                    return await self.query_gtex_expression_multi(gene, tissues)
                elif query_type == "hgnc":
                    return await self.query_hgnc_symbol(gene)
                elif query_type == "gnomad_constraint":
                    return await self.query_gnomad_gene_constraint(gene)
                else:
                    raise ValueError(f"Unknown query_type: {query_type}")
        
        results = {}
        total = len(gene_symbols)
        
        for i in range(0, total, CHUNK_SIZE):
            chunk = gene_symbols[i:i + CHUNK_SIZE]
            # Create tasks for this chunk
            chunk_tasks = {gene: asyncio.create_task(_query_one(gene)) for gene in chunk}
            # Wait for all in this chunk
            chunk_results = await asyncio.gather(*chunk_tasks.values(), return_exceptions=True)
            # Store results
            for gene, result in zip(chunk_tasks.keys(), chunk_results):
                if isinstance(result, Exception):
                    results[gene] = {
                        "gene": gene,
                        "source": "failed",
                        "confidence": "low",
                        "error": str(result),
                    }
                else:
                    results[gene] = result
            
            # Brief pause between chunks to be polite to APIs
            if i + CHUNK_SIZE < total:
                await asyncio.sleep(0.5)
        
        return results
    
    # =====================================================================
    # v0.9.4 P2: Runtime API Health Check
    # =====================================================================
    
    async def probe_api_health(self) -> Dict[str, Dict[str, Any]]:
        """
        Probe all configured APIs with lightweight queries to verify availability.
        
        Returns dict mapping api_name → {status, latency_ms, error_msg, details}.
        
        Status values:
        - "OK": API responded successfully within timeout
        - "SLOW": API responded but exceeded latency threshold
        - "ERROR": API returned error or timeout
        - "SKIPPED": API not probed (offline mode or no probe query defined)
        """
        health = {}
        
        # Probe queries for each API (lightweight, fast, unlikely to be rate-limited)
        PROBE_QUERIES = {
            "ensembl": {
                "method": "GET",
                "endpoint": "/lookup/symbol/homo_sapiens/BRCA1?content-type=application/json",
                "timeout": 10,
                "latency_threshold_ms": 3000,
            },
            "uniprot": {
                "method": "GET", 
                "endpoint": "/uniprotkb/P04637?format=json&fields=accession,gene_names",
                "timeout": 15,
                "latency_threshold_ms": 5000,
            },
            "gtex": {
                "method": "GET",
                "endpoint": "/expression/medianGeneExpression",
                "params": {
                    "gencodeId": "ENSG00000012048.23",
                    "datasetId": "gtex_v8",
                    "tissueSiteDetailIds": "Whole_Blood",
                    "format": "json",
                },
                "timeout": 15,
                "latency_threshold_ms": 5000,
            },
            "gnomad": {
                "method": "POST",
                "endpoint": "/",
                "json_body": {
                    "query": "query { gene(gene_symbol: \"BRCA1\", reference_genome: GRCh38) { gene_id } }",
                    "variables": {},
                },
                "timeout": 10,
                "latency_threshold_ms": 3000,
            },
            "ncbi_eutils": {
                "method": "GET",
                "endpoint": "/esearch.fcgi?db=gene&term=BRCA1[Gene]&retmax=1&retmode=json",
                "timeout": 10,
                "latency_threshold_ms": 3000,
            },
            "hgnc": {
                "method": "GET",
                "endpoint": "/search/symbol/BRCA1",
                "timeout": 10,
                "latency_threshold_ms": 2000,
            },
        }
        
        for api_name, probe in PROBE_QUERIES.items():
            if api_name not in self.config.apis:
                health[api_name] = {"status": "SKIPPED", "reason": "API not configured"}
                continue
            
            cfg = self.config.apis[api_name]
            url = f"{cfg.base_url}/{probe['endpoint'].lstrip('/')}"
            timeout = probe.get("timeout", 10)
            latency_threshold = probe.get("latency_threshold_ms", 3000)
            
            try:
                await self._rate_limit(api_name)
                start = time.time()
                
                async with self._session.request(
                    method=probe["method"],
                    url=url,
                    params=probe.get("params"),
                    json=probe.get("json_body"),
                    proxy=self._proxy_url,  # v0.10.5: auto-detected proxy
                    timeout=aiohttp.ClientTimeout(total=timeout),
                ) as resp:
                    elapsed_ms = (time.time() - start) * 1000
                    status = resp.status
                    
                    if status == 200:
                        rating = "OK" if elapsed_ms < latency_threshold else "SLOW"
                        health[api_name] = {
                            "status": rating,
                            "latency_ms": round(elapsed_ms, 1),
                            "http_status": status,
                            "details": f"Responded in {elapsed_ms:.0f}ms",
                        }
                    elif status == 429:
                        health[api_name] = {
                            "status": "ERROR",
                            "latency_ms": round(elapsed_ms, 1),
                            "http_status": 429,
                            "error_msg": "Rate limited (429) — API may be overloaded",
                            "details": "Consider reducing rate_limit_per_sec or adding delays",
                        }
                    else:
                        health[api_name] = {
                            "status": "ERROR",
                            "latency_ms": round(elapsed_ms, 1),
                            "http_status": status,
                            "error_msg": f"HTTP {status}",
                            "details": "Unexpected response status",
                        }
                        
            except asyncio.TimeoutError:
                health[api_name] = {
                    "status": "ERROR",
                    "latency_ms": timeout * 1000,
                    "http_status": None,
                    "error_msg": f"Timeout after {timeout}s",
                    "details": "API may be unreachable or too slow",
                }
            except aiohttp.ClientError as e:
                health[api_name] = {
                    "status": "ERROR",
                    "latency_ms": 0,
                    "http_status": None,
                    "error_msg": f"Connection error: {str(e)[:100]}",
                    "details": "Network or DNS issue — check connectivity and proxy settings",
                }
            except Exception as e:
                health[api_name] = {
                    "status": "ERROR",
                    "latency_ms": 0,
                    "http_status": None,
                    "error_msg": f"Unexpected: {type(e).__name__}: {str(e)[:100]}",
                    "details": "Unknown error during health check",
                }
        
        return health
    
    @staticmethod
    def format_health_report(health: Dict[str, Dict[str, Any]]) -> str:
        """Format health check results into a human-readable report."""
        lines = ["=== GPA API Health Check ===\n"]
        
        ok_count = 0
        slow_count = 0
        error_count = 0
        
        for api_name, result in sorted(health.items()):
            status = result.get("status", "UNKNOWN")
            latency = result.get("latency_ms", 0)
            
            if status == "OK":
                icon = "✅"
                ok_count += 1
            elif status == "SLOW":
                icon = "⚠️"
                slow_count += 1
            elif status == "ERROR":
                icon = "❌"
                error_count += 1
            else:
                icon = "⏭️"
            
            detail = result.get("error_msg") or result.get("details", "")
            lines.append(f"  {icon} {api_name:12s} {status:6s} ({latency:.0f}ms)  {detail}")
        
        lines.append(f"\n  Summary: {ok_count} OK, {slow_count} SLOW, {error_count} ERROR\n")
        return "\n".join(lines)


# =============================================================================
# Standalone convenience functions for non-async contexts
# =============================================================================

def run_async(coro):
    """Run an async coroutine from sync code."""
    return asyncio.run(coro)


async def demo():
    """Demo: query a few genes across multiple APIs."""
    from dgra_config import DGRAGlobalConfig
    
    config = DGRAGlobalConfig.from_env()
    cache = DGRACache(config.cache_db_path, default_ttl_days=config.cache_ttl_days)
    
    async with DGRAAPIClient(config, cache) as client:
        # Test Ensembl
        print("=== Ensembl: MYH11 ===")
        result = await client.query_ensembl_gene("MYH11")
        print(f"  Canonical: {result.get('canonical_transcript')}")
        print(f"  Biotype: {result.get('biotype')}")
        print(f"  Source: {result.get('source')}")
        
        # Test UniProt
        print("\n=== UniProt: MYH11 ===")
        result = await client.query_uniprot_by_gene("MYH11")
        print(f"  ID: {result.get('uniprot_id')}")
        print(f"  Length: {result.get('sequence_length')}")
        print(f"  Domains: {len(result.get('domains', []))}")
        print(f"  GO terms: {len(result.get('go_terms', []))}")
        print(f"  Source: {result.get('source')}")
        
        # Test cache stats
        print("\n=== Cache Stats ===")
        stats = cache.get_stats()
        for api, s in stats.items():
            print(f"  {api}: hits={s['hits']}, misses={s['misses']}, rate={s['hit_rate']:.1%}")


if __name__ == "__main__":
    run_async(demo())
