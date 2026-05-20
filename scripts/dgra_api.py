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
    """
    
    def __init__(self, config: DGRAGlobalConfig, cache: DGRACache):
        self.config = config
        self.cache = cache
        self._session: Optional[aiohttp.ClientSession] = None
        self._last_request_time: Dict[str, float] = {}  # api_name -> timestamp
    
    async def __aenter__(self):
        self._session = aiohttp.ClientSession(
            connector=aiohttp.TCPConnector(limit=50, limit_per_host=20),
            timeout=aiohttp.ClientTimeout(total=120),
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
        cache_key_params = {"url": url, **(params or {})}
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
                
                async with self._session.request(
                    method=method,
                    url=url,
                    params=params,
                    json=json_body,
                    headers=headers,
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
                    
                    elif http_status >= 500:
                        # Server error - retry
                        last_error = DGRAAPIError(api_name, f"Server error {http_status}", http_status)
                        await asyncio.sleep(cfg.retry_delay * (2 ** attempt))
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
                await asyncio.sleep(cfg.retry_delay * (2 ** attempt))
            
            except aiohttp.ClientError as e:
                last_error = DGRAAPIError(api_name, f"Connection error: {e}")
                await asyncio.sleep(cfg.retry_delay * (2 ** attempt))
        
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
    
    # =====================================================================
    # gnomAD GraphQL API
    # =====================================================================
    
    async def query_gnomad_variant(self, chrom: str, pos: int, ref: str, alt: str,
                                    dataset: str = "gnomad_r4") -> Dict[str, Any]:
        """
        Query gnomAD for variant allele frequency and constraint metrics.
        
        Uses GraphQL API. Returns structured frequency data.
        
        Returns:
        {
            "variant_id": "1-12345-A-G",
            "af": 0.00123,
            "af_popmax": 0.00234,
            "an": 152000,
            "hom_count": 2,
            "gene_constraint": {"lof_z": 2.5, "pLI": 0.99},
            "source": "gnomad|cache",
            "confidence": "medium",
        }
        """
        query = """
        query VariantQuery($variantId: String!, $datasetId: DatasetId!) {
            variant(variantId: $variantId, dataset: $datasetId) {
                variantId
                exome {
                    af
                    an
                    ac
                    homozygote_count
                }
                genome {
                    af
                    an
                    ac
                    homozygote_count
                }
            }
        }
        """
        variant_id = f"{chrom}-{pos}-{ref}-{alt}"
        
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
        )
        
        if result["data"] and result["http_status"] == 200:
            data = result["data"]
            variant_data = data.get("data", {}).get("variant", {})
            
            if not variant_data:
                return {
                    "variant_id": variant_id,
                    "af": None,
                    "status": "NOT_CAPTURED",
                    "source": "gnomad",
                    "confidence": "medium",
                    "note": "Variant not in gnomAD dataset",
                    "raw": data,
                }
            
            # Combine exome + genome AF
            exome = variant_data.get("exome", {}) or {}
            genome = variant_data.get("genome", {}) or {}
            
            exome_af = exome.get("af")
            genome_af = genome.get("af")
            
            # Use whichever is available, prefer combined
            if exome_af is not None and genome_af is not None:
                combined_af = (exome_af * exome.get("an", 0) + genome_af * genome.get("an", 0)) / \
                              (exome.get("an", 0) + genome.get("an", 1))
            elif exome_af is not None:
                combined_af = exome_af
            elif genome_af is not None:
                combined_af = genome_af
            else:
                combined_af = None
            
            return {
                "variant_id": variant_id,
                "af": combined_af,
                "af_exome": exome_af,
                "af_genome": genome_af,
                "an_exome": exome.get("an"),
                "an_genome": genome.get("an"),
                "hom_count": (exome.get("homozygote_count") or 0) + (genome.get("homozygote_count") or 0),
                "status": "CAPTURED",
                "source": "cache" if result["from_cache"] else "gnomad",
                "confidence": result["confidence"],
                "raw": data,
            }
        
        return {
            "variant_id": variant_id,
            "af": None,
            "status": "QUERY_FAILED",
            "source": "failed",
            "confidence": "low",
            "error": result.get("error"),
        }
    
    # =====================================================================
    # NCBI E-utilities (ClinVar, Gene)
    # =====================================================================
    
    async def query_ncbi_clinvar(self, gene: str, hgvs: Optional[str] = None) -> Dict[str, Any]:
        """
        Query ClinVar via NCBI E-utilities for clinical significance.
        
        Uses esearch to find ClinVar records, then efetch for details.
        
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
        # Step 1: esearch
        search_term = f"{gene}[Gene] AND ClinVar[Title]"
        if hgvs:
            search_term += f" AND {hgvs}"
        
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
            return {
                "gene": gene,
                "clinvar_id": None,
                "clinical_significance": None,
                "source": "failed",
                "confidence": "low",
                "error": search_result.get("error"),
            }
        
        search_data = search_result["data"]
        idlist = search_data.get("esearchresult", {}).get("idlist", [])
        
        if not idlist:
            return {
                "gene": gene,
                "clinvar_id": None,
                "clinical_significance": None,
                "source": "clinvar",
                "confidence": "medium",
                "note": "No ClinVar records found",
            }
        
        # Step 2: efetch first record
        clinvar_id = idlist[0]
        fetch_result = await self._request_with_retry(
            api_name="clinvar_eutils",
            endpoint="/efetch.fcgi",
            params={
                "db": "clinvar",
                "id": clinvar_id,
                "retmode": "json",
            },
        )
        
        if fetch_result["data"] and fetch_result["http_status"] == 200:
            # Parse ClinVar JSON (structure is complex, Phase 2 will do full parsing)
            # Skeleton: just return the raw data for now
            return {
                "gene": gene,
                "clinvar_id": clinvar_id,
                "clinical_significance": "pending_parsing",  # Phase 2
                "review_status": "pending_parsing",
                "source": "cache" if fetch_result["from_cache"] else "clinvar",
                "confidence": fetch_result["confidence"],
                "raw": fetch_result["data"],
            }
        
        return {
            "gene": gene,
            "clinvar_id": clinvar_id,
            "clinical_significance": None,
            "source": "failed",
            "confidence": "low",
            "error": fetch_result.get("error"),
        }
    
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
        
        Strategy: semaphore (max 20 concurrent) + chunked batches (30 per batch)
        to avoid overwhelming public APIs while maintaining throughput.
        """
        CHUNK_SIZE = 30
        MAX_CONCURRENT = 20
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
