#!/usr/bin/env python3
"""
NCBI ClinVar Direct Query Module for GPA
v0.10.5 - 2026-06-06

Replaces MyVariant.info's unreliable ClinVar aggregation with direct
NCBI E-utilities queries. Queries ClinVar by gene + position and parses
clinical significance from XML responses.

新增：本地 ClinVar VCF 查询作为优先路径（离线可用，零延迟）。
本地 VCF 路径: /Users/zhaorongli/WorkBuddy/2026-05-24-17-27-51/tools/clinvar.vcf.gz

Rate limit: NCBI recommends ≤3 requests/second. We use 1 request/second
to be conservative and avoid IP blocking.

API Docs:
- ESearch: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.ESearch
- EFetch: https://www.ncbi.nlm.nih.gov/books/NBK25499/#chapter4.EFetch
"""

import asyncio
import aiohttp
import json
import xml.etree.ElementTree as ET
from typing import List, Dict, Optional, Tuple, Any
from dataclasses import dataclass
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Local ClinVar VCF integration
# ---------------------------------------------------------------------------
_LOCAL_CLINVAR_MODULE = Path("/Users/zhaorongli/.workbuddy/scripts/clinvar_vcf_local.py")


def _try_local_clinvar_query(chrom: str, pos: int, ref: str, alt: str) -> Optional[Dict[str, Any]]:
    """Try querying the local ClinVar VCF before falling back to NCBI API.
    
    Returns parsed dict if found, None if not found or VCF unavailable.
    """
    if not _LOCAL_CLINVAR_MODULE.exists():
        return None
    
    try:
        # Import via exec to avoid hard dependency
        import importlib.util
        spec = importlib.util.spec_from_file_location("clinvar_vcf_local", _LOCAL_CLINVAR_MODULE)
        if spec is None or spec.loader is None:
            return None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        
        result = mod.query_variant(chrom, pos, ref, alt)
        return result
    except Exception:
        # Any error -> fallback to API
        return None


@dataclass
class ClinVarResult:
    """Structured ClinVar query result."""
    variant_id: str                    # chrom:pos_ref>alt
    gene: str
    
    # ClinVar data
    clinvar_accession: Optional[str] = None      # VCV or RCV accession
    clinvar_significance: Optional[str] = None   # Pathogenic / Likely pathogenic / etc.
    clinvar_review_status: Optional[str] = None  # criteria provided, single submitter / etc.
    clinvar_variant_name: Optional[str] = None   # e.g. "NM_001885.2(CRYAB):c.149_163del"
    clinvar_variant_type: Optional[str] = None   # Deletion / SNV / etc.
    
    # Position match info
    position_match: bool = False       # Does the ClinVar record's position match our query?
    match_quality: str = "unknown"     # exact | nearby | no_match
    
    status: str = "unknown"            # success | not_found | error | skipped
    error: Optional[str] = None


def _local_result_to_clinvar_result(
    local: Dict[str, Any],
    variant_id: str,
    gene: str,
) -> "ClinVarResult":
    """Convert local VCF query result to ClinVarResult dataclass."""
    return ClinVarResult(
        variant_id=variant_id,
        gene=gene,
        clinvar_accession=local.get("ID"),
        clinvar_significance=local.get("CLNSIG"),
        clinvar_review_status=local.get("CLNREVSTAT"),
        clinvar_variant_name=local.get("CLNHGVS"),
        clinvar_variant_type=local.get("MC"),
        position_match=True,
        match_quality="exact",
        status="success",
    )


# ---------------------------------------------------------------------------
# Consequence-aware filtering
# ---------------------------------------------------------------------------

_CLINVAR_RELEVANT_CONSEQUENCES = {
    "missense_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "inframe_deletion",
    "inframe_insertion",
    "frameshift_variant",           # Some frameshifts ARE in ClinVar
    "stop_gained",                   # Some nonsense variants ARE in ClinVar
    "splice_donor_variant",
    "splice_acceptor_variant",
}

_CLINVAR_IRRELEVANT_CONSEQUENCES = {
    "synonymous_variant",
    "intron_variant",
    "upstream_gene_variant",
    "downstream_gene_variant",
    "3_prime_UTR_variant",
    "5_prime_UTR_variant",
    "non_coding_transcript_exon_variant",
    "intergenic_variant",
    "regulatory_region_variant",
}


def variant_needs_clinvar(consequence: str) -> bool:
    """Determine if a variant's consequence warrants a ClinVar query.
    
    ClinVar is most informative for missense and splice variants where
    computational prediction alone is insufficient. For PVS1-level
    consequences (stop_gained, frameshift), ClinVar provides supporting
    evidence but is not critical for tiering.
    """
    if not consequence or consequence == "_UNKNOWN_":
        return False
    
    # Split compound consequences
    terms = {t.strip() for t in consequence.lower().split(",")}
    
    # If any relevant term is present, query ClinVar
    if terms & {t.lower() for t in _CLINVAR_RELEVANT_CONSEQUENCES}:
        return True
    
    # If only irrelevant terms, skip
    if terms & {t.lower() for t in _CLINVAR_IRRELEVANT_CONSEQUENCES}:
        return False
    
    # Default: query (unknown consequence)
    return True


# ---------------------------------------------------------------------------
# NCBI E-utilities query functions
# ---------------------------------------------------------------------------

async def _ncbi_esearch_clinvar(
    gene: str,
    chrom: str,
    pos: int,
    session: aiohttp.ClientSession,
    timeout: int = 30,
) -> Tuple[List[str], Optional[str]]:
    """Search ClinVar by gene + position using NCBI ESearch.
    
    Returns:
        (list of ClinVar IDs, error message or None)
    """
    chrom_std = chrom.replace("chr", "").replace("CHR", "")
    # NCBI CHRPOS field uses zero-padded 9-digit format
    pos_padded = str(pos).zfill(9)
    
    # Build query: gene AND chrom AND position
    # Use GRCh38 assembly qualifier when possible
    query = f'{gene}[Gene] AND {chrom_std}[CHR] AND {pos_padded}[CHRPOS]'
    
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
    params = {
        "db": "clinvar",
        "term": query,
        "retmode": "json",
        "retmax": 10,  # Max 10 results per position
    }
    
    try:
        async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=timeout)) as resp:
            if resp.status != 200:
                text = await resp.text()
                return [], f"HTTP {resp.status}: {text[:200]}"
            
            data = await resp.json()
            idlist = data.get("esearchresult", {}).get("idlist", [])
            return idlist, None
            
    except asyncio.TimeoutError:
        return [], "ESearch timeout"
    except Exception as e:
        return [], f"ESearch error: {type(e).__name__}: {e}"


async def _ncbi_efetch_clinvar(
    clinvar_id: str,
    session: aiohttp.ClientSession,
    timeout: int = 30,
) -> Optional[Dict]:
    """Fetch ClinVar record details via EFetch.
    
    Returns parsed dict with clinical significance, review status, etc.
    Returns None if record not found or error.
    """
    url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    params = {
        "db": "clinvar",
        "id": clinvar_id,
        "rettype": "vcv",
        "is_variationid": "true",
    }
    
    try:
        async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=timeout)) as resp:
            if resp.status != 200:
                return None
            
            xml_text = await resp.text()
            return _parse_clinvar_xml(xml_text)
            
    except Exception:
        return None


def _parse_clinvar_xml(xml_text: str) -> Optional[Dict]:
    """Parse ClinVar XML (VCV format) to extract clinical significance.
    
    Returns dict with:
        - accession: VCV accession
        - variant_name: VariationName attribute
        - variant_type: VariationType attribute
        - clinical_significance: Description text
        - review_status: ReviewStatus text
        - positions: list of {chr, start, stop, ref, alt} from SequenceLocation
    """
    try:
        root = ET.fromstring(xml_text)
    except ET.ParseError:
        return None
    
    result = {
        "accession": None,
        "variant_name": None,
        "variant_type": None,
        "clinical_significance": None,
        "review_status": None,
        "positions": [],
    }
    
    # Find VariationArchive element
    var_arch = root.find(".//VariationArchive")
    if var_arch is not None:
        result["accession"] = var_arch.get("Accession")
        result["variant_name"] = var_arch.get("VariationName")
        result["variant_type"] = var_arch.get("VariationType")
    
    # Find ClinicalSignificance
    # VCV format: <ClassifiedRecord><ClinicalAssertionList><ClinicalAssertion>...<Interpretation><Description>...
    # Also: <ClassifiedRecord><ReviewStatus>...
    for review in root.iter("ReviewStatus"):
        text = review.text
        if text:
            result["review_status"] = text
            break
    
    for desc in root.iter("Description"):
        text = desc.text
        if text and text.lower() not in ("not provided", ""):
            # Prioritize descriptions with DateLastEvaluated (more reliable)
            if desc.get("DateLastEvaluated"):
                result["clinical_significance"] = text
                break
            elif result["clinical_significance"] is None:
                result["clinical_significance"] = text
    
    # Extract positions from SequenceLocation
    for seq_loc in root.iter("SequenceLocation"):
        if seq_loc.get("Assembly") == "GRCh38":
            result["positions"].append({
                "chr": seq_loc.get("Chr"),
                "start": seq_loc.get("start"),
                "stop": seq_loc.get("stop"),
                "ref": seq_loc.get("referenceAlleleVCF"),
                "alt": seq_loc.get("alternateAlleleVCF"),
                "position_vcf": seq_loc.get("positionVCF"),
            })
    
    return result


def _check_position_match(
    query_chrom: str,
    query_pos: int,
    clinvar_positions: List[Dict],
) -> Tuple[bool, str]:
    """Check if ClinVar record's position matches our query variant.
    
    Returns: (is_match, match_quality)
        match_quality: "exact" | "nearby" | "overlap" | "no_match"
    """
    chrom_std = query_chrom.replace("chr", "").replace("CHR", "")
    
    for pos in clinvar_positions:
        cv_chr = pos.get("chr", "")
        cv_start = pos.get("start")
        cv_stop = pos.get("stop")
        cv_pos_vcf = pos.get("position_vcf")
        
        # Check chromosome match
        if str(cv_chr) != str(chrom_std):
            continue
        
        # Exact position match (VCF position)
        if cv_pos_vcf and int(cv_pos_vcf) == query_pos:
            return True, "exact"
        
        # Position within the ClinVar variant's span
        if cv_start and cv_stop:
            start = int(cv_start)
            stop = int(cv_stop)
            if start <= query_pos <= stop:
                return True, "overlap"
        
        # Nearby (within 5bp) - could be same variant with different representation
        if cv_pos_vcf:
            diff = abs(int(cv_pos_vcf) - query_pos)
            if diff <= 5:
                return True, "nearby"
    
    return False, "no_match"


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

async def query_clinvar_variant(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    gene: str,
    consequence: str,
    session: aiohttp.ClientSession,
) -> ClinVarResult:
    """Query ClinVar for a single variant using local VCF first, then NCBI E-utilities.
    
    Flow:
        1. Try local ClinVar VCF query (offline, zero latency)
        2. If not found locally -> ESearch: find ClinVar IDs by gene + position
        3. EFetch: get details for each ID
        4. Parse XML: extract clinical significance
        5. Check position match
        6. Return best-matching record
    
    Args:
        chrom: Chromosome (with or without 'chr' prefix)
        pos: Position (1-based)
        ref: Reference allele
        alt: Alternate allele
        gene: Gene symbol
        consequence: VEP consequence (for filtering)
        session: aiohttp ClientSession
    
    Returns:
        ClinVarResult object
    """
    variant_id = f"{chrom}:{pos}_{ref}>{alt}"
    
    # Step 0: Check if consequence warrants a ClinVar query
    if not variant_needs_clinvar(consequence):
        return ClinVarResult(
            variant_id=variant_id,
            gene=gene,
            status="skipped",
            error=f"Consequence '{consequence}' not relevant for ClinVar lookup",
        )
    
    # Step 0.5: Try local ClinVar VCF query first (offline, fast)
    local_result = _try_local_clinvar_query(chrom, pos, ref, alt)
    if local_result:
        return _local_result_to_clinvar_result(local_result, variant_id, gene)
    
    # Step 1: ESearch (NCBI API fallback)
    clinvar_ids, error = await _ncbi_esearch_clinvar(gene, chrom, pos, session)
    if error:
        return ClinVarResult(
            variant_id=variant_id,
            gene=gene,
            status="error",
            error=f"ESearch failed: {error}",
        )
    
    if not clinvar_ids:
        return ClinVarResult(
            variant_id=variant_id,
            gene=gene,
            status="not_found",
            error="No ClinVar records found for this gene+position",
        )
    
    # Step 2: EFetch each ID and find best match
    best_match = None
    best_quality = "no_match"
    
    for cv_id in clinvar_ids:
        record = await _ncbi_efetch_clinvar(cv_id, session)
        if not record:
            continue
        
        is_match, quality = _check_position_match(chrom, pos, record.get("positions", []))
        
        if is_match and quality in ("exact", "overlap"):
            # Best possible match - use immediately
            return ClinVarResult(
                variant_id=variant_id,
                gene=gene,
                clinvar_accession=record.get("accession"),
                clinvar_significance=record.get("clinical_significance"),
                clinvar_review_status=record.get("review_status"),
                clinvar_variant_name=record.get("variant_name"),
                clinvar_variant_type=record.get("variant_type"),
                position_match=True,
                match_quality=quality,
                status="success",
            )
        elif is_match and quality == "nearby" and best_quality != "exact":
            # Nearby match - keep as fallback
            if best_quality in ("no_match",):
                best_match = record
                best_quality = quality
    
    if best_match:
        return ClinVarResult(
            variant_id=variant_id,
            gene=gene,
            clinvar_accession=best_match.get("accession"),
            clinvar_significance=best_match.get("clinical_significance"),
            clinvar_review_status=best_match.get("review_status"),
            clinvar_variant_name=best_match.get("variant_name"),
            clinvar_variant_type=best_match.get("variant_type"),
            position_match=True,
            match_quality=best_quality,
            status="success",
        )
    
    # No position match found - ClinVar has records for this gene but not this exact variant
    return ClinVarResult(
        variant_id=variant_id,
        gene=gene,
        status="not_found",
        error=f"Found {len(clinvar_ids)} ClinVar record(s) for {gene} at this position, but none match the exact variant",
    )


async def query_clinvar_batch(
    variants: List[Tuple[str, int, str, str, str, str]],  # (chrom, pos, ref, alt, gene, consequence)
    session: aiohttp.ClientSession,
    semaphore: Optional[asyncio.Semaphore] = None,
    rate_limit_delay: float = 1.0,  # NCBI: 1 request per second
) -> Dict[str, ClinVarResult]:
    """Query ClinVar for multiple variants with rate limiting.
    
    Args:
        variants: List of (chrom, pos, ref, alt, gene, consequence) tuples
        session: aiohttp ClientSession
        semaphore: Optional semaphore for concurrency control
        rate_limit_delay: Delay between requests in seconds (default 1.0)
    
    Returns:
        Dict mapping variant_id -> ClinVarResult
    """
    sem = semaphore or asyncio.Semaphore(1)  # Default: 1 concurrent request
    results: Dict[str, ClinVarResult] = {}
    
    async def _query_one(var_data):
        chrom, pos, ref, alt, gene, consequence = var_data
        async with sem:
            result = await query_clinvar_variant(chrom, pos, ref, alt, gene, consequence, session)
            await asyncio.sleep(rate_limit_delay)
            return result
    
    tasks = [_query_one(v) for v in variants]
    batch_results = await asyncio.gather(*tasks, return_exceptions=True)
    
    for var_data, result in zip(variants, batch_results):
        chrom, pos, ref, alt, gene, consequence = var_data
        variant_id = f"{chrom}:{pos}_{ref}>{alt}"
        
        if isinstance(result, Exception):
            results[variant_id] = ClinVarResult(
                variant_id=variant_id,
                gene=gene,
                status="error",
                error=f"Exception: {type(result).__name__}: {result}",
            )
        else:
            results[variant_id] = result
    
    return results


def apply_clinvar_results(
    variants: List[Any],  # List[Variant] - avoid circular import
    clinvar_results: Dict[str, ClinVarResult],
) -> Dict[str, int]:
    """Apply ClinVar query results to GPA Variant objects.
    
    Only fills in missing or UNKNOWN clinvar fields.
    
    Args:
        variants: List of GPA Variant objects
        clinvar_results: Dict from query_clinvar_batch
    
    Returns:
        Stats dict: {filled, not_found, skipped, errors}
    """
    stats = {"filled": 0, "not_found": 0, "skipped": 0, "errors": 0, "total": 0}
    
    _UNKNOWN_MARKER = "_UNKNOWN_"
    
    for v in variants:
        variant_id = f"{v.chrom}:{v.pos}_{v.ref}>{v.alt}"
        result = clinvar_results.get(variant_id)
        if not result:
            continue
        
        stats["total"] += 1
        
        if result.status == "skipped":
            stats["skipped"] += 1
            continue
        if result.status == "not_found":
            stats["not_found"] += 1
            continue
        if result.status == "error":
            stats["errors"] += 1
            continue
        
        # Fill ClinVar if missing or UNKNOWN
        clinvar_current = getattr(v, 'clinvar', '') or ''
        is_clinvar_missing = not clinvar_current or clinvar_current == _UNKNOWN_MARKER or clinvar_current == "UNKNOWN"
        
        if is_clinvar_missing and result.clinvar_significance:
            v.clinvar = result.clinvar_significance
            if result.clinvar_review_status:
                v.clinvar_review_status = result.clinvar_review_status
            stats["filled"] += 1
    
    return stats


# ---------------------------------------------------------------------------
# Self-test
# ---------------------------------------------------------------------------

async def _test_clinvar():
    """Self-test with known variants."""
    async with aiohttp.ClientSession(trust_env=False) as session:
        # Test 1: Known pathogenic variant (BRCA1)
        print("=== Test 1: BRCA1 known pathogenic ===")
        r1 = await query_clinvar_variant(
            "chr17", 43094692, "G", "A", "BRCA1", "missense_variant", session
        )
        print(f"  status={r1.status}, sig={r1.clinvar_significance}, review={r1.clinvar_review_status}")
        
        # Test 2: CRYAB position (should find c.149_163del nearby)
        print("\n=== Test 2: CRYAB 11:111911576 ===")
        r2 = await query_clinvar_variant(
            "chr11", 111911576, "CGA", "C", "CRYAB", "frameshift_variant", session
        )
        print(f"  status={r2.status}, sig={r2.clinvar_significance}, match={r2.match_quality}")
        
        # Test 3: Intron variant (should be skipped)
        print("\n=== Test 3: Intron variant (should skip) ===")
        r3 = await query_clinvar_variant(
            "chr1", 1000000, "A", "G", "GENE1", "intron_variant", session
        )
        print(f"  status={r3.status}, error={r3.error}")


if __name__ == "__main__":
    asyncio.run(_test_clinvar())
