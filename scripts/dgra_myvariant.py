#!/usr/bin/env python3
"""
MyVariant.info Integration Module for GPA
v0.9.2 - 2026-05-23

Aggregates variant annotation from multiple databases via a single API call:
- gnomAD (genome AF, population frequencies)
- ClinVar (clinical significance, review status)
- CADD (phred score)
- dbSNP (rsID)

Benefits over individual gnomAD GraphQL + ClinVar E-utilities queries:
1. Single API call per variant (vs 2+ calls)
2. Better ClinVar coverage (position-based lookup)
3. CADD scores included without separate query
4. Batch endpoint: up to 1000 variants per POST

API Docs: https://docs.myvariant.info/en/latest/
"""

import asyncio
import aiohttp
import json
import time
from typing import List, Dict, Any, Optional, Tuple
from dataclasses import dataclass
from pathlib import Path


@dataclass
class MyVariantResult:
    """Structured result from MyVariant.info query."""
    variant_id: str                    # chr1:g.12345A>G
    queried: bool                      # Was this variant actually queried?
    
    # gnomAD
    gnomad_af: Optional[float] = None
    gnomad_ac: Optional[int] = None
    gnomad_an: Optional[int] = None
    gnomad_populations: Optional[Dict[str, Dict]] = None
    gnomad_source: str = ""
    
    # ClinVar
    clinvar_significance: Optional[str] = None
    clinvar_review_status: Optional[str] = None
    clinvar_rcv: Optional[str] = None
    clinvar_source: str = ""
    
    # CADD
    cadd_phred: Optional[float] = None
    cadd_source: str = ""
    
    # dbSNP
    dbsnp_rs: Optional[str] = None
    
    # Metadata
    status: str = "unknown"           # success | not_found | error | skipped
    error: Optional[str] = None
    raw_response: Optional[Dict] = None


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _build_variant_id(chrom: str, pos: int, ref: str, alt: str) -> str:
    """Build MyVariant.info variant ID: chr1:g.12345A>G
    
    Handles chrom prefix stripping (chr1 → 1, chrX → X).
    """
    chrom_std = chrom.replace("chr", "").replace("CHR", "")
    if chrom_std.upper().startswith("CHR"):
        chrom_std = chrom_std[3:]
    # MyVariant uses chr prefix in the ID
    return f"chr{chrom_std}:g.{pos}{ref}>{alt}"


def _parse_gnomad_from_myvariant(raw: Dict) -> Tuple[Optional[float], Optional[int], Optional[int], Optional[Dict]]:
    """Extract gnomAD AF/AC/AN/populations from MyVariant response."""
    gnomad = raw.get("gnomad_genome") or raw.get("gnomad_exome")
    if not gnomad:
        return None, None, None, None
    
    af = None
    ac = None
    an = None
    pops = {}
    
    # gnomAD genome data
    if isinstance(gnomad, dict):
        af_data = gnomad.get("af", {})
        if isinstance(af_data, dict):
            af = af_data.get("af")
            ac = af_data.get("ac")
            an = af_data.get("an")
        
        # Population frequencies
        # MyVariant stores per-population ACs, we can compute AFs if AN is known
        # Structure: gnomad_genome.ac.ac_eas = count, gnomad_genome.an.an_eas = total
        ac_data = gnomad.get("ac", {})
        an_data = gnomad.get("an", {})
        if isinstance(ac_data, dict) and isinstance(an_data, dict):
            pop_codes = ["afr", "amr", "asj", "eas", "fin", "nfe", "mid", "oth", "sas"]
            for pop in pop_codes:
                pop_ac = ac_data.get(f"ac_{pop}")
                pop_an = an_data.get(f"an_{pop}")
                if pop_ac is not None and pop_an and pop_an > 0:
                    pops[pop.upper()] = {
                        "af": pop_ac / pop_an,
                        "ac": pop_ac,
                        "an": pop_an,
                    }
    
    return af, ac, an, pops if pops else None


def _parse_clinvar_from_myvariant(raw: Dict) -> Tuple[Optional[str], Optional[str], Optional[str]]:
    """Extract ClinVar significance from MyVariant response.
    
    Returns: (significance, review_status, rcv_accession)
    """
    clinvar = raw.get("clinvar")
    if not clinvar or not isinstance(clinvar, dict):
        return None, None, None
    
    # Try rcv list first
    rcv_list = clinvar.get("rcv", [])
    if isinstance(rcv_list, list) and rcv_list:
        # Take first RCV entry
        rcv = rcv_list[0]
        if isinstance(rcv, dict):
            sig = rcv.get("clinical_significance")
            review = rcv.get("review_status")
            acc = rcv.get("accession")
            return sig, review, acc
    
    # Fallback: direct fields
    sig = clinvar.get("clinical_significance")
    review = clinvar.get("review_status")
    acc = clinvar.get("rcv_accession")
    return sig, review, acc


def _parse_cadd_from_myvariant(raw: Dict) -> Optional[float]:
    """Extract CADD phred score from MyVariant response."""
    cadd = raw.get("cadd")
    if not cadd or not isinstance(cadd, dict):
        return None
    return cadd.get("phred")


def _parse_dbsnp_from_myvariant(raw: Dict) -> Optional[str]:
    """Extract dbSNP rsID from MyVariant response."""
    dbsnp = raw.get("dbsnp")
    if not dbsnp or not isinstance(dbsnp, dict):
        return None
    rsid = dbsnp.get("rsid")
    return f"rs{rsid}" if isinstance(rsid, int) else str(rsid) if rsid else None


# ---------------------------------------------------------------------------
# Core query functions
# ---------------------------------------------------------------------------

async def query_myvariant_single(
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    session: aiohttp.ClientSession,
    fields: Optional[List[str]] = None,
) -> MyVariantResult:
    """Query MyVariant.info for a single variant.
    
    Args:
        chrom: Chromosome (with or without 'chr' prefix)
        pos: Position (1-based)
        ref: Reference allele
        alt: Alternate allele
        session: aiohttp ClientSession
        fields: Optional list of fields to request (default: all)
    
    Returns:
        MyVariantResult object
    """
    variant_id = _build_variant_id(chrom, pos, ref, alt)
    url = f"https://myvariant.info/v1/variant/{variant_id}"
    
    params = {}
    if fields:
        params["fields"] = ",".join(fields)
    
    try:
        async with session.get(url, params=params, timeout=aiohttp.ClientTimeout(total=30)) as resp:
            if resp.status == 404:
                return MyVariantResult(
                    variant_id=variant_id,
                    queried=True,
                    status="not_found",
                    error="Variant not in MyVariant.info database",
                )
            if resp.status != 200:
                text = await resp.text()
                return MyVariantResult(
                    variant_id=variant_id,
                    queried=True,
                    status="error",
                    error=f"HTTP {resp.status}: {text[:200]}",
                )
            
            raw = await resp.json()
            
            # Parse gnomAD
            gnomad_af, gnomad_ac, gnomad_an, gnomad_pops = _parse_gnomad_from_myvariant(raw)
            
            # Parse ClinVar
            clinvar_sig, clinvar_review, clinvar_rcv = _parse_clinvar_from_myvariant(raw)
            
            # Parse CADD
            cadd_phred = _parse_cadd_from_myvariant(raw)
            
            # Parse dbSNP
            dbsnp_rs = _parse_dbsnp_from_myvariant(raw)
            
            return MyVariantResult(
                variant_id=variant_id,
                queried=True,
                status="success",
                gnomad_af=gnomad_af,
                gnomad_ac=gnomad_ac,
                gnomad_an=gnomad_an,
                gnomad_populations=gnomad_pops,
                gnomad_source="myvariant_gnomad",
                clinvar_significance=clinvar_sig,
                clinvar_review_status=clinvar_review,
                clinvar_rcv=clinvar_rcv,
                clinvar_source="myvariant_clinvar",
                cadd_phred=cadd_phred,
                cadd_source="myvariant_cadd",
                dbsnp_rs=dbsnp_rs,
                raw_response=raw,
            )
    
    except asyncio.TimeoutError:
        return MyVariantResult(
            variant_id=variant_id,
            queried=True,
            status="error",
            error="Timeout after 30s",
        )
    except aiohttp.ClientError as e:
        return MyVariantResult(
            variant_id=variant_id,
            queried=True,
            status="error",
            error=f"Client error: {e}",
        )
    except Exception as e:
        return MyVariantResult(
            variant_id=variant_id,
            queried=True,
            status="error",
            error=f"Unexpected: {type(e).__name__}: {e}",
        )


async def query_myvariant_batch(
    variants: List[Tuple[str, int, str, str]],
    session: aiohttp.ClientSession,
    semaphore: Optional[asyncio.Semaphore] = None,
    fields: Optional[List[str]] = None,
    batch_size: int = 100,
) -> Dict[str, MyVariantResult]:
    """Query MyVariant.info for multiple variants using batch endpoint.
    
    Args:
        variants: List of (chrom, pos, ref, alt) tuples
        session: aiohttp ClientSession
        semaphore: Optional semaphore for concurrency control
        fields: Optional fields to request
        batch_size: Max variants per batch (MyVariant supports up to 1000)
    
    Returns:
        Dict mapping variant_id → MyVariantResult
    """
    if not variants:
        return {}
    
    results: Dict[str, MyVariantResult] = {}
    variant_ids = [_build_variant_id(c, p, r, a) for c, p, r, a in variants]
    
    # Build lookup from variant_id to (chrom, pos, ref, alt)
    id_to_variant = {vid: var for vid, var in zip(variant_ids, variants)}
    
    sem = semaphore or asyncio.Semaphore(10)
    
    # MyVariant batch endpoint supports up to 1000 variants per POST
    # We chunk into batches
    for i in range(0, len(variant_ids), batch_size):
        batch_ids = variant_ids[i:i + batch_size]
        
        async with sem:
            url = "https://myvariant.info/v1/variant"
            headers = {"Content-Type": "application/json"}
            payload = {
                "ids": batch_ids,
            }
            if fields:
                payload["fields"] = ",".join(fields)
            
            try:
                async with session.post(
                    url,
                    headers=headers,
                    json=payload,
                    timeout=aiohttp.ClientTimeout(total=60),
                ) as resp:
                    if resp.status != 200:
                        text = await resp.text()
                        # Mark all in batch as error
                        for vid in batch_ids:
                            results[vid] = MyVariantResult(
                                variant_id=vid,
                                queried=True,
                                status="error",
                                error=f"Batch HTTP {resp.status}: {text[:200]}",
                            )
                        continue
                    
                    batch_results = await resp.json()
                    if not isinstance(batch_results, list):
                        batch_results = [batch_results]
                    
                    for raw in batch_results:
                        if not isinstance(raw, dict):
                            continue
                        
                        vid = raw.get("_id", raw.get("query", ""))
                        if not vid:
                            continue
                        
                        # Check for not_found
                        if raw.get("notfound", False) or "_id" not in raw:
                            results[vid] = MyVariantResult(
                                variant_id=vid,
                                queried=True,
                                status="not_found",
                                error="Variant not in MyVariant.info database",
                            )
                            continue
                        
                        # Parse successful result
                        gnomad_af, gnomad_ac, gnomad_an, gnomad_pops = _parse_gnomad_from_myvariant(raw)
                        clinvar_sig, clinvar_review, clinvar_rcv = _parse_clinvar_from_myvariant(raw)
                        cadd_phred = _parse_cadd_from_myvariant(raw)
                        dbsnp_rs = _parse_dbsnp_from_myvariant(raw)
                        
                        results[vid] = MyVariantResult(
                            variant_id=vid,
                            queried=True,
                            status="success",
                            gnomad_af=gnomad_af,
                            gnomad_ac=gnomad_ac,
                            gnomad_an=gnomad_an,
                            gnomad_populations=gnomad_pops,
                            gnomad_source="myvariant_gnomad",
                            clinvar_significance=clinvar_sig,
                            clinvar_review_status=clinvar_review,
                            clinvar_rcv=clinvar_rcv,
                            clinvar_source="myvariant_clinvar",
                            cadd_phred=cadd_phred,
                            cadd_source="myvariant_cadd",
                            dbsnp_rs=dbsnp_rs,
                            raw_response=raw,
                        )
            
            except asyncio.TimeoutError:
                for vid in batch_ids:
                    results[vid] = MyVariantResult(
                        variant_id=vid,
                        queried=True,
                        status="error",
                        error="Batch timeout after 60s",
                    )
            except aiohttp.ClientError as e:
                for vid in batch_ids:
                    results[vid] = MyVariantResult(
                        variant_id=vid,
                        queried=True,
                        status="error",
                        error=f"Batch client error: {e}",
                    )
            except Exception as e:
                for vid in batch_ids:
                    results[vid] = MyVariantResult(
                        variant_id=vid,
                        queried=True,
                        status="error",
                        error=f"Batch unexpected: {type(e).__name__}: {e}",
                    )
    
    # Ensure all requested variants have a result (even if not returned by API)
    for vid in variant_ids:
        if vid not in results:
            results[vid] = MyVariantResult(
                variant_id=vid,
                queried=False,
                status="skipped",
                error="Not returned by API",
            )
    
    return results


# ---------------------------------------------------------------------------
# Convenience: apply MyVariant results to GPA Variant objects
# ---------------------------------------------------------------------------

def apply_myvariant_results(
    variants: List[Any],  # List[Variant] - avoid circular import
    myvariant_results: Dict[str, MyVariantResult],
    fill_clinvar: bool = True,
    fill_cadd: bool = True,
) -> Dict[str, int]:
    """Apply MyVariant.info results to GPA Variant objects.
    
    Only fills in missing fields (gnomad_af, clinvar, cadd_phred) —
    does NOT overwrite existing data from VEP annotation.
    
    Args:
        variants: List of GPA Variant objects
        myvariant_results: Dict from query_myvariant_batch
        fill_clinvar: Whether to fill ClinVar data (False for variants where
                     functional impact is already certain, e.g., stop_gained)
        fill_cadd: Whether to fill CADD scores
    
    Returns:
        Stats dict: {gnomad_filled, clinvar_filled, cadd_filled, not_found, errors}
    """
    stats = {
        "gnomad_filled": 0,
        "clinvar_filled": 0,
        "cadd_filled": 0,
        "not_found": 0,
        "errors": 0,
        "total_queried": 0,
    }
    
    # Build lookup: variant key → result
    # Variant key format: "chr1:12345_A>G"
    result_by_key: Dict[str, MyVariantResult] = {}
    for vid, result in myvariant_results.items():
        # vid format: "chr1:g.12345A>G"
        # Convert to our key format
        try:
            parts = vid.replace("chr", "").replace(":g.", ":").split(">")
            if len(parts) == 2:
                left = parts[0]
                alt = parts[1]
                # left = "1:12345A"
                chrom_pos_ref = left
                key = f"chr{chrom_pos_ref}>{alt}"
                result_by_key[key] = result
                # Also store without chr prefix
                result_by_key[chrom_pos_ref + ">" + alt] = result
        except Exception:
            pass
    
    for v in variants:
        # Build key for this variant
        key = f"{v.chrom}:{v.pos}_{v.ref}>{v.alt}"
        key2 = f"{v.chrom}:{v.pos}{v.ref}>{v.alt}"
        key3 = f"{v.chrom.replace('chr', '')}:{v.pos}{v.ref}>{v.alt}"
        
        result = result_by_key.get(key) or result_by_key.get(key2) or result_by_key.get(key3)
        if not result:
            continue
        
        stats["total_queried"] += 1
        
        if result.status == "not_found":
            stats["not_found"] += 1
            continue
        if result.status == "error":
            stats["errors"] += 1
            continue
        
        # Fill gnomAD AF if missing
        if v.gnomad_af is None and result.gnomad_af is not None:
            v.gnomad_af = result.gnomad_af
            if result.gnomad_populations:
                v.gnomad_populations = result.gnomad_populations
            stats["gnomad_filled"] += 1
        
        # Fill ClinVar if missing or UNKNOWN AND fill_clinvar is enabled
        _UNKNOWN_MARKER = "_UNKNOWN_"
        if fill_clinvar:
            clinvar_current = getattr(v, 'clinvar', '') or ''
            is_clinvar_missing = not clinvar_current or clinvar_current == _UNKNOWN_MARKER or clinvar_current == "UNKNOWN"
            if is_clinvar_missing and result.clinvar_significance:
                v.clinvar = result.clinvar_significance
                if result.clinvar_review_status:
                    v.clinvar_review_status = result.clinvar_review_status
                stats["clinvar_filled"] += 1
        
        # Store CADD phred in vcf_info if enabled
        if fill_cadd and result.cadd_phred is not None:
            if not hasattr(v, 'vcf_info') or v.vcf_info is None:
                v.vcf_info = {}
            v.vcf_info["cadd_phred"] = result.cadd_phred
            stats["cadd_filled"] += 1
    
    return stats


# ---------------------------------------------------------------------------
# Test / validation
# ---------------------------------------------------------------------------

async def _test_myvariant():
    """Self-test with a known common variant."""
    async with aiohttp.ClientSession(trust_env=False) as session:
        # Known variant: chr1:218631822 G>A (common in AFR)
        result = await query_myvariant_single("chr1", 218631822, "G", "A", session)
        print(f"Single query test:")
        print(f"  status: {result.status}")
        print(f"  gnomad_af: {result.gnomad_af}")
        print(f"  clinvar: {result.clinvar_significance}")
        print(f"  cadd_phred: {result.cadd_phred}")
        
        # Batch test
        batch = [
            ("chr1", 218631822, "G", "A"),
            ("chr12", 859434, "A", "C"),      # WNK1 — likely not found
            ("chr2", 232480764, "C", "T"),    # ECEL1 — likely not found
        ]
        results = await query_myvariant_batch(batch, session)
        print(f"\nBatch test ({len(batch)} variants):")
        for vid, r in results.items():
            print(f"  {vid}: status={r.status}, gnomad_af={r.gnomad_af}, clinvar={r.clinvar_significance}")


if __name__ == "__main__":
    asyncio.run(_test_myvariant())
