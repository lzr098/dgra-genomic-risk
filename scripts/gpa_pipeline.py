#!/usr/bin/env python3
"""
GPA Pipeline Module

Main analysis pipeline: async orchestration of API calls, tier classification,
and report generation.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
v0.11.1: Re-architected into phase-based workflow engine with JSON checkpoints.
"""

import asyncio
import json
import logging
import re

import aiohttp  # v0.10.3-fix: required for gnomad query exception handling

from dataclasses import asdict
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Any, Tuple

from dgra_core import (
    Variant, GPAConfig, _UNKNOWN,
    _COMMON_TS_GENES, _KNOWN_AML_DRIVERS, _GENE_FAMILY_REDUNDANCY,
    OFFLINE_ARCHIVE_DIR,
    _save_offline_archive, _load_offline_archive,
    assess_tissue_relevance,
    classify_gnomad_frequency, normalize_gene_symbols,
    map_variant_to_domain, evaluate_gene_constraint,
    predict_nmd, evaluate_missense_tier,
    aggregate_gtex_expression,
    _x_linked_female_adjustment,
    detect_pseudogene_artifact,
    _variant_has_pathogenic_evidence,
    correct_transcript_priority,
)
from gpa_tier_classifier import classify_variant_tier
from gpa_phaser import PhaseResult
from gpa_multi_hit import detect_multi_hit_genes
from gpa_qc import _run_qc_checks
from gpa_report import generate_tier_report, generate_json_report
from dgra_cache import DGRACache
from dgra_api import DGRAAPIClient
from api_hub import APIHub
from dgra_hgnc_local import LocalHGNC, local_hgnc_available
from gpa_workflow_engine import WorkflowEngine, WorkflowError
from gpa_audit_trail import AuditTrail


# =============================================================================
# Phase functions for the workflow engine
# =============================================================================

async def _phase_0_preflight_parse(context: Dict[str, Any]) -> Dict[str, Any]:
    """Phase 0: environment preflight + input parsing."""
    variants_data = context["variants_data"]
    user_phenotypes = context.get("user_phenotypes")
    config = context.get("config")
    if config is None:
        config = GPAConfig()

    global_config = config.to_global()

    # Preflight health check
    from gpa_preflight import run_preflight_check, suggest_action
    preflight, route_map = await run_preflight_check(global_config)
    if not preflight.is_ready():
        action = suggest_action(preflight)
        if action == "abort":
            print("[GPA Preflight] 环境检查未通过，中止分析。")
            print(preflight.to_markdown())
            context["output"] = {"error": "Preflight failed", "report": preflight.to_dict()}
            raise WorkflowError("Preflight failed")
        elif action == "offline":
            print("[GPA Preflight] 环境检查未通过，建议切换到离线模式（跳过所有 API 调用）。")
            print("  如需继续离线模式，请显式设置 config.offline_mode=True 后重试。")
            print("  当前默认行为：中止任务，保持在线优先。")
            context["output"] = {
                "error": "Preflight failed — online mode required. Set offline_mode=True explicitly to proceed.",
                "report": preflight.to_dict(),
            }
            raise WorkflowError("Preflight failed")

    tissue_profile = config.get_tissue_profile()
    profile_name = tissue_profile.get("display_name", config.tissue_profile)

    # Parse variants
    variants: List[Variant] = []
    for vd in variants_data:
        missing = []

        raw_impact = str(vd.get("IMPACT", "")).strip()
        _IMPACT_CN_MAP = {"高": "HIGH", "中等": "MODERATE", "低": "LOW", "修饰": "MODIFIER"}
        if raw_impact in _IMPACT_CN_MAP:
            raw_impact = _IMPACT_CN_MAP[raw_impact]
        if not raw_impact:
            raw_impact = _UNKNOWN
            missing.append("IMPACT")

        raw_consequence = str(vd.get("Consequence", "")).strip()
        _CONSEQUENCE_CN_MAP = {
            "错义变异": "missense_variant",
            "无义变异": "stop_gained",
            "获得终止密码子": "stop_gained",
            "移码变异": "frameshift_variant",
            "框内插入": "inframe_insertion",
            "剪接位点变异": "splice_donor_variant",
            "剪接区域变异": "splice_region_variant",
            "剪接供体区域变异": "splice_donor_variant",
            "剪接供体第5位碱基变异": "splice_donor_variant",
            "剪接多嘧啶束变异": "splice_polypyrimidine_tract_variant",
            "内含子变异": "intron_variant",
            "基因上游变异": "upstream_gene_variant",
            "基因下游变异": "downstream_gene_variant",
            "同义变异": "synonymous_variant",
            "非翻译区变异": "UTR_variant",
            "3'非翻译区变异": "3_prime_UTR_variant",
            "5'非翻译区变异": "5_prime_UTR_variant",
            "非编码转录本外显子变异": "non_coding_transcript_exon_variant",
        }
        if raw_consequence in _CONSEQUENCE_CN_MAP:
            raw_consequence = _CONSEQUENCE_CN_MAP[raw_consequence]
        if not raw_consequence:
            raw_consequence = _UNKNOWN
            missing.append("Consequence")

        raw_clinvar = str(vd.get("CLIN_SIG", "")).strip()
        if not raw_clinvar:
            raw_clinvar = _UNKNOWN
            missing.append("CLIN_SIG")

        raw_clinvar_review = str(vd.get("CLNREVSTAT", "")).strip()
        if not raw_clinvar_review or raw_clinvar_review == "nan":
            raw_clinvar_review = None

        raw_dp = str(vd.get("DP", "")).strip()
        try:
            dp_val = int(float(raw_dp)) if raw_dp and raw_dp != _UNKNOWN else 0
        except ValueError:
            dp_val = 0
        if not raw_dp:
            missing.append("DP")

        raw_gq = str(vd.get("GQ", "")).strip()
        try:
            gq_val = float(raw_gq) if raw_gq and raw_gq not in ("", _UNKNOWN, "None") else 0.0
        except (ValueError, TypeError):
            gq_val = 0.0
        if not raw_gq:
            missing.append("GQ")

        raw_vaf = str(vd.get("VAF", "")).strip()
        if not raw_vaf or raw_vaf == ".":
            vaf_val = None
        else:
            try:
                vaf_val = float(raw_vaf)
            except (ValueError, TypeError):
                vaf_val = None
        if not raw_vaf or raw_vaf == ".":
            missing.append("VAF")

        raw_gnomad = str(vd.get("gnomAD_AF", "")).strip()
        gnomad_val = float(raw_gnomad) if raw_gnomad and raw_gnomad != "N/A" else None
        if not raw_gnomad:
            missing.append("gnomAD_AF")

        quality_confidence = "high"
        if missing:
            n_critical = sum(1 for f in missing if f in ("IMPACT", "Consequence", "VAF", "CLIN_SIG"))
            if n_critical >= 3:
                quality_confidence = "unknown"
            elif n_critical >= 1:
                quality_confidence = "low"
            else:
                quality_confidence = "medium"

        v = Variant(
            chrom=vd.get("CHROM", ""),
            pos=int(vd.get("POS", 0) or 0),
            ref=vd.get("REF", ""),
            alt=vd.get("ALT", ""),
            gene=vd.get("GENE", "") or vd.get("Gene", ""),
            transcript=vd.get("Feature", ""),
            exon=vd.get("EXON", ""),
            impact=raw_impact,
            consequence=raw_consequence,
            hgvsp=vd.get("HGVSp", ""),
            hgvsc=vd.get("HGVSc", ""),
            clinvar=raw_clinvar,
            dp=dp_val,
            gq=gq_val,
            gt=vd.get("GT", ""),
            vaf=vaf_val,
            gnomad_af=gnomad_val,
            classification=vd.get("classification", ""),
            is_tsg=vd.get("is_tsg", "") == "Yes" or vd.get("is_tsg", False) == True,
            is_oncogene=vd.get("is_oncogene", "") == "Yes" or vd.get("is_oncogene", False) == True,
            clinvar_review_status=raw_clinvar_review,
            quality_confidence=quality_confidence,
            missing_fields=missing,
        )
        variants.append(v)

    # Annotation quality gate
    if variants:
        total = len(variants)
        n_no_gene = sum(1 for v in variants if not v.gene or v.gene == _UNKNOWN)
        n_no_impact = sum(1 for v in variants if v.impact == _UNKNOWN)
        n_no_consequence = sum(1 for v in variants if v.consequence == _UNKNOWN)
        if n_no_gene >= total * 0.8 and n_no_impact >= total * 0.8 and n_no_consequence >= total * 0.8:
            raise ValueError(
                f"CRITICAL: Input appears to be an unannotated raw VCF. "
                f"{n_no_gene}/{total} variants missing Gene, "
                f"{n_no_impact}/{total} missing IMPACT, "
                f"{n_no_consequence}/{total} missing Consequence. "
                f"Raw VCF must be annotated via VCFAnnotator before pipeline entry. "
                f"Use run_gpa_from_file() or dgra_core.py --input for VCF inputs."
            )

    unique_genes = list({v.gene for v in variants})
    print(f"[GPA] {len(variants)} variants across {len(unique_genes)} unique genes")
    print(f"[GPA] Tissue profile: {profile_name} | Offline: {config.offline_mode}")

    context.update({
        "config": config,
        "global_config": global_config,
        "preflight_report": preflight.to_dict(),
        "proxy_route_map": route_map,
        "tissue_profile": tissue_profile,
        "profile_name": profile_name,
        "variants": variants,
        "unique_genes": unique_genes,
        "user_phenotypes": user_phenotypes,
    })
    return context


async def _phase_1_api_queries(context: Dict[str, Any]) -> Dict[str, Any]:
    """Phase 1: gene-level and variant-level API enrichment."""
    config = context["config"]
    global_config = context["global_config"]
    variants = context["variants"]
    unique_genes = context["unique_genes"]
    audit_trail = context.get("audit_trail")

    ensembl_data: Dict[str, Any] = {}
    uniprot_data: Dict[str, Any] = {}
    gtex_data: Dict[str, Any] = {}
    hgnc_data: Dict[str, Any] = {}
    gnomad_constraint_data: Dict[str, Any] = {}

    if not config.offline_mode and unique_genes:
        cache = DGRACache(global_config.cache_db_path)
        async with DGRAAPIClient(global_config, cache, audit_trail=audit_trail) as client:
            gtex_data = {}
            try:
                import sys
                gtex_local_path = str(Path.home() / ".workbuddy" / "scripts")
                if gtex_local_path not in sys.path:
                    sys.path.insert(0, gtex_local_path)
                from gtex_local import query_gtex_local
                for gene in unique_genes:
                    local_expr = query_gtex_local(gene, tissues=None)
                    if local_expr:
                        all_expressing = sorted(
                            [(t, v) for t, v in local_expr.items() if v > 0],
                            key=lambda x: x[1], reverse=True
                        )
                        max_tpm = max(local_expr.values()) if local_expr else 0.0
                        gtex_data[gene] = {
                            "median_tpm": max_tpm,
                            "max_tpm": max_tpm,
                            "all_tissues": [{"tissue": t, "tpm": v} for t, v in all_expressing],
                            "expressing_tissues": len(all_expressing),
                            "source": "gtex_local_db",
                            "confidence": "medium",
                        }
                print(f"[GPA] GTEx local DB: expression data for {len(gtex_data)}/{len(unique_genes)} genes")
            except Exception as e:
                print(f"[GPA] GTEx local DB query failed: {e}")

            ensembl_raw, uniprot_raw, hgnc_raw, gnomad_constraint_raw = await asyncio.gather(
                client.batch_query_genes(unique_genes, "ensembl"),
                client.batch_query_genes(unique_genes, "uniprot"),
                client.batch_query_genes(unique_genes, "hgnc"),
                client.batch_query_genes(unique_genes, "gnomad_constraint"),
            )
            ensembl_data = {g: ensembl_raw.get(g, {}) for g in unique_genes}
            uniprot_data = {g: uniprot_raw.get(g, {}) for g in unique_genes}
            hgnc_data = {g: hgnc_raw.get(g, {}) for g in unique_genes}
            gnomad_constraint_data = {g: gnomad_constraint_raw.get(g, {}) for g in unique_genes}

            # MyVariant.info batch query for gnomAD/CADD
            variants_needing_enrichment = [
                v for v in variants
                if v.gnomad_af is None and v.chrom and v.pos and v.ref and v.alt
            ]
            if variants_needing_enrichment:
                print(f"[GPA] MyVariant.info: batch querying {len(variants_needing_enrichment)} variants for gnomAD/CADD")
                try:
                    from dgra_myvariant import query_myvariant_batch, apply_myvariant_results
                    mv_sem = asyncio.Semaphore(10)
                    timeout_obj = aiohttp.ClientTimeout(total=300)
                    mv_variants = [(v.chrom, v.pos, v.ref, v.alt) for v in variants_needing_enrichment]
                    mv_results = await query_myvariant_batch(mv_variants, client.hub.session, semaphore=mv_sem, batch_size=1000)
                    mv_stats = apply_myvariant_results(variants, mv_results)
                    print(f"[GPA] MyVariant.info: {mv_stats['gnomad_filled']} gnomAD, {mv_stats['clinvar_filled']} ClinVar, {mv_stats['cadd_filled']} CADD filled | {mv_stats['not_found']} not_found, {mv_stats['errors']} errors")
                except Exception as e:
                    print(f"[GPA] MyVariant.info batch query failed (non-critical, falling back): {type(e).__name__}: {e}")

            # gnomAD variant frequency batch query
            variants_without_af = [v for v in variants if v.gnomad_af is None and v.chrom and v.pos and v.ref and v.alt]
            if variants_without_af:
                print(f"[GPA] gnomAD: querying {len(variants_without_af)} variants without AF data")
                gnomad_sem = asyncio.Semaphore(2)
                async def _query_one_gnomad(v):
                    async with gnomad_sem:
                        try:
                            return await client.query_gnomad_variant(v.chrom, v.pos, v.ref, v.alt)
                        except asyncio.TimeoutError as e:
                            print(f"[GPA] gnomAD query TIMEOUT for {v.gene} {v.chrom}:{v.pos}: {e}")
                            return {"status": "API_FAILED", "error": f"timeout: {e}", "source": "failed"}
                        except aiohttp.ClientError as e:
                            print(f"[GPA] gnomAD query CLIENT_ERROR for {v.gene} {v.chrom}:{v.pos}: {e}")
                            return {"status": "API_FAILED", "error": f"client_error: {e}", "source": "failed"}
                        except Exception as e:
                            print(f"[GPA] gnomAD query FAILED for {v.gene} {v.chrom}:{v.pos}: {e}")
                            return {"status": "API_FAILED", "error": str(e), "source": "failed"}

                gnomad_results = await asyncio.gather(*[_query_one_gnomad(v) for v in variants_without_af])
                n_success = n_failed = n_not_captured = n_query_error = n_myvariant_fallback = 0
                _myvariant_fallback_variants = []

                for v, result in zip(variants_without_af, gnomad_results):
                    if result and result.get("source") in ("gnomad", "cache", "failed"):
                        af = result.get("af")
                        gnomad_status = result.get("status", "NOT_CAPTURED")
                        if af is not None:
                            v.gnomad_af = af
                            v.gnomad_populations = result.get("af_populations", {})
                            v.gnomad_status = "SUCCESS"
                            n_success += 1
                            print(f"[GPA] gnomAD: {v.gene} {v.chrom}:{v.pos} AF={v.gnomad_af}")
                        elif gnomad_status == "QUERY_ERROR":
                            v.gnomad_status = "QUERY_ERROR"
                            v.gnomad_error_msg = result.get("note", "gnomAD GraphQL query error")
                            v.gnomad_af_warning = True
                            v.gnomad_populations = result.get("af_populations", {})
                            n_query_error += 1
                            _myvariant_fallback_variants.append(v)
                            print(f"[GPA] gnomAD: {v.gene} {v.chrom}:{v.pos} QUERY_ERROR → will try MyVariant fallback")
                        else:
                            v.gnomad_populations = {}
                            v.gnomad_status = result.get("status", "NOT_CAPTURED")
                            n_not_captured += 1
                            _myvariant_fallback_variants.append(v)
                            print(f"[GPA] gnomAD: {v.gene} {v.chrom}:{v.pos} NOT_CAPTURED → will try MyVariant fallback")
                    elif result and result.get("status") == "API_FAILED":
                        v.gnomad_status = "API_FAILED"
                        v.gnomad_error_msg = result.get("error", "unknown")
                        v.gnomad_af_warning = True
                        v.gnomad_populations = {}
                        n_failed += 1
                        print(f"[GPA] gnomAD: {v.gene} {v.chrom}:{v.pos} API_FAILED ({v.gnomad_error_msg})")

                # MyVariant.info fallback for NOT_CAPTURED/QUERY_ERROR
                if _myvariant_fallback_variants:
                    print(f"[GPA] MyVariant.info fallback: querying {len(_myvariant_fallback_variants)} variants not found in gnomAD GraphQL")
                    try:
                        from dgra_myvariant import query_myvariant_batch, apply_myvariant_results
                        mv_fb_sem = asyncio.Semaphore(10)
                        mv_fb_variants = [(v.chrom, v.pos, v.ref, v.alt) for v in _myvariant_fallback_variants]
                        for v in _myvariant_fallback_variants:
                            v.gnomad_af = None
                        mv_fb_results = await query_myvariant_batch(mv_fb_variants, client.hub.session, semaphore=mv_fb_sem, batch_size=1000)
                        mv_fb_stats = apply_myvariant_results(_myvariant_fallback_variants, mv_fb_results)
                        n_myvariant_fallback = mv_fb_stats.get("gnomad_filled", 0)
                        for v in _myvariant_fallback_variants:
                            if v.gnomad_af is not None:
                                v.gnomad_status = "MYVARIANT_FALLBACK"
                                v.gnomad_af_warning = False
                        print(f"[GPA] MyVariant.info fallback: {n_myvariant_fallback} variants filled with gnomAD data")
                    except Exception as e:
                        print(f"[GPA] MyVariant.info fallback failed (non-critical): {type(e).__name__}: {e}")

                print(f"[GPA] gnomAD results: {n_success} success, {n_not_captured} not in dataset, {n_query_error} query errors, {n_failed} API failures | {n_myvariant_fallback} recovered via MyVariant fallback")

                # Retry API_FAILED variants once with longer timeout
                _retry_variants = [v for v in variants_without_af if v.gnomad_status == "API_FAILED"]
                if _retry_variants:
                    print(f"[GPA] gnomAD retry: re-querying {len(_retry_variants)} variants with API_FAILED status (longer timeout)")
                    gnomad_retry_sem = asyncio.Semaphore(1)
                    async def _retry_one_gnomad(v):
                        async with gnomad_retry_sem:
                            try:
                                return await asyncio.wait_for(
                                    client.query_gnomad_variant(v.chrom, v.pos, v.ref, v.alt),
                                    timeout=60
                                )
                            except asyncio.TimeoutError as e:
                                print(f"[GPA] gnomAD retry TIMEOUT for {v.gene} {v.chrom}:{v.pos}: {e}")
                                return {"status": "API_FAILED", "error": f"retry_timeout: {e}", "source": "failed"}
                            except Exception as e:
                                print(f"[GPA] gnomAD retry FAILED for {v.gene} {v.chrom}:{v.pos}: {e}")
                                return {"status": "API_FAILED", "error": f"retry_failed: {e}", "source": "failed"}

                    retry_results = await asyncio.gather(*[_retry_one_gnomad(v) for v in _retry_variants])
                    n_retry_success = n_retry_failed = 0
                    for v, result in zip(_retry_variants, retry_results):
                        if result and result.get("source") in ("gnomad", "cache", "failed"):
                            af = result.get("af")
                            if af is not None:
                                v.gnomad_af = af
                                v.gnomad_populations = result.get("af_populations", {})
                                v.gnomad_status = "SUCCESS"
                                v.gnomad_af_warning = False
                                n_retry_success += 1
                                print(f"[GPA] gnomAD retry SUCCESS: {v.gene} {v.chrom}:{v.pos} AF={v.gnomad_af}")
                            else:
                                v.gnomad_status = result.get("status", "NOT_CAPTURED")
                                n_retry_failed += 1
                        else:
                            v.gnomad_status = "API_FAILED"
                            v.gnomad_error_msg = result.get("error", "retry_failed")
                            n_retry_failed += 1
                    print(f"[GPA] gnomAD retry results: {n_retry_success} success, {n_retry_failed} still failed")

                    def _is_tier_candidate(v):
                        impact = str(v.impact or "").upper()
                        clinvar = str(v.clinvar or "").lower()
                        is_high = impact == "HIGH" or impact == "UNKNOWN"
                        is_pathogenic = "pathogenic" in clinvar and "conflicting" not in clinvar
                        return is_high or is_pathogenic

                    _still_failed_tier_candidates = [
                        v for v in _retry_variants
                        if v.gnomad_status == "API_FAILED" and _is_tier_candidate(v)
                    ]
                    if _still_failed_tier_candidates:
                        print(f"\n{'='*60}")
                        print("⚠️  USER ALERT: Tier 1/2 candidate variants with gnomAD query failure")
                        print(f"   {len(_still_failed_tier_candidates)} variants could not be verified in gnomAD")
                        print("   These variants have been conservatively downgraded to Tier 2")
                        print("   MANUAL ACTION REQUIRED: Please verify gnomAD/ClinVar status externally")
                        print(f"{'='*60}\n")

            # Direct NCBI ClinVar query
            variants_needing_clinvar = [
                v for v in variants
                if v.clinvar == _UNKNOWN or v.clinvar == "UNKNOWN"
            ]
            if variants_needing_clinvar:
                from dgra_clinvar import variant_needs_clinvar, query_clinvar_batch, apply_clinvar_results
                clinvar_candidates = [v for v in variants_needing_clinvar if variant_needs_clinvar(v.consequence)]
                if clinvar_candidates:
                    n_cv = len(clinvar_candidates)
                    if n_cv <= 1000:
                        print(f"[GPA] ClinVar (NCBI): querying {n_cv} variants (1 req/s, consequence-filtered from {len(variants_needing_clinvar)})")
                        try:
                            cv_variants = [(v.chrom, v.pos, v.ref, v.alt, v.gene, v.consequence) for v in clinvar_candidates]
                            cv_sem = asyncio.Semaphore(1)
                            cv_results = await query_clinvar_batch(cv_variants, client.hub.session, semaphore=cv_sem, rate_limit_delay=1.0)
                            cv_stats = apply_clinvar_results(variants, cv_results)
                            print(f"[GPA] ClinVar (NCBI): {cv_stats['filled']} filled, {cv_stats['not_found']} not_found, {cv_stats['skipped']} skipped, {cv_stats['errors']} errors")
                        except Exception as e:
                            print(f"[GPA] ClinVar (NCBI) query failed (non-critical): {type(e).__name__}: {e}")
                    else:
                        print(f"[GPA] ClinVar: {n_cv} candidates > 1000 threshold — skipping direct NCBI query. "
                              f"ClinVar data will come from MyVariant.info batch (if variant also lacks gnomAD AF). "
                              f"For more accurate ClinVar annotation, re-run with offline_mode=True and provide a pre-annotated VCF.")

        print(f"[GPA] API batch query complete: Ensembl={len(ensembl_data)}, UniProt={len(uniprot_data)}, GTEx={len(gtex_data)}, HGNC={len(hgnc_data)}, gnomAD_constraint={len(gnomad_constraint_data)}")
        for gene in unique_genes:
            _save_offline_archive(gene, ensembl_data, uniprot_data, gtex_data, config.tissue_profile, gnomad_constraint_data)
        print(f"[GPA] Offline archive saved for {len(unique_genes)} genes to {OFFLINE_ARCHIVE_DIR}")
    else:
        loaded = 0
        for gene in unique_genes:
            archive = _load_offline_archive(gene)
            if archive:
                ensembl_data[gene] = archive.get("ensembl", {})
                uniprot_data[gene] = archive.get("uniprot", {})
                gtex_data[gene] = archive.get("gtex", {})
                gc = archive.get("gnomad_constraint")
                if gc and gc.get("status") == "CAPTURED":
                    gnomad_constraint_data[gene] = gc
                loaded += 1
        print(f"[GPA] Offline mode: loaded archived data for {loaded}/{len(unique_genes)} genes from {OFFLINE_ARCHIVE_DIR}")
        if loaded == 0:
            print("[GPA] Offline mode: no archive found, using local fallbacks only (conservative)")

        hgnc_data = {}
        try:
            if local_hgnc_available():
                local_hgnc = LocalHGNC()
                hgnc_data = local_hgnc.batch_resolve(list(unique_genes))
                meta = local_hgnc.metadata()
                print(
                    f"[GPA] Offline mode: loaded local HGNC lookup "
                    f"({meta.get('total_symbols', 0):,} symbols, "
                    f"{meta.get('approved_count', 0):,} approved)"
                )
            else:
                print("[GPA] Offline mode: local HGNC lookup not found, using conservative known-symbol fallback")
        except Exception as e:
            print(f"[GPA] Offline mode local HGNC lookup failed (non-critical): {e}")

        try:
            from dgra_clinvar import _try_local_clinvar_query
            variants_needing_clinvar = [v for v in variants if v.clinvar == _UNKNOWN or v.clinvar == "UNKNOWN"]
            if variants_needing_clinvar:
                print(f"[GPA] Offline mode: querying local ClinVar VCF for {len(variants_needing_clinvar)} variants...")
                filled = 0
                not_found = 0
                for v in variants_needing_clinvar:
                    local_result = _try_local_clinvar_query(v.chrom, v.pos, v.ref, v.alt)
                    if local_result:
                        v.clinvar = local_result.get("CLNSIG", _UNKNOWN)
                        v.clinvar_review_status = local_result.get("CLNREVSTAT")
                        v.clinvar_accession = local_result.get("ID")
                        filled += 1
                    else:
                        not_found += 1
                print(f"[GPA] Offline mode local ClinVar: {filled} filled, {not_found} not_found")
        except Exception as e:
            print(f"[GPA] Offline mode local ClinVar query failed (non-critical): {e}")

    context.update({
        "ensembl_data": ensembl_data,
        "uniprot_data": uniprot_data,
        "gtex_data": gtex_data,
        "hgnc_data": hgnc_data,
        "gnomad_constraint_data": gnomad_constraint_data,
    })
    return context


async def _phase_2_post_processing(context: Dict[str, Any]) -> Dict[str, Any]:
    """Phase 2: HGNC normalization, NMD, transcript correction, domain/tissue mapping, multi-hit, QC."""
    config = context["config"]
    variants = context["variants"]
    ensembl_data = context["ensembl_data"]
    uniprot_data = context["uniprot_data"]
    gtex_data = context["gtex_data"]
    hgnc_data = context["hgnc_data"]
    gnomad_constraint_data = context["gnomad_constraint_data"]
    tissue_profile = context["tissue_profile"]

    # HGNC normalization
    hgnc_warnings = normalize_gene_symbols(variants, hgnc_data, offline_mode=config.offline_mode)
    if hgnc_warnings:
        print(f"[GPA] HGNC normalization: {len(hgnc_warnings)} warnings")
        for w in hgnc_warnings[:5]:
            print(f"  - {w}")
        if len(hgnc_warnings) > 5:
            print(f"  ... and {len(hgnc_warnings) - 5} more")

    # Gene constraint
    for v in variants:
        gc = gnomad_constraint_data.get(v.gene, {})
        if gc and gc.get("status") == "CAPTURED":
            v.gene_constraint = {
                "pLI": gc.get("pLI"),
                "loeuf": gc.get("loeuf"),
                "lof_z": gc.get("lof_z"),
                "mis_z": gc.get("mis_z"),
                "oe_lof": gc.get("oe_lof"),
                "source": gc.get("source", "gnomad"),
            }

    # NMD prediction
    nmd_count = 0
    for v in variants:
        lof_terms = {"frameshift", "nonsense", "stop_gained", "start_lost"}
        is_truncating = any(term in v.consequence.lower() for term in lof_terms)
        if is_truncating:
            nmd_result = predict_nmd(v, ensembl_data.get(v.gene, {}) if ensembl_data else None)
            v.nmd_prediction = nmd_result
            if v.gene_constraint is None:
                v.gene_constraint = {}
            v.gene_constraint["nmd_prediction"] = nmd_result
            nmd_count += 1
    if nmd_count > 0:
        print(f"[GPA] NMD prediction computed for {nmd_count} truncating variants")

    # Transcript correction
    for v in variants:
        v, warning = await correct_transcript_priority(v, ensembl_data)
        if warning:
            v.transcript_warning = json.dumps(warning)

    # VEP reannotation for TRANSCRIPT_DISCREPANCY variants
    discrepancy_variants = []
    for v in variants:
        if v.transcript_warning:
            try:
                tw = json.loads(v.transcript_warning)
            except (json.JSONDecodeError, TypeError):
                tw = {}
            if tw.get("type") == "TRANSCRIPT_DISCREPANCY":
                discrepancy_variants.append(v)

    if discrepancy_variants and not config.offline_mode:
        print(f"[GPA] VEP reannotation: {len(discrepancy_variants)} variants with TRANSCRIPT_DISCREPANCY")
        vep_inputs = []
        for v in discrepancy_variants:
            vep_inputs.append({
                "chrom": v.chrom,
                "pos": v.pos,
                "ref": v.ref,
                "alt": v.alt,
                "key": f"{v.chrom}:{v.pos}_{v.ref}>{v.alt}",
            })

        global_config = config.to_global()
        cache = DGRACache(global_config.cache_db_path)
        async with DGRAAPIClient(global_config, cache) as client:
            vep_results = await client.batch_query_vep_region(vep_inputs)

        updated_count = 0
        failed_count = 0
        for v in discrepancy_variants:
            key = f"{v.chrom}:{v.pos}_{v.ref}>{v.alt}"
            vep_result = vep_results.get(key, {})

            if vep_result.get("error"):
                failed_count += 1
                try:
                    tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
                except (json.JSONDecodeError, TypeError):
                    tw = {}
                tw["vep_reannotation"] = {
                    "status": "failed",
                    "error": vep_result.get("error"),
                    "fallback": "Using annotator's original annotation",
                }
                tw["vep_reannotation_failed"] = True
                v.transcript_warning = json.dumps(tw)
                v.quality_confidence = "LOW"
                v.tier_confidence = "LOW"
                continue

            original = {
                "consequence": v.consequence,
                "impact": v.impact,
                "hgvsc": v.hgvsc,
                "hgvsp": v.hgvsp,
                "transcript": v.transcript,
            }

            consequence_terms = vep_result.get("consequence_terms", [])
            if consequence_terms:
                v.consequence = consequence_terms[0].replace("_", " ")
            if vep_result.get("impact"):
                v.impact = vep_result["impact"]
            if vep_result.get("hgvsc"):
                v.hgvsc = vep_result["hgvsc"]
            if vep_result.get("hgvsp"):
                v.hgvsp = vep_result["hgvsp"]
            if vep_result.get("transcript_id"):
                v.transcript = vep_result["transcript_id"]

            try:
                tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
            except (json.JSONDecodeError, TypeError):
                tw = {}
            tw["vep_reannotation"] = {
                "status": "success",
                "original": original,
                "canonical": {
                    "consequence": v.consequence,
                    "impact": v.impact,
                    "hgvsc": v.hgvsc,
                    "hgvsp": v.hgvsp,
                    "transcript": v.transcript,
                    "transcript_id": vep_result.get("transcript_id"),
                    "protein_domains": vep_result.get("protein_domains", []),
                },
                "source": vep_result.get("source", "ensembl"),
                "confidence": vep_result.get("confidence", "medium"),
            }
            v.transcript_warning = json.dumps(tw)
            updated_count += 1

        print(f"[GPA] VEP reannotation complete: {updated_count} updated, {failed_count} failed")

    elif discrepancy_variants and config.offline_mode:
        print(f"[GPA] VEP reannotation skipped: offline mode ({len(discrepancy_variants)} discrepancies)")
        for v in discrepancy_variants:
            try:
                tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
            except (json.JSONDecodeError, TypeError):
                tw = {}
            tw["vep_reannotation"] = {
                "status": "skipped",
                "reason": "offline_mode",
                "fallback": "Using annotator's original annotation",
            }
            tw["vep_reannotation_failed"] = True
            v.transcript_warning = json.dumps(tw)
            v.quality_confidence = "LOW"
            v.tier_confidence = "LOW"

    # Pseudogene detection
    for v in variants:
        pg_warning = detect_pseudogene_artifact(v)
        if pg_warning:
            v.pseudogene_warning = json.dumps(pg_warning)

    # gnomAD classification
    for v in variants:
        if v.gnomad_status == "API_FAILED":
            continue
        gnomad_info = classify_gnomad_frequency(
            v.gnomad_af, v.gene,
            af_by_population=v.gnomad_populations,
            target_population=getattr(config, 'target_population', None)
        )
        v.gnomad_status = gnomad_info["status"]

    # Protein domain mapping
    for v in variants:
        v.domain_info = map_variant_to_domain(v, uniprot_data)

    # Tissue relevance
    tissue_assessments = {}
    for v in variants:
        tissue = assess_tissue_relevance(v, tissue_profile, gtex_data)
        tissue_assessments[v.gene] = tissue
        v.tissue_relevance = tissue

    # Multi-hit detection
    multi_hits = detect_multi_hit_genes(variants, gtex_data)

    # QC checks
    qc_summary = _run_qc_checks(variants)
    if qc_summary["flagged"] > 0:
        print(f"[GPA] QC: {qc_summary['flagged']}/{qc_summary['total']} variants flagged: {qc_summary['by_flag']}")

    context.update({
        "tissue_assessments": tissue_assessments,
        "multi_hits": multi_hits,
        "qc_summary": qc_summary,
    })
    return context


async def _phase_3_classification_report(context: Dict[str, Any]) -> Dict[str, Any]:
    """Phase 3: phenotype matching, SpliceAI, tier classification, adjustments, report generation."""
    config = context["config"]
    variants = context["variants"]
    audit_trail = context.get("audit_trail")
    tissue_profile = context["tissue_profile"]
    tissue_assessments = context["tissue_assessments"]
    multi_hits = context["multi_hits"]
    qc_summary = context["qc_summary"]
    gtex_data = context["gtex_data"]
    user_phenotypes = context.get("user_phenotypes")
    profile_name = context["profile_name"]
    unique_genes = context["unique_genes"]
    ensembl_data = context["ensembl_data"]
    uniprot_data = context["uniprot_data"]

    # Phenotype association
    if user_phenotypes:
        from gpa_phenotype_match import PhenotypeMatcher
        matcher = PhenotypeMatcher(offline_mode=config.offline_mode)
        gene_symbols = [v.gene for v in variants]
        match_results = await matcher.match_batch(gene_symbols, user_phenotypes)
        for v, mr in zip(variants, match_results):
            v.phenotype_match_score = mr.get("score")
            v.phenotype_match_explanation = mr.get("explanation", "")
            v.phenotype_match_confidence = mr.get("confidence", "")
            v.phenotype_matched_pairs = mr.get("matched_pairs", [])
            v.phenotype_known_list = mr.get("known_phenotypes", [])

    # SpliceAI (optional)
    if getattr(config, 'spliceai_enabled', False):
        from dgra_splice_predictor import (
            query_spliceai_batch, should_query_spliceai, reset_spliceai_cache
        )
        reset_spliceai_cache()
        spliceai_sem = asyncio.Semaphore(getattr(config, 'spliceai_concurrency', 5))
        spliceai_candidates = [v for v in variants if should_query_spliceai(v.consequence)]
        if spliceai_candidates:
            print(f"[GPA] SpliceAI: querying {len(spliceai_candidates)} splice variants (concurrency={getattr(config, 'spliceai_concurrency', 5)})")
            spliceai_results = await query_spliceai_batch(
                spliceai_candidates, spliceai_sem,
                timeout=getattr(config, 'spliceai_timeout', 45),
            )
            for v in variants:
                from dgra_splice_predictor import _cache_key as _splice_key
                key = _splice_key(v.chrom, v.pos, v.ref, v.alt)
                if key in spliceai_results:
                    v.spliceai_result = spliceai_results[key]
                elif should_query_spliceai(v.consequence):
                    v.spliceai_result = {"source": "not_in_db", "delta_score": None, "predicted_impact": None}
            print("[GPA] SpliceAI: batch complete")

    # Tier classification
    for v in variants:
        tissue = tissue_assessments[v.gene]
        gnomad_info = classify_gnomad_frequency(
            v.gnomad_af, v.gene,
            af_by_population=v.gnomad_populations,
            target_population=getattr(config, 'target_population', None)
        )
        tw = json.loads(v.transcript_warning) if v.transcript_warning else None
        pw = json.loads(v.pseudogene_warning) if v.pseudogene_warning else None

        tier, reason, actions = classify_variant_tier(
            v, v.domain_info, tissue, gnomad_info, tw, pw, tissue_profile, config
        )
        v.tier = tier
        v.tier_reason = reason
        v.tier_actions = actions

        if "HLA_REGION_SHORT_READ_UNRELIABLE" in v.qc_flags and v.tier <= 2:
            v.tier_actions.append(
                "⚠️ HLA_REGION: 该变异位于 HLA 高多态区域（chr6:29.7M-33.1M），"
                "短读测序数据在此区域可能不可靠（比对伪影、参考基因组偏差）。"
                "建议通过 HLA 分型工具（HLA-HD、OptiType）或长读测序验证。"
            )

    # HLA multi-hit filtering
    _HLA_GENES = {
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
        "HLA-DPA1", "HLA-DPB1", "HLA-E", "HLA-F", "HLA-G", "HLA-H",
        "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S",
        "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
        "MICA", "MICB", "TAP1", "TAP2",
    }
    multi_hit_genes = {mh["gene"] for mh in multi_hits}
    hla_multi_hits = {g for g in multi_hit_genes if g in _HLA_GENES}
    multi_hit_genes -= hla_multi_hits

    # X-linked female adjustment
    for v in variants:
        adj_tier, adj_reason = _x_linked_female_adjustment(
            v.tier, v.chrom, v.gt, v.gene_constraint
        )
        if adj_reason:
            v.tier = adj_tier
            v.tier_reason += f" | {adj_reason}"

    # Gene family redundancy
    for v in variants:
        if v.gene in _GENE_FAMILY_REDUNDANCY and v.tier == 1:
            redundancy = _GENE_FAMILY_REDUNDANCY[v.gene]
            if redundancy.get("compensation_level") == "complete":
                v.tier = 2
                v.tier_reason += f" | REDUCED: {redundancy['reason']} - complete paralog compensation"
            elif redundancy.get("compensation_level") == "partial":
                v.tier_reason += f" | NOTE: {redundancy['reason']} - partial compensation may mitigate risk"

    # Report generation
    report_md = generate_tier_report(variants, config, tissue_profile, multi_hits, gtex_data=gtex_data)
    json_report = generate_json_report(variants, config, tissue_profile, multi_hits, report_md, qc_summary)

    if audit_trail:
        audit_trail.record_metric("total_variants", len(variants))
        audit_trail.record_metric("tier1_count", len([v for v in variants if v.tier == 1]))
        audit_trail.record_metric("tier2_count", len([v for v in variants if v.tier == 2]))
        audit_trail.record_metric("tier3_count", len([v for v in variants if v.tier == 3]))
        audit_trail.record_metric("multi_hit_genes", len(multi_hits))
        for v in variants:
            if v.tier <= 2:
                audit_trail.record_decision(
                    decision=f"Tier {v.tier}",
                    reason=v.tier_reason or "",
                    subject=f"{v.gene} {v.chrom}:{v.pos} {v.ref}>{v.alt}",
                )

    def _count_valid(data_dict):
        return sum(1 for d in data_dict.values() if d and d.get("source") not in ("failed", None))

    output = {
        "meta": {
            "tissue_profile": config.tissue_profile,
            "profile_display_name": profile_name,
            "analysis_date": datetime.now().isoformat(),
            "total_variants": len(variants),
            "offline_mode": config.offline_mode,
            "api_coverage": {
                "genes_queried": len(unique_genes),
                "ensembl_success": _count_valid(ensembl_data),
                "uniprot_success": _count_valid(uniprot_data),
                "gtex_success": _count_valid(gtex_data),
            }
        },
        "summary": {
            "tier1_gene_count": len(set(v.gene for v in variants if v.tier == 1)),
            "tier1_variant_count": len([v for v in variants if v.tier == 1]),
            "tier2_gene_count": len(set(v.gene for v in variants if v.tier == 2)),
            "tier2_variant_count": len([v for v in variants if v.tier == 2]),
            "tier3_gene_count": len(set(v.gene for v in variants if v.tier == 3)),
            "tier3_variant_count": len([v for v in variants if v.tier == 3]),
            "multi_hit_genes": [mh["gene"] for mh in multi_hits],
        },
        "tier1_variants": [asdict(v) for v in variants if v.tier == 1],
        "tier2_variants": [asdict(v) for v in variants if v.tier == 2],
        "tier3_variants": [asdict(v) for v in variants if v.tier == 3],
        "multi_hit_details": multi_hits,
        "report_markdown": report_md,
        "json_report": json_report,
        "qc_summary": qc_summary,
    }

    context["output"] = output
    return context


async def run_dgra_pipeline(variants_data: List[Dict],
                      user_phenotypes: Optional[str] = None,
                      config: Optional[GPAConfig] = None,
                      enable_audit: Optional[bool] = None) -> Dict:
    """
    Main GPA analysis pipeline with dynamic tissue context.
    v0.4: Async, batch API queries with cache.
    v0.7: Optional phenotype association for Tier 1/2 variants.
    v0.11.1: Phase-based workflow engine with JSON checkpoints.

    Args:
        variants_data: List of variant dicts from VCF annotation
        user_phenotypes: Optional clinical phenotype description
        config: GPA configuration

    Returns:
        Dict with report and structured results
    """
    if config is None:
        config = GPAConfig()

    global_config = config.to_global()
    checkpoint_dir = Path(global_config.cache_db_path).parent / "gpa_checkpoints"
    resume = getattr(config, "resume", False)

    # Audit trail is opt-in via enable_audit argument or config.enable_audit
    if enable_audit is None:
        enable_audit = getattr(config, "enable_audit", False)
    audit_output_dir = Path(global_config.cache_db_path).parent / "gpa_audit"
    audit_trail = AuditTrail(output_dir=audit_output_dir, enabled=enable_audit)

    phases = [
        ("phase_0_preflight_parse", _phase_0_preflight_parse),
        ("phase_1_api_queries", _phase_1_api_queries),
        ("phase_2_post_processing", _phase_2_post_processing),
        ("phase_3_classification_report", _phase_3_classification_report),
    ]

    engine = WorkflowEngine(
        phases=phases,
        checkpoint_dir=checkpoint_dir,
        resume=resume,
        keep_checkpoints=False,
        audit_trail=audit_trail,
    )

    initial_context = {
        "variants_data": variants_data,
        "user_phenotypes": user_phenotypes,
        "config": config,
        "audit_trail": audit_trail,
    }

    final_ctx = await engine.run(initial_context)

    if "output" not in final_ctx:
        return {"error": "Pipeline completed but no output was produced"}

    output = final_ctx["output"]
    workflow_meta = final_ctx.get("_workflow_meta", {})
    if "trace_path" in workflow_meta:
        output["trace_path"] = workflow_meta["trace_path"]
    if "audit_log_path" in workflow_meta:
        output["audit_log_path"] = workflow_meta["audit_log_path"]
    return output


# =============================================================================
# Multi-Organ Assessment  (v0.5 P1-7)
# =============================================================================

async def run_multi_organ_assessment(variants_data: List[Dict],
                                      user_phenotypes: Optional[str] = None,
                                      config: Optional[GPAConfig] = None) -> Dict:
    """
    v0.5 P1-7: Multi-organ joint assessment.

    Runs GPA for each profile in config.multi_organ_profiles, then generates
    a joint risk matrix taking the MAX tier across profiles per variant.

    API queries are performed per-profile (GTEx tissue-specific), but cached
    responses minimize redundant calls.

    Args:
        variants_data: List of variant dicts from VCF annotation
        config: GPA configuration with multi_organ_profiles set

    Returns:
        Dict with per-profile results + joint report
    """
    if config is None:
        config = GPAConfig()

    if not config.multi_organ_profiles or len(config.multi_organ_profiles) == 0:
        raise ValueError("multi_organ_profiles must be set for multi-organ assessment")

    profiles = config.multi_organ_profiles
    print(f"[GPA] Multi-organ assessment: {len(profiles)} profiles - {', '.join(profiles)}")

    # Run each profile independently
    profile_results = {}
    for profile_name in profiles:
        profile_config = GPAConfig(
            tissue_profile=profile_name,
            offline_mode=config.offline_mode,
            somatic_mode=config.somatic_mode,
            target_population=config.target_population,
            min_dp=config.min_dp,
            min_gq=config.min_gq,
            common_af_threshold=config.common_af_threshold,
            low_af_threshold=config.low_af_threshold,
            vaf_deviation_threshold=config.vaf_deviation_threshold,
            force_sync=config.force_sync,
        )
        print(f"\n[GPA] === Running profile: {profile_name} ===")
        result = await run_dgra_pipeline(variants_data, user_phenotypes=user_phenotypes, config=profile_config)
        profile_results[profile_name] = result

    # Build joint risk matrix: max tier across profiles per variant
    variant_tiers_by_profile = {}  # gene_pos -> {profile: tier}
    variant_details = {}  # gene_pos -> variant dict

    for profile_name, result in profile_results.items():
        for tier_field in ["tier1_variants", "tier2_variants", "tier3_variants"]:
            tier_num = int(tier_field.replace("tier", "").replace("_variants", ""))
            for v_dict in result.get(tier_field, []):
                key = (v_dict.get("gene", ""), v_dict.get("chrom", ""), v_dict.get("pos", 0))
                if key not in variant_tiers_by_profile:
                    variant_tiers_by_profile[key] = {}
                    variant_details[key] = v_dict
                variant_tiers_by_profile[key][profile_name] = tier_num

    # Compute max tier per variant
    joint_tiers = {}
    for key, tiers in variant_tiers_by_profile.items():
        max_tier = max(tiers.values())
        joint_tiers[key] = {
            "gene": key[0],
            "chrom": key[1],
            "pos": key[2],
            "max_tier": max_tier,
            "per_profile": tiers,
            "details": variant_details[key],
        }

    # Generate joint report
    joint_report = generate_multi_organ_report(profile_results, joint_tiers, profiles)

    # Summary
    tier1_joint = [k for k, v in joint_tiers.items() if v["max_tier"] == 1]
    tier2_joint = [k for k, v in joint_tiers.items() if v["max_tier"] == 2]
    tier3_joint = [k for k, v in joint_tiers.items() if v["max_tier"] == 3]

    # Cross-profile high-concern variants (Tier 1 in any profile)
    high_concern = []
    for key, info in joint_tiers.items():
        if info["max_tier"] == 1:
            profiles_t1 = [p for p, t in info["per_profile"].items() if t == 1]
            high_concern.append({
                "gene": info["gene"],
                "chrom": info["chrom"],
                "pos": info["pos"],
                "tier_1_in": profiles_t1,
                "all_profiles": info["per_profile"],
            })

    return {
        "meta": {
            "multi_organ": True,
            "profiles": profiles,
            "analysis_date": datetime.now().isoformat(),
            "total_variants": len(variants_data),
        },
        "profile_results": profile_results,
        "joint_summary": {
            "tier1_count": len(tier1_joint),
            "tier2_count": len(tier2_joint),
            "tier3_count": len(tier3_joint),
            "high_concern_variants": high_concern,
        },
        "joint_risk_matrix": joint_tiers,
        "joint_report_markdown": joint_report,
    }


def generate_multi_organ_report(profile_results: Dict[str, Dict],
                                 joint_tiers: Dict,
                                 profiles: List[str]) -> str:
    """
    v0.5 P1-7: Generate a joint multi-organ risk assessment report.

    Includes:
    - Joint summary (Tier 1/2/3 counts across all profiles)
    - Risk matrix table (variant x profile)
    - High-concern variants (Tier 1 in any profile)
    - Per-profile reports (appended)
    """
    report = []

    report.append("# GPA 多器官联合关联分析报告\n")
    report.append(f"**评估器官**: {', '.join(profiles)}\n")
    report.append(f"**分析日期**: {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    report.append(f"**联合策略**: 跨器官取最高 Tier(最保守)\n\n")

    # Joint summary
    tier1_count = sum(1 for k, v in joint_tiers.items() if v["max_tier"] == 1)
    tier2_count = sum(1 for k, v in joint_tiers.items() if v["max_tier"] == 2)
    tier3_count = sum(1 for k, v in joint_tiers.items() if v["max_tier"] == 3)

    report.append("## 联合风险摘要\n\n")
    report.append(f"- **Tier 1 (需干预)**: {tier1_count} 个变异\n")
    report.append(f"- **Tier 2 (需知情)**: {tier2_count} 个变异\n")
    report.append(f"- **Tier 3 (无需担忧)**: {tier3_count} 个变异\n\n")

    # Risk matrix
    report.append("## 联合风险矩阵\n\n")
    report.append("| 基因 | 位置 | " + " | ".join(profiles) + " | **最高 Tier** |\n")
    report.append("|------|------|" + "|".join(["------"] * len(profiles)) + "|----------|\n")

    for key, info in sorted(joint_tiers.items(), key=lambda x: (x[1]["gene"], x[1]["pos"])):
        gene = info["gene"]
        pos = f"{info['chrom']}:{info['pos']}"
        tier_cells = []
        for p in profiles:
            t = info["per_profile"].get(p, "-")
            tier_cells.append(str(t))
        max_tier = info["max_tier"]
        report.append(f"| {gene} | {pos} | " + " | ".join(tier_cells) + f" | **{max_tier}** |\n")

    report.append("\n")

    # High-concern variants
    high_concern = [(k, v) for k, v in joint_tiers.items() if v["max_tier"] == 1]
    if high_concern:
        report.append("## 高关注变异(任一器官 Tier 1)\n\n")
        for key, info in high_concern:
            gene = info["gene"]
            pos = f"{info['chrom']}:{info['pos']}"
            t1_profiles = [p for p, t in info["per_profile"].items() if t == 1]
            report.append(f"- **{gene}** ({pos}): Tier 1 于 {', '.join(t1_profiles)}\n")
            # Show all profile tiers
            tier_detail = ", ".join([f"{p}: Tier {t}" for p, t in info["per_profile"].items()])
            report.append(f"  - 全器官评估: {tier_detail}\n")
        report.append("\n")

    # Per-profile reports
    report.append("---\n\n")
    report.append("# 各器官独立评估报告\n\n")
    for profile_name in profiles:
        result = profile_results[profile_name]
        report.append(f"## [{profile_name}] {result['meta']['profile_display_name']}\n\n")
        report.append(result.get("report_markdown", "(无报告)\n"))
        report.append("\n---\n\n")

    return "\n".join(report)

# =============================================================================
# CLI Interface  (v0.4: --offline + asyncio.run)
# =============================================================================
