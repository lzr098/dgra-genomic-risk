#!/usr/bin/env python3
"""
GPA Two-Phase Pipeline (v0.10.1)

Optimization for large VCF datasets: splits analysis into
  Phase 1: Fast local triage using only VEP annotation + local gene lists
  Phase 2: API enrichment (gnomAD, SpliceAI, phenotype LLM) ONLY for Tier 1/2 candidates

This reduces API calls by 50-200x for typical germline VCFs where >95% of variants
are common SNPs or low-impact.

Usage:
    from gpa_two_phase import run_two_phase_pipeline
    results = await run_two_phase_pipeline(variants_data, config, user_phenotypes)
"""

import asyncio
import json
import time
from typing import List, Dict, Optional, Any, Tuple, Set

# Import from core (lazy to avoid circular)
from dgra_core import (
    Variant, GPAConfig, _UNKNOWN, classify_gnomad_frequency,
    normalize_gene_symbols, map_variant_to_domain, assess_tissue_relevance,
    predict_nmd, correct_transcript_priority, detect_pseudogene_artifact,
    _save_offline_archive, _load_offline_archive,
)
from dgra_cache import DGRACache
from dgra_api import DGRAAPIClient
from api_hub import APIHub, AdaptiveRateLimiter


# ─── Phase 1: Fast Local Triage ────────────────────────────────────────────

# Impact categories for pre-filtering
_HIGH_IMPACT_TERMS = {"HIGH"}
_MODERATE_IMPACT_TERMS = {"MODERATE"}

# Known neuromuscular disease genes (supplement to tissue_context.json lists)
# These are genes where MODERATE impact variants warrant Phase 2 enrichment
# even without prior ClinVar pathogenicity
_NEUROMUSCULAR_DISEASE_GENES: Set[str] = {
    # Muscular Dystrophies
    "DMD", "DYSF", "CAPN3", "SGCA", "SGCB", "SGCG", "SGCD",
    "FKRP", "FKTN", "POMT1", "POMT2", "POMGNT1", "LARGE1",
    "DAG1", "LMNA", "EMD", "SYNE1", "SYNE2", "TMEM43",
    "TTN", "MYOT", "FLNC", "BAG3", "DES", "CRYAB", "LDB3",
    "MYH7", "ACTA1", "TPM2", "TPM3", "TNNT1", "NEB",
    "RYR1", "CACNA1S", "MTM1", "DNM2", "BIN1", "AMPH",
    "SPEG", "CNTN1", "TRIM32", "TCAP", "ANO5", "PLEC",
    "COL6A1", "COL6A2", "COL6A3", "LAMA2", "ITGA7",
    "PABPN1", "GNE", "MYH2", "MATR3", "VCP", "HNRNPA1",
    "HNRNPA2B1", "SQSTM1", "TIA1",

    # Distal Myopathies
    "MYH7", "TIA1", "DNAJB6", "CRYAB",

    # Myotonic Disorders
    "CLCN1", "SCN4A", "DMPK", "CNBP",

    # Congenital Myopathies
    "ACTA1", "NEB", "TPM2", "TPM3", "TNNT1", "CFL2",
    "RYR1", "SEPN1", "MTM1", "DNM2", "BIN1", "TTN",
    "MYH7", "KLHL40", "KLHL41", "LMOD3",

    # Metabolic Myopathies
    "GAA", "AGL", "PYGM", "PFKM", "PGAM2", "LDHA",
    "CPT2", "ETFDH", "ETFA", "ETFB",

    # Motor Neuron / ALS
    "SOD1", "TARDBP", "FUS", "C9orf72", "UBQLN2",
    "OPTN", "TBK1", "CHCHD10", "KIF5A", "NEK1",

    # Charcot-Marie-Tooth / Neuropathy
    "PMP22", "MPZ", "GJB1", "MFN2", "GDAP1", "NEFL",
    "HSPB1", "HSPB8", "LITAF", "EGR2",

    # Mitochondrial
    "POLG", "TK2", "RRM2B", "SUCLA2", "SUCLG1",
    "OPA1", "MFN2",

    # Other Neuromuscular
    "SMN1", "SMN2", "IGHMBP2", "DYNC1H1", "BICD2",
    "CHRNA1", "CHRNB1", "CHRND", "CHRNE", "RAPSN",
    "DOK7", "MUSK", "AGRN", "COLQ",
}


def _parse_variants_phase1(variants_data: List[Dict]) -> List[Variant]:
    """Parse variant dicts into Variant objects (shared with main pipeline)."""
    # Chinese impact/consequence mapping (copied from gpa_pipeline.py)
    _IMPACT_CN_MAP = {"高": "HIGH", "中等": "MODERATE", "低": "LOW", "修饰": "MODIFIER"}
    _CONSEQUENCE_CN_MAP = {
        "错义变异": "missense_variant", "无义变异": "stop_gained",
        "获得终止密码子": "stop_gained", "移码变异": "frameshift_variant",
        "框内插入": "inframe_insertion", "剪接位点变异": "splice_donor_variant",
        "剪接区域变异": "splice_region_variant",
        "剪接供体区域变异": "splice_donor_variant",
        "剪接供体第5位碱基变异": "splice_donor_variant",
        "剪接多嘧啶束变异": "splice_polypyrimidine_tract_variant",
        "内含子变异": "intron_variant", "基因上游变异": "upstream_gene_variant",
        "基因下游变异": "downstream_gene_variant", "同义变异": "synonymous_variant",
        "非翻译区变异": "UTR_variant",
        "3'非翻译区变异": "3_prime_UTR_variant",
        "5'非翻译区变异": "5_prime_UTR_variant",
        "非编码转录本外显子变异": "non_coding_transcript_exon_variant",
    }

    variants = []
    for vd in variants_data:
        raw_impact = str(vd.get("IMPACT", "")).strip()
        if raw_impact in _IMPACT_CN_MAP:
            raw_impact = _IMPACT_CN_MAP[raw_impact]
        if not raw_impact:
            raw_impact = _UNKNOWN

        raw_consequence = str(vd.get("Consequence", "")).strip()
        if raw_consequence in _CONSEQUENCE_CN_MAP:
            raw_consequence = _CONSEQUENCE_CN_MAP[raw_consequence]
        if not raw_consequence:
            raw_consequence = _UNKNOWN

        raw_clinvar = str(vd.get("CLIN_SIG", "")).strip()
        if not raw_clinvar:
            raw_clinvar = _UNKNOWN

        raw_dp = str(vd.get("DP", "")).strip()
        try:
            dp_val = int(float(raw_dp)) if raw_dp and raw_dp != _UNKNOWN else 0
        except ValueError:
            dp_val = 0

        raw_gq = str(vd.get("GQ", "")).strip()
        try:
            gq_val = float(raw_gq) if raw_gq and raw_gq not in ("", _UNKNOWN, "None") else 0.0
        except (ValueError, TypeError):
            gq_val = 0.0

        raw_vaf = str(vd.get("VAF", "")).strip()
        if not raw_vaf or raw_vaf == ".":
            vaf_val = None
        else:
            try:
                vaf_val = float(raw_vaf)
            except (ValueError, TypeError):
                vaf_val = None

        raw_gnomad = str(vd.get("gnomAD_AF", "")).strip()
        gnomad_val = float(raw_gnomad) if raw_gnomad and raw_gnomad != "N/A" else None

        gene = vd.get("GENE", "") or vd.get("Gene", "")

        v = Variant(
            chrom=vd.get("chrom", vd.get("CHROM", "")),
            pos=int(vd.get("pos", vd.get("POS", 0)) or 0),
            ref=vd.get("ref", vd.get("REF", "")),
            alt=vd.get("alt", vd.get("ALT", "")),
            gene=gene,
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
        )
        variants.append(v)
    return variants


# v0.10.2: Consequence terms that are LOW impact in VEP but potentially pathogenic
_SPLICE_RELEVANT_CONSEQUENCES = {
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_acceptor_variant",
    "incomplete_terminal_codon_variant",
}


def _is_potentially_pathogenic(v: Variant) -> bool:
    """Phase 1 fast check: does this variant warrant deeper investigation?

    v0.10.2 enhancements:
    - Consequence-aware: splice_region variants (LOW impact) now pass through
    - VEP gnomAD AF pre-filter: common missense (EAS AF > 1%) excluded
    """
    # ── EAS AF pre-filter (applies to MODERATE impact missense) ──
    # Use VEP-provided gnomAD AF if available (already in v.gnomad_af)
    # This avoids passing large numbers of common polymorphisms to Phase 2
    if v.impact == "MODERATE" and v.gene not in _NEUROMUSCULAR_DISEASE_GENES:
        if v.gnomad_af is not None and v.gnomad_af > 0.01:
            # Common variant (AF > 1%) in a non-disease gene → skip
            return False

    # 1. HIGH impact (stop_gained, frameshift, splice donor/acceptor)
    if v.impact == "HIGH":
        # Skip common homozygous (AF=1.00 ref mismatches)
        if v.gnomad_af and v.gnomad_af > 0.95 and v.gt == "1/1":
            return False
        return True

    # 2. MODERATE impact
    if v.impact == "MODERATE":
        # 2a. ClinVar Pathogenic/Likely_pathogenic
        clinvar_lower = (v.clinvar or "").lower()
        if any(kw in clinvar_lower for kw in ("pathogenic", "致病", "likely_pathogenic", "可能致病")):
            return True
        # 2b. Known neuromuscular disease gene
        if v.gene in _NEUROMUSCULAR_DISEASE_GENES:
            # For disease genes, still exclude if AF > 5% (clearly benign)
            if v.gnomad_af is not None and v.gnomad_af > 0.05:
                return False
            return True
        # 2c. Very rare (AF < 0.001 if known, or unknown)
        if v.gnomad_af is not None and v.gnomad_af < 0.001:
            return True
        if v.gnomad_af is None:
            # Unknown frequency → conservative: include for enrichment
            return True

    # 3. LOW/MODIFIER impact: consequence-aware filtering (v0.10.2)
    if v.impact in ("LOW", "MODIFIER"):
        # 3a. Splice-relevant consequences → pass through (may disrupt splicing)
        if v.consequence in _SPLICE_RELEVANT_CONSEQUENCES:
            return True
        # 3b. ClinVar Pathogenic (rare but possible for regulatory variants)
        if v.clinvar not in (_UNKNOWN, ""):
            clinvar_lower = (v.clinvar or "").lower()
            if any(kw in clinvar_lower for kw in ("pathogenic", "致病")):
                return True

    return False


def _fast_tier_classification(v: Variant) -> Tuple[int, str]:
    """
    Phase 1: Fast local-only tier assignment.
    Conservative: errs on the side of higher tier.
    Phase 2 will refine with API data.
    """
    if not _is_potentially_pathogenic(v):
        return 3, "Pre-filtered: low impact, no pathogenic evidence, not in disease gene list"

    # Candidate — assign preliminary tier
    if v.impact == "HIGH":
        return 1, "PRELIMINARY: HIGH impact variant (Phase 2 enrichment pending)"
    elif v.clinvar and ("pathogenic" in (v.clinvar or "").lower() or "致病" in (v.clinvar or "")):
        return 1, "PRELIMINARY: ClinVar Pathogenic (Phase 2 enrichment pending)"
    elif v.consequence in _SPLICE_RELEVANT_CONSEQUENCES:
        return 1, "PRELIMINARY: Splice-relevant consequence (Phase 2 SpliceAI pending)"
    else:
        return 2, "PRELIMINARY: Candidate variant (Phase 2 enrichment pending)"


def _candidate_priority_score(v: Variant) -> tuple:
    """
    v0.10.4: Priority score for candidate ranking when max_candidates is enforced.
    Lower tuple = higher priority. Used to truncate excess candidates.
    """
    # 1. Tier priority (1 > 2)
    tier_score = 0 if getattr(v, 'tier', 2) == 1 else 1

    # 2. Impact priority
    impact_order = {"HIGH": 0, "MODERATE": 1, "LOW": 2, "MODIFIER": 3}
    impact_score = impact_order.get(v.impact, 4)

    # 3. ClinVar priority
    clinvar_lower = (v.clinvar or "").lower()
    has_pathogenic = any(kw in clinvar_lower for kw in ("pathogenic", "致病", "likely_pathogenic", "可能致病"))
    clinvar_score = 0 if has_pathogenic else 1

    # 4. Frequency priority: unknown (None) is most interesting
    if v.gnomad_af is None:
        af_score = 0
    elif v.gnomad_af < 0.001:
        af_score = 1
    elif v.gnomad_af < 0.01:
        af_score = 2
    else:
        af_score = 3

    return (tier_score, impact_score, clinvar_score, af_score)


# ─── Phase 2: Deep Enrichment for Candidates ───────────────────────────────

async def _enrich_candidate_genes(
    variants: List[Variant],
    candidate_indices: List[int],
    config: GPAConfig,
    tissue_profile: Dict,
    tracker: Optional[Any] = None,
) -> Tuple[Dict, Dict, Dict]:
    """
    Phase 2: API enrichment ONLY for candidate genes.
    Returns (ensembl_data, uniprot_data, gnomad_constraint_data).
    """
    if not candidate_indices:
        return {}, {}, {}

    candidate_genes = list({variants[i].gene for i in candidate_indices})
    n_genes = len(candidate_genes)
    global_config = config.to_global()

    if config.offline_mode:
        ensembl_data, uniprot_data, gnomad_constraint_data = {}, {}, {}
        for idx, gene in enumerate(candidate_genes):
            archive = _load_offline_archive(gene)
            if archive:
                ensembl_data[gene] = archive.get("ensembl", {})
                uniprot_data[gene] = archive.get("uniprot", {})
                gc = archive.get("gnomad_constraint")
                if gc and gc.get("status") == "CAPTURED":
                    gnomad_constraint_data[gene] = gc
            if tracker and (idx + 1) % max(1, n_genes // 10) == 0:
                tracker.step_progress("phase2", "gene_apis", idx + 1, n_genes)
        return ensembl_data, uniprot_data, gnomad_constraint_data

    cache = DGRACache(global_config.cache_db_path)
    async with DGRAAPIClient(global_config, cache) as client:
        if tracker:
            tracker.api_call("phase2", "gene_apis", "ensembl", "started", {"genes": n_genes})
            tracker.api_call("phase2", "gene_apis", "uniprot", "started", {"genes": n_genes})
        # v0.12.1 FIX: gnomAD GraphQL is 429-rate-limited in batch mode;
        # skip gnomad_constraint here — gnomAD AF is provided by MyVariant.info.
        ensembl_raw, uniprot_raw = await asyncio.gather(
            client.batch_query_genes(candidate_genes, "ensembl"),
            client.batch_query_genes(candidate_genes, "uniprot"),
        )
        if tracker:
            tracker.api_call("phase2", "gene_apis", "ensembl", "completed",
                             {"filled": len([g for g in candidate_genes if ensembl_raw.get(g)])})
            tracker.api_call("phase2", "gene_apis", "uniprot", "completed",
                             {"filled": len([g for g in candidate_genes if uniprot_raw.get(g)])})
        ensembl_data = {g: ensembl_raw.get(g, {}) for g in candidate_genes}
        uniprot_data = {g: uniprot_raw.get(g, {}) for g in candidate_genes}
        gnomad_constraint_data = {}  # skipped — gnomAD GraphQL is rate-limited

    return ensembl_data, uniprot_data, gnomad_constraint_data


async def _enrich_gtex(
    variants: List[Variant],
    candidate_indices: List[int],
    config: GPAConfig,
    tissue_profile: Dict,
    user_phenotypes: Optional[str] = None,
    tracker: Optional[Any] = None,
) -> Dict[str, Dict]:
    """
    Phase 2: GTEx expression data for candidate genes (v0.10.8).
    
    v0.10.8 ARCHITECTURE CHANGE:
    - Queries ALL 54 GTEx tissues (not limited to profile tissues)
    - Performs phenotype-tissue association analysis
    - GTEx data is used for phenotype-relevance scoring, NOT as a hard tier gate
    
    Returns {gene: {
        "median_tpm": float, "max_tpm": float, "all_tissues": [...],
        "source": "gtex_multi", "expressing_tissues": int,
        "phenotype_max_tpm": float, "phenotype_tissues": [...],
        "global_max_tpm": float,
    }}
    """
    if not candidate_indices:
        return {}
    
    # v0.10.8: Query ALL GTEx tissues, not just profile-specific ones
    ALL_GTEX_TISSUES = [
        "Adipose - Subcutaneous", "Adipose - Visceral (Omentum)", "Adrenal Gland",
        "Artery - Aorta", "Artery - Coronary", "Artery - Tibial", "Bladder",
        "Brain - Amygdala", "Brain - Anterior cingulate cortex (BA24)", "Brain - Caudate (basal ganglia)",
        "Brain - Cerebellar Hemisphere", "Brain - Cerebellum", "Brain - Cortex",
        "Brain - Frontal Cortex (BA9)", "Brain - Hippocampus", "Brain - Hypothalamus",
        "Brain - Nucleus accumbens (basal ganglia)", "Brain - Putamen (basal ganglia)",
        "Brain - Spinal cord (cervical c-1)", "Brain - Substantia nigra",
        "Breast - Mammary Tissue", "Cells - EBV-transformed lymphocytes",
        "Cells - Cultured fibroblasts", "Cervix - Ectocervix", "Cervix - Endocervix",
        "Colon - Sigmoid", "Colon - Transverse", "Esophagus - Gastroesophageal Junction",
        "Esophagus - Mucosa", "Esophagus - Muscularis", "Fallopian Tube",
        "Heart - Atrial Appendage", "Heart - Left Ventricle", "Kidney - Cortex",
        "Kidney - Medulla", "Liver", "Lung", "Minor Salivary Gland",
        "Muscle - Skeletal", "Nerve - Tibial", "Ovary", "Pancreas", "Pituitary",
        "Prostate", "Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)",
        "Small Intestine - Terminal Ileum", "Spleen", "Stomach", "Testis",
        "Thyroid", "Uterus", "Vagina", "Whole Blood",
    ]
    
    # Phenotype-to-tissue keyword mapping (Chinese + English)
    PHENOTYPE_TISSUE_MAP = {
        # Eye / Vision (GTEx has no retina; use brain as proxy for vision-related genes)
        "retina": ["Brain - Cortex", "Brain - Cerebellum", "Brain - Hippocampus"],
        "eye": ["Brain - Cortex", "Brain - Cerebellum", "Brain - Hippocampus"],
        "vision": ["Brain - Cortex", "Brain - Cerebellum"],
        "optic": ["Brain - Cortex", "Brain - Cerebellum"],
        " ophthalm": ["Brain - Cortex", "Brain - Cerebellum"],
        "视网膜": ["Brain - Cortex", "Brain - Cerebellum", "Brain - Hippocampus"],
        "眼": ["Brain - Cortex", "Brain - Cerebellum"],
        "视觉": ["Brain - Cortex", "Brain - Cerebellum"],
        "牵牛花": ["Brain - Cortex", "Brain - Cerebellum"],
        # Brain / Neurological
        "brain": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "neuro": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "neural": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "cognitive": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "autism": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "脑": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "神经": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "智力": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "自闭症": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        "认知": [t for t in ALL_GTEX_TISSUES if t.startswith("Brain -")],
        # Heart
        "heart": ["Heart - Atrial Appendage", "Heart - Left Ventricle", "Artery - Aorta", "Artery - Coronary"],
        "cardiac": ["Heart - Atrial Appendage", "Heart - Left Ventricle", "Artery - Aorta", "Artery - Coronary"],
        "cardio": ["Heart - Atrial Appendage", "Heart - Left Ventricle", "Artery - Aorta", "Artery - Coronary"],
        "心": ["Heart - Atrial Appendage", "Heart - Left Ventricle"],
        # Muscle
        "muscle": ["Muscle - Skeletal", "Heart - Left Ventricle"],
        "myopathy": ["Muscle - Skeletal", "Heart - Left Ventricle"],
        "肌肉": ["Muscle - Skeletal"],
        # Kidney
        "kidney": ["Kidney - Cortex", "Kidney - Medulla"],
        "renal": ["Kidney - Cortex", "Kidney - Medulla"],
        "肾": ["Kidney - Cortex", "Kidney - Medulla"],
        # Liver
        "liver": ["Liver"],
        "hepatic": ["Liver"],
        "肝": ["Liver"],
        # Lung
        "lung": ["Lung"],
        "pulmonary": ["Lung"],
        "肺": ["Lung"],
        # Blood
        "blood": ["Whole Blood", "Cells - EBV-transformed lymphocytes"],
        "hemato": ["Whole Blood", "Cells - EBV-transformed lymphocytes", "Spleen"],
        "leukemia": ["Whole Blood", "Cells - EBV-transformed lymphocytes", "Spleen"],
        "血液": ["Whole Blood", "Cells - EBV-transformed lymphocytes"],
        # Skin
        "skin": ["Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)"],
        "皮": ["Skin - Not Sun Exposed (Suprapubic)", "Skin - Sun Exposed (Lower leg)"],
        # Adrenal
        "adrenal": ["Adrenal Gland"],
        "肾上腺": ["Adrenal Gland"],
        # Thyroid
        "thyroid": ["Thyroid"],
        "甲状腺": ["Thyroid"],
        # Nerve
        "nerve": ["Nerve - Tibial", "Brain - Spinal cord (cervical c-1)"],
        "神经": ["Nerve - Tibial", "Brain - Spinal cord (cervical c-1)"],
    }
    
    # Determine phenotype-relevant tissues from user input
    phenotype_tissues: set = set()
    phenotype_matched_keywords: list = []
    if user_phenotypes:
        pheno_lower = user_phenotypes.lower()
        for keyword, tissues in PHENOTYPE_TISSUE_MAP.items():
            if keyword.lower() in pheno_lower:
                phenotype_tissues.update(tissues)
                phenotype_matched_keywords.append(keyword)
    
    # Fallback: if no phenotype match, use profile tissues
    if not phenotype_tissues:
        profile_tissues = tissue_profile.get("gtex_tissues", [])
        if not profile_tissues:
            single = tissue_profile.get("gtex_tissue")
            if single:
                profile_tissues = [single]
        phenotype_tissues.update(profile_tissues)
    
    phenotype_tissues = sorted(phenotype_tissues)
    if phenotype_matched_keywords:
        print(f"[GPA Phase 2] GTEx phenotype-tissue: matched keywords={phenotype_matched_keywords}")
        print(f"[GPA Phase 2] GTEx phenotype-relevant tissues: {len(phenotype_tissues)} tissues")
    
    candidate_genes = sorted({variants[i].gene for i in candidate_indices})
    gtex_data: Dict[str, Dict] = {}

    # v0.10.16: Query local GTEx DB directly (zero external API calls)
    try:
        import sys
        from pathlib import Path
        gtex_local_path = str(Path.home() / ".workbuddy" / "scripts")
        if gtex_local_path not in sys.path:
            sys.path.insert(0, gtex_local_path)
        from gtex_local import query_gtex_local

        n_genes = len(candidate_genes)
        progress_interval = max(1, n_genes // 10)
        for idx, gene in enumerate(candidate_genes):
            local_expr = query_gtex_local(gene, tissues=None)
            if not local_expr:
                continue

            tissue_tpm = dict(local_expr)
            global_max = max(tissue_tpm.values()) if tissue_tpm else 0.0

            # Phenotype-relevant max
            phenotype_values = [tissue_tpm.get(t) for t in phenotype_tissues if t in tissue_tpm]
            phenotype_max = max(phenotype_values) if phenotype_values else 0.0

            all_expressing = sorted(
                [(t, tpm) for t, tpm in tissue_tpm.items() if tpm > 0],
                key=lambda x: x[1], reverse=True
            )
            expressing_count = len(all_expressing)

            gtex_data[gene] = {
                "median_tpm": global_max,
                "max_tpm": global_max,
                "all_tissues": [{"tissue": t, "tpm": tpm} for t, tpm in all_expressing],
                "expressing_tissues": expressing_count,
                "source": "gtex_local_db",
                "confidence": "medium",
                "phenotype_max_tpm": phenotype_max,
                "phenotype_tissues": phenotype_tissues,
                "phenotype_matched_keywords": phenotype_matched_keywords,
                "global_max_tpm": global_max,
            }
            if tracker and (idx + 1) % progress_interval == 0:
                tracker.step_progress("phase2", "gtex", idx + 1, n_genes)
        print(f"[GPA Phase 2] GTEx local DB: expression data for {len(gtex_data)}/{len(candidate_genes)} candidate genes")
    except Exception as e:
        print(f"[GPA Phase 2] GTEx local DB query failed: {type(e).__name__}: {e}")
        if tracker:
            tracker.error("phase2", "gtex", f"GTEx local DB query failed: {e}")

    return gtex_data


# ─── v0.10.3: Consequence-aware API enrichment rules ───────────────────────

# ClinVar only needed for variants where functional impact is uncertain
_CLINVAR_RELEVANT_CONSEQUENCES = {
    "missense_variant",
    "splice_region_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "inframe_deletion",
    "inframe_insertion",
    "frameshift_variant",
    "incomplete_terminal_codon_variant",
}

# SpliceAI only needed for splice-site-proximal variants
_SPLICEAI_RELEVANT_CONSEQUENCES = {
    "splice_acceptor_variant",
    "splice_donor_variant",
    "splice_donor_5th_base_variant",
    "splice_donor_region_variant",
    "splice_polypyrimidine_tract_variant",
    "splice_region_variant",
}

def _variant_needs_clinvar(v: Variant) -> bool:
    """Check if variant's consequence warrants ClinVar lookup."""
    cons = v.consequence or ""
    # Check if any relevant consequence term is present (VEP consequences are comma-separated)
    cons_terms = set(c.strip() for c in cons.split(","))
    return bool(cons_terms & _CLINVAR_RELEVANT_CONSEQUENCES)

def _variant_needs_spliceai(v: Variant) -> bool:
    """Check if variant's consequence warrants SpliceAI lookup."""
    cons = v.consequence or ""
    cons_terms = set(c.strip() for c in cons.split(","))
    # Direct splice consequence match
    if bool(cons_terms & _SPLICEAI_RELEVANT_CONSEQUENCES):
        return True
    # intron_variant within 50bp of splice site is also relevant
    # (we check distance if exon boundary info available in vcf_info)
    if "intron_variant" in cons_terms:
        vcf_info = getattr(v, 'vcf_info', None) or {}
        intron_dist = vcf_info.get('SPLICE_DIST')
        if intron_dist is not None and abs(int(intron_dist)) <= 50:
            return True
    return False


async def _enrich_variant_frequencies(
    variants: List[Variant],
    candidate_indices: List[int],
    config: GPAConfig,
    tracker: Optional[Any] = None,
) -> None:
    """
    Phase 2: gnomAD frequency lookup ONLY for Tier 1/2 candidate variants.

    v0.10.13 OPTIMIZATION: Skip gnomAD queries for Tier 3 variants entirely.
    Phase 1 _fast_tier_classification() already assigns preliminary tiers.
    We only enrich variants with preliminary tier 1 or 2, since Tier 3 variants
    (common SNPs, low impact, high AF) do not need gnomAD frequency data for
    final classification.
    """
    if not candidate_indices:
        return

    import aiohttp
    global_config = config.to_global()

    # v0.10.13: Filter to Tier 1/2 candidates only before querying gnomAD
    tier12_candidates_no_af = [
        variants[i] for i in candidate_indices
        if variants[i].gnomad_af is None and getattr(variants[i], 'tier', 3) in (1, 2)
    ]

    # Count how many were skipped due to Tier 3 status
    tier3_skipped = [
        variants[i] for i in candidate_indices
        if variants[i].gnomad_af is None and getattr(variants[i], 'tier', 3) == 3
    ]
    if tier3_skipped:
        print(f"[GPA Phase 2] gnomAD: skipping {len(tier3_skipped)} Tier 3 variants (no AF enrichment needed)")

    if not tier12_candidates_no_af:
        print(f"[GPA Phase 2] All Tier 1/2 candidates already have AF data — skipping gnomAD queries")
        return

    print(f"[GPA Phase 2] gnomAD: querying {len(tier12_candidates_no_af)} Tier 1/2 candidate variants for AF")
    if tracker:
        tracker.api_call("phase2", "gnomad_freq", "myvariant", "started",
                         {"variants": len(tier12_candidates_no_af), "batch_size": 100})

    # v0.12.2 FIX (E2): unify aiohttp proxy behavior with DGRAAPIClient.
    # MyVariant session must use the same proxy as DGRAAPIClient to avoid
    # inconsistent network behavior (some APIs via proxy, others direct).
    mv_proxy = None
    if getattr(global_config, 'proxy', None) and global_config.proxy != "__DIRECT__":
        mv_proxy = global_config.proxy
    elif getattr(global_config, '_proxy_route_map', None) is not None:
        try:
            mv_proxy = global_config._proxy_route_map.get_proxy("myvariant")
        except Exception:
            pass

    # v0.12.2: Respect skip_gnomad config (replaces phase2_analysis.py monkey-patch)
    if getattr(config, 'skip_gnomad', False):
        print("[GPA Phase 2] skip_gnomad=True — skipping MyVariant/gnomAD frequency queries")
    else:
        # Try MyVariant.info batch first — query Tier 1/2 candidates for gnomAD AF + ClinVar + CADD.
        # v0.13.0 FIX: batch_size 1000 -> 100, semaphore 10 -> 5, add inter-batch delay.
        try:
            from dgra_myvariant import query_myvariant_batch, apply_myvariant_results
            mv_sem = asyncio.Semaphore(5)
            timeout_obj = aiohttp.ClientTimeout(total=120)
            mv_variants = [(v.chrom, v.pos, v.ref, v.alt) for v in tier12_candidates_no_af]
            n_mv = len(mv_variants)
            print(f"[GPA Phase 2] MyVariant.info: {n_mv} variants, batch=100, sem=5")
            # ponytail: use APIHub session instead of ad-hoc ClientSession
            async with APIHub(global_config, None, detect_proxy=False) as hub:
                mv_results = await query_myvariant_batch(
                    mv_variants, hub.session, semaphore=mv_sem, batch_size=100,
                    proxy=mv_proxy,
                )
            mv_stats = apply_myvariant_results(tier12_candidates_no_af, mv_results)
            print(f"[GPA Phase 2] MyVariant.info: {mv_stats['gnomad_filled']} gnomAD, "
                  f"{mv_stats['clinvar_filled']} ClinVar, "
                  f"{mv_stats['cadd_filled']} CADD filled, "
                  f"{mv_stats['not_found']} not_found, "
                  f"{mv_stats['errors']} errors")
            if tracker:
                tracker.api_call("phase2", "gnomad_freq", "myvariant", "completed", mv_stats)
            # v0.12.2 FIX (B3): diagnostic warning when MyVariant returns all empty
            if mv_stats['gnomad_filled'] == 0 and mv_stats['errors'] > 0:
                warn_msg = (f"MyVariant.info returned zero gnomAD AF for all "
                            f"{len(tier12_candidates_no_af)} variants ({mv_stats['errors']} errors)")
                print(f"[GPA Phase 2] WARNING: {warn_msg}")
                if tracker:
                    tracker.warning("phase2", "gnomad_freq", warn_msg, mv_stats)
        except Exception as e:
            print(f"[GPA Phase 2] MyVariant.info batch query failed (non-critical): {type(e).__name__}: {e}")
            if tracker:
                tracker.error("phase2", "gnomad_freq", f"MyVariant.info failed: {e}")

        # v0.12.2 FIX (B1): restore gnomAD REST fallback when MyVariant fails.
        # v0.13.0: Replace fixed Semaphore(2) with AdaptiveRateLimiter for gnomAD GraphQL.
        still_no_af = [v for v in tier12_candidates_no_af if v.gnomad_af is None]
        if still_no_af:
            print(f"[GPA Phase 2] gnomAD REST fallback: querying {len(still_no_af)} Tier 1/2 candidates without AF")
            if tracker:
                tracker.api_call("phase2", "gnomad_freq", "gnomad_rest_fallback", "started",
                                 {"variants": len(still_no_af)})
            cache = DGRACache(global_config.cache_db_path)
            async with DGRAAPIClient(global_config, cache) as client:
                # v0.13.0: Adaptive rate limiter — starts at 2 req/s,
                # halves on 429, boosts +20% after 5 consecutive successes.
                gnomad_limiter = AdaptiveRateLimiter(
                    initial_rate=2.0, min_rate=0.2, max_rate=4.0,
                    success_threshold=5, rate_boost=1.2, rate_cut=0.5,
                )

                progress_counter = {"done": 0}
                n_fallback = len(still_no_af)
                progress_interval = max(1, n_fallback // 10)

                async def _query_one(v):
                    await gnomad_limiter.acquire()
                    try:
                        result = await client.query_gnomad_variant(v.chrom, v.pos, v.ref, v.alt)
                        # Check if result indicates 429/rate-limit at GraphQL level
                        if result.get("status") == "QUERY_ERROR" and "rate" in str(result.get("note", "")).lower():
                            gnomad_limiter.report_429()
                        else:
                            gnomad_limiter.report_success()
                        return result
                    except Exception as e:
                        err_str = str(e).lower()
                        if "429" in err_str or "rate" in err_str:
                            gnomad_limiter.report_429()
                        else:
                            gnomad_limiter.report_error(is_429=False)
                        return {"status": "API_FAILED", "error": str(e), "source": "failed"}
                    finally:
                        if tracker:
                            progress_counter["done"] += 1
                            if progress_counter["done"] % progress_interval == 0:
                                tracker.step_progress("phase2", "gnomad_freq",
                                                      progress_counter["done"], n_fallback)

                results = await asyncio.gather(*[_query_one(v) for v in still_no_af])
                n_filled = 0
                for v, result in zip(still_no_af, results):
                    if result and result.get("source") in ("gnomad", "cache", "failed"):
                        af = result.get("af")
                        if af is not None:
                            v.gnomad_af = af
                            v.gnomad_populations = result.get("af_populations", {})
                            v.gnomad_status = "SUCCESS"
                            n_filled += 1
                stats = gnomad_limiter.stats
                print(f"[GPA Phase 2] gnomAD REST fallback: {n_filled}/{len(still_no_af)} filled, "
                      f"rate={stats['current_rate']} req/s, 429s={stats['429_count']}")
                if tracker:
                    tracker.api_call("phase2", "gnomad_freq", "gnomad_rest_fallback", "completed",
                                     {"filled": n_filled, "total": len(still_no_af), **stats})

        # Final warning if AF is still missing for many candidates
        still_no_af_final = [v for v in tier12_candidates_no_af if v.gnomad_af is None]
        if still_no_af_final:
            print(f"[GPA Phase 2] WARNING: {len(still_no_af_final)}/{len(tier12_candidates_no_af)} Tier 1/2 variants "
                  f"still have no gnomAD AF data. Frequency-based filtering (EAS AF>1%) will not apply. "
                  f"Candidate set may be inflated with common polymorphisms.")

    # v0.10.6 FIX: NCBI ClinVar direct query for variants still UNKNOWN after MyVariant.info.
    # v0.12.2: Respect skip_clinvar config (replaces phase2_analysis.py monkey-patch)
    if getattr(config, 'skip_clinvar', False):
        print("[GPA Phase 2] skip_clinvar=True — skipping NCBI ClinVar queries")
    else:
        clinvar_unknown = [variants[i] for i in candidate_indices
                          if variants[i].clinvar in (_UNKNOWN, "UNKNOWN", "") and getattr(variants[i], 'tier', 3) in (1, 2)]
        if clinvar_unknown:
            # Consequence-aware filtering: skip intron/UTR where ClinVar is rarely informative
            cv_candidates = [v for v in clinvar_unknown if _variant_needs_clinvar(v)]
            if cv_candidates:
                n_cv = len(cv_candidates)
                print(f"[GPA Phase 2] ClinVar (NCBI): querying {n_cv} variants (serial, 3 req/s)")
                cache = DGRACache(global_config.cache_db_path)
                async with DGRAAPIClient(global_config, cache) as client:
                    # v0.12.2 FIX: NCBI E-utilities true serial query.
                    # asyncio.gather + Semaphore(1) does NOT serialize because
                    # coroutines are created before acquire; HTTP connections
                    # may be established concurrently. Use for-loop + sleep.
                    cv_results: List[Dict[str, Any]] = []
                    for i, v in enumerate(cv_candidates):
                        try:
                            result = await client.query_ncbi_clinvar(
                                gene=v.gene, chrom=v.chrom, pos=v.pos
                            )
                            cv_results.append(result)
                        except Exception as e:
                            cv_results.append({"clinical_significance": None, "error": str(e), "source": "failed"})
                        # 0.34s delay = ~3 req/s (NCBI safe limit without API key)
                        if i < n_cv - 1:
                            await asyncio.sleep(0.34)
                    n_found = 0
                    for v, cv_result in zip(cv_candidates, cv_results):
                        sig = cv_result.get("clinical_significance")
                        if sig:
                            v.clinvar = sig
                            v.clinvar_review_status = cv_result.get("review_status")
                            n_found += 1
                    print(f"[GPA Phase 2] ClinVar (NCBI): {n_found}/{n_cv} filled")


def _sa_get(sa, key, default=None):
    """Uniform accessor for spliceai_result (supports dict and SpliceAIResult dataclass)."""
    if isinstance(sa, dict):
        return sa.get(key, default)
    return getattr(sa, key, default)


async def _enrich_spliceai(
    variants: List[Variant],
    candidate_indices: List[int],
    config: GPAConfig,
    tracker: Optional[Any] = None,
) -> None:
    """Phase 2: SpliceAI for ALL Tier 1/2 candidate variants — AUTO-ENABLED in v0.10.7.

    v0.10.7 CHANGE: Removed manual spliceai_enabled gate. SpliceAI is now
    automatically queried for any candidate variant with splice-relevant
    consequences. This ensures HIGH-impact splice variants (e.g. SAG c.72_75+15del)
    are always assessed without requiring user to pass --spliceai.

    v0.10.13 OPTIMIZATION: SpliceAI is now queried for ALL Tier 1/2 variants,
    not just those with splice-relevant consequences. Variants without splice
    changes will return delta=0, which is valid evidence (confirms no splicing
    impact). This provides comprehensive splicing evidence for all candidate
    variants and enables splice-aware classification for missense and other
    non-canonical variants that may still affect splicing.
    """
    if not candidate_indices:
        return

    from dgra_splice_predictor import (
        query_spliceai_batch, _cache_key as _splice_key
    )

    # v0.10.13: Query SpliceAI for ALL Tier 1/2 candidates, not just splice-relevant ones.
    # The SpliceAI API handles non-splice variants gracefully (returns delta=0 or not_in_db).
    # This ensures we have complete splicing evidence for every candidate variant.
    tier12_candidates = [
        variants[i] for i in candidate_indices
        if getattr(variants[i], 'tier', 3) in (1, 2)
    ]

    if not tier12_candidates:
        print(f"[GPA Phase 2] SpliceAI: no Tier 1/2 variants — skipping")
        return

    print(f"[GPA Phase 2] SpliceAI: querying {len(tier12_candidates)} Tier 1/2 candidate variants "
          f"(concurrency={getattr(config, 'spliceai_concurrency', 5)})")
    if tracker:
        tracker.api_call("phase2", "spliceai", "spliceai_batch", "started",
                         {"variants": len(tier12_candidates)})

    spliceai_sem = asyncio.Semaphore(getattr(config, 'spliceai_concurrency', 5))
    spliceai_results = await query_spliceai_batch(
        tier12_candidates, spliceai_sem,
        timeout=getattr(config, 'spliceai_timeout', 45),
    )

    # Attach SpliceAI results to all Tier 1/2 candidates.
    # For variants without splice changes, the API returns delta=0 or not_in_db,
    # both of which are valid evidence entries.
    n_spliceai_found = 0
    for v in tier12_candidates:
        key = _splice_key(v.chrom, v.pos, v.ref, v.alt)
        if key in spliceai_results:
            v.spliceai_result = spliceai_results[key]
            if _sa_get(spliceai_results[key], "delta_score") is not None:
                n_spliceai_found += 1
        else:
            # Not queried or not in results — mark as not_in_db with null scores.
            # delta_score=None is used by the tier classifier to distinguish
            # "not queried" from "queried and delta=0".
            v.spliceai_result = {"source": "not_in_db", "delta_score": None, "predicted_impact": None}

    print(f"[GPA Phase 2] SpliceAI: batch complete for {len(tier12_candidates)} variants")
    if tracker:
        tracker.api_call("phase2", "spliceai", "spliceai_batch", "completed",
                         {"variants": len(tier12_candidates), "with_scores": n_spliceai_found})


async def _enrich_phenotype(
    variants: List[Variant],
    candidate_indices: List[int],
    user_phenotypes: Optional[str],
    tracker: Optional[Any] = None,
    config: Optional[Any] = None,
) -> None:
    """Phase 2: Phenotype LLM matching ONLY for candidate genes with local DB data."""
    if not user_phenotypes or not candidate_indices:
        return

    candidate_genes = sorted({variants[i].gene for i in candidate_indices})

    # v0.10.3: Pre-filter — only query genes that have known phenotype data locally
    # v0.10.17: Pass offline_mode so matcher can also use OMIM-backed keyword matching.
    offline_mode = bool(getattr(config, "offline_mode", False))
    try:
        from gpa_phenotype_match import PhenotypeMatcher
        matcher = PhenotypeMatcher(offline_mode=offline_mode)
        # Check which genes have local phenotype entries by reading _local_db directly
        genes_with_data = []
        genes_without_data = []
        for gene in candidate_genes:
            known = []
            if matcher._local_db and gene in matcher._local_db:
                entries = matcher._local_db[gene].get("phenotypes", [])
                known = [e.get("name", "") for e in entries if e.get("name")]
            if known:
                genes_with_data.append(gene)
            else:
                genes_without_data.append(gene)
        if genes_without_data:
            print(f"[GPA Phase 2] Phenotype: {len(genes_without_data)} genes skipped (no local phenotype data): {', '.join(genes_without_data[:5])}{'...' if len(genes_without_data) > 5 else ''}")
        if not genes_with_data:
            print("[GPA Phase 2] Phenotype: no candidate genes with known phenotype data — skipping LLM calls")
            if tracker:
                tracker.api_call("phase2", "phenotype", "phenotype_match", "skipped",
                                 {"reason": "no_local_data", "genes_checked": len(candidate_genes)})
            for i in candidate_indices:
                v = variants[i]
                v.phenotype_match_score = 0.0
                v.phenotype_match_explanation = "No known phenotypes found for this gene in local database."
                v.phenotype_match_confidence = "low"
            return
    except Exception as e:
        print(f"[GPA Phase 2] Phenotype pre-filter error (proceeding with all): {e}")
        if tracker:
            tracker.warning("phase2", "phenotype", f"Phenotype pre-filter error: {e}")
        genes_with_data = candidate_genes
        genes_without_data = []

    print(f"[GPA Phase 2] Phenotype matching: {len(genes_with_data)} candidate genes (with local data)")

    try:
        from gpa_phenotype_match import PhenotypeMatcher
        matcher = PhenotypeMatcher(offline_mode=offline_mode)
        if tracker:
            tracker.api_call("phase2", "phenotype", "phenotype_match", "started",
                             {"genes": len(genes_with_data)})
        match_results = await matcher.match_batch(genes_with_data, user_phenotypes)

        # Map results back by gene
        gene_results = {gene: mr for gene, mr in zip(genes_with_data, match_results)}
        # Add zero scores for genes without data
        for gene in genes_without_data:
            gene_results[gene] = {
                "score": 0.0,
                "explanation": "No known phenotypes found for this gene in local database.",
                "confidence": "low",
                "matched_pairs": [],
                "known_phenotypes": [],
            }

        for i in candidate_indices:
            v = variants[i]
            mr = gene_results.get(v.gene, {})
            v.phenotype_match_score = mr.get("score")
            v.phenotype_match_explanation = mr.get("explanation", "")
            v.phenotype_match_confidence = mr.get("confidence", "")
            v.phenotype_matched_pairs = mr.get("matched_pairs", [])
            v.phenotype_known_list = mr.get("known_phenotypes", [])
        if tracker:
            tracker.api_call("phase2", "phenotype", "phenotype_match", "completed",
                             {"genes": len(genes_with_data)})
    except Exception as e:
        print(f"[GPA Phase 2] Phenotype matching failed (non-critical): {type(e).__name__}: {e}")
        if tracker:
            tracker.error("phase2", "phenotype", f"Phenotype matching failed: {e}")


# ─── Main Two-Phase Entry Point ────────────────────────────────────────────

async def run_two_phase_pipeline(
    variants_data: List[Dict],
    config: Optional[GPAConfig] = None,
    user_phenotypes: Optional[str] = None,
    max_candidates: int = 150,
    progress_log_path: Optional[str] = None,
    tier1_only: bool = False,
) -> Dict[str, Any]:
    """
    Two-phase GPA pipeline optimized for large VCF datasets.

    Phase 1: Fast local triage (VEP data + local gene lists → preliminary tiers)
    Phase 2: API enrichment + final classification ONLY for Tier 1/2 candidates
        (or ONLY Tier 1 candidates when tier1_only=True)

    Args:
        progress_log_path: Optional path to a JSON Lines progress log file.
            If provided, fine-grained progress events are written here and
            can be polled by an external monitor.
        tier1_only: If True, restrict Phase 2 API enrichment and final
            classification to Tier 1 candidates only. This significantly
            reduces runtime when the user wants a focused high-confidence
            report first.

    Returns dict with report_markdown, summary, tier1/2/3_variants, etc.
    """
    try:
        return await _run_two_phase_pipeline_impl(
            variants_data, config, user_phenotypes, max_candidates,
            progress_log_path=progress_log_path,
            tier1_only=tier1_only,
        )
    except Exception as e:  # noqa: BROAD_EXCEPT — outer process-level guard: never let pipeline crash the host process
        import traceback
        err_msg = f"{type(e).__name__}: {e}"
        print(f"[GPA Two-Phase] PIPELINE CRASH: {err_msg}")
        traceback.print_exc()
        return {
            "error": err_msg,
            "report_markdown": f"# Analysis Failed\\n\\nPipeline crashed: `{err_msg}`\\n",
            "summary": {
                "tier1_variant_count": 0,
                "tier2_variant_count": 0,
                "tier3_variant_count": 0,
                "total_variants": len(variants_data),
                "error": err_msg,
            },
            "tier1_variants": [],
            "tier2_variants": [],
            "tier3_variants": [],
            "multi_hit_details": [],
            "meta": {"error": err_msg, "offline_mode": getattr(config, "offline_mode", False)},
            "json_report": {},
        }


async def _run_two_phase_pipeline_impl(
    variants_data: List[Dict],
    config: Optional[GPAConfig] = None,
    user_phenotypes: Optional[str] = None,
    max_candidates: int = 150,
    progress_log_path: Optional[str] = None,
    tier1_only: bool = False,
) -> Dict[str, Any]:
    """Internal implementation — protected by run_two_phase_pipeline wrapper."""
    if config is None:
        config = GPAConfig()

    # v0.10.16: Initialize fine-grained progress tracker
    from gpa_progress import GPAProgressTracker
    tracker = GPAProgressTracker(
        log_path=progress_log_path,
        enabled=progress_log_path is not None,
    )

    t0 = time.time()
    global_config = config.to_global()

    # v0.10.12: Preflight health check + per-API proxy route discovery
    from gpa_preflight import run_preflight_check, suggest_action
    preflight, route_map = await run_preflight_check(global_config)
    print(route_map.to_markdown())
    # Attach route map to global_config so DGRAAPIClient instances can pick it up
    global_config._proxy_route_map = route_map
    if not preflight.is_ready():
        action = suggest_action(preflight)
        if action == "abort":
            print("[GPA Preflight] 环境检查未通过，中止分析。")
            print(preflight.to_markdown())
            return {"error": "Preflight failed", "report": preflight.to_dict()}
        elif action == "offline":
            print("=" * 60)
            print("[GPA Preflight] ⚠️  WARNING: 切换到离线模式（跳过所有 API 调用）")
            print("[GPA Preflight] 离线模式下 VEP 注释不可用，大量变异将无法 tier 分级。")
            print("[GPA Preflight] 39,000+ 变异可能仅有极少数能被评估。")
            print("=" * 60)
            config.offline_mode = True
            global_config.offline_mode = True

    tissue_profile = config.get_tissue_profile()
    profile_name = tissue_profile.get("display_name", config.tissue_profile)

    # ── Phase 0: Parse all variants ──
    tracker.phase_start("phase0", "Parsing variants", {"input_variants": len(variants_data)})
    variants = _parse_variants_phase1(variants_data)
    unique_genes = list({v.gene for v in variants})
    n_total = len(variants)
    print(f"[GPA Two-Phase] {n_total} variants across {len(unique_genes)} unique genes")
    print(f"[GPA Two-Phase] Tissue profile: {profile_name} | Offline: {config.offline_mode}")
    tracker.phase_end("phase0", "Variants parsed", {
        "total_variants": n_total,
        "unique_genes": len(unique_genes),
        "tissue_profile": profile_name,
        "offline": config.offline_mode,
    })

    # ── Phase 1: Fast Local Triage ──
    tracker.phase_start("phase1", "Fast local triage", {"total_variants": n_total})
    phase1_start = time.time()

    # Phase 1.1: Pre-filter — identify candidate variants
    tracker.step_start("phase1", "prefilter", "Identifying candidate variants")
    candidate_indices = []
    progress_interval = max(1, n_total // 20)
    for i, v in enumerate(variants):
        if _is_potentially_pathogenic(v):
            candidate_indices.append(i)
        if (i + 1) % progress_interval == 0 or i == n_total - 1:
            tracker.step_progress("phase1", "prefilter", i + 1, n_total)

    n_candidates = len(candidate_indices)
    reduction = (1 - n_candidates / n_total) * 100 if n_total > 0 else 0
    tracker.step_end("phase1", "prefilter", "Candidate identification complete", {
        "candidates": n_candidates,
        "total": n_total,
        "reduction_percent": round(reduction, 2),
    })

    # Phase 1.2: Assign preliminary tiers to ALL variants
    tracker.step_start("phase1", "tier_assignment", "Assigning preliminary tiers")
    candidate_tier1 = []
    candidate_tier2 = []
    candidate_gene_set = set()
    for i in candidate_indices:
        tier, reason = _fast_tier_classification(variants[i])
        variants[i].tier = tier
        variants[i].tier_reason = reason
        if tier == 1:
            candidate_tier1.append(i)
        else:
            candidate_tier2.append(i)
        candidate_gene_set.add(variants[i].gene)
    tracker.step_end("phase1", "tier_assignment", "Preliminary tiers assigned", {
        "tier1": len(candidate_tier1),
        "tier2": len(candidate_tier2),
        "tier3": n_total - n_candidates,
        "candidate_genes": len(candidate_gene_set),
    })

    # v0.10.4: Warn if candidate count exceeds threshold, but do NOT truncate.
    # User wants to be notified rather than silently dropping variants.
    all_candidate_indices = list(set(candidate_tier1 + candidate_tier2))
    if len(all_candidate_indices) > max_candidates:
        tracker.warning("phase1", "candidate_threshold",
                        f"Candidate count {len(all_candidate_indices)} exceeds threshold {max_candidates}",
                        {"candidates": len(all_candidate_indices), "threshold": max_candidates,
                         "tier1": len(candidate_tier1), "tier2": len(candidate_tier2)})
        print("=" * 60)
        print(f"[GPA WARNING] Candidate count ({len(all_candidate_indices)}) exceeds "
              f"threshold={max_candidates}")
        print(f"              Tier 1: {len(candidate_tier1)}, Tier 2: {len(candidate_tier2)}")
        print(f"              All {len(all_candidate_indices)} candidates will proceed to "
              f"Phase 2 API enrichment (no truncation).")
        print("              This may take longer due to gnomAD rate limits.")
        print("=" * 60)
    else:
        all_candidate_indices = list(set(all_candidate_indices))

    # Non-candidates → Tier 3
    # v0.10.16 FIX: candidate_indices was a list, making the membership test
    # O(n*m) and causing hangs on large VCFs. Convert to set first.
    candidate_index_set = set(candidate_indices)
    non_candidate_indices = [i for i in range(n_total) if i not in candidate_index_set]
    for i in non_candidate_indices:
        variants[i].tier = 3
        variants[i].tier_reason = "Pre-filtered: low impact, no ClinVar pathogenicity, not in disease gene list"
        variants[i].tier_actions = ["Archive only"]

    print(f"[GPA Phase 1] Preliminary tiers: {len(candidate_tier1)} Tier 1, "
          f"{len(candidate_tier2)} Tier 2, {len(non_candidate_indices)} Tier 3 "
          f"({len(candidate_gene_set)} candidate genes)")
    tracker.phase_end("phase1", "Fast local triage complete", {
        "candidates": n_candidates,
        "tier1": len(candidate_tier1),
        "tier2": len(candidate_tier2),
        "tier3": len(non_candidate_indices),
        "candidate_genes": len(candidate_gene_set),
        "reduction_percent": round(reduction, 2),
    })

    # v0.10.16: Optional Tier-1-only mode for rapid focused reporting.
    # When enabled, drop Tier 2 candidates from Phase 2 enrichment and
    # mark them as preliminary Tier 3 (they will be re-classified if the
    # user later runs the full pipeline).
    if tier1_only:
        dropped_tier2 = len(candidate_tier2)
        if dropped_tier2 > 0:
            print(f"[GPA Two-Phase] Tier-1-only mode: dropping {dropped_tier2} Tier 2 candidates from Phase 2")
            for i in candidate_tier2:
                variants[i].tier = 3
                variants[i].tier_reason = "[TIER1_ONLY] Preliminary Tier 2 candidate excluded from Phase 2 enrichment"
                variants[i].tier_actions = ["Re-run full pipeline if broader review needed"]
            all_candidate_indices = list(candidate_tier1)
            candidate_gene_set = {variants[i].gene for i in all_candidate_indices}
            tracker.warning("phase1", "tier1_only",
                            f"Tier-1-only mode enabled: {dropped_tier2} Tier 2 candidates excluded",
                            {"tier1": len(candidate_tier1), "tier2_excluded": dropped_tier2})

    # If no candidates, skip Phase 2
    if not all_candidate_indices:
        print(f"[GPA Two-Phase] No candidates — skipping Phase 2 enrichment. "
              f"Total: {time.time() - t0:.1f}s")
        return _build_output(
            variants, [], tissue_profile, config, profile_name,
            "Phase 1 only (no candidates found)"
        )

    # ── Phase 2: Deep Enrichment for Candidates ONLY ──
    phase2_start = time.time()
    tracker.phase_start("phase2", "API enrichment for candidates", {
        "candidates": len(all_candidate_indices),
        "candidate_genes": len(candidate_gene_set),
        "tier1": len(candidate_tier1),
        "tier2": len(candidate_tier2),
        "tier1_only": tier1_only,
    })

    print(f"[GPA Phase 2] Starting enrichment for {len(all_candidate_indices)} candidates "
          f"({len(candidate_gene_set)} genes)...")

    # Phase 2.1: Per-gene API enrichment (Ensembl/UniProt/gnomAD_constraint)
    tracker.step_start("phase2", "gene_apis", "Ensembl/UniProt/gnomAD constraint enrichment",
                       {"candidate_genes": len(candidate_gene_set)})
    ensembl_data, uniprot_data, gnomad_constraint_data = await _enrich_candidate_genes(
        variants, all_candidate_indices, config, tissue_profile, tracker=tracker
    )
    tracker.step_end("phase2", "gene_apis", "Gene API enrichment complete", {
        "ensembl_genes": len(ensembl_data),
        "uniprot_genes": len(uniprot_data),
        "gnomad_constraint_genes": len(gnomad_constraint_data),
    })
    print(f"[GPA Phase 2] Gene APIs: Ensembl={len(ensembl_data)}, "
          f"UniProt={len(uniprot_data)}, gnomAD_constraint={len(gnomad_constraint_data)}")

    # Phase 2.1b: GTEx expression for candidate genes (v0.10.8)
    tracker.step_start("phase2", "gtex", "GTEx expression enrichment",
                       {"candidate_genes": len(candidate_gene_set)})
    gtex_data = await _enrich_gtex(
        variants, all_candidate_indices, config, tissue_profile, user_phenotypes, tracker=tracker
    )
    tracker.step_end("phase2", "gtex", "GTEx enrichment complete", {
        "gtex_genes": len(gtex_data),
    })
    print(f"[GPA Phase 2] GTEx: expression data for {len(gtex_data)} candidate genes")

    # Phase 2.2: gnomAD frequency enrichment
    tracker.step_start("phase2", "gnomad_freq", "gnomAD frequency enrichment",
                       {"candidates": len(all_candidate_indices)})
    await _enrich_variant_frequencies(variants, all_candidate_indices, config, tracker=tracker)
    tracker.step_end("phase2", "gnomad_freq", "gnomAD frequency enrichment complete")

    # Phase 2.3: SpliceAI (auto-enabled for Tier 1/2 candidates in v0.10.7/0.10.13)
    tracker.step_start("phase2", "spliceai", "SpliceAI enrichment",
                       {"candidates": len(all_candidate_indices)})
    if all_candidate_indices and config.spliceai_enabled:
        await _enrich_spliceai(variants, all_candidate_indices, config, tracker=tracker)
    else:
        tracker.step_end("phase2", "spliceai", "SpliceAI skipped (disabled)")
    tracker.step_end("phase2", "spliceai", "SpliceAI enrichment complete")

    # Phase 2.4: Phenotype LLM matching
    tracker.step_start("phase2", "phenotype", "Phenotype matching",
                       {"candidates": len(all_candidate_indices)})
    await _enrich_phenotype(variants, all_candidate_indices, user_phenotypes, tracker=tracker, config=config)
    tracker.step_end("phase2", "phenotype", "Phenotype matching complete")

    # Phase 2.5: Process candidate genes with enriched data
    tracker.step_start("phase2", "post_process", "Post-processing candidates",
                       {"candidates": len(all_candidate_indices)})
    # NMD prediction for truncating variants
    lof_terms = {"frameshift", "nonsense", "stop_gained", "start_lost"}
    for i in all_candidate_indices:
        v = variants[i]
        if any(term in v.consequence.lower() for term in lof_terms):
            v.nmd_prediction = predict_nmd(v, ensembl_data.get(v.gene, {}))
            v.gene_constraint = {"nmd_prediction": v.nmd_prediction}

    # Transcript correction
    for i in all_candidate_indices:
        v = variants[i]
        v, warning = await correct_transcript_priority(v, ensembl_data)
        if warning:
            v.transcript_warning = json.dumps(warning)

    # Pseudogene detection
    for i in all_candidate_indices:
        pg_warning = detect_pseudogene_artifact(variants[i])
        if pg_warning:
            variants[i].pseudogene_warning = json.dumps(pg_warning)

    # gnomAD classification
    for i in all_candidate_indices:
        v = variants[i]
        gnomad_info = classify_gnomad_frequency(
            v.gnomad_af, v.gene,
            af_by_population=getattr(v, 'gnomad_populations', None),
            target_population=getattr(config, 'target_population', None)
        )
        v.gnomad_status = gnomad_info.get("status", v.gnomad_status or "UNKNOWN")

    # Domain mapping
    for i in all_candidate_indices:
        variants[i].domain_info = map_variant_to_domain(variants[i], uniprot_data)

    # Tissue relevance (v0.10.7: GTEx data now available from Phase 2.1b)
    tissue_assessments = {}
    for i in all_candidate_indices:
        v = variants[i]
        tissue = assess_tissue_relevance(v, tissue_profile, gtex_data)
        tissue_assessments[v.gene] = tissue
        v.tissue_relevance = tissue
    tracker.step_end("phase2", "post_process", "Post-processing complete")

    # Phase 2.6: Final tier classification for candidates
    tracker.step_start("phase2", "final_classification", "Final tier classification",
                       {"candidates": len(all_candidate_indices)})
    from gpa_tier_classifier import classify_variant_tier
    n_classified = 0
    progress_interval = max(1, len(all_candidate_indices) // 10)
    for i in all_candidate_indices:
        v = variants[i]
        tissue = tissue_assessments.get(v.gene, {})
        gnomad_info = classify_gnomad_frequency(
            v.gnomad_af, v.gene,
            af_by_population=getattr(v, 'gnomad_populations', None),
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
        n_classified += 1
        if tracker and n_classified % progress_interval == 0:
            tracker.step_progress("phase2", "final_classification", n_classified, len(all_candidate_indices))
    tracker.step_end("phase2", "final_classification", "Final classification complete",
                     {"classified": n_classified})

    phase2_duration = time.time() - phase2_start
    total_duration = time.time() - t0
    print(f"[GPA Phase 2] Enrichment complete in {phase2_duration:.1f}s")
    print(f"[GPA Two-Phase] Total: {total_duration:.1f}s")
    tracker.phase_end("phase2", "API enrichment complete", {
        "duration_sec": round(phase2_duration, 2),
    })

    tracker.phase_start("report", "Report generation", {"variants": len(variants)})
    output = _build_output(
        variants, all_candidate_indices, tissue_profile, config, profile_name,
        f"Two-phase: Phase 1 filtered {n_total}→{n_candidates} candidates, "
        f"Phase 2 enriched {len(candidate_gene_set)} genes",
        gtex_data=gtex_data,
    )
    tracker.phase_end("report", "Report generation complete", {
        "tier1": output["summary"]["tier1_variant_count"],
        "tier2": output["summary"]["tier2_variant_count"],
        "tier3": output["summary"]["tier3_variant_count"],
        "total_duration_sec": round(total_duration, 2),
    })
    tracker.finish()
    return output


def _build_output(
    variants: List[Variant],
    enriched_indices: List[int],
    tissue_profile: Dict,
    config: GPAConfig,
    profile_name: str,
    method_note: str,
    gtex_data: Optional[Dict[str, Dict]] = None,
) -> Dict[str, Any]:
    """Build standardized output dict for report generation."""
    from gpa_report import generate_tier_report, generate_json_report
    from gpa_qc import _run_qc_checks

    enriched_set = set(enriched_indices)

    tier1 = [v for v in variants if v.tier == 1]
    tier2 = [v for v in variants if v.tier == 2]
    tier3 = [v for v in variants if v.tier == 3]

    # Multi-hit detection on all variants
    from gpa_multi_hit import detect_multi_hit_genes
    multi_hits = detect_multi_hit_genes(variants, {})

    # QC checks
    qc_summary = _run_qc_checks(variants)

    # Mark enrichment status
    # v0.10.4: Use enumerate() to avoid O(n²) variants.index(v) lookup
    for idx, v in enumerate(variants):
        if idx not in enriched_set and v.tier == 3:
            v.tier_reason = f"[PRELIMINARY] {v.tier_reason}"

    report_md = generate_tier_report(
        variants, config, tissue_profile, multi_hits, gtex_data=gtex_data
    )

    json_report = generate_json_report(
        variants, config, tissue_profile, multi_hits, report_md, qc_summary
    )

    # Build summary
    summary = {
        "tier1_gene_count": len({v.gene for v in tier1}),
        "tier1_variant_count": len(tier1),
        "tier2_gene_count": len({v.gene for v in tier2}),
        "tier2_variant_count": len(tier2),
        "tier3_gene_count": len({v.gene for v in tier3}),
        "tier3_variant_count": len(tier3),
        "multi_hit_genes": [mh["gene"] for mh in multi_hits],
        "total_variants": len(variants),
        "enriched_candidates": len(enriched_indices),
        "method": method_note,
    }

    return {
        "report_markdown": report_md,
        "summary": summary,
        "tier1_variants": [v.__dict__ for v in tier1],
        "tier2_variants": [v.__dict__ for v in tier2],
        "tier3_variants": [v.__dict__ for v in tier3],
        "multi_hit_details": multi_hits,
        "meta": {
            "tissue_profile": config.tissue_profile,
            "profile_display_name": profile_name,
            "total_variants": len(variants),
            "offline_mode": config.offline_mode,
            "method": method_note,
        },
        "json_report": json_report,
    }
