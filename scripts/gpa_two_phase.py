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
            chrom=vd.get("CHROM", ""),
            pos=int(vd.get("POS", 0) or 0),
            ref=vd.get("REF", ""),
            alt=vd.get("ALT", ""),
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
) -> Tuple[Dict, Dict, Dict]:
    """
    Phase 2: API enrichment ONLY for candidate genes.
    Returns (ensembl_data, uniprot_data, gnomad_constraint_data).
    """
    if not candidate_indices:
        return {}, {}, {}

    candidate_genes = list({variants[i].gene for i in candidate_indices})
    global_config = config.to_global()

    if config.offline_mode:
        ensembl_data, uniprot_data, gnomad_constraint_data = {}, {}, {}
        for gene in candidate_genes:
            archive = _load_offline_archive(gene)
            if archive:
                ensembl_data[gene] = archive.get("ensembl", {})
                uniprot_data[gene] = archive.get("uniprot", {})
                gc = archive.get("gnomad_constraint")
                if gc and gc.get("status") == "CAPTURED":
                    gnomad_constraint_data[gene] = gc
        return ensembl_data, uniprot_data, gnomad_constraint_data

    cache = DGRACache(global_config.cache_db_path)
    async with DGRAAPIClient(global_config, cache) as client:
        ensembl_raw, uniprot_raw, gnomad_constraint_raw = await asyncio.gather(
            client.batch_query_genes(candidate_genes, "ensembl"),
            client.batch_query_genes(candidate_genes, "uniprot"),
            client.batch_query_genes(candidate_genes, "gnomad_constraint"),
        )
        ensembl_data = {g: ensembl_raw.get(g, {}) for g in candidate_genes}
        uniprot_data = {g: uniprot_raw.get(g, {}) for g in candidate_genes}
        gnomad_constraint_data = {g: gnomad_constraint_raw.get(g, {}) for g in candidate_genes}

    return ensembl_data, uniprot_data, gnomad_constraint_data


async def _enrich_gtex(
    variants: List[Variant],
    candidate_indices: List[int],
    config: GPAConfig,
    tissue_profile: Dict,
    user_phenotypes: Optional[str] = None,
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
    global_config = config.to_global()
    
    cache = DGRACache(global_config.cache_db_path)
    gtex_data: Dict[str, Dict] = {}
    
    async with DGRAAPIClient(global_config, cache) as client:
        # Semaphore to be polite to GTEx API
        gtex_sem = asyncio.Semaphore(5)
        
        async def _query_one_gene(gene: str) -> Tuple[str, Optional[Dict]]:
            async with gtex_sem:
                try:
                    # v0.10.8: Query ALL tissues in a single batch call
                    results = await client.query_gtex_expression_multi(gene, ALL_GTEX_TISSUES)
                    if not results:
                        return gene, None
                    
                    # Build tissue→TPM map
                    tissue_tpm = {}
                    for r in results:
                        t = r.get("tissue")
                        tpm = r.get("median_tpm")
                        if t is not None and tpm is not None:
                            tissue_tpm[t] = tpm
                    
                    if not tissue_tpm:
                        return gene, None
                    
                    # Global max across all tissues
                    global_max = max(tissue_tpm.values())
                    
                    # Phenotype-relevant max
                    phenotype_values = [tissue_tpm.get(t) for t in phenotype_tissues if t in tissue_tpm]
                    phenotype_max = max(phenotype_values) if phenotype_values else 0.0
                    
                    # All tissues with TPM > 0, sorted by TPM desc
                    all_expressing = sorted(
                        [(t, tpm) for t, tpm in tissue_tpm.items() if tpm > 0],
                        key=lambda x: x[1], reverse=True
                    )
                    
                    expressing_count = len(all_expressing)
                    
                    return gene, {
                        "median_tpm": global_max,  # Use global max as representative
                        "max_tpm": global_max,
                        "all_tissues": [{"tissue": t, "tpm": tpm} for t, tpm in all_expressing],
                        "expressing_tissues": expressing_count,
                        "source": "gtex_multi",
                        "confidence": "medium",
                        # v0.10.8: Phenotype-relevance data
                        "phenotype_max_tpm": phenotype_max,
                        "phenotype_tissues": phenotype_tissues,
                        "phenotype_matched_keywords": phenotype_matched_keywords,
                        "global_max_tpm": global_max,
                    }
                except Exception as e:
                    print(f"[GPA Phase 2] GTEx warning for {gene}: {type(e).__name__}: {e}")
                    return gene, None
        
        gtex_results = await asyncio.gather(*[_query_one_gene(g) for g in candidate_genes])
        for gene, data in gtex_results:
            if data:
                gtex_data[gene] = data
    
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

    # Try MyVariant.info batch first — query Tier 1/2 candidates for gnomAD AF + ClinVar + CADD.
    # v0.10.4: All candidate variants get ClinVar/CADD regardless of consequence.
    # Frameshift/stop_gained may still have ClinVar conflict annotations or
    # additional phenotypic evidence that affects tier weighting.
    # v0.10.13: Restricted to Tier 1/2 only to reduce API calls.
    try:
        from dgra_myvariant import query_myvariant_batch, apply_myvariant_results
        mv_sem = asyncio.Semaphore(10)
        timeout_obj = aiohttp.ClientTimeout(total=120)
        mv_variants = [(v.chrom, v.pos, v.ref, v.alt) for v in tier12_candidates_no_af]
        # v0.10.4: Always fill ClinVar and CADD for all candidates
        # v0.10.13: Only for Tier 1/2 candidates
        async with aiohttp.ClientSession(timeout=timeout_obj, trust_env=False) as mv_session:
            mv_results = await query_myvariant_batch(mv_variants, mv_session, semaphore=mv_sem, batch_size=1000)
        mv_stats = apply_myvariant_results(tier12_candidates_no_af, mv_results)
        print(f"[GPA Phase 2] MyVariant.info: {mv_stats['gnomad_filled']} gnomAD, "
              f"{mv_stats['clinvar_filled']} ClinVar, "
              f"{mv_stats['cadd_filled']} CADD filled")
    except Exception as e:
        print(f"[GPA Phase 2] MyVariant.info batch query failed (non-critical): {type(e).__name__}: {e}")

    # Fallback: gnomAD GraphQL ONLY for Tier 1/2 candidates still without AF
    still_no_af = [v for v in tier12_candidates_no_af if v.gnomad_af is None]
    if not still_no_af:
        return

    print(f"[GPA Phase 2] gnomAD GraphQL: querying {len(still_no_af)} Tier 1/2 candidates without AF")
    cache = DGRACache(global_config.cache_db_path)
    async with DGRAAPIClient(global_config, cache) as client:
        gnomad_sem = asyncio.Semaphore(2)
        async def _query_one(v):
            async with gnomad_sem:
                try:
                    return await client.query_gnomad_variant(v.chrom, v.pos, v.ref, v.alt)
                except Exception as e:
                    return {"status": "API_FAILED", "error": str(e), "source": "failed"}
        results = await asyncio.gather(*[_query_one(v) for v in still_no_af])
        for v, result in zip(still_no_af, results):
            if result and result.get("source") in ("gnomad", "cache", "failed"):
                af = result.get("af")
                if af is not None:
                    v.gnomad_af = af
                    v.gnomad_populations = result.get("af_populations", {})
                    v.gnomad_status = "SUCCESS"

        # v0.10.6 FIX: NCBI ClinVar direct query for variants still UNKNOWN after MyVariant.info.
        # MyVariant.info ClinVar coverage is incomplete; NCBI ESummary provides accurate data.
        # v0.10.13: Only query Tier 1/2 variants for ClinVar.
        clinvar_unknown = [variants[i] for i in candidate_indices
                          if variants[i].clinvar in (_UNKNOWN, "UNKNOWN", "") and getattr(variants[i], 'tier', 3) in (1, 2)]
        if clinvar_unknown:
            # Consequence-aware filtering: skip intron/UTR where ClinVar is rarely informative
            cv_candidates = [v for v in clinvar_unknown if _variant_needs_clinvar(v)]
            if cv_candidates:
                n_cv = len(cv_candidates)
                print(f"[GPA Phase 2] ClinVar (NCBI): querying {n_cv} variants (1 req/s)")
                cv_sem = asyncio.Semaphore(1)  # NCBI: 1 req/s
                async def _query_one_clinvar(v):
                    async with cv_sem:
                        try:
                            return await client.query_ncbi_clinvar(
                                gene=v.gene, chrom=v.chrom, pos=v.pos
                            )
                        except Exception as e:
                            return {"clinical_significance": None, "error": str(e), "source": "failed"}
                cv_results = await asyncio.gather(*[_query_one_clinvar(v) for v in cv_candidates])
                n_found = 0
                for v, cv_result in zip(cv_candidates, cv_results):
                    sig = cv_result.get("clinical_significance")
                    if sig:
                        v.clinvar = sig
                        v.clinvar_review_status = cv_result.get("review_status")
                        n_found += 1
                print(f"[GPA Phase 2] ClinVar (NCBI): {n_found}/{n_cv} filled")


async def _enrich_spliceai(
    variants: List[Variant],
    candidate_indices: List[int],
    config: GPAConfig,
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

    spliceai_sem = asyncio.Semaphore(getattr(config, 'spliceai_concurrency', 5))
    spliceai_results = await query_spliceai_batch(
        tier12_candidates, spliceai_sem,
        timeout=getattr(config, 'spliceai_timeout', 45),
    )

    # Attach SpliceAI results to all Tier 1/2 candidates.
    # For variants without splice changes, the API returns delta=0 or not_in_db,
    # both of which are valid evidence entries.
    for v in tier12_candidates:
        key = _splice_key(v.chrom, v.pos, v.ref, v.alt)
        if key in spliceai_results:
            v.spliceai_result = spliceai_results[key]
        else:
            # Not queried or not in results — mark as not_in_db with null scores.
            # delta_score=None is used by the tier classifier to distinguish
            # "not queried" from "queried and delta=0".
            v.spliceai_result = {"source": "not_in_db", "delta_score": None, "predicted_impact": None}

    print(f"[GPA Phase 2] SpliceAI: batch complete for {len(tier12_candidates)} variants")


async def _enrich_phenotype(
    variants: List[Variant],
    candidate_indices: List[int],
    user_phenotypes: Optional[str],
) -> None:
    """Phase 2: Phenotype LLM matching ONLY for candidate genes with local DB data."""
    if not user_phenotypes or not candidate_indices:
        return

    candidate_genes = sorted({variants[i].gene for i in candidate_indices})

    # v0.10.3: Pre-filter — only query genes that have known phenotype data locally
    try:
        from gpa_phenotype_match import PhenotypeMatcher
        matcher = PhenotypeMatcher()
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
            for i in candidate_indices:
                v = variants[i]
                v.phenotype_match_score = 0.0
                v.phenotype_match_explanation = "No known phenotypes found for this gene in local database."
                v.phenotype_match_confidence = "low"
            return
    except Exception as e:
        print(f"[GPA Phase 2] Phenotype pre-filter error (proceeding with all): {e}")
        genes_with_data = candidate_genes
        genes_without_data = []

    print(f"[GPA Phase 2] Phenotype matching: {len(genes_with_data)} candidate genes (with local data)")

    try:
        from gpa_phenotype_match import PhenotypeMatcher
        matcher = PhenotypeMatcher()
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
    except Exception as e:
        print(f"[GPA Phase 2] Phenotype matching failed (non-critical): {type(e).__name__}: {e}")


# ─── Main Two-Phase Entry Point ────────────────────────────────────────────

async def run_two_phase_pipeline(
    variants_data: List[Dict],
    config: Optional[GPAConfig] = None,
    user_phenotypes: Optional[str] = None,
    max_candidates: int = 150,
) -> Dict[str, Any]:
    """
    Two-phase GPA pipeline optimized for large VCF datasets.

    Phase 1: Fast local triage (VEP data + local gene lists → preliminary tiers)
    Phase 2: API enrichment + final classification ONLY for Tier 1/2 candidates

    Returns dict with report_markdown, summary, tier1/2/3_variants, etc.
    """
    try:
        return await _run_two_phase_pipeline_impl(
            variants_data, config, user_phenotypes, max_candidates
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
) -> Dict[str, Any]:
    """Internal implementation — protected by run_two_phase_pipeline wrapper."""
    if config is None:
        config = GPAConfig()

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
    variants = _parse_variants_phase1(variants_data)
    unique_genes = list({v.gene for v in variants})
    n_total = len(variants)
    print(f"[GPA Two-Phase] {n_total} variants across {len(unique_genes)} unique genes")
    print(f"[GPA Two-Phase] Tissue profile: {profile_name} | Offline: {config.offline_mode}")

    # ── Phase 1: Fast Local Triage ──
    phase1_start = time.time()

    # Phase 1.1: Pre-filter — identify candidate variants
    candidate_indices = []
    for i, v in enumerate(variants):
        if _is_potentially_pathogenic(v):
            candidate_indices.append(i)

    n_candidates = len(candidate_indices)
    reduction = (1 - n_candidates / n_total) * 100 if n_total > 0 else 0
    print(f"[GPA Phase 1] Fast triage complete: {n_candidates}/{n_total} candidates "
          f"({reduction:.1f}% reduction) in {time.time() - phase1_start:.1f}s")

    # Phase 1.2: Assign preliminary tiers to ALL variants
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

    # v0.10.4: Warn if candidate count exceeds threshold, but do NOT truncate.
    # User wants to be notified rather than silently dropping variants.
    all_candidate_indices = list(set(candidate_tier1 + candidate_tier2))
    if len(all_candidate_indices) > max_candidates:
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
    non_candidate_indices = [i for i in range(n_total) if i not in candidate_indices]
    for i in non_candidate_indices:
        variants[i].tier = 3
        variants[i].tier_reason = "Pre-filtered: low impact, no ClinVar pathogenicity, not in disease gene list"
        variants[i].tier_actions = ["Archive only"]

    print(f"[GPA Phase 1] Preliminary tiers: {len(candidate_tier1)} Tier 1, "
          f"{len(candidate_tier2)} Tier 2, {len(non_candidate_indices)} Tier 3 "
          f"({len(candidate_gene_set)} candidate genes)")

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

    print(f"[GPA Phase 2] Starting enrichment for {len(all_candidate_indices)} candidates "
          f"({len(candidate_gene_set)} genes)...")

    # Phase 2.1: Per-gene API enrichment (Ensembl/UniProt/gnomAD_constraint)
    ensembl_data, uniprot_data, gnomad_constraint_data = await _enrich_candidate_genes(
        variants, all_candidate_indices, config, tissue_profile
    )
    print(f"[GPA Phase 2] Gene APIs: Ensembl={len(ensembl_data)}, "
          f"UniProt={len(uniprot_data)}, gnomAD_constraint={len(gnomad_constraint_data)}")

    # Phase 2.1b: GTEx expression for candidate genes (v0.10.8)
    # v0.10.8: Queries ALL GTEx tissues + phenotype-tissue association
    gtex_data = await _enrich_gtex(
        variants, all_candidate_indices, config, tissue_profile, user_phenotypes
    )
    print(f"[GPA Phase 2] GTEx: expression data for {len(gtex_data)} candidate genes")

    # Phase 2.2: gnomAD frequency enrichment
    await _enrich_variant_frequencies(variants, all_candidate_indices, config)

    # Phase 2.3: SpliceAI (auto-enabled for Tier 1/2 candidates in v0.10.7/0.10.13)
    # v0.10.13: Enable SpliceAI flag so that classify_variant_tier uses the results.
    # The two-phase pipeline auto-queries SpliceAI for all Tier 1/2 candidates,
    # so we set spliceai_enabled=True to ensure the classifier evaluates upgrades/downgrades.
    if all_candidate_indices:
        config.spliceai_enabled = True
    await _enrich_spliceai(variants, all_candidate_indices, config)

    # Phase 2.4: Phenotype LLM matching
    await _enrich_phenotype(variants, all_candidate_indices, user_phenotypes)

    # Phase 2.5: Process candidate genes with enriched data
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

    # Phase 2.6: Final tier classification for candidates
    from gpa_tier_classifier import classify_variant_tier
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

    phase2_duration = time.time() - phase2_start
    total_duration = time.time() - t0
    print(f"[GPA Phase 2] Enrichment complete in {phase2_duration:.1f}s")
    print(f"[GPA Two-Phase] Total: {total_duration:.1f}s")

    return _build_output(
        variants, all_candidate_indices, tissue_profile, config, profile_name,
        f"Two-phase: Phase 1 filtered {n_total}→{n_candidates} candidates, "
        f"Phase 2 enriched {len(candidate_gene_set)} genes"
    )


def _build_output(
    variants: List[Variant],
    enriched_indices: List[int],
    tissue_profile: Dict,
    config: GPAConfig,
    profile_name: str,
    method_note: str,
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
        variants, config, tissue_profile, multi_hits
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
