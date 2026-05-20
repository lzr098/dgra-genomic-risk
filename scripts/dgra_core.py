#!/usr/bin/env python3
"""
DGRA Core Engine - Donor Genomic Risk Assessment
v0.4 - 2026-05-19

API-first pipeline with cache layer. Replaces hardcoded gene dictionaries
with live queries to Ensembl, UniProt, GTEx, and gnomAD.

Phase 2 changes:
- MANE_SELECT / PROTEIN_DOMAINS / tissue gene lists → API queries
- run_dgra_pipeline is now async, with batch concurrent API calls
- --offline mode skips APIs, uses cache + local overrides only
"""

import json
import sys
import csv
import re
import asyncio
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple
import argparse

# v0.4 API layer imports
from dgra_config import DGRAGlobalConfig
from dgra_cache import DGRACache
from dgra_api import DGRAAPIClient

# =============================================================================
# Offline Archive Persistence (v0.4)
# =============================================================================

OFFLINE_ARCHIVE_DIR = Path(__file__).parent.parent / "references" / "offline_data"

def _save_offline_archive(gene: str, ensembl_data: dict, uniprot_data: dict,
                          gtex_data: dict, tissue_profile: str):
    """Persist API results for future offline use."""
    OFFLINE_ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    archive = {
        "gene": gene,
        "tissue_profile": tissue_profile,
        "saved_at": datetime.utcnow().isoformat() + "Z",
        "ensembl": ensembl_data.get(gene, {}),
        "uniprot": uniprot_data.get(gene, {}),
        "gtex": gtex_data.get(gene, {}),
    }
    path = OFFLINE_ARCHIVE_DIR / f"{gene}.json"
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(archive, f, indent=2, ensure_ascii=False, default=str)

def _load_offline_archive(gene: str) -> Optional[Dict]:
    """Load previously saved API results for offline mode."""
    path = OFFLINE_ARCHIVE_DIR / f"{gene}.json"
    if not path.exists():
        return None
    try:
        with open(path, 'r', encoding='utf-8') as f:
            return json.load(f)
    except Exception:
        return None

# =============================================================================
# Data Structures
# =============================================================================

@dataclass
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    gene: str
    transcript: str
    exon: str  # e.g., "E5/15"
    impact: str  # HIGH, MODERATE, LOW
    consequence: str  # frameshift, missense, splice_donor, etc.
    hgvsp: str  # protein change, e.g., "p.Thr111SerfsTer22"
    hgvsc: str  # cDNA change
    clinvar: str  # Pathogenic, Benign, etc.
    gnomad_af: Optional[float] = None
    dp: int = 0
    gq: float = 0.0
    gt: str = ""  # 0/1, 1/1, etc.
    vaf: Optional[float] = None
    
    # Computed fields
    tier: Optional[int] = None
    tier_reason: str = ""
    tier_actions: List[str] = field(default_factory=list)
    domain_info: Optional[Dict] = None
    transcript_warning: Optional[str] = None
    pseudogene_warning: Optional[str] = None
    gnomad_status: Optional[str] = None
    tissue_relevance: Optional[Dict] = None

@dataclass
class DGRAConfig:
    """User-facing config (kept simple). Maps to DGRAGlobalConfig internally.
    v0.4-fix: tissue_profile has NO default — must be specified by user.
    """
    min_dp: int = 20
    min_gq: float = 90.0
    common_af_threshold: float = 0.01
    low_af_threshold: float = 0.001
    vaf_deviation_threshold: float = 0.20
    tissue_profile: Optional[str] = None  # NO default — caller must specify
    offline_mode: bool = False

    def to_global(self) -> DGRAGlobalConfig:
        gc = DGRAGlobalConfig()
        gc.tissue_profile = self.tissue_profile or ""
        gc.offline_mode = self.offline_mode
        gc.min_dp = self.min_dp
        gc.min_gq = self.min_gq
        gc.common_af_threshold = self.common_af_threshold
        gc.low_af_threshold = self.low_af_threshold
        gc.vaf_deviation_threshold = self.vaf_deviation_threshold
        return gc

    def get_tissue_profile(self) -> Dict:
        """Load tissue profile from references/tissue_context.json.
        Raises if tissue_profile is not set.
        """
        if not self.tissue_profile:
            raise ValueError(
                "tissue_profile is required. Available profiles: hematopoietic, cardiovascular, "
                "hepatic, renal, neurological. Specify via --tissue or config.tissue_profile."
            )
        ref_path = Path(__file__).parent.parent / "references" / "tissue_context.json"
        with open(ref_path, 'r') as f:
            data = json.load(f)
        profiles = data.get("profiles", {})
        if self.tissue_profile not in profiles:
            available = ", ".join(profiles.keys())
            raise ValueError(f"Unknown tissue profile '{self.tissue_profile}'. Available: {available}")
        return profiles[self.tissue_profile]

# =============================================================================
# Module A: Transcript Priority Correction  (v0.4: Ensembl API)
# =============================================================================

async def correct_transcript_priority(variant: Variant,
                                       ensembl_data: Dict[str, Dict]) -> Tuple[Variant, Optional[Dict]]:
    """
    Check if annotator selected non-canonical isoform.
    Uses Ensembl API canonical transcript; falls back to input VCF transcript.
    
    Args:
        variant: Variant object
        ensembl_data: Pre-fetched Ensembl gene data {gene: {canonical_transcript, ...}}
    
    Returns: (corrected_variant, warning_dict or None)
    """
    gene = variant.gene
    selected_tx = variant.transcript.split()[0] if variant.transcript else ""
    
    # Get canonical from Ensembl API
    ens = ensembl_data.get(gene, {})
    canonical = ens.get("canonical_transcript", "")
    
    if not canonical:
        # No API data: fall back to trusting annotator's choice
        return variant, None
    
    # If already canonical
    if selected_tx.startswith(canonical):
        return variant, None
    
    # Annotator selected non-canonical isoform
    warning = {
        "type": "TRANSCRIPT_DISCREPANCY",
        "gene": gene,
        "annotator_selected": selected_tx,
        "canonical": canonical,
        "message": f"Annotator selected {selected_tx}, but Ensembl canonical {canonical} should be used.",
        "recommendation": "Re-assess impact using canonical transcript.",
        "source": ens.get("source", "unknown"),
        "confidence": ens.get("confidence", "medium"),
    }
    
    variant.transcript_warning = json.dumps(warning)
    return variant, warning

# =============================================================================
# Module B: Pseudogene Interference Detection  (unchanged in v0.4)
# =============================================================================

PSEUDOGENE_CONFIG = {
    "SETBP1": {"pseudogene": "VWFP1", "vaf_expected": 0.5, "vaf_min": 0.30, "notes": "VWFP1 on chr1 shares homology"},
    "VWF": {"pseudogene": "VWFP1", "vaf_expected": 0.5, "vaf_min": 0.25, "notes": "VWFP1 exon 23-34 homology"},
    "GBA": {"pseudogene": "GBAP1", "vaf_expected": 0.5, "vaf_min": 0.20, "notes": "GBAP1 on chr1"},
    "PMS2": {"pseudogene": "PMS2CL", "vaf_expected": 0.5, "vaf_min": 0.25, "notes": "PMS2CL on chr7"},
}

def detect_pseudogene_artifact(variant: Variant) -> Optional[Dict]:
    """
    Detect if VAF deviation suggests pseudogene interference.
    """
    gene = variant.gene
    if gene not in PSEUDOGENE_CONFIG:
        return None
    
    config = PSEUDOGENE_CONFIG[gene]
    
    # Determine expected VAF
    if variant.gt in ["1/1", "1|1"]:
        expected_vaf = 1.0
    else:
        expected_vaf = config["vaf_expected"]
    
    observed_vaf = variant.vaf
    if observed_vaf is None:
        # Try to estimate from DP if available
        return None
    
    deviation = abs(observed_vaf - expected_vaf)
    
    # For heterozygous, check if VAF is too low
    if expected_vaf == 0.5 and observed_vaf < config["vaf_min"]:
        return {
            "type": "PSEUDOGENE_SUSPECTED",
            "gene": gene,
            "pseudogene": config["pseudogene"],
            "expected_vaf": expected_vaf,
            "observed_vaf": observed_vaf,
            "deviation": deviation,
            "notes": config["notes"],
            "recommendation": "Consider long-read sequencing or Sanger validation."
        }
    
    return None

# =============================================================================
# Module C: gnomAD Frequency Handler
# =============================================================================

def classify_gnomad_frequency(af: Optional[float], gene: str) -> Dict:
    """
    Classify gnomAD allele frequency with zero-frequency handling.
    """
    # Genes with germline warnings (clonal hematopoiesis filtering)
    GERMLINE_WARNING_GENES = {"ASXL1", "DNMT3A", "TET2", "TP53", "PPM1D", "JAK2", "CBL", "IDH1", "IDH2"}

    if af is None:
        return {
            "af": "NOT_CAPTURED",
            "status": "gnomAD database does not capture this locus",
            "interpretation": "Population frequency data unavailable. Cannot judge benign based on gnomAD.",
            "action": "Continue with other modules (domain, ClinVar, zygosity).",
            "risk_adjustment": "Do NOT downgrade risk due to gnomAD absence."
        }

    if gene in GERMLINE_WARNING_GENES:
        return {
            "af": af,
            "status": "GERMLINE_WARNING",
            "interpretation": f"{gene} is filtered in gnomAD germline due to clonal hematopoiesis.",
            "action": "Use alternative databases (ExAC pre-filtered, ClinVar, LOVD).",
            "risk_adjustment": "Evaluate independently of gnomAD."
        }

    if af > 0.01:
        return {
            "af": af,
            "status": "common_polymorphism",
            "interpretation": "High population frequency (>1%), likely benign.",
            "action": "Require strong evidence to upgrade risk.",
            "risk_adjustment": "Default Tier 3 unless other evidence is strong."
        }
    elif af > 0.001:
        return {
            "af": af,
            "status": "low_frequency",
            "interpretation": "Moderate frequency (0.1-1%), needs functional assessment.",
            "action": "Include in comprehensive evaluation."
        }
    else:
        return {
            "af": af if af else "extremely_rare",
            "status": "rare_variant",
            "interpretation": "Very rare in population (<0.1%), likely under negative selection.",
            "action": "Focus attention, requires literature search.",
            "risk_adjustment": "Do NOT downgrade risk by default."
        }

# =============================================================================
# Module D: Protein Domain Mapping
# =============================================================================

# =============================================================================
# Module D: Protein Domain Mapping  (v0.4: UniProt API)
# =============================================================================

# v0.3 hardcoded PROTEIN_DOMAINS dict removed — now fetched from UniProt API
# Protein domains are loaded via dgra_api.query_uniprot_by_gene() and passed
# as uniprot_data dict to this function.

def parse_protein_position(hgvsp: str) -> Optional[int]:
    """Extract amino acid position from p. string. Handles NP_ prefix."""
    if not hgvsp or hgvsp == "":
        return None
    # Strip NP_ prefix if present: NP_000543.3:p.Val1565Leu -> p.Val1565Leu
    if ':p.' in hgvsp:
        hgvsp = hgvsp.split(':p.', 1)[1]
        hgvsp = 'p.' + hgvsp
    # Match p.XXX123 or p.Ala123 or p.123
    match = re.search(r'p\.[A-Za-z]+(\d+)', hgvsp)
    if match:
        return int(match.group(1))
    match = re.search(r'p\.(\d+)', hgvsp)
    if match:
        return int(match.group(1))
    return None

def map_variant_to_domain(variant: Variant, uniprot_data: Dict[str, Dict]) -> Dict:
    """
    Map variant to protein functional domain using UniProt API data.
    
    Args:
        variant: Variant object
        uniprot_data: Pre-fetched UniProt data {gene: {domains: [...], sequence_length: N}}
    
    Returns:
        Domain mapping dict with damage assessment.
    """
    gene = variant.gene
    up = uniprot_data.get(gene, {})
    domains = up.get("domains", [])
    seq_length = up.get("sequence_length")
    aa_pos = parse_protein_position(variant.hgvsp)
    
    # No API data available
    if not domains:
        if seq_length and aa_pos and aa_pos > seq_length:
            return {
                "domain": "unknown",
                "note": f"Parsed aa{aa_pos} exceeds UniProt sequence length ({seq_length}). "
                        "Possible transcript mismatch.",
                "gene": gene,
                "hgvsp": variant.hgvsp,
                "source": up.get("source", "failed"),
                "confidence": up.get("confidence", "low"),
                "interpro_id": None,
                "interpro_url": None,
            }
        
        return {
            "domain": "unknown",
            "note": "No UniProt domain annotation available for this gene.",
            "gene": gene,
            "hgvsp": variant.hgvsp,
            "source": up.get("source", "failed"),
            "confidence": up.get("confidence", "low"),
            "interpro_id": None,
            "interpro_url": None,
        }
    
    if not aa_pos:
        return {
            "domain": "unknown",
            "note": "Could not parse protein position from hgvsp",
            "gene": gene,
            "hgvsp": variant.hgvsp,
            "source": up.get("source", "unknown"),
            "interpro_id": None,
            "interpro_url": None,
        }
    
    # Search for domain overlap
    for d in domains:
        start = d.get("start")
        end = d.get("end")
        if start is None or end is None:
            continue
        
        if start <= aa_pos <= end:
            # Determine damage type
            if "frameshift" in variant.consequence or "nonsense" in variant.consequence:
                relative_pos = (aa_pos - start) / (end - start + 1)
                if relative_pos < 0.3:
                    damage = "N-terminal destruction, domain completely lost"
                    integrity = "completely_destroyed"
                elif relative_pos < 0.7:
                    damage = "Mid-domain truncation"
                    integrity = "partially_destroyed"
                else:
                    damage = "C-terminal truncation, partial domain retention"
                    integrity = "partially_retained"
            elif "splice" in variant.consequence:
                damage = "Splice site disruption, may cause exon skipping"
                integrity = "splicing_disrupted"
            else:
                damage = "Point mutation, assess specific amino acid change"
                integrity = "point_mutation"
            
            return {
                "domain": d.get("name", "unnamed"),
                "domain_range": f"aa{start}-{end}",
                "position_in_domain": f"aa{aa_pos}",
                "relative_position": round(relative_pos if "relative_pos" in locals() else (aa_pos - start) / (end - start + 1), 2),
                "function": d.get("type", "unknown"),
                "damage_type": damage,
                "domain_integrity": integrity,
                "gene": gene,
                "source": up.get("source", "unknown"),
                "confidence": up.get("confidence", "medium"),
                "interpro_id": (up.get("interpro_ids") or [None])[0],
                "interpro_url": f"https://www.ebi.ac.uk/interpro/entry/InterPro/{(up.get('interpro_ids') or [None])[0]}/" if (up.get("interpro_ids") or [None])[0] else None,
            }
    
    # Position outside all annotated domains
    return {
        "domain": "inter-domain / unannotated",
        "position": f"aa{aa_pos}",
        "note": "Position outside known functional domains",
        "gene": gene,
        "source": up.get("source", "unknown"),
        "confidence": up.get("confidence", "medium"),
        "interpro_id": (up.get("interpro_ids") or [None])[0],
        "interpro_url": f"https://www.ebi.ac.uk/interpro/entry/InterPro/{(up.get('interpro_ids') or [None])[0]}/" if (up.get("interpro_ids") or [None])[0] else None,
    }

# =============================================================================
# Module E: Tissue Relevance Assessment  (v0.4: GTEx API + fallback)
# =============================================================================

def assess_tissue_relevance(variant: Variant, tissue_profile: Dict,
                            gtex_data: Dict[str, Dict]) -> Dict:
    """
    Assess if gene is relevant to target tissue/organ context.
    
    Priority:
    1. GTEx API expression data (if available) — auto-classify by TPM thresholds
    2. Local tissue_context.json profile (fallback for API failures)
    3. Unknown if neither available — conservative, do NOT fast-track
    
    Args:
        variant: Variant object
        tissue_profile: Loaded tissue profile (tier_rules + special_gene_lists)
        gtex_data: Pre-fetched GTEx data {gene: {median_tpm, ...}}
    """
    gene = variant.gene
    profile_name = tissue_profile.get("display_name", "target tissue")
    
    # --- Priority 1: GTEx API data ---
    gtex = gtex_data.get(gene, {})
    tpm = gtex.get("median_tpm")
    
    if tpm is not None:
        # Auto-classify based on expression + ClinVar status
        if tpm >= 10.0:
            relevance = "primary"
            rationale = f"High {profile_name} expression (TPM={tpm:.1f}) per GTEx."
        elif tpm >= 1.0:
            relevance = "secondary"
            rationale = f"Moderate {profile_name} expression (TPM={tpm:.1f}) per GTEx."
        elif tpm > 0:
            relevance = "none"
            rationale = f"Low {profile_name} expression (TPM={tpm:.2f}) per GTEx."
        else:
            relevance = "none"
            rationale = f"No detectable {profile_name} expression per GTEx."
        
        # Fast track for none + benign
        if relevance == "none":
            if variant.clinvar and "Pathogenic" in variant.clinvar:
                return {
                    "tier_suggestion": 2,
                    "relevance": relevance,
                    "reason": f"{gene} is not {profile_name}-relevant (GTEx TPM={tpm:.2f}) but ClinVar pathogenic.",
                    "clinical_note": "Inform donor, record in medical history. Does not affect decision for this context.",
                    "fast_track": False,
                    "rationale": rationale,
                    "gtex_tpm": tpm,
                    "source": gtex.get("source", "gtex"),
                }
            
            return {
                "tier_suggestion": 3,
                "relevance": relevance,
                "reason": f"{gene} has no {profile_name} relevance (GTEx TPM={tpm:.2f}).",
                "clinical_note": f"No impact on {profile_name} function or safety.",
                "fast_track": True,
                "rationale": rationale,
                "gtex_tpm": tpm,
                "source": gtex.get("source", "gtex"),
            }
        
        # Primary or secondary: standard pipeline
        return {
            "tier_suggestion": "assess_via_standard_pipeline",
            "relevance": relevance,
            "gtex_tpm": tpm,
            "rationale": rationale,
            "fast_track": False,
            "clinical_note": f"{gene} is {relevance}-relevant to {profile_name} (GTEx TPM={tpm:.1f}).",
            "source": gtex.get("source", "gtex"),
        }
    
    # --- Priority 2: Local tissue_context.json fallback ---
    genes_local = tissue_profile.get("genes", {})
    relevance_info = genes_local.get(gene)
    
    if relevance_info:
        relevance = relevance_info.get("relevance", "unknown")
        gtex_rpkm_local = relevance_info.get("gtex_rpkm", None)
        rationale = relevance_info.get("rationale", "")
        
        if relevance == "none":
            if variant.clinvar and "Pathogenic" in variant.clinvar:
                return {
                    "tier_suggestion": 2,
                    "relevance": relevance,
                    "reason": f"{gene} is not {profile_name}-relevant but ClinVar pathogenic.",
                    "clinical_note": "Inform donor, record in medical history.",
                    "fast_track": False,
                    "rationale": rationale,
                    "gtex_rpkm": gtex_rpkm_local,
                    "source": "local_fallback",
                }
            
            return {
                "tier_suggestion": 3,
                "relevance": relevance,
                "reason": f"{gene} has no {profile_name} relevance.",
                "clinical_note": f"No impact on {profile_name} function or safety.",
                "fast_track": True,
                "rationale": rationale,
                "gtex_rpkm": gtex_rpkm_local,
                "source": "local_fallback",
            }
        
        return {
            "tier_suggestion": "assess_via_standard_pipeline",
            "relevance": relevance,
            "gtex_rpkm": gtex_rpkm_local,
            "rationale": rationale,
            "fast_track": False,
            "clinical_note": f"{gene} is {relevance}-relevant to {profile_name}.",
            "source": "local_fallback",
        }
    
    # --- Priority 3: Special gene lists (irreplaceable clinical rules) ---
    special_lists = tissue_profile.get("special_gene_lists", {})
    _SPECIAL_LIST_RELEVANCE = {
        # hematopoietic
        "coagulation": "primary",
        "fa_dna_repair": "primary",
        "drug_metabolism": "secondary",
        "kir_cluster": "secondary",
        # cardiovascular
        "cardiomyopathy": "primary",
        "ion_channel": "primary",
        "aortopathy": "primary",
        "arrhythmia": "primary",
        # hepatic
        "bilirubin_metabolism": "secondary",
        "cyp450": "secondary",
        "cholestatic": "primary",
        "hemochromatosis": "primary",
        # renal
        "renal_ciliopathy": "primary",
        "tubulopathy": "primary",
        "nephrotic": "primary",
        "congenital_abnormalities": "primary",
        # neurological
        "neurodegeneration": "primary",
        "leukodystrophy": "primary",
        "epilepsy": "primary",
        "movement_disorder": "primary",
    }
    for list_name, genes_in_list in special_lists.items():
        if gene in genes_in_list:
            relevance = _SPECIAL_LIST_RELEVANCE.get(list_name, "secondary")
            rationale = f"{gene} is in '{list_name}' special list for {profile_name}."
            
            return {
                "tier_suggestion": "assess_via_standard_pipeline",
                "relevance": relevance,
                "rationale": rationale,
                "fast_track": False,
                "action": "Proceed with standard domain + ClinVar + gnomAD assessment.",
                "clinical_note": f"{gene} is {relevance}-relevant to {profile_name} ({list_name} list).",
                "gtex_rpkm": None,
                "source": f"special_list:{list_name}",
            }
    
    # --- Priority 4: Completely unknown ---
    return {
        "tier_suggestion": "assess_via_standard_pipeline",
        "relevance": "unknown",
        "rationale": f"{gene} not in GTEx or local tissue profile '{profile_name}'.",
        "fast_track": False,
        "action": "Proceed with standard domain + ClinVar + gnomAD assessment. Do NOT fast-track.",
        "clinical_note": "Gene relevance unknown — conservative assessment.",
        "gtex_rpkm": None,
        "source": "unknown",
    }

# =============================================================================
# Module F: Three-Tier Risk Classification
# =============================================================================

TIER1_ACTION_GENES = {
    "VWF": {"reason": "Coagulation disorder affects collection safety", "condition": "ClinVar_Pathogenic"},
}

def classify_variant_tier(variant: Variant, domain_info: Dict, tissue_assessment: Dict,
                          gnomad_info: Dict, transcript_warning: Optional[Dict],
                          pseudogene_warning: Optional[Dict], tissue_profile: Dict) -> Tuple[int, str, List[str]]:
    """
    Three-tier classification with dynamic tissue context.
    Returns: (tier, reason, actions)
    """
    gene = variant.gene
    actions = []
    profile_name = tissue_profile.get("display_name", "target tissue")
    tier_rules = tissue_profile.get("tier_rules", {})
    special_lists = tissue_profile.get("special_gene_lists", {})

    # Priority 0: Fast track for non-relevant tissue genes
    if tissue_assessment.get("fast_track") and tissue_assessment.get("tier_suggestion") == 3:
        if not (variant.clinvar and "Pathogenic" in variant.clinvar):
            return 3, tissue_assessment["reason"], []

    # Priority 1: Tier 1 checks
    # 1a. Known high-risk special gene lists with pathogenic variant
    for list_name, gene_list in special_lists.items():
        if gene in gene_list:
            if "coagulation" in list_name.lower() and variant.clinvar and "Pathogenic" in variant.clinvar:
                actions.append("Assess bleeding history and coagulation function before collection")
                actions.append("Consider PBSC over BM to minimize bleeding risk")
                return 1, f"{gene} pathogenic variant affects collection safety (coagulation gene)", actions
            if "fa_dna_repair" in list_name.lower() and variant.clinvar and "Pathogenic" in variant.clinvar:
                actions.append("Assess if donor has Fanconi anemia phenotype")
                actions.append("Biallelic = ineligible donor; heterozygous = acceptable but monitor engraftment")
                return 1, f"{gene} pathogenic variant in FA pathway - marrow failure risk", actions

    # 1b. Homozygous truncating in primary tissue gene
    if variant.gt in ["1/1", "1|1"] and variant.impact == "HIGH":
        if tissue_assessment.get("relevance") == "primary":
            actions.append("Confirm homozygosity via secondary method")
            actions.append("Assess if phenotype is consistent with expected tissue function")
            return 1, f"Homozygous truncating variant in primary tissue gene {gene}", actions

    # Priority 2: Tier 2 checks
    # 2a. Primary tissue gene, heterozygous, function affected
    if tissue_assessment.get("relevance") in ["primary", "secondary"] and variant.gt in ["0/1", "0|1"]:
        if variant.impact == "HIGH":
            reason = f"Heterozygous {variant.consequence} in tissue-relevant gene {gene}"
            if domain_info and domain_info.get("domain_integrity") in ["completely_destroyed", "partially_destroyed"]:
                reason += f", {domain_info['domain']} domain disrupted"
            actions.append("Inform donor of carrier status")
            actions.append("Monitor post-intervention recovery/function")
            return 2, reason, actions

    # 2b. Non-primary but ClinVar pathogenic
    if variant.clinvar and "Pathogenic" in variant.clinvar and tissue_assessment.get("relevance") == "none":
        actions.append("Inform donor of genetic finding")
        actions.append("Refer for relevant specialist evaluation if indicated")
        return 2, f"ClinVar pathogenic variant in {gene} - donor's own health may be affected", actions

    # 2c. Drug metabolism genes (if applicable to this tissue context)
    drug_genes = special_lists.get("drug_metabolism", [])
    if gene in drug_genes:
        actions.append(f"Monitor post-intervention drug levels if relevant medications used")
        return 2, f"Drug metabolism variant may affect pharmacokinetics", actions

    # Priority 3: Tier 3 - everything else
    reason_parts = []
    if gnomad_info.get("status") == "common_polymorphism":
        reason_parts.append(f"Common polymorphism (AF={gnomad_info.get('af')})")
    if variant.clinvar and "Benign" in variant.clinvar:
        reason_parts.append("ClinVar benign")
    if tissue_assessment.get("relevance") == "none":
        reason_parts.append("No tissue relevance")

    reason = "; ".join(reason_parts) if reason_parts else "Low risk based on combined assessment"
    return 3, reason, []

# =============================================================================
# Multi-hit Gene Detection
# =============================================================================

def _variant_has_pathogenic_evidence(v: Variant, gtex_data: Optional[Dict] = None) -> bool:
    """
    Check if a variant has evidence suggesting pathogenicity.
    
    Criteria (OR):
      1. Affects protein domain (has specific domain mapping, not unknown/inter-domain)
         AND gene is expressed in target tissue (GTEx TPM >= 1.0)
      2. ClinVar pathogenic/likely pathogenic or HIGH impact or rare gnomAD (<0.001)
      3. Splice site change (consequence contains 'splice')
    
    Note: Domain impact in non-expressed genes (TPM < 1.0) is irrelevant for
    the target tissue context.
    """
    # Quick check: splice site changes are always considered
    consequence = str(v.consequence or "").lower()
    if "splice" in consequence:
        return True
    
    # Check tissue expression for domain relevance
    # If gene is not expressed in target tissue (TPM < 1.0), domain mapping
    # does not constitute pathogenic evidence for this context
    tissue_tpm = None
    if gtex_data and v.gene in gtex_data:
        tissue_tpm = gtex_data[v.gene].get("median_tpm")
    
    # Check ClinVar status for benign exclusion
    clinvar = str(v.clinvar or "").lower()
    is_benign = "benign" in clinvar and "conflicting" not in clinvar
    
    # 1. Domain impact — only counts if gene is expressed in target tissue
    # AND variant is not already classified as benign
    di = v.domain_info
    if di and not is_benign:
        domain = di.get("domain", "")
        if domain and domain not in ("unknown", "N/A", "inter-domain / unannotated"):
            # Domain impact is only relevant if tissue-expressed
            if tissue_tpm is not None and tissue_tpm < 1.0:
                pass  # Low expression: domain not relevant for this tissue
            else:
                return True
    
    # 2. Pathogenic evidence — always counts regardless of tissue expression
    if "pathogenic" in clinvar and "conflicting" not in clinvar:
        return True
    if v.impact == "HIGH":
        return True
    if v.gnomad_af is not None and v.gnomad_af < 0.001:
        return True
    
    # 3. Splice site changes — always considered
    consequence = str(v.consequence or "").lower()
    if "splice" in consequence:
        return True
    
    return False


def detect_multi_hit_genes(variants: List[Variant], gtex_data: Optional[Dict] = None) -> List[Dict]:
    """
    Detect genes with multiple pathogenic variants that may require phase analysis.
    
    Only counts variants with evidence of pathogenicity:
      - Domain impact, or
      - ClinVar pathogenic / HIGH impact / rare gnomAD, or  
      - Splice site change
    
    Normal polymorphisms (benign / common / no domain impact) are excluded.
    """
    # Group variants by gene
    gene_variants = {}
    for v in variants:
        gene_variants.setdefault(v.gene, []).append(v)
    
    multi_hits = []
    for gene, var_list in gene_variants.items():
        # Count only variants with pathogenic evidence
        pathogenic_vars = [v for v in var_list if _variant_has_pathogenic_evidence(v, gtex_data)]
        
        if len(pathogenic_vars) >= 2:
            # Collect details for each pathogenic variant
            var_details = []
            for v in pathogenic_vars:
                detail = {
                    "hgvsp": v.hgvsp,
                    "hgvsc": v.hgvsc,
                    "chrom": v.chrom,
                    "pos": v.pos,
                    "impact": v.impact,
                    "clinvar": v.clinvar,
                    "gnomad_af": v.gnomad_af,
                    "consequence": v.consequence,
                }
                if v.domain_info:
                    detail["domain"] = v.domain_info.get("domain")
                    detail["domain_range"] = v.domain_info.get("domain_range")
                var_details.append(detail)
            
            multi_hits.append({
                "gene": gene,
                "variant_count": len(var_list),           # total variants in gene
                "pathogenic_count": len(pathogenic_vars),  # variants with evidence
                "warning": "MULTI_HIT_GENE",
                "pathogenic_variants": var_details,
                "phases": {
                    "cis": "Both variants on same allele → other allele normal → heterozygous function retained",
                    "trans": "Variants on different alleles → compound heterozygous → function may be severely impaired"
                },
                "required_evidence": [
                    "Trio/family segregation analysis",
                    "Long-read sequencing (PacBio/Nanopore)",
                    "Allele-specific expression analysis"
                ],
                "action": "Priority P1: Must confirm phase before final assessment",
                "impact": "If trans: may elevate to Tier 1 regardless of individual variant assessment"
            })

    return multi_hits

# =============================================================================
# Patient-Donor Cross-check
# =============================================================================

def cross_check_patient_donor(patient_mutations: List[Dict], donor_variants: List[Variant]) -> List[Dict]:
    """
    Check if patient somatic driver mutations are present in donor germline.
    """
    results = []
    donor_genes = {v.gene: v for v in donor_variants}

    for pm in patient_mutations:
        gene = pm["gene"]
        if gene in donor_genes:
            dv = donor_genes[gene]
            results.append({
                "patient_mutation": pm,
                "donor_status": "PRESENT",
                "donor_variant": asdict(dv),
                "clinical_significance": f"HIGH RISK: Patient's somatic driver mutation is in donor germline. "
                                        f"If inherited, transplant may reintroduce leukemic clone.",
                "action": "URGENT: Verify if exact same variant. If yes, consider alternative donor."
            })
        else:
            results.append({
                "patient_mutation": pm,
                "donor_status": "NOT_PRESENT",
                "clinical_significance": "Favorable: Donor did not inherit this driver mutation.",
                "action": "No special action needed for this mutation."
            })

    return results

# =============================================================================
# Report Generation
# =============================================================================

def generate_tier_report(variants: List[Variant], config: DGRAConfig,
                        tissue_profile: Dict, multi_hits: List[Dict],
                        cross_check: List[Dict]) -> str:
    """
    Generate Markdown report with three-tier structure and dynamic tissue context.
    """
    # Sort by tier
    tier1 = [v for v in variants if v.tier == 1]
    tier2 = [v for v in variants if v.tier == 2]
    tier3 = [v for v in variants if v.tier == 3]

    profile_name = tissue_profile.get("display_name", config.tissue_profile)

    report = []
    report.append("# DGRA Report - Donor Genomic Risk Assessment v0.4\n")
    report.append(f"**Analysis Context**: {profile_name}\n")
    report.append(f"**Tissue Profile**: `{config.tissue_profile}`\n")
    report.append(f"**Offline Mode**: {'Yes' if config.offline_mode else 'No'}\n")
    report.append(f"**Analysis Date**: {datetime.now().isoformat()}\n")
    report.append(f"**Total Variants Assessed**: {len(variants)}\n")
    report.append(f"**Tier 1 (Action Required)**: {len(tier1)}\n")
    report.append(f"**Tier 2 (Inform & Monitor)**: {len(tier2)}\n")
    report.append(f"**Tier 3 (No Concern)**: {len(tier3)}\n\n")

    # Multi-hit warnings
    if multi_hits:
        report.append("## ⚠️ Multi-Hit Gene Warnings\n")
        for mh in multi_hits:
            report.append(f"### {mh['gene']} - {mh['variant_count']} variants detected\n")
            report.append(f"- **Warning**: {mh['warning']}\n")
            report.append(f"- **Cis hypothesis**: {mh['phases']['cis']}\n")
            report.append(f"- **Trans hypothesis**: {mh['phases']['trans']}\n")
            report.append(f"- **Required evidence**: {', '.join(mh['required_evidence'])}\n")
            report.append(f"- **Action**: {mh['action']}\n\n")

    # Patient-donor cross-check
    if cross_check:
        report.append("## Patient-Donor Cross-Check\n")
        for cc in cross_check:
            status_icon = "🔴" if cc['donor_status'] == "PRESENT" else "🟢"
            report.append(f"### {status_icon} {cc['patient_mutation']['gene']}\n")
            report.append(f"- **Donor status**: {cc['donor_status']}\n")
            report.append(f"- **Significance**: {cc['clinical_significance']}\n")
            report.append(f"- **Action**: {cc['action']}\n\n")

    # Tier 1
    if tier1:
        report.append("---\n\n## 🔴 Tier 1: Action Required\n")
        report.append(f"*Variants requiring intervention for {profile_name} context*\n\n")
        
        # Group by gene
        from collections import OrderedDict
        gene_groups = OrderedDict()
        for v in tier1:
            gene_groups.setdefault(v.gene, []).append(v)
        
        for gene, var_list in gene_groups.items():
            report.append(f"### {gene}\n")
            
            # Gene summary
            report.append(f"**基因**: {gene} | **变异数**: {len(var_list)}\n\n")
            
            # Variant table
            report.append("| # | 染色体位置 | 转录本 | 变异名称 | 功能域 | 合子型 | ClinVar | 说明 |\n")
            report.append("|---|-----------|--------|---------|--------|--------|---------|------|\n")
            
            for i, v in enumerate(var_list, 1):
                # Position
                pos = f"{v.chrom}:{v.pos}"
                
                # Transcript
                tx = v.transcript or "N/A"
                
                # Variant name
                var_name = v.hgvsp or v.hgvsc or "N/A"
                
                # Domain
                di = v.domain_info
                if di:
                    domain = f"{di.get('domain', 'N/A')} ({di.get('domain_range', 'N/A')})"
                else:
                    domain = "N/A"
                
                # Zygosity
                zyg = v.gt or "N/A"
                
                # ClinVar
                clin = v.clinvar or "N/A"
                
                # Reason (shortened)
                reason = v.tier_reason[:80] + "..." if len(v.tier_reason) > 80 else v.tier_reason
                reason = reason.replace("|", "/")  # avoid markdown table break
                
                report.append(f"| {i} | {pos} | {tx} | {var_name} | {domain} | {zyg} | {clin} | {reason} |\n")
            
            report.append(f"\n**详细说明**:\n")
            for i, v in enumerate(var_list, 1):
                report.append(f"{i}. **{v.hgvsp or v.hgvsc}** ({v.chrom}:{v.pos}):\n")
                report.append(f"   - 影响程度: {v.impact} | 后果: {v.consequence}\n")
                if v.domain_info:
                    di = v.domain_info
                    report.append(f"   - 功能域: {di.get('domain', 'N/A')} {di.get('domain_range', 'N/A')}\n")
                rp = di.get('relative_position')
                if isinstance(rp, (int, float)):
                    report.append(f"   - 域内位置: {di.get('position_in_domain', 'N/A')} (相对: {rp:.2f})\n")
                else:
                    report.append(f"   - 域内位置: {di.get('position_in_domain', 'N/A')} (相对: {rp})\n")
                    report.append(f"   - 损伤评估: {di.get('damage_type', 'N/A')}\n")
                if v.tissue_relevance:
                    tr = v.tissue_relevance
                    report.append(f"   - 组织相关性: {tr.get('relevance', 'N/A')} | GTEx TPM: {tr.get('gtex_tpm', 'N/A')}\n")
                report.append(f"   - 分级原因: {v.tier_reason}\n")
                if v.tier_actions:
                    report.append(f"   - 建议措施: {'; '.join(v.tier_actions)}\n")
                report.append(f"\n")

    # Tier 2
    if tier2:
        report.append("---\n\n## 🟡 Tier 2: Inform & Monitor\n")
        report.append(f"*Variants donors should be informed of for {profile_name} context*\n\n")
        
        # Group by gene
        from collections import OrderedDict
        gene_groups_t2 = OrderedDict()
        for v in tier2:
            gene_groups_t2.setdefault(v.gene, []).append(v)
        
        for gene, var_list in gene_groups_t2.items():
            report.append(f"### {gene}\n")
            report.append(f"**基因**: {gene} | **变异数**: {len(var_list)}\n\n")
            
            # Variant table
            report.append("| # | 染色体位置 | 转录本 | 变异名称 | 功能域 | 合子型 | ClinVar | 说明 |\n")
            report.append("|---|-----------|--------|---------|--------|--------|---------|------|\n")
            
            for i, v in enumerate(var_list, 1):
                pos = f"{v.chrom}:{v.pos}"
                tx = v.transcript or "N/A"
                var_name = v.hgvsp or v.hgvsc or "N/A"
                di = v.domain_info
                if di:
                    domain = f"{di.get('domain', 'N/A')} ({di.get('domain_range', 'N/A')})"
                else:
                    domain = "N/A"
                zyg = v.gt or "N/A"
                clin = v.clinvar or "N/A"
                reason = v.tier_reason[:80] + "..." if len(v.tier_reason) > 80 else v.tier_reason
                reason = reason.replace("|", "/")
                report.append(f"| {i} | {pos} | {tx} | {var_name} | {domain} | {zyg} | {clin} | {reason} |\n")
            
            report.append(f"\n**详细说明**:\n")
            for i, v in enumerate(var_list, 1):
                report.append(f"{i}. **{v.hgvsp or v.hgvsc}** ({v.chrom}:{v.pos}):\n")
                report.append(f"   - 影响程度: {v.impact} | 后果: {v.consequence}\n")
                if v.domain_info:
                    di = v.domain_info
                    report.append(f"   - 功能域: {di.get('domain', 'N/A')} {di.get('domain_range', 'N/A')}\n")
                    rp = di.get('relative_position')
                    if isinstance(rp, (int, float)):
                        report.append(f"   - 域内位置: {di.get('position_in_domain', 'N/A')} (相对: {rp:.2f})\n")
                    else:
                        report.append(f"   - 域内位置: {di.get('position_in_domain', 'N/A')} (相对: {rp})\n")
                    report.append(f"   - 损伤评估: {di.get('damage_type', 'N/A')}\n")
                if v.tissue_relevance:
                    tr = v.tissue_relevance
                    report.append(f"   - 组织相关性: {tr.get('relevance', 'N/A')} | GTEx TPM: {tr.get('gtex_tpm', 'N/A')}\n")
                report.append(f"   - 分级原因: {v.tier_reason}\n")
                if v.tier_actions:
                    report.append(f"   - 建议措施: {'; '.join(v.tier_actions)}\n")
                report.append(f"\n")

    # Tier 3
    if tier3:
        report.append("---\n\n## 🟢 Tier 3: No Concern\n")
        report.append(f"*Variants with no {profile_name} relevance*\n\n")
        
        # Group by gene
        gene_groups_t3 = {}
        for v in tier3:
            gene_groups_t3.setdefault(v.gene, []).append(v)
        
        for gene, var_list in gene_groups_t3.items():
            if len(var_list) <= 3:
                # Short list: inline table
                report.append(f"**{gene}** ({len(var_list)} variants):\n")
                report.append("| 位置 | 变异 | 功能域 | 原因 |\n")
                report.append("|------|------|--------|------|\n")
                for v in var_list:
                    pos = f"{v.chrom}:{v.pos}"
                    var_name = v.hgvsp or v.hgvsc or "N/A"
                    di = v.domain_info
                    domain = f"{di.get('domain', 'N/A')}" if di else "N/A"
                    reason = v.tier_reason[:50] + "..." if len(v.tier_reason) > 50 else v.tier_reason
                    reason = reason.replace("|", "/")
                    report.append(f"| {pos} | {var_name} | {domain} | {reason} |\n")
                report.append(f"\n")
            else:
                # Many variants: just count
                report.append(f"**{gene}**: {len(var_list)} variants — 详见原始数据\n\n")
        
        report.append(f"\n")

    # Methodology
    report.append("---\n\n## 方法学附录\n")
    report.append(f"### 组织背景: {profile_name}\n")
    report.append(f"- **GTEx 参考组织**: {tissue_profile.get('gtex_tissue', 'N/A')}\n")
    report.append(f"- **快速排除规则**: {tissue_profile.get('fast_track_rule', 'N/A')}\n\n")
    report.append("### 分析流程\n")
    report.append("1. **转录本校正**: Ensembl REST API → canonical transcript → 本地回退\n")
    report.append("2. **假基因检测**: VAF 偏差分析识别已知假基因对\n")
    report.append("3. **gnomAD 整合**: AF>1% 常见; AF<0.1% 罕见; NOT_CAPTURED 明确标注\n")
    report.append("4. **蛋白功能域映射**: UniProt REST API → DOMAIN/REGION 特征 → 本地回退\n")
    report.append("5. **组织相关性评估**: GTEx API → median TPM → 自动分级 + 本地回退\n")
    report.append("6. **三级分类**: Action (Tier 1) → Inform (Tier 2) → No concern (Tier 3)\n")
    report.append("7. **患者-供者交叉核对**: 体细胞驱动突变的遗传检测\n")
    report.append("8. **缓存**: 所有 API 响应缓存 30 天 (SQLite); 离线模式仅用缓存\n")

    return "\n".join(report)

# =============================================================================
# Main Pipeline  (v0.4: async + batch API queries)
# =============================================================================

async def run_dgra_pipeline(variants_data: List[Dict], patient_mutations: List[Dict] = None,
                      config: Optional[DGRAConfig] = None) -> Dict:
    """
    Main DGRA analysis pipeline with dynamic tissue context.
    v0.4: Async, batch API queries with cache.

    Args:
        variants_data: List of variant dicts from VCF annotation
        patient_mutations: List of patient somatic driver mutations
        config: DGRA configuration (includes tissue_profile + offline_mode)

    Returns:
        Dict with report and structured results
    """
    if config is None:
        config = DGRAConfig()

    # Convert user config to global config
    global_config = config.to_global()

    # Load tissue profile (keeps tier_rules + special_gene_lists)
    tissue_profile = config.get_tissue_profile()
    profile_name = tissue_profile.get("display_name", config.tissue_profile)

    # Parse variants
    variants = []
    for vd in variants_data:
        v = Variant(
            chrom=vd.get("CHROM", ""),
            pos=int(vd.get("POS", 0)),
            ref=vd.get("REF", ""),
            alt=vd.get("ALT", ""),
            gene=vd.get("GENE", ""),
            transcript=vd.get("Feature", ""),
            exon=vd.get("EXON", ""),
            impact=vd.get("IMPACT", ""),
            consequence=vd.get("Consequence", ""),
            hgvsp=vd.get("HGVSp", ""),
            hgvsc=vd.get("HGVSc", ""),
            clinvar=vd.get("CLIN_SIG", ""),
            dp=int(vd.get("DP", 0)) if vd.get("DP") else 0,
            gq=float(vd.get("GQ", 0)) if vd.get("GQ") else 0,
            gt=vd.get("GT", ""),
            vaf=float(vd.get("VAF", 0)) if vd.get("VAF") else None,
            gnomad_af=float(vd.get("gnomAD_AF", 0)) if vd.get("gnomAD_AF") else None,
        )
        variants.append(v)

    # Collect all unique genes for batch API queries
    unique_genes = list({v.gene for v in variants})
    print(f"[DGRA] {len(variants)} variants across {len(unique_genes)} unique genes")
    print(f"[DGRA] Tissue profile: {profile_name} | Offline: {config.offline_mode}")

    # ------------------------------------------------------------------
    # Batch API queries (concurrent)
    # ------------------------------------------------------------------
    ensembl_data = {}
    uniprot_data = {}
    gtex_data = {}

    if not config.offline_mode and unique_genes:
        cache = DGRACache(global_config.cache_db_path)
        async with DGRAAPIClient(global_config, cache) as client:
            # Batch query all three APIs concurrently
            ensembl_raw, uniprot_raw, gtex_raw = await asyncio.gather(
                client.batch_query_genes(unique_genes, "ensembl"),
                client.batch_query_genes(unique_genes, "uniprot"),
                client.batch_query_genes(
                    unique_genes, "gtex",
                    tissue=tissue_profile.get("gtex_tissue", "Whole Blood")
                ),
            )
            ensembl_data = {g: ensembl_raw.get(g, {}) for g in unique_genes}
            uniprot_data = {g: uniprot_raw.get(g, {}) for g in unique_genes}
            gtex_data = {g: gtex_raw.get(g, {}) for g in unique_genes}
        print(f"[DGRA] API batch query complete: Ensembl={len(ensembl_data)}, UniProt={len(uniprot_data)}, GTEx={len(gtex_data)}")
        # Persist successful API results for future offline use
        for gene in unique_genes:
            _save_offline_archive(gene, ensembl_data, uniprot_data, gtex_data, config.tissue_profile)
        print(f"[DGRA] Offline archive saved for {len(unique_genes)} genes to {OFFLINE_ARCHIVE_DIR}")
    else:
        # Offline mode: try to load archived data first, then fall back to local rules
        loaded = 0
        for gene in unique_genes:
            archive = _load_offline_archive(gene)
            if archive:
                ensembl_data[gene] = archive.get("ensembl", {})
                uniprot_data[gene] = archive.get("uniprot", {})
                gtex_data[gene] = archive.get("gtex", {})
                loaded += 1
        print(f"[DGRA] Offline mode: loaded archived data for {loaded}/{len(unique_genes)} genes from {OFFLINE_ARCHIVE_DIR}")
        if loaded == 0:
            print("[DGRA] Offline mode: no archive found, using local fallbacks only (conservative)")

    # ------------------------------------------------------------------
    # Step 1: Transcript correction (v0.4: Ensembl API)
    # ------------------------------------------------------------------
    for v in variants:
        v, warning = await correct_transcript_priority(v, ensembl_data)
        if warning:
            v.transcript_warning = json.dumps(warning)

    # Step 2: Pseudogene detection (unchanged)
    for v in variants:
        pg_warning = detect_pseudogene_artifact(v)
        if pg_warning:
            v.pseudogene_warning = json.dumps(pg_warning)

    # Step 3: gnomAD classification (unchanged)
    for v in variants:
        gnomad_info = classify_gnomad_frequency(v.gnomad_af, v.gene)
        v.gnomad_status = gnomad_info["status"]

    # Step 4: Protein domain mapping (v0.4: UniProt API)
    for v in variants:
        v.domain_info = map_variant_to_domain(v, uniprot_data)

    # Step 5: Tissue relevance assessment (v0.4: GTEx API + local fallback)
    tissue_assessments = {}
    for v in variants:
        tissue = assess_tissue_relevance(v, tissue_profile, gtex_data)
        tissue_assessments[v.gene] = tissue
        v.tissue_relevance = tissue

    # Step 6: Multi-hit detection (unchanged)
    multi_hits = detect_multi_hit_genes(variants, gtex_data)

    # Step 7: Three-tier classification (with tissue context)
    for v in variants:
        tissue = tissue_assessments[v.gene]
        gnomad_info = classify_gnomad_frequency(v.gnomad_af, v.gene)
        tw = json.loads(v.transcript_warning) if v.transcript_warning else None
        pw = json.loads(v.pseudogene_warning) if v.pseudogene_warning else None

        tier, reason, actions = classify_variant_tier(
            v, v.domain_info, tissue, gnomad_info, tw, pw, tissue_profile
        )
        v.tier = tier
        v.tier_reason = reason
        v.tier_actions = actions

    # Handle multi-hit elevation
    # Phase 1 v0.4.3: Exclude HLA genes from elevation — high polymorphism is normal
    _HLA_GENES = {
        "HLA-A", "HLA-B", "HLA-C", "HLA-DRB1", "HLA-DQA1", "HLA-DQB1",
        "HLA-DPA1", "HLA-DPB1", "HLA-E", "HLA-F", "HLA-G", "HLA-H",
        "HLA-J", "HLA-K", "HLA-L", "HLA-N", "HLA-P", "HLA-S",
        "HLA-DMA", "HLA-DMB", "HLA-DOA", "HLA-DOB",
        "MICA", "MICB", "TAP1", "TAP2",
    }
    multi_hit_genes = {mh["gene"] for mh in multi_hits}
    # Filter out HLA genes: natural polymorphism, not pathogenic multi-hit
    hla_multi_hits = {g for g in multi_hit_genes if g in _HLA_GENES}
    multi_hit_genes -= hla_multi_hits  # Exclude from elevation
    
    # v0.4.4: Only elevate variants that themselves have pathogenic evidence.
    # A gene may have 6 variants but only 3 pass _variant_has_pathogenic_evidence;
    # only those 3 should be elevated, not all 6.
    for v in variants:
        if v.gene in multi_hit_genes and v.tier > 1:
            if _variant_has_pathogenic_evidence(v, gtex_data):
                v.tier = 1
                v.tier_reason += " | ELEVATED: Multi-hit gene, phase unknown. Must confirm cis/trans."
                actions = v.tier_actions or []
                actions.append("URGENT: Confirm phase before final assessment")
                v.tier_actions = actions
            # Variants without pathogenic evidence remain at their original tier
            v.tier_reason += " | ELEVATED: Multi-hit gene, phase unknown. Must confirm cis/trans."
            v.tier_actions.append("URGENT: Confirm phase before final assessment")
    
    # Note: HLA multi-hits are logged but NOT elevated — they are normal polymorphism
    # The hla_multi_hits set can be used for reporting if needed

    # Step 8: Patient-donor cross-check (unchanged)
    cross_check = []
    if patient_mutations:
        cross_check = cross_check_patient_donor(patient_mutations, variants)

    # Step 9: Generate report (with tissue context)
    report_md = generate_tier_report(variants, config, tissue_profile, multi_hits, cross_check)

    # Compile structured output
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
            "tier1_count": len([v for v in variants if v.tier == 1]),
            "tier2_count": len([v for v in variants if v.tier == 2]),
            "tier3_count": len([v for v in variants if v.tier == 3]),
            "multi_hit_genes": [mh["gene"] for mh in multi_hits],
            "patient_inherited_mutations": [cc["patient_mutation"]["gene"] for cc in cross_check if cc["donor_status"] == "PRESENT"]
        },
        "tier1_variants": [asdict(v) for v in variants if v.tier == 1],
        "tier2_variants": [asdict(v) for v in variants if v.tier == 2],
        "tier3_variants": [asdict(v) for v in variants if v.tier == 3],
        "multi_hit_details": multi_hits,
        "patient_donor_cross_check": cross_check,
        "report_markdown": report_md
    }

    return output

# =============================================================================
# CLI Interface  (v0.4: --offline + asyncio.run)
# =============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="DGRA - Donor Genomic Risk Assessment (v0.4 API-first with cache)",
        epilog="Available tissue profiles: hematopoietic (default), cardiovascular, hepatic, renal, neurological. "
               "Define custom profiles in references/tissue_context.json"
    )
    parser.add_argument("--input", "-i", required=True, help="Input CSV/TSV with annotated variants")
    parser.add_argument("--patient-mutations", "-p", help="JSON file with patient somatic mutations")
    parser.add_argument("--output", "-o", default="dgra_report.md", help="Output report file")
    parser.add_argument("--json", "-j", help="Output JSON file with structured results")
    parser.add_argument("--tissue", "-t", required=True,
                        help="Tissue/organ context profile (REQUIRED). "
                             "Controls which genes are considered relevant for tier classification. "
                             "Available: hematopoietic, cardiovascular, hepatic, renal, neurological")
    parser.add_argument("--offline", action="store_true",
                        help="Offline mode: skip all API calls, use cache + local references only")

    args = parser.parse_args()

    # Read variants from CSV/TSV
    variants_data = []
    with open(args.input, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t' if args.input.endswith('.tsv') else ',')
        for row in reader:
            variants_data.append(dict(row))

    # Read patient mutations if provided
    patient_mutations = None
    if args.patient_mutations:
        with open(args.patient_mutations, 'r') as f:
            patient_mutations = json.load(f)

    # Run async pipeline with tissue context
    config = DGRAConfig(tissue_profile=args.tissue, offline_mode=args.offline)
    results = asyncio.run(run_dgra_pipeline(variants_data, patient_mutations, config))

    # Write report
    with open(args.output, 'w') as f:
        f.write(results["report_markdown"])

    profile_name = results["meta"]["profile_display_name"]
    print(f"DGRA Report Generated: {args.output}")
    print(f"Tissue Context: {profile_name} ({args.tissue})")
    print(f"Summary: Tier 1={results['summary']['tier1_count']}, "
          f"Tier 2={results['summary']['tier2_count']}, "
          f"Tier 3={results['summary']['tier3_count']}")

    if results['summary']['multi_hit_genes']:
        print(f"Multi-hit genes: {', '.join(results['summary']['multi_hit_genes'])}")

    # Write JSON if requested
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"Structured output written to: {args.json}")

if __name__ == "__main__":
    main()
