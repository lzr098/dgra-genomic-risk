#!/usr/bin/env python3
"""
GPA Core Engine - Genomic Phenotype Association

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
import os
import asyncio
from datetime import datetime
from pathlib import Path
from dataclasses import dataclass, field, asdict
from typing import List, Dict, Optional, Tuple, Any
import argparse

from version import __version__

import aiohttp

# v0.8.0: SpliceAI requires aiohttp for async HTTP queries
try:
    import aiohttp
except ImportError:
    aiohttp = None

# v0.5 P2-3: YAML config support
from dgra_config import DGRAGlobalConfig as GPAGlobalConfig, DGRAFileConfig as GPAFileConfig, DEFAULT_CONFIG_PATH
from dgra_cache import DGRACache
from dgra_api import DGRAAPIClient
from dgra_splice_predictor import SpliceAIPredictor, should_query_spliceai

# =============================================================================
# Offline Archive Persistence (v0.4)
# =============================================================================

OFFLINE_ARCHIVE_DIR = Path(__file__).parent.parent / "references" / "offline_data"

def _save_offline_archive(gene: str, ensembl_data: dict, uniprot_data: dict,
                          gtex_data: dict, tissue_profile: str,
                          gnomad_constraint_data: dict = None):
    """Persist API results for future offline use.
    v0.5 P1-4: Now also saves gnomAD gene constraint data."""
    OFFLINE_ARCHIVE_DIR.mkdir(parents=True, exist_ok=True)
    archive = {
        "gene": gene,
        "tissue_profile": tissue_profile,
        "saved_at": datetime.utcnow().isoformat() + "Z",
        "ensembl": ensembl_data.get(gene, {}),
        "uniprot": uniprot_data.get(gene, {}),
        "gtex": gtex_data.get(gene, {}),
    }
    if gnomad_constraint_data is not None:
        archive["gnomad_constraint"] = gnomad_constraint_data.get(gene, {})
    path = OFFLINE_ARCHIVE_DIR / f"{gene}.json"
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(archive, f, indent=2, ensure_ascii=False, default=str)

def _load_offline_archive(gene: str) -> Optional[Dict]:
    """Load previously saved API results for offline mode.
    v0.5 P1-4: Returns gnomad_constraint field if present."""
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

# v0.5 P0-7: Sentinel for missing/unknown field values. Used instead of injecting
# false defaults like IMPACT="MODERATE" or VAF=0.5, which systematically
# underestimate risk.
_UNKNOWN = "UNKNOWN"


def _is_unknown(val) -> bool:
    """Check if a value represents an unknown/missing field.

    Handles the three sentinel values used across the codebase:
    - _UNKNOWN ("UNKNOWN") — explicit unknown sentinel
    - "" — empty string
    - None — Python None
    """
    return val == _UNKNOWN or val == "" or val is None


# v0.4.5: Common cancer gene lists for somatic mode
_COMMON_TS_GENES = {
    "TP53", "RB1", "CDKN2A", "CDKN2B", "PTEN", "NF1", "NF2", "APC", "BRCA1", "BRCA2",
    "ATM", "CHEK2", "MLH1", "MSH2", "MSH6", "PMS2", "VHL", "WT1", "TSC1", "TSC2",
    "PHF6", "BCOR", "BCORL1", "ASXL1", "RUNX1", "CEBPA", "GATA2", "ETV6", "DDX41",
    "SAMD9", "SAMD9L", "TP53BP1", "BRCC3", "RAD51", "RAD51C", "RAD51D", "PALB2",
    "BARD1", "NBN", "ATM", "CHEK2",
}

_KNOWN_AML_DRIVERS = {
    "FLT3", "NPM1", "IDH1", "IDH2", "DNMT3A", "TET2", "ASXL1", "RUNX1", "CEBPA",
    "TP53", "KIT", "NRAS", "KRAS", "PTPN11", "CBL", "JAK2", "JAK3", "SH2B3",
    "BCOR", "BCORL1", "PHF6", "GATA2", "ETV6", "DDX41", "SAMD9", "SAMD9L",
    "KDM6A", "KDM5C", "KMT2A", "KMT2D", "EZH2", "STAG2", "RAD21", "SMC1A",
    "SMC3", "ZRSR2", "SRSF2", "SF3B1", "U2AF1", "U2AF2",
}

# v0.5.1 OPT-P2-2: Gene family redundancy - reduce multi-hit false positives
# Genes with functional paralogs/isoforms that can compensate for LOF
_GENE_FAMILY_REDUNDANCY = {
    # Mitochondrial ADP/ATP translocases - 4 paralogs with overlapping function
    "SLC25A5": {
        "paralogs": ["SLC25A4", "SLC25A6", "SLC25A31"],  # ANT1, ANT3, ANT4
        "compensation_level": "partial",  # Not complete, but significant
        "reason": "ANT family has 4 paralogs; SLC25A5 (ANT2) loss partially compensated",
    },
    # CYP family - extensive redundancy for drug metabolism
    "CYP2D6": {
        "paralogs": ["CYP2C19", "CYP3A4", "CYP1A2", "CYP2C9"],
        "compensation_level": "partial",
        "reason": "CYP450 family redundancy; other isoforms handle most substrates",
    },
    # SLC drug transporters
    "SLC22A1": {
        "paralogs": ["SLC22A2", "SLC22A3"],  # OCT2, OCT3
        "compensation_level": "partial",
        "reason": "OCT family redundancy for cation transport",
    },
    # HLA class I - extensive polymorphism is normal, null alleles common
    "HLA-A": {"paralogs": [], "compensation_level": "complete", "reason": "HLA null alleles are normal polymorphism"},
    "HLA-B": {"paralogs": [], "compensation_level": "complete", "reason": "HLA null alleles are normal polymorphism"},
    "HLA-C": {"paralogs": [], "compensation_level": "complete", "reason": "HLA null alleles are normal polymorphism"},
}

@dataclass
class Evidence:
    """Structured evidence entry for variant tier classification.
    v0.5 P1-9: Replaces free-text tier_reason with traceable evidence chain.
    """
    source: str           # "ClinVar", "gnomAD", "GTEx", "ACMG", "GeneConstraint", "NMD", "DomainMapping", "Somatic", "SpecialList", "MultiHit", "Phase", "FastTrack"
    rule: str             # Human-readable rule/condition, e.g. "PM2: gnomAD_AF=0.0001 < 0.001"
    weight: float = 1.0   # Contribution weight to final tier (default 1.0)
    confidence: str = "high"  # "high" | "moderate" | "low"
    raw_data: Optional[Dict] = None  # Key API response fields for audit


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
    gnomad_populations: Optional[Dict[str, Dict]] = None  # v0.5 P1-1: per-pop AFs
    dp: int = 0
    gq: float = 0.0
    gt: str = ""  # 0/1, 1/1, etc.
    vaf: Optional[float] = None

    # v0.4.5: Somatic annotation fields
    classification: str = ""  # OncoKB: Oncogenic, Likely Oncogenic, VUS, etc.
    is_tsg: bool = False
    is_oncogene: bool = False

    # v0.9.1: gnomAD query status (hotfix DDX3X misclassification)
    gnomad_status: str = "UNKNOWN"  # SUCCESS | API_FAILED | NOT_CAPTURED | ERROR
    gnomad_error_msg: Optional[str] = None
    gnomad_af_warning: bool = False

    # Computed fields
    tier: Optional[int] = None
    tier_reason: str = ""
    tier_actions: List[str] = field(default_factory=list)
    domain_info: Optional[Dict] = None
    transcript_warning: Optional[str] = None
    pseudogene_warning: Optional[str] = None
    tissue_relevance: Optional[Dict] = None

    # v0.5 P0-7: Quality confidence when key fields are missing.
    quality_confidence: str = "high"  # "high" | "medium" | "low" | "unknown"
    missing_fields: List[str] = field(default_factory=list)

    # v0.5 P1-2: Original gene symbol before HGNC normalization
    gene_original: str = ""

    # v0.5 P1-9: Structured evidence chain for traceable tier classification
    evidence_chain: List[Evidence] = field(default_factory=list)

    # v0.5 P1-13: Quality control flags for input validation
    qc_flags: List[str] = field(default_factory=list)

    # v0.5 P1-11: Upgrade conditions - forward-looking evidence gap descriptions
    upgrade_conditions: List[str] = field(default_factory=list)

    # v0.5 P1-4: Gene constraint metrics (pLI, LOEUF, lof_z, mis_z)
    gene_constraint: Optional[Dict] = None

    # v0.5 P1-10: Tier confidence classification
    tier_confidence: str = "LOW"  # "HIGH" | "MEDIUM" | "LOW"

    # v0.7: Phenotype association fields
    phenotype_match_score: Optional[float] = None
    phenotype_match_explanation: str = ""
    phenotype_match_confidence: str = ""
    phenotype_matched_pairs: List = field(default_factory=list)
    phenotype_known_list: List[str] = field(default_factory=list)

    # v0.7.2: ClinVar review status (CLNREVSTAT)
    clinvar_review_status: Optional[str] = None
    # v0.8.0: SpliceAI lookup result (pre-computed in pipeline)
    spliceai_result: Optional[Dict[str, Any]] = None

    # v0.9.0: Transcript selection (disease-aware, multi-transcript保留)
    primary_transcript: Optional[str] = None
    primary_consequence: Optional[str] = None
    primary_hgvsc: Optional[str] = None
    primary_hgvsp: Optional[str] = None
    primary_impact: Optional[str] = None
    alternative_transcripts: List[Dict[str, Any]] = field(default_factory=list)
    transcript_selection_method: str = "canonical"  # canonical / tissue_expression / llm_disease_match / ambiguous / user_specified
    transcript_ambiguity_flag: bool = False
    transcript_selection_log: str = ""

    # VCF原始信息保留
    vcf_filter: Optional[str] = None
    vcf_info: Dict[str, Any] = field(default_factory=dict)
    vcf_format: Optional[str] = None
    vcf_sample: Dict[str, str] = field(default_factory=dict)

@dataclass
class GPAConfig:
    """User-facing config (kept simple). Maps to GPAGlobalConfig internally.
    v0.4-fix: tissue_profile has NO default - must be specified by user.
    """
    min_dp: int = 20
    min_gq: float = 90.0
    common_af_threshold: float = 0.01
    low_af_threshold: float = 0.001
    vaf_deviation_threshold: float = 0.20
    tissue_profile: Optional[str] = None  # NO default - caller must specify
    target_population: Optional[str] = None  # v0.5 P1-1: EAS, AMR, AFR, NFE, SAS, etc.
    offline_mode: bool = False
    somatic_mode: bool = False  # v0.4.5: tumor/somatic driver analysis mode
    multi_organ_profiles: Optional[List[str]] = None  # v0.5 P1-7: multi-organ assessment
    gene_sync_enabled: bool = True  # v0.5 P1-8: auto-sync special_gene_lists
    filter_stats: Optional[Dict[str, Any]] = None  # v0.7.1: variant pre-filtering statistics
    filter_preset: Optional[str] = None  # v0.7.1: filter preset name used
    gene_sync_ttl_days: int = 7  # v0.5 P1-8: sync cache TTL
    force_sync: bool = False  # v0.5 P1-8: force sync gene lists (bypass cache)
    evidence_detail: str = "brief"  # v0.5 P1-9: "brief" | "full" - evidence chain detail level in report
    database_version: Optional[str] = None  # v0.5 P1-15: freeze analysis DB version for reproducibility
    # v0.8.0: SpliceAI splice-prediction integration (default OFF - must be explicitly enabled)
    spliceai_enabled: bool = False
    spliceai_concurrency: int = 5
    # v0.9.0: VCF annotation + disease-aware transcript selection
    disease_description: Optional[str] = None
    annotator: str = "auto"
    vep_cache: Optional[str] = None
    # v0.10.1: Two-phase pipeline for large VCF optimization
    two_phase: bool = False

    def to_global(self) -> GPAGlobalConfig:
        gc = GPAGlobalConfig()
        gc.tissue_profile = self.tissue_profile or ""
        gc.target_population = self.target_population or ""
        gc.offline_mode = self.offline_mode
        gc.somatic_mode = self.somatic_mode
        gc.gene_sync_enabled = self.gene_sync_enabled
        gc.gene_sync_ttl_days = self.gene_sync_ttl_days
        gc.min_dp = self.min_dp
        gc.min_gq = self.min_gq
        gc.common_af_threshold = self.common_af_threshold
        gc.low_af_threshold = self.low_af_threshold
        gc.vaf_deviation_threshold = self.vaf_deviation_threshold
        # v0.10.3: Allow cache DB override for sandboxed environments
        import os
        cache_override = os.environ.get("DGRA_CACHE_DB_PATH")
        if cache_override:
            gc.cache_db_path = Path(cache_override)
        return gc

    def get_tissue_profile(self, force_sync: bool = False) -> Dict:
        """Load tissue profile from references/tissue_context.json.
        v0.5 P1-8: Merged with external sync + user extensions for special_gene_lists.

        Priority of special_gene_lists:
          1. Hardcoded CORE (safety-critical, immutable)
          2. User add/remove extensions (user_gene_lists.json)
          3. Auto-synced from Orphanet / OMIM
          4. Static JSON from tissue_context.json

        Raises if tissue_profile is not set.
        Args:
            force_sync: Force a fresh sync even if cache is not expired.
        """
        if not self.tissue_profile:
            raise ValueError(
                "tissue_profile is required. Available profiles: general, hematopoietic, cardiovascular, "
                "hepatic, renal, neurological. Specify via --tissue or config.tissue_profile."
            )
        ref_path = Path(__file__).parent.parent / "references" / "tissue_context.json"
        with open(ref_path, 'r') as f:
            data = json.load(f)
        profiles = data.get("profiles", {})
        if self.tissue_profile not in profiles:
            available = ", ".join(profiles.keys())
            raise ValueError(f"Unknown tissue profile '{self.tissue_profile}'. Available: {available}")

        profile = dict(profiles[self.tissue_profile])

        # v0.5 P1-8: Merge special_gene_lists with sync + user extensions
        try:
            from dgra_gene_sync import get_merged_gene_lists_sync
            try:
                merged_lists = get_merged_gene_lists_sync(
                    tissue_profile=self.tissue_profile,
                    offline_mode=self.offline_mode,
                    sync_enabled=self.gene_sync_enabled,
                    ttl_days=self.gene_sync_ttl_days,
                    force_sync=force_sync,
                )
                if merged_lists:
                    profile["special_gene_lists"] = merged_lists
            except RuntimeError as e:
                # asyncio.run() fails when already in a running event loop
                print(f"[GPA] Gene list sync skipped (event loop conflict) - using static lists")
        except Exception as e:
            # Non-blocking: if merge fails, keep static lists
            print(f"[GPA] Gene list sync warning: {e} - using static lists")

        return profile

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

# Load pseudogene database from JSON (v0.5.3)
_PSEUDOGENE_DB_PATH = Path(__file__).resolve().parent.parent / "references" / "pseudogene_database.json"
_PSEUDOGENE_CONFIG = {}

def _load_pseudogene_database():
    """Load pseudogene database from JSON. Falls back to empty if missing."""
    global _PSEUDOGENE_CONFIG
    if _PSEUDOGENE_CONFIG:
        return _PSEUDOGENE_CONFIG
    try:
        with open(_PSEUDOGENE_DB_PATH, "r") as f:
            db = json.load(f)
        _PSEUDOGENE_CONFIG = db.get("genes", {})
        print(f"[GPA] Loaded pseudogene database: {len(_PSEUDOGENE_CONFIG)} genes")
        return _PSEUDOGENE_CONFIG
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"[GPA] WARNING: pseudogene_database.json not found or invalid ({e}). Using empty config.")
        _PSEUDOGENE_CONFIG = {}
        return _PSEUDOGENE_CONFIG


# Legacy hardcoded config - kept as fallback if JSON missing
_PSEUDOGENE_CONFIG_LEGACY = {
    "SETBP1": {"pseudogene": "VWFP1", "vaf_expected": 0.5, "vaf_min": 0.30, "notes": "VWFP1 on chr1 shares homology"},
    "VWF": {"pseudogene": "VWFP1", "vaf_expected": 0.5, "vaf_min": 0.25, "notes": "VWFP1 exon 23-34 homology"},
    "GBA": {"pseudogene": "GBAP1", "vaf_expected": 0.5, "vaf_min": 0.20, "notes": "GBAP1 on chr1"},
    "PMS2": {"pseudogene": "PMS2CL", "vaf_expected": 0.5, "vaf_min": 0.25, "notes": "PMS2CL on chr7"},
}


# =============================================================================
# v0.6: Lightweight Pseudogene Lookup (manual curation + Ensembl REST)
# =============================================================================

_PSEUDOGENE_LOOKUP_PATH = Path(__file__).resolve().parent.parent / "references" / "pseudogene_lookup.json"
_PSEUDOGENE_LOOKUP: Optional[Dict] = None


def _load_pseudogene_lookup() -> Dict:
    """
    Load lightweight pseudogene lookup from JSON.
    Returns dict keyed by functional gene symbol.
    Falls back to empty dict if file missing.
    """
    global _PSEUDOGENE_LOOKUP
    if _PSEUDOGENE_LOOKUP is not None:
        return _PSEUDOGENE_LOOKUP

    try:
        with open(_PSEUDOGENE_LOOKUP_PATH, "r", encoding="utf-8") as f:
            data = json.load(f)
        _PSEUDOGENE_LOOKUP = data.get("pairs", {})
        print(f"[GPA] Loaded pseudogene lookup: {len(_PSEUDOGENE_LOOKUP)} genes")
        return _PSEUDOGENE_LOOKUP
    except (FileNotFoundError, json.JSONDecodeError) as e:
        print(f"[GPA] INFO: pseudogene_lookup.json not found ({e}). Using empty lookup.")
        _PSEUDOGENE_LOOKUP = {}
        return _PSEUDOGENE_LOOKUP


def get_pseudogenes_for_gene(gene: str, offline_mode: bool = True) -> List[str]:
    """
    Return list of known pseudogenes for a given functional gene.

    Resolution order:
      1. Local pseudogene_lookup.json (manual curation)
      2. Legacy pseudogene_database.json
      3. (Future) Ensembl REST API - not yet implemented

    Args:
        gene: Functional gene symbol (e.g., "VWF")
        offline_mode: If True, skip API calls (default for reliability)

    Returns:
        List of pseudogene names (deduplicated). Empty if none known.
    """
    pseudogenes: List[str] = []

    # 1. v0.6 lookup
    lookup = _load_pseudogene_lookup()
    entry = lookup.get(gene)
    if entry:
        pseudogenes.extend(entry.get("pseudogenes", []))

    # 2. Legacy database (v0.5.3)
    legacy = _load_pseudogene_database()
    legacy_cfg = legacy.get(gene) or _PSEUDOGENE_CONFIG_LEGACY.get(gene)
    if legacy_cfg:
        legacy_pgs = legacy_cfg.get("pseudogenes", [])
        if isinstance(legacy_pgs, str):
            legacy_pgs = [legacy_pgs]
        pseudogenes.extend(legacy_pgs)

    # 3. (Future) Ensembl REST query when offline_mode=False

    # Deduplicate while preserving order
    seen = set()
    result = []
    for pg in pseudogenes:
        if pg and pg not in seen:
            seen.add(pg)
            result.append(pg)

    return result


def _calculate_pseudogene_score(
    observed_vaf: float,
    genotype: str,
    pseudogenes: List[str],
    gene: str,
) -> Dict:
    """
    Calculate pseudogene interference score (0.0-1.0).

    Scoring rules:
      - Homozygous (1/1): score=0 (no pseudogene interference expected)
      - Heterozygous (0/1) with VAF ~0.50: score=0
      - Heterozygous with VAF <0.30: score=0.8+ (strong interference)
      - Heterozygous with VAF 0.30-0.40: score=0.4-0.7 (suspected)
      - Heterozygous with VAF >0.65: score=0.3 (possible read bias)

    Returns dict with score, level, and recommendation.
    """
    if genotype in ["1/1", "1|1"]:
        return {
            "score": 0.0,
            "level": "none",
            "reason": "Homozygous variant - pseudogene interference unlikely",
            "recommendation": None,
        }

    if genotype not in ["0/1", "0|1"]:
        return {
            "score": 0.0,
            "level": "unknown_genotype",
            "reason": f"Genotype {genotype} not evaluable",
            "recommendation": None,
        }

    # Heterozygous scoring
    if observed_vaf < 0.20:
        return {
            "score": 0.9,
            "level": "strong_interference",
            "reason": f"VAF={observed_vaf:.2f} far below expected 0.50",
            "recommendation": "Strong evidence of pseudogene interference. Sanger validation or long-read sequencing strongly recommended.",
        }
    elif observed_vaf < 0.30:
        return {
            "score": 0.75,
            "level": "interference",
            "reason": f"VAF={observed_vaf:.2f} significantly below expected 0.50",
            "recommendation": "Probable pseudogene interference. Consider Sanger validation.",
        }
    elif observed_vaf < 0.40:
        return {
            "score": 0.55,
            "level": "suspected",
            "reason": f"VAF={observed_vaf:.2f} below expected 0.50",
            "recommendation": "Possible pseudogene interference. Caution advised.",
        }
    elif observed_vaf > 0.65:
        return {
            "score": 0.30,
            "level": "bias_suspected",
            "reason": f"VAF={observed_vaf:.2f} above expected 0.50",
            "recommendation": "Possible read bias or allele-specific expression. Review alignment.",
        }
    else:
        return {
            "score": 0.0,
            "level": "none",
            "reason": f"VAF={observed_vaf:.2f} within expected heterozygous range",
            "recommendation": None,
        }


def detect_pseudogene_artifact(variant: Variant) -> Optional[Dict]:
    """
    Detect pseudogene interference using v0.6 lookup system.

    v0.6: Unified detection using pseudogene_lookup.json + legacy DB.
    - Checks if gene has known pseudogenes
    - Evaluates VAF against expected heterozygous ratio
    - Returns structured score instead of binary flag

    Does NOT modify tier directly. Only provides evidence for confidence adjustment.
    """
    gene = variant.gene
    pseudogenes = get_pseudogenes_for_gene(gene)
    if not pseudogenes:
        return None

    observed_vaf = variant.vaf
    genotype = variant.gt

    if observed_vaf is None or genotype is None:
        return None

    score_info = _calculate_pseudogene_score(observed_vaf, genotype, pseudogenes, gene)

    # Build evidence dict
    if score_info["score"] >= 0.75:
        artifact_type = "PSEUDOGENE_INTERFERENCE"
    elif score_info["score"] >= 0.40:
        artifact_type = "PSEUDOGENE_SUSPECTED"
    else:
        return None  # No meaningful interference detected

    # Fetch lookup entry for metadata
    lookup = _load_pseudogene_lookup()
    entry = lookup.get(gene, {})
    strategy = entry.get("detection_strategy", "vaf_mismatch")
    confidence = entry.get("confidence", "medium")
    notes = entry.get("notes", "")

    return {
        "type": artifact_type,
        "gene": gene,
        "pseudogenes": pseudogenes,
        "strategy": strategy,
        "score": score_info["score"],
        "level": score_info["level"],
        "observed_vaf": observed_vaf,
        "expected_vaf": 0.5,
        "notes": notes,
        "recommendation": score_info["recommendation"],
        "confidence": confidence,
    }

# =============================================================================
# Module C: gnomAD Frequency Handler
# =============================================================================

def classify_gnomad_frequency(af: Optional[float], gene: str,
                               af_by_population: Optional[Dict[str, Dict]] = None,
                               target_population: Optional[str] = None) -> Dict:
    """
    Classify gnomAD allele frequency with zero-frequency handling.
    v0.5 P1-1: Population subgroup frequencies (EAS, AMR, AFR, NFE, SAS, etc.)

    If target_population is specified and available in af_by_population,
    uses that population's AF for classification instead of overall AF.
    """
    # Genes with germline warnings (clonal hematopoiesis filtering)
    GERMLINE_WARNING_GENES = {"ASXL1", "DNMT3A", "TET2", "TP53", "PPM1D", "JAK2", "CBL", "IDH1", "IDH2"}

    # Determine effective AF for classification
    effective_af = af
    pop_note = ""

    if target_population and af_by_population and target_population in af_by_population:
        pop_data = af_by_population[target_population]
        pop_af = pop_data.get("af")
        if pop_af is not None:
            effective_af = pop_af
            pop_note = f" Using {target_population} AF={pop_af:.6f} (overall AF={af})."

    # If overall AF is available but target population is missing, note it
    if target_population and (not af_by_population or target_population not in (af_by_population or {})):
        if af is not None:
            pop_note = f" {target_population} data unavailable; using overall AF={af}."
        else:
            pop_note = f" {target_population} data unavailable."

    if effective_af is None:
        result = {
            "af": "NOT_CAPTURED",
            "status": "NOT_CAPTURED",
            "interpretation": "Population frequency data unavailable. Cannot judge benign based on gnomAD." + pop_note,
            "action": "Continue with other modules (domain, ClinVar, zygosity).",
            "risk_adjustment": "Do NOT downgrade risk due to gnomAD absence."
        }
        if af_by_population:
            result["af_populations"] = af_by_population
        return result

    if gene in GERMLINE_WARNING_GENES:
        result = {
            "af": effective_af,
            "status": "GERMLINE_WARNING",
            "interpretation": f"{gene} is filtered in gnomAD germline due to clonal hematopoiesis." + pop_note,
            "action": "Use alternative databases (ExAC pre-filtered, ClinVar, LOVD).",
            "risk_adjustment": "Evaluate independently of gnomAD."
        }
        if af_by_population:
            result["af_populations"] = af_by_population
        return result

    if effective_af > 0.01:
        result = {
            "af": effective_af,
            "status": "common_polymorphism",
            "interpretation": "High population frequency (>1%), likely benign." + pop_note,
            "action": "Require strong evidence to upgrade risk.",
            "risk_adjustment": "Default Tier 3 unless other evidence is strong."
        }
    elif effective_af > 0.001:
        result = {
            "af": effective_af,
            "status": "low_frequency",
            "interpretation": "Moderate frequency (0.1-1%), needs functional assessment." + pop_note,
            "action": "Include in comprehensive evaluation."
        }
    else:
        result = {
            "af": effective_af if effective_af else "extremely_rare",
            "status": "rare_variant",
            "interpretation": "Very rare in population (<0.1%), likely under negative selection." + pop_note,
            "action": "Focus attention, requires literature search.",
            "risk_adjustment": "Do NOT downgrade risk by default."
        }

    if af_by_population:
        result["af_populations"] = af_by_population
    return result

# =============================================================================
# Module C.5: HGNC Gene Symbol Normalization (v0.5 P1-2)
# =============================================================================

def normalize_gene_symbols(
    variants: List[Variant],
    hgnc_results: Dict[str, Dict[str, Any]],
    offline_mode: bool = False,
) -> List[str]:
    """
    Normalize gene symbols using HGNC API results.

    Mutates Variant objects in-place:
    - Sets gene_original to the original symbol
    - Updates gene to the approved HGNC symbol
    - Sets transcript_warning for outdated/withdrawn/unknown symbols

    Returns list of WARNING strings for the report.

    Args:
        variants: Variant objects to normalize
        hgnc_results: Dict mapping original symbol -> HGNC query result
        offline_mode: If True, skip validation for known symbols; warn on unknown
    """
    warnings = []

    for v in variants:
        original = v.gene
        v.gene_original = original

        hgnc = hgnc_results.get(original, {})
        status = hgnc.get("status", "query_failed")
        approved = hgnc.get("approved_symbol", original)

        if status == "approved":
            # Symbol is valid, may need case correction
            if original != approved:
                v.gene = approved
                warnings.append(
                    f"HGNC: Gene symbol corrected \"{original}\" → \"{approved}\" (approved)"
                )
            # else: exact match, no change needed

        elif status in ("previous", "alias"):
            # Outdated symbol - auto-replace with approved
            v.gene = approved
            prev_syms = ", ".join(hgnc.get("previous_symbols", [])) or "N/A"
            alias_syms = ", ".join(hgnc.get("alias_symbols", [])) or "N/A"
            warnings.append(
                f"HGNC WARNING: Outdated symbol \"{original}\" → \"{approved}\" "
                f"(status={status}, previous={prev_syms}, alias={alias_syms})"
            )
            # Set transcript_warning for visibility
            tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
            tw["hgnc_warning"] = {
                "original": original,
                "approved": approved,
                "status": status,
                "action": "Symbol auto-corrected. Verify variant coordinates match approved symbol.",
            }
            v.transcript_warning = json.dumps(tw)

        elif status == "withdrawn":
            # Symbol withdrawn - cannot map, keep original, mark warning
            warnings.append(
                f"HGNC WARNING: Symbol \"{original}\" is WITHDRAWN. "
                f"Approved symbol: \"{approved}\". Cannot auto-correct - manual review required."
            )
            tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
            tw["hgnc_warning"] = {
                "original": original,
                "approved": approved,
                "status": "withdrawn",
                "action": "Symbol withdrawn. Manual mapping required. Verify coordinates independently.",
            }
            v.transcript_warning = json.dumps(tw)

        elif status in ("not_found", "query_failed"):
            # Unknown symbol - keep original, mark warning
            if offline_mode:
                # In offline mode, we can't validate. Only warn if symbol looks suspicious.
                # Known symbols from our common lists are accepted without warning.
                known_symbols = {
                    "VWF", "F8", "F9", "F7", "F10", "F11", "F12", "F13A1", "F13B",
                    "PROC", "PROS1", "SERPINC1", "SERPIND1", "PLG", "THBD",
                    "BRCA1", "BRCA2", "TP53", "PTEN", "APC", "MLH1", "MSH2", "MSH6", "PMS2",
                    "ATM", "CHEK2", "PALB2", "CDH1", "STK11",
                    "RUNX1", "CEBPA", "GATA2", "ETV6", "DDX41", "SAMD9", "SAMD9L",
                    "ASXL1", "BCOR", "BCORL1", "PHF6", "IDH1", "IDH2", "DNMT3A", "TET2",
                    "FLT3", "NPM1", "NRAS", "KRAS", "KIT",
                    "CYP2D6", "CYP2C19", "CYP3A4", "CYP3A5", "ABCB1", "TPMT", "DPYD", "UGT1A1",
                    "SCN5A", "KCNQ1", "KCNH2", "RYR2", "DSP", "PKP2", "LMNA", "TTN",
                }
                if original not in known_symbols:
                    warnings.append(
                        f"HGNC WARNING: Symbol \"{original}\" could not be validated (offline mode). "
                        f"Not in known gene lists - may be outdated or misspelled."
                    )
                    tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
                    tw["hgnc_warning"] = {
                        "original": original,
                        "status": "unvalidated_offline",
                        "action": "Offline mode: symbol not in known lists. Manual verification recommended.",
                    }
                    v.transcript_warning = json.dumps(tw)
                # else: known symbol, silently accept in offline mode
            else:
                # Online mode: API returned not_found - definitely warn
                warnings.append(
                    f"HGNC WARNING: Symbol \"{original}\" not found in HGNC database. "
                    f"May be misspelled, outdated, or non-standard. Keeping original symbol."
                )
                tw = json.loads(v.transcript_warning) if v.transcript_warning else {}
                tw["hgnc_warning"] = {
                    "original": original,
                    "status": status,
                    "action": "Symbol not found in HGNC. Verify spelling or check if gene has been renamed.",
                }
                v.transcript_warning = json.dumps(tw)

    return warnings


# =============================================================================
# Module D: Protein Domain Mapping
# =============================================================================

# =============================================================================
# Module D: Protein Domain Mapping  (v0.4: UniProt API)
# =============================================================================

# v0.3 hardcoded PROTEIN_DOMAINS dict removed - now fetched from UniProt API
# Protein domains are loaded via dgra_api.query_uniprot_by_gene() and passed
# as uniprot_data dict to this function.

def parse_protein_position(hgvsp: str) -> Optional[int]:
    """Extract amino acid position from p. string. Handles NP_ prefix."""
    hgvsp = str(hgvsp) if hgvsp is not None else ""
    if not hgvsp or hgvsp == "" or hgvsp == "nan":
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
# Module D.5: Gene Constraint Evaluation (v0.5 P1-4: pLI/LOEUF)
# =============================================================================

def evaluate_gene_constraint(variant: Variant) -> Dict[str, Any]:
    """
    Evaluate gnomAD gene constraint for variant tier modulation.

    Only applies to LOF variants (frameshift, nonsense, splice, start_lost, stop_gained).
    Non-LOF variants get constraint info for report but no tier upgrade.

    Returns:
    {
        "constraint_level": "strong" | "moderate" | "tolerant" | "unknown",
        "pLI": 0.99,
        "loeuf": 0.12,
        "lof_z": 5.2,
        "mis_z": 3.1,
        "is_lof": True,
        "tier_adjustment": 1,  # -1, 0, or +1 (only +1 for strong LOF constraint)
        "reason": "LOF-intolerant gene (pLI=0.99, LOEUF=0.12)",
    }
    """
    gc = variant.gene_constraint
    if not gc or gc.get("status") in ("NO_CONSTRAINT_DATA", "QUERY_FAILED"):
        return {
            "constraint_level": "unknown",
            "pLI": None,
            "loeuf": None,
            "lof_z": None,
            "mis_z": None,
            "is_lof": False,
            "tier_adjustment": 0,
            "reason": "No constraint data available",
        }

    pLI = gc.get("pLI")
    loeuf = gc.get("loeuf")
    lof_z = gc.get("lof_z")
    mis_z = gc.get("mis_z")

    # Determine if variant is LOF
    lof_consequences = {
        "frameshift", "nonsense", "stop_gained", "start_lost",
        "splice_donor", "splice_acceptor", "splice_site",
    }
    is_lof = any(lof_term in variant.consequence.lower() for lof_term in lof_consequences)

    # Determine constraint level
    constraint_level = "unknown"
    tier_adjustment = 0
    reason = ""

    if pLI is not None and loeuf is not None:
        if pLI >= 0.9 or loeuf <= 0.35:
            constraint_level = "strong"
            if is_lof:
                tier_adjustment = 1
                reason = f"LOF-intolerant gene (pLI={pLI:.2f}, LOEUF={loeuf:.2f}) - heterozygous LOF likely pathogenic"
            else:
                reason = f"LOF-intolerant gene (pLI={pLI:.2f}, LOEUF={loeuf:.2f}) - but variant is non-LOF, no tier upgrade"
        elif (pLI >= 0.5 and pLI < 0.9) or (loeuf >= 0.35 and loeuf <= 0.8):
            constraint_level = "moderate"
            reason = f"Moderate LOF constraint (pLI={pLI:.2f}, LOEUF={loeuf:.2f})"
        elif pLI < 0.5 and loeuf > 0.8:
            constraint_level = "tolerant"
            reason = f"LOF-tolerant gene (pLI={pLI:.2f}, LOEUF={loeuf:.2f}) - background LOF tolerated"
    elif pLI is not None:
        # Only pLI available
        if pLI >= 0.9:
            constraint_level = "strong"
            if is_lof:
                tier_adjustment = 1
                reason = f"LOF-intolerant gene (pLI={pLI:.2f}) - heterozygous LOF likely pathogenic"
            else:
                reason = f"LOF-intolerant gene (pLI={pLI:.2f}) - but variant is non-LOF"
        elif pLI >= 0.5:
            constraint_level = "moderate"
            reason = f"Moderate LOF constraint (pLI={pLI:.2f})"
        else:
            constraint_level = "tolerant"
            reason = f"LOF-tolerant gene (pLI={pLI:.2f})"
    elif loeuf is not None:
        # Only LOEUF available
        if loeuf <= 0.35:
            constraint_level = "strong"
            if is_lof:
                tier_adjustment = 1
                reason = f"LOF-intolerant gene (LOEUF={loeuf:.2f}) - heterozygous LOF likely pathogenic"
            else:
                reason = f"LOF-intolerant gene (LOEUF={loeuf:.2f}) - but variant is non-LOF"
        elif loeuf <= 0.8:
            constraint_level = "moderate"
            reason = f"Moderate LOF constraint (LOEUF={loeuf:.2f})"
        else:
            constraint_level = "tolerant"
            reason = f"LOF-tolerant gene (LOEUF={loeuf:.2f})"
    else:
        reason = "Constraint metrics incomplete"

    return {
        "constraint_level": constraint_level,
        "pLI": pLI,
        "loeuf": loeuf,
        "lof_z": lof_z,
        "mis_z": mis_z,
        "is_lof": is_lof,
        "tier_adjustment": tier_adjustment,
        "reason": reason,
    }


# =============================================================================
# Module D.6: NMD Prediction (v0.5 P1-5: PVS1 refinement)
# =============================================================================

def predict_nmd(variant: Variant, ensembl_data: Optional[Dict] = None) -> Dict[str, Any]:
    """
    Predict nonsense-mediated decay (NMD) sensitivity based on truncation position.

    ClinGen PVS1 guidance:
    - Internal exons (up to penultimate exon, excluding last 50-55bp) → NMD sensitive → PVS1 applies
    - Last 50-55bp of penultimate exon → possible NMD escape → PVS1 downgraded to PM/PP
    - Last exon → NMD escape → PVS1 does NOT apply
    - UTR regions → NMD does NOT apply

    Returns:
    {
        "status": "sensitive" | "escape" | "possible_escape" | "not_applicable" | "unknown",
        "reason": "Truncation in internal exon → classic NMD",
        "confidence": "high" | "moderate" | "low",
        "pvs1_applicable": True/False,
        "pvs1_strength": "strong" | "moderate" | "weak" | "not_applicable",
    }
    """
    consequence = str(variant.consequence or "").lower()

    # Check if variant is a truncating variant (LOF)
    lof_terms = {"frameshift", "nonsense", "stop_gained", "start_lost"}
    is_truncating = any(term in consequence for term in lof_terms)

    # Splice variants are not handled by simple exon position rules
    if "splice" in consequence:
        return {
            "status": "unknown",
            "reason": "Splice variant - NMD prediction requires transcript-level analysis",
            "confidence": "low",
            "pvs1_applicable": True,  # Conservative: assume PVS1 applies for splice
            "pvs1_strength": "strong",
        }

    if not is_truncating:
        return {
            "status": "not_applicable",
            "reason": "Not a truncating variant - NMD prediction not applicable",
            "confidence": "high",
            "pvs1_applicable": False,
            "pvs1_strength": "not_applicable",
        }

    # Parse exon field (e.g., "2/15" or "15/15")
    exon_str = str(variant.exon or "").strip()
    if _is_unknown(exon_str):
        # No exon info - conservative assumption: NMD sensitive
        return {
            "status": "unknown",
            "reason": "NMD status unknown, assuming sensitive (conservative)",
            "confidence": "low",
            "pvs1_applicable": True,
            "pvs1_strength": "strong",
        }

    # Parse exon position
    try:
        if "/" in exon_str:
            current_exon, total_exons = exon_str.split("/", 1)
            current = int(current_exon.strip())
            total = int(total_exons.strip())
        elif " of " in exon_str.lower():
            # Format: "2 of 15"
            parts = exon_str.lower().split(" of ")
            current = int(parts[0].strip())
            total = int(parts[1].strip())
        else:
            # Single number - can't determine position
            return {
                "status": "unknown",
                "reason": f"Cannot determine exon position from '{exon_str}' - assuming sensitive",
                "confidence": "low",
                "pvs1_applicable": True,
                "pvs1_strength": "strong",
            }
    except (ValueError, IndexError):
        return {
            "status": "unknown",
            "reason": f"Cannot parse exon '{exon_str}' - assuming sensitive",
            "confidence": "low",
            "pvs1_applicable": True,
            "pvs1_strength": "strong",
        }

    if total <= 0:
        return {
            "status": "unknown",
            "reason": "Invalid total exon count - assuming sensitive",
            "confidence": "low",
            "pvs1_applicable": True,
            "pvs1_strength": "strong",
        }

    # Determine NMD status based on exon position
    if current == total:
        # Last exon - NMD escape
        return {
            "status": "escape",
            "reason": f"Truncation in last exon ({current}/{total}) - NMD escape",
            "confidence": "high",
            "pvs1_applicable": False,
            "pvs1_strength": "not_applicable",
        }
    elif current == total - 1:
        # Penultimate exon - possible escape if within last 50-55bp of CDS
        # Without exact CDS position, we use a conservative estimate:
        # If we have transcript length info from Ensembl, we could be more precise
        return {
            "status": "possible_escape",
            "reason": f"Truncation in penultimate exon ({current}/{total}) - possible NMD escape if within last 50-55bp",
            "confidence": "moderate",
            "pvs1_applicable": False,  # Conservative: don't apply PVS1 if uncertain
            "pvs1_strength": "moderate",  # Downgraded to PM/PP level
        }
    else:
        # Internal exon - classic NMD
        return {
            "status": "sensitive",
            "reason": f"Truncation in internal exon ({current}/{total}) - classic NMD",
            "confidence": "high",
            "pvs1_applicable": True,
            "pvs1_strength": "strong",
        }


# =============================================================================
# Module D.7: Missense Stratification (v0.5 P1-5)
# =============================================================================

def evaluate_missense_tier(variant: Variant, domain_info: Optional[Dict],
                           gene_constraint: Optional[Dict]) -> Dict[str, Any]:
    """
    Stratify missense variants by domain impact and conservation.

    Returns a score (0-1) and tier recommendation for missense variants.

    Scoring:
    - Domain-critical: completely_destroyed + mis_z > 3.09 → score=0.9
    - Likely damaging: partially_destroyed + mis_z 2-3.09 → score=0.6
    - Tolerated: tolerated/inter-domain + mis_z < 2 → score=0.1
    - Unknown: no domain info → score=0.3 (conservative)

    Returns:
    {
        "score": 0.9,
        "tier_recommendation": 2,  # or 3
        "category": "domain_critical" | "likely_damaging" | "tolerated" | "unknown",
        "reason": "Missense in critical domain residue (mis_z=4.2)",
    }
    """
    consequence = str(variant.consequence or "").lower()

    # Only applies to missense variants
    if "missense" not in consequence:
        return {
            "score": 0.0,
            "tier_recommendation": None,
            "category": "not_missense",
            "reason": "Not a missense variant",
        }

    # Get mis_z score
    mis_z = None
    if gene_constraint and gene_constraint.get("mis_z") is not None:
        mis_z = float(gene_constraint["mis_z"])

    # Get domain integrity
    domain_integrity = None
    if domain_info:
        domain_integrity = domain_info.get("domain_integrity")

    # Stratify
    if domain_integrity == "completely_destroyed":
        if mis_z is not None and mis_z > 3.09:
            return {
                "score": 0.9,
                "tier_recommendation": 2,
                "category": "domain_critical",
                "reason": f"Missense in critical domain residue (mis_z={mis_z:.2f}) - high pathogenic potential",
            }
        else:
            return {
                "score": 0.7,
                "tier_recommendation": 2,
                "category": "likely_damaging",
                "reason": f"Missense destroys domain structure (mis_z={mis_z:.2f if mis_z else 'N/A'})",
            }
    elif domain_integrity == "partially_destroyed":
        if mis_z is not None and mis_z >= 2.0:
            return {
                "score": 0.6,
                "tier_recommendation": 2,
                "category": "likely_damaging",
                "reason": f"Missense in conserved domain (mis_z={mis_z:.2f}) - likely damaging",
            }
        else:
            mis_z_fmt = f"{mis_z:.2f}" if mis_z is not None else "N/A"
            return {
                "score": 0.4,
                "tier_recommendation": 2,
                "category": "possibly_damaging",
                "reason": f"Missense partially disrupts domain (mis_z={mis_z_fmt})",
            }
    elif domain_integrity == "tolerated":
        if mis_z is not None and mis_z < 2.0:
            return {
                "score": 0.1,
                "tier_recommendation": 3,
                "category": "tolerated",
                "reason": f"Missense in non-critical region (mis_z={mis_z:.2f}) - likely tolerated",
            }
        else:
            return {
                "score": 0.3,
                "tier_recommendation": 2,
                "category": "uncertain",
                "reason": f"Missense tolerated by structure but in conserved region (mis_z={mis_z:.2f})",
            }
    else:
        # No domain info or unknown
        if mis_z is not None and mis_z > 3.09:
            return {
                "score": 0.5,
                "tier_recommendation": 2,
                "category": "conservation_concern",
                "reason": f"Missense in highly constrained gene (mis_z={mis_z:.2f}) - domain info unavailable",
            }
        else:
            return {
                "score": 0.3,
                "tier_recommendation": 2,
                "category": "unknown",
                "reason": "Missense with unknown domain impact - conservative Tier 2",
            }


# =============================================================================
# Module E: Tissue Relevance Assessment  (v0.4: GTEx API + fallback)
# =============================================================================

def aggregate_gtex_expression(results: List[Dict[str, Any]],
                               strategy: str = "max") -> Dict[str, Any]:
    """
    Aggregate multi-tissue GTEx expression results.

    v0.5 P1-6: Multi-tissue GTEx aggregation. Takes a list of per-tissue
    GTEx results and returns a single aggregated result dict compatible
    with the single-tissue format.

    Args:
        results: List of GTEx result dicts (from query_gtex_expression_multi)
        strategy: Aggregation strategy - "max" (default, conservative),
                  "mean", or "median"

    Returns:
        Aggregated result dict with:
        - median_tpm: aggregated TPM value
        - max_tpm: maximum TPM across tissues
        - mean_tpm: mean TPM across tissues
        - expressing_tissues: count of tissues with TPM >= 1.0
        - primary_tissues: list of tissues with TPM >= 10.0
        - all_tissues: list of (tissue, tpm) tuples
        - source: "gtex_multi"
    """
    if not results:
        return {
            "median_tpm": None,
            "max_tpm": None,
            "mean_tpm": None,
            "expressing_tissues": 0,
            "primary_tissues": [],
            "all_tissues": [],
            "source": "gtex_multi",
            "unit": "TPM",
        }

    # Extract valid TPM values
    valid = []
    all_tissues = []
    for r in results:
        tpm = r.get("median_tpm")
        tissue = r.get("tissue", "unknown")
        if tpm is not None and tpm >= 0:
            valid.append(tpm)
            all_tissues.append((tissue, tpm))

    if not valid:
        return {
            "median_tpm": None,
            "max_tpm": None,
            "mean_tpm": None,
            "expressing_tissues": 0,
            "primary_tissues": [],
            "all_tissues": all_tissues,
            "source": "gtex_multi",
            "unit": "TPM",
        }

    max_tpm = max(valid)
    mean_tpm = sum(valid) / len(valid)

    # Determine aggregated TPM based on strategy
    if strategy == "max":
        aggregated_tpm = max_tpm
    elif strategy == "mean":
        aggregated_tpm = mean_tpm
    elif strategy == "median":
        sorted_tpms = sorted(valid)
        mid = len(sorted_tpms) // 2
        aggregated_tpm = sorted_tpms[mid] if len(sorted_tpms) % 2 == 1 else (sorted_tpms[mid - 1] + sorted_tpms[mid]) / 2
    else:
        aggregated_tpm = max_tpm

    expressing_tissues = sum(1 for t in valid if t >= 1.0)
    primary_tissues = [tissue for tissue, tpm in all_tissues if tpm >= 10.0]

    return {
        "median_tpm": aggregated_tpm,  # Named "median_tpm" for backward compatibility
        "max_tpm": max_tpm,
        "mean_tpm": mean_tpm,
        "expressing_tissues": expressing_tissues,
        "primary_tissues": primary_tissues,
        "all_tissues": all_tissues,
        "source": "gtex_multi",
        "unit": "TPM",
    }


def assess_tissue_relevance(variant: Variant, tissue_profile: Dict,
                            gtex_data: Dict[str, Dict]) -> Dict:
    """
    Assess if gene is relevant to target tissue/organ context.

    Priority:
    1. GTEx API expression data (if available) - auto-classify by TPM thresholds
    2. Local tissue_context.json profile (fallback for API failures)
    3. Unknown if neither available - conservative, do NOT fast-track

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

    # v0.5 P1-6: Check if this is multi-tissue aggregated data
    is_multi_tissue = gtex.get("source") == "gtex_multi"
    all_tissues = gtex.get("all_tissues", [])
    max_tpm = gtex.get("max_tpm")
    mean_tpm = gtex.get("mean_tpm")
    primary_tissues = gtex.get("primary_tissues", [])
    expressing_tissues = gtex.get("expressing_tissues", 0)

    # v0.10.8: Use phenotype-aware TPM when available (Phase 2 full GTEx query)
    phenotype_max_tpm = gtex.get("phenotype_max_tpm")
    global_max_tpm = gtex.get("global_max_tpm")
    phenotype_tissues = gtex.get("phenotype_tissues", [])
    
    # For relevance assessment, prefer phenotype-matched tissues
    # If phenotype tissues are not in GTEx (e.g. retina not in GTEx v8),
    # phenotype_max_tpm will be 0 even if the gene is tissue-specific
    assess_tpm = phenotype_max_tpm if phenotype_max_tpm is not None else tpm
    
    # v0.10.8: GTEx fast-track REMOVED. Expression data is now used only for
    # phenotype-tissue association, not as a hard tier gate.
    # The fast_track field is kept for backward compatibility but ignored
    # by the tier classifier (see gpa_tier_classifier.py v0.10.8).
    
    # v0.5 P0-6: General profile - skip GTEx tissue-specific fast-track.
    # When gtex_tissue is null, expression-based fast-track is disabled;
    # assessment relies on special gene lists and ClinVar/gnomAD instead.
    if assess_tpm is not None and tissue_profile.get("gtex_tissue") is not None:
        # v0.5 P1-6: Multi-tissue rationale
        if is_multi_tissue and all_tissues:
            tissue_count = len(all_tissues)
            if assess_tpm >= 10.0:
                relevance = "primary"
                rationale = f"High phenotype-relevant expression (max TPM={assess_tpm:.1f} across matched tissues) per GTEx."
            elif assess_tpm >= 1.0:
                relevance = "secondary"
                rationale = f"Moderate phenotype-relevant expression (max TPM={assess_tpm:.1f} across matched tissues) per GTEx."
            elif assess_tpm > 0:
                relevance = "none"
                rationale = f"Low phenotype-relevant expression (max TPM={assess_tpm:.2f} across matched tissues) per GTEx."
            else:
                relevance = "none"
                # v0.10.8: Note when GTEx lacks the phenotype-relevant tissue
                if phenotype_tissues and not any("retina" in t.lower() for t in phenotype_tissues):
                    missing = ", ".join(phenotype_tissues[:3])
                    rationale = f"No detectable expression in GTEx-matched tissues ({missing}). Note: GTEx v8 lacks many specialized tissues (e.g. retina)."
                else:
                    rationale = f"No detectable {profile_name} expression across {tissue_count} tissues per GTEx."
        else:
            # Single tissue (backward compatible)
            if assess_tpm >= 10.0:
                relevance = "primary"
                rationale = f"High {profile_name} expression (TPM={assess_tpm:.1f}) per GTEx."
            elif assess_tpm >= 1.0:
                relevance = "secondary"
                rationale = f"Moderate {profile_name} expression (TPM={assess_tpm:.1f}) per GTEx."
            elif assess_tpm > 0:
                relevance = "none"
                rationale = f"Low {profile_name} expression (TPM={assess_tpm:.2f}) per GTEx."
            else:
                relevance = "none"
                rationale = f"No detectable {profile_name} expression per GTEx."

        # v0.10.8: Fast track REMOVED. Always return standard pipeline suggestion.
        if relevance == "none":
            clinical_note = f"{gene} has low/no GTEx expression in phenotype-matched tissues."
            if global_max_tpm and global_max_tpm > assess_tpm:
                clinical_note += f" However, global max TPM={global_max_tpm:.1f} in other tissues suggests tissue-specific expression."
            if phenotype_tissues and not any("retina" in t.lower() for t in phenotype_tissues):
                clinical_note += " Note: GTEx v8 lacks retinal/eye tissues; tissue-specific genes may be underrepresented."
            
            return {
                "tier_suggestion": "assess_via_standard_pipeline",
                "relevance": relevance,
                "reason": f"{gene} GTEx expression low in matched tissues (TPM={assess_tpm:.2f}).",
                "clinical_note": clinical_note,
                "fast_track": False,  # v0.10.8: disabled
                "rationale": rationale,
                "gtex_tpm": assess_tpm,
                "global_max_tpm": global_max_tpm,
                "source": gtex.get("source", "gtex"),
            }

        # Primary or secondary: standard pipeline
        # v0.5 P1-6: Enhanced clinical note for multi-tissue
        if is_multi_tissue and all_tissues:
            tissue_detail = "; ".join([f"{t}:{v:.1f}" for t, v in all_tissues])
            clinical_note = f"{gene} is {relevance}-relevant to phenotype-matched tissues (max TPM={assess_tpm:.1f} across {len(all_tissues)} tissues: {tissue_detail})."
        else:
            clinical_note = f"{gene} is {relevance}-relevant to {profile_name} (GTEx TPM={assess_tpm:.1f})."

        return {
            "tier_suggestion": "assess_via_standard_pipeline",
            "relevance": relevance,
            "gtex_tpm": assess_tpm,
            "global_max_tpm": global_max_tpm,
            "rationale": rationale,
            "fast_track": False,
            "clinical_note": clinical_note,
            "source": gtex.get("source", "gtex"),
        }

    # --- Priority 1b: GTEx data available but tissue is null (general profile) ---
    if tpm is not None and tissue_profile.get("gtex_tissue") is None:
        # For general profile, use max TPM across all queried tissues as a proxy
        # for "expressed somewhere important". Still don't fast-track based on
        # expression alone; rely on special lists.
        if tpm >= 10.0:
            relevance = "primary"
            rationale = f"High expression (TPM={tpm:.1f}) in at least one tissue. Relevant to general health."
        elif assess_tpm >= 1.0:
            relevance = "secondary"
            rationale = f"Moderate expression (TPM={assess_tpm:.1f}) in at least one tissue."
        else:
            relevance = "none"
            rationale = f"Low expression (TPM={assess_tpm:.2f}) - not prominently expressed."
        return {
            "tier_suggestion": "assess_via_standard_pipeline",
            "relevance": relevance,
            "gtex_tpm": assess_tpm,
            "global_max_tpm": global_max_tpm,
            "rationale": rationale,
            "fast_track": False,
            "clinical_note": f"{gene} general relevance: {relevance} (TPM={assess_tpm:.1f}). No tissue-specific fast-track.",
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
            # v0.10.8: Fast track REMOVED even for local fallback
            return {
                "tier_suggestion": "assess_via_standard_pipeline",
                "relevance": relevance,
                "reason": f"{gene} has no {profile_name} relevance per local database.",
                "clinical_note": f"No impact on {profile_name} function or safety per local database.",
                "fast_track": False,
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
        # general + hematopoietic (overlapping keys share same values)
        "cancer_predisposition": "primary",
        "cardiac_safety": "primary",
        "coagulation": "primary",
        "fa_dna_repair": "primary",
        "drug_metabolism": "secondary",
        "kir_cluster": "secondary",
        "immunodeficiency": "primary",
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
        "clinical_note": "Gene relevance unknown - conservative assessment.",
        "gtex_rpkm": None,
        "source": "unknown",
    }

# =============================================================================
# Module F: Three-Tier Risk Classification
# =============================================================================

TIER1_ACTION_GENES = {
    "VWF": {"reason": "Coagulation disorder - vWD risk in patient", "condition": "ClinVar_Pathogenic"},
}

# v0.5.1 OPT-P0-2: X-linked gene female heterozygous risk adjustment
def _x_linked_female_adjustment(tier: int, chrom: str, gt: str,
                                gene_constraint: Optional[Dict] = None) -> Tuple[int, str]:
    """
    Adjust tier for X-linked genes in female heterozygous carriers.

    Biological basis:
    - Female XX: random X-inactivation (lyonization) ~50% cells express wild-type X
    - If gene is haplosufficient (pLI < 0.5 or LOEUF > 0.35): 50% wild-type sufficient
    - If gene is haploinsufficient (pLI > 0.9): maintain tier

    Returns: (adjusted_tier, reason)
    """
    if chrom not in ('X', 'chrX'):
        return tier, ""
    if gt not in ('0/1', '0|1'):
        return tier, ""

    pLI = 0.0
    loeuf = 1.0
    if gene_constraint:
        pLI = gene_constraint.get("pLI", 0.0) or 0.0
        loeuf = gene_constraint.get("loeuf", 1.0) or 1.0

    is_haplosufficient = (pLI < 0.5) or (loeuf > 0.35)
    is_haploinsufficient = (pLI > 0.9) and (loeuf < 0.35)

    if is_haplosufficient and tier == 1:
        return 2, (f"X-linked female het + haplosufficient "
                   f"(pLI={pLI:.2f}, LOEUF={loeuf:.2f}) - "
                   f"50% wild-type via X-inactivation sufficient")
    elif is_haplosufficient and tier == 2:
        return 3, (f"X-linked female het + haplosufficient - no concern")
    elif is_haploinsufficient:
        return tier, (f"X-linked female het but haploinsufficient "
                      f"(pLI={pLI:.2f}, LOEUF={loeuf:.2f}) - maintaining tier {tier}")
    return tier, ""


# v0.5.1 OPT-P1-3: C-terminal truncation severity assessment
def _is_minimal_c_terminal_truncation(hgvsp: Optional[str]) -> bool:
    """
    Assess if a HIGH-impact variant involves minimal C-terminal truncation.
    Such truncations (<5% of protein from C-terminus) are often benign
    due to NMD escape or non-critical tail domains.

    Returns True if the variant appears to be a minimal C-terminal truncation.
    """
    if not hgvsp:
        return False
    hgvsp_str = str(hgvsp)
    import re
    # Match nonsense: p.Glu293Ter
    ter_match = re.search(r'p\.[A-Za-z]+(\d+)Ter', hgvsp_str)
    if ter_match:
        stop_pos = int(ter_match.group(1))
        # Conservative heuristic: stop position >=280 aa is likely near C-terminus
        # for most proteins (median ~400 aa). Combined with ClinVar benign in caller.
        if stop_pos >= 280:
            return True
    # Match frameshift: p.Ile249LeufsTer3
    fs_match = re.search(r'p\.[A-Za-z]+(\d+)[A-Za-z]*fsTer(\d+)', hgvsp_str)
    if fs_match:
        original_pos = int(fs_match.group(1))
        fs_stop_count = int(fs_match.group(2))
        # If frameshift occurs late in protein AND early termination
        # e.g., position >=250 with fsTer within 10 aa -> minimal impact
        if original_pos >= 250 and fs_stop_count <= 10:
            return True
    return False


# =============================================================================
# Multi-hit Gene Detection (v0.5.1 OPT: ClinVar benign + synonymous exclusion)
# =============================================================================

def _variant_has_pathogenic_evidence(v: Variant, gtex_data: Optional[Dict] = None) -> bool:
    """
    Check if a variant has evidence suggesting pathogenicity.
    v0.5 P0-7: UNKNOWN fields treated conservatively - UNKNOWN impact is treated as HIGH,
    UNKNOWN clinvar does not trigger benign exclusion.

    Criteria (OR):
      1. Affects protein domain (has specific domain mapping, not unknown/inter-domain)
         AND gene is expressed in target tissue (GTEx TPM >= 1.0)
      2. ClinVar pathogenic/likely pathogenic or HIGH impact or rare gnomAD (<0.001)
      3. Splice site change (consequence contains 'splice')
    """
    # === OPT-P1-4: Synonymous (LOW) variants NEVER have pathogenic evidence ===
    if v.impact == "LOW":
        cons = str(v.consequence or "").lower()
        if "synonymous" in cons or "同义" in cons:
            return False

    # === OPT-P0-1: ClinVar benign exclusion for LOW/MODERATE impact ===
    clinvar_lower = str(v.clinvar or "").lower()
    is_benign_cv = (("benign" in clinvar_lower and "conflicting" not in clinvar_lower)
                    or "likely_benign" in clinvar_lower)
    if is_benign_cv and v.impact in ("LOW", "MODERATE"):
        return False

    # === OPT-P1-3: Minimal C-terminal truncation + ClinVar benign -> not pathogenic ===
    if v.impact == "HIGH" and is_benign_cv:
        if _is_minimal_c_terminal_truncation(v.hgvsp):
            return False

    # Quick check: splice site changes are always considered
    consequence_lower = str(v.consequence or "").lower()
    if "splice" in consequence_lower:
        return True

    # Check tissue expression for domain relevance
    tissue_tpm = None
    if gtex_data and v.gene in gtex_data:
        tissue_tpm = gtex_data[v.gene].get("median_tpm")

    # 1. Domain impact - only counts if gene is expressed in target tissue
    di = v.domain_info
    if di:
        domain = di.get("domain", "")
        if domain and domain not in ("unknown", "N/A", "inter-domain / unannotated"):
            if tissue_tpm is not None and tissue_tpm < 1.0:
                pass  # Low expression: domain not relevant for this tissue
            else:
                return True

    # 2. Pathogenic evidence
    if v.clinvar and not _is_unknown(v.clinvar) and "pathogenic" in clinvar_lower and "conflicting" not in clinvar_lower:
        return True
    # HIGH impact - but C-terminal truncation + ClinVar benign already filtered above
    if v.impact == "HIGH" or _is_unknown(v.impact):
        return True
    if v.gnomad_af is not None and v.gnomad_af < 0.001:
        return True

    # 3. Splice site changes - always considered
    if "splice" in consequence_lower:
        return True

    return False


# =============================================================================
# Phase Analysis (v0.4.5)
# v0.10.0: Extracted to gpa_phaser.py
# =============================================================================

# NOTE: Lazy imports at bottom of file moved into main() to avoid circular imports

# =============================================================================
# Main Pipeline  (v0.4: async + batch API queries)
# =============================================================================

# =============================================================================
# Pipeline
# v0.10.0: Extracted to gpa_pipeline.py
# =============================================================================

def main():
    # v0.10.0: Lazy imports to avoid circular dependency with gpa_* modules
    from gpa_phaser import PhaseResult, determine_phase
    from gpa_multi_hit import detect_multi_hit_genes
    from gpa_report import _get_version_info, generate_tier_report, generate_json_report
    from gpa_qc import _run_qc_checks
    from gpa_input import InputType, detect_input_type, variants_from_vep_annotation, parse_annotated_vcf
    from gpa_pipeline import run_dgra_pipeline, run_multi_organ_assessment

    parser = argparse.ArgumentParser(
        description="GPA - Genomic Phenotype Association (v0.4 API-first with cache)",
        epilog="Available tissue profiles: hematopoietic (default), cardiovascular, hepatic, renal, neurological. "
               "Define custom profiles in references/tissue_context.json"
    )
    parser.add_argument("--input", "-i", required=True, help="Input CSV/TSV with annotated variants")
    parser.add_argument("--output", "-o", default="dgra_report.md", help="Output report file")
    parser.add_argument("--json", "-j", help="Output JSON file with full structured results (backward compatible)")
    parser.add_argument("--output-json", dest="output_json",
                        help="Output P1-12 structured JSON report for downstream systems (v0.5 P1-12)")
    parser.add_argument("--tissue", "-t", default="general",
                        help="Tissue/organ context profile. "
                             "Controls which genes are considered relevant for tier classification. "
                             "Available: general, hematopoietic, cardiovascular, hepatic, renal, neurological. "
                             "Default: general. Mutually exclusive with --multi-organ.")
    parser.add_argument("--multi-organ", default=None,
                        help="Multi-organ assessment: comma-separated profiles, e.g. 'hematopoietic,cardiovascular'. "
                             "Runs independent assessment for each profile and generates a joint risk matrix. "
                             "Takes max tier across profiles. Mutually exclusive with --tissue. (v0.5 P1-7)")
    parser.add_argument("--offline", action="store_true",
                        help="Offline mode: skip all API calls, use cache + local references only")
    parser.add_argument("--somatic", action="store_true",
                        help="Somatic mode: tumor driver mutation analysis. "
                             "TSG truncating + oncogene hotspots = Tier 1")
    parser.add_argument("--sync-gene-lists", action="store_true",
                        help="Force sync special_gene_lists from external sources (Orphanet, OMIM) before analysis. "
                             "Bypasses cache TTL. (v0.5 P1-8)")
    parser.add_argument("--evidence-detail", choices=["brief", "full"], default="brief",
                        help="Evidence chain detail level in report: brief (top 3) or full (all). (v0.5 P1-9)")
    parser.add_argument("--target-population", "--population",
                        choices=["EAS", "AMR", "AFR", "NFE", "SAS", "ASJ", "FIN", "MID", "OTH"],
                        help="Target population for gnomAD subgroup AF classification (v0.5 P1-1). "
                             "Uses that population's AF instead of overall AF for frequency-based tiering.")
    parser.add_argument("--database-version",
                        help="Freeze analysis to a specific database version for reproducibility "
                             "(e.g., 'gnomAD v4.1'). Recorded in output meta. (v0.5 P1-15)")
    # v0.5 P2-3: YAML config file support
    parser.add_argument("--config", "-c", type=Path, default=None,
                        help="Path to dgra.yaml configuration file. Overrides built-in defaults. "
                             "If not specified, uses references/dgra.yaml if it exists, "
                             "otherwise falls back to built-in defaults. (v0.5 P2-3)")

    # v0.7: Phenotype association
    parser.add_argument("--phenotypes", default=None,
                        help="Clinical phenotype description for phenotype-gene association analysis. "
                             "e.g. 'distal muscle weakness, myopathic damage, slow progression'. "
                             "Only applied to Tier 1/2 variants. Requires LLM API key for best accuracy.")
    parser.add_argument("--filter-preset", default=None, choices=["strict", "clinical", "broad"],
                        help="Filter preset name used for pre-filtering. Displayed in report header. (v0.7.1)")
    parser.add_argument("--filter-stats", default=None,
                        help="JSON string of filter statistics. Displayed in report header. (v0.7.1)")
    parser.add_argument("--llm-model", default="gpt-4o-mini",
                        help="LLM model for phenotype semantic matching. Default: gpt-4o-mini. "
                             "Alternative: gpt-4o, claude-3-haiku.")
    # v0.8.0: SpliceAI splice-prediction integration (default OFF)
    parser.add_argument("--spliceai", action="store_true",
                        help="Enable SpliceAI splice-prediction lookup for splice variants. "
                             "Default OFF — only applies to canonical splice (acceptor/donor) and splice_region. "
                             "(v0.8.0)")
    parser.add_argument("--spliceai-concurrency", type=int, default=5,
                        help="Max concurrent SpliceAI API requests (default: 5). (v0.8.0)")

    # v0.9.0: VCF annotation + transcript selection
    parser.add_argument("--disease-description", default=None,
                        help="Clinical disease description for disease-aware transcript selection. "
                             "e.g. 'limb-girdle muscular dystrophy, proximal muscle weakness'. "
                             "Only used for raw VCF input; optional — falls back to canonical/MANE if not provided.")
    parser.add_argument("--annotator", default="auto", choices=["auto", "vep_api", "vep_local"],
                        help="Variant annotator for raw VCF: auto (default, zero-config VEP API), "
                             "vep_api (Ensembl REST), vep_local (local VEP command). (v0.9.0)")
    parser.add_argument("--vep-cache", default=None,
                        help="Path to local VEP cache directory. Required for --annotator vep_local. (v0.9.0)")

    # v0.10.1: Two-phase pipeline for large VCF datasets
    parser.add_argument("--two-phase", action="store_true",
                        help="Enable two-phase pipeline: fast local triage first, then API enrichment only for "
                             "Tier 1/2 candidates. Reduces API calls by 50-200x for large VCFs. "
                             "Recommended for VCF input with > 1000 variants. (v0.10.1)")

    args = parser.parse_args()

    # v0.5 P2-3: Load YAML config if provided or default exists
    file_config = None
    config_path = args.config
    if config_path is None and DEFAULT_CONFIG_PATH.exists():
        config_path = DEFAULT_CONFIG_PATH

    if config_path:
        try:
            file_config = GPAFileConfig.from_yaml(config_path)
            print(f"[GPA] Loaded configuration from {config_path}")
        except FileNotFoundError:
            print(f"[GPA] Config file not found: {config_path}, using built-in defaults")
        except Exception as e:
            print(f"[GPA] Warning: Failed to load config {config_path}: {e}")

    # v0.5 P1-7: Validate --multi-organ vs --tissue mutual exclusion
    multi_organ = None
    if args.multi_organ:
        if args.tissue != "general":
            # tissue was explicitly set (not default)
            print("Error: --tissue and --multi-organ are mutually exclusive. Use one or the other.")
            sys.exit(1)
        multi_organ = [p.strip() for p in args.multi_organ.split(",") if p.strip()]
        valid_tissues = {"general", "hematopoietic", "cardiovascular", "hepatic", "renal", "neurological"}
        invalid = [p for p in multi_organ if p not in valid_tissues]
        if invalid:
            print(f"Error: Invalid multi-organ profile(s): {', '.join(invalid)}. Valid: {', '.join(sorted(valid_tissues))}")
            sys.exit(1)
        if len(multi_organ) < 1 or len(multi_organ) > 3:
            print("Error: --multi-organ requires 1-3 profiles.")
            sys.exit(1)
        print(f"Multi-organ assessment: {', '.join(multi_organ)}")

    # v0.9.0: Input type detection
    input_type = detect_input_type(args.input)
    print(f"[GPA] Input type detected: {input_type.value}")

    variants_data = []
    if input_type == InputType.RAW_VCF:
        # v0.9.0: Annotate raw VCF
        print("[GPA] Raw VCF detected — starting annotation pipeline...")
        from gpa_vcf_annotator import VCFAnnotator
        from gpa_transcript_selector import TranscriptSelector

        annotator_name = args.annotator if hasattr(args, 'annotator') else "auto"
        vep_cache_path = args.vep_cache if hasattr(args, 'vep_cache') else None
        annotator = VCFAnnotator(
            annotator=annotator_name,
            genome="auto",
            max_concurrency=5,
            timeout=30,
            vep_cache=vep_cache_path,
            interactive=False,
        )
        async def _annotate_and_close():
            try:
                annotated = await annotator.annotate(args.input)
                return annotated
            except Exception as e:  # noqa: BROAD_EXCEPT — process-level guard around annotator
                print(f"[GPA] VCF annotation failed: {type(e).__name__}: {e}", file=sys.stderr)
                raise
        try:
            annotated = asyncio.run(_annotate_and_close())
        except Exception:  # noqa: BROAD_EXCEPT — outer guard to prevent asyncio.run crash from killing process
            # Graceful degradation: return empty list so pipeline can continue
            # or re-raise depending on user preference. For now, re-raise with context.
            print("[GPA] ERROR: VCF annotation failed. Check network/proxy settings.")
            raise

        # Disease-aware transcript selection
        selector = None
        if args.disease_description:
            selector = TranscriptSelector(
                tissue_profile=args.tissue,
                disease_description=args.disease_description,
            )

        variants_data = variants_from_vep_annotation(annotated, selector)
        print(f"[GPA] Annotation complete: {len(variants_data)} variant-gene entries from VCF")

    elif input_type == InputType.ANNOTATED_VCF:
        # v0.10.0: Parse annotated VCF (CSQ in INFO) using VCFParser
        print("[GPA] Annotated VCF detected — parsing CSQ fields...")
        variants_data = parse_annotated_vcf(args.input, sample_idx=0)
        print(f"[GPA] Parsed {len(variants_data)} variants from annotated VCF")

    else:
        # Default: CSV/TSV (existing behavior)
        with open(args.input, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t' if args.input.endswith('.tsv') else ',')
            for row in reader:
                variants_data.append(dict(row))

    # v0.7: Set LLM model env var if specified
    if args.llm_model:
        os.environ.setdefault("GPA_LLM_MODEL", args.llm_model)

    # Run async pipeline with tissue context
    filter_stats = None
    if args.filter_stats:
        try:
            filter_stats = json.loads(args.filter_stats)
        except json.JSONDecodeError:
            print(f"[GPA] Warning: Invalid --filter-stats JSON, ignoring")

    config = GPAConfig(
        tissue_profile=args.tissue,
        offline_mode=args.offline,
        somatic_mode=args.somatic,
        target_population=args.target_population,
        multi_organ_profiles=multi_organ,
        force_sync=args.sync_gene_lists,
        evidence_detail=args.evidence_detail,
        database_version=args.database_version,
        filter_stats=filter_stats,
        filter_preset=args.filter_preset,
        # v0.8.0: SpliceAI
        spliceai_enabled=getattr(args, 'spliceai', False),
        spliceai_concurrency=getattr(args, 'spliceai_concurrency', 5),
        # v0.9.0: VCF annotation + transcript selection
        disease_description=args.disease_description,
        annotator=args.annotator,
        vep_cache=args.vep_cache,
        # v0.10.1: Two-phase pipeline
        two_phase=getattr(args, 'two_phase', False),
    )

    # v0.5 P2-3: Apply YAML config overrides to user config
    if file_config:
        file_config.apply_to_user_config(config)

    # v0.5 P2-3: Also build global config with file overrides (for API layer)
    global_config = config.to_global()
    if file_config:
        base_dir = config_path.parent if config_path else Path(__file__).parent.parent
        file_config.apply_to_global(global_config, base_dir)

    # v0.10.1: Preflight health check — verify all dependencies before starting analysis
    from gpa_preflight import run_preflight_check, suggest_action
    preflight, _route_map = asyncio.run(run_preflight_check(global_config))
    if not preflight.is_ready():
        action = suggest_action(preflight)
        if action == "abort":
            print("[GPA Preflight] 环境检查未通过，中止分析。")
            print(preflight.to_markdown())
            sys.exit(1)
        elif action == "offline":
            # v0.11.1: offline mode requires explicit user confirmation
            print("[GPA Preflight] 环境检查未通过，建议切换到离线模式（跳过所有 API 调用）。")
            print("  如需继续离线模式，请显式设置 --offline 参数后重试。")
            print("  当前默认行为：中止任务，保持在线优先。")
            sys.exit(1)

    # v0.5 P1-7: Multi-organ path
    if multi_organ:
        results = asyncio.run(run_multi_organ_assessment(variants_data, user_phenotypes=args.phenotypes, config=config))

        # Write joint report
        with open(args.output, 'w') as f:
            f.write(results["joint_report_markdown"])

        print(f"GPA Multi-Organ Report Generated: {args.output}")
        print(f"Profiles assessed: {', '.join(multi_organ)}")
        print(f"Joint Summary: Tier 1={results['joint_summary']['tier1_gene_count']} genes / {results['joint_summary']['tier1_variant_count']} variants, "
              f"Tier 2={results['joint_summary']['tier2_gene_count']} genes / {results['joint_summary']['tier2_variant_count']} variants, "
              f"Tier 3={results['joint_summary']['tier3_gene_count']} genes / {results['joint_summary']['tier3_variant_count']} variants")

        if results['joint_summary']['high_concern_variants']:
            genes = [v['gene'] for v in results['joint_summary']['high_concern_variants']]
            print(f"High-concern variants (Tier 1 in any profile): {', '.join(genes)}")

        # Write JSON if requested (backward compatible)
        if args.json:
            with open(args.json, 'w') as f:
                json.dump(results, f, indent=2, default=str)
            print(f"Structured output written to: {args.json}")

        # v0.5 P1-12: Write P1-12 structured JSON report if requested
        if args.output_json:
            with open(args.output_json, 'w') as f:
                json.dump(results.get("json_report", {}), f, indent=2, default=str, ensure_ascii=False)
            print(f"P1-12 JSON report written to: {args.output_json}")

        return

    # Single-organ path (original behavior)
    # v0.10.1: Two-phase pipeline for large VCF datasets
    # v0.10.2: Auto-enable two-phase for raw VCF unless explicitly overridden
    two_phase_explicit = getattr(args, 'two_phase', False)
    two_phase_enabled = two_phase_explicit or getattr(config, 'two_phase', False)
    if not two_phase_enabled and input_type == InputType.RAW_VCF:
        two_phase_enabled = True
        print("[GPA] Auto-enabling two-phase pipeline for raw VCF (fast local triage + API enrichment for candidates)")
    if two_phase_enabled:
        print("[GPA] Two-phase pipeline enabled — Phase 1: fast local triage, Phase 2: API enrichment for candidates only")
        from gpa_two_phase import run_two_phase_pipeline
        results = asyncio.run(run_two_phase_pipeline(variants_data, config=config, user_phenotypes=args.phenotypes, max_candidates=150))
    else:
        results = asyncio.run(run_dgra_pipeline(variants_data, user_phenotypes=args.phenotypes, config=config))

    # Write report
    with open(args.output, 'w') as f:
        f.write(results["report_markdown"])

    profile_name = results["meta"]["profile_display_name"]
    print(f"GPA Report Generated: {args.output}")
    print(f"Tissue Context: {profile_name} ({args.tissue})")
    print(f"Summary: Tier 1={results['summary']['tier1_gene_count']} genes / {results['summary']['tier1_variant_count']} variants, "
          f"Tier 2={results['summary']['tier2_gene_count']} genes / {results['summary']['tier2_variant_count']} variants, "
          f"Tier 3={results['summary']['tier3_gene_count']} genes / {results['summary']['tier3_variant_count']} variants")

    if results['summary']['multi_hit_genes']:
        print(f"Multi-hit genes: {', '.join(results['summary']['multi_hit_genes'])}")

    # Write JSON if requested (backward compatible)
    if args.json:
        with open(args.json, 'w') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"Structured output written to: {args.json}")

    # v0.5 P1-12: Write P1-12 structured JSON report if requested
    if args.output_json:
        with open(args.output_json, 'w') as f:
            json.dump(results.get("json_report", {}), f, indent=2, default=str, ensure_ascii=False)
        print(f"P1-12 JSON report written to: {args.output_json}")

# Backward-compatibility lazy re-exports to avoid circular imports
# (moved to gpa_tier_classifier.py / gpa_pipeline.py in v0.10.0)
def __getattr__(name):
    if name == "classify_variant_tier":
        from gpa_tier_classifier import classify_variant_tier as _cvt
        return _cvt
    if name == "run_dgra_pipeline":
        from gpa_pipeline import run_dgra_pipeline as _rgp
        return _rgp
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


if __name__ == "__main__":
    main()
