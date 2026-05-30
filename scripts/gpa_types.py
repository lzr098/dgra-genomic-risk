#!/usr/bin/env python3
"""
gpa_types.py - Shared types and constants for GPA pipeline.

Extracted from dgra_core.py to eliminate circular imports.
v0.10.11
"""

import json
from pathlib import Path
from dataclasses import dataclass, field
from typing import List, Dict, Optional, Tuple, Any

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
    "BARD1", "NBN",
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

    def to_global(self):
        from dgra_config import DGRAGlobalConfig as GPAGlobalConfig
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
        with open(ref_path, 'r', encoding='utf-8') as f:
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
        except (RuntimeError, ValueError) as e:
            # Non-blocking: if merge fails, keep static lists
            print(f"[GPA] Gene list sync warning: {e} - using static lists")

        return profile
