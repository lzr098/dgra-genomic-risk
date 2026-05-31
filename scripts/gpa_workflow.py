#!/usr/bin/env python3
"""
GPA Workflow Definition (v0.11.0)
=================================
Workflow-as-Code: declarative, executable pipeline steps for GPA analysis.

This file defines the canonical workflow that ALL GPA executions must follow.
It is the single source of truth for:
  - Which APIs/modules are called
  - In what order
  - Under what conditions steps may be skipped
  - What happens on failure

MODIFICATION POLICY:
  - This file may ONLY be modified in "optimize" mode.
  - Every modification must be reviewed and confirmed by the user.
  - NEVER modify this file during task execution ("run" mode).

Usage:
    from gpa_workflow import STANDARD_WORKFLOW, WorkflowStep
    for step in STANDARD_WORKFLOW:
        print(step.name, step.required)
"""

from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any
from enum import Enum


class FailureAction(Enum):
    """What to do when a workflow step fails."""
    ABORT = "abort"           # Stop the entire pipeline
    WARN = "warn"             # Log warning, continue
    FALLBACK = "fallback"     # Use fallback logic, continue
    SKIP = "skip"             # Skip this step, continue


@dataclass(frozen=True)
class WorkflowStep:
    """
    A single step in the GPA analysis pipeline.

    Attributes:
        name: Human-readable step name (snake_case, used as key)
        module: Python module that implements this step
        function: Function name within the module (optional)
        api_name: External API being called (e.g., "ensembl", "gnomAD").
                  None for internal processing steps.
        required: If True, this step MUST execute. Skipping is an error.
        skip_condition: Human-readable condition under which this step may be
                        automatically skipped. Empty string = never skip.
        timeout_sec: Maximum time allowed for this step (0 = no limit)
        on_failure: Action to take if this step fails
        description: Brief explanation of what this step does
        produces: List of data fields this step produces
        consumes: List of data fields this step requires as input
    """
    name: str
    module: str
    function: Optional[str] = None
    api_name: Optional[str] = None
    required: bool = True
    skip_condition: str = ""
    timeout_sec: int = 0
    on_failure: FailureAction = FailureAction.ABORT
    description: str = ""
    produces: List[str] = field(default_factory=list)
    consumes: List[str] = field(default_factory=list)

    def can_be_skipped(self, context: Dict[str, Any]) -> tuple[bool, str]:
        """
        Determine if this step should be skipped based on execution context.

        Returns:
            (should_skip, reason) — reason is empty if should_skip is False
        """
        if self.required:
            return False, ""

        # Evaluate skip conditions based on context
        if self.name == "gtex_expression":
            tissue_profile = context.get("tissue_profile", {})
            if not tissue_profile.get("gtex_tissue"):
                return True, "No GTEx tissue mapping configured for this profile"

        elif self.name == "myvariant_enrichment":
            variants = context.get("variants", [])
            needing = [v for v in variants if v.get("gnomAD_AF") is None]
            if not needing:
                return True, "All variants already have gnomAD AF data"

        elif self.name == "clinvar_ncbi":
            variants = context.get("variants", [])
            # v0.11.0: ClinVar is now required — only skip if truly no variants to query
            needing_clinvar = [v for v in variants if v.get("CLIN_SIG") in ("", "UNKNOWN", "_UNKNOWN")]
            if not needing_clinvar:
                return True, "No variants need ClinVar annotation (all have known status)"

        elif self.name == "spliceai_prediction":
            variants = context.get("variants", [])
            # v0.11.0: SpliceAI is required for critical positions
            splice_variants = [v for v in variants if any(
                term in v.get("Consequence", "").lower()
                for term in ["splice", "donor", "acceptor", "intron"]
            )]
            if not splice_variants:
                return True, "No splice-related or critical region variants to analyze"

        elif self.name == "vep_annotation":
            input_type = context.get("input_type", "")
            if input_type != "RAW_VCF":
                return True, f"Input type is {input_type}, not RAW_VCF — VEP annotation not needed"

        elif self.name == "transcript_selection":
            input_type = context.get("input_type", "")
            if input_type != "RAW_VCF":
                return True, f"Input type is {input_type} — transcript selection already done"

        elif self.name == "vep_reannotation":
            discrepancy_count = context.get("discrepancy_count", 0)
            if discrepancy_count == 0:
                return True, "No TRANSCRIPT_DISCREPANCY variants to reannotate"
            if context.get("offline_mode", False):
                return True, "Offline mode — VEP reannotation requires API access"

        elif self.name == "phenotype_matching":
            if not context.get("user_phenotypes"):
                return True, "No user phenotypes provided"

        elif self.name == "gnomad_variant_frequency":
            variants = context.get("variants", [])
            # v0.11.0: gnomAD variant frequency is now required
            without_af = [v for v in variants if v.get("gnomAD_AF") is None]
            if not without_af:
                return True, "All variants already have gnomAD AF data (fresh query not needed)"

        return False, ""


# =============================================================================
# STANDARD WORKFLOW — Canonical GPA Analysis Pipeline
# =============================================================================
# This is the single source of truth. All executions must follow this order.
# Steps marked required=False MAY be skipped with user notification.
# Steps marked required=True MUST NOT be skipped.

STANDARD_WORKFLOW: List[WorkflowStep] = [
    # -------------------------------------------------------------------------
    # Phase 0: Environment & Validation
    # -------------------------------------------------------------------------
    WorkflowStep(
        name="preflight_check",
        module="gpa_preflight",
        function="run_preflight_check",
        required=True,
        timeout_sec=60,
        on_failure=FailureAction.ABORT,
        description="Verify all API dependencies, proxy routes, and network connectivity before analysis",
        produces=["preflight_report", "proxy_route_map"],
    ),

    # -------------------------------------------------------------------------
    # Phase 1: Input Processing & VEP Annotation (for raw VCF)
    # -------------------------------------------------------------------------
    WorkflowStep(
        name="input_parsing",
        module="dgra_input_parsers",
        function="parse_input",
        required=True,
        timeout_sec=120,
        on_failure=FailureAction.ABORT,
        description="Parse input file (VCF, TSV, CSV, Excel) and detect format",
        produces=["variants", "input_type", "input_stats"],
    ),

    WorkflowStep(
        name="vep_annotation",
        module="gpa_vcf_annotator",
        function="annotate",
        api_name="ensembl",
        required=False,
        skip_condition="Input is not RAW_VCF (already annotated)",
        timeout_sec=1800,  # 30 min for large VCFs
        on_failure=FailureAction.WARN,
        description="Annotate raw VCF via Ensembl VEP REST API (batch processing with checkpoint)",
        produces=["annotated_variants", "vep_stats"],
        consumes=["raw_vcf_path"],
    ),

    WorkflowStep(
        name="transcript_selection",
        module="gpa_transcript_selector",
        function="select_transcripts",
        required=False,
        skip_condition="Input is not RAW_VCF",
        timeout_sec=300,
        on_failure=FailureAction.WARN,
        description="Disease-aware transcript selection from VEP output (canonical/MANE/LLM)",
        produces=["selected_transcripts", "transcript_ambiguity_flags"],
        consumes=["annotated_variants"],
    ),

    # -------------------------------------------------------------------------
    # Phase 2: Gene-Level Batch API Queries
    # -------------------------------------------------------------------------
    WorkflowStep(
        name="ensembl_gene_batch",
        module="dgra_api",
        function="batch_query_genes",
        api_name="ensembl",
        required=True,
        timeout_sec=300,
        on_failure=FailureAction.FALLBACK,
        description="Batch query Ensembl for gene info, transcripts, and functional data",
        produces=["ensembl_data"],
        consumes=["unique_genes"],
    ),

    WorkflowStep(
        name="uniprot_gene_batch",
        module="dgra_api",
        function="batch_query_genes",
        api_name="uniprot",
        required=True,
        timeout_sec=300,
        on_failure=FailureAction.FALLBACK,
        description="Batch query UniProt for protein domains and functional annotations",
        produces=["uniprot_data"],
        consumes=["unique_genes"],
    ),

    WorkflowStep(
        name="hgnc_gene_batch",
        module="dgra_api",
        function="batch_query_genes",
        api_name="hgnc",
        required=True,
        timeout_sec=300,
        on_failure=FailureAction.FALLBACK,
        description="Batch query HGNC for gene symbol normalization and aliases",
        produces=["hgnc_data"],
        consumes=["unique_genes"],
    ),

    WorkflowStep(
        name="gnomad_constraint_batch",
        module="dgra_api",
        function="batch_query_genes",
        api_name="gnomad_constraint",
        required=True,
        timeout_sec=300,
        on_failure=FailureAction.FALLBACK,
        description="Batch query gnomAD for gene constraint metrics (pLI, LOEUF)",
        produces=["gnomad_constraint_data"],
        consumes=["unique_genes"],
    ),

    WorkflowStep(
        name="gtex_expression",
        module="dgra_api",
        function="batch_query_genes",
        api_name="gtex",
        required=True,
        skip_condition="No GTEx tissue mapping configured",
        timeout_sec=600,
        on_failure=FailureAction.FALLBACK,
        description="Query GTEx for tissue-specific gene expression levels (required for transcript selection confirmation)",
        produces=["gtex_data"],
        consumes=["unique_genes", "tissue_profile"],
    ),

    # -------------------------------------------------------------------------
    # Phase 3: Variant-Level Enrichment
    # -------------------------------------------------------------------------
    WorkflowStep(
        name="myvariant_enrichment",
        module="dgra_myvariant",
        function="query_myvariant_batch",
        api_name="myvariant",
        required=False,
        skip_condition="All variants already have gnomAD AF data",
        timeout_sec=300,
        on_failure=FailureAction.WARN,
        description="Batch query MyVariant.info for gnomAD/CADD/ClinVar aggregation",
        produces=["myvariant_results"],
        consumes=["variants_needing_af"],
    ),

    WorkflowStep(
        name="gnomad_variant_frequency",
        module="dgra_api",
        function="query_gnomad_variant",
        api_name="gnomad",
        required=True,
        skip_condition="No variants to query",
        timeout_sec=600,
        on_failure=FailureAction.WARN,
        description="Query gnomAD GraphQL for variant-level allele frequencies (required for population frequency assessment)",
        produces=["gnomad_variant_results"],
        consumes=["variants_without_af"],
    ),

    WorkflowStep(
        name="clinvar_ncbi",
        module="dgra_clinvar",
        function="query_clinvar_batch",
        api_name="ncbi_eutils",
        required=True,
        skip_condition="No variants to query",
        timeout_sec=1200,  # 1 req/sec, up to 1000 variants = ~17 min
        on_failure=FailureAction.WARN,
        description="Direct NCBI E-utilities query for accurate ClinVar annotations (required for clinical interpretation)",
        produces=["clinvar_results"],
        consumes=["variants_needing_clinvar"],
    ),

    # -------------------------------------------------------------------------
    # Phase 4: Post-Processing
    # -------------------------------------------------------------------------
    WorkflowStep(
        name="hgnc_normalization",
        module="dgra_core",
        function="normalize_gene_symbols",
        required=True,
        timeout_sec=60,
        on_failure=FailureAction.WARN,
        description="Normalize gene symbols using HGNC data (handle aliases, outdated symbols)",
        produces=["normalized_variants", "hgnc_warnings"],
        consumes=["variants", "hgnc_data"],
    ),

    WorkflowStep(
        name="gene_constraint_population",
        module="dgra_core",
        function="evaluate_gene_constraint",
        required=True,
        timeout_sec=60,
        on_failure=FailureAction.WARN,
        description="Populate gene constraint metrics (pLI, LOEUF) onto variants",
        produces=["variants_with_constraint"],
        consumes=["variants", "gnomad_constraint_data"],
    ),

    WorkflowStep(
        name="nmd_prediction",
        module="dgra_core",
        function="predict_nmd",
        required=True,
        timeout_sec=120,
        on_failure=FailureAction.WARN,
        description="Predict nonsense-mediated decay for truncating variants",
        produces=["variants_with_nmd"],
        consumes=["variants", "ensembl_data"],
    ),

    WorkflowStep(
        name="transcript_correction",
        module="dgra_core",
        function="correct_transcript_priority",
        required=True,
        timeout_sec=300,
        on_failure=FailureAction.WARN,
        description="Correct transcript priority using Ensembl canonical/MANE data",
        produces=["variants_with_corrected_tx", "transcript_warnings"],
        consumes=["variants", "ensembl_data"],
    ),

    WorkflowStep(
        name="vep_reannotation",
        module="dgra_api",
        function="batch_query_vep_region",
        api_name="ensembl",
        required=False,
        skip_condition="No TRANSCRIPT_DISCREPANCY variants OR offline mode",
        timeout_sec=300,
        on_failure=FailureAction.WARN,
        description="Reannotate variants with TRANSCRIPT_DISCREPANCY using canonical transcript",
        produces=["vep_reannotation_results"],
        consumes=["discrepancy_variants"],
    ),

    WorkflowStep(
        name="tissue_relevance",
        module="dgra_core",
        function="assess_tissue_relevance",
        required=True,
        timeout_sec=120,
        on_failure=FailureAction.WARN,
        description="Assess tissue relevance for each variant using GTEx expression data",
        produces=["variants_with_tissue_relevance"],
        consumes=["variants", "gtex_data", "tissue_profile"],
    ),

    WorkflowStep(
        name="spliceai_prediction",
        module="dgra_splice_predictor",
        function="query_spliceai_batch",
        api_name="spliceai",
        required=True,
        skip_condition="No splice-related or critical region variants to analyze",
        timeout_sec=300,
        on_failure=FailureAction.WARN,
        description="SpliceAI prediction for splice donor/acceptor/region variants (required for critical positions)",
        produces=["spliceai_results"],
        consumes=["splice_variants"],
    ),

    # -------------------------------------------------------------------------
    # Phase 5: Classification & Reporting
    # -------------------------------------------------------------------------
    WorkflowStep(
        name="tier_classification",
        module="gpa_tier_classifier",
        function="classify_variant_tier",
        required=True,
        timeout_sec=300,
        on_failure=FailureAction.ABORT,
        description="Three-tier risk classification (Tier 1/2/3) with confidence scoring",
        produces=["classified_variants", "tier_summary"],
        consumes=["variants"],
    ),

    WorkflowStep(
        name="multi_hit_detection",
        module="gpa_multi_hit",
        function="detect_multi_hit_genes",
        required=True,
        timeout_sec=120,
        on_failure=FailureAction.WARN,
        description="Detect genes with multiple hits and assess cis/trans phase",
        produces=["multi_hit_genes", "phase_results"],
        consumes=["classified_variants"],
    ),

    WorkflowStep(
        name="phenotype_matching",
        module="gpa_phenotype_match",
        function="match_phenotypes",
        required=False,
        skip_condition="No user phenotypes provided",
        timeout_sec=120,
        on_failure=FailureAction.WARN,
        description="LLM-based phenotype-gene association matching",
        produces=["phenotype_matches"],
        consumes=["tier12_variants", "user_phenotypes"],
    ),

    WorkflowStep(
        name="qc_checks",
        module="gpa_qc",
        function="_run_qc_checks",
        required=True,
        timeout_sec=60,
        on_failure=FailureAction.WARN,
        description="Quality control checks (annotation coverage, conflicting ClinVar, etc.)",
        produces=["qc_flags", "qc_warnings"],
        consumes=["classified_variants"],
    ),

    WorkflowStep(
        name="report_generation",
        module="gpa_report",
        function="generate_tier_report",
        required=True,
        timeout_sec=120,
        on_failure=FailureAction.ABORT,
        description="Generate Markdown report with Tier 1/2/3 details and transcript validation",
        produces=["report_markdown", "report_json"],
        consumes=["classified_variants", "multi_hit_genes", "qc_flags"],
    ),
]


# =============================================================================
# Workflow Metadata
# =============================================================================

WORKFLOW_VERSION = "0.11.0-required"
WORKFLOW_NAME = "GPA Standard Analysis Pipeline"
WORKFLOW_DESCRIPTION = """
Canonical GPA analysis workflow. All steps must be executed in order.
Required steps cannot be skipped. Optional steps may be auto-skipped with
user notification if skip conditions are met.

v0.11.0-required changes:
- GTEx expression: required=True (for transcript selection confirmation)
- gnomAD variant frequency: required=True (population frequency assessment)
- ClinVar NCBI: required=True (clinical interpretation)
- SpliceAI prediction: required=True (critical position analysis)
"""


def get_step_by_name(name: str) -> Optional[WorkflowStep]:
    """Look up a workflow step by name."""
    for step in STANDARD_WORKFLOW:
        if step.name == name:
            return step
    return None


def get_required_steps() -> List[WorkflowStep]:
    """Return all required (non-skippable) steps."""
    return [s for s in STANDARD_WORKFLOW if s.required]


def get_optional_steps() -> List[WorkflowStep]:
    """Return all optional (skippable) steps."""
    return [s for s in STANDARD_WORKFLOW if not s.required]


def get_api_steps() -> List[WorkflowStep]:
    """Return all steps that call external APIs."""
    return [s for s in STANDARD_WORKFLOW if s.api_name is not None]


def validate_workflow() -> List[str]:
    """
    Validate the workflow definition for consistency.
    Returns list of error messages (empty if valid).
    """
    errors = []
    names = set()
    for i, step in enumerate(STANDARD_WORKFLOW):
        if step.name in names:
            errors.append(f"Duplicate step name: {step.name}")
        names.add(step.name)
        if not step.module:
            errors.append(f"Step {step.name}: module is required")
        # v0.11.0: Required steps MAY have skip_condition for "no input data" cases
        # (e.g., required ClinVar step skips when all variants already have ClinVar status)
        if not step.required and not step.skip_condition:
            errors.append(f"Step {step.name}: optional steps must have skip_condition")
        if step.timeout_sec < 0:
            errors.append(f"Step {step.name}: timeout must be >= 0")
    return errors


# Run validation on import
_validation_errors = validate_workflow()
if _validation_errors:
    import warnings
    for err in _validation_errors:
        warnings.warn(f"Workflow validation error: {err}", RuntimeWarning)
