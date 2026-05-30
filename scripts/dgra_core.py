#!/usr/bin/env python3
"""
GPA Core Engine - Genomic Phenotype Association

Orchestrator module. Types moved to gpa_types.py, analysis functions to gpa_analysis.py.
v0.10.11: Eliminated circular imports through module decomposition.
"""

import json
import sys
import csv
import asyncio
from pathlib import Path
import argparse

from version import __version__

from gpa_types import GPAConfig
from dgra_config import DGRAFileConfig as GPAFileConfig, DEFAULT_CONFIG_PATH

# v0.8.0: SpliceAI requires aiohttp for async HTTP queries
try:
    import aiohttp
except ImportError:
    aiohttp = None

def main():
    # v0.10.0: Lazy imports to avoid circular dependency with gpa_* modules
    from gpa_phaser import PhaseResult, determine_phase
    from gpa_multi_hit import detect_multi_hit_genes
    from gpa_report import _get_version_info, generate_tier_report, generate_json_report
    from gpa_qc import _run_qc_checks
    from gpa_input import InputType, detect_input_type, variants_from_vep_annotation
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
        except (RuntimeError, ValueError) as e:
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
        )
        async def _annotate_and_close():
            annotated = await annotator.annotate(args.input)
            await annotator.close()
            return annotated
        annotated = asyncio.run(_annotate_and_close())

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
        # v0.9.0: Parse annotated VCF (CSQ in INFO)
        print("[GPA] Annotated VCF detected — parsing CSQ fields...")
        # For v0.9.0, annotated VCF parsing is simplified; full support in future
        raise NotImplementedError("Annotated VCF input parsing not yet implemented in v0.9.0. "
                                  "Please convert to TSV/CSV or use raw VCF.")

    else:
        # Default: CSV/TSV (existing behavior)
        with open(args.input, 'r', encoding='utf-8') as f:
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

    # v0.5 P1-7: Multi-organ path
    if multi_organ:
        results = asyncio.run(run_multi_organ_assessment(variants_data, user_phenotypes=args.phenotypes, config=config))

        # Write joint report
        with open(args.output, 'w', encoding='utf-8') as f:
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
            with open(args.json, 'w', encoding='utf-8') as f:
                json.dump(results, f, indent=2, default=str)
            print(f"Structured output written to: {args.json}")

        # v0.5 P1-12: Write P1-12 structured JSON report if requested
        if args.output_json:
            with open(args.output_json, 'w', encoding='utf-8') as f:
                json.dump(results.get("json_report", {}), f, indent=2, default=str, ensure_ascii=False)
            print(f"P1-12 JSON report written to: {args.output_json}")

        return

    # Single-organ path (original behavior)
    # v0.10.1: Two-phase pipeline for large VCF datasets
    if getattr(args, 'two_phase', False) or getattr(config, 'two_phase', False):
        print("[GPA] Two-phase pipeline enabled — Phase 1: fast local triage, Phase 2: API enrichment for candidates only")
        from gpa_two_phase import run_two_phase_pipeline
        results = asyncio.run(run_two_phase_pipeline(variants_data, config=config, user_phenotypes=args.phenotypes, max_candidates=150))
    else:
        results = asyncio.run(run_dgra_pipeline(variants_data, user_phenotypes=args.phenotypes, config=config))

    # Write report
    with open(args.output, 'w', encoding='utf-8') as f:
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
        with open(args.json, 'w', encoding='utf-8') as f:
            json.dump(results, f, indent=2, default=str)
        print(f"Structured output written to: {args.json}")

    # v0.5 P1-12: Write P1-12 structured JSON report if requested
    if args.output_json:
        with open(args.output_json, 'w', encoding='utf-8') as f:
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
