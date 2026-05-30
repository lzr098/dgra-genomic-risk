#!/usr/bin/env python3
"""
GPA Multi-hit Gene Detection Module

Detects genes with multiple potentially pathogenic variants that may
constitute compound heterozygosity. Includes pairwise phase analysis
for close variant pairs.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
"""

from __future__ import annotations

from typing import List, Dict, Optional, TYPE_CHECKING

from gpa_types import Variant
from gpa_analysis import (
    _variant_has_pathogenic_evidence,
)


def detect_multi_hit_genes(variants: List[Variant], gtex_data: Optional[Dict] = None) -> List[Dict]:
    """
    Detect genes with multiple pathogenic variants that may require phase analysis.

    v0.9.4 FIX: Added pairwise phase analysis for close variant pairs.
    Previously, all pathogenic variants in a gene were passed to determine_phase(),
    causing distant benign variants to inflate max_gap_bp and incorrectly mark
    close, clinically-relevant pairs as "infeasible_short_reads".

    Now: Pairwise analysis for pairs <= 1000bp in addition to whole-gene assessment.

    Only counts variants with evidence of pathogenicity:
      - Domain impact, or
      - ClinVar pathogenic / HIGH impact / rare gnomAD, or
      - Splice site change

    Normal polymorphisms (benign / common / no domain impact) are excluded.
    """

    # v0.10.11: Lazy import to avoid circular dependency
    from gpa_phaser import determine_phase

    # Group variants by gene
    gene_variants = {}
    for v in variants:
        if not v.gene or v.gene.strip() == "":
            continue
        gene_variants.setdefault(v.gene, []).append(v)

    multi_hits = []
    for gene, var_list in gene_variants.items():
        # Count only variants with pathogenic evidence
        pathogenic_vars = [v for v in var_list if _variant_has_pathogenic_evidence(v, gtex_data)]

        if len(pathogenic_vars) >= 2:
            # v0.4.5: Whole-gene phase analysis (kept for overall assessment)
            phase_result = determine_phase(pathogenic_vars)

            # v0.9.4 FIX: Pairwise phase analysis for close variant pairs
            pairwise_phases = []
            for i in range(len(pathogenic_vars)):
                for j in range(i + 1, len(pathogenic_vars)):
                    v1, v2 = pathogenic_vars[i], pathogenic_vars[j]
                    dist = abs(v1.pos - v2.pos)
                    if dist <= 1000:
                        pair_phase = determine_phase([v1, v2])
                        pairwise_phases.append({
                            "variant1": {"chrom": v1.chrom, "pos": v1.pos, "hgvsp": v1.hgvsp, "impact": v1.impact},
                            "variant2": {"chrom": v2.chrom, "pos": v2.pos, "hgvsp": v2.hgvsp, "impact": v2.impact},
                            "distance_bp": dist,
                            "phase_result": {
                                "status": pair_phase.phase_status,
                                "confidence": pair_phase.confidence,
                                "method": pair_phase.method,
                                "evidence": pair_phase.evidence,
                            }
                        })

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

            # v0.4.5: 相位状态临床解读
            phase_clinical = {
                "cis": "两个变异位于同一单倍型 → 另一单倍型正常 → 保留 50% 功能",
                "trans": "两个变异位于不同单倍型 → 复合杂合 → 功能可能完全丧失",
                "cis_both": "两条单倍型均携带变异 → 纯合/复合 → 功能严重受损",
                "ambiguous": "相位关系不确定 → 需进一步验证",
                "unphased": "超出短 reads 相位范围 → 需 trio 或长读长",
                "cis_likely": "高概率 cis,但未 100% 确认 → 建议验证",
                "trans_likely": "高概率 trans,但未 100% 确认 → 建议验证"
            }

            multi_hits.append({
                "gene": gene,
                "variant_count": len(var_list),           # total variants in gene
                "pathogenic_count": len(pathogenic_vars),  # variants with evidence
                "warning": "MULTI_HIT_GENE",
                "pathogenic_variants": var_details,
                "phase_result": {
                    "status": phase_result.phase_status,
                    "confidence": phase_result.confidence,
                    "method": phase_result.method,
                    "evidence": phase_result.evidence,
                    "max_gap_bp": phase_result.max_gap_bp,
                    "min_gap_bp": phase_result.min_gap_bp,
                    "n_variants": phase_result.n_variants
                },
                "phase_clinical_significance": phase_clinical.get(phase_result.phase_status, "未知"),
                # v0.9.4 FIX: Add pairwise phase analysis
                "pairwise_phase_analysis": pairwise_phases,
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
