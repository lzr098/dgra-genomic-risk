#!/usr/bin/env python3
"""
GPA Tier Classification Module

Three-tier variant classification with dynamic tissue context,
ACMG-like evidence chains, and upgrade condition generation.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
"""

import json
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
from gpa_types import (
    Variant,
    GPAConfig,
    Evidence,
    _UNKNOWN,
    _is_unknown,
    _COMMON_TS_GENES,
    _KNOWN_AML_DRIVERS,
    _GENE_FAMILY_REDUNDANCY,
)
from gpa_analysis import (
    evaluate_gene_constraint,
    evaluate_missense_tier,
    predict_nmd,
)




# v0.7 Phase 3: Rare disease gene list (from gene_phenotype_map.json)
_RARE_DISEASE_GENES: Optional[set] = None


def _load_rare_disease_genes() -> set:
    """Load rare disease gene list from gene_phenotype_map.json.
    Genes with OMIM/ClinVar phenotypes are considered rare disease-related.
    """
    global _RARE_DISEASE_GENES
    if _RARE_DISEASE_GENES is not None:
        return _RARE_DISEASE_GENES

    map_path = Path(__file__).parent.parent / "references" / "gene_phenotype_map.json"
    try:
        with open(map_path, 'r', encoding='utf-8') as f:
            data = json.load(f)
        _RARE_DISEASE_GENES = set(data.keys())
        return _RARE_DISEASE_GENES
    except (FileNotFoundError, IsADirectoryError, PermissionError, ValueError, json.JSONDecodeError):
        _RARE_DISEASE_GENES = set()
        return _RARE_DISEASE_GENES


def _is_rare_disease_gene(gene: str) -> bool:
    """Check if gene is in the rare disease gene list."""
    return gene in _load_rare_disease_genes()


def classify_variant_tier(variant: Variant, domain_info: Dict, tissue_assessment: Dict,
                          gnomad_info: Dict, transcript_warning: Optional[Dict],
                          pseudogene_warning: Optional[Dict], tissue_profile: Dict,
                          config: Optional[GPAConfig] = None) -> Tuple[int, str, List[str]]:
    """
    Three-tier classification with dynamic tissue context.
    v0.4.5: Added somatic_mode support for tumor driver analysis.
    Returns: (tier, reason, actions)
    """
    gene = variant.gene
    actions = []
    # v0.5 P1-9: Initialize evidence chain for structured traceability
    # If variant already has evidence (e.g., from previous analysis or testing), preserve it
    evidence_chain = list(variant.evidence_chain)
    def _add_evidence(source, rule, weight=1.0, confidence="high", raw_data=None):
        evidence_chain.append(Evidence(source=source, rule=rule, weight=weight, confidence=confidence, raw_data=raw_data))

    def _clinvar_is_conflicting(clinvar):
        """Detect ClinVar conflicting interpretations.
        v0.7.1: Conflicting = pathogenic AND benign/VUS keywords present simultaneously.
        Standard composite ratings like 'Pathogenic/Likely_pathogenic' are NOT conflicting.
        """
        if _is_unknown(clinvar):
            return False
        clinvar_lower = clinvar.lower()
        if "conflicting" in clinvar_lower:
            return True
        pathogenic_keywords = ["pathogenic", "致病", "likely_pathogenic", "可能致病"]
        benign_or_vus_keywords = ["benign", "良性", "likely_benign", "可能良性", "vus", "意义不明", "uncertain"]
        has_pathogenic = any(kw in clinvar_lower for kw in pathogenic_keywords)
        has_benign_or_vus = any(kw in clinvar_lower for kw in benign_or_vus_keywords)
        if has_pathogenic and has_benign_or_vus:
            return True
        if "/" in clinvar and clinvar.count("/") == 1:
            return False
        return False

    def _clinvar_pathogenic(clinvar):
        """ClinVar pathogenic check - UNKNOWN does NOT trigger this.
        v0.5.2: Support both English 'Pathogenic' and Chinese '致病'.
        v0.7.1: Conflicting interpretations return False."""
        if _clinvar_is_conflicting(clinvar):
            return False
        if _is_unknown(clinvar):
            return False
        clinvar_lower = clinvar.lower()
        return ("pathogenic" in clinvar_lower or
                "致病" in clinvar or
                "likely_pathogenic" in clinvar_lower or
                "可能致病" in clinvar)

    def _clinvar_benign(clinvar):
        """ClinVar benign check - UNKNOWN does NOT trigger this.
        v0.5.2: Support both English 'Benign' and Chinese '良性'.
        v0.7.1: Conflicting interpretations return False."""
        if _clinvar_is_conflicting(clinvar):
            return False
        if _is_unknown(clinvar):
            return False
        clinvar_lower = clinvar.lower()
        return (("benign" in clinvar_lower or "良性" in clinvar)
                and "conflicting" not in clinvar_lower)

    def _parse_clinvar_confidence(clnrevstat):
        """v0.7.2: Map ClinVar CLNREVSTAT review status to confidence weight (0.30~0.95).

        ClinVar review status is text, not numeric stars:
            - practice_guideline → 0.95 (★★★★)
            - reviewed_by_expert_panel → 0.80 (★★★☆)
            - multiple_submitters_no_conflict → 0.55 (★★☆☆)
            - single_submitter → 0.40 (★☆☆☆)
            - no_assertion / no_criteria / conflicting / missing → 0.30
        """
        if not clnrevstat:
            return 0.30
        cs = clnrevstat.lower()
        if "practice_guideline" in cs:
            return 0.95
        if "reviewed_by_expert_panel" in cs:
            return 0.80
        if "multiple_submitters" in cs and "noconflict" in cs.replace("_", ""):
            return 0.55
        if "noconflicts" in cs.replace("_", ""):
            return 0.50
        if "single_submitter" in cs:
            return 0.40
        if "conflicting" in cs:
            return 0.30
        if "no_assertion" in cs or "no_criteria" in cs:
            return 0.30
        return 0.30

    # v0.7.2: Pre-compute ClinVar confidence once for this variant
    clinvar_conf = _parse_clinvar_confidence(variant.clinvar_review_status)
    if _clinvar_is_conflicting(variant.clinvar):
        if "CLINVAR_CONFLICTING" not in variant.qc_flags:
            variant.qc_flags.append("CLINVAR_CONFLICTING")
        _add_evidence(
            source="ClinVar",
            rule=f"Conflicting ClinVar interpretation: '{variant.clinvar}' - NOT used for tier upgrade",
            weight=0.0,
            confidence="low",
            raw_data={"clinvar": variant.clinvar},
        )
    if pseudogene_warning:
        pg_score = pseudogene_warning.get("score", 0)
        pg_level = pseudogene_warning.get("level", "unknown")
        pg_conf = "low" if pg_score >= 0.75 else ("moderate" if pg_score >= 0.40 else "high")
        _add_evidence(
            source="PseudogeneDetection",
            rule=f"Pseudogene_{pg_level}_score={pg_score:.2f}",
            weight=0.0,  # Does not affect tier directly
            confidence=pg_conf,
            raw_data=pseudogene_warning,
        )

    def _confidence_from_data():
        """Determine confidence based on data quality."""
        if variant.missing_fields:
            return "low" if len(variant.missing_fields) >= 3 else "moderate"
        if getattr(config, 'offline_mode', False):
            return "low"
        return "high"

    def _calculate_tier_confidence(chain):
        """v0.5 P1-10: Calculate tier confidence based on evidence chain.

        HIGH: >=3 independent high-confidence sources, no conflict.
        MEDIUM: 2 independent sources, or mixed quality, or minor conflict.
        LOW: <=1 source, or conflicting evidence, or many UNKNOWN fields.
        """
        if not chain:
            return "LOW"

        # Count unique high-confidence sources
        unique_sources = set()
        high_conf_count = 0
        has_conflict = False

        for ev in chain:
            if ev.confidence in ("high", "HIGH"):
                unique_sources.add(ev.source)
                high_conf_count += 1

        # Check for conflicting evidence
        clinvar_path = any("ClinVar" in ev.source and "Pathogenic" in ev.rule for ev in chain)
        gnomad_common = any(("gnomAD" in ev.source or "Frequency" in ev.source) and "common" in ev.rule.lower() for ev in chain)
        clinvar_benign = any("ClinVar" in ev.source and "Benign" in ev.rule for ev in chain)

        if (clinvar_path and gnomad_common) or (clinvar_path and clinvar_benign):
            has_conflict = True

        # Check for many UNKNOWN fields
        many_unknown = len(variant.missing_fields) >= 3 if variant.missing_fields else False

        if has_conflict or many_unknown:
            return "LOW"

        # v0.6: Pseudogene interference scoring - does NOT change tier, only confidence
        if pseudogene_warning:
            pg_score = pseudogene_warning.get("score", 0)
            if pg_score >= 0.75:
                return "LOW"  # Strong interference: confidence drops
            elif pg_score >= 0.40:
                return "MEDIUM"  # Suspected interference: confidence downgraded
            elif pg_score > 0:
                return "MEDIUM"  # Minor bias: slight downgrade
            # v0.5.3 legacy fallback (backward compatibility with old-format warnings)
            elif pseudogene_warning.get("type") == "PSEUDOGENE_INTERFERENCE":
                return "LOW"

        # v0.5.3: QC flags force LOW confidence
        if "VAF_GT_MISMATCH" in variant.qc_flags:
            return "LOW"

        if len(unique_sources) >= 3 and high_conf_count >= 3:
            return "HIGH"
        elif len(unique_sources) >= 2:
            return "MEDIUM"
        else:
            return "LOW"

    def _generate_upgrade_conditions(variant, tier, tissue_assessment, gnomad_info):
        """v0.5 P1-11: Generate forward-looking upgrade conditions based on evidence gaps.

        Tier 2 → Tier 1: What would make this variant actionable?
        Tier 3 → Tier 2: What would make this variant worth monitoring?
        Tier 1: No upgrade conditions (already highest tier).
        """
        conditions = []
        gene = variant.gene
        clinvar = variant.clinvar
        impact = variant.impact
        consequence = variant.consequence
        gnomad_status = gnomad_info.get("status", "")
        gnomad_af = gnomad_info.get("af")
        if gnomad_af is not None:
            try:
                gnomad_af = float(gnomad_af)
            except (ValueError, TypeError):
                gnomad_af = None

        if tier == 2:
            # Tier 2 → Tier 1 upgrade paths
            # Condition 1: ClinVar upgrade
            if clinvar and "Pathogenic" not in clinvar and "pathogenic" not in clinvar.lower():
                conditions.append(f"若 ClinVar 收录为 Pathogenic/Likely_pathogenic 则升级为 Tier 1")

            # Condition 2: Functional evidence
            if impact not in ("HIGH", ""):
                conditions.append(f"若功能实验证实 {consequence} 有害(如蛋白稳定性下降)则升级为 Tier 1")

            # Condition 3: Zygosity upgrade
            if variant.gt == "0/1":
                conditions.append(f"若后续验证为纯合变异 (1/1) 且基因对 {tissue_assessment.get('relevance', 'target')} 组织关键则升级为 Tier 1")

            # Condition 4: gnomAD AF near threshold
            if gnomad_af and gnomad_af > 0.001:
                conditions.append(f"若东亚人群 AF < 0.001% 或该位点在患者中富集则升级为 Tier 1")

            # Condition 5: Domain info upgrade
            if not variant.domain_info:
                conditions.append(f"若位于关键功能域或保守残基(如 ATP结合位点)则升级为 Tier 1")

        elif tier == 3:
            # Tier 3 → Tier 2 upgrade paths
            # Condition 1: de novo validation
            conditions.append(f"若后续家系验证为 de novo(非遗传)或患者表型与该基因高度匹配则升级为 Tier 2")

            # Condition 2: ClinVar upgrade from benign
            if clinvar and ("Benign" in clinvar or "benign" in clinvar.lower()):
                conditions.append(f"若 ClinVar 重新评级为 VUS 或以上,或新功能证据出现则升级为 Tier 2")

            # Condition 3: Common polymorphism but in special domain
            if gnomad_status == "common_polymorphism":
                conditions.append(f"若功能域分析显示该位点位于关键结构域或影响剪接则升级为 Tier 2")

            # Condition 4: Missing tissue relevance
            if tissue_assessment.get("relevance") == "none":
                conditions.append(f"若 GTEx 或其他数据显示该基因在 {tissue_assessment.get('relevance', 'target')} 组织中高表达则升级为 Tier 2")

            # Condition 5: Domain info gap
            if not variant.domain_info:
                conditions.append(f"若后续实验证实该变异影响蛋白功能或相互作用则升级为 Tier 2")

        return conditions

    profile_name = tissue_profile.get("display_name", "target tissue")
    tier_rules = tissue_profile.get("tier_rules", {})
    special_lists = tissue_profile.get("special_gene_lists", {})

    # v0.5 P0-7: Helpers for conservative UNKNOWN handling
    def _impact_high(impact):
        """Impact HIGH check - UNKNOWN is treated as HIGH (conservative, no downgrade).
        v0.5.2: TRANSCRIPT_DISCREPANCY with non-coding annotator transcript (NR_/XM_/XR_)
        but canonical protein-coding (ENST/ENSG) → downgrade HIGH to prevent false Tier 1.
        """
        if _is_unknown(impact):
            return True  # Missing impact data → assume worst case

        # v0.5.2: Transcript discrepancy check
        if (transcript_warning and
            transcript_warning.get("type") == "TRANSCRIPT_DISCREPANCY"):
            annotator_tx = transcript_warning.get("annotator_selected", "")
            canonical_tx = transcript_warning.get("canonical", "")
            is_annotator_noncoding = annotator_tx.startswith(("NR_", "XM_", "XR_"))
            is_canonical_protein = canonical_tx.startswith(("ENST", "ENSG"))
            if (is_annotator_noncoding and is_canonical_protein and
                impact == "HIGH"):
                actions.append(f"WARNING: Annotator used non-coding transcript {annotator_tx} "
                              f"but canonical {canonical_tx} is protein-coding. "
                              f"HIGH impact downgraded to MODERATE for tier classification.")
                _add_evidence("TranscriptWarning",
                    f"Non-coding annotator tx {annotator_tx} vs canonical {canonical_tx} → impact downgraded",
                    weight=0.5, confidence="medium",
                    raw_data={"annotator_tx": annotator_tx, "canonical_tx": canonical_tx})
                return False  # Treat as non-HIGH for tier classification

        return impact == "HIGH"

    # Report missing fields in actions
    if variant.missing_fields:
        actions.append(f"Missing fields: {', '.join(variant.missing_fields)} - conservative assessment applied")

    # v0.10.8: GTEx fast-track REMOVED from Phase 1 tier classification.
    # Previously, low GTEx TPM in the target profile could force Tier 3,
    # causing false negatives for tissue-specific genes (e.g. SAG in retina,
    # FOXE3 in eye) because GTEx v8 does not include eye tissues.
    # 
    # GTEx expression data is now used ONLY in Phase 2 for phenotype-tissue
    # association analysis, not as a hard tier gate.
    #
    # If tissue_assessment indicates low expression, we add it as soft
    # evidence but do NOT downgrade the tier here.
    if tissue_assessment.get("gtex_tpm") is not None and tissue_assessment.get("relevance") == "none":
        _add_evidence("TissueContext",
            f"Low GTEx expression in {profile_name} profile (TPM={tissue_assessment.get('gtex_tpm', 0):.2f}); "
            f"phenotype-tissue association assessed in Phase 2",
            weight=0.15, confidence="low",
            raw_data={"relevance": tissue_assessment.get("relevance"), "gtex_tpm": tissue_assessment.get("gtex_tpm")})

    # v0.4.5: Somatic mode overrides for tumor driver analysis
    # In somatic mode, tier classification prioritizes driver mutation evidence
    # over germline carrier-state logic
    if getattr(config, 'somatic_mode', False):
        # v0.9.0 fix: VAF > 0.5 check only meaningful in somatic context.
        # Germline homozygous variants (VAF≈1.0) are normal and MUST NOT be demoted.
        if hasattr(variant, 'vaf') and variant.vaf is not None and variant.vaf > 0.5:
            # Likely germline polymorphism contamination in somatic sample
            actions.append("VAF > 0.5 suggests germline contamination - verify if intended somatic analysis")
            variant.evidence_chain = evidence_chain
            upgrade_conditions = _generate_upgrade_conditions(variant, 3, tissue_assessment, gnomad_info)
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            return 3, f"VAF={variant.vaf:.3f} > 0.5 - likely germline polymorphism, not somatic driver", actions

        # Somatic mode Tier 1: Core driver mutations
        # 1a. TSG loss-of-function in tissue-relevant gene = Tier 1 (core driver)
        if tissue_assessment.get("relevance") in ["primary", "secondary"] and _impact_high(variant.impact):
            # Check if gene is known TSG (from OncoKB annotation or common TSG list)
            is_tsg = getattr(variant, 'is_tsg', False) or gene in _COMMON_TS_GENES
            if is_tsg:
                reason = f"Somatic TSG loss-of-function: {variant.consequence} in {gene}"
                if domain_info and domain_info.get("domain_integrity") in ["completely_destroyed", "partially_destroyed"]:
                    reason += f", {domain_info['domain']} domain disrupted"
                actions.append("Confirm somatic origin (VAF < 0.5, tumor-normal pair)")
                actions.append("Assess as core leukemic driver - target for MRD monitoring")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = []  # Tier 1: no upgrade
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 1, reason, actions

        # 1b. Oncogene hotspot / functional domain mutation = Tier 1
        is_oncogene = getattr(variant, 'is_oncogene', False)
        oncokb_class = getattr(variant, 'classification', '')
        if is_oncogene or oncokb_class in ("Oncogenic", "Likely Oncogenic"):
            if tissue_assessment.get("relevance") in ["primary", "secondary"]:
                reason = f"Somatic oncogene driver: {gene} {variant.hgvsp or variant.hgvsc}"
                if oncokb_class:
                    reason += f" (OncoKB: {oncokb_class})"
                actions.append("Confirm somatic origin")
                actions.append("Assess as core leukemic driver - potential therapeutic target")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = []  # Tier 1: no upgrade
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 1, reason, actions

        # 1c. Known AML driver genes with HIGH impact = Tier 1
        if gene in _KNOWN_AML_DRIVERS and _impact_high(variant.impact):
            actions.append("Known AML driver gene with truncating mutation")
            actions.append("Assess for therapeutic targeting or MRD monitoring")
            variant.evidence_chain = evidence_chain
            upgrade_conditions = []  # Tier 1: no upgrade
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            return 1, f"Known AML driver {gene} with {variant.consequence} - core somatic driver", actions

    # =====================================================================
    # v0.9.4 P0 FIX 2026-05-24: Population frequency guard — EAS AF > 50% auto-Tier 3
    # Prevents false Tier 1 calls for common polymorphisms (e.g., MAD2L2 rs2233004
    # with EAS AF=93.5%, OR2B11 p.Phe8SerfsTer2 with AF=48%).
    # 
    # Applied BEFORE Tier 1 checks. A variant with population frequency >50% in
    # the target population (default: EAS) is extremely unlikely to be pathogenic.
    # 
    # Exceptions:
    # - Autosomal recessive carrier screening: homozygous common LOF in AR genes
    #   may still be relevant if ClinVar pathogenic AND rare disease gene
    # - Pharmacogenomic variants: common PGx alleles are not pathogenic but clinically
    #   actionable (handled by separate PGx module)
    # =====================================================================
    _EAS_COMMON_AF_THRESHOLD = 0.50   # EAS AF > 50% → Tier 3
    _GLOBAL_COMMON_AF_THRESHOLD = 0.80  # Global AF > 80% → Tier 3 regardless of population
    
    gnomad_af = gnomad_info.get("af")
    try:
        if gnomad_af is not None:
            gnomad_af = float(gnomad_af)
    except (ValueError, TypeError):
        gnomad_af = None
    gnomad_af_populations = gnomad_info.get("af_populations", {})
    
    # Check EAS population frequency first (most relevant for Chinese/Asian cohorts)
    eas_af = None
    if gnomad_af_populations:
        for pop_code in ("EAS", "eas"):
            if pop_code in gnomad_af_populations:
                eas_af = gnomad_af_populations[pop_code].get("af")
                if eas_af is not None:
                    try:
                        eas_af = float(eas_af)
                    except (ValueError, TypeError):
                        eas_af = None
                break
    
    # Determine if this is a common polymorphism that should be auto-Tier 3
    force_tier3_frequency = False
    frequency_override_reason = ""
    
    if eas_af is not None and eas_af > _EAS_COMMON_AF_THRESHOLD:
        force_tier3_frequency = True
        frequency_override_reason = (
            f"EAS population frequency {eas_af:.1%} > {_EAS_COMMON_AF_THRESHOLD:.0%} threshold "
            f"— extremely common in Asian populations, cannot be pathogenic"
        )
    elif gnomad_af is not None and gnomad_af > _GLOBAL_COMMON_AF_THRESHOLD:
        # Global AF > 80% but EAS may not be available
        force_tier3_frequency = True
        frequency_override_reason = (
            f"Global AF {gnomad_af:.1%} > {_GLOBAL_COMMON_AF_THRESHOLD:.0%} — "
            f"near-fixation variant, cannot be pathogenic"
        )
    
    if force_tier3_frequency:
        _add_evidence("Frequency", frequency_override_reason,
                     weight=-1.0, confidence="high",
                     raw_data={"gnomad_af": gnomad_af, "eas_af": eas_af,
                              "gnomad_status": gnomad_info.get("status")})
        actions.append(f"⚠️ {frequency_override_reason}")
        actions.append("Auto Tier 3 — population frequency incompatible with monogenic disease")
        # Clear any pre-existing ClinVar pathogenic evidence since it's overridden by frequency
        if _clinvar_pathogenic(variant.clinvar):
            actions.append(f"ClinVar pathogenic call ({variant.clinvar}) overridden by population frequency")
        variant.evidence_chain = evidence_chain
        upgrade_conditions = []  # Tier 3: no upgrade
        variant.upgrade_conditions = upgrade_conditions
        variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
        variant.qc_flags.append("POPULATION_FREQUENCY_OVERRIDE")
        return 3, f"Auto Tier 3: {frequency_override_reason}", actions

    # =====================================================================
    # v0.9.4 P1 FIX 2026-05-25: Moderate frequency guard — AF > 1% → Tier 3
    # Post-mortem analysis revealed that many Tier S/Tier 1 candidates were
    # actually common benign polymorphisms missed because VEP API calls lacked
    # gnomAD/check_existing parameters.
    #
    # A variant at >1% AF in any population is very unlikely to be a high-
    # penetrance monogenic disease mutation (ACMG PM2/BS1). 
    #
    # Exception: ClinVar Pathogenic/Likely_pathogenic with expert_panel or
    # practice_guideline review status — these known founder/polymorphic
    # variants can be pathogenic despite higher AF.
    # =====================================================================
    _MODERATE_AF_THRESHOLD = 0.01  # AF > 1% → Tier 3 (benign polymorphism)

    moderate_af_exceeded = False
    moderate_af_value = None
    moderate_af_pop = ""

    if eas_af is not None and eas_af > _MODERATE_AF_THRESHOLD:
        moderate_af_exceeded = True
        moderate_af_value = eas_af
        moderate_af_pop = "EAS"
    elif gnomad_af is not None and gnomad_af > _MODERATE_AF_THRESHOLD:
        moderate_af_exceeded = True
        moderate_af_value = gnomad_af
        moderate_af_pop = "global"

    if moderate_af_exceeded:
        # Check for ClinVar exception
        clinvar_conf = _parse_clinvar_confidence(variant.clinvar_review_status)
        has_high_confidence_clinvar = (
            _clinvar_pathogenic(variant.clinvar) and clinvar_conf >= 0.80
        )

        if not has_high_confidence_clinvar:
            reason = (
                f"{moderate_af_pop} AF {moderate_af_value:.4f} > {_MODERATE_AF_THRESHOLD:.0%} "
                f"— common benign polymorphism"
            )
            _add_evidence("Frequency", reason,
                         weight=-0.8, confidence="high",
                         raw_data={"gnomad_af": gnomad_af, "eas_af": eas_af,
                                  "threshold": _MODERATE_AF_THRESHOLD})
            if variant.clinvar and not _is_unknown(variant.clinvar):
                actions.append(f"ClinVar {variant.clinvar} overridden by population frequency >1%")
            actions.append(f"Auto Tier 3: {reason}")
            variant.evidence_chain = evidence_chain
            variant.upgrade_conditions = []
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            variant.qc_flags.append("BENIGN_POLYMORPHISM_FREQUENCY")
            return 3, f"Auto Tier 3: {reason}", actions
        else:
            actions.append(
                f"⚠️ {moderate_af_pop} AF {moderate_af_value:.4f} > 1% but ClinVar "
                f"Pathogenic (review: {variant.clinvar_review_status}, "
                f"confidence: {clinvar_conf:.0%}) — keeping for review"
            )

    # Priority 1: Tier 1 checks (germline disease risk logic)
    # 1a. Known high-risk special gene lists with pathogenic variant
    for list_name, gene_list in special_lists.items():
        if gene in gene_list:
            if "coagulation" in list_name.lower() and _clinvar_pathogenic(variant.clinvar):
                # v0.7.2: weight scaled by ClinVar review status confidence
                _add_evidence("ClinVar", f"Pathogenic in coagulation gene {gene} → Tier 1 (review_status={variant.clinvar_review_status}, conf={clinvar_conf:.2f})", weight=1.0*clinvar_conf, confidence="high" if clinvar_conf >= 0.8 else "medium" if clinvar_conf >= 0.5 else "low", raw_data={"clinvar": variant.clinvar, "gene_list": list_name, "review_status": variant.clinvar_review_status, "clinvar_conf": clinvar_conf})
                actions.append("Assess bleeding history and coagulation function")
                actions.append("Consider peripheral blood stem cell over bone marrow if applicable")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = []  # Tier 1: no upgrade
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 1, f"{gene} pathogenic variant in coagulation gene - bleeding risk", actions
            if "fa_dna_repair" in list_name.lower() and _clinvar_pathogenic(variant.clinvar):
                actions.append("Assess if patient has Fanconi anemia phenotype")
                actions.append("Biallelic: high personal health risk; heterozygous: carrier status")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = []  # Tier 1: no upgrade
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 1, f"{gene} pathogenic variant in FA pathway - marrow failure risk", actions

    # 1b. Homozygous truncating in primary tissue gene
    if variant.gt in ["1/1", "1|1"] and _impact_high(variant.impact):
        if tissue_assessment.get("relevance") == "primary":
            # === v0.9.1: gnomAD frequency guard (DDX3X hotfix) ===
            if variant.gnomad_status == "API_FAILED":
                _add_evidence("gnomAD", f"gnomAD API FAILED ({variant.gnomad_error_msg}) — cannot confirm rarity. Downgrading from Tier 1 to Tier 2.", weight=-0.8, confidence="low", raw_data={"gnomad_status": variant.gnomad_status, "error": variant.gnomad_error_msg})
                actions.append("⚠️ gnomAD query failed — frequency unverified. Downgraded to Tier 2 pending external verification.")
                variant.gnomad_af_warning = True
                variant.evidence_chain = evidence_chain
                upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                upgrade_conditions.append("若gnomAD查询恢复正常且确认AF<1%,可升级为Tier 1")
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 2, f"Priority 1b candidate (homozygous HIGH in primary tissue gene {gene}), but gnomAD query FAILED ({variant.gnomad_error_msg}). Downgraded to Tier 2 pending frequency verification.", actions
            elif variant.gnomad_status == "NOT_CAPTURED":
                _add_evidence("Zygosity", f"Homozygous LOF in primary tissue gene {gene} → Tier 1 (gnomAD NOT_CAPTURED, confidence=MEDIUM)", weight=1.0, confidence="medium", raw_data={"gt": variant.gt, "impact": variant.impact, "relevance": "primary", "gnomad_status": "NOT_CAPTURED"})
                actions.append("Confirm homozygosity via secondary method")
                actions.append("Assess if phenotype is consistent with expected tissue function")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = []  # Tier 1: no upgrade
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 1, f"Homozygous truncating variant in primary tissue gene {gene} (gnomAD not captured — may be rare/indel)", actions
            # === v0.9.1: known common polymorphism guard ===
            elif variant.gnomad_af is not None and variant.gnomad_af > 0.01:
                _add_evidence("gnomAD", f"AF={variant.gnomad_af:.3f} > 1% — common polymorphism, not Tier 1", weight=-1.0, confidence="high", raw_data={"gnomad_af": variant.gnomad_af})
                actions.append("Common polymorphism — no clinical action needed for homozygous state")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = []  # Tier 3: no upgrade
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 3, f"Homozygous HIGH in {gene} but AF={variant.gnomad_af:.3f} > 1% — common polymorphism", actions
            # === original logic (SUCCESS with confirmed rare AF) ===
            _add_evidence("Zygosity", f"Homozygous LOF in primary tissue gene {gene} → Tier 1", weight=1.0, confidence="high", raw_data={"gt": variant.gt, "impact": variant.impact, "relevance": "primary", "gnomad_af": variant.gnomad_af})
            actions.append("Confirm homozygosity via secondary method")
            actions.append("Assess if phenotype is consistent with expected tissue function")
            variant.evidence_chain = evidence_chain
            upgrade_conditions = []  # Tier 1: no upgrade
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            return 1, f"Homozygous truncating variant in primary tissue gene {gene}", actions

    # Priority 1c: ClinVar Pathogenic + HIGH impact + primary/secondary tissue
    # v0.5.2 FIX: Heterozygous pathogenic truncating variants in tissue-relevant genes
    # were incorrectly falling to Tier 2. ClinVar Pathogenic + HIGH + relevant tissue
    # should be Tier 1 regardless of zygosity (heterozygous pathogenic = actionable).
    # v0.7 Phase 3: If phenotype_match_score is provided and < 0.6 → Tier 2 (phenotype mismatch)
    if _clinvar_pathogenic(variant.clinvar) and _impact_high(variant.impact):
        if tissue_assessment.get("relevance") in ["primary", "secondary"]:
            # v0.7 Phase 3: Check phenotype match score
            pms = getattr(variant, 'phenotype_match_score', None)
            if pms is not None and pms < 0.6:
                # Phenotype mismatch → Tier 2 (not Tier 3)
                _add_evidence("ClinVar", f"Pathogenic + HIGH + tissue-relevant but phenotype mismatch (score={pms:.2f}) → Tier 2 (review_status={variant.clinvar_review_status}, conf={clinvar_conf:.2f})", weight=0.30*clinvar_conf, confidence="medium" if clinvar_conf >= 0.5 else "low", raw_data={"clinvar": variant.clinvar, "impact": variant.impact, "phenotype_match_score": pms, "review_status": variant.clinvar_review_status, "clinvar_conf": clinvar_conf})
                actions.append("Confirm variant via secondary method")
                actions.append("Assess phenotypic severity and clinical relevance")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                upgrade_conditions.append(f"若表型验证与 {gene} 已知疾病匹配,可升级为 Tier 1")
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 2, f"ClinVar pathogenic {variant.consequence} in tissue-relevant gene {gene} - phenotype mismatch (score={pms:.2f})", actions
            # Phenotype match or no phenotype provided → Tier 1 (original logic)
            _add_evidence("ClinVar", f"Pathogenic + HIGH + tissue-relevant → Tier 1 for {gene} (review_status={variant.clinvar_review_status}, conf={clinvar_conf:.2f})", weight=1.0*clinvar_conf, confidence="high" if clinvar_conf >= 0.8 else "medium" if clinvar_conf >= 0.5 else "low", raw_data={"clinvar": variant.clinvar, "impact": variant.impact, "relevance": tissue_assessment.get("relevance"), "review_status": variant.clinvar_review_status, "clinvar_conf": clinvar_conf})
            actions.append("Confirm variant via secondary method")
            actions.append("Assess phenotypic severity and clinical relevance")
            variant.evidence_chain = evidence_chain
            upgrade_conditions = []
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            return 1, f"ClinVar pathogenic {variant.consequence} in tissue-relevant gene {gene}", actions

    # v0.5 P1-4 + P1-5: Gene constraint tier adjustment with NMD refinement
    # Only applies to LOF variants in strongly constrained genes
    # Placed BEFORE Priority 2 so that NMD-sensitive LOF-intolerant variants
    # go directly to Tier 1, while NMD-escape variants fall through to Priority 2/3
    constraint_eval = evaluate_gene_constraint(variant)
    if constraint_eval.get("tier_adjustment") == 1:
        # P1-5: Check NMD prediction before applying PVS1
        nmd = variant.nmd_prediction or predict_nmd(variant)
        nmd_status = nmd.get("status", "unknown")

        if nmd_status == "escape":
            # Last exon - NMD escape, PVS1 does NOT apply
            # Do NOT upgrade to Tier 1, continue to Priority 2/3
            pass  # Fall through to Priority 2
        elif nmd_status == "possible_escape":
            # Penultimate exon - possible escape, PVS1 downgraded to PM/PP
            if variant.gt in ["0/1", "0|1"] and _impact_high(variant.impact):
                if tissue_assessment.get("relevance") in ["primary", "secondary"]:
                    reason = f"Heterozygous {variant.consequence} in LOF-intolerant gene {gene}"
                    reason += f" ({constraint_eval['reason']})"
                    reason += " - PVS1_Strong→PM: possible NMD escape in penultimate exon"
                    actions.append(f"Gene constraint: pLI={constraint_eval['pLI']:.2f}, LOEUF={constraint_eval['loeuf']:.2f}")
                    actions.append("Possible NMD escape - PVS1 downgraded to moderate evidence")
                    actions.append("Haploinsufficiency possible but uncertain - consider functional validation")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 2, reason, actions  # Tier 2, not Tier 1
                else:
                    reason = f"ClinVar pathogenic variant in LOF-intolerant gene {gene}"
                    reason += f" ({constraint_eval['reason']})"
                    reason += " - possible NMD escape, PVS1_Strong→PM"
                    actions.append(f"Gene constraint: pLI={constraint_eval['pLI']:.2f}, LOEUF={constraint_eval['loeuf']:.2f}")
                    actions.append("Possible NMD escape - functional assessment needed")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 2, reason, actions
        elif nmd_status == "unknown":
            # NMD uncertain - conservative: apply PVS1 but annotate uncertainty
            if variant.gt in ["0/1", "0|1"] and _impact_high(variant.impact):
                if tissue_assessment.get("relevance") in ["primary", "secondary"]:
                    reason = f"Heterozygous {variant.consequence} in LOF-intolerant gene {gene}"
                    reason += f" ({constraint_eval['reason']})"
                    reason += f" - NMD status unknown ({nmd.get('reason', 'assuming sensitive')}), PVS1 applied conservatively"
                    actions.append(f"Gene constraint: pLI={constraint_eval['pLI']:.2f}, LOEUF={constraint_eval['loeuf']:.2f}")
                    actions.append("NMD prediction uncertain - assumed sensitive, functional validation recommended")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = []  # Tier 1: no upgrade
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 1, reason, actions
                else:
                    reason = f"ClinVar pathogenic variant in LOF-intolerant gene {gene}"
                    reason += f" ({constraint_eval['reason']})"
                    reason += " - NMD uncertain, assumed sensitive"
                    actions.append(f"Gene constraint: pLI={constraint_eval['pLI']:.2f}, LOEUF={constraint_eval['loeuf']:.2f}")
                    actions.append("Non-tissue-relevant but LOF-intolerant - patient's own health risk")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 2, reason, actions
        else:
            # NMD sensitive - classic PVS1 applies
            if variant.gt in ["0/1", "0|1"] and _impact_high(variant.impact):
                if tissue_assessment.get("relevance") in ["primary", "secondary"]:
                    _add_evidence("GeneConstraint", f"LOF-intolerant + NMD-sensitive → Tier 1 (PVS1) for {gene}", weight=0.8, confidence="high", raw_data={"pLI": constraint_eval.get('pLI'), "loeuf": constraint_eval.get('loeuf'), "nmd_status": "sensitive"})
                    reason = f"Heterozygous {variant.consequence} in LOF-intolerant gene {gene}"
                    reason += f" ({constraint_eval['reason']})"
                    reason += " — NMD sensitive, PVS1 fully applicable"
                    actions.append(f"Gene constraint: pLI={constraint_eval['pLI']:.2f}, LOEUF={constraint_eval['loeuf']:.2f}")
                    actions.append("Haploinsufficiency likely — consider functional validation")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = []  # Tier 1: no upgrade
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)

                    # v0.9.5: SpliceAI evidence for Tier 1 splice variants.
                    # Strong SpliceAI (delta >= 0.5) strengthens PP3 evidence.
                    # Delta=0 suggests VEP HIGH may be overcalled → downgrade to Tier 2.
                    _VALID_SA_SOURCES = {"spliceai", "spliceai_lookup", "vep_rest"}
                    if getattr(config, 'spliceai_enabled', False) and variant.spliceai_result:
                        sa = variant.spliceai_result
                        sa_source = getattr(sa, "source", "")
                        sa_delta = getattr(sa, "delta_score", 0.0) or 0.0
                        if sa_source in _VALID_SA_SOURCES:
                            if sa_delta >= 0.5:
                                _add_evidence("SpliceAI", f"SpliceAI strong (delta={sa_delta:.2f}) for {variant.consequence} — strengthens Tier 1 splice evidence", weight=0.6, confidence="high", raw_data={"delta_score": sa_delta, "details": getattr(sa, 'raw_response', None)})
                                actions.append(f"SpliceAI predicts strong splice disruption (delta={sa_delta:.2f}) — confirm via RNA-seq")
                            elif sa_delta == 0.0:
                                _add_evidence("SpliceAI", f"SpliceAI delta=0 — no splice change predicted for {variant.consequence}, downgrading from Tier 1", weight=-0.5, confidence="high", raw_data={"delta_score": 0.0, "predicted_impact": "none"})
                                actions.append("SpliceAI predicts no splice disruption — VEP HIGH may be overcalled; consider RNA-seq validation")
                                variant.evidence_chain = evidence_chain
                                upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                                variant.upgrade_conditions = upgrade_conditions
                                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                                return 2, f"SpliceAI delta=0 — no splice change for {gene} {variant.consequence}, downgraded from Tier 1", actions

                    variant.evidence_chain = evidence_chain
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 1, reason, actions
                else:
                    reason = f"ClinVar pathogenic variant in LOF-intolerant gene {gene}"
                    reason += f" ({constraint_eval['reason']})"
                    reason += " - NMD sensitive, patient's own health risk"
                    actions.append(f"Gene constraint: pLI={constraint_eval['pLI']:.2f}, LOEUF={constraint_eval['loeuf']:.2f}")
                    actions.append("Non-tissue-relevant but LOF-intolerant - assess phenotypic impact")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 2, reason, actions

    # Priority 2: Tier 2 checks
    # 2a. Primary/secondary tissue gene, heterozygous, function affected
    # v0.10.3: Include "unknown" relevance for HIGH-impact variants.
    # When GTEx data is unavailable (offline mode), tissue relevance defaults to
    # "unknown" rather than "none". HIGH-impact truncating variants should still
    # be flagged for monitoring (Tier 2) rather than silently relegated to Tier 3.
    _relevance = tissue_assessment.get("relevance")
    if _relevance in ["primary", "secondary", "unknown"] and variant.gt in ["0/1", "0|1"]:
        if _impact_high(variant.impact):
            if _relevance == "unknown":
                _add_evidence("TissueRelevance", f"Heterozygous LOF in {gene} → Tier 2 (tissue relevance unknown, GTEx offline)", weight=0.5, confidence="low", raw_data={"relevance": _relevance, "gt": variant.gt, "impact": variant.impact, "domain": domain_info.get("domain") if domain_info else None})
                reason = f"Heterozygous {variant.consequence} in {gene} — tissue relevance unknown (GTEx offline), conservative Tier 2"
            else:
                _add_evidence("TissueRelevance", f"Heterozygous LOF in tissue-relevant {gene} → Tier 2", weight=0.6, confidence=_confidence_from_data(), raw_data={"relevance": _relevance, "gt": variant.gt, "impact": variant.impact, "domain": domain_info.get("domain") if domain_info else None})
                reason = f"Heterozygous {variant.consequence} in tissue-relevant gene {gene}"
            if domain_info and domain_info.get("domain_integrity") in ["completely_destroyed", "partially_destroyed"]:
                reason += f", {domain_info['domain']} domain disrupted"
            actions.append("Inform patient of carrier status")
            actions.append("Monitor post-intervention recovery/function")
            variant.evidence_chain = evidence_chain
            upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)

            # v0.8.0: SpliceAI downgrade check for HIGH-impact variants ascribed to canonical splice in Tier 2
            # v0.9.5: Accept both "spliceai" and "vep_rest" as valid SpliceAI sources.
            _VALID_SA_SOURCES = {"spliceai", "spliceai_lookup", "vep_rest"}
            if getattr(config, 'spliceai_enabled', False) and variant.spliceai_result:
                sa = variant.spliceai_result
                if getattr(sa, "source", "") in _VALID_SA_SOURCES and getattr(sa, "delta_score", 0.0) == 0.0:
                    _add_evidence("SpliceAI", f"SpliceAI delta=0 — no splice change predicted for {variant.consequence}, downgrading from Tier 2", weight=-0.5, confidence="high", raw_data={"delta_score": 0.0, "predicted_impact": "none"})
                    actions.append("SpliceAI predicts no splice disruption — VEP HIGH may be overcalled; consider RNA-seq validation")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = _generate_upgrade_conditions(variant, 3, tissue_assessment, gnomad_info)
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 3, f"SpliceAI delta=0 — no splice change for {gene} {variant.consequence}, downgraded from Tier 2", actions

            return 2, reason, actions

        # v0.5 P1-5: Missense stratification (impact is MODERATE, not HIGH)
        if "missense" in variant.consequence.lower():
            missense_eval = evaluate_missense_tier(variant, domain_info, variant.gene_constraint)
            if missense_eval.get("tier_recommendation") == 2:
                reason = f"Heterozygous missense in tissue-relevant gene {gene}"
                reason += f" - {missense_eval['reason']}"
                actions.append("Inform patient of carrier status")
                actions.append("Monitor post-intervention recovery/function")
                variant.evidence_chain = evidence_chain
                upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                variant.upgrade_conditions = upgrade_conditions
                variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                return 2, reason, actions
            elif missense_eval.get("tier_recommendation") == 3:
                # Missense is tolerated - continue to Priority 3
                pass  # Fall through to Tier 3 logic

    # 2b. Non-primary but ClinVar pathogenic
    if _clinvar_pathogenic(variant.clinvar) and tissue_assessment.get("relevance") == "none":
        # v0.7.2: weight scaled by ClinVar review status confidence
        _add_evidence("ClinVar", f"Pathogenic but non-tissue-relevant {gene} → Tier 2 (review_status={variant.clinvar_review_status}, conf={clinvar_conf:.2f})", weight=0.7*clinvar_conf, confidence="high" if clinvar_conf >= 0.8 else "medium" if clinvar_conf >= 0.5 else "low", raw_data={"clinvar": variant.clinvar, "relevance": "none", "review_status": variant.clinvar_review_status, "clinvar_conf": clinvar_conf})
        actions.append("Inform patient of genetic finding")
        actions.append("Refer for relevant specialist evaluation if indicated")
        variant.evidence_chain = evidence_chain
        upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
        variant.upgrade_conditions = upgrade_conditions
        variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
        return 2, f"ClinVar pathogenic variant in {gene} - patient's own health may be affected", actions

    # 2c. Drug metabolism genes (if applicable to this tissue context)
    drug_genes = special_lists.get("drug_metabolism", [])
    if gene in drug_genes:
        actions.append(f"Monitor post-intervention drug levels if relevant medications used")
        variant.evidence_chain = evidence_chain
        upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
        variant.upgrade_conditions = upgrade_conditions
        variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
        return 2, f"Drug metabolism variant may affect pharmacokinetics", actions

    # Priority 2d: Truncating variants with very low/absent population frequency
    # When all annotation layers (tissue mapping, gene constraint, ClinVar) are
    # incomplete, the combination of HIGH impact + rare frequency is itself
    # sufficient biological evidence for Tier 2 (questionable). This is a
    # catch-all safety net — NOT a gene-specific override.
    if _impact_high(variant.impact):
        af_very_rare = gnomad_af is None or (gnomad_af is not None and gnomad_af < 0.001)  # < 0.1%
        clinvar_not_benign = not _clinvar_benign(variant.clinvar)
        if af_very_rare and clinvar_not_benign:
            _add_evidence(
                "RareTruncating",
                f"Truncating variant with AF={'None' if gnomad_af is None else f'{gnomad_af:.5f}'} "
                f"— Tier 2 (questionable) due to biological severity despite limited annotation",
                weight=0.6,
                confidence="low",
                raw_data={"gnomad_af": gnomad_af, "impact": variant.impact, "consequence": variant.consequence},
            )
            actions.append(
                f"Rare truncating variant in {gene} — Tier 2 (questionable); "
                f"recommend Sanger validation and phenotype correlation"
            )
            variant.evidence_chain = evidence_chain
            upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            return 2, f"Rare truncating variant in {gene} — Tier 2 (questionable)", actions

    # Priority 3: Tier 3 - everything else
    reason_parts = []
    if gnomad_info.get("status") == "common_polymorphism":
        # v0.7 Phase 3: Rare disease genes with AF>1% - do NOT auto Tier 3
        if _is_rare_disease_gene(gene) and not _clinvar_benign(variant.clinvar):
            _add_evidence("Frequency", f"AF>1% but rare disease gene {gene} → Tier 2 (not auto Tier 3)", weight=0.2, confidence="medium", raw_data={"gnomad_status": gnomad_info.get("status"), "af": gnomad_info.get("af"), "gene": gene, "rare_disease": True})
            variant.qc_flags.append("COMMON_POLYMORPHISM_BUT_RARE_DISEASE_GENE")
            actions.append(f"Rare disease gene {gene} with common polymorphism - monitor for phenotype correlation")
            variant.evidence_chain = evidence_chain
            upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
            variant.upgrade_conditions = upgrade_conditions
            variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
            return 2, f"Common polymorphism in rare disease gene {gene}", actions
        reason_parts.append(f"Common polymorphism (AF={gnomad_info.get('af')})")
    if _clinvar_benign(variant.clinvar):
        reason_parts.append("ClinVar benign")
        # v0.7.2: Add separate ClinVar benign evidence with negative weight scaled by review status
        _add_evidence("ClinVar", f"ClinVar Benign (review_status={variant.clinvar_review_status}, conf={clinvar_conf:.2f})", weight=-0.5*clinvar_conf, confidence="high" if clinvar_conf >= 0.8 else "medium" if clinvar_conf >= 0.5 else "low", raw_data={"clinvar": variant.clinvar, "review_status": variant.clinvar_review_status, "clinvar_conf": clinvar_conf})
    if tissue_assessment.get("relevance") == "none":
        reason_parts.append("No tissue relevance")

    reason = "; ".join(reason_parts) if reason_parts else "Low risk based on combined assessment"
    _add_evidence("Frequency", f"Common polymorphism / benign / no tissue relevance → Tier 3", weight=0.2, confidence="high", raw_data={"gnomad_status": gnomad_info.get("status"), "clinvar": variant.clinvar, "relevance": tissue_assessment.get("relevance")})

    # v0.8.0: SpliceAI evidence for Tier 3 splice variants (default OFF)
    # v0.9.5: Accept both "spliceai" and "vep_rest" as valid SpliceAI sources.
    # If SpliceAI is enabled and pre-computed result exists, evaluate upgrade/downgrade.
    _VALID_SA_SOURCES = {"spliceai", "spliceai_lookup", "vep_rest"}
    if getattr(config, 'spliceai_enabled', False) and variant.spliceai_result:
        sa = variant.spliceai_result
        source = getattr(sa, 'source', 'unknown')
        if source in _VALID_SA_SOURCES:
            delta = getattr(sa, 'delta_score', None)
            if delta is not None:
                impact = getattr(sa, 'predicted_impact', 'none')
                if impact == "strong" and delta >= 0.5:
                    # Strong splice prediction for a Tier 3 variant → upgrade to Tier 2
                    _add_evidence("SpliceAI", f"SpliceAI strong (delta={delta:.2f}) for {variant.consequence} -> upgrade to Tier 2", weight=0.8, confidence="high", raw_data={"delta_score": delta, "predicted_impact": impact, "details": getattr(sa, 'raw_response', None)})
                    actions.append("SpliceAI predicts strong splice disruption - confirm via RNA-seq or functional assay")
                    variant.evidence_chain = evidence_chain
                    upgrade_conditions = _generate_upgrade_conditions(variant, 2, tissue_assessment, gnomad_info)
                    variant.upgrade_conditions = upgrade_conditions
                    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
                    return 2, f"SpliceAI strong splice prediction (delta={delta:.2f}) for {gene} {variant.consequence} - upgraded from Tier 3", actions
                elif impact == "moderate" and delta >= 0.2:
                    _add_evidence("SpliceAI", f"SpliceAI moderate (delta={delta:.2f}) for {variant.consequence}", weight=0.4, confidence="medium", raw_data={"delta_score": delta, "predicted_impact": impact, "details": getattr(sa, 'raw_response', None)})
                elif impact in ("weak", "none") and delta == 0.0:
                    # No predicted splice change - supports VEP overcall, keep Tier 3
                    _add_evidence("SpliceAI", f"SpliceAI delta=0 - no splice change predicted for {variant.consequence}", weight=-0.5, confidence="high", raw_data={"delta_score": delta, "predicted_impact": impact, "details": getattr(sa, 'raw_response', None)})
        elif source == "api_error":
            variant.qc_flags.append("SPLICEAI_API_ERROR")
        elif source == "not_in_db":
            _add_evidence("SpliceAI", "Not in SpliceAI database - no splice prediction available", weight=0.0, confidence="low", raw_data={"source": "not_in_db"})

    # v0.5 P1-11: Generate upgrade conditions before final tier assignment
    upgrade_conditions = _generate_upgrade_conditions(variant, 3, tissue_assessment, gnomad_info)
    variant.upgrade_conditions = upgrade_conditions
    variant.evidence_chain = evidence_chain
    variant.tier_confidence = _calculate_tier_confidence(evidence_chain)
    return 3, reason, []
