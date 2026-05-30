#!/usr/bin/env python3
"""
GPA Report Generation Module

Markdown and JSON report builders for GPA analysis results.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
"""

from __future__ import annotations

import json
import re
from datetime import datetime
from pathlib import Path
from typing import List, Dict, Optional, Any, TYPE_CHECKING

try:
    from jinja2 import Template
    _HAS_JINJA2 = True
except ImportError:
    _HAS_JINJA2 = False
    Template = None  # type: ignore

if TYPE_CHECKING:
    from dgra_core import Variant, GPAConfig

from version import __version__

# Offline archive path (shared with dgra_core)
_OFFLINE_ARCHIVE_DIR = Path(__file__).resolve().parent.parent / "references" / "offline_data"

# v0.10.0: Jinja2 report templates — header/summary separated from per-variant logic
if _HAS_JINJA2:
    _REPORT_HEADER_TEMPLATE = Template("""
# GPA Report - Genomic Phenotype Association v{{ version }}

**Analysis Context**: {{ profile_name }}
**Tissue Profile**: `{{ tissue_profile }}`
**Offline Mode**: {{ 'Yes' if offline_mode else 'No' }}
**Analysis Date**: {{ analysis_date }}
{% if filter_stats %}
**Input Variants**: {{ filter_stats.input_count | default('?') }} → **Assessed**: {{ filter_stats.output_count | default('?') }} (excluded: {{ filter_stats.excluded | default('?') }})
{% if impact_distribution %}**Impact Distribution (input)**: {{ impact_distribution }}
{% endif %}
{% if filter_retention %}**Filter Retention**: {{ filter_retention }}
{% endif %}
{% if filter_preset %}**Filter Preset**: `{{ filter_preset }}`
{% endif %}
{% else %}
**Total Variants Assessed**: {{ total_variants }}
{% endif %}
{% if version_info %}
**GPA Version**: {{ version_info.dgra_version | default('0.9.0') }}
{% if version_info.cache_version %}**Cache Version**: {{ version_info.cache_version }}
{% endif %}
{% if version_info.offline_archive_date and version_info.offline_archive_date not in ('no_archive', 'empty', 'unknown') %}**Offline Archive Date**: {{ version_info.offline_archive_date }}
{% endif %}
{% if version_info.dgra_core_commit and version_info.dgra_core_commit != 'unknown' %}**Code Commit**: `{{ version_info.dgra_core_commit }}`
{% endif %}
{% if version_info.database_version %}**Database Version**: {{ version_info.database_version }}
{% endif %}
{% endif %}

**Tier 1 基因**: {{ tier1_genes }} 个 | **Tier 1 突变**: {{ tier1_variants }} 个
**Tier 2 基因**: {{ tier2_genes }} 个 | **Tier 2 突变**: {{ tier2_variants }} 个
**Tier 3 基因**: {{ tier3_genes }} 个 | **Tier 3 突变**: {{ tier3_variants }} 个
""")
else:
    _REPORT_HEADER_TEMPLATE = None  # type: ignore


def _build_report_header(
    profile_name: str,
    config: GPAConfig,
    variants: List[Variant],
    tier1: List[Variant],
    tier2: List[Variant],
    tier3: List[Variant],
) -> str:
    """Fallback header builder when Jinja2 is not available."""
    report_lines = []
    report_lines.append(f"# GPA Report - Genomic Phenotype Association v{__version__}\n")
    report_lines.append(f"**Analysis Context**: {profile_name}\n")
    report_lines.append(f"**Tissue Profile**: `{config.tissue_profile}`\n")
    report_lines.append(f"**Offline Mode**: {'Yes' if config.offline_mode else 'No'}\n")
    report_lines.append(f"**Analysis Date**: {datetime.now().isoformat()}\n")

    if config.filter_stats:
        fs = config.filter_stats
        report_lines.append(
            f"**Input Variants**: {fs.get('input_count', '?')} → **Assessed**: "
            f"{fs.get('output_count', '?')} (excluded: {fs.get('excluded', '?')})\n"
        )
        impact_parts = []
        for imp in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']:
            cnt = fs.get('by_impact', {}).get(imp, 0)
            if cnt > 0:
                impact_parts.append(f"{imp}: {cnt}")
        if impact_parts:
            report_lines.append(f"**Impact Distribution (input)**: {' | '.join(impact_parts)}\n")
        retention_parts = []
        if fs.get('splice_retained', 0) > 0:
            retention_parts.append(f"splice retained: {fs['splice_retained']}")
        if fs.get('synonymous_tissue_retained', 0) > 0:
            retention_parts.append(f"synonymous tissue: {fs['synonymous_tissue_retained']}")
        if fs.get('clinvar_conflicting_retained', 0) > 0:
            retention_parts.append(f"ClinVar conflict: {fs['clinvar_conflicting_retained']}")
        if retention_parts:
            report_lines.append(f"**Filter Retention**: {' | '.join(retention_parts)}\n")
        if config.filter_preset:
            report_lines.append(f"**Filter Preset**: `{config.filter_preset}`\n")
    else:
        report_lines.append(f"**Total Variants Assessed**: {len(variants)}\n")

    version_info = _get_version_info(config)
    report_lines.append(f"**GPA Version**: {version_info.get('dgra_version', '0.9.0')}\n")
    if version_info.get('cache_version'):
        report_lines.append(f"**Cache Version**: {version_info['cache_version']}\n")
    if version_info.get('offline_archive_date') and version_info['offline_archive_date'] not in ('no_archive', 'empty', 'unknown'):
        report_lines.append(f"**Offline Archive Date**: {version_info['offline_archive_date']}\n")
    if version_info.get('dgra_core_commit') and version_info['dgra_core_commit'] != 'unknown':
        report_lines.append(f"**Code Commit**: `{version_info['dgra_core_commit']}`\n")
    if version_info.get('database_version'):
        report_lines.append(f"**Database Version**: {version_info['database_version']}\n")
    report_lines.append("\n")

    report_lines.append(f"**Tier 1 基因**: {len(set(v.gene for v in tier1))} 个 | **Tier 1 突变**: {len(tier1)} 个\n")
    report_lines.append(f"**Tier 2 基因**: {len(set(v.gene for v in tier2))} 个 | **Tier 2 突变**: {len(tier2)} 个\n")
    report_lines.append(f"**Tier 3 基因**: {len(set(v.gene for v in tier3))} 个 | **Tier 3 突变**: {len(tier3)} 个\n\n")

    return "".join(report_lines)


def _get_version_info(config: GPAConfig) -> Dict:
    """
    Gather analysis version and provenance metadata.
    v0.5 P1-15: Full version tracking for reproducibility.
    """
    import hashlib
    import subprocess

    version_info = {
        "dgra_version": "0.9.0",
        "analysis_date": datetime.now().isoformat(),
    }

    # Cache version: hash of SQLite cache file
    cache_path = getattr(config, 'cache_db_path', None)
    if cache_path and Path(cache_path).exists():
        try:
            with open(cache_path, 'rb') as f:
                cache_hash = hashlib.sha256(f.read()).hexdigest()[:16]
            version_info["cache_version"] = cache_hash
        except Exception:
            version_info["cache_version"] = "unknown"
    else:
        version_info["cache_version"] = "no_cache"

    # Offline archive: latest file modification time
    if _OFFLINE_ARCHIVE_DIR.exists():
        try:
            mtimes = [p.stat().st_mtime for p in _OFFLINE_ARCHIVE_DIR.iterdir() if p.is_file()]
            if mtimes:
                latest_mtime = max(mtimes)
                version_info["offline_archive_date"] = datetime.fromtimestamp(latest_mtime).isoformat()
            else:
                version_info["offline_archive_date"] = "empty"
        except Exception:
            version_info["offline_archive_date"] = "unknown"
    else:
        version_info["offline_archive_date"] = "no_archive"

    # Git commit hash of dgra_core.py
    try:
        script_dir = Path(__file__).resolve().parent
        result = subprocess.run(
            ["git", "-C", str(script_dir), "rev-parse", "HEAD"],
            capture_output=True, text=True, timeout=5
        )
        if result.returncode == 0:
            version_info["dgra_core_commit"] = result.stdout.strip()[:12]
        else:
            version_info["dgra_core_commit"] = "unknown"
    except Exception:
        version_info["dgra_core_commit"] = "unknown"

    # Database version override (CLI --database-version)
    db_version = getattr(config, 'database_version', None)
    if db_version:
        version_info["database_version"] = db_version

    return version_info

# =============================================================================
# =============================================================================
# Report Generation
# =============================================================================

def _format_vep_reannotation_note(v: Variant) -> Optional[str]:
    """Format a human-readable VEP reannotation note for report display.

    Returns None if the variant was not VEP-reannotated.
    """
    if not v.transcript_warning:
        return None
    try:
        tw = json.loads(v.transcript_warning)
    except (json.JSONDecodeError, TypeError):
        return None

    vep = tw.get("vep_reannotation")
    if not vep or vep.get("status") != "success":
        return None

    original = vep.get("original", {})
    canonical = vep.get("canonical", {})
    orig_tx = original.get("transcript", "N/A")
    orig_consq = original.get("consequence", "N/A")
    orig_impact = original.get("impact", "N/A")
    new_tx = canonical.get("transcript_id", v.transcript or "N/A")
    new_consq = v.consequence or canonical.get("consequence", "N/A")
    new_impact = v.impact or canonical.get("impact", "N/A")

    note = (
        f"⚠️ 后果已按 canonical transcript ({new_tx}) 重新注释:"
        f"原 {orig_tx} {orig_consq}/{orig_impact} → "
        f"VEP 结果 {new_consq}/{new_impact}"
    )
    return note


def _generate_pseudogene_assessment_section(variants: List[Variant]) -> Optional[str]:
    """
    v0.6: Generate standalone pseudogene interference assessment section.

    Returns Markdown string with:
      - Summary table of all variants with pseudogene warnings
      - Score-based classification (0.0-1.0)
      - Per-gene pseudogene list from lookup
      - Recommendations (Sanger, long-read, etc.)

    Does NOT modify tier - only provides analytical assessment.
    """
    # Collect variants with pseudogene warnings
    pg_variants = []
    for v in variants:
        if v.pseudogene_warning:
            try:
                pw = json.loads(v.pseudogene_warning)
                if pw.get("score", 0) > 0:
                    pg_variants.append((v, pw))
            except (json.JSONDecodeError, TypeError):
                continue

    if not pg_variants:
        return None

    lines = []
    lines.append("## 🧬 假基因干扰评估\n")
    lines.append(f"*基于 v0.6 轻量版假基因数据库({len(_load_pseudogene_lookup())} 个临床相关基因对)*\n")
    lines.append(f"**检测到 {len(pg_variants)} 个变异存在假基因干扰风险**\n\n")

    # Summary table
    lines.append("| 基因 | 位点 | 基因型 | 观察VAF | 预期VAF | 干扰评分 | 等级 | 假基因 | 建议 |\n")
    lines.append("|------|------|--------|---------|---------|----------|------|--------|------|\n")

    # Severity classification
    severity_map = {
        "strong_interference": ("🔴", "高度干扰"),
        "interference": ("🟠", "中度干扰"),
        "suspected": ("🟡", "疑似干扰"),
        "bias_suspected": ("🟡", "疑似偏倚"),
    }

    for v, pw in pg_variants:
        level = pw.get("level", "unknown")
        icon, desc = severity_map.get(level, ("⚪", "未知"))

        pos = f"{v.chrom}:{v.pos}"
        gt = v.gt or "N/A"
        obs_vaf = pw.get("observed_vaf", v.vaf or "N/A")
        exp_vaf = pw.get("expected_vaf", 0.5)
        score = pw.get("score", 0)
        pgs = ", ".join(pw.get("pseudogenes", []))
        rec = pw.get("recommendation", "建议验证")

        # Truncate recommendation for table
        rec_short = rec[:40] + "..." if len(rec) > 40 else rec

        lines.append(f"| {v.gene} | {pos} | {gt} | {obs_vaf} | {exp_vaf} | {score:.2f} | {icon} {desc} | {pgs} | {rec_short} |\n")

    lines.append("\n")

    # Detailed per-variant analysis
    lines.append("### 详细分析\n\n")
    for i, (v, pw) in enumerate(pg_variants, 1):
        level = pw.get("level", "unknown")
        icon, desc = severity_map.get(level, ("⚪", "未知"))

        lines.append(f"**{i}. {v.gene} ({v.chrom}:{v.pos})** {icon} {desc}\n\n")
        lines.append(f"- **变异**: {v.hgvsp or v.hgvsc or 'N/A'}\n")
        lines.append(f"- **基因型**: {v.gt or 'N/A'} | **观察VAF**: {pw.get('observed_vaf', v.vaf)} | **预期VAF**: {pw.get('expected_vaf', 0.5)}\n")
        lines.append(f"- **干扰评分**: {pw.get('score', 0):.2f}/1.0\n")
        lines.append(f"- **检测策略**: {pw.get('strategy', 'vaf_mismatch')}\n")

        pgs = pw.get("pseudogenes", [])
        if pgs:
            lines.append(f"- **相关假基因**: {', '.join(pgs)}\n")
            # Add notes from lookup if available
            lookup = _load_pseudogene_lookup()
            entry = lookup.get(v.gene, {})
            if entry.get("notes"):
                lines.append(f"- **注释**: {entry['notes']}\n")

        lines.append(f"- **建议**: {pw.get('recommendation', '建议验证')}\n")

        # Tier impact note
        lines.append(f"- **分级影响**: 本评估**不直接修改 Tier**,仅下调置信度。当前 Tier {v.tier or 'N/A'},置信度 {v.tier_confidence or 'N/A'}\n")

        lines.append("\n")

    # Overall recommendation
    strong_count = sum(1 for _, pw in pg_variants if pw.get("score", 0) >= 0.75)
    if strong_count > 0:
        lines.append(f"### ⚠️ 重点关注\n\n")
        lines.append(f"检测到 **{strong_count} 个**高度假基因干扰变异(评分≥0.75)。")
        lines.append("强烈建议使用 Sanger 测序或长读长测序验证以下位点:\n\n")
        for v, pw in pg_variants:
            if pw.get("score", 0) >= 0.75:
                lines.append(f"- {v.gene}: {v.chrom}:{v.pos} ({v.hgvsp or v.hgvsc})\n")
        lines.append("\n")

    return "".join(lines)


# =============================================================================
# Phenotype Association Assessment Section (v0.7 Phase 4)
# =============================================================================

def _generate_transcript_selection_section(variants: List[Variant]) -> Optional[str]:
    """
    v0.9.0: Generate transcript selection assessment section.

    Only appears when there are variants with transcript_ambiguity_flag=True
    or variants with alternative_transcripts.

    Lists primary transcript, selection method, and alternatives.
    Flags ambiguous cases with ⚠️.
    """
    # Collect variants with transcript selection info
    tx_variants = []
    for v in variants:
        if v.transcript_ambiguity_flag or v.alternative_transcripts:
            tx_variants.append(v)

    if not tx_variants:
        return None

    lines = []
    lines.append("## 🧬 转录本选择评估\n")
    lines.append(f"*检测到 {len(tx_variants)} 个变异存在多转录本注释，展示主转录本选择结果*\n\n")

    # Summary table
    lines.append("### 汇总表\n\n")
    lines.append("| 基因 | 位点 | 主转录本 | 选择方法 | 后果 | 影响 | 歧义标记 | 备选转录本数 |\n")
    lines.append("|------|------|----------|----------|------|------|----------|--------------|\n")

    for v in tx_variants:
        pos = f"{v.chrom}:{v.pos}"
        primary_tx = v.primary_transcript or v.transcript or "N/A"
        method = v.transcript_selection_method or "canonical"
        primary_consequence = v.primary_consequence or v.consequence or "N/A"
        primary_impact = v.primary_impact or v.impact or "N/A"

        # Ambiguity flag
        if v.transcript_ambiguity_flag:
            ambiguity = "⚠️ 歧义"
        else:
            ambiguity = "✅ 明确"

        alt_count = len(v.alternative_transcripts) if v.alternative_transcripts else 0

        lines.append(
            f"| {v.gene} | {pos} | {primary_tx} | {method} | {primary_consequence} | {primary_impact} | {ambiguity} | {alt_count} |\n"
        )

    lines.append("\n")

    # Detailed per-variant analysis
    lines.append("### 逐变异详细分析\n\n")
    for i, v in enumerate(tx_variants, 1):
        pos = f"{v.chrom}:{v.pos}"
        primary_tx = v.primary_transcript or v.transcript or "N/A"
        method = v.transcript_selection_method or "canonical"

        if v.transcript_ambiguity_flag:
            lines.append(f"**{i}. {v.gene} ({pos})** ⚠️ **转录本选择存在歧义**\n\n")
        else:
            lines.append(f"**{i}. {v.gene} ({pos})**\n\n")

        lines.append(f"- **主转录本**: {primary_tx}\n")
        lines.append(f"- **选择方法**: {method}\n")

        if v.primary_consequence:
            lines.append(f"- **主转录本后果**: {v.primary_consequence}\n")
        if v.primary_hgvsc:
            lines.append(f"- **cDNA 变化 (HGVS.c)**: {v.primary_hgvsc}\n")
        if v.primary_hgvsp:
            lines.append(f"- **蛋白变化 (HGVS.p)**: {v.primary_hgvsp}\n")
        if v.primary_impact:
            lines.append(f"- **主转录本影响等级**: {v.primary_impact}\n")

        if v.transcript_selection_log:
            lines.append(f"- **选择日志**: {v.transcript_selection_log}\n")

        # List alternative transcripts
        if v.alternative_transcripts:
            lines.append(f"- **备选转录本** ({len(v.alternative_transcripts)} 个):\n")
            for alt in v.alternative_transcripts[:5]:
                tx_id = alt.get("transcript_id", "N/A")
                tx_consequence = ", ".join(alt.get("consequence_terms", [])) if isinstance(alt.get("consequence_terms"), list) else alt.get("consequence_terms", "N/A")
                tx_impact = alt.get("impact", "N/A")
                is_canonical = "canonical" if alt.get("canonical") else ""
                is_mane = "MANE" if alt.get("mane_select") else ""
                flags = ", ".join(filter(None, [is_canonical, is_mane]))
                flag_str = f" ({flags})" if flags else ""
                lines.append(f"  - {tx_id}{flag_str}: {tx_consequence} | {tx_impact}\n")
            if len(v.alternative_transcripts) > 5:
                lines.append(f"  - ... 共 {len(v.alternative_transcripts)} 个备选转录本\n")

        # Flag ambiguous cases
        if v.transcript_ambiguity_flag:
            lines.append(
                f"- ⚠️ **歧义警告**: 该变异的主转录本选择存在不确定性，"
                f"top candidates 分数差距小于阈值。"
            )
            if v.transcript_selection_method == "llm_disease_match":
                lines.append("已通过 LLM 疾病描述辅助选择，但仍建议人工复核。\n")
            else:
                lines.append("建议结合临床表现和文献进一步验证。\n")

        lines.append(f"- **当前分级**: Tier {v.tier or 'N/A'} | 置信度: {v.tier_confidence or 'N/A'}\n")
        lines.append("\n")

    return "".join(lines)


# =============================================================================
# Tier 1/2 Mandatory Transcript Validation (v0.10.14)
# =============================================================================

def _generate_tier12_transcript_validation_section(variants: List[Variant]) -> Optional[str]:
    """
    v0.10.14: Mandatory transcript validation for ALL Tier 1/2 variants.

    Compares the selected (primary) transcript against canonical and MANE Select
    transcripts. Flags discrepancies with warnings. Includes isoform count,
    canonical HGVS, protein domain impact, tissue expression, and literature
    concordance assessment.

    This section ALWAYS appears for Tier 1/2 variants regardless of whether
    transcript_ambiguity_flag is set.
    """
    tier12 = [v for v in variants if v.tier in (1, 2)]
    if not tier12:
        return None

    lines = []
    lines.append("## 🧬 转录本验证 (Tier 1/2 强制)\n")
    lines.append(f"*对全部 **{len(tier12)}** 个 Tier 1/2 变异进行转录本校正验证*\n\n")

    # Summary table
    lines.append("### 验证汇总表\n\n")
    lines.append("| 基因 | 位点 | 选中转录本 | 选择方法 | Canonical | MANE Select | Isoform数 | 验证状态 |\n")
    lines.append("|------|------|-----------|----------|-----------|-------------|-----------|----------|\n")

    discrepancy_count = 0
    for v in tier12:
        # Find canonical and MANE Select from alternative_transcripts
        canonical_tx: Optional[str] = None
        mane_tx: Optional[str] = None
        isoform_count = len(v.alternative_transcripts)

        for tx in v.alternative_transcripts:
            if tx.get("canonical"):
                canonical_tx = tx.get("transcript_id", "")
            if tx.get("mane_select"):
                mane_tx = tx.get("transcript_id", "")

        selected_tx = v.primary_transcript or v.transcript or "N/A"

        # Determine validation status
        is_canonical = (selected_tx == canonical_tx) if canonical_tx else False
        is_mane = (selected_tx == mane_tx) if mane_tx else False

        if is_canonical or is_mane:
            status = "一致"
            status_icon = "✅"
        else:
            status = "差异"
            status_icon = "⚠️"
            discrepancy_count += 1

        tx_type = v.transcript_selection_method or "canonical"

        lines.append(
            f"| {v.gene} | {v.chrom}:{v.pos} | {selected_tx} | {tx_type} | "
            f"{canonical_tx or 'N/A'} | {mane_tx or 'N/A'} | {isoform_count} | "
            f"{status_icon} {status} |\n"
        )

    lines.append("\n")

    if discrepancy_count > 0:
        lines.append(
            f"> ⚠️ **警告**: {discrepancy_count}/{len(tier12)} 个 Tier 1/2 变异的选中转录本 "
            f"与 Canonical/MANE Select 不一致，详见下方逐变异分析。\n\n"
        )

    # Detailed per-variant validation
    lines.append("### 逐变异详细验证\n\n")
    for i, v in enumerate(tier12, 1):
        pos = f"{v.chrom}:{v.pos}"
        selected_tx = v.primary_transcript or v.transcript or "N/A"

        # Re-find canonical/MANE for this variant
        canonical_tx = None
        mane_tx = None
        canonical_hgvsc = ""
        canonical_hgvsp = ""
        canonical_domains = []
        selected_domains = []

        for tx in v.alternative_transcripts:
            tx_id = tx.get("transcript_id", "")
            if tx.get("canonical"):
                canonical_tx = tx_id
                canonical_hgvsc = tx.get("hgvsc", "")
                canonical_hgvsp = tx.get("hgvsp", "")
                canonical_domains = tx.get("protein_domains", [])
            if tx.get("mane_select"):
                mane_tx = tx_id
            if tx_id == selected_tx:
                selected_domains = tx.get("protein_domains", [])

        is_canonical = (selected_tx == canonical_tx) if canonical_tx else False
        is_mane = (selected_tx == mane_tx) if mane_tx else False

        if is_canonical or is_mane:
            lines.append(f"**{i}. {v.gene} ({pos})** ✅ 转录本选择验证通过\n\n")
        else:
            lines.append(f"**{i}. {v.gene} ({pos})** ⚠️ **转录本选择存在差异**\n\n")

        lines.append(f"- **选中转录本**: `{selected_tx}` (方法: {v.transcript_selection_method or 'canonical'})\n")
        lines.append(f"- **Canonical 转录本**: `{canonical_tx or 'N/A'}`\n")
        lines.append(f"- **MANE Select 转录本**: `{mane_tx or 'N/A'}`\n")
        lines.append(f"- **Isoform 总数**: {len(v.alternative_transcripts)}\n")

        # HGVS comparison
        if canonical_hgvsc or canonical_hgvsp:
            lines.append(f"- **Canonical HGVS**: c.{canonical_hgvsc or 'N/A'} | p.{canonical_hgvsp or 'N/A'}\n")
        if v.primary_hgvsc or v.primary_hgvsp:
            lines.append(f"- **选中转录本 HGVS**: c.{v.primary_hgvsc or 'N/A'} | p.{v.primary_hgvsp or 'N/A'}\n")

        # Protein domain impact
        if selected_domains or canonical_domains:
            lines.append(f"- **蛋白功能域影响**:\n")
            if selected_domains:
                lines.append(f"  - 选中转录本: {', '.join(str(d) for d in selected_domains[:5])}{'...' if len(selected_domains) > 5 else ''}\n")
            if canonical_domains and canonical_domains != selected_domains:
                lines.append(f"  - Canonical 转录本: {', '.join(str(d) for d in canonical_domains[:5])}{'...' if len(canonical_domains) > 5 else ''}\n")
        else:
            lines.append(f"- **蛋白功能域影响**: 无明确功能域注释\n")

        # Tissue expression (from tissue_relevance)
        if v.tissue_relevance:
            tr = v.tissue_relevance
            tpm = tr.get("gtex_tpm")
            relevance = tr.get("relevance", "N/A")
            if tpm is not None:
                lines.append(f"- **组织表达**: {relevance} (GTEx TPM={tpm:.2f})\n")
            else:
                lines.append(f"- **组织表达**: {relevance} (GTEx 数据不可用)\n")
        else:
            lines.append(f"- **组织表达**: 未评估\n")

        # Discrepancy warning details
        if not is_canonical and canonical_tx:
            lines.append(
                f"- ⚠️ **差异说明**: 选中转录本 `{selected_tx}` 与 Canonical 转录本 "
                f"`{canonical_tx}` 不同。若 Canonical 转录本后果更轻（如同义/内含子），"
                f"当前 Tier 分级可能需下调。\n"
            )
        elif not is_mane and mane_tx and is_canonical:
            lines.append(
                f"- ℹ️ **说明**: 选中转录本为 Canonical 但非 MANE Select。"
                f"MANE Select `{mane_tx}` 是临床首选参考转录本，建议复核。\n"
            )

        # Transcript selection log
        if v.transcript_selection_log:
            lines.append(f"- **选择日志**: {v.transcript_selection_log}\n")

        # Literature concordance
        clinvar_lit = ""
        if v.clinvar and v.clinvar not in ("UNKNOWN", ""):
            clinvar_lit = f"ClinVar: {v.clinvar}"
        if clinvar_lit:
            lines.append(f"- **文献/数据库一致性**: {clinvar_lit}\n")
        else:
            lines.append(f"- **文献/数据库一致性**: 无明确文献记录\n")

        lines.append(f"- **当前分级**: Tier {v.tier} | 置信度: {v.tier_confidence or 'N/A'}\n")
        lines.append("\n")

    # Overall recommendation
    if discrepancy_count > 0:
        lines.append("### ⚠️ 转录本差异处理建议\n\n")
        lines.append(
            f"检测到 **{discrepancy_count}** 个 Tier 1/2 变异的选中转录本与 Canonical/MANE Select 不一致。\n\n"
        )
        lines.append("建议采取以下措施:\n")
        lines.append("1. **人工复核**: 确认选中转录本是否为患者组织中最主要的表达 isoform\n")
        lines.append("2. **功能验证**: 若 Canonical 转录本后果更轻，考虑使用 Canonical 后果重新分级\n")
        lines.append("3. **文献检索**: 查询该基因疾病相关文献，确认致病变异通常发生在哪个转录本上\n")
        lines.append("4. **RNA-seq 验证**: 若条件允许，通过 RNA-seq 确认患者组织中实际表达的主要 isoform\n")
        lines.append("\n")

    return "".join(lines)


# =============================================================================
# Phenotype Association Assessment Section (v0.7 Phase 4)
# =============================================================================

def _generate_phenotype_assessment_section(variants: List[Variant]) -> Optional[str]:
    """
    v0.7 Phase 4: Generate standalone phenotype association assessment section.

    Core principle: Only Tier 1/2 variants undergo phenotype association analysis.
    This saves ~95% computation (typically 5-20 variants out of 500).

    Returns Markdown string with:
      - Summary table of all Tier 1/2 variants with phenotype match results
      - Per-variant detailed analysis (score, explanation, matched pairs, known phenotypes)
      - High-score (≥0.75) variants → validation recommendation
      - ClinVar Pathogenic + low match score → explicit mismatch warning
    """
    # Collect Tier 1/2 variants with phenotype data
    pheno_variants = []
    for v in variants:
        if v.tier in (1, 2) and v.phenotype_match_score is not None:
            pheno_variants.append(v)

    if not pheno_variants:
        return None

    lines = []
    lines.append("## 🧬 表型关联评估\n")
    lines.append(f"*仅对 Tier 1/2 变异执行表型关联分析(共 {len(pheno_variants)} 个变异)*\n\n")

    # Summary table
    lines.append("### 汇总表\n\n")
    lines.append("| 基因 | 位点 | 合子型 | VAF | 匹配评分 | 关联等级 | 假基因状态 | 建议 |\n")
    lines.append("|------|------|--------|-----|----------|----------|------------|------|\n")

    for v in pheno_variants:
        pos = f"{v.chrom}:{v.pos}"
        gt = v.gt or "N/A"
        vaf = f"{v.vaf:.3f}" if v.vaf is not None else "N/A"
        score = v.phenotype_match_score

        # Association level
        conf = v.phenotype_match_confidence or "low"
        if score >= 0.75:
            level_icon = "🟢"
            level_text = "高度关联"
        elif score >= 0.40:
            level_icon = "🟡"
            level_text = "中度关联"
        elif score > 0:
            level_icon = "🔴"
            level_text = "低度关联"
        else:
            level_icon = "⚪"
            level_text = "无关联"

        # Pseudogene status
        pg_status = "无"
        if v.pseudogene_warning:
            try:
                pw = json.loads(v.pseudogene_warning)
                if pw.get("score", 0) > 0:
                    pg_level = pw.get("level", "unknown")
                    pg_status = f"{pw.get('score', 0):.2f} ({pg_level})"
            except (json.JSONDecodeError, TypeError):
                pass

        # Recommendation
        if score >= 0.75:
            rec = "✅ 表型高度匹配,建议确认"
        elif score >= 0.40:
            rec = "⚠️ 部分匹配,建议结合临床评估"
        else:
            rec = "⚠️ 表型匹配度低,建议复核"

        lines.append(f"| {v.gene} | {pos} | {gt} | {vaf} | {score:.2f} | {level_icon} {level_text} | {pg_status} | {rec} |\n")

    lines.append("\n")

    # Detailed per-variant analysis
    lines.append("### 逐变异详细分析\n\n")
    for i, v in enumerate(pheno_variants, 1):
        pos = f"{v.chrom}:{v.pos}"
        score = v.phenotype_match_score
        conf = v.phenotype_match_confidence or "low"

        if score >= 0.75:
            level_icon = "🟢"
            level_text = "高度关联"
        elif score >= 0.40:
            level_icon = "🟡"
            level_text = "中度关联"
        elif score > 0:
            level_icon = "🔴"
            level_text = "低度关联"
        else:
            level_icon = "⚪"
            level_text = "无关联"

        lines.append(f"**{i}. {v.gene} ({pos})** {level_icon} {level_text}\n\n")
        lines.append(f"- **变异**: {v.hgvsp or v.hgvsc or 'N/A'}\n")
        lines.append(f"- **表型匹配评分**: {score:.2f}/1.0 (置信度: {conf})\n")

        if v.phenotype_match_explanation:
            lines.append(f"- **匹配解释**: {v.phenotype_match_explanation}\n")

        if v.phenotype_matched_pairs:
            pairs_str = ", ".join([f"'{u}' → '{k}'" for u, k in v.phenotype_matched_pairs[:5]])
            ellipsis = "..." if len(v.phenotype_matched_pairs) > 5 else ""
            lines.append(f"- **匹配对**: {pairs_str}{ellipsis}\n")

        if v.phenotype_known_list:
            known_str = ", ".join(v.phenotype_known_list[:5])
            ellipsis = "..." if len(v.phenotype_known_list) > 5 else ""
            lines.append(f"- **基因已知表型**: {known_str}{ellipsis}\n")

        # ClinVar Pathogenic + low score → explicit warning
        clinvar_is_pathogenic = v.clinvar and any(kw in v.clinvar.lower() for kw in ("pathogenic", "致病"))
        if clinvar_is_pathogenic and score < 0.6:
            lines.append(f"- ⚠️ **ClinVar 致病性提示**: 该变异在 ClinVar 中标注为致病/可能致病,但与输入表型的匹配度较低({score:.2f})。")
            lines.append("建议结合患者具体临床表现进一步验证,不排除该变异与其他未报告表型相关。\n")

        # Tier context
        lines.append(f"- **当前分级**: Tier {v.tier} | 置信度: {v.tier_confidence or 'N/A'}\n")
        lines.append(f"- **分级原因**: {v.tier_reason}\n")

        lines.append("\n")

    # High-score variants validation recommendation
    high_score = [v for v in pheno_variants if v.phenotype_match_score >= 0.75]
    if high_score:
        lines.append("### ⚠️ 高分匹配变异验证建议\n\n")
        lines.append(f"检测到 **{len(high_score)}** 个表型高度匹配变异(评分≥0.75)。")
        lines.append("这些变异与输入的临床表型高度吻合,强烈建议进行验证:\n\n")
        for v in high_score:
            pos = f"{v.chrom}:{v.pos}"
            lines.append(f"- **{v.gene}**: {pos} ({v.hgvsp or v.hgvsc or 'N/A'}) - 评分 {v.phenotype_match_score:.2f}\n")
        lines.append("\n**建议验证方法**:\n")
        lines.append("- Sanger 测序验证\n")
        lines.append("- 长读长测序(如 PacBio / Oxford Nanopore)验证结构变异\n")
        lines.append("- 家系共分离分析(如可获得父母/子女样本)\n")
        lines.append("- 功能实验(如可行)\n\n")

    return "".join(lines)



# =============================================================================
# Report Generation
# =============================================================================

def generate_tier_report(variants: List[Variant], config: GPAConfig,
                        tissue_profile: Dict, multi_hits: List[Dict]) -> str:
    """
    Generate Markdown report with three-tier structure and dynamic tissue context.
    """
    # Sort by tier
    tier1 = [v for v in variants if v.tier == 1]
    tier2 = [v for v in variants if v.tier == 2]
    tier3 = [v for v in variants if v.tier == 3]

    profile_name = tissue_profile.get("display_name", config.tissue_profile)

    # v0.10.0: Build header via Jinja2 template
    version_info = _get_version_info(config)
    fs = config.filter_stats
    impact_parts = []
    if fs:
        for imp in ['HIGH', 'MODERATE', 'LOW', 'MODIFIER']:
            cnt = fs.get('by_impact', {}).get(imp, 0)
            if cnt > 0:
                impact_parts.append(f"{imp}: {cnt}")
    retention_parts = []
    if fs:
        if fs.get('splice_retained', 0) > 0:
            retention_parts.append(f"splice retained: {fs['splice_retained']}")
        if fs.get('synonymous_tissue_retained', 0) > 0:
            retention_parts.append(f"synonymous tissue: {fs['synonymous_tissue_retained']}")
        if fs.get('clinvar_conflicting_retained', 0) > 0:
            retention_parts.append(f"ClinVar conflict: {fs['clinvar_conflicting_retained']}")

    if _HAS_JINJA2 and _REPORT_HEADER_TEMPLATE is not None:
        header_md = _REPORT_HEADER_TEMPLATE.render(
            version=__version__,
            profile_name=profile_name,
            tissue_profile=config.tissue_profile,
            offline_mode=config.offline_mode,
            analysis_date=datetime.now().isoformat(),
            filter_stats=fs,
            impact_distribution=" | ".join(impact_parts) if impact_parts else None,
            filter_retention=" | ".join(retention_parts) if retention_parts else None,
            filter_preset=config.filter_preset,
            total_variants=len(variants),
            version_info=version_info,
            tier1_genes=len(set(v.gene for v in tier1)),
            tier1_variants=len(tier1),
            tier2_genes=len(set(v.gene for v in tier2)),
            tier2_variants=len(tier2),
            tier3_genes=len(set(v.gene for v in tier3)),
            tier3_variants=len(tier3),
        )
    else:
        header_md = _build_report_header(
            profile_name=profile_name,
            config=config,
            variants=variants,
            tier1=tier1,
            tier2=tier2,
            tier3=tier3,
        )

    report = [header_md]

    # v0.5 P1-13: QC summary table
    # Collect QC flags from all variants
    all_qc_flags = []
    for v in variants:
        all_qc_flags.extend(v.qc_flags)

    if all_qc_flags:
        from collections import Counter
        flag_counts = Counter(all_qc_flags)
        report.append("## ⚠️ 输入 QC 异常汇总\n")
        report.append(f"**总异常数**: {len(all_qc_flags)} 条(涉及 {len(set((v.chrom, v.pos) for v in variants if v.qc_flags))} 个变异)\n\n")
        report.append("| QC 标志 | 计数 | 说明 |\n")
        report.append("|---------|------|------|\n")

        flag_descriptions = {
            "INVALID_VAF": "VAF 超出 [0,1] 范围",
            "LOW_DEPTH": "测序深度 < 10x",
            "LOW_COMPLEXITY_REGION": "位于低复杂度/重复区域",
            "INVALID_GENE_SYMBOL": "基因名格式不符合 HGNC 规范",
            "VAF_GT_MISMATCH": "VAF 与基因型不一致(杂合 35-65%,纯合 85-100%,野生型 0-5%)",
        }

        # Show top 5 most frequent flags
        for flag, count in flag_counts.most_common(5):
            desc = flag_descriptions.get(flag, "未知标志")
            report.append(f"| {flag} | {count} | {desc} |\n")

        if len(flag_counts) > 5:
            report.append(f"| ... | ... | 共 {len(flag_counts)} 种异常 |\n")

        # Show first 5 flagged variant details
        flagged = [v for v in variants if v.qc_flags][:5]
        if flagged:
            report.append(f"\n**前 5 条异常变异**:\n")
            report.append("| 位置 | 基因 | 异常标志 |\n")
            report.append("|------|------|----------|\n")
            for v in flagged:
                report.append(f"| {v.chrom}:{v.pos} | {v.gene} | {', '.join(v.qc_flags)} |\n")

        if len([v for v in variants if v.qc_flags]) > 5:
            report.append(f"\n*... 共 {len([v for v in variants if v.qc_flags])} 个变异有异常,详见 JSON 输出*\n")

        report.append("\n")

    # v0.6: Pseudogene interference assessment - independent section
    pseudogene_assessment = _generate_pseudogene_assessment_section(variants)
    if pseudogene_assessment:
        report.append(pseudogene_assessment)
        report.append("\n")

    # Multi-hit warnings
    if multi_hits:
        report.append("## ⚠️ Multi-Hit Gene Warnings\n")
        for mh in multi_hits:
            report.append(f"### {mh['gene']} - {mh['variant_count']} variants detected\n")
            report.append(f"- **Warning**: {mh['warning']}\n")

            # v0.10.3: Phase analysis — concise display
            phase = mh.get('phase_result', {})
            if phase:
                status = phase.get('status', 'unknown')
                confidence = phase.get('confidence', 'unknown')
                method = phase.get('method', 'unknown')
                evidence = phase.get('evidence', 'N/A')
                max_gap = phase.get('max_gap_bp', 'N/A')

                # Phase status emoji
                if status in ('cis', 'cis_both', 'cis_likely'):
                    status_icon = '🟢'
                elif status == 'trans':
                    status_icon = '🔴'
                elif status == 'ambiguous':
                    status_icon = '🟡'
                elif status == 'unphased':
                    status_icon = '⚪'
                else:
                    status_icon = '❓'

                # v0.10.3: Only show detailed phase analysis when there is a real result.
                # UNPHASED / infeasible_short_reads means no phasing was possible —
                # show a one-liner instead of filling the report with boilerplate.
                _has_real_phase = status not in ('unphased', 'unknown') and 'infeasible' not in method

                if _has_real_phase:
                    report.append(f"\n**相位分析**: {status_icon} **{status.upper()}** (置信度: {confidence})\n")
                    report.append(f"- **判定方法**: {method}\n")
                    report.append(f"- **间距**: {max_gap}bp\n")
                    report.append(f"- **证据**: {evidence}\n")
                    phase_clinical = mh.get('phase_clinical_significance', '')
                    if phase_clinical:
                        report.append(f"- **临床意义**: {phase_clinical}\n")
                    report.append(f"\n- **Cis hypothesis**: {mh['phases']['cis']}\n")
                    report.append(f"- **Trans hypothesis**: {mh['phases']['trans']}\n")
                    report.append(f"- **Required evidence**: {', '.join(mh['required_evidence'])}\n")
                    report.append(f"- **Action**: {mh['action']}\n\n")
                else:
                    # Concise one-liner for unphased / infeasible
                    report.append(f"\n**相位分析**: {status_icon} **{status.upper()}** — {evidence}\n\n")

    # v0.7 Phase 4: Phenotype association assessment - independent section
    phenotype_assessment = _generate_phenotype_assessment_section(variants)
    if phenotype_assessment:
        report.append(phenotype_assessment)
        report.append("\n")

    # v0.9.0: Transcript selection assessment - independent section (ambiguous cases)
    tx_selection = _generate_transcript_selection_section(variants)
    if tx_selection:
        report.append(tx_selection)
        report.append("\n")

    # v0.10.14: Mandatory transcript validation for ALL Tier 1/2 variants
    tx_validation = _generate_tier12_transcript_validation_section(variants)
    if tx_validation:
        report.append(tx_validation)
        report.append("\n")

    # v0.10.3: Tier 1 + Tier 2 merged into a gene-centric view
    tier12 = tier1 + tier2
    if tier12:
        report.append("---\n\n## 重要变异分析（Tier 1/2）\n")
        report.append(f"*基因-重要突变聚合视图 for {profile_name} context*\n\n")

        tier12_genes = set(v.gene for v in tier12)
        report.append(f"**涉及基因总数**: {len(tier12_genes)} 个 | **变异总数**: {len(tier12)} 个")
        if tier1:
            report.append(f" (Tier 1: {len(tier1)}, Tier 2: {len(tier2)})")
        else:
            report.append(f" (全为 Tier 2)")
        report.append("\n\n")

        # List multi-hit genes among Tier 1/2
        multi_hit_tier12_genes = tier12_genes.intersection(set(mh['gene'] for mh in multi_hits))
        if multi_hit_tier12_genes:
            report.append(f"**其中 Multi-hit 基因** ({len(multi_hit_tier12_genes)} 个): {', '.join(sorted(multi_hit_tier12_genes))}\n")
            report.append(f"*注: Multi-hit 基因因检测到多个变异被标记关注,但各变异保持独立分级*\n\n")

        # Group by gene (merge tier1 + tier2 within each gene)
        from collections import OrderedDict
        gene_groups = OrderedDict()
        for v in tier12:
            gene_groups.setdefault(v.gene, []).append(v)

        for gene, var_list in gene_groups.items():
            # Count per tier within this gene
            t1_count = sum(1 for v in var_list if v.tier == 1)
            t2_count = sum(1 for v in var_list if v.tier == 2)
            tier_badge = ""
            if t1_count > 0:
                tier_badge += f" 🔴T1×{t1_count}"
            if t2_count > 0:
                tier_badge += f" 🟡T2×{t2_count}"

            report.append(f"### {gene}{tier_badge}\n")

            # Multi-hit indicator
            is_multi_hit = gene in [mh['gene'] for mh in multi_hits]
            if is_multi_hit:
                report.append(f"**[Multi-hit 基因]** | **变异数**: {len(var_list)}\n\n")
            else:
                report.append(f"**基因**: {gene} | **变异数**: {len(var_list)}\n\n")

            # Variant table with tier column
            report.append("| # | Tier | 染色体位置 | 转录本 | 变异名称 | 功能域 | 合子型 | ClinVar | 基因约束 | 置信度 | 说明 |\n")
            report.append("|---|------|-----------|--------|---------|--------|--------|---------|----------|--------|------|\n")

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

                # Tier badge
                tier_badge = "🔴" if v.tier == 1 else "🟡" if v.tier == 2 else "⚪"

                # Gene constraint
                gc_info = ""
                if v.gene_constraint:
                    pLI = v.gene_constraint.get("pLI")
                    loeuf = v.gene_constraint.get("loeuf")
                    parts = []
                    if pLI is not None:
                        parts.append(f"pLI={pLI:.2f}")
                    if loeuf is not None:
                        parts.append(f"LOEUF={loeuf:.2f}")
                    if parts:
                        gc_info = " | ".join(parts)

                conf = v.tier_confidence or "UNKNOWN"
                conf_icon = "⚠️" if conf == "LOW" else ""

                report.append(f"| {i} | {tier_badge} | {pos} | {tx} | {var_name} | {domain} | {zyg} | {clin} | {gc_info} | {conf_icon} {conf} | {reason} |\n")

            report.append(f"\n**详细说明**:\n")
            for i, v in enumerate(var_list, 1):
                tier_label = "🔴 Tier 1" if v.tier == 1 else "🟡 Tier 2"
                report.append(f"{i}. **{v.hgvsp or v.hgvsc}** ({v.chrom}:{v.pos}) — {tier_label}:\n")
                vep_note = _format_vep_reannotation_note(v)
                if vep_note:
                    report.append(f"   - **{vep_note}**\n")
                if "VAF_GT_MISMATCH" in v.qc_flags:
                    report.append(f"   - ⚠️ **VAF 与基因型不一致**(GT={v.gt}, VAF={v.vaf:.2f}),提示可能存在假基因干扰、CNV 或比对错误\n")
                if v.pseudogene_warning:
                    try:
                        pw = json.loads(v.pseudogene_warning)
                        if pw.get("type") == "PSEUDOGENE_INTERFERENCE":
                            report.append(f"   - ⚠️ **假基因干扰**: {pw.get('gene')} 观察 VAF={pw.get('observed_vaf', 'N/A')},远低于预期杂合 0.5,疑似 {', '.join(pw.get('pseudogenes', []))} 假基因读取\n")
                    except (json.JSONDecodeError, TypeError):
                        pass
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
                if v.gene_constraint:
                    gc = v.gene_constraint
                    parts = []
                    if gc.get("pLI") is not None:
                        parts.append(f"pLI={gc['pLI']:.2f}")
                    if gc.get("loeuf") is not None:
                        parts.append(f"LOEUF={gc['loeuf']:.2f}")
                    if gc.get("lof_z") is not None:
                        parts.append(f"lof_z={gc['lof_z']:.2f}")
                    if gc.get("mis_z") is not None:
                        parts.append(f"mis_z={gc['mis_z']:.2f}")
                    if parts:
                        report.append(f"   - 基因约束: {' | '.join(parts)} (来源: {gc.get('source', 'N/A')})\n")
                report.append(f"   - 分级原因: {v.tier_reason}\n")
                if v.phenotype_match_score is not None:
                    score = v.phenotype_match_score
                    conf = v.phenotype_match_confidence
                    conf_icon = "🟢" if conf == "high" else "🟡" if conf == "medium" else "🔴"
                    report.append(f"   - **表型关联**: {conf_icon} Score={score:.2f} (置信度: {conf})\n")
                    if v.phenotype_match_explanation:
                        report.append(f"     - 解释: {v.phenotype_match_explanation}\n")
                    if v.phenotype_matched_pairs:
                        pairs_str = ", ".join([f"'{u}'→'{k}'" for u, k in v.phenotype_matched_pairs[:3]])
                        report.append(f"     - 匹配对: {pairs_str}{'...' if len(v.phenotype_matched_pairs) > 3 else ''}\n")
                    if v.phenotype_known_list:
                        known_str = ", ".join(v.phenotype_known_list[:3])
                        report.append(f"     - 基因已知表型: {known_str}{'...' if len(v.phenotype_known_list) > 3 else ''}\n")
                if v.evidence_chain:
                    report.append(f"   - **证据链** ({len(v.evidence_chain)} 条):\n")
                    chain = v.evidence_chain
                    evidence_detail = getattr(config, 'evidence_detail', 'brief')
                    if evidence_detail == 'brief' and len(chain) > 3:
                        chain = chain[:3]
                        report.append(f"     *(brief mode: 显示前 3 条关键证据)*\n")
                    for ev in chain:
                        report.append(f"     - **{ev.source}**: {ev.rule} (权重={ev.weight:.2f}, 置信度={ev.confidence})\n")
                        if evidence_detail == 'full' and ev.raw_data:
                            raw_summary = {k: v for k, v in ev.raw_data.items() if k in ('af', 'tpm', 'pLI', 'loeuf', 'status', 'domain')}
                            if raw_summary:
                                report.append(f"       原始数据: {raw_summary}\n")
                if v.tier_actions:
                    report.append(f"   - 建议措施: {'; '.join(v.tier_actions)}\n")
                if v.upgrade_conditions:
                    report.append(f"   - **升级条件**:\n")
                    for uc in v.upgrade_conditions:
                        report.append(f"     → {uc}\n")
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
                report.append("| 位置 | 变异 | 功能域 | 置信度 | 原因 |\n")
                report.append("|------|------|--------|--------|------|\n")
                for v in var_list:
                    pos = f"{v.chrom}:{v.pos}"
                    var_name = v.hgvsp or v.hgvsc or "N/A"
                    di = v.domain_info
                    domain = f"{di.get('domain', 'N/A')}" if di else "N/A"
                    # v0.5 P1-10: Confidence
                    conf = v.tier_confidence or "UNKNOWN"
                    conf_icon = "⚠️" if conf == "LOW" else ""
                    reason = v.tier_reason[:50] + "..." if len(v.tier_reason) > 50 else v.tier_reason
                    reason = reason.replace("|", "/")
                    # v0.5.2: VEP reannotation note in table
                    vep_note = _format_vep_reannotation_note(v)
                    if vep_note:
                        reason = f"⚠️ VEP reannotated | {reason}"
                    report.append(f"| {pos} | {var_name} | {domain} | {conf_icon} {conf} | {reason} |\n")
                report.append(f"\n")
                # v0.5.2: VEP reannotation detail for Tier 3 short-list
                for v in var_list:
                    vep_note = _format_vep_reannotation_note(v)
                    if vep_note:
                        report.append(f"   *{v.hgvsp or v.hgvsc}*: {vep_note}\n")
                report.append(f"\n")
            else:
                # Many variants: just count
                report.append(f"**{gene}**: {len(var_list)} variants - 详见原始数据\n\n")

            # v0.5 P1-11: Show upgrade conditions for Tier 3 short-list genes
            if len(var_list) <= 3:
                for v in var_list:
                    if v.upgrade_conditions:
                        report.append(f"   *{v.hgvsp or v.hgvsc} 升级条件*: {' / '.join(v.upgrade_conditions[:2])}\n")
                report.append(f"\n")

    # Methodology
    report.append("---\n\n## 方法学附录\n")
    report.append(f"### 分析背景: {profile_name}\n")
    report.append(f"- **GTEx 参考组织**: {tissue_profile.get('gtex_tissue', 'N/A')}\n")
    report.append(f"- **快速排除规则**: {tissue_profile.get('fast_track_rule', 'N/A')}\n\n")
    report.append("### 分析流程\n")
    report.append("1. **转录本校正**: Ensembl REST API → canonical transcript → 本地回退\n")
    report.append("2. **假基因检测**: VAF 偏差分析识别已知假基因对\n")
    report.append("3. **gnomAD 整合**: AF>1% 常见; AF<0.1% 罕见; NOT_CAPTURED 明确标注\n")
    report.append("4. **蛋白功能域映射**: UniProt REST API → DOMAIN/REGION 特征 → 本地回退\n")
    report.append("5. **组织相关性评估**: GTEx API → median TPM → 自动分级 + 本地回退\n")
    report.append("6. **表型关联分析** (v0.7): LLM 语义匹配 → 基因已知表型 vs 用户输入表型 → 仅 Tier 1/2 执行\n")
    report.append("7. **三级分类**: Action (Tier 1) → Inform (Tier 2) → No concern (Tier 3)\n")
    report.append("8. **基因检测历史记录**: 保留变异信息用于后续分析\n")
    report.append("9. **缓存**: 所有 API 响应缓存 30 天 (SQLite); 离线模式仅用缓存\n")

    return "\n".join(report)

# =============================================================================
# JSON Structured Report Generation (v0.5 P1-12)
# =============================================================================

def generate_json_report(variants: List[Variant], config: GPAConfig,
                         tissue_profile: Dict, multi_hits: List[Dict],
                         report_md: str,
                         qc_summary: Optional[Dict] = None) -> Dict:
    """
    Generate structured JSON report for downstream system consumption.
    v0.5 P1-12: Complete structured output alongside Markdown report.
    """
    from dgra_core import _UNKNOWN

    profile_name = tissue_profile.get("display_name", config.tissue_profile)

    # Meta section - v0.5 P1-15: include full version metadata
    meta = {
        "dgra_version": "0.9.0",
        "analysis_date": datetime.now().isoformat(),
        "input_format": "vcf",
        "tissue_profile": config.tissue_profile,
        "target_population": getattr(config, 'target_population', None),
        "scoring_model": "weighted",
        "evidence_detail": getattr(config, 'evidence_detail', 'brief'),
        "offline_mode": config.offline_mode,
        "somatic_mode": config.somatic_mode,
        "multi_organ_profiles": config.multi_organ_profiles,
    }
    # Merge version/provenance info (P1-15)
    meta.update(_get_version_info(config))

    # Summary section - v0.5.2: gene-level and variant-level counts
    tier1_genes = set(v.gene for v in variants if v.tier == 1)
    tier2_genes = set(v.gene for v in variants if v.tier == 2)
    tier3_genes = set(v.gene for v in variants if v.tier == 3)
    summary = {
        "total_variants": len(variants),
        "tier1_gene_count": len(tier1_genes),
        "tier1_variant_count": len([v for v in variants if v.tier == 1]),
        "tier2_gene_count": len(tier2_genes),
        "tier2_variant_count": len([v for v in variants if v.tier == 2]),
        "tier3_gene_count": len(tier3_genes),
        "tier3_variant_count": len([v for v in variants if v.tier == 3]),
        "multi_hit_genes": [mh["gene"] for mh in multi_hits],
    }

    # Variants array - structured per variant
    variants_json = []
    for v in variants:
        # Build evidence chain JSON
        evidence_chain = []
        for ev in v.evidence_chain:
            evidence_chain.append({
                "source": ev.source,
                "rule": ev.rule,
                "weight": ev.weight,
                "confidence": ev.confidence,
                "raw_data": ev.raw_data,
            })

        # Build gnomAD section
        gnomad_section = {
            "overall_af": v.gnomad_af,
            "popmax_af": None,
            "eas_af": None,
        }
        if v.gnomad_populations:
            for pop_code, pop_data in v.gnomad_populations.items():
                if pop_code == "EAS":
                    gnomad_section["eas_af"] = pop_data.get("af")
                # Track popmax
                pop_af = pop_data.get("af")
                if pop_af is not None:
                    current_max = gnomad_section["popmax_af"] or 0
                    if pop_af > current_max:
                        gnomad_section["popmax_af"] = pop_af

        # Build GTEx section
        gtex_section = {"tissue": None, "tpm_median": None}
        if v.tissue_relevance:
            gtex_section["tissue"] = v.tissue_relevance.get("gtex_tissue")
            gtex_section["tpm_median"] = v.tissue_relevance.get("gtex_tpm")

        # Build NMD prediction
        nmd_section = {"status": "not_applicable"}
        if v.gene_constraint and "nmd_prediction" in v.gene_constraint:
            nmd_section = v.gene_constraint["nmd_prediction"]

        # Build ClinVar section
        clinvar_section = {
            "clinical_significance": v.clinvar if v.clinvar != _UNKNOWN else None,
            "review_status": None,
        }

        variant_json = {
            "gene": v.gene,
            "gene_original": v.gene_original or v.gene,
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "hgvsc": v.hgvsc or None,
            "hgvsp": v.hgvsp or None,
            "transcript": v.transcript or None,
            "exon": v.exon or None,
            "impact": v.impact if v.impact != _UNKNOWN else None,
            "consequence": v.consequence if v.consequence != _UNKNOWN else None,
            "zygosity": v.gt or None,
            "vaf": v.vaf,
            "dp": v.dp,
            "gq": v.gq,
            "tier": v.tier,
            "tier_confidence": v.tier_confidence,
            "tier_reason": v.tier_reason or None,
            "tier_actions": v.tier_actions,
            "evidence_chain": evidence_chain,
            "upgrade_conditions": v.upgrade_conditions,
            "gene_constraint": v.gene_constraint,
            "clinvar": clinvar_section,
            "gnomAD": gnomad_section,
            "gtex": gtex_section,
            "nmd_prediction": nmd_section,
            "domain_info": v.domain_info,
            "tissue_relevance": v.tissue_relevance,
            "quality_confidence": v.quality_confidence,
            "missing_fields": v.missing_fields,
            "pseudogene_warning": json.loads(v.pseudogene_warning) if v.pseudogene_warning else None,
            "transcript_warning": json.loads(v.transcript_warning) if v.transcript_warning else None,
        }
        variants_json.append(variant_json)

    # v0.7 Phase 4: Phenotype association data for JSON output
    phenotype_data = []
    pheno_variants = [v for v in variants if v.tier in (1, 2) and v.phenotype_match_score is not None]
    for v in pheno_variants:
        phenotype_data.append({
            "gene": v.gene,
            "chrom": v.chrom,
            "pos": v.pos,
            "ref": v.ref,
            "alt": v.alt,
            "hgvsp": v.hgvsp or None,
            "hgvsc": v.hgvsc or None,
            "tier": v.tier,
            "phenotype_match_score": v.phenotype_match_score,
            "phenotype_match_confidence": v.phenotype_match_confidence or None,
            "phenotype_match_explanation": v.phenotype_match_explanation or None,
            "phenotype_matched_pairs": v.phenotype_matched_pairs,
            "phenotype_known_list": v.phenotype_known_list,
        })

    phenotype_association = {
        "total_tier12_with_phenotype": len(pheno_variants),
        "high_match_count": len([v for v in pheno_variants if v.phenotype_match_score >= 0.75]),
        "variants": phenotype_data,
    }

    # Assemble final JSON
    json_report = {
        "meta": meta,
        "summary": summary,
        "variants": variants_json,
        "multi_hit_details": multi_hits,
        "qc_summary": qc_summary,  # v0.5 P1-13: input QC flags
        "phenotype_association": phenotype_association,
        "report_md": report_md,
    }

    return json_report


