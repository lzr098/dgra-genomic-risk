#!/usr/bin/env python3
"""
GPA Phase Analysis Module

Determines cis/trans phase relationships between multiple variants in the same gene.
Used by multi-hit detection to assess compound heterozygosity.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import List, Dict, Optional, TYPE_CHECKING

if TYPE_CHECKING:
    from dgra_core import Variant


@dataclass
class PhaseResult:
    phase_status: str       # cis / trans / cis_both / ambiguous / unphased / cis_likely / trans_likely
    confidence: str         # high / medium / low / none
    method: str             # gatk_phased_gt / short_reads_overlap / paired_end / reads_direct / trio_segregation / ld_inference / infeasible_short_reads
    evidence: str           # 详细证据描述
    max_gap_bp: int = 0
    min_gap_bp: int = 0
    n_variants: int = 0


def _parse_gt_field(gt_str: str) -> Dict:
    """解析 VCF GT 字段,返回 {is_phased, allele_0, allele_1}"""
    gt_str = str(gt_str) if gt_str is not None else ""
    if not gt_str or gt_str in ('.', './.', '.|.', 'nan'):
        return {"is_phased": False, "allele_0": -1, "allele_1": -1}

    if '|' in gt_str:
        parts = gt_str.split('|')
        return {"is_phased": True, "allele_0": int(parts[0]), "allele_1": int(parts[1])}
    elif '/' in gt_str:
        parts = gt_str.split('/')
        return {"is_phased": False, "allele_0": int(parts[0]), "allele_1": int(parts[1])}
    else:
        # Single allele (haploid)
        val = int(gt_str)
        return {"is_phased": False, "allele_0": val, "allele_1": val}


def _level1_gatk_phase(variants: List[Variant]) -> Optional[PhaseResult]:
    """Level 1: 基于 GATK phased GT 判断相位"""
    parsed = [_parse_gt_field(v.gt) for v in variants]

    # 检查是否全部 phased
    all_phased = all(p["is_phased"] for p in parsed)
    if not all_phased:
        return None

    hap0_alleles = [p["allele_0"] for p in parsed]
    hap1_alleles = [p["allele_1"] for p in parsed]

    # 情况 1: 所有变异都是 1|1 → 两条单倍型都携带
    if set(hap0_alleles) == {1} and set(hap1_alleles) == {1}:
        return PhaseResult(
            phase_status="cis_both",
            confidence="high",
            method="gatk_phased_gt",
            evidence=f"所有 {len(variants)} 个变异 GT=1|1,GATK local assembly 确认两条单倍型均携带"
        )

    # 情况 2: Hap0 全为 ALT, Hap1 全为 REF → cis(杂合)
    if set(hap0_alleles) == {1} and set(hap1_alleles) == {0}:
        return PhaseResult(
            phase_status="cis",
            confidence="high",
            method="gatk_phased_gt",
            evidence="所有 ALT 等位基因位于同一单倍型 (Hap0),REF 位于另一单倍型"
        )

    # 情况 3: Hap0 全为 REF, Hap1 全为 ALT → cis(对称)
    if set(hap0_alleles) == {0} and set(hap1_alleles) == {1}:
        return PhaseResult(
            phase_status="cis",
            confidence="high",
            method="gatk_phased_gt",
            evidence="所有 ALT 等位基因位于同一单倍型 (Hap1),REF 位于另一单倍型"
        )

    # 情况 4: Hap0 上同时存在 REF 和 ALT → trans
    if 0 in hap0_alleles and 1 in hap0_alleles:
        return PhaseResult(
            phase_status="trans",
            confidence="high",
            method="gatk_phased_gt",
            evidence=f"单倍型 Hap0 上同时存在 REF 和 ALT ({hap0_alleles}),确认 trans 关系"
        )

    # 情况 5: Hap1 上同时存在 REF 和 ALT → trans
    if 0 in hap1_alleles and 1 in hap1_alleles:
        return PhaseResult(
            phase_status="trans",
            confidence="high",
            method="gatk_phased_gt",
            evidence=f"单倍型 Hap1 上同时存在 REF 和 ALT ({hap1_alleles}),确认 trans 关系"
        )

    # 其他情况(如包含缺失 -1)
    return None


def _level2_distance_assessment(variants: List[Variant]) -> Dict:
    """Level 2: 基于变异间距判断相位可行性"""
    positions = sorted([v.pos for v in variants])
    gaps = [positions[i+1] - positions[i] for i in range(len(positions) - 1)]
    max_gap = max(gaps) if gaps else 0
    min_gap = min(gaps) if gaps else 0

    if max_gap < 50:
        return {
            "feasible": True,
            "confidence": "high",
            "method": "short_reads_overlap",
            "evidence": f"间距 {min_gap}-{max_gap}bp,同一 150bp read 必然覆盖所有变异"
        }
    elif max_gap < 150:
        return {
            "feasible": True,
            "confidence": "high",
            "method": "short_reads_overlap_or_paired_end",
            "evidence": f"间距 {min_gap}-{max_gap}bp,同一 read (靠近 3' 端) 或 pair-end 覆盖"
        }
    elif max_gap < 500:
        return {
            "feasible": True,
            "confidence": "medium",
            "method": "paired_end_only",
            "evidence": f"间距 {min_gap}-{max_gap}bp,依赖 pair-end insert size (通常 300-500bp)"
        }
    else:
        return {
            "feasible": False,
            "confidence": "none",
            "method": "infeasible_short_reads",
            "evidence": f"最大间距 {max_gap}bp 超出 short-read 相位范围"
        }


def determine_phase(variants: List[Variant]) -> PhaseResult:
    """
    主函数:分层决策判断 multi-hit 变异的相位关系

    优先级:
    1. GATK phased GT(最可靠)
    2. 间距可行性判断(短 reads 范围评估)
    3. 标记为需进一步验证(trio / 长读长)
    """
    positions = sorted([v.pos for v in variants])
    max_gap = max(positions[i+1] - positions[i] for i in range(len(positions) - 1)) if len(positions) > 1 else 0
    min_gap = min(positions[i+1] - positions[i] for i in range(len(positions) - 1)) if len(positions) > 1 else 0

    # Level 1: GATK Phased GT
    result = _level1_gatk_phase(variants)
    if result:
        result.max_gap_bp = max_gap
        result.min_gap_bp = min_gap
        result.n_variants = len(variants)
        return result

    # Level 2: 间距可行性判断
    distance = _level2_distance_assessment(variants)

    if not distance["feasible"]:
        # 短 reads 不可行
        return PhaseResult(
            phase_status="unphased",
            confidence="none",
            method=distance["method"],
            evidence=f"{distance['evidence']}。建议: trio 测序 或 PacBio/Nanopore 长读长",
            max_gap_bp=max_gap,
            min_gap_bp=min_gap,
            n_variants=len(variants)
        )

    # 间距可行但未 phased
    # 根据间距范围给出 cis 可能性评估
    if distance["method"] == "short_reads_overlap":
        # <50bp: 同一 read 必然覆盖 → 高置信度 cis(如果都是杂合)
        # 但如果是 0/1 (unphased),我们无法确认是 cis 还是 trans
        # 只能标记为"技术上可行,需 reads 分析确认"
        return PhaseResult(
            phase_status="cis_likely" if all(_parse_gt_field(v.gt)["allele_0"] == _parse_gt_field(v.gt)["allele_1"] for v in variants) else "ambiguous",
            confidence="high",
            method=distance["method"],
            evidence=f"{distance['evidence']}。GATK 未输出 phased GT,但物理距离保证 reads 重叠。建议 IGV 验证 reads 直接比对",
            max_gap_bp=max_gap,
            min_gap_bp=min_gap,
            n_variants=len(variants)
        )
    elif distance["method"] == "short_reads_overlap_or_paired_end":
        return PhaseResult(
            phase_status="ambiguous",
            confidence="medium",
            method=distance["method"],
            evidence=f"{distance['evidence']}。短 reads 可能 phase,需 reads 分析或 trio 确认",
            max_gap_bp=max_gap,
            min_gap_bp=min_gap,
            n_variants=len(variants)
        )
    else:  # paired_end_only
        return PhaseResult(
            phase_status="ambiguous",
            confidence="low",
            method=distance["method"],
            evidence=f"{distance['evidence']}。pair-end phase 可靠性低,建议 trio 或长读长",
            max_gap_bp=max_gap,
            min_gap_bp=min_gap,
            n_variants=len(variants)
        )
