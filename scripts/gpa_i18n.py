#!/usr/bin/env python3
"""
gpa_i18n.py — GPA Internationalization / Normalization Module (v0.7.1)

Centralized bilingual (Chinese/English) consequence normalization.
Used by adapters, filters, and downstream analysis modules.
"""

from typing import List

# =============================================================================
# Chinese VEP CSV Column Header Mapping (v0.9.1 hotfix)
# =============================================================================

CHINESE_COLUMN_MAP = {
    # 位置信息
    "Uploaded_variation": "Uploaded_variation",
    "位置": "Location",
    "基因位置": "Gene",
    "基因": "Gene",
    "特征": "Feature",
    "转录本": "Feature",
    "转录本ID": "Feature",
    "转录本编号": "Feature",

    # 后果
    "变异后果": "Consequence",
    "后果": "Consequence",
    "变异类型": "Consequence",
    "影响程度": "IMPACT",
    "影响": "IMPACT",

    # cDNA/蛋白
    "CDNA位置": "cDNA_position",
    "CDS位置": "CDS_position",
    "蛋白位置": "Protein_position",
    "氨基酸改变": "Amino_acids",
    "密码子": "Codons",
    "HGVC": "HGVSc",
    "HGVSc": "HGVSc",
    "HGVSp": "HGVSp",
    "hgvsc": "HGVSc",
    "hgvsp": "HGVSp",

    # 等位基因
    "现有等位基因": "Existing_variation",
    "rs号": "Existing_variation",
    "参考等位基因": "REF",
    "替代等位基因": "ALT",

    # 频率
    "gnomAD频率": "gnomAD_AF",
    "gnomad频率": "gnomAD_AF",
    "gnomad_af": "gnomAD_AF",
    "GnomAD_AF": "gnomAD_AF",
    "千人基因组频率": "AFR_AF",
    "最大人群频率": "MAX_AF",

    # 致病性
    "ClinVar": "CLIN_SIG",
    "临床意义": "CLIN_SIG",
    "clinvar": "CLIN_SIG",

    # 样本
    "样本": "SAMPLE",
    "基因型": "GT",
    "测序深度": "DP",
    "质量值": "GQ",
    "等位基因频率": "VAF",

    # 其他
    "距离": "DISTANCE",
    "链": "STRAND",
    "突变频谱": "VARIANT_CLASS",
    "最小等位基因": "Allele",
    "SYMBOL": "SYMBOL",
    "基因符号": "SYMBOL",
}


def _normalize_header_key(key: str) -> str:
    """Normalize header key for fuzzy matching."""
    return key.strip().replace(" ", "").replace("_", "").replace("-", "").lower()


def translate_chinese_header(headers: List[str]) -> List[str]:
    """将中文VEP CSV表头翻译为英文标准列名。"""
    translated = []
    # Build normalized lookup
    normalized_map = {_normalize_header_key(k): v for k, v in CHINESE_COLUMN_MAP.items()}
    for h in headers:
        h_stripped = h.strip()
        # Exact match first
        if h_stripped in CHINESE_COLUMN_MAP:
            translated.append(CHINESE_COLUMN_MAP[h_stripped])
            continue
        # Fuzzy match
        h_norm = _normalize_header_key(h_stripped)
        if h_norm in normalized_map:
            translated.append(normalized_map[h_norm])
            continue
        # Keep original if no match
        translated.append(h_stripped)
    return translated


def is_chinese_header(headers: List[str]) -> bool:
    """检测表头是否包含中文字符。"""
    for h in headers:
        if any("\u4e00" <= ch <= "\u9fff" for ch in h):
            return True
    return False


# =============================================================================
# Consequence Term Mapping: Chinese VEP/SnpEff/ANNOVAR → English SO terms
# =============================================================================

CONSEQUENCE_MAP = {
    # High impact
    "移码变异": "frameshift_variant",
    "获得终止密码子": "stop_gained",
    "失去终止密码子": "stop_lost",
    "起始密码子缺失": "start_lost",
    "剪接受体位点变异": "splice_acceptor_variant",
    "剪接供体位点变异": "splice_donor_variant",
    "剪接位点变异": "splice_variant",
    "转录本缺失": "transcript_ablation",
    "转录本扩增": "transcript_amplification",
    # Moderate impact
    "错义变异": "missense_variant",
    "框内插入": "inframe_insertion",
    "框内缺失": "inframe_deletion",
    "蛋白质改变变异": "protein_altering_variant",
    "剪接区域变异": "splice_region_variant",
    "剪接区变异": "splice_region_variant",
    "剪接多聚嘧啶束变异": "splice_polypyrimidine_tract_variant",
    "剪接供体第5碱基变异": "splice_donor_5th_base_variant",
    "起始密码子保留变异": "start_retained_variant",
    "终止密码子保留变异": "stop_retained_variant",
    "不完全末端密码子变异": "incomplete_terminal_codon_variant",
    # Low impact
    "同义变异": "synonymous_variant",
    "5'UTR变异": "5_prime_UTR_variant",
    "3'UTR变异": "3_prime_UTR_variant",
    "内含子变异": "intron_variant",
    "上游基因变异": "upstream_gene_variant",
    "下游基因变异": "downstream_gene_variant",
    "非编码转录本外显子变异": "non_coding_transcript_exon_variant",
    "编码序列变异": "coding_sequence_variant",
    "成熟miRNA变异": "mature_miRNA_variant",
    "NMD转录本变异": "NMD_transcript_variant",
    "非编码转录本变异": "non_coding_transcript_variant",
    # Modifier
    "基因间变异": "intergenic_variant",
    "调控区域变异": "regulatory_region_variant",
    "TFBS缺失": "TFBS_ablation",
    "TFBS扩增": "TFBS_amplification",
    "TF结合位点变异": "TF_binding_site_variant",
    "调控区域缺失": "regulatory_region_ablation",
    "调控区域扩增": "regulatory_region_amplification",
    "特征延伸": "feature_elongation",
    "特征截短": "feature_truncation",
}


# =============================================================================
# Normalization Functions
# =============================================================================

def normalize_consequence(cons_str: str) -> List[str]:
    """
    Normalize a mixed Chinese/English consequence string into English SO terms.

    Input:  "错义变异,剪接区域变异"  or  "missense_variant,splice_region_variant"
    Output: ["missense_variant", "splice_region_variant"]
    """
    if not cons_str or cons_str.strip() in ("", ".", "UNKNOWN", "N/A"):
        return []

    # Split by common delimiters (VEP uses ',', some tools use '&' or ';')
    raw_terms = [t.strip() for t in cons_str.replace("&", ",").replace(";", ",").split(",")]

    normalized = []
    for t in raw_terms:
        if not t:
            continue
        # Direct Chinese → English mapping
        if t in CONSEQUENCE_MAP:
            normalized.append(CONSEQUENCE_MAP[t])
        else:
            # Already English or unknown — pass through as-is
            normalized.append(t)
    return normalized


def infer_impact_from_consequence(cons_str: str) -> str:
    """
    Unified IMPACT inference from consequence string.
    Supports both Chinese and English input.

    Returns: "HIGH" | "MODERATE" | "LOW" | "MODIFIER" | "" (unknown)
    """
    if not cons_str or cons_str.strip() in ("", ".", "UNKNOWN", "N/A"):
        return ""  # core.py maps empty → UNKNOWN → conservative HIGH

    # Normalize to English SO terms first
    terms = normalize_consequence(cons_str)
    if not terms:
        return ""

    # Combine all terms into a single lowercase string for matching
    cons_combined = ",".join(terms).lower()

    # HIGH impact consequences
    high_keywords = {
        "stop_gained", "stop_lost", "frameshift_variant", "splice_acceptor_variant",
        "splice_donor_variant", "start_lost", "transcript_ablation", "transcript_amplification",
        "splice_variant",
    }
    # MODERATE impact consequences
    moderate_keywords = {
        "missense_variant", "inframe_deletion", "inframe_insertion",
        "protein_altering_variant", "splice_region_variant",
        "splice_polypyrimidine_tract_variant", "splice_donor_5th_base_variant",
        "incomplete_terminal_codon_variant", "stop_retained_variant", "start_retained_variant",
    }
    # LOW impact consequences
    low_keywords = {
        "synonymous_variant", "5_prime_utr_variant", "3_prime_utr_variant",
        "intron_variant", "upstream_gene_variant", "downstream_gene_variant",
        "non_coding_transcript_exon_variant", "coding_sequence_variant",
        "mature_mirna_variant", "nmd_transcript_variant", "non_coding_transcript_variant",
    }

    if any(h in cons_combined for h in high_keywords):
        return "HIGH"
    if any(m in cons_combined for m in moderate_keywords):
        return "MODERATE"
    if any(l in cons_combined for l in low_keywords):
        return "LOW"

    # Everything else → MODIFIER (or pass through as unknown for conservative handling)
    return "MODIFIER"


# =============================================================================
# ClinVar text normalization (for downstream matching)
# =============================================================================

def normalize_clinvar(clinvar_str: str) -> str:
    """
    Normalize ClinVar text: standardize delimiters, strip extra whitespace.
    Does NOT change semantic meaning — purely formatting.
    """
    if not clinvar_str:
        return ""
    # Replace & with / (VEP style)
    normalized = clinvar_str.replace("_", " ").replace("&", "/")
    return normalized.strip()
