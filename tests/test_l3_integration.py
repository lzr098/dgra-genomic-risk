"""
L3 集成测试 — 端到端核心流程验证
原始 VCF → 注释 → 过滤 → 选择 → 分级
中文 CSV → 全 pipeline
gnomAD API 失败级联
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import run_tests, MockGnomAD, MockTissueProfile, MockTissueAssessment, make_variant


def test_integration_raw_vcf_to_tier_ddx3x_tier3():
    """原始 VCF 流程：DDX3X 常见变异 → Tier 3。"""
    from dgra_core import classify_variant_tier, GPAConfig, Variant
    from dgra_variant_filter import filter_variants

    # 模拟原始 VCF 解析后的变异 dict
    raw = {
        "CHROM": "X", "POS": "41357831", "REF": "A", "ALT": "T",
        "GENE": "DDX3X", "Feature": "ENST00000373383",
        "EXON": "5/15", "IMPACT": "HIGH", "Consequence": "stop_gained",
        "HGVSp": "p.Glu233Ter", "HGVSc": "c.697G>A",
        "CLIN_SIG": "", "GT": "0/1", "DP": "60", "GQ": "99", "VAF": "0.45",
        "gnomAD_AF": "0.45",
    }

    # 过滤
    filtered, stats = filter_variants([raw], preset="clinical")
    assert stats["output_count"] == 1

    # 构建 Variant (模拟 VCF annotator 后的状态)
    v = Variant(
        chrom=raw["CHROM"], pos=int(raw["POS"]), ref=raw["REF"], alt=raw["ALT"],
        gene=raw["GENE"], transcript=raw["Feature"], exon=raw["EXON"],
        impact=raw["IMPACT"], consequence=raw["Consequence"],
        hgvsp=raw["HGVSp"], hgvsc=raw["HGVSc"], clinvar=raw["CLIN_SIG"],
        gt=raw["GT"], dp=int(raw["DP"]), gq=float(raw["GQ"]),
        vaf=float(raw["VAF"]), gnomad_af=float(raw["gnomAD_AF"]),
        gnomad_status="SUCCESS",
    )

    # DDX3X 不是造血核心基因，用 none assessment 触发 fast-track
    tissue = MockTissueAssessment.none()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.success_common("X", 41357831, "A", "T", af=0.45)

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 3, f"Expected Tier 3 for DDX3X common variant, got {tier}. Reason: {reason}"


def test_integration_chinese_csv_tp53_tier1():
    """中文 CSV 流程：TP53 Pathogenic → Tier 1。"""
    from dgra_core import classify_variant_tier, GPAConfig, Variant
    from dgra_adapters import VEPAdapter
    from gpa_i18n import translate_chinese_headers

    # 先翻译表头
    cn_headers = [
        "位置", "基因", "转录本", "变异后果", "影响程度",
        "HGVSc", "HGVSp", "ClinVar", "gnomAD频率",
        "样本", "基因型", "测序深度", "质量值",
    ]
    en_headers = translate_chinese_headers(cn_headers)
    assert "CLIN_SIG" in en_headers

    # 模拟中文 CSV 一行（英文键，已翻译）
    raw_row = {
        "Location": "17_7579472_C_T",
        "Gene": "TP53",
        "Feature": "ENST00000269305",
        "Consequence": "stop_gained",
        "IMPACT": "HIGH",
        "HGVSc": "c.818C>T",
        "HGVSp": "p.Arg273Ter",
        "CLIN_SIG": "致病",
        "gnomAD_AF": "0.00001",
        "GT": "0/1",
        "DP": "80",
        "GQ": "99",
    }

    # VEPAdapter 适配（英文键输入）
    adapter = VEPAdapter()
    adapted = adapter.adapt(raw_row)
    assert adapted["IMPACT"] == "HIGH"
    assert adapted["Consequence"] == "stop_gained"

    # 构建 Variant
    v = Variant(
        chrom="17", pos=7579472, ref="C", alt="T",
        gene=adapted["GENE"], transcript=adapted["Feature"],
        exon="5/11", impact=adapted["IMPACT"],
        consequence=adapted["Consequence"],
        hgvsp=adapted["HGVSp"], hgvsc=adapted["HGVSc"],
        clinvar=adapted["CLIN_SIG"],
        gt="0/1", dp=80, gq=99.0, vaf=0.45,
        gnomad_af=0.00001, gnomad_status="SUCCESS",
    )

    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.success_rare("17", 7579472, "C", "T", af=0.00001)

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 1, f"Expected Tier 1 for TP53 Pathogenic, got {tier}. Reason: {reason}"


def test_integration_chinese_csv_ddx3x_tier3():
    """中文 CSV 流程：DDX3X 常见变异 → Tier 3。"""
    from dgra_core import classify_variant_tier, GPAConfig, Variant
    from dgra_adapters import VEPAdapter
    from gpa_i18n import translate_chinese_headers

    # 英文键输入
    raw_row = {
        "Location": "X_41357831_A_T",
        "Gene": "DDX3X",
        "Feature": "ENST00000373383",
        "Consequence": "stop_gained",
        "IMPACT": "HIGH",
        "HGVSc": "c.697G>A",
        "HGVSp": "p.Glu233Ter",
        "CLIN_SIG": "",
        "gnomAD_AF": "0.45",
        "GT": "0/1",
        "DP": "60",
        "GQ": "99",
    }

    adapter = VEPAdapter()
    adapted = adapter.adapt(raw_row)
    assert adapted["IMPACT"] == "HIGH"

    v = Variant(
        chrom="X", pos=41357831, ref="A", alt="T",
        gene="DDX3X", transcript="ENST00000373383",
        exon="5/15", impact="HIGH", consequence="stop_gained",
        hgvsp="p.Glu233Ter", hgvsc="c.697G>A",
        clinvar="", gt="0/1", dp=60, gq=99.0, vaf=0.45,
        gnomad_af=0.45, gnomad_status="SUCCESS",
    )

    # DDX3X 不是造血核心基因
    tissue = MockTissueAssessment.none()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.success_common("X", 41357831, "A", "T", af=0.45)

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 3, f"Expected Tier 3 for DDX3X common, got {tier}. Reason: {reason}"


def test_integration_gnomad_api_failed_cascade_tier2():
    """gnomAD API 失败级联：候选 Tier 1 降级为 Tier 2。"""
    from dgra_core import classify_variant_tier, GPAConfig

    v = make_variant(
        gene="RUNX1", gt="1/1", impact="HIGH",
        gnomad_status="API_FAILED",
        gnomad_error_msg="HTTP 500",
    )
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.api_failed("1", 100, "A", "G")

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 2, f"Expected Tier 2 for API_FAILED cascade, got {tier}. Reason: {reason}"
    assert any("API_FAILED" in str(a) or "Downgraded" in str(a) for a in actions)


def test_integration_full_pipeline_3_variants():
    """3 变异完整 pipeline：过滤 → 分级。"""
    from dgra_core import classify_variant_tier, GPAConfig, Variant
    from dgra_variant_filter import filter_variants

    raw_variants = [
        {"CHROM": "17", "POS": "7579472", "REF": "C", "ALT": "T", "GENE": "TP53",
         "Feature": "ENST00000269305", "EXON": "5/11", "IMPACT": "HIGH",
         "Consequence": "stop_gained", "HGVSp": "p.Arg273Ter", "HGVSc": "c.818C>T",
         "CLIN_SIG": "Pathogenic", "GT": "0/1", "DP": "80", "GQ": "99", "VAF": "0.45",
         "gnomAD_AF": "0.00001"},
        {"CHROM": "X", "POS": "41357831", "REF": "A", "ALT": "T", "GENE": "DDX3X",
         "Feature": "ENST00000373383", "EXON": "5/15", "IMPACT": "HIGH",
         "Consequence": "stop_gained", "HGVSp": "p.Glu233Ter", "HGVSc": "c.697G>A",
         "CLIN_SIG": "", "GT": "0/1", "DP": "60", "GQ": "99", "VAF": "0.45",
         "gnomAD_AF": "0.45"},
        {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "BRCA1",
         "Feature": "ENST00000346300", "EXON": "10/24", "IMPACT": "MODERATE",
         "Consequence": "missense", "HGVSp": "p.Val1Leu", "HGVSc": "c.1A>G",
         "CLIN_SIG": "", "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.48",
         "gnomAD_AF": "0.02"},
    ]

    # 过滤
    filtered, stats = filter_variants(raw_variants, preset="clinical")
    assert stats["output_count"] == 3

    # 分级（TP53用primary，DDX3X/BRCA1用none触发fast-track）
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")

    results = []
    for raw in filtered:
        v = Variant(
            chrom=raw["CHROM"], pos=int(raw["POS"]), ref=raw["REF"], alt=raw["ALT"],
            gene=raw["GENE"], transcript=raw["Feature"], exon=raw["EXON"],
            impact=raw["IMPACT"], consequence=raw["Consequence"],
            hgvsp=raw["HGVSp"], hgvsc=raw["HGVSc"], clinvar=raw["CLIN_SIG"],
            gt=raw["GT"], dp=int(raw["DP"]), gq=float(raw["GQ"]),
            vaf=float(raw["VAF"]),
            gnomad_af=float(raw["gnomAD_AF"]) if raw["gnomAD_AF"] else None,
            gnomad_status="SUCCESS",
        )
        if v.gene == "TP53":
            tissue = MockTissueAssessment.primary()
            gnomad_info = MockGnomAD.success_rare("17", 7579472, "C", "T", af=0.00001)
        else:
            tissue = MockTissueAssessment.none()
            if v.gene == "DDX3X":
                gnomad_info = MockGnomAD.success_common("X", 41357831, "A", "T", af=0.45)
            else:
                gnomad_info = MockGnomAD.success_common("1", 100, "A", "G", af=0.02)

        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, gnomad_info, {}, None, profile, config
        )
        results.append((v.gene, tier, reason))

    # 验证
    tp53 = next(r for r in results if r[0] == "TP53")
    ddx3x = next(r for r in results if r[0] == "DDX3X")
    brca1 = next(r for r in results if r[0] == "BRCA1")

    assert tp53[1] == 1, f"TP53 should be Tier 1, got {tp53[1]}"
    assert ddx3x[1] == 3, f"DDX3X should be Tier 3, got {ddx3x[1]}"
    assert brca1[1] == 3, f"BRCA1 should be Tier 3, got {brca1[1]}"


if __name__ == "__main__":
    print("=" * 60)
    print("L3 Integration Tests")
    print("=" * 60)
    tests = [
        test_integration_raw_vcf_to_tier_ddx3x_tier3,
        test_integration_chinese_csv_tp53_tier1,
        test_integration_chinese_csv_ddx3x_tier3,
        test_integration_gnomad_api_failed_cascade_tier2,
        test_integration_full_pipeline_3_variants,
    ]
    run_tests("L3", tests)
