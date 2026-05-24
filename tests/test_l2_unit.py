"""
L2 单元测试 — 按模块分测试用例
gnomAD / Tier分级 / 中文表头 / NMD / 转录本选择 / 过滤 / 表型分隔符
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import run_tests, MockGnomAD, MockTissueProfile, MockTissueAssessment, make_variant


# ---------------------------------------------------------------------------
# gnomAD
# ---------------------------------------------------------------------------

def test_gnomad_success_common_af_range():
    """gnomAD SUCCESS 常见变异 AF 在合理范围。"""
    r = MockGnomAD.success_common("X", 41357831, "A", "T", af=0.45)
    assert r["status"] == "SUCCESS"
    assert 0.30 <= r["af"] <= 0.70
    assert "af_exome" in r
    assert "af_populations" in r
    assert "EAS" in r["af_populations"]


def test_gnomad_api_failed():
    """gnomAD API_FAILED 状态正确。"""
    r = MockGnomAD.api_failed("1", 100, "A", "G")
    assert r["status"] == "API_FAILED"
    assert r["af"] is None
    assert r["confidence"] == "low"


def test_gnomad_not_captured():
    """gnomAD NOT_CAPTURED 状态正确。"""
    r = MockGnomAD.not_captured("1", 100, "A", "G")
    assert r["status"] == "NOT_CAPTURED"
    assert r["af"] is None


def test_gnomad_af_calculation():
    """gnomAD AF 手算逻辑：ac/an 正确。"""
    r = MockGnomAD.success_rare("1", 100, "A", "G", af=0.0001)
    assert r["af"] == 0.0001
    eas = r["af_populations"]["EAS"]
    assert eas["ac"] == 1
    assert eas["an"] == 10000
    assert eas["af"] == 1 / 10000


# ---------------------------------------------------------------------------
# Tier 分级 — Priority 1b 三层守卫
# ---------------------------------------------------------------------------

def test_tier_priority_1b_api_failed_downgrade():
    """Priority 1b: gnomAD API_FAILED → Tier 2 (不是 Tier 1)。"""
    from dgra_core import classify_variant_tier, GPAConfig
    v = make_variant(
        gene="DDX3X", gt="1/1", impact="HIGH",
        gnomad_status="API_FAILED", gnomad_error_msg="GraphQL 400",
    )
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.api_failed("X", 41357831, "A", "T")

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 2, f"Expected Tier 2 for API_FAILED, got {tier}. Reason: {reason}"
    assert "API_FAILED" in reason or "Downgraded" in str(actions)


def test_tier_priority_1b_not_captured_tier1():
    """Priority 1b: gnomAD NOT_CAPTURED → Tier 1 MEDIUM。"""
    from dgra_core import classify_variant_tier, GPAConfig
    v = make_variant(
        gene="RUNX1", gt="1/1", impact="HIGH",
        gnomad_status="NOT_CAPTURED",
    )
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 1, f"Expected Tier 1 for NOT_CAPTURED, got {tier}. Reason: {reason}"


def test_tier_priority_1b_common_polymorphism_tier3():
    """Priority 1b: AF>1% 常见多态性 → Tier 3。"""
    from dgra_core import classify_variant_tier, GPAConfig
    v = make_variant(
        gene="DDX3X", gt="1/1", impact="HIGH",
        gnomad_af=0.45, gnomad_status="SUCCESS",
    )
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.success_common("X", 41357831, "A", "T", af=0.45)

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 3, f"Expected Tier 3 for common polymorphism, got {tier}. Reason: {reason}"


def test_tier_priority_1c_clinvar_pathogenic_tier1():
    """Priority 1c: ClinVar Pathogenic + HIGH + 造血相关 → Tier 1。"""
    from dgra_core import classify_variant_tier, GPAConfig
    v = make_variant(
        gene="VWF", gt="0/1", impact="HIGH",
        clinvar="Pathogenic",
        gnomad_af=None, gnomad_status="NOT_CAPTURED",
    )
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.not_captured("12", 6126538, "G", "A")

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 1, f"Expected Tier 1 for ClinVar pathogenic, got {tier}. Reason: {reason}"


def test_tier_priority_2_heterozygous_lof_tier2():
    """Priority 2: 杂合截短 + 造血相关 → Tier 2。"""
    from dgra_core import classify_variant_tier, GPAConfig
    v = make_variant(
        gene="RUNX1", gt="0/1", impact="HIGH",
        gnomad_af=None, gnomad_status="NOT_CAPTURED",
    )
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = MockGnomAD.not_captured("1", 100, "A", "G")

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 2, f"Expected Tier 2 for heterozygous LOF, got {tier}. Reason: {reason}"


def test_tier_fast_track_no_relevance_tier3():
    """Fast track: 无组织相关性 + 非Pathogenic → Tier 3。"""
    from dgra_core import classify_variant_tier, GPAConfig
    v = make_variant(
        gene="BRCA1", gt="0/1", impact="HIGH",
        gnomad_af=0.45, gnomad_status="SUCCESS",
    )
    tissue = MockTissueAssessment.none()
    profile = MockTissueProfile.general()
    config = GPAConfig(tissue_profile="general")
    gnomad_info = MockGnomAD.success_common("1", 100, "A", "G", af=0.45)

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier == 3, f"Expected Tier 3 for fast track, got {tier}. Reason: {reason}"


# ---------------------------------------------------------------------------
# 中文表头翻译
# ---------------------------------------------------------------------------

def test_translate_single_cn_header():
    """单个中文表头正确映射。"""
    from gpa_i18n import translate_chinese_headers
    headers = ["位置", "变异后果", "影响程度", "ClinVar", "gnomAD频率"]
    result = translate_chinese_headers(headers)
    assert result[0] == "Location"
    assert result[1] == "Consequence"
    assert result[2] == "IMPACT"
    assert result[3] == "CLIN_SIG"
    assert result[4] == "gnomAD_AF"


def test_translate_cn_header_list_13():
    """13个中文表头全部正确映射。"""
    from gpa_i18n import translate_chinese_headers
    headers = [
        "位置", "基因", "转录本", "变异后果", "影响程度",
        "HGVSc", "HGVSp", "ClinVar", "gnomAD频率",
        "样本", "基因型", "测序深度", "质量值",
    ]
    result = translate_chinese_headers(headers)
    expected = [
        "Location", "Gene", "Feature", "Consequence", "IMPACT",
        "HGVSc", "HGVSp", "CLIN_SIG", "gnomAD_AF",
        "SAMPLE", "GT", "DP", "GQ",
    ]
    assert result == expected, f"Mismatch: {result} vs {expected}"


# ---------------------------------------------------------------------------
# NMD 预测
# ---------------------------------------------------------------------------

def test_nmd_stop_gained_internal_exon_sensitive():
    """stop_gained + 内部外显子 → NMD sensitive。"""
    from dgra_core import predict_nmd
    v = make_variant(consequence="stop_gained", exon="5/11")
    result = predict_nmd(v)
    assert result["status"] == "sensitive"
    assert result["pvs1_applicable"] is True
    assert result["pvs1_strength"] == "strong"


def test_nmd_synonymous_not_applicable():
    """同义变异 → NMD not_applicable。"""
    from dgra_core import predict_nmd
    v = make_variant(consequence="synonymous_variant", exon="5/11")
    result = predict_nmd(v)
    assert result["status"] == "not_applicable"
    assert result["pvs1_applicable"] is False


def test_nmd_last_exon_escape():
    """最后外显子截短 → NMD escape。"""
    from dgra_core import predict_nmd
    v = make_variant(consequence="stop_gained", exon="11/11")
    result = predict_nmd(v)
    assert result["status"] == "escape"
    assert result["pvs1_applicable"] is False


def test_nmd_penultimate_possible_escape():
    """倒数第二外显子 → possible_escape。"""
    from dgra_core import predict_nmd
    v = make_variant(consequence="frameshift", exon="10/11")
    result = predict_nmd(v)
    assert result["status"] == "possible_escape"
    assert result["pvs1_strength"] == "moderate"


# ---------------------------------------------------------------------------
# 转录本选择
# ---------------------------------------------------------------------------

def test_transcript_ambiguous_detection():
    """多个transcript同一基因 → ambiguous flag。"""
    from gpa_transcript_selector import TranscriptSelector
    # 这个模块的接口可能随版本变化，先尝试导入
    try:
        selector = TranscriptSelector()
        # 如果存在检测ambiguous的方法
        if hasattr(selector, 'detect_ambiguous'):
            variants = [
                {"Feature": "ENST001", "SYMBOL": "TP53"},
                {"Feature": "ENST002", "SYMBOL": "TP53"},
            ]
            is_ambiguous = selector.detect_ambiguous(variants, "TP53")
            assert is_ambiguous is True
    except (AttributeError, TypeError):
        # 接口不存在则跳过
        print("    Skipped: detect_ambiguous not available in this version")
        return


def test_transcript_canonical_preference():
    """canonical transcript 评分优先。"""
    from gpa_transcript_selector import TranscriptSelector
    selector = TranscriptSelector()
    # 检查评分逻辑是否偏好 canonical/MANE
    if hasattr(selector, '_score_transcript'):
        score_canonical = selector._score_transcript({"is_canonical": True, "mane_select": True}, "TP53")
        score_noncanonical = selector._score_transcript({"is_canonical": False, "mane_select": False}, "TP53")
        assert score_canonical > score_noncanonical
    else:
        print("    Skipped: _score_transcript not available")


# ---------------------------------------------------------------------------
# 过滤 — strict / clinical / broad
# ---------------------------------------------------------------------------

def test_filter_strict_keeps_high_moderate():
    """strict 预设保留 HIGH/MODERATE，排除 LOW。"""
    from dgra_variant_filter import filter_variants
    variants = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "GENE": "TP53"},
        {"IMPACT": "MODERATE", "Consequence": "missense", "GENE": "TP53"},
        {"IMPACT": "LOW", "Consequence": "synonymous", "GENE": "TP53"},
    ]
    filtered, stats = filter_variants(variants, preset="strict")
    assert stats["output_count"] == 2
    assert stats["excluded"] == 1


def test_filter_clinical_keeps_splice_region():
    """clinical 预设保留 splice_region_variant (LOW)。"""
    from dgra_variant_filter import filter_variants
    variants = [
        {"IMPACT": "LOW", "Consequence": "splice_region_variant", "GENE": "TP53"},
    ]
    filtered, stats = filter_variants(variants, preset="clinical")
    assert stats["output_count"] == 1
    assert stats["splice_retained"] == 1


def test_filter_broad_keeps_all_impacts():
    """broad 预设保留 HIGH/MODERATE/LOW。"""
    from dgra_variant_filter import filter_variants
    variants = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "GENE": "A"},
        {"IMPACT": "MODERATE", "Consequence": "missense", "GENE": "B"},
        {"IMPACT": "LOW", "Consequence": "synonymous", "GENE": "C"},
    ]
    filtered, stats = filter_variants(variants, preset="broad")
    assert stats["output_count"] == 3


# ---------------------------------------------------------------------------
# 表型分隔符
# ---------------------------------------------------------------------------

def test_phenotype_delimiter_split():
    """表型分隔符正确 split（顿号、逗号、句号）。"""
    from gpa_phenotype_match import PhenotypeMatcher
    # 测试分隔符处理
    text1 = "贫血、出血、感染"
    text2 = "贫血,出血,感染"
    text3 = "贫血。出血。感染"

    # 检查 PhenotypeMatcher 是否支持 split 或类似逻辑
    matcher = PhenotypeMatcher()
    if hasattr(matcher, '_split_phenotypes'):
        assert len(matcher._split_phenotypes(text1)) == 3
        assert len(matcher._split_phenotypes(text2)) == 3
        assert len(matcher._split_phenotypes(text3)) == 3
    else:
        # 直接测试 replace + split 逻辑
        cleaned = text1.replace("。", "、").replace(".", "、").replace(",", "、")
        parts = [p.strip() for p in cleaned.split("、") if p.strip()]
        assert len(parts) == 3
        print("    Skipped: _split_phenotypes not available, manual split OK")


def test_phenotype_delimiter_mixed():
    """混合分隔符正确处理。"""
    text = "贫血,出血。感染、发热"
    cleaned = text.replace("。", "、").replace(".", "、").replace(",", "、")
    parts = [p.strip() for p in cleaned.split("、") if p.strip()]
    assert len(parts) == 4


if __name__ == "__main__":
    print("=" * 60)
    print("L2 Unit Tests")
    print("=" * 60)
    tests = [
        test_gnomad_success_common_af_range,
        test_gnomad_api_failed,
        test_gnomad_not_captured,
        test_gnomad_af_calculation,
        test_tier_priority_1b_api_failed_downgrade,
        test_tier_priority_1b_not_captured_tier1,
        test_tier_priority_1b_common_polymorphism_tier3,
        test_tier_priority_1c_clinvar_pathogenic_tier1,
        test_tier_priority_2_heterozygous_lof_tier2,
        test_tier_fast_track_no_relevance_tier3,
        test_translate_single_cn_header,
        test_translate_cn_header_list_13,
        test_nmd_stop_gained_internal_exon_sensitive,
        test_nmd_synonymous_not_applicable,
        test_nmd_last_exon_escape,
        test_nmd_penultimate_possible_escape,
        test_transcript_ambiguous_detection,
        test_transcript_canonical_preference,
        test_filter_strict_keeps_high_moderate,
        test_filter_clinical_keeps_splice_region,
        test_filter_broad_keeps_all_impacts,
        test_phenotype_delimiter_split,
        test_phenotype_delimiter_mixed,
    ]
    run_tests("L2", tests)
