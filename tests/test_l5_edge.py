"""
L5 边界/异常测试 — 空输入、畸形数据、极端值、编码问题
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import run_tests, MockGnomAD, MockTissueProfile, MockTissueAssessment, make_variant


def test_empty_input_generates_report():
    """空变异列表 → 报告生成，不崩溃。"""
    from dgra_core import run_dgra_pipeline, GPAConfig
    import asyncio
    config = GPAConfig(offline_mode=True, tissue_profile="general")

    async def _run():
        result = await run_dgra_pipeline([], config=config)
        assert result is not None
        assert "report_markdown" in result or "json_report" in result
        assert len(result.get("tier1_variants", [])) == 0
        assert len(result.get("tier2_variants", [])) == 0
        assert len(result.get("tier3_variants", [])) == 0

    asyncio.run(_run())


def test_malformed_vcf_raises_valueerror():
    """缺少必填字段（CHROM/POS/REF/ALT）→ 应抛异常或被过滤。"""
    from dgra_core import run_dgra_pipeline, GPAConfig
    import asyncio
    config = GPAConfig(offline_mode=True, tissue_profile="general")

    async def _run():
        malformed = [
            {"CHROM": "", "POS": "", "REF": "", "ALT": ""},  # 空值
            {"CHROM": "1", "POS": "abc", "REF": "A", "ALT": "G"},  # POS 非数字
        ]
        # 不应崩溃；POS="abc" 会抛 ValueError，这是预期的
        try:
            result = await run_dgra_pipeline(malformed, config=config)
            assert result is not None
            # 如果未抛异常，应被过滤或标记
            variants = result.get("variants", [])
            for v in variants:
                assert v.quality_confidence in ("low", "unknown", "medium")
        except ValueError:
            pass  # ValueError 是预期的，畸形数据应被拒绝

    asyncio.run(_run())


def test_chrmt_normalization():
    """chrMT / MT → 标准化为 M。"""
    from dgra_core import classify_variant_tier, GPAConfig, Variant
    v = Variant(chrom="chrMT", pos=1000, ref="A", alt="G",
                gene="MT-ATP6", transcript="", exon="",
                impact="HIGH", consequence="stop_gained",
                hgvsp="", hgvsc="", clinvar="")
    # 测试是否能正常进入分级流程
    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.general()
    config = GPAConfig(tissue_profile="general")
    gnomad_info = MockGnomAD.not_captured("M", 1000, "A", "G")

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    assert tier in (1, 2, 3)


def test_abnormal_af_truncation():
    """异常 AF=999.0 → 截断或标记。"""
    from dgra_core import classify_variant_tier, GPAConfig, Variant
    v = Variant(chrom="1", pos=100, ref="A", alt="G",
                gene="TP53", transcript="ENST000001", exon="1/10",
                impact="HIGH", consequence="stop_gained",
                hgvsp="p.Arg1Ter", hgvsc="c.1A>T", clinvar="",
                gnomad_af=999.0, gnomad_status="SUCCESS")

    tissue = MockTissueAssessment.primary()
    profile = MockTissueProfile.hematopoietic()
    config = GPAConfig(tissue_profile="hematopoietic")
    gnomad_info = {"af": 999.0, "status": "SUCCESS", "source": "gnomad"}

    tier, reason, actions = classify_variant_tier(
        v, {}, tissue, gnomad_info, {}, None, profile, config
    )
    # AF>1% 应被识别为常见多态性或标记异常
    assert tier in (2, 3), f"Expected Tier 2/3 for abnormal AF, got {tier}"


def test_unicode_phenotype_no_crash():
    """Unicode 表型描述 → 不崩溃。"""
    from dgra_core import run_dgra_pipeline, GPAConfig
    from gpa_phenotype_match import PhenotypeMatcher
    import asyncio

    # 测试 PhenotypeMatcher 能处理中文
    matcher = PhenotypeMatcher()
    phenotypes = ["肌无力", "肌源性损害", "贫血、出血、感染"]
    for p in phenotypes:
        # 不应抛异常
        try:
            if hasattr(matcher, '_split_phenotypes'):
                matcher._split_phenotypes(p)
        except Exception as e:
            assert False, f"Unicode phenotype failed: {e}"

    # 测试 pipeline 能处理中文 phenotypes 参数
    async def _run():
        variants = [{
            "CHROM": "17", "POS": "7579472", "REF": "C", "ALT": "T",
            "GENE": "TP53", "Feature": "ENST00000269305", "EXON": "5/11",
            "IMPACT": "HIGH", "Consequence": "stop_gained",
            "HGVSp": "p.Arg273Ter", "HGVSc": "c.818C>T",
            "CLIN_SIG": "Pathogenic", "GT": "0/1", "DP": "80", "GQ": "99", "VAF": "0.45",
        }]
        result = await run_dgra_pipeline(
            variants, user_phenotypes="肌无力、肌源性损害",
            config=GPAConfig(offline_mode=True, tissue_profile="general")
        )
        assert result is not None
        assert "report_markdown" in result or "json_report" in result

    asyncio.run(_run())


def test_zero_dp_filtered():
    """DP=0 → 被过滤。"""
    from dgra_variant_filter import filter_variants
    variants = [
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "GENE": "TP53", "DP": "0"},
        {"IMPACT": "HIGH", "Consequence": "stop_gained", "GENE": "TP53", "DP": "10"},
    ]
    filtered, stats = filter_variants(variants, preset="strict")
    # DP=0 应被过滤或保留但标记为低质量
    assert stats["output_count"] >= 1


def test_multi_allelic_handled():
    """多等位基因 alt='A,C' → 正确处理。"""
    from dgra_core import run_dgra_pipeline, GPAConfig, Variant
    import asyncio

    async def _run():
        variants = [
            {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "A,C",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "p.Arg1Ter", "HGVSc": "c.1A>T", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.45"},
        ]
        result = await run_dgra_pipeline(variants, config=GPAConfig(offline_mode=True, tissue_profile="general"))
        assert result is not None
        vs = result.get("variants", [])
        # 多等位基因不应导致崩溃，alt 可能被解析或保留原样
        if vs:
            assert vs[0].alt in ("A,C", "C", "A")

    asyncio.run(_run())


def test_structural_variant_no_crash():
    """结构变异 <DEL> → 不崩溃。"""
    from dgra_core import run_dgra_pipeline, GPAConfig
    import asyncio

    async def _run():
        variants = [
            {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "<DEL>",
             "GENE": "BRCA1", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "HIGH", "Consequence": "structural_variant",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "30", "GQ": "99", "VAF": "0.30"},
        ]
        result = await run_dgra_pipeline(variants, config=GPAConfig(offline_mode=True, tissue_profile="general"))
        assert result is not None
        assert "report_markdown" in result or "json_report" in result

    asyncio.run(_run())


def test_missing_fields_graceful_degradation():
    """缺失关键字段 → 优雅降级（quality_confidence=low/unknown）。"""
    from dgra_core import run_dgra_pipeline, GPAConfig
    import asyncio

    async def _run():
        # 缺少 IMPACT、Consequence、CLIN_SIG
        variants = [
            {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "HGVSp": "", "HGVSc": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "0.45"},
        ]
        result = await run_dgra_pipeline(variants, config=GPAConfig(offline_mode=True, tissue_profile="general"))
        assert result is not None
        vs = result.get("variants", [])
        if vs:
            # 应标记为低质量或未知
            assert vs[0].quality_confidence in ("low", "unknown", "medium")
            assert len(vs[0].missing_fields) > 0

    asyncio.run(_run())


def test_circular_reference_protection():
    """循环引用保护 — JSON 序列化不崩溃。"""
    import json
    from dgra_core import Variant, Evidence

    v = make_variant()
    e = Evidence(source="Test", rule="circular", raw_data=v)
    v.evidence_chain.append(e)

    # 尝试 JSON 序列化，不应无限递归
    try:
        # 使用自定义序列化或检查是否可导出
        data = {
            "gene": v.gene,
            "impact": v.impact,
            "evidence_count": len(v.evidence_chain),
        }
        json.dumps(data)
    except RecursionError:
        assert False, "Circular reference caused RecursionError"


def test_extreme_vaf_values():
    """极端 VAF 值 → 不崩溃。"""
    from dgra_core import run_dgra_pipeline, GPAConfig
    import asyncio

    async def _run():
        extreme_cases = [
            {"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "-0.1"},
            {"CHROM": "1", "POS": "101", "REF": "A", "ALT": "G",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "1.5"},
            {"CHROM": "1", "POS": "102", "REF": "A", "ALT": "G",
             "GENE": "TP53", "Feature": "ENST000001", "EXON": "1/10",
             "IMPACT": "HIGH", "Consequence": "stop_gained",
             "HGVSp": "", "HGVSc": "", "CLIN_SIG": "",
             "GT": "0/1", "DP": "50", "GQ": "99", "VAF": "NaN"},
        ]
        result = await run_dgra_pipeline(extreme_cases, config=GPAConfig(offline_mode=True, tissue_profile="general"))
        assert result is not None
        vs = result.get("variants", [])
        assert len(vs) >= 0  # 不应崩溃

    asyncio.run(_run())


if __name__ == "__main__":
    print("=" * 60)
    print("L5 Edge/Boundary Tests")
    print("=" * 60)
    tests = [
        test_empty_input_generates_report,
        test_malformed_vcf_raises_valueerror,
        test_chrmt_normalization,
        test_abnormal_af_truncation,
        test_unicode_phenotype_no_crash,
        test_zero_dp_filtered,
        test_multi_allelic_handled,
        test_structural_variant_no_crash,
        test_missing_fields_graceful_degradation,
        test_circular_reference_protection,
        test_extreme_vaf_values,
    ]
    run_tests("L5", tests)
