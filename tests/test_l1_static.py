"""
L1 静态测试 — 模块可导入性、数据结构默认值、循环依赖检查
"""

import sys
from pathlib import Path

SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

from conftest import run_tests


def test_import_all_modules():
    """所有核心模块可成功导入，无ImportError。"""
    modules = [
        "dgra_core",
        "dgra_api",
        "dgra_adapters",
        "dgra_variant_filter",
        "dgra_config",
        "dgra_cache",
        "dgra_input_parsers",
        "gpa_i18n",
        "gpa_phenotype_match",
        "gpa_transcript_selector",
        "dgra_splice_predictor",
        "dgra_myvariant",
        "dgra_batch_runner",
        "dgra_pseudogene_sync",
    ]
    for mod in modules:
        __import__(mod)
        print(f"    imported {mod}")


def test_variant_dataclass_defaults():
    """Variant dataclass字段默认值正确。"""
    from dgra_core import Variant
    v = Variant(
        chrom="1", pos=100, ref="A", alt="G",
        gene="TP53", transcript="ENST000001", exon="1/10",
        impact="HIGH", consequence="stop_gained",
        hgvsp="p.Arg1Ter", hgvsc="c.1A>T", clinvar="",
    )
    assert v.gnomad_af is None
    assert v.dp == 0
    assert v.gq == 0.0
    assert v.gt == ""
    assert v.vaf is None
    assert v.tier is None
    assert v.tier_reason == ""
    assert v.tier_actions == []
    assert v.quality_confidence == "high"
    assert v.missing_fields == []
    assert v.evidence_chain == []
    assert v.qc_flags == []
    assert v.upgrade_conditions == []
    assert v.gene_constraint is None
    assert v.tier_confidence == "LOW"
    assert v.gnomad_status == "UNKNOWN"
    assert v.gnomad_af_warning is False
    assert v.gnomad_error_msg is None
    assert v.phenotype_match_score is None
    assert v.transcript_selection_method == "canonical"
    assert v.transcript_ambiguity_flag is False


def test_evidence_dataclass_defaults():
    """Evidence dataclass默认值正确。"""
    from dgra_core import Evidence
    e = Evidence(source="ClinVar", rule="test")
    assert e.weight == 1.0
    assert e.confidence == "high"
    assert e.raw_data is None


def test_gpa_config_defaults():
    """GPAConfig默认值正确。"""
    from dgra_core import GPAConfig
    c = GPAConfig()
    assert c.min_dp == 20
    assert c.min_gq == 90.0
    assert c.common_af_threshold == 0.01
    assert c.low_af_threshold == 0.001
    assert c.vaf_deviation_threshold == 0.20
    assert c.tissue_profile is None
    assert c.offline_mode is False
    assert c.somatic_mode is False
    assert c.spliceai_enabled is False


def test_no_circular_imports():
    """核心模块间无循环依赖导致的运行时错误。"""
    # 已经在 test_import_all_modules 中验证了可导入
    # 额外检查：重新导入不会触发 RecursionError 或 AttributeError
    import importlib
    for mod_name in ["dgra_core", "dgra_api", "dgra_adapters", "gpa_i18n"]:
        mod = importlib.import_module(mod_name)
        importlib.reload(mod)
        assert mod is not None


def test_filter_presets_exist():
    """过滤预设存在且格式正确。"""
    from dgra_variant_filter import FILTER_PRESETS
    assert "strict" in FILTER_PRESETS
    assert "clinical" in FILTER_PRESETS
    assert "broad" in FILTER_PRESETS
    for name, preset in FILTER_PRESETS.items():
        assert "impacts" in preset
        assert isinstance(preset["impacts"], set)


def test_chinese_column_map_populated():
    """中文表头映射表非空且包含关键列。"""
    from gpa_i18n import CHINESE_COLUMN_MAP
    assert len(CHINESE_COLUMN_MAP) > 10
    assert "位置" in CHINESE_COLUMN_MAP
    assert "变异后果" in CHINESE_COLUMN_MAP
    assert "影响程度" in CHINESE_COLUMN_MAP
    assert "ClinVar" in CHINESE_COLUMN_MAP


if __name__ == "__main__":
    print("=" * 60)
    print("L1 Static Tests")
    print("=" * 60)
    tests = [
        test_import_all_modules,
        test_variant_dataclass_defaults,
        test_evidence_dataclass_defaults,
        test_gpa_config_defaults,
        test_no_circular_imports,
        test_filter_presets_exist,
        test_chinese_column_map_populated,
    ]
    run_tests("L1", tests)
