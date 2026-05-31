"""
L1 Static Tests — Module importability, dataclass defaults, circular dependency checks.

Run: pytest -m l1 tests/l1_static/
"""

import importlib
import sys
from pathlib import Path

import pytest

SCRIPT_DIR = Path(__file__).parent.parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


# =============================================================================
# Module imports
# =============================================================================

CORE_MODULES = [
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
    "gpa_pipeline",
    "gpa_tier_classifier",
    "gpa_report",
    "gpa_phaser",
    "gpa_multi_hit",
    "gpa_qc",
    "gpa_vcf_annotator",
    "gpa_preflight",
    "gpa_two_phase",
    "gpa_workflow",
    "gpa_workflow_runner",
    "gpa_workflow_pm",
    "gpa_review_gate",
    "dgra_cli_wrapper",
    "dgra_gene_sync",
    "convert_csv_to_gpa",
]


@pytest.mark.l1
@pytest.mark.smoke
@pytest.mark.parametrize("module_name", CORE_MODULES)
def test_module_imports(module_name: str):
    """All core modules import successfully without ImportError."""
    mod = importlib.import_module(module_name)
    assert mod is not None


# =============================================================================
# Dataclass defaults
# =============================================================================

@pytest.mark.l1
@pytest.mark.smoke
def test_variant_dataclass_defaults():
    """Variant dataclass fields have correct default values."""
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


@pytest.mark.l1
def test_evidence_dataclass_defaults():
    """Evidence dataclass defaults are correct."""
    from dgra_core import Evidence
    e = Evidence(source="ClinVar", rule="test")
    assert e.weight == 1.0
    assert e.confidence == "high"
    assert e.raw_data is None


@pytest.mark.l1
def test_gpa_config_defaults():
    """GPAConfig default values are correct."""
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


# =============================================================================
# Circular dependency check
# =============================================================================

@pytest.mark.l1
@pytest.mark.smoke
def test_no_circular_imports():
    """Core modules have no circular import errors."""
    for mod_name in ["dgra_core", "dgra_api", "dgra_adapters", "gpa_i18n"]:
        mod = importlib.import_module(mod_name)
        # Save original dict snapshot to restore after reload
        original_dict = dict(mod.__dict__)
        importlib.reload(mod)
        assert mod is not None
        # Restore original dict to prevent class identity issues in downstream tests
        mod.__dict__.clear()
        mod.__dict__.update(original_dict)


# =============================================================================
# Configuration files
# =============================================================================

@pytest.mark.l1
def test_filter_presets_exist():
    """Filter presets exist and are well-formed."""
    from dgra_variant_filter import FILTER_PRESETS
    assert "strict" in FILTER_PRESETS
    assert "clinical" in FILTER_PRESETS
    assert "broad" in FILTER_PRESETS
    for name, preset in FILTER_PRESETS.items():
        assert "impacts" in preset
        assert isinstance(preset["impacts"], set)


@pytest.mark.l1
def test_chinese_column_map_populated():
    """Chinese column map is non-empty and contains key columns."""
    from gpa_i18n import CHINESE_COLUMN_MAP
    assert len(CHINESE_COLUMN_MAP) > 10
    assert "位置" in CHINESE_COLUMN_MAP
    assert "变异后果" in CHINESE_COLUMN_MAP
    assert "影响程度" in CHINESE_COLUMN_MAP
    assert "ClinVar" in CHINESE_COLUMN_MAP


@pytest.mark.l1
def test_tissue_context_json_exists():
    """tissue_context.json exists and is valid JSON."""
    import json
    path = Path(__file__).parent.parent.parent / "references" / "tissue_context.json"
    assert path.exists(), f"{path} not found"
    data = json.loads(path.read_text(encoding="utf-8"))
    assert isinstance(data, dict)
    assert len(data) > 0
