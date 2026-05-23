"""GPA v0.9.1 Hotfix Tests — DDX3X misclassification bug fixes"""
import unittest, sys, os, asyncio
from pathlib import Path

sys.path.insert(0, os.path.dirname(__file__) + "/../scripts")

from dgra_core import Variant, classify_variant_tier, GPAConfig
from gpa_i18n import is_chinese_header, translate_chinese_header

# _impact_high is nested inside classify_variant_tier; define locally for tests
def _impact_high(impact):
    if impact is None or str(impact).strip().upper() in ('UNKNOWN', 'NA', 'N/A', '', 'NONE'):
        return True
    return str(impact).strip().upper() == 'HIGH'


def make_variant(**kwargs):
    """Helper to create a Variant with all required fields filled."""
    defaults = dict(
        chrom="1", pos=12345, ref="G", alt="A",
        gene="TEST", transcript="ENST000001", exon="E1/5",
        impact="HIGH", consequence="missense_variant",
        hgvsp="p.Ala1Val", hgvsc="c.2C>T", clinvar="VUS",
    )
    defaults.update(kwargs)
    return Variant(**defaults)


class TestGnomADStatus(unittest.TestCase):
    """Phase 1+2: gnomAD API status marking"""

    def test_api_failed_status(self):
        v = make_variant(
            chrom="X", pos=41357831, ref="A", alt="T",
            gene="DDX3X", consequence="splice_acceptor_variant",
            hgvsp="", hgvsc="c.3131-2A>T", clinvar="Benign",
            gt="1/1", gnomad_status="API_FAILED",
            gnomad_error_msg="timeout",
            gnomad_af_warning=True,
        )
        self.assertEqual(v.gnomad_status, "API_FAILED")
        self.assertTrue(v.gnomad_af_warning)
        self.assertEqual(v.gnomad_error_msg, "timeout")

    def test_not_captured_status(self):
        v = make_variant(
            chrom="1", pos=12345, ref="G", alt="A",
            gene="TEST", gt="1/1", gnomad_status="NOT_CAPTURED",
        )
        self.assertEqual(v.gnomad_status, "NOT_CAPTURED")
        self.assertIsNone(v.gnomad_error_msg)

    def test_success_status(self):
        v = make_variant(
            chrom="1", pos=12345, ref="G", alt="A",
            gene="TEST", gt="0/1",
            gnomad_status="SUCCESS", gnomad_af=0.0001,
        )
        self.assertEqual(v.gnomad_status, "SUCCESS")
        self.assertEqual(v.gnomad_af, 0.0001)


class TestTier1Guard(unittest.TestCase):
    """Phase 3: Priority 1b gnomAD guard logic"""

    def test_api_failed_downgrades_to_tier2(self):
        """DDX3X scenario: API_FAILED + homozygous HIGH → Tier 2"""
        v = make_variant(
            chrom="X", pos=41357831, ref="A", alt="T",
            gene="DDX3X", consequence="splice_acceptor_variant",
            hgvsp="", hgvsc="c.3131-2A>T", clinvar="Benign",
            gt="1/1", gnomad_status="API_FAILED",
            gnomad_error_msg="gnomAD API timeout after 3 retries",
            gnomad_af_warning=True,
        )
        tissue = {"relevance": "primary"}
        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, {}, None, None, tissue, GPAConfig()
        )
        self.assertEqual(tier, 2, f"Expected Tier 2, got {tier}. Reason: {reason}")
        self.assertIn("gnomAD query FAILED", reason)
        self.assertTrue(v.gnomad_af_warning)
        self.assertTrue(any("Downgraded to Tier 2" in a for a in actions))

    def test_not_captured_keeps_tier1_medium_conf(self):
        """NOT_CAPTURED + homozygous HIGH → Tier 1, confidence=MEDIUM"""
        v = make_variant(
            chrom="1", pos=99999, ref="C", alt="T",
            gene="RARE", consequence="frameshift_variant",
            hgvsp="p.Arg100Ter", hgvsc="c.298C>T", clinvar="VUS",
            gt="1/1", gnomad_status="NOT_CAPTURED",
        )
        tissue = {"relevance": "primary"}
        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, {}, None, None, tissue, GPAConfig()
        )
        self.assertEqual(tier, 1)
        self.assertIn("gnomAD not captured", reason)
        # v0.9.1: NOT_CAPTURED keeps Tier 1 but with note about missing frequency data

    def test_common_polymorphism_excluded(self):
        """SUCCESS + AF>1% → Tier 3"""
        v = make_variant(
            chrom="X", pos=41357831, ref="A", alt="T",
            gene="DDX3X", consequence="splice_acceptor_variant",
            hgvsp="", hgvsc="c.3131-2A>T", clinvar="Benign",
            gt="1/1", gnomad_af=0.60,
            gnomad_status="SUCCESS",
        )
        tissue = {"relevance": "primary"}
        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, {}, None, None, tissue, GPAConfig()
        )
        self.assertEqual(tier, 3, f"Expected Tier 3 for AF=60%, got {tier}")
        self.assertIn("common polymorphism", reason.lower())

    def test_rare_af_tier1_high_conf(self):
        """SUCCESS + AF<1% → Tier 1 HIGH confidence"""
        v = make_variant(
            chrom="1", pos=12345, ref="G", alt="A",
            gene="VWF", consequence="stop_gained",
            hgvsp="p.Gln1311Ter", hgvsc="c.3931C>T", clinvar="Pathogenic",
            gt="1/1", gnomad_af=0.00001,
            gnomad_status="SUCCESS",
        )
        tissue = {"relevance": "primary"}
        tier, reason, actions = classify_variant_tier(
            v, {}, tissue, {}, None, None, tissue, GPAConfig()
        )
        self.assertEqual(tier, 1)


class TestChineseHeaderTranslation(unittest.TestCase):
    """Phase 4: Chinese VEP CSV header translation"""

    def test_detect_chinese_header(self):
        headers = ["位置", "基因", "变异后果", "影响程度", "gnomAD频率", "ClinVar"]
        self.assertTrue(is_chinese_header(headers))

    def test_detect_english_header(self):
        headers = ["Location", "Gene", "Consequence", "IMPACT", "gnomAD_AF", "CLIN_SIG"]
        self.assertFalse(is_chinese_header(headers))

    def test_translate_full_chinese_header(self):
        headers = [
            "Uploaded_variation", "位置", "基因", "转录本", "变异后果",
            "影响程度", "CDNA位置", "CDS位置", "蛋白位置", "氨基酸改变",
            "密码子", "HGVC", "HGVSp", "现有等位基因", "rs号",
            "参考等位基因", "替代等位基因", "gnomAD频率", "ClinVar",
            "样本", "基因型", "测序深度", "质量值", "距离", "链", "突变频谱",
        ]
        translated = translate_chinese_header(headers)
        expected = {
            "位置": "Location",
            "基因": "Gene",
            "转录本": "Feature",
            "变异后果": "Consequence",
            "影响程度": "IMPACT",
            "gnomAD频率": "gnomAD_AF",
            "ClinVar": "CLIN_SIG",
            "基因型": "GT",
            "测序深度": "DP",
        }
        for cn, en in expected.items():
            idx = headers.index(cn)
            self.assertEqual(translated[idx], en, f"Header '{cn}' should map to '{en}', got '{translated[idx]}'")

    def test_preserve_unmapped_headers(self):
        headers = ["未知列", "位置", "AnotherUnknown"]
        translated = translate_chinese_header(headers)
        self.assertEqual(translated[0], "未知列")
        self.assertEqual(translated[1], "Location")
        self.assertEqual(translated[2], "AnotherUnknown")


if __name__ == "__main__":
    unittest.main(verbosity=2)
