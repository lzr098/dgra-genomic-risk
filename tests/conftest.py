"""
GPA Test Fixtures — Mock utilities for L1-L5 test suite.
No external API calls; pure Python mocks.
"""

from dataclasses import dataclass, field
from typing import Dict, Any, List, Optional
import sys
from pathlib import Path

# Ensure scripts dir is on path
SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))


class MockGnomAD:
    """Mock gnomAD API responses."""

    @staticmethod
    def success_common(chrom: str, pos: int, ref: str, alt: str, af: float = 0.45):
        """Return a SUCCESS response for a common variant (DDX3X-like)."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": af,
            "af_popmax": af,
            "af_populations": {
                "EAS": {"af": af, "ac": 100, "an": 222},
                "NFE": {"af": 0.35, "ac": 80, "an": 228},
            },
            "af_exome": af,
            "af_genome": None,
            "an_exome": 222,
            "an_genome": None,
            "hom_count": 12,
            "status": "SUCCESS",
            "source": "gnomad",
            "confidence": "medium",
            "raw": {},
        }

    @staticmethod
    def success_rare(chrom: str, pos: int, ref: str, alt: str, af: float = 0.0001):
        """Return a SUCCESS response for a rare variant."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": af,
            "af_popmax": af,
            "af_populations": {
                "EAS": {"af": af, "ac": 1, "an": 10000},
            },
            "af_exome": af,
            "af_genome": None,
            "an_exome": 10000,
            "an_genome": None,
            "hom_count": 0,
            "status": "SUCCESS",
            "source": "gnomad",
            "confidence": "medium",
            "raw": {},
        }

    @staticmethod
    def api_failed(chrom: str, pos: int, ref: str, alt: str):
        """Return an API_FAILED response."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": None,
            "af_popmax": None,
            "af_populations": {},
            "status": "API_FAILED",
            "source": "failed",
            "confidence": "low",
            "error": "GraphQL 400: unknown error",
        }

    @staticmethod
    def not_captured(chrom: str, pos: int, ref: str, alt: str):
        """Return a NOT_CAPTURED response."""
        return {
            "variant_id": f"{chrom}-{pos}-{ref}-{alt}",
            "af": None,
            "af_popmax": None,
            "af_populations": {},
            "status": "NOT_CAPTURED",
            "source": "gnomad",
            "confidence": "medium",
            "note": "Variant not in gnomAD dataset",
            "raw": {},
        }


class MockEnsemblVEP:
    """Mock Ensembl VEP responses."""

    @staticmethod
    def vep_result(variant_id: str, gene: str, consequence: str, impact: str,
                   hgvsp: str = "", hgvsc: str = "", exon: str = ""):
        """Return a minimal VEP result dict."""
        return {
            "input": variant_id,
            "transcript_consequences": [{
                "gene_symbol": gene,
                "consequence_terms": [consequence],
                "impact": impact,
                "hgvsp": hgvsp,
                "hgvsc": hgvsc,
                "exon": exon,
            }]
        }

    @staticmethod
    def batch_results(variants: List[Dict]) -> List[Dict]:
        """Return mocked VEP results for a batch of variants."""
        return [MockEnsemblVEP.vep_result(**v) for v in variants]


class MockTissueProfile:
    """Mock tissue profile for testing tier classification."""

    @staticmethod
    def hematopoietic():
        return {
            "display_name": "hematopoietic",
            "tier_rules": {},
            "special_gene_lists": {
                "coagulation": {"VWF", "F8", "F9"},
                "fa_dna_repair": {"BRCA1", "BRCA2", "FANCA"},
                "drug_metabolism": {"CYP2D6", "TPMT", "DPYD"},
            },
            "tissue_genes": {
                "RUNX1", "CEBPA", "GATA2", "ASXL1", "BCOR", "BCORL1", "PHF6",
                "FLT3", "NPM1", "IDH1", "IDH2", "DNMT3A", "TET2",
                "TP53", "KIT", "NRAS", "KRAS", "PTPN11",
            },
        }

    @staticmethod
    def general():
        return {
            "display_name": "general",
            "tier_rules": {},
            "special_gene_lists": {
                "coagulation": {"VWF", "F8", "F9"},
                "fa_dna_repair": {"BRCA1", "BRCA2"},
            },
            "tissue_genes": set(),
        }


class MockTissueAssessment:
    """Build tissue_assessment dicts for classify_variant_tier."""

    @staticmethod
    def primary(gtex_tpm: float = 50.0):
        return {
            "relevance": "primary",
            "gtex_tpm": gtex_tpm,
            "fast_track": False,
            "tier_suggestion": None,
            "reason": "Primary tissue gene",
        }

    @staticmethod
    def secondary(gtex_tpm: float = 10.0):
        return {
            "relevance": "secondary",
            "gtex_tpm": gtex_tpm,
            "fast_track": False,
            "tier_suggestion": None,
            "reason": "Secondary tissue gene",
        }

    @staticmethod
    def none(gtex_tpm: float = 0.0):
        return {
            "relevance": "none",
            "gtex_tpm": gtex_tpm,
            "fast_track": True,
            "tier_suggestion": 3,
            "reason": "No tissue relevance",
        }


def make_variant(**kwargs) -> Any:
    """Factory: create a dgra_core.Variant with sensible defaults."""
    from dgra_core import Variant
    defaults = {
        "chrom": "1",
        "pos": 100000,
        "ref": "A",
        "alt": "G",
        "gene": "TP53",
        "transcript": "ENST00000269305",
        "exon": "5/11",
        "impact": "HIGH",
        "consequence": "stop_gained",
        "hgvsp": "p.Arg273Ter",
        "hgvsc": "c.818C>T",
        "clinvar": "",
        "gnomad_af": None,
        "dp": 50,
        "gq": 99.0,
        "gt": "0/1",
        "vaf": 0.45,
    }
    defaults.update(kwargs)
    return Variant(**defaults)


def run_tests(module_name: str, test_funcs: List):
    """Simple test runner — no pytest dependency."""
    passed = 0
    failed = 0
    errors = []
    for fn in test_funcs:
        name = getattr(fn, "__name__", str(fn))
        try:
            # Detect async
            import asyncio
            if asyncio.iscoroutinefunction(fn):
                asyncio.run(fn())
            else:
                fn()
            passed += 1
            print(f"  ✅ {name}")
        except Exception as e:
            failed += 1
            errors.append((name, str(e)))
            print(f"  ❌ {name}: {e}")
    print(f"\n{module_name}: {passed} passed, {failed} failed")
    if errors:
        for name, err in errors:
            print(f"    {name}: {err}")
    return passed, failed
