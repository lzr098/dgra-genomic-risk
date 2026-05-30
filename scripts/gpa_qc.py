#!/usr/bin/env python3
"""
GPA Quality Control Module

Input validation and quality flagging for variants.
Conservative approach: flags anomalies but does NOT reject analysis.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
"""

import json
import re
from pathlib import Path
from typing import List, Dict
from gpa_types import Variant




# =============================================================================
# RepeatMasker (lazy-loaded)
# =============================================================================

_REPEATMASKER_DATA = None  # Lazy-loaded cache


def _load_repeatmasker():
    """Load repeatmasker regions from references."""
    global _REPEATMASKER_DATA
    if _REPEATMASKER_DATA is not None:
        return _REPEATMASKER_DATA
    rm_path = Path(__file__).resolve().parent.parent / "references" / "repeatmasker_regions.json"
    if rm_path.exists():
        with open(rm_path, 'r', encoding='utf-8') as f:
            _REPEATMASKER_DATA = json.load(f)
    else:
        _REPEATMASKER_DATA = []
    return _REPEATMASKER_DATA


def _is_in_repeat_region(chrom, pos):
    """Check if position falls within any repeatmasker region."""
    regions = _load_repeatmasker()
    for region in regions:
        if str(region.get("chrom", "")) == str(chrom) and region.get("start", 0) <= pos <= region.get("end", 0):
            return True
    return False


# =============================================================================
# QC Checks
# =============================================================================

def _run_qc_checks(variants: List[Variant]) -> Dict:
    """
    v0.5 P1-13: Input quality control checks.

    Flags anomalies but does NOT reject analysis - conservative approach.
    Returns QC summary dict for report rendering.
    """
    qc_summary = {
        "total": len(variants),
        "flagged": 0,
        "by_flag": {},
        "flagged_variants": [],
    }

    for v in variants:
        flags = []

        # 1. VAF range check
        if v.vaf is not None and (v.vaf < 0 or v.vaf > 1):
            flags.append("INVALID_VAF")

        # 2. DP depth check
        if v.dp < 10:
            flags.append("LOW_DEPTH")
        elif v.dp >= 10:
            # Check allele support (AD/DP ratio) if AD is available
            # AD not in standard Variant dataclass, skip if unavailable
            pass  # AD field not in core Variant; could be extended later

        # 3. Low complexity / repeat region check
        if _is_in_repeat_region(v.chrom, v.pos):
            flags.append("LOW_COMPLEXITY_REGION")

        # 4. Gene symbol format check (P1-2 HGNC rules)
        gene = v.gene
        if gene:
            # Invalid if starts with digit, too long, or contains illegal chars
            if gene[0].isdigit():
                flags.append("INVALID_GENE_SYMBOL")
            elif len(gene) > 50:
                flags.append("INVALID_GENE_SYMBOL")
            elif not re.match(r'^[A-Za-z][A-Za-z0-9\-]*$', gene):
                flags.append("INVALID_GENE_SYMBOL")

        # v0.5 P1-14: HGNC validation - check transcript_warning for HGNC status
        if v.transcript_warning:
            try:
                tw = json.loads(v.transcript_warning)
                hgnc_warn = tw.get("hgnc_warning", {})
                hgnc_status = hgnc_warn.get("status", "")
                # Withdrawn, not_found, query_failed, or unvalidated_offline = invalid
                if hgnc_status in ("withdrawn", "not_found", "query_failed", "unvalidated_offline"):
                    if "INVALID_GENE_SYMBOL" not in flags:
                        flags.append("INVALID_GENE_SYMBOL")
            except (json.JSONDecodeError, AttributeError):
                pass

        # 5. VAF-GT consistency check (v0.5.3)
        gt_raw = str(v.gt) if v.gt is not None else ""
        if gt_raw in ('.', './.', '.|.', 'nan', 'None', ''):
            gt_raw = ""
        if v.vaf is not None and gt_raw:
            gt = gt_raw.replace("|", "/")
            vaf = v.vaf
            if gt == "0/1":
                if vaf < 0.20 or vaf > 0.80:
                    flags.append("VAF_GT_MISMATCH")
            elif gt == "1/1":
                if vaf < 0.70:
                    flags.append("VAF_GT_MISMATCH")
            elif gt == "0/0":
                if vaf > 0.10:
                    flags.append("VAF_GT_MISMATCH")

        v.qc_flags = flags

        if flags:
            qc_summary["flagged"] += 1
            qc_summary["flagged_variants"].append({
                "gene": v.gene,
                "chrom": v.chrom,
                "pos": v.pos,
                "flags": flags,
            })
            for f in flags:
                qc_summary["by_flag"][f] = qc_summary["by_flag"].get(f, 0) + 1

    return qc_summary
