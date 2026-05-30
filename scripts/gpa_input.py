#!/usr/bin/env python3
"""
GPA Input Handling Module

Detects input file types and parses VEP-annotated variants into
dgra_core's internal format.

Extracted from dgra_core.py in v0.10.0 God Module refactoring.
"""

import json
from enum import Enum
from typing import List, Dict, Any, Optional

from dgra_input_parsers import VEP_ANNOTATION_FIELD


class InputType(Enum):
    RAW_VCF = "raw_vcf"
    ANNOTATED_VCF = "annotated_vcf"
    ANNOTATED_TABLE = "annotated_table"
    FREE_TEXT = "free_text"
    UNKNOWN = "unknown"


def detect_input_type(input_path: str) -> InputType:
    """Detect input file type: raw VCF, annotated VCF, or annotated table."""
    path_lower = input_path.lower()

    if path_lower.endswith(('.vcf', '.vcf.gz', '.bcf')):
        if _has_vcf_annotation(input_path):
            return InputType.ANNOTATED_VCF
        else:
            return InputType.RAW_VCF
    elif path_lower.endswith(('.tsv', '.csv', '.xlsx', '.xls', '.xlsm')):
        return InputType.ANNOTATED_TABLE
    elif path_lower.endswith(('.txt', '.md')):
        return InputType.FREE_TEXT
    else:
        # Try reading content to determine
        try:
            with open(input_path, 'r', encoding='utf-8') as f:
                first_lines = [f.readline() for _ in range(5)]
                content = ''.join(first_lines)
                if '##fileformat=VCF' in content or '#CHROM\tPOS' in content:
                    if _has_vcf_annotation_from_content(content):
                        return InputType.ANNOTATED_VCF
                    else:
                        return InputType.RAW_VCF
                elif 'CHROM' in content or 'Gene' in content or 'Consequence' in content:
                    return InputType.ANNOTATED_TABLE
                else:
                    return InputType.FREE_TEXT
        except (IndexError, ValueError):
            return InputType.UNKNOWN


def _has_vcf_annotation(vcf_path: str) -> bool:
    """Check if VCF has CSQ or ANN annotation in INFO."""
    import gzip
    opener = gzip.open if vcf_path.endswith('.gz') else open
    try:
        with opener(vcf_path, 'rt', encoding='utf-8') as f:
            for i, line in enumerate(f):
                if i > 1000:
                    break
                if line.startswith(f'##INFO=<ID={VEP_ANNOTATION_FIELD}'):
                    return True
                if line.startswith('##INFO=<ID=ANN'):
                    return True
                if not line.startswith('#') and f'{VEP_ANNOTATION_FIELD}=' in line:
                    return True
    except (IndexError, ValueError):
        pass
    return False


def _has_vcf_annotation_from_content(content: str) -> bool:
    """Check VCF annotation from content string."""
    return f'##INFO=<ID={VEP_ANNOTATION_FIELD}' in content or '##INFO=<ID=ANN' in content


def variants_from_vep_annotation(
    annotated_variants: List[Dict[str, Any]],
    selector: Optional[Any] = None,
) -> List[Dict[str, Any]]:
    """
    Convert VCFAnnotator output (with transcript_consequences) into dgra_core variants_data format.
    If selector is provided, perform disease-aware transcript selection.
    """
    variants_data = []
    for v in annotated_variants:
        tx_consequences = v.get("transcript_consequences", [])
        if not tx_consequences:
            # No transcript consequences — still include with minimal info
            variants_data.append({
                "CHROM": v.get("chrom", ""),
                "POS": str(v.get("pos", "")),
                "REF": v.get("ref", ""),
                "ALT": v.get("alt", ""),
                "Gene": "",
                "Consequence": "",
                "IMPACT": "",
                "HGVSc": "",
                "HGVSp": "",
                "DP": str(v.get("dp", "")),
                "QUAL": str(v.get("qual", "")),
                "GT": v.get("gt", ""),
                "GQ": str(v.get("gq", "")) if v.get("gq") is not None else "",
                "VAF": str(v.get("vaf", "")) if v.get("vaf") is not None else "",
            })
            continue

        # Group by gene
        gene_txs: Dict[str, List[Dict]] = {}
        for tx in tx_consequences:
            gene = tx.get("gene_symbol", "")
            if not gene:
                continue
            gene_txs.setdefault(gene, []).append(tx)

        # For each gene, select primary transcript
        for gene, txs in gene_txs.items():
            if selector:
                result = selector.select(gene, txs)
                primary = result.primary
                alternatives = result.alternatives
            else:
                # Fallback: pick canonical or first
                primary = next((t for t in txs if t.get("canonical")), txs[0])
                alternatives = [t for t in txs if t != primary]

            # Build variant dict in dgra_core expected format
            vd = {
                "CHROM": v.get("chrom", ""),
                "POS": str(v.get("pos", "")),
                "REF": v.get("ref", ""),
                "ALT": v.get("alt", ""),
                "Gene": gene,
                "Feature": primary.get("transcript_id", ""),
                "Consequence": ",".join(primary.get("consequence_terms", [])),
                "IMPACT": primary.get("impact", ""),
                "HGVSc": primary.get("hgvsc", ""),
                "HGVSp": primary.get("hgvsp", ""),
                # v0.10.1: Extract ClinVar from VEP annotation (saves MyVariant re-query)
                "CLIN_SIG": primary.get("clin_sig", ""),
                "DP": str(v.get("dp", "")),
                "QUAL": str(v.get("qual", "")),
                "GT": v.get("gt", ""),
                "GQ": str(v.get("gq", "")) if v.get("gq") is not None else "",
                "VAF": str(v.get("vaf", "")) if v.get("vaf") is not None else "",
                # v0.9.0 transcript selection metadata
                "primary_transcript": primary.get("transcript_id", ""),
                "primary_hgvsc": primary.get("hgvsc", ""),
                "primary_hgvsp": primary.get("hgvsp", ""),
                "primary_consequence": ",".join(primary.get("consequence_terms", [])),
                "primary_impact": primary.get("impact", ""),
                "alternative_transcripts": json.dumps(alternatives) if alternatives else "",
                "transcript_selection_method": getattr(selector, "method", "canonical") if selector else "canonical",
                "transcript_ambiguity_flag": "",
            }
            variants_data.append(vd)

    return variants_data
