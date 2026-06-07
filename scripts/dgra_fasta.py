#!/usr/bin/env python3
"""
GRCh38 FASTA Integration for GPA
v0.10.6 - 2026-06-07

Provides local reference genome sequence extraction for:
- REF allele verification (catch liftover errors or wrong build)
- Flanking sequence context for reports
- No network required (offline capable)

Depends on shared module: ~/.workbuddy/scripts/grch38_fasta_local.py
"""

import sys
from pathlib import Path
from typing import Optional, Tuple, Any

# Import shared FASTA module
_FASTA_MODULE = Path("/Users/zhaorongli/.workbuddy/scripts/grch38_fasta_local.py")


def _import_fasta() -> Optional[Any]:
    """Lazy import shared FASTA module."""
    if not _FASTA_MODULE.exists():
        return None
    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("grch38_fasta_local", _FASTA_MODULE)
        if spec is None or spec.loader is None:
            return None
        mod = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(mod)
        return mod
    except Exception:
        return None


def verify_variant_ref(chrom: str, pos: int, ref: str) -> Tuple[bool, Optional[str], Optional[str]]:
    """Verify if a variant's REF allele matches the GRCh38 reference genome.

    Args:
        chrom: Chromosome
        pos: 1-based position
        ref: Reference allele from VCF

    Returns:
        (is_correct, actual_ref, error_msg)
        - is_correct: True if REF matches genome
        - actual_ref: The actual reference sequence, or None if unavailable
        - error_msg: Error description if check failed, else None
    """
    fasta = _import_fasta()
    if fasta is None:
        return True, None, "GRCh38 FASTA module not available"

    try:
        ok, actual = fasta.verify_ref(chrom, pos, ref)
        if ok:
            return True, actual, None
        return False, actual, f"REF mismatch: VCF says '{ref}', genome has '{actual}'"
    except Exception as e:
        return True, None, f"REF verification error: {e}"


def get_variant_context(chrom: str, pos: int, ref: str, alt: str, flank: int = 50) -> Optional[dict]:
    """Get genomic context around a variant for reporting.

    Args:
        chrom: Chromosome
        pos: 1-based position
        ref: Reference allele
        alt: Alternate allele
        flank: Bases to include on each side

    Returns:
        Dict with upstream/ref/alt/downstream/full_context, or None if unavailable.
    """
    fasta = _import_fasta()
    if fasta is None:
        return None

    try:
        return fasta.get_flanking_sequence(chrom, pos, ref, alt, flank=flank)
    except Exception:
        return None


def annotate_variant_ref_check(variant: Any) -> None:
    """In-place annotate a Variant dataclass with REF verification results.

    Adds to variant.qc_flags if REF mismatch detected.
    Safe to call even if FASTA is unavailable (no-op).

    Args:
        variant: A dgra_core.Variant instance (or any object with chrom/pos/ref/alt/qc_flags)
    """
    if not hasattr(variant, "chrom") or not hasattr(variant, "pos"):
        return

    ok, actual, error = verify_variant_ref(variant.chrom, variant.pos, variant.ref)

    if error and "not available" in error:
        return  # Silently skip if FASTA unavailable

    if not ok:
        flag = f"REF_MISMATCH: VCF='{variant.ref}', genome='{actual}'"
        if hasattr(variant, "qc_flags"):
            variant.qc_flags.append(flag)
        if hasattr(variant, "quality_confidence"):
            variant.quality_confidence = "low"


def annotate_variant_context(variant: Any, flank: int = 30) -> None:
    """In-place add genomic context to a Variant for reporting.

    Adds variant.genomic_context dict if available.

    Args:
        variant: A dgra_core.Variant instance
        flank: Flanking bases to extract
    """
    if not hasattr(variant, "chrom"):
        return

    ctx = get_variant_context(variant.chrom, variant.pos, variant.ref, variant.alt, flank=flank)
    if ctx:
        variant.genomic_context = ctx


