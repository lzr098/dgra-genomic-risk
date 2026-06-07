#!/usr/bin/env python3
"""OMIM integration module for GPA (gpa-genomic-phenotype).

v0.10.6 - 2026-06-07

Provides OMIM gene-disease annotations for tier classification.
Wraps the shared ~/.workbuddy/scripts/omim_local.py module.

Key additions:
- is_rare_disease_gene(gene): supersedes gene_phenotype_map.json lookup
- get_gene_disease_info(gene): full gene→disease annotation for tier evidence
- get_inheritance(gene): inheritance pattern for variant interpretation
"""

import sys
from pathlib import Path
from typing import Any

# Load shared OMIM module
_SHARED_OMIM = Path.home() / ".workbuddy/scripts/omim_local.py"

_omim = None
_omim_available = False


def _load_omim() -> Any:
    """Lazy-load the shared OMIM module."""
    global _omim, _omim_available
    if _omim is not None:
        return _omim

    if not _SHARED_OMIM.exists():
        _omim_available = False
        _omim = None
        return None

    try:
        import importlib.util
        spec = importlib.util.spec_from_file_location("omim_local", _SHARED_OMIM)
        if spec is None or spec.loader is None:
            _omim_available = False
            return None
        _omim = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(_omim)
        _omim_available = True
        return _omim
    except Exception:
        _omim_available = False
        _omim = None
        return None


def is_rare_disease_gene(gene: str) -> bool:
    """Check if gene has Mendelian disease association in OMIM.

    This supersedes the static gene_phenotype_map.json lookup.
    Queries the local OMIM SQLite database (29K records).

    Args:
        gene: HGNC gene symbol

    Returns:
        True if gene is associated with Mendelian disease(s)
    """
    omim = _load_omim()
    if omim is None:
        return False

    try:
        return omim.is_mendelian_gene(gene)
    except Exception:
        return False


def get_gene_disease_info(gene: str) -> dict[str, Any] | None:
    """Get comprehensive gene-disease annotation from OMIM.

    Args:
        gene: HGNC gene symbol

    Returns:
        Dict with gene info + associated phenotypes, or None if not found
        {
            "gene": {mim_number, title, chromosome, start, end},
            "phenotypes": [{mim_number, title, inheritance, clinical_synopsis}],
            "total_records": int
        }
    """
    omim = _load_omim()
    if omim is None:
        return None

    try:
        return omim.get_gene_phenotype(gene)
    except Exception:
        return None


def get_inheritance(gene: str) -> str | None:
    """Get primary inheritance pattern for a gene.

    Args:
        gene: HGNC gene symbol

    Returns:
        Inheritance pattern string (e.g. 'AD', 'AR', 'XLD') or None
    """
    info = get_gene_disease_info(gene)
    if not info or not info.get("phenotypes"):
        return None

    # Return the inheritance of the first phenotype entry
    inheritances = [
        p.get("inheritance")
        for p in info["phenotypes"]
        if p.get("inheritance")
    ]
    return inheritances[0] if inheritances else None


def is_available() -> bool:
    """Check if OMIM database is available."""
    _load_omim()
    return _omim_available
