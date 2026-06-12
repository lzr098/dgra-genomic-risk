"""
Local HGNC symbol lookup for GPA offline mode.

Reads ~/.workbuddy/data/hgnc/hgnc_lookup.json built by
~/.workbuddy/scripts/build_hgnc_local.py and provides an API-compatible
dict mapping original symbol -> HGNC query result.

This allows GPA offline mode to validate/normalize gene symbols with the
same coverage as the online HGNC REST API, without network access.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Dict, Optional


DEFAULT_HGNC_LOOKUP_PATH = Path.home() / ".workbuddy/data/hgnc/hgnc_lookup.json"


class LocalHGNC:
    """Local HGNC symbol resolver."""

    def __init__(self, lookup_path: Path | str | None = None):
        self.lookup_path = Path(lookup_path or os.getenv(
            "HGNC_LOOKUP_PATH", DEFAULT_HGNC_LOOKUP_PATH
        ))
        self._lookup: Dict[str, Dict[str, Any]] | None = None
        self._metadata: Dict[str, Any] | None = None

    def _load(self) -> None:
        if self._lookup is not None:
            return
        if not self.lookup_path.exists():
            self._lookup = {}
            self._metadata = {"missing": True}
            return
        try:
            with self.lookup_path.open("r", encoding="utf-8") as fh:
                data = json.load(fh)
            self._lookup = data.get("lookup", {})
            self._metadata = data.get("metadata", {})
        except Exception:
            self._lookup = {}
            self._metadata = {"load_error": True}

    def is_available(self) -> bool:
        self._load()
        return bool(self._lookup)

    def metadata(self) -> Dict[str, Any]:
        self._load()
        return dict(self._metadata or {})

    def resolve(self, symbol: str) -> Optional[Dict[str, Any]]:
        """Return HGNC result dict for a symbol, or None if not found."""
        self._load()
        if not symbol:
            return None
        return self._lookup.get(symbol.upper())

    def batch_resolve(self, symbols: list[str]) -> Dict[str, Dict[str, Any]]:
        """Resolve multiple symbols into a dict compatible with HGNC API results."""
        self._load()
        results: Dict[str, Dict[str, Any]] = {}
        for symbol in symbols:
            result = self.resolve(symbol)
            if result:
                # Return a copy so callers can mutate safely
                results[symbol] = dict(result)
            else:
                results[symbol] = {
                    "input": symbol,
                    "approved_symbol": symbol,
                    "hgnc_id": None,
                    "status": "not_found",
                    "previous_symbols": [],
                    "alias_symbols": [],
                    "locus_type": None,
                    "source": "local_hgnc",
                    "confidence": "high",
                }
        return results


def load_local_hgnc(lookup_path: Path | str | None = None) -> LocalHGNC:
    """Factory function returning a loaded LocalHGNC instance."""
    return LocalHGNC(lookup_path)


def local_hgnc_available(lookup_path: Path | str | None = None) -> bool:
    """Quick check whether a local HGNC lookup file exists and is non-empty."""
    return LocalHGNC(lookup_path).is_available()
