#!/usr/bin/env python3
"""
GnomAD Local Frequency Archive Module (v0.10.15)

Provides local disk-based allele frequency storage as a fallback for gnomAD API.
Architecture:
  - SQLite cache (30d TTL) → local TSV archive → gnomAD GraphQL API
  - Local archive accumulates queried variants over time, no pre-download required.
  - Format: bgzip + tabix indexed TSV for fast point lookups.

Usage:
    from dgra_gnomad_local import GnomADLocalArchive
    archive = GnomADLocalArchive()
    result = archive.query(chrom, pos, ref, alt)  # returns dict or None
    archive.save(chrom, pos, ref, alt, af, af_eas, af_nfe)
"""

import gzip
import os
import tempfile
from pathlib import Path
from typing import Dict, Optional

# Optional pysam for tabix queries; fallback to linear scan if unavailable
try:
    import pysam
    _HAS_PYSAM = True
except ImportError:
    _HAS_PYSAM = False


class GnomADLocalArchive:
    """
    Local disk archive for gnomAD allele frequencies.

    File format (bgzipped TSV, 1-based POS):
        CHROM  POS  REF  ALT  AF  AF_eas  AF_nfe  queried_at

    Index: tabix (if pysam available) or linear scan fallback.
    """

    _HEADER = "CHROM\tPOS\tREF\tALT\tAF\tAF_eas\tAF_nfe\tqueried_at\n"

    def __init__(
        self,
        archive_path: Optional[Path] = None,
        enable_write: bool = True,
    ):
        """
        Args:
            archive_path: Path to local archive file (default: refs/gnomad_freq_local.tsv.bgz)
            enable_write: Allow appending new entries (default True)
        """
        if archive_path is None:
            self.archive_path = (
                Path(__file__).resolve().parent.parent
                / "references"
                / "gnomad_freq_local.tsv.bgz"
            )
        else:
            self.archive_path = Path(archive_path)

        self.enable_write = enable_write
        self._tabix = None

        # Ensure parent directory exists
        self.archive_path.parent.mkdir(parents=True, exist_ok=True)

        # Initialize file with header if it doesn't exist
        if not self.archive_path.exists():
            self._init_archive()
        elif self.archive_path.stat().st_size == 0:
            self._init_archive()

        # Try to open tabix index
        if _HAS_PYSAM and self._has_tabix_index():
            try:
                self._tabix = pysam.TabixFile(str(self.archive_path))
            except Exception:
                self._tabix = None

    def _init_archive(self):
        """Create archive file with header."""
        with gzip.open(self.archive_path, "wt", compresslevel=6) as f:
            f.write(self._HEADER)

    def _has_tabix_index(self) -> bool:
        """Check if tabix index file exists."""
        return Path(str(self.archive_path) + ".tbi").exists()

    def _build_index(self):
        """Build tabix index using pysam (requires pysam and bgzip format)."""
        if not _HAS_PYSAM:
            return False
        try:
            # pysam.tabix_index requires bgzip format
            pysam.tabix_index(
                str(self.archive_path),
                preset="generic",
                force=True,
            )
            # Reopen tabix
            if self._tabix:
                self._tabix.close()
            self._tabix = pysam.TabixFile(str(self.archive_path))
            return True
        except Exception as e:
            print(f"[GnomADLocal] Failed to build tabix index: {e}")
            return False

    def query(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
    ) -> Optional[Dict]:
        """
        Query local archive for a variant.

        Args:
            chrom: Chromosome (with or without 'chr' prefix)
            pos: 1-based position
            ref: Reference allele
            alt: Alternate allele

        Returns:
            Dict with frequency data if found, None otherwise.
        """
        if not self.archive_path.exists():
            return None

        # Normalize chromosome
        chrom_std = chrom.replace("chr", "").replace("CHR", "") if chrom.upper().startswith("CHR") else chrom
        chrom_with_chr = f"chr{chrom_std}" if not chrom.upper().startswith("CHR") else chrom

        # Try tabix query first (fast)
        if self._tabix:
            try:
                for row in self._tabix.fetch(chrom_std, pos - 1, pos):
                    parts = row.strip().split("\t")
                    if len(parts) >= 8:
                        if (
                            int(parts[1]) == pos
                            and parts[2] == ref
                            and parts[3] == alt
                        ):
                            return {
                                "chrom": parts[0],
                                "pos": int(parts[1]),
                                "ref": parts[2],
                                "alt": parts[3],
                                "af": self._parse_float(parts[4]),
                                "af_eas": self._parse_float(parts[5]),
                                "af_nfe": self._parse_float(parts[6]),
                                "queried_at": parts[7],
                                "source": "gnomad_local_archive",
                                "confidence": "medium",
                            }
            except Exception:
                pass  # Fallback to linear scan

        # Linear scan fallback (slower but always works)
        try:
            with gzip.open(self.archive_path, "rt") as f:
                next(f)  # Skip header
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) < 8:
                        continue
                    if (
                        parts[0] in (chrom_std, chrom_with_chr)
                        and int(parts[1]) == pos
                        and parts[2] == ref
                        and parts[3] == alt
                    ):
                        return {
                            "chrom": parts[0],
                            "pos": int(parts[1]),
                            "ref": parts[2],
                            "alt": parts[3],
                            "af": self._parse_float(parts[4]),
                            "af_eas": self._parse_float(parts[5]),
                            "af_nfe": self._parse_float(parts[6]),
                            "queried_at": parts[7],
                            "source": "gnomad_local_archive",
                            "confidence": "medium",
                        }
        except Exception:
            return None

        return None

    @staticmethod
    def _parse_float(val: str):
        """Parse float or return None."""
        if not val or val == "N/A":
            return None
        try:
            return float(val)
        except (ValueError, TypeError):
            return None

    def save(
        self,
        chrom: str,
        pos: int,
        ref: str,
        alt: str,
        af: Optional[float],
        af_eas: Optional[float],
        af_nfe: Optional[float],
        queried_at: Optional[str] = None,
    ) -> bool:
        """
        Append a variant frequency record to the local archive.

        Args:
            chrom: Chromosome (with or without 'chr' prefix)
            pos: 1-based position
            ref: Reference allele
            alt: Alternate allele
            af: Global allele frequency
            af_eas: EAS population AF
            af_nfe: NFE population AF
            queried_at: ISO timestamp (default: now)

        Returns:
            True if saved successfully.
        """
        if not self.enable_write:
            return False

        import datetime

        if queried_at is None:
            queried_at = datetime.datetime.now().isoformat()

        # Normalize chromosome (store without 'chr' prefix for consistency)
        chrom_std = chrom.replace("chr", "").replace("CHR", "") if chrom.upper().startswith("CHR") else chrom

        line = (
            f"{chrom_std}\t{pos}\t{ref}\t{alt}\t"
            f"{af if af is not None else 'N/A'}\t"
            f"{af_eas if af_eas is not None else 'N/A'}\t"
            f"{af_nfe if af_nfe is not None else 'N/A'}\t"
            f"{queried_at}\n"
        )

        try:
            # Append to bgzip file
            with gzip.open(self.archive_path, "at", compresslevel=6) as f:
                f.write(line)

            # Rebuild tabix index periodically (every 100 writes)
            # Simple heuristic: check file size mod
            if _HAS_PYSAM and self.archive_path.stat().st_size % (100 * 200) < 200:
                self._build_index()

            return True
        except Exception as e:
            print(f"[GnomADLocal] Failed to save record: {e}")
            return False

    def stats(self) -> Dict:
        """Return archive statistics."""
        if not self.archive_path.exists():
            return {"records": 0, "size_mb": 0, "indexed": False}

        try:
            size_mb = self.archive_path.stat().st_size / (1024 * 1024)
            records = 0
            with gzip.open(self.archive_path, "rt") as f:
                next(f)  # Skip header
                for _ in f:
                    records += 1
            return {
                "records": records,
                "size_mb": round(size_mb, 2),
                "indexed": self._has_tabix_index(),
                "path": str(self.archive_path),
            }
        except Exception:
            return {"records": 0, "size_mb": 0, "indexed": False}

    def close(self):
        """Close tabix file handle if open."""
        if self._tabix:
            self._tabix.close()
            self._tabix = None
