#!/usr/bin/env python3
"""
GPA VCF Annotator Module (v0.10.3)

Handles raw (unannotated) VCF input by:
1. Detecting whether a VCF is already annotated (CSQ/INFO fields)
2. Detecting genome build from VCF header
3. Lightweight pre-filtering (QUAL<20 or DP<10)
4. Annotating via Docker VEP (preferred), local VEP, or REST API (fallback)
5. Returning all transcript consequences per variant

Auto-detection priority:
  1. Docker VEP (offline cache or database mode) — handles any dataset size
  2. System `vep` command — for small datasets only
  3. VEP REST API — zero-config fallback, rate-limited

v0.10.3: Fixed incorrect 5000-variant threshold that forced VEP API for large
  datasets, bypassing Docker VEP. Docker VEP now correctly prioritized.
v0.10.3: Docker VEP supports --database mode when no local cache exists,
  auto-mounting user's GRCh38 FASTA if available.
"""

import asyncio
import gzip
import json
import logging
import os
import re
import subprocess
import tempfile
from datetime import datetime
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional, Tuple

import aiohttp

try:
    from api_hub import APIHub
except Exception:
    APIHub = None  # type: ignore[misc,assignment]

from dgra_input_parsers import VEP_ANNOTATION_FIELD

logger = logging.getLogger(__name__)

# Ensembl VEP REST API
VEP_API_URL = "https://rest.ensembl.org/vep/human/region"

# Default parameters to fetch all transcripts with rich annotation
VEP_DEFAULT_PARAMS = {
    "canonical": "1",
    "mane_select": "1",
    "mane_plus_clinical": "1",
    "domains": "1",
    "protein": "1",
    "hgvs": "1",
    "numbers": "1",
    "pick": "0",  # Do NOT pick one transcript — return all
}


class VEPBatchFailureError(Exception):
    """Raised when one or more VEP batches fail and user chooses to abort."""

    def __init__(self, failed_batches: List[Dict[str, Any]], message: str = ""):
        self.failed_batches = failed_batches
        self.message = message or f"{len(failed_batches)} VEP batch(es) failed"
        super().__init__(self.message)


class VCFAnnotator:
    """Annotate raw VCF files using VEP API or local VEP."""

    def __init__(
        self,
        annotator: str = "auto",
        genome: str = "auto",
        batch_size: int = 200,
        max_concurrency: int = 5,
        timeout: int = 30,
        vep_cache: Optional[str] = None,
        proxy: Optional[str] = None,
        proxy_route_map: Optional[Any] = None,
        vep_params: Optional[Dict[str, str]] = None,
        shard_dir: Optional[str] = None,
        resume: bool = False,
        checkpoint_path: Optional[str] = None,
        interactive: bool = True,
        hub: Optional["APIHub"] = None,
    ):
        """
        Args:
            annotator: "auto", "vep_api", "vep_local", "vep_docker", "annovar", "snpeff"
            genome: "auto", "GRCh37", "GRCh38"
            batch_size: variants per batch (VEP API chunk size)
            max_concurrency: max concurrent API requests
            timeout: request timeout in seconds
            vep_cache: path to local VEP cache (for vep_local and vep_docker)
            proxy: None = use system proxy, "__DIRECT__" = disable proxy
            proxy_route_map: gpa_proxy_routes.ProxyRouteMap — per-API proxy routing
            vep_params: extra VEP API parameters merged with defaults,
                e.g. {"check_existing": "1", "SIFT": "1", "PolyPhen": "1", "CADD": "1"}
            shard_dir: directory for shard-based incremental annotation storage.
                When set, VEP results are saved per-shard (1k variants each).
            resume: when True and shard_dir is set, skip already-annotated shards.
            checkpoint_path: path to JSON checkpoint file. If exists and non-empty,
                annotate() loads from it instead of calling VEP. After successful
                annotation, results are saved to this path for resume on restart.
            interactive: when True and VEP batches fail, pause and prompt user
                for action (retry / skip / abort). When False, failed batches are
                silently marked with error dicts (legacy behavior).
        """
        self.annotator = annotator
        self.genome = genome
        self.batch_size = batch_size
        self.max_concurrency = max_concurrency
        self.timeout = timeout
        self.vep_cache = vep_cache
        self.vep_docker_cache = vep_cache or os.path.expanduser("~/.workbuddy/tools/vep/cache")
        self.proxy = proxy
        self.proxy_route_map = proxy_route_map
        self.vep_params = vep_params or {}
        self.shard_dir = shard_dir
        self.resume = resume
        self.checkpoint_path = checkpoint_path
        self.interactive = interactive
        self.hub = hub
        self._own_hub = hub is None
        self._session: Optional[aiohttp.ClientSession] = None

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    async def annotate(
        self,
        vcf_path: str,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Dict[str, Any]]:
        """
        Main entry: detect input type, pre-filter, annotate, return enriched variants.
        If checkpoint_path is set and the file exists, loads from checkpoint instead
        of calling VEP. After successful annotation, saves results to checkpoint.

        Returns:
            List of variant dicts, each containing:
                chrom, pos, ref, alt, qual, filter, dp, gt,
                and transcript_consequences (list of all transcript annotations).
        """
        # v0.10.13: Checkpoint resume — skip VEP if checkpoint exists
        if self.checkpoint_path and os.path.exists(self.checkpoint_path):
            size = os.path.getsize(self.checkpoint_path)
            if size > 100:
                try:
                    with open(self.checkpoint_path, 'r', encoding='utf-8') as f:
                        annotated = json.load(f)
                    logger.info(
                        f"[VCFAnnotator] Checkpoint loaded: {len(annotated)} variants "
                        f"from {self.checkpoint_path}"
                    )
                    return annotated
                except Exception as e:  # noqa: BROAD_EXCEPT — checkpoint corruption is non-fatal, re-run VEP
                    logger.warning(
                        f"[VCFAnnotator] Checkpoint corrupt ({e}), re-running VEP..."
                    )

        try:
            annotated = await self._annotate_internal(vcf_path, progress_callback)
        except Exception as e:  # noqa: BROAD_EXCEPT — process-level guard: ensures session cleanup on any failure
            logger.error(f"[VCFAnnotator] Annotation failed: {type(e).__name__}: {e}")
            raise
        finally:
            await self.close()

        # Save checkpoint on success
        if self.checkpoint_path and annotated:
            try:
                with open(self.checkpoint_path, 'w', encoding='utf-8') as f:
                    json.dump(annotated, f, ensure_ascii=False, indent=1)
                logger.info(
                    f"[VCFAnnotator] Checkpoint saved: {len(annotated)} variants "
                    f"to {self.checkpoint_path}"
                )
            except Exception as e:  # noqa: BROAD_EXCEPT — checkpoint save failure is non-fatal
                logger.warning(f"[VCFAnnotator] Failed to save checkpoint: {e}")

        return annotated

    async def _annotate_internal(
        self,
        vcf_path: str,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Dict[str, Any]]:
        """Internal annotate implementation."""
        vcf_path = Path(vcf_path)
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF not found: {vcf_path}")

        # 1. Detect genome build
        genome_build = self._detect_genome(vcf_path) if self.genome == "auto" else self.genome
        logger.info(f"[VCFAnnotator] Genome build detected: {genome_build}")

        # 2. Parse raw variants + pre-filter
        raw_variants = self._parse_vcf(vcf_path)
        total_raw = len(raw_variants)
        filtered_variants = self._prefilter(raw_variants)
        logger.info(
            f"[VCFAnnotator] Parsed {total_raw} variants, "
            f"retained {len(filtered_variants)} after pre-filter (QUAL≥20, DP≥10)"
        )

        if not filtered_variants:
            return []

        # 3. Resolve annotator
        resolved = self._resolve_annotator(len(filtered_variants))
        logger.info(f"[VCFAnnotator] Using annotator: {resolved}")

        # 4. Annotate
        if resolved == "vep_api":
            annotated = await self._annotate_vep_api(
                filtered_variants, genome_build, progress_callback
            )
        elif resolved == "vep_local":
            annotated = await self._annotate_vep_local(
                filtered_variants, genome_build, progress_callback
            )
        elif resolved == "vep_docker":
            annotated = await self._annotate_vep_docker(
                filtered_variants, genome_build, progress_callback
            )
        else:
            raise NotImplementedError(f"Annotator '{resolved}' not yet implemented in v0.9.0")

        # 5. Merge annotation back into variant dicts
        # v0.12.0: Merge ALL fields from annotation (including CSQ-style fields)
        # Old code only copied transcript_consequences + vep_summary, dropping GENE/IMPACT/etc.
        results = []
        for v, ann in zip(filtered_variants, annotated):
            v.update(ann)
            results.append(v)

        logger.info(f"[VCFAnnotator] Annotation complete for {len(results)} variants")
        return results

    def is_annotated_vcf(self, vcf_path: str) -> bool:
        """Check if VCF already contains annotation (CSQ in INFO)."""
        path = Path(vcf_path)
        opener = self._vcf_opener(path)
        try:
            with opener(path, "rt") as fh:
                for line in fh:
                    if line.startswith("#"):
                        if VEP_ANNOTATION_FIELD in line or "Consequence" in line:
                            return True
                    else:
                        # Only scan first 50 non-header lines
                        break
        except Exception as e:
            logger.warning(f"[VCFAnnotator] Cannot check annotation status: {e}")
        return False

    # ------------------------------------------------------------------
    # Genome detection
    # ------------------------------------------------------------------

    def _detect_genome(self, vcf_path: Path) -> str:
        """Infer genome build from VCF header."""
        opener = self._vcf_opener(vcf_path)
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if not line.startswith("##"):
                    break
                low = line.lower()
                if "grch38" in low or "hg38" in low or "b38" in low:
                    return "GRCh38"
                if "grch37" in low or "hg19" in low or "b37" in low:
                    return "GRCh37"
        # Fallback: count variants on chrM length as heuristic (crude)
        return "GRCh38"

    def _vcf_opener(self, vcf_path: Path):
        """Return appropriate file opener based on actual file content, not just extension."""
        try:
            with open(vcf_path, "rb") as fh:
                magic = fh.read(2)
                if magic == b'\x1f\x8b':
                    return gzip.open
        except Exception:
            pass
        return open

    # ------------------------------------------------------------------
    # VCF parsing
    # ------------------------------------------------------------------

    def _parse_vcf(self, vcf_path: Path) -> List[Dict[str, Any]]:
        """Extract CHROM, POS, REF, ALT, QUAL, FILTER, DP, GT from VCF.

        v0.12.2: Added boundary handling for symbolic alleles and low-quality records.
        """
        opener = self._vcf_opener(vcf_path)
        variants: List[Dict[str, Any]] = []
        with opener(vcf_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                chrom, pos, _id, ref, alt, qual, filt, info = parts[:8]

                # v0.12.2: Skip records with missing or failed FILTER
                if filt not in (".", "PASS", ""):
                    continue

                # v0.12.2: Skip symbolic/non-reference alleles and SV breakpoints
                # <NON_REF>, *, and breakends (]chr:pos[) are not supported
                alts = alt.split(",")
                valid_alts = []
                for a in alts:
                    if a == ".":
                        continue
                    if a.startswith("<") or a == "*":
                        continue
                    if any(c in a for c in ("[", "]")):
                        continue
                    valid_alts.append(a)
                if not valid_alts:
                    continue

                for a in valid_alts:
                    v = {
                        "chrom": self._normalize_chrom(chrom),
                        "pos": int(pos),
                        "ref": ref,
                        "alt": a,
                        "qual": float(qual) if qual != "." else 0.0,
                        "filter": filt,
                    }
                    # Extract DP from INFO
                    v["dp"] = self._extract_dp(info)
                    # Extract GT from FORMAT/SAMPLE if available
                    if len(parts) >= 10:
                        v["gt"] = self._extract_gt(parts[8], parts[9])
                    else:
                        v["gt"] = "./."
                    variants.append(v)
        return variants

    @staticmethod
    def _normalize_chrom(chrom: str) -> str:
        """Strip 'chr' prefix for Ensembl compatibility."""
        chrom = chrom.upper()
        if chrom.startswith("CHR"):
            chrom = chrom[3:]
        # Ensembl uses MT not M
        if chrom == "M":
            chrom = "MT"
        return chrom

    @staticmethod
    def _extract_dp(info: str) -> Optional[int]:
        """Extract DP from INFO field."""
        for field in info.split(";"):
            if field.startswith("DP="):
                try:
                    return int(field.split("=", 1)[1])
                except ValueError:
                    return None
        return None

    @staticmethod
    def _extract_gt(format_col: str, sample_col: str) -> str:
        """Extract GT from FORMAT/SAMPLE columns."""
        fmt = format_col.split(":")
        if "GT" not in fmt:
            return "./."
        gt_idx = fmt.index("GT")
        vals = sample_col.split(":")
        if gt_idx < len(vals):
            return vals[gt_idx]
        return "./."

    # ------------------------------------------------------------------
    # Pre-filtering
    # ------------------------------------------------------------------

    @staticmethod
    def _prefilter(variants: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
        """Lightweight filter: exclude QUAL<20 or DP<10."""
        kept = []
        for v in variants:
            qual = v.get("qual", 0)
            dp = v.get("dp")
            if qual < 20:
                continue
            if dp is not None and dp < 10:
                continue
            kept.append(v)
        return kept

    # ------------------------------------------------------------------
    # Annotator resolution
    # ------------------------------------------------------------------

    def _resolve_annotator(self, n_variants: int = 0) -> str:
        """Resolve auto → concrete annotator.

        Priority in auto mode:
        1. Docker VEP (vep_docker) — offline cache or database mode, no rate limits,
           handles any dataset size efficiently.
        2. System vep command (vep_local) — for small datasets only (<5000 variants).
        3. VEP REST API (vep_api) — zero-config fallback, may be rate-limited.

        v0.10.3 FIX: Removed the erroneous 5000-variant threshold that forced VEP API
        for large datasets, completely bypassing Docker VEP. Docker VEP with local
        cache annotates ~100-500 variants/sec — 50k variants take minutes, not hours.
        Even database mode (no cache) is faster than REST API for large batches due
        to no rate-limiting and single-roundtrip batch processing.

        Local `vep` command still capped at 5000 variants because it runs inline
        and blocks the event loop; Docker VEP has no such limitation.
        """
        if self.annotator != "auto":
            return self.annotator

        # 1. Docker VEP — ALWAYS check first, regardless of variant count
        docker_available, docker_mode = self._vep_docker_available()
        if docker_available:
            logger.info(
                f"[VCFAnnotator] Auto-detected Docker VEP ({docker_mode} mode) → "
                f"using vep_docker ({n_variants} variants)"
            )
            return "vep_docker"

        # 2. System vep command — only for small datasets (blocks event loop)
        LOCAL_VEP_MAX_VARIANTS = 5000
        if n_variants <= LOCAL_VEP_MAX_VARIANTS and self._vep_local_available():
            logger.info(
                f"[VCFAnnotator] Auto-detected local VEP → using vep_local "
                f"({n_variants} variants <= {LOCAL_VEP_MAX_VARIANTS} threshold)"
            )
            return "vep_local"

        # 3. REST API — last resort
        logger.info(
            f"[VCFAnnotator] No local annotators available → using vep_api (REST) "
            f"for {n_variants} variants"
        )
        return "vep_api"

    @staticmethod
    def _vep_local_available() -> bool:
        """Check if `vep` command exists."""
        try:
            subprocess.run(
                ["vep", "--help"],
                capture_output=True,
                timeout=5,
                check=False,
            )
            return True
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False

    @staticmethod
    def _vep_docker_available(cache_dir: str = "~/.workbuddy/tools/vep/cache") -> Tuple[bool, str]:
        """Check if Docker VEP is available and determine best run mode.

        Returns:
            Tuple of (available: bool, mode: str) where mode is:
            - "offline"  : local cache exists → --cache --offline (fastest)
            - "database" : no cache but Docker + image available → --database
                           (connects to Ensembl DB; slower but no rate limits)
            - ""         : not available

        v0.10.3: Now returns mode string instead of bool, allowing callers to
        choose between --cache/offline and --database modes. Also detects user's
        local GRCh38 FASTA for use with --database --fasta.
        """
        # 1. Docker daemon
        try:
            result = subprocess.run(
                ["docker", "info"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if result.returncode != 0:
                return False, ""
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False, ""

        # 2. VEP image
        try:
            result = subprocess.run(
                ["docker", "image", "inspect", "ensemblorg/ensembl-vep:latest"],
                capture_output=True,
                text=True,
                timeout=10,
            )
            if result.returncode != 0:
                return False, ""
        except (FileNotFoundError, subprocess.TimeoutExpired):
            return False, ""

        # 3. Check for local cache (offline mode — preferred)
        expanded = os.path.expanduser(cache_dir)
        species_dir = os.path.join(expanded, "homo_sapiens")
        if os.path.isdir(species_dir):
            for root, _, files in os.walk(species_dir):
                for f in files:
                    if f.endswith(".gz"):
                        return True, "offline"

        # 4. No cache — still usable in database mode (no rate limits vs REST API)
        return True, "database"

    # ------------------------------------------------------------------
    # Shard-based incremental annotation (P1 fix: architecture + resume)
    # ------------------------------------------------------------------

    SHARD_SIZE = 1000  # variants per shard

    def _shard_path(self, shard_idx: int) -> Path:
        """Get path for a shard's annotation JSON file."""
        return Path(self.shard_dir) / f"shard_{shard_idx:05d}.json"

    def _index_path(self) -> Path:
        """Get path for the missing-by-shard index file."""
        return Path(self.shard_dir) / "missing_by_shard.json"

    def _get_completed_shards(self) -> set:
        """Return set of shard indices that already have complete annotation files."""
        if not self.shard_dir:
            return set()
        index_path = self._index_path()
        if not index_path.exists():
            return set()
        try:
            with open(index_path, "r") as f:
                missing = json.load(f)
            total_shards = missing.get("total_shards", 0)
            done = missing.get("completed_shards", [])
            return set(done)
        except Exception:
            return set()

    def _save_shard_atomic(self, shard_idx: int, data: list) -> None:
        """Atomically write a shard annotation file (tmp → rename)."""
        import os as _os
        shard_path = self._shard_path(shard_idx)
        tmp_path = Path(str(shard_path) + ".tmp")
        with open(tmp_path, "w") as f:
            json.dump(data, f, ensure_ascii=False)
        _os.replace(str(tmp_path), str(shard_path))

    def _mark_shard_complete(self, shard_idx: int) -> None:
        """Update missing-by-shard index after completing a shard."""
        index_path = self._index_path()
        tmp_path = Path(str(index_path) + ".tmp")
        try:
            if index_path.exists():
                with open(index_path, "r") as f:
                    index = json.load(f)
            else:
                index = {}
            completed = set(index.get("completed_shards", []))
            completed.add(shard_idx)
            index["completed_shards"] = sorted(completed)
            index["updated_at"] = datetime.now().isoformat()
            with open(tmp_path, "w") as f:
                json.dump(index, f, ensure_ascii=False)
            import os as _os
            _os.replace(str(tmp_path), str(index_path))
        except Exception:
            pass

    def _init_shard_dir(self, total_variants: int) -> None:
        """Initialize shard directory with missing_by_shard.json."""
        if not self.shard_dir:
            return
        Path(self.shard_dir).mkdir(parents=True, exist_ok=True)
        total_shards = (total_variants + self.SHARD_SIZE - 1) // self.SHARD_SIZE
        index_path = self._index_path()
        if not index_path.exists() or not self.resume:
            index = {
                "total_shards": total_shards,
                "shard_size": self.SHARD_SIZE,
                "completed_shards": [],
                "created_at": datetime.now().isoformat(),
            }
            with open(index_path, "w") as f:
                json.dump(index, f, ensure_ascii=False)

    # ------------------------------------------------------------------
    # VEP API annotation
    # ------------------------------------------------------------------

    async def _annotate_vep_api(
        self,
        variants: List[Dict[str, Any]],
        genome: str,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Dict[str, Any]]:
        """Annotate via Ensembl VEP REST API with controlled concurrency.

        v0.10.14: If batches fail after all retries and interactive=True,
        pause and prompt user for action (retry / skip / abort).
        """
        # v0.10.12: per-API proxy routing
        vep_proxy: Optional[str] = None
        if self.proxy_route_map is not None:
            vep_proxy = self.proxy_route_map.get_proxy("ensembl")
            if vep_proxy:
                logger.info(f"[VCFAnnotator] VEP API using proxy: {vep_proxy}")
            else:
                logger.info("[VCFAnnotator] VEP API using direct connection")
        elif self.proxy and self.proxy != "__DIRECT__":
            vep_proxy = self.proxy

        if not self._session:
            if self.hub is not None:
                self._session = self.hub.session
            else:
                timeout = aiohttp.ClientTimeout(total=self.timeout)
                connector = aiohttp.TCPConnector(force_close=True)
                trust_env = self.proxy != "__DIRECT__" and self.proxy_route_map is None
                self._session = aiohttp.ClientSession(
                    connector=connector, timeout=timeout, trust_env=trust_env
                )

        semaphore = asyncio.Semaphore(self.max_concurrency)
        total = len(variants)
        annotated: List[Dict[str, Any]] = [None] * total
        failed_batches: List[Dict[str, Any]] = []

        async def _query_batch(start_idx: int, batch: List[Dict[str, Any]]) -> bool:
            """Query one batch. Returns True on success, False on failure."""
            async with semaphore:
                body = {
                    "variants": [
                        f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ."
                        for v in batch
                    ]
                }
                params = dict(VEP_DEFAULT_PARAMS)
                params.update(self.vep_params)
                if genome == "GRCh37":
                    params["refseq"] = "1"

                backoff = [1, 2, 4]
                last_error = ""
                for attempt, delay in enumerate(backoff + [None]):
                    try:
                        async with self._session.post(
                            VEP_API_URL,
                            params=params,
                            json=body,
                            headers={"Content-Type": "application/json"},
                            proxy=vep_proxy,
                        ) as resp:
                            if resp.status == 200:
                                data = await resp.json()
                                results = self._parse_vep_response(data, batch)
                                for i, r in zip(
                                    range(start_idx, start_idx + len(batch)), results
                                ):
                                    annotated[i] = r
                                return True
                            elif resp.status == 429:
                                retry_after = int(resp.headers.get("Retry-After", delay or 2))
                                logger.warning(
                                    f"VEP API 429, retry after {retry_after}s (attempt {attempt + 1})"
                                )
                                if delay is not None:
                                    await asyncio.sleep(retry_after)
                                    continue
                                last_error = f"HTTP 429 (rate limited)"
                            else:
                                logger.warning(
                                    f"VEP API HTTP {resp.status}, attempt {attempt + 1}"
                                )
                                if delay is not None:
                                    await asyncio.sleep(delay)
                                    continue
                                last_error = f"HTTP {resp.status}"
                    except asyncio.TimeoutError:
                        logger.warning(f"VEP API timeout, attempt {attempt + 1}")
                        if delay is not None:
                            await asyncio.sleep(delay)
                            continue
                        last_error = "TimeoutError"
                    except aiohttp.ClientError as e:
                        logger.warning(f"VEP API client error: {e}, attempt {attempt + 1}")
                        if delay is not None:
                            await asyncio.sleep(delay)
                            continue
                        last_error = f"ClientError: {e}"
                    except OSError as e:
                        logger.warning(f"VEP API network error: {e}, attempt {attempt + 1}")
                        if delay is not None:
                            await asyncio.sleep(delay)
                            continue
                        last_error = f"OSError: {e}"
                    except Exception as e:  # noqa: BROAD_EXCEPT
                        logger.error(f"VEP API unexpected error: {type(e).__name__}: {e}, attempt {attempt + 1}")
                        if delay is not None:
                            await asyncio.sleep(delay)
                            continue
                        last_error = f"{type(e).__name__}: {e}"

                # All retries failed — record for later handling
                failed_batches.append({
                    "start_idx": start_idx,
                    "batch": batch,
                    "error": last_error or "VEP API failed after all retries",
                })
                return False

        # Launch all batches
        tasks = []
        for i in range(0, total, self.batch_size):
            batch = variants[i : i + self.batch_size]
            task = _query_batch(i, batch)
            tasks.append(task)

        # Progress tracking
        completed = 0
        for coro in asyncio.as_completed(tasks):
            await coro
            completed += self.batch_size
            if progress_callback:
                progress_callback(min(completed, total), total)

        # v0.10.14: Handle failed batches interactively
        if failed_batches:
            annotated = await self._handle_failed_batches(
                annotated, failed_batches, variants, genome, semaphore, vep_proxy,
            )

        return annotated

    async def _handle_failed_batches(
        self,
        annotated: List[Dict[str, Any]],
        failed_batches: List[Dict[str, Any]],
        variants: List[Dict[str, Any]],
        genome: str,
        semaphore: asyncio.Semaphore,
        vep_proxy: Optional[str],
    ) -> List[Dict[str, Any]]:
        """Pause and handle failed VEP batches. May retry, skip, or abort."""
        n_failed = len(failed_batches)
        n_variants = sum(len(b["batch"]) for b in failed_batches)

        # Build failure report
        print("\n" + "=" * 70)
        print("[VCFAnnotator] ⚠️  VEP BATCH FAILURE REPORT")
        print("=" * 70)
        print(f"Failed batches: {n_failed}  |  Affected variants: {n_variants}")
        print("-" * 70)
        for fb in failed_batches:
            start = fb["start_idx"]
            batch = fb["batch"]
            error = fb["error"]
            positions = [f"{v['chrom']}:{v['pos']}" for v in batch[:3]]
            ellipsis = " ..." if len(batch) > 3 else ""
            print(f"  Batch {start:>6}-{start + len(batch):<6}  {error:<30}  variants: {', '.join(positions)}{ellipsis}")
        print("-" * 70)

        if self.interactive:
            print("\nOptions:")
            print("  [R]etry  — re-submit failed batches (may resolve transient errors)")
            print("  [S]kip   — mark failed variants with error annotations and continue")
            print("  [A]bort  — stop the analysis and raise VEPBatchFailureError")
            print("-" * 70)

            # Read user input (run in thread to avoid blocking event loop)
            def _read_choice() -> str:
                while True:
                    choice = input("Your choice [R/s/a]: ").strip().lower()
                    if choice in ("", "r", "s", "a"):
                        return choice if choice else "r"
                    print("Invalid choice. Please enter R, s, or a.")

            choice = await asyncio.to_thread(_read_choice)
        else:
            choice = "s"  # Non-interactive: default to skip
            print(f"[VCFAnnotator] interactive=False — auto-skipping {n_failed} failed batches")

        if choice == "r":
            # Retry failed batches
            print(f"[VCFAnnotator] Retrying {n_failed} failed batches...")
            retry_failed: List[Dict[str, Any]] = []
            for fb in failed_batches:
                success = await self._query_single_batch(
                    fb["start_idx"], fb["batch"], genome, semaphore, vep_proxy, annotated
                )
                if not success:
                    retry_failed.append(fb)

            if retry_failed:
                print(f"[VCFAnnotator] ⚠️  {len(retry_failed)} batches still failed after retry.")
                # Second chance: prompt again or auto-skip
                if self.interactive:
                    def _read_second_choice() -> str:
                        while True:
                            c = input("Some batches still failed. [S]kip remaining / [A]bort: ").strip().lower()
                            if c in ("", "s", "a"):
                                return c if c else "s"
                    second_choice = await asyncio.to_thread(_read_second_choice)
                else:
                    second_choice = "s"

                if second_choice == "a":
                    raise VEPBatchFailureError(retry_failed, "Batches failed after retry")
                # Skip remaining
                for fb in retry_failed:
                    for i in range(fb["start_idx"], fb["start_idx"] + len(fb["batch"])):
                        annotated[i] = {
                            "transcript_consequences": [],
                            "vep_summary": {"error": f"VEP API failed: {fb['error']}"},
                        }
            else:
                print("[VCFAnnotator] All failed batches succeeded on retry.")

        elif choice == "a":
            raise VEPBatchFailureError(failed_batches)

        else:  # skip
            for fb in failed_batches:
                for i in range(fb["start_idx"], fb["start_idx"] + len(fb["batch"])):
                    annotated[i] = {
                        "transcript_consequences": [],
                        "vep_summary": {"error": f"VEP API failed: {fb['error']}"},
                    }
            print(f"[VCFAnnotator] Skipped {n_failed} failed batches, marked {n_variants} variants with error annotations.")

        print("=" * 70 + "\n")
        return annotated

    async def _query_single_batch(
        self,
        start_idx: int,
        batch: List[Dict[str, Any]],
        genome: str,
        semaphore: asyncio.Semaphore,
        vep_proxy: Optional[str],
        annotated: List[Dict[str, Any]],
    ) -> bool:
        """Query a single batch (used for retry). Returns True on success."""
        async with semaphore:
            body = {
                "variants": [
                    f"{v['chrom']} {v['pos']} . {v['ref']} {v['alt']} . . ."
                    for v in batch
                ]
            }
            params = dict(VEP_DEFAULT_PARAMS)
            params.update(self.vep_params)
            if genome == "GRCh37":
                params["refseq"] = "1"

            try:
                async with self._session.post(
                    VEP_API_URL,
                    params=params,
                    json=body,
                    headers={"Content-Type": "application/json"},
                    proxy=vep_proxy,
                ) as resp:
                    if resp.status == 200:
                        data = await resp.json()
                        results = self._parse_vep_response(data, batch)
                        for i, r in zip(
                            range(start_idx, start_idx + len(batch)), results
                        ):
                            annotated[i] = r
                        return True
                    else:
                        logger.warning(f"VEP retry batch {start_idx}: HTTP {resp.status}")
                        return False
            except Exception as e:
                logger.warning(f"VEP retry batch {start_idx}: {type(e).__name__}: {e}")
                return False

    @staticmethod
    def _parse_vep_response(
        data: List[Dict[str, Any]], batch: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Parse VEP REST API response into per-variant annotation dicts.

        Guarantees output length == len(batch) even if VEP returns fewer/more results.
        """
        results = []
        # VEP returns one entry per input line
        for entry, v in zip(data, batch):
            if not isinstance(entry, dict):
                results.append({
                    "transcript_consequences": [],
                    "vep_summary": {"error": "Invalid VEP response format"},
                })
                continue

            # Summary info
            summary = {
                "most_severe_consequence": entry.get("most_severe_consequence", ""),
                "variant_class": entry.get("variant_class", ""),
            }
            if "colocated_variants" in entry:
                summary["colocated_variants"] = entry["colocated_variants"]

            # Transcript consequences
            tx_list = []
            for tc in entry.get("transcript_consequences", []):
                tx = {
                    "transcript_id": tc.get("transcript_id", ""),
                    "gene_symbol": tc.get("gene_symbol", ""),
                    "gene_id": tc.get("gene_id", ""),
                    "consequence_terms": tc.get("consequence_terms", []),
                    "impact": tc.get("impact", ""),
                    "hgvsc": tc.get("hgvsc", ""),
                    "hgvsp": tc.get("hgvsp", ""),
                    "canonical": tc.get("canonical", 0),
                    "mane_select": tc.get("mane_select", 0),
                    "mane_plus_clinical": tc.get("mane_plus_clinical", 0),
                    "exon": tc.get("exon", ""),
                    "intron": tc.get("intron", ""),
                    "protein_domains": tc.get("domains", []),
                    # v0.11.5: biotype is critical for transcript selection filtering
                    "biotype": tc.get("biotype", ""),
                }
                tx_list.append(tx)

            # ── v0.12.0: Also write top-level CSQ-style fields for _parse_variants_phase1() ──
            # Pick the best transcript for CSQ fields
            best_tx = None
            for tc in entry.get("transcript_consequences", []):
                if tc.get("canonical") or tc.get("mane_select"):
                    best_tx = tc
                    break
            if best_tx is None and tx_list:
                best_tx = tx_list[0]

            csq_fields = {}
            if best_tx:
                gene_sym = best_tx.get("gene_symbol", "")
                csq_fields["GENE"] = gene_sym
                csq_fields["IMPACT"] = best_tx.get("impact", "UNKNOWN")
                ct = best_tx.get("consequence_terms", [])
                csq_fields["Consequence"] = ",".join(ct) if ct else ""
                csq_fields["HGVSc"] = best_tx.get("hgvsc", "")
                csq_fields["HGVSp"] = best_tx.get("hgvsp", "")
            else:
                csq_fields["GENE"] = ""
                csq_fields["IMPACT"] = "UNKNOWN"
                csq_fields["Consequence"] = ""
                csq_fields["HGVSc"] = ""
                csq_fields["HGVSp"] = ""

            # ClinVar + gnomAD from colocated_variants
            clin_sig_parts = []
            gnomad_afs = []
            if "colocated_variants" in entry:
                for cv in entry["colocated_variants"]:
                    # ClinVar — clin_sig_allele format: "ALLELE:significance"
                    # For indels, ALLELE can be very long; only keep significance part
                    cs = cv.get("clin_sig_allele", "")
                    if cs:
                        # Extract significance after the last ":" (handles long ALT sequences)
                        if ":" in cs:
                            sig = cs.split(":", 1)[1].strip()
                        else:
                            sig = cs.strip()
                        if sig and sig not in clin_sig_parts:
                            clin_sig_parts.append(sig)
                    # gnomAD AF
                    # VEP frequencies structure: {allele: {pop: af}}
                    freqs = cv.get("frequencies", {})
                    if isinstance(freqs, dict):
                        for allele_freqs in freqs.values():
                            if isinstance(allele_freqs, dict):
                                af = allele_freqs.get("gnomade_af", None)
                                if af is not None:
                                    try:
                                        gnomad_afs.append(float(af))
                                    except (ValueError, TypeError):
                                        pass
            csq_fields["CLIN_SIG"] = ",".join(clin_sig_parts) if clin_sig_parts else ""
            csq_fields["gnomad_af"] = min(gnomad_afs) if gnomad_afs else None

            # Merge CSQ fields into result dict
            result = {
                "transcript_consequences": tx_list,
                "vep_summary": summary,
            }
            result.update(csq_fields)
            results.append(result)

        # v0.10.0 P2-2: VEP may return fewer results than input batch (filtered rows)
        # Pad remaining slots with error-marked entries so caller gets exactly len(batch) items
        missing = len(batch) - len(results)
        if missing > 0:
            for _ in range(missing):
                results.append({
                    "transcript_consequences": [],
                    "vep_summary": {"error": "VEP response shorter than input batch"},
                })

        return results

    # ------------------------------------------------------------------
    # Local VEP annotation (async wrapper around subprocess)
    # ------------------------------------------------------------------

    async def _annotate_vep_local(
        self,
        variants: List[Dict[str, Any]],
        genome: str,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Dict[str, Any]]:
        """Annotate using local VEP command."""
        # Write variants to temp VCF
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as tmp_in:
            tmp_in.write("##fileformat=VCFv4.2\n")
            tmp_in.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for v in variants:
                tmp_in.write(
                    f"{v['chrom']}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}"
                    f"\t{v.get('qual', '.')}\t.\t.\n"
                )
            tmp_in_path = tmp_in.name

        out_path = tmp_in_path.replace(".vcf", "_vep.vcf")
        cache_flag = ["--cache", "--cache_version", "108"] if not self.vep_cache else [
            "--cache", "--dir_cache", self.vep_cache
        ]
        assembly = "GRCh38" if genome == "GRCh38" else "GRCh37"

        cmd = [
            "vep",
            "--input_file", tmp_in_path,
            "--output_file", out_path,
            "--vcf",
            "--assembly", assembly,
            "--canonical",
            "--mane",
            "--domains",
            "--protein",
            "--hgvs",
            "--numbers",
            # v0.10.17: Add population AF and existing variant annotation.
            # Fields that are absent in the cache are left empty by VEP rather
            # than filled with placeholder values.
            # "--af_gnomad",  # v0.10.18: VEP 115 Docker lacks gnomAD plugin
            # "--af_gnomad_exome",
            # "--af_gnomad_genome",
            "--check_existing",
            "--pubmed",
            "--fork", "4",
            *cache_flag,
            "--offline" if not self.vep_cache else "",
        ]
        cmd = [c for c in cmd if c]  # Remove empty

        logger.info(f"[VCFAnnotator] Running local VEP: {' '.join(cmd[:8])} ...")
        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()
            if proc.returncode != 0:
                logger.error(f"Local VEP failed: {stderr.decode()[:500]}")
                return [
                    {
                        "transcript_consequences": [],
                        "vep_summary": {"error": f"Local VEP exit {proc.returncode}"},
                    }
                ] * len(variants)
        except Exception as e:
            logger.error(f"Local VEP execution error: {e}")
            return [
                {
                    "transcript_consequences": [],
                    "vep_summary": {"error": str(e)},
                }
            ] * len(variants)

        # Parse VEP output VCF
        # (Simplified: in practice would parse CSQ field from output VCF)
        # For v0.9.0, fall back to parsing the VEP output
        return self._parse_vep_local_output(out_path, variants)

    async def _annotate_vep_docker(
        self,
        variants: List[Dict[str, Any]],
        genome: str,
        progress_callback: Optional[Callable[[int, int], None]] = None,
    ) -> List[Dict[str, Any]]:
        """Annotate using Docker VEP (ensemblorg/ensembl-vep).

        v0.10.3: Supports both --cache/--offline (fastest, requires local cache)
        and --database modes (no cache required, auto-mounts user's FASTA).
        """
        # Write variants to temp VCF
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".vcf", delete=False
        ) as tmp_in:
            tmp_in.write("##fileformat=VCFv4.2\n")
            tmp_in.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for v in variants:
                tmp_in.write(
                    f"{v['chrom']}\t{v['pos']}\t.\t{v['ref']}\t{v['alt']}"
                    f"\t{v.get('qual', '.')}\t.\t.\n"
                )
            tmp_in_path = tmp_in.name

        out_path = tmp_in_path.replace(".vcf", "_vep.vcf")
        assembly = "GRCh38" if genome == "GRCh38" else "GRCh37"

        # Determine run mode and build volumes + args
        _, docker_mode = self._vep_docker_available()
        volumes: List[str] = []
        vep_args: List[str] = []

        if docker_mode == "offline":
            cache_dir = os.path.expanduser(self.vep_docker_cache)
            container_cache = "/data/vep_cache"
            volumes = [
                "-v", f"{cache_dir}:{container_cache}:ro",
            ]
            vep_args = [
                "--cache",
                "--dir_cache", container_cache,
                "--cache_version", "115",
                "--offline",
            ]
        else:
            # database mode — connect to Ensembl DB; no local cache needed
            vep_args = ["--database"]
            # Auto-detect and mount user's local GRCh38 FASTA
            fasta_candidates = [
                os.path.expanduser(
                    "~/.workbuddy/data/genome/"
                    "Homo_sapiens.GRCh38.dna.primary_assembly.fa"
                ),
                os.path.expanduser("~/.workbuddy/data/genome/GRCh38.fa"),
                os.path.expanduser("~/data/genome/Homo_sapiens.GRCh38.dna.primary_assembly.fa"),
            ]
            for fp in fasta_candidates:
                if os.path.isfile(fp):
                    fasta_dir = os.path.dirname(fp)
                    fasta_name = os.path.basename(fp)
                    volumes.extend(["-v", f"{fasta_dir}:/data/fasta:ro"])
                    vep_args.extend(["--fasta", f"/data/fasta/{fasta_name}"])
                    logger.info(
                        f"[VCFAnnotator] Auto-mounted FASTA for database mode: {fp}"
                    )
                    break

        cmd = (
            ["docker", "run", "--rm"]
            + volumes
            + [
                "-v", f"{os.path.dirname(tmp_in_path)}:/data/input:ro",
                "-v", f"{os.path.dirname(out_path)}:/data/output",
                "ensemblorg/ensembl-vep:latest",
                "vep",
                "--input_file", f"/data/input/{os.path.basename(tmp_in_path)}",
                "--output_file", f"/data/output/{os.path.basename(out_path)}",
                "--vcf",
                "--assembly", assembly,
                "--canonical",
                "--mane",
                "--domains",
                "--protein",
                "--hgvs",
                "--numbers",
                # v0.10.17: Add population AF and existing variant annotation.
                # "--af_gnomad",  # v0.10.18: VEP 115 Docker lacks gnomAD plugin
                # "--af_gnomad_exome",
                # "--af_gnomad_genome",
                "--check_existing",
                "--pubmed",
                "--fork", "4",
            ]
            + vep_args
        )

        logger.info(
            f"[VCFAnnotator] Running Docker VEP ({docker_mode}): "
            f"{' '.join(cmd[:6])} ... ({len(variants)} variants)"
        )
        try:
            proc = await asyncio.create_subprocess_exec(
                *cmd,
                stdout=asyncio.subprocess.PIPE,
                stderr=asyncio.subprocess.PIPE,
            )
            stdout, stderr = await proc.communicate()
            if proc.returncode != 0:
                err = stderr.decode()[:800]
                logger.error(f"Docker VEP failed: {err}")
                return [
                    {
                        "transcript_consequences": [],
                        "vep_summary": {"error": f"Docker VEP exit {proc.returncode}: {err}"},
                    }
                ] * len(variants)
        except Exception as e:
            logger.error(f"Docker VEP execution error: {e}")
            return [
                {
                    "transcript_consequences": [],
                    "vep_summary": {"error": str(e)},
                }
            ] * len(variants)

        # Parse VEP output VCF (same format as local VEP)
        return self._parse_vep_local_output(out_path, variants)

    def _parse_vep_local_output(
        self, out_path: str, variants: List[Dict[str, Any]]
    ) -> List[Dict[str, Any]]:
        """Parse VEP-local output VCF (with CSQ in INFO)."""
        results = []
        opener = self._vcf_opener(Path(out_path))
        with opener(out_path, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) < 8:
                    continue
                info = parts[7]
                # Extract VEP annotation field
                csq_match = re.search(rf"{VEP_ANNOTATION_FIELD}=([^;]+)", info)
                if not csq_match:
                    results.append({
                        "transcript_consequences": [],
                        "vep_summary": {},
                    })
                    continue
                csq_str = csq_match.group(1)
                # Parse CSQ (simplified)
                tx_list = []
                for csq in csq_str.split(","):
                    fields = csq.split("|")
                    # VEP CSQ order varies — we need header to know positions
                    # For simplicity, map by position if we know the order
                    # In production, parse ##INFO=<ID=CSQ,...> header
                    tx_list.append({
                        "raw_csq": csq,
                        "consequence_terms": [fields[1]] if len(fields) > 1 else [],
                        "impact": fields[2] if len(fields) > 2 else "",
                        "transcript_id": fields[6] if len(fields) > 6 else "",
                        "gene_symbol": fields[3] if len(fields) > 3 else "",
                        "hgvsc": fields[9] if len(fields) > 9 else "",
                        "hgvsp": fields[10] if len(fields) > 10 else "",
                    })
                results.append({
                    "transcript_consequences": tx_list,
                    "vep_summary": {},
                })
        return results

    # ------------------------------------------------------------------
    # Cleanup
    # ------------------------------------------------------------------

    async def close(self):
        """Close aiohttp session.

        IMPORTANT: Callers MUST explicitly await close() when done.
        Do NOT rely on garbage collection — aiohttp.ClientSession.close()
        is async and cannot be safely triggered from __del__ without
        RuntimeWarning: coroutine 'ClientSession.close' was never awaited.

        ponytail: if a shared APIHub owns the session, leave it alone.
        """
        if self._session and self._own_hub:
            await self._session.close()
            self._session = None
