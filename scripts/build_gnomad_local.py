#!/usr/bin/env python3
"""
GnomAD v4.1 Local Frequency Builder (v0.10.15)

Downloads gnomAD v4.1 exome sites VCF by chromosome, extracts high-frequency
variants (AF > 1%), and builds a lightweight local archive for fast lookups.

Requirements:
    - bcftools (installed and in PATH)
    - ~50GB temporary storage (input VCFs are deleted after processing)
    - ~200-300MB final output

Usage:
    python build_gnomad_local.py [--af-threshold 0.01] [--output refs/gnomad_freq_local.tsv.bgz]

Output:
    - gnomad_freq_local.tsv.bgz (bgzip + tabix indexed)
    - Columns: CHROM, POS, REF, ALT, AF, AF_eas, AF_nfe
"""

import argparse
import os
import subprocess
import tempfile
from pathlib import Path
from typing import List, Optional

# gnomAD v4.1 exome sites on Google Cloud Storage
GNOMAD_V4_EXOME_BASE = (
    "gs://gcp-public-data--gnomad/release/4.1/vcf/exomes/"
    "gnomad.exomes.v4.1.sites.chr{chrom}.vcf.bgz"
)

CHROMOSOMES = [str(i) for i in range(1, 23)] + ["X", "Y"]


def run_cmd(cmd: List[str], cwd: Optional[Path] = None) -> None:
    """Run shell command and check return code."""
    print(f"[build] {' '.join(cmd)}")
    result = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
    if result.returncode != 0:
        raise RuntimeError(f"Command failed: {' '.join(cmd)}\nstderr: {result.stderr}")


def download_chromosome(chrom: str, out_dir: Path) -> Path:
    """Download gnomAD v4.1 exome sites VCF for one chromosome."""
    url = GNOMAD_V4_EXOME_BASE.format(chrom=chrom)
    out_path = out_dir / f"gnomad.chr{chrom}.vcf.bgz"

    if out_path.exists():
        print(f"[build] chr{chrom}: already downloaded, skipping")
        return out_path

    print(f"[build] Downloading chr{chrom} from gnomAD...")
    # Use gsutil if available, otherwise curl via public HTTP URL
    gsutil_ok = subprocess.run(["which", "gsutil"], capture_output=True).returncode == 0
    if gsutil_ok:
        run_cmd(["gsutil", "cp", url, str(out_path)])
    else:
        # Fallback: public HTTPS URL (may be slower)
        https_url = url.replace("gs://", "https://storage.googleapis.com/")
        run_cmd(["curl", "-L", "-o", str(out_path), https_url])

    return out_path


def extract_frequencies(
    vcf_path: Path,
    output_tsv: Path,
    af_threshold: float = 0.01,
) -> None:
    """Extract AF, AF_eas, AF_nfe from gnomAD VCF to TSV."""
    # bcftools query format: extract INFO/AF, INFO/AF_eas, INFO/AF_nfe
    # gnomAD v4.1 INFO fields: AF, AF_afr, AF_amr, AF_eas, AF_fin, AF_mid, AF_nfe, AF_sas, etc.
    fmt = "%CHROM\t%POS\t%REF\t%ALT\t%INFO/AF\t%INFO/AF_eas\t%INFO/AF_nfe\n"

    # Filter: overall AF > threshold OR any population AF > threshold
    expr = f"INFO/AF>{af_threshold} || INFO/AF_eas>{af_threshold} || INFO/AF_nfe>{af_threshold}"

    print(f"[build] Extracting AF>{af_threshold} from {vcf_path.name}...")

    with open(output_tsv, "a") as out_f:
        proc = subprocess.run(
            [
                "bcftools", "query",
                "-f", fmt,
                "-i", expr,
                str(vcf_path),
            ],
            capture_output=True,
            text=True,
        )
        if proc.returncode != 0:
            print(f"[build] Warning: bcftools failed for {vcf_path.name}: {proc.stderr}")
            return
        out_f.write(proc.stdout)


def build_index(tsv_path: Path) -> None:
    """Sort, bgzip, and tabix index the TSV."""
    print(f"[build] Sorting and indexing {tsv_path.name}...")

    sorted_tsv = tsv_path.with_suffix(".sorted.tsv")
    run_cmd(["sort", "-k1,1", "-k2,2n", "-o", str(sorted_tsv), str(tsv_path)])

    # Add header
    header_tsv = tsv_path.with_suffix(".header.tsv")
    with open(header_tsv, "w") as f:
        f.write("CHROM\tPOS\tREF\tALT\tAF\tAF_eas\tAF_nfe\n")

    final_tsv = tsv_path.with_suffix(".final.tsv")
    run_cmd(["cat", str(header_tsv), str(sorted_tsv), "-o", str(final_tsv)])

    # bgzip
    bgz_path = tsv_path.with_suffix(".tsv.bgz")
    run_cmd(["bgzip", "-f", str(final_tsv)])
    bgz_final = final_tsv.with_suffix(".tsv.bgz")
    if bgz_final.exists():
        bgz_final.rename(bgz_path)

    # tabix index
    run_cmd(["tabix", "-s", "1", "-b", "2", "-e", "2", str(bgz_path)])

    # Cleanup intermediates
    for p in [sorted_tsv, header_tsv, final_tsv]:
        if p.exists():
            p.unlink()

    print(f"[build] Indexed archive: {bgz_path}")
    print(f"[build] Index: {bgz_path}.tbi")


def main():
    parser = argparse.ArgumentParser(description="Build local gnomAD frequency archive")
    parser.add_argument(
        "--af-threshold",
        type=float,
        default=0.01,
        help="Minimum allele frequency threshold (default: 0.01 = 1%%)",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("references/gnomad_freq_local.tsv.bgz"),
        help="Output path for local archive",
    )
    parser.add_argument(
        "--tmp-dir",
        type=Path,
        default=None,
        help="Temporary directory for downloads (default: system temp)",
    )
    parser.add_argument(
        "--chromosomes",
        type=str,
        default=",".join(CHROMOSOMES),
        help="Comma-separated chromosome list (default: 1-22,X,Y)",
    )
    parser.add_argument(
        "--keep-vcf",
        action="store_true",
        help="Keep downloaded VCF files (default: delete after processing)",
    )
    args = parser.parse_args()

    chromosomes = args.chromosomes.split(",")
    out_dir = args.output.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    tmpdir = Path(args.tmp_dir) if args.tmp_dir else Path(tempfile.mkdtemp(prefix="gnomad_build_"))
    tmpdir.mkdir(parents=True, exist_ok=True)

    # Intermediate TSV (unsorted)
    raw_tsv = tmpdir / "gnomad_freq_raw.tsv"
    if raw_tsv.exists():
        raw_tsv.unlink()

    print(f"[build] Building gnomAD local archive")
    print(f"[build] AF threshold: {args.af_threshold}")
    print(f"[build] Chromosomes: {chromosomes}")
    print(f"[build] Temp dir: {tmpdir}")
    print(f"[build] Output: {args.output}")

    try:
        for chrom in chromosomes:
            vcf_path = download_chromosome(chrom, tmpdir)
            extract_frequencies(vcf_path, raw_tsv, af_threshold=args.af_threshold)

            if not args.keep_vcf:
                vcf_path.unlink()
                print(f"[build] Deleted {vcf_path.name}")

        # Sort, bgzip, index
        build_index(raw_tsv)

        # Move final archive to output path
        final_bgz = raw_tsv.with_suffix(".tsv.bgz")
        final_tbi = raw_tsv.with_suffix(".tsv.bgz.tbi")
        if final_bgz.exists():
            final_bgz.rename(args.output)
        if final_tbi.exists():
            final_tbi.rename(Path(str(args.output) + ".tbi"))

        # Stats
        size_mb = args.output.stat().st_size / (1024 * 1024)
        print(f"[build] Done! Archive size: {size_mb:.1f} MB")
        print(f"[build] Location: {args.output}")

    finally:
        # Cleanup temp dir
        if not args.keep_vcf and tmpdir.exists():
            for f in tmpdir.iterdir():
                f.unlink()
            tmpdir.rmdir()
            print(f"[build] Cleaned up temp dir: {tmpdir}")


if __name__ == "__main__":
    main()
