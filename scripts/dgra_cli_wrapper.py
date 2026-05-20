#!/usr/bin/env python3
"""
DGRA CLI Wrapper — 简化供者基因组风险评估的调用接口
供 OpenClaw agent 直接使用。

用法:
    python3 scripts/dgra_cli_wrapper.py --variants '[{...}]' --tissue hematopoietic
    python3 scripts/dgra_cli_wrapper.py --input-file donor_variants.tsv --tissue hematopoietic

功能:
    1. 将 variant list (JSON) 写入临时 TSV 输入文件
    2. 调用 dgra_core.py 执行分析
    3. 解析 JSON 输出并返回结构化 dict
    4. 失败时返回 error dict，不抛异常
"""

import json
import sys
import tempfile
import csv
import subprocess
from pathlib import Path
from typing import List, Dict, Any, Optional


# DGRA core script path
SCRIPT_DIR = Path(__file__).resolve().parent
DGRA_CORE = SCRIPT_DIR / "dgra_core.py"
REFS_DIR = SCRIPT_DIR.parent / "references"

# Required TSV columns
REQUIRED_COLS = [
    "CHROM", "POS", "REF", "ALT", "GENE", "Feature", "EXON",
    "IMPACT", "Consequence", "HGVSp", "HGVSc", "CLIN_SIG",
    "GT", "DP", "GQ", "VAF", "gnomAD_AF"
]

# Optional columns with defaults
OPTIONAL_DEFAULTS = {
    "Feature": "",
    "EXON": "",
    "IMPACT": "MODERATE",
    "Consequence": "missense_variant",
    "HGVSp": "",
    "HGVSc": "",
    "CLIN_SIG": "",
    "GT": "0/1",
    "DP": "30",
    "GQ": "99",
    "VAF": "0.5",
    "gnomAD_AF": "",
}


def _write_tsv(variants: List[Dict[str, Any]], tsv_path: Path) -> None:
    """将 variant dict list 写入 TSV，补全缺失列。"""
    with open(tsv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=REQUIRED_COLS, delimiter="\t")
        writer.writeheader()
        for v in variants:
            row = {}
            for col in REQUIRED_COLS:
                val = v.get(col)
                if val is None or val == "":
                    val = OPTIONAL_DEFAULTS.get(col, "")
                row[col] = str(val)
            writer.writerow(row)


def _write_patient_mutations(mutations: List[Dict[str, Any]], json_path: Path) -> None:
    """将患者突变列表写入 JSON 文件。"""
    with open(json_path, "w", encoding="utf-8") as f:
        json.dump(mutations, f, indent=2, ensure_ascii=False)


def run_dgra(
    variants: List[Dict[str, Any]],
    tissue: str = "hematopoietic",
    patient_mutations: Optional[List[Dict[str, Any]]] = None,
    offline: bool = False,
) -> Dict[str, Any]:
    """
    运行 DGRA 分析管道。

    Args:
        variants: variant dict 列表，每个 dict 至少包含 CHROM, POS, REF, ALT, GENE
        tissue: 组织类型，默认 hematopoietic
        patient_mutations: 可选，患者体细胞突变列表用于交叉比对
        offline: 是否离线模式（跳过 API）

    Returns:
        dict: {"success": True, "results": {...}, "report_md": "..."}
        或 {"success": False, "error": "..."}
    """
    if not variants:
        return {"success": False, "error": "variants list is empty"}

    # Validate tissue
    valid_tissues = {"hematopoietic", "cardiovascular", "hepatic", "renal", "neurological"}
    if tissue not in valid_tissues:
        return {
            "success": False,
            "error": f"Invalid tissue '{tissue}'. Valid: {', '.join(sorted(valid_tissues))}",
        }

    # 创建临时文件
    with tempfile.TemporaryDirectory(prefix="dgra_wrapper_") as tmpdir:
        tmp = Path(tmpdir)
        tsv_path = tmp / "variants.tsv"
        json_out = tmp / "results.json"
        md_out = tmp / "report.md"

        # 写输入
        try:
            _write_tsv(variants, tsv_path)
        except Exception as e:
            return {"success": False, "error": f"Failed to write TSV: {e}"}

        # 构造命令行
        cmd = [
            sys.executable,
            str(DGRA_CORE),
            "--input", str(tsv_path),
            "--tissue", tissue,
            "--output", str(md_out),
            "--json", str(json_out),
        ]
        if offline:
            cmd.append("--offline")

        # 患者突变（可选）
        patient_json = None
        if patient_mutations:
            patient_json = tmp / "patient_mutations.json"
            _write_patient_mutations(patient_mutations, patient_json)
            cmd.extend(["--patient-mutations", str(patient_json)])

        # 执行
        try:
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                timeout=300,  # 5 分钟上限
                cwd=str(SCRIPT_DIR),
            )
        except subprocess.TimeoutExpired:
            return {"success": False, "error": "DGRA analysis timed out after 5 minutes"}
        except Exception as e:
            return {"success": False, "error": f"Subprocess failed: {e}"}

        if result.returncode != 0:
            return {
                "success": False,
                "error": f"dgra_core.py exited with code {result.returncode}",
                "stderr": result.stderr,
                "stdout": result.stdout,
            }

        # 解析输出
        try:
            with open(json_out, "r", encoding="utf-8") as f:
                results = json.load(f)
        except Exception as e:
            return {
                "success": False,
                "error": f"Failed to parse JSON output: {e}",
                "stdout": result.stdout,
            }

        try:
            with open(md_out, "r", encoding="utf-8") as f:
                report_md = f.read()
        except Exception:
            report_md = ""

        return {
            "success": True,
            "results": results,
            "report_md": report_md,
            "stdout": result.stdout,
        }


def main():
    import argparse

    parser = argparse.ArgumentParser(description="DGRA CLI Wrapper")
    parser.add_argument(
        "--variants",
        required=True,
        help='JSON array of variant dicts, e.g. \'[{"CHROM":"1","POS":12345,"REF":"A","ALT":"G","GENE":"VWF"}]\'',
    )
    parser.add_argument("--tissue", default="hematopoietic", help="Tissue profile")
    parser.add_argument("--patient-mutations", help="JSON array of patient mutations")
    parser.add_argument("--offline", action="store_true", help="Offline mode")
    parser.add_argument("--output-json", help="Write result JSON to this file")

    args = parser.parse_args()

    try:
        variants = json.loads(args.variants)
    except json.JSONDecodeError as e:
        print(json.dumps({"success": False, "error": f"Invalid variants JSON: {e}"}, indent=2))
        sys.exit(1)

    patient_mutations = None
    if args.patient_mutations:
        try:
            patient_mutations = json.loads(args.patient_mutations)
        except json.JSONDecodeError as e:
            print(json.dumps({"success": False, "error": f"Invalid patient_mutations JSON: {e}"}, indent=2))
            sys.exit(1)

    result = run_dgra(
        variants=variants,
        tissue=args.tissue,
        patient_mutations=patient_mutations,
        offline=args.offline,
    )

    output = json.dumps(result, indent=2, ensure_ascii=False, default=str)

    if args.output_json:
        with open(args.output_json, "w", encoding="utf-8") as f:
            f.write(output)
    else:
        print(output)


if __name__ == "__main__":
    main()
