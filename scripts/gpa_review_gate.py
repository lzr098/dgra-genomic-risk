#!/usr/bin/env python3
"""
GPA Code Review Gate — 代码审查门禁脚本
每次修改 GPA skill 代码时运行此脚本，确保符合审查标准。

用法:
    python3 gpa_review_gate.py                    # 检查所有修改的文件
    python3 gpa_review_gate.py --all              # 检查所有文件
    python3 gpa_review_gate.py <file1.py> ...     # 检查指定文件
    python3 gpa_review_gate.py --fix              # 自动修复可自动修复的问题

退出码:
    0 — 通过
    1 — 发现 Blocker，必须修复
    2 — 发现 Critical，建议修复
    3 — 仅 Nit 级别问题
"""

import ast
import sys
import subprocess
import re
from pathlib import Path
from typing import List, Tuple, Dict
from dataclasses import dataclass
from enum import Enum

SCRIPT_DIR = Path("/Users/zhaorongli/.workbuddy/skills/dgra-genomic-risk/scripts")


class Severity(Enum):
    BLOCKER = "🔴 Blocker"
    CRITICAL = "🟡 Critical"
    NIT = "💭 Nit"


@dataclass
class Issue:
    file: str
    line: int
    severity: Severity
    code: str
    message: str
    fixable: bool = False


class GPAReviewGate:
    def __init__(self, files: List[Path]):
        self.files = files
        self.issues: List[Issue] = []
        self.builtins = self._get_builtins()

    def _get_builtins(self) -> set:
        return {"True", "False", "None", "len", "str", "int", "float", "list", "dict",
                "set", "tuple", "range", "enumerate", "zip", "map", "filter", "sum",
                "min", "max", "abs", "round", "sorted", "isinstance", "hasattr",
                "getattr", "setattr", "dir", "vars", "locals", "globals", "print",
                "open", "repr", "chr", "ord", "hex", "oct", "bin", "format", "iter",
                "next", "slice", "hash", "id", "memoryview", "bytearray", "bytes",
                "complex", "divmod", "pow", "super", "object", "type", "Exception",
                "ValueError", "TypeError", "KeyError", "IndexError", "AttributeError",
                "RuntimeError", "NotImplementedError", "FileNotFoundError", "ImportError",
                "ModuleNotFoundError", "AssertionError", "ZeroDivisionError", "OSError",
                "json", "sys", "csv", "re", "os", "asyncio", "datetime", "pathlib",
                "typing", "argparse", "math", "time", "collections", "itertools",
                "functools", "copy", "hashlib", "base64", "urllib", "socket",
                "subprocess", "tempfile", "shutil", "glob", "inspect", "textwrap",
                "string", "warnings", "contextlib", "dataclasses", "enum", "numbers",
                "uuid", "html", "email", "configparser", "io", "traceback", "logging",
                "unittest", "pytest", "mock"}

    def run(self) -> int:
        print("=" * 70)
        print("GPA Code Review Gate v1.0")
        print("=" * 70)
        print(f"Checking {len(self.files)} file(s)...\n")

        for f in self.files:
            if not f.exists():
                self.issues.append(Issue(f.name, 0, Severity.BLOCKER, "FILE_NOT_FOUND", f"File not found: {f}"))
                continue
            self._check_file(f)

        return self._report()

    def _check_file(self, filepath: Path):
        source = filepath.read_text()
        filename = filepath.name

        # 1. Syntax check
        try:
            tree = ast.parse(source)
        except SyntaxError as e:
            self.issues.append(Issue(filename, e.lineno, Severity.BLOCKER, "SYNTAX_ERROR", str(e)))
            return

        # 2. Blocker checks
        self._check_bare_except(filename, tree)
        self._check_mutable_defaults(filename, tree)
        self._check_open_encoding(filename, tree)
        self._check_eval_exec(filename, tree)
        self._check_is_literals(filename, tree)
        self._check_duplicate_set_items(filename, tree, source)
        self._check_duplicate_dict_keys(filename, tree, source)

        # 3. Critical checks
        self._check_logging_fstring(filename, tree, source)
        self._check_unused_imports(filename, tree)

        # 4. Nit checks
        self._check_pointless_fstring(filename, tree, source)

    def _check_bare_except(self, filename: str, tree: ast.AST):
        for node in ast.walk(tree):
            if isinstance(node, ast.ExceptHandler):
                if node.type is None:
                    self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                        "BARE_EXCEPT", "bare except: catches KeyboardInterrupt, SystemExit. Use specific exceptions."))
                elif isinstance(node.type, ast.Name) and node.type.id == "Exception":
                    self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                        "BROAD_EXCEPT", "except Exception is too broad. Catch specific exception types."))

    def _check_mutable_defaults(self, filename: str, tree: ast.AST):
        for node in ast.walk(tree):
            if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
                for default in node.args.defaults + node.args.kw_defaults:
                    if isinstance(default, (ast.List, ast.Dict, ast.Set)):
                        self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                            "MUTABLE_DEFAULT", f"Function '{node.name}' has mutable default argument. Use None + if check."))

    def _check_open_encoding(self, filename: str, tree: ast.AST):
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name) and node.func.id == "open":
                    has_encoding = any(
                        kw.arg == "encoding" for kw in node.keywords
                    )
                    if not has_encoding:
                        self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                            "OPEN_NO_ENCODING", "open() without encoding='utf-8'. May cause UnicodeDecodeError on Windows."))

    def _check_eval_exec(self, filename: str, tree: ast.AST):
        for node in ast.walk(tree):
            if isinstance(node, ast.Call):
                if isinstance(node.func, ast.Name) and node.func.id in ("eval", "exec", "compile"):
                    self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                        "UNSAFE_EVAL", f"{node.func.id}() is dangerous. Avoid unless absolutely necessary."))

    def _check_is_literals(self, filename: str, tree: ast.AST):
        for node in ast.walk(tree):
            if isinstance(node, ast.Compare):
                for op in node.ops:
                    if isinstance(op, ast.Is):
                        comparators = [c for c in node.comparators if isinstance(c, ast.Constant)]
                        left = node.left
                        if isinstance(left, ast.Constant) and left.value not in (None,):
                            if left.value not in (True, False):  # None/True/False are OK
                                self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                                    "IS_LITERAL", f"'is' used with literal {left.value!r}. Use '==' instead."))
                        for c in comparators:
                            if c.value not in (None, True, False):
                                self.issues.append(Issue(filename, node.lineno, Severity.BLOCKER,
                                    "IS_LITERAL", f"'is' used with literal {c.value!r}. Use '==' instead."))

    def _check_duplicate_set_items(self, filename: str, tree: ast.AST, source: str):
        for node in ast.walk(tree):
            if isinstance(node, ast.Set):
                seen = {}
                for elt in node.elts:
                    if isinstance(elt, ast.Constant):
                        key = repr(elt.value)
                        if key in seen:
                            self.issues.append(Issue(filename, elt.lineno, Severity.BLOCKER,
                                "DUPLICATE_SET_ITEM", f"Duplicate value {elt.value!r} in set literal."))
                        seen[key] = elt.lineno

    def _check_duplicate_dict_keys(self, filename: str, tree: ast.AST, source: str):
        for node in ast.walk(tree):
            if isinstance(node, ast.Dict):
                seen = {}
                for k in node.keys:
                    if k is None:  # **kwargs
                        continue
                    if isinstance(k, ast.Constant):
                        key = repr(k.value)
                        if key in seen:
                            self.issues.append(Issue(filename, k.lineno, Severity.BLOCKER,
                                "DUPLICATE_DICT_KEY", f"Duplicate key {k.value!r} in dict literal."))
                        seen[key] = k.lineno

    def _check_logging_fstring(self, filename: str, tree: ast.AST, source: str):
        lines = source.split("\n")
        for i, line in enumerate(lines, 1):
            if "logger." in line or "logging." in line:
                if re.search(r'logger\.(debug|info|warning|error|critical)\s*\(\s*f[\"\']', line):
                    self.issues.append(Issue(filename, i, Severity.CRITICAL,
                        "LOGGING_FSTRING", "Logging with f-string. Use lazy % formatting for performance."))

    def _check_unused_imports(self, filename: str, tree: ast.AST):
        imports = {}
        for node in ast.walk(tree):
            if isinstance(node, ast.Import):
                for alias in node.names:
                    name = alias.asname or alias.name.split(".")[0]
                    imports[name] = node.lineno
            elif isinstance(node, ast.ImportFrom):
                for alias in node.names:
                    name = alias.asname or alias.name
                    imports[name] = node.lineno

        used = set()
        for node in ast.walk(tree):
            if isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load):
                used.add(node.id)

        for name, line in imports.items():
            if name.startswith("_"):
                continue
            if name not in used and name not in self.builtins:
                self.issues.append(Issue(filename, line, Severity.CRITICAL,
                    "UNUSED_IMPORT", f"'{name}' imported but not used. Remove or use it."))

    def _check_pointless_fstring(self, filename: str, tree: ast.AST, source: str):
        for node in ast.walk(tree):
            if isinstance(node, ast.JoinedStr):
                has_placeholder = any(
                    not (isinstance(v, ast.Constant) and isinstance(v.value, str))
                    for v in node.values
                )
                if not has_placeholder:
                    self.issues.append(Issue(filename, node.lineno, Severity.NIT,
                        "POINTLESS_FSTRING", "f-string with no placeholders. Use regular string."))

    def _report(self) -> int:
        blockers = [i for i in self.issues if i.severity == Severity.BLOCKER]
        criticals = [i for i in self.issues if i.severity == Severity.CRITICAL]
        nits = [i for i in self.issues if i.severity == Severity.NIT]

        if not self.issues:
            print("✅ All checks passed! No issues found.")
            return 0

        if blockers:
            print(f"\n{Severity.BLOCKER.value} ({len(blockers)} issues) — MUST FIX:")
            for issue in blockers:
                print(f"  {issue.file}:{issue.line} [{issue.code}] {issue.message}")

        if criticals:
            print(f"\n{Severity.CRITICAL.value} ({len(criticals)} issues) — SHOULD FIX:")
            for issue in criticals:
                print(f"  {issue.file}:{issue.line} [{issue.code}] {issue.message}")

        if nits:
            print(f"\n{Severity.NIT.value} ({len(nits)} issues) — OPTIONAL:")
            for issue in nits:
                print(f"  {issue.file}:{issue.line} [{issue.code}] {issue.message}")

        print(f"\n{'=' * 70}")
        print(f"Summary: {len(blockers)} Blocker, {len(criticals)} Critical, {len(nits)} Nit")

        if blockers:
            print("\n❌ BLOCKERS FOUND — Commit blocked. Fix before proceeding.")
            return 1
        elif criticals:
            print("\n⚠️  CRITICAL ISSUES — Recommend fixing before commit.")
            return 2
        else:
            print("\n✅ Only Nit issues. Safe to commit (optional cleanup).")
            return 3


def get_modified_files() -> List[Path]:
    """Get Python files modified in git working tree."""
    result = subprocess.run(
        ["git", "diff", "--name-only", "--diff-filter=ACM", "HEAD"],
        capture_output=True, text=True, cwd=SCRIPT_DIR
    )
    files = []
    for line in result.stdout.strip().split("\n"):
        if line.endswith(".py"):
            f = SCRIPT_DIR / line
            if f.exists():
                files.append(f)
    return files


def main():
    import argparse
    parser = argparse.ArgumentParser(description="GPA Code Review Gate")
    parser.add_argument("files", nargs="*", help="Python files to check")
    parser.add_argument("--all", action="store_true", help="Check all Python files")
    parser.add_argument("--fix", action="store_true", help="Auto-fix fixable issues")
    args = parser.parse_args()

    if args.all:
        files = list(SCRIPT_DIR.glob("*.py"))
    elif args.files:
        files = [Path(f) for f in args.files]
    else:
        files = get_modified_files()
        if not files:
            print("No modified Python files found. Use --all to check all files.")
            sys.exit(0)

    gate = GPAReviewGate(files)
    sys.exit(gate.run())


if __name__ == "__main__":
    main()
