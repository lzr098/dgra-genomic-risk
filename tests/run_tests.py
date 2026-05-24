#!/usr/bin/env python3
"""
GPA v0.9.4 L1-L5 Test Suite Runner
No pytest dependency — pure Python assert-based runner.
Usage: python3 tests/run_tests.py
"""

import sys
from pathlib import Path

# Ensure paths
SCRIPT_DIR = Path(__file__).parent.parent / "scripts"
if str(SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(SCRIPT_DIR))

TESTS_DIR = Path(__file__).parent
if str(TESTS_DIR) not in sys.path:
    sys.path.insert(0, str(TESTS_DIR))


from conftest import run_tests


def main():
    modules = [
        ("L1 Static", "test_l1_static.py"),
        ("L2 Unit", "test_l2_unit.py"),
        ("L3 Integration", "test_l3_integration.py"),
    ]

    total_passed = 0
    total_failed = 0

    print("=" * 70)
    print("GPA v0.9.4 Test Suite — L1 / L2 / L3")
    print("=" * 70)

    for name, filename in modules:
        print(f"\n{'='*70}")
        print(f"Running {name} ({filename})")
        print("=" * 70)
        filepath = TESTS_DIR / filename
        if not filepath.exists():
            print(f"  SKIP: {filepath} not found")
            continue

        # Execute the test module in a clean namespace
        namespace = {"__name__": "__test_runner__", "__file__": str(filepath)}
        exec(compile(filepath.read_text(), str(filepath), "exec"), namespace)

        # Extract test functions (matching test_* pattern)
        test_funcs = [
            namespace[k] for k in sorted(namespace.keys())
            if k.startswith("test_") and callable(namespace[k])
        ]

        p, f = run_tests(name, test_funcs)
        total_passed += p
        total_failed += f

    print(f"\n{'='*70}")
    print(f"OVERALL: {total_passed} passed, {total_failed} failed")
    print("=" * 70)

    return 0 if total_failed == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
