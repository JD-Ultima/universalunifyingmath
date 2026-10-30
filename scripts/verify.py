#!/usr/bin/env python3
"""Run selected assertions and regressions without rewriting recorded outputs.

This suite excludes long searches, optional SDP/element checks, manuscript
compilation, and fresh Lean builds. --upstream adds text/record consistency
checks against the exact recorded Chain Theory commit.
"""
import argparse
import importlib.metadata
import os
from pathlib import Path
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
SOURCE_COMMIT = "a28db8998bce5bd5ef11195c05f13c86cd9dac5c"
DEPENDENCIES = {"numpy": "2.5.3", "sympy": "1.14.0", "mpmath": "1.3.0",
                "scipy": "1.18.1", "z3-solver": "5.3.3.0"}
SELECTED = [
    "chain-theory/alpha/records_independent.py",
    "chain-theory/alpha/config_enumeration.py",
    "chain-theory/alpha/algebraic.py",
    "chain-theory/alpha/mutations_classified.py",
    "chain-theory/alpha/freedom.py",
    "chain-theory/alpha/xcheck_operator_form.py",
    "geometric-unity/A_spin_obstruction.py",
    "geometric-unity/xcheck_round5.py",
    "geometric-unity/projection_regression.py",
    "howard/arith_candidates.py",
    "howard/axiom_lattice.py",
    "howard/semirings.py",
]


def source_commands(path):
    compendium = path / "compendium"
    if not (compendium / "index/statements.tsv").is_file():
        raise ValueError("--upstream must point to a Chain Theory checkout")
    head = subprocess.run(["git", "-C", str(path), "rev-parse", "HEAD"],
                          capture_output=True, text=True, encoding="utf-8", timeout=30)
    if head.returncode or head.stdout.strip() != SOURCE_COMMIT:
        raise ValueError(f"--upstream must be checked out at {SOURCE_COMMIT}")
    dirty = subprocess.run(["git", "-C", str(path), "status", "--porcelain",
                            "--untracked-files=all", "--ignored"],
                           capture_output=True, text=True, encoding="utf-8", timeout=30)
    if dirty.returncode or dirty.stdout.strip():
        raise ValueError("--upstream must be pristine, with no tracked changes, untracked files, or ignored files; "
                         "use a fresh clone checked out at the recorded commit")
    prefix = ROOT / "computation/chain-theory/compendium"
    logs = prefix / "outputs"
    return [
        ("source census", [str(prefix / "census.py"), str(compendium / "index/statements.tsv")]),
        ("source premises", [str(prefix / "premises.py"), str(compendium)]),
        ("source hydrogen registers", [str(prefix / "hydrogen_registers.py"), str(compendium)]),
        ("Lean coverage using recorded logs", [
            str(prefix / "lean_tables.py"), str(compendium),
            "--axioms", str(logs / "lean_axioms.txt"), "--checks", str(logs / "lean_table_checks.txt"),
            "--round6-tsv", str(logs / "lean_resolve.tsv")]),
    ]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--upstream", type=Path, help="optional pristine Chain Theory checkout (no extra files)")
    parser.add_argument("--timeout", type=float, default=180, help="seconds per check (default: 180)")
    parser.add_argument("--list", action="store_true", help="list selected checks without running them")
    args = parser.parse_args()
    if args.timeout <= 0:
        parser.error("--timeout must be positive")
    checks = [("focused regressions", ["-m", "unittest", "discover", "-s", "tests", "-v"])]
    checks += [(item, [str(ROOT / "computation" / item)]) for item in SELECTED]
    if args.upstream:
        try:
            checks += source_commands(args.upstream.resolve())
        except (ValueError, OSError, subprocess.TimeoutExpired) as error:
            parser.error(str(error))
    if args.list:
        for label, _ in checks:
            print(label)
        return 0
    missing = []
    for package, expected in DEPENDENCIES.items():
        try:
            actual = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            actual = "missing"
        if actual != expected:
            missing.append(f"{package}: expected {expected}, found {actual}")
    if missing:
        print("\n".join(missing), file=sys.stderr)
        print("Install dependencies: python -m pip install -r requirements.txt", file=sys.stderr)
        return 2
    environment = dict(os.environ, PYTHONHASHSEED="0", PYTHONIOENCODING="utf-8")
    failed = 0
    for label, arguments in checks:
        started = time.monotonic()
        try:
            result = subprocess.run([sys.executable, "-X", "utf8", *arguments],
                                    cwd=ROOT, env=environment, capture_output=True,
                                    text=True, encoding="utf-8", errors="replace", timeout=args.timeout)
        except (OSError, subprocess.TimeoutExpired) as error:
            print(f"FAIL {label}: {error}", flush=True)
            failed += 1
            continue
        seconds = time.monotonic() - started
        if result.returncode:
            failed += 1
            print(f"FAIL {label} ({seconds:.1f}s; exit {result.returncode})", flush=True)
            print(result.stdout)
            print(result.stderr)
        else:
            print(f"PASS {label} ({seconds:.1f}s)", flush=True)
    print(f"\n{len(checks) - failed}/{len(checks)} selected checks passed.")
    print("Long searches, optional checks, manuscript compilation, and fresh Lean builds are outside this suite.")
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())
