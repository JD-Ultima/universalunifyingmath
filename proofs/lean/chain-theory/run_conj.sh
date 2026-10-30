#!/usr/bin/env bash
# Usage: LEAN=/path/to/lean-4.35.0/bin/lean ./run_conj.sh > Conjugation.out
# Standalone: checks Conjugation.lean with Lean 4.35.0, core library only (no upstream file).
set -euo pipefail
"$LEAN" --version | grep -q "version 4.35.0" || { echo "expected Lean 4.35.0" >&2; exit 1; }
if grep -n "sorry" "$(dirname "$0")/Conjugation.lean"; then echo "sorry found" >&2; exit 1; fi
"$LEAN" "$(dirname "$0")/Conjugation.lean"
echo "exit 0"
