#!/usr/bin/env bash
# Usage: LEAN=/path/to/lean-4.35.0/bin/lean ./run_gauge.sh > GaugeSelection.out
# Standalone: checks GaugeSelection.lean with Lean 4.35.0, core library only (no upstream file).
set -euo pipefail
"$LEAN" --version | grep -q "version 4.35.0" || { echo "expected Lean 4.35.0" >&2; exit 1; }
"$LEAN" "$(dirname "$0")/GaugeSelection.lean"
