#!/usr/bin/env bash
# Usage: LEAN=/path/to/lean-4.35.0/bin/lean UPSTREAM=/path/to/chain-theory ./mutation_spotcheck.sh > mutation_spotcheck.out
# Lean check of the mechanical mutation space of computation/chain-theory/alpha/mechanical_mutations.py.
# Every single mutant and a seeded sample of 300 pairs (seed 20261030) is written as Parts 0-1 of the pinned
# ChainAlpha.lean with the mutation's text edits, followed by `#eval readout objectWorld objectTwoPoint`. Lean 4.35.0
# checks each file (2 at a time, 120 s timeout); its result is compared with the Python transcription's prediction
# (a readout, "fails to elaborate", or "infeasible": a timeout or resource abort). The upstream file is never
# modified. Also prints the Lean integer conventions (division, remainder, subtraction) that the transcription follows.
set -euo pipefail
"$LEAN" --version | grep -q "version 4.35.0" || { echo "expected Lean 4.35.0" >&2; exit 1; }
up="$UPSTREAM/papers/alpha/lean/ChainAlpha.lean"
echo "522844e911cf5f6091b603db0d1d042f50544d2ec8d10ed9336b722945c331da  $up" | sha256sum -c - | sed "s|$UPSTREAM|\$UPSTREAM|"
conv=$(mktemp -d)
cat > "$conv/Conventions.lean" <<'LEAN'
#eval ((-7 : Int) / 2, (7 : Int) / (-2), (-7 : Int) / (-2), (5 : Int) / 0, (-7 : Int) % 2, (5 : Int) % 0)
#eval ((5 : Nat) / 0, (5 : Nat) % 0, (3 : Nat) - 5)
LEAN
echo "Lean arithmetic the transcription follows: Int (-7)/2, 7/(-2), (-7)/(-2), 5/0, (-7)%2, 5%0; Nat 5/0, 5%0, 3-5:"
"$LEAN" "$conv/Conventions.lean"
here=$(cd "$(dirname "$0")" && pwd)
py="$here/../../../computation/chain-theory/alpha/mechanical_mutations.py"
work=$(mktemp -d)
python3 "$py" --lean "$up" --emit-lean "$work" >&2
find "$work" -name '*.lean' | sort | xargs -P 2 -I{} sh -c \
  'f="$1"; timeout 120 "$0" "$f" > "${f%.lean}.out" 2>&1; echo "exit $?" >> "${f%.lean}.out"' "$LEAN" {}
python3 "$py" --compare "$work"
