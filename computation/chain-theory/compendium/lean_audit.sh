#!/usr/bin/env bash
# Round 6: build the compendium's Lean with core Lean 4.35.0 in a scratch directory, print the axioms of every
# tied declaration, run our own checks, and resolve the restated / kernel declarations.
#
# Usage: LEAN=/path/to/lean-4.35.0/bin/lean SRC=/path/to/chain-theory/compendium WORK=/fresh/scratch/dir ./lean_audit.sh
#   SRC is the compendium directory at commit a28db899 (read as data: no lake, lakefile or script of it is run;
#   lean_build.py copies the library sources into WORK and compiles the copies). WORK must not exist.
#   Build B needs about 16 GB of memory at --jobs 3 (one AlphaDerivation module peaks at 7.7 GB under 4.35.0).
#   To continue an interrupted Build B, rerun its line with --resume.
set -euo pipefail
cd "$(dirname "$0")"
"$LEAN" --version | grep -q "version 4.35.0" || { echo "expected Lean 4.35.0" >&2; exit 1; }
mkdir -p "$WORK" outputs
# Build A: the published sources, unchanged, AxiomDeclaration library (the first module that fails under 4.35.0)
python3 -I lean_build.py --lean "$LEAN" --src "$SRC/lean" --work "$WORK/A" --targets AxiomDeclaration \
    --record outputs/lean_build_A.tsv > outputs/lean_build_A.txt || true
# Build B: all 493 modules of the source's `lake build`, with the five proofs of lean_build_patch.tsv re-tacticked
python3 -I lean_build.py --lean "$LEAN" --src "$SRC/lean" --work "$WORK/B" --patch lean_build_patch.tsv \
    --jobs 3 --mem-mb 13500 --mem-scale 1.35 --record outputs/lean_build_B.tsv > outputs/lean_build_B.txt
# why the five proofs fail under 4.35.0 (decide finds no instance; rfl proves the same statements)
LEAN_PATH="$WORK/B/out" "$LEAN" lean_audit/Diagnose.lean > outputs/lean_diagnose.txt
# #print axioms for every tied name (and the source's own AxiomsCheck lists), in a file of ours
python3 -I lean_resolve.py "$SRC" --emit-axioms-file "$WORK/AuditAxioms.lean"
LEAN_PATH="$WORK/B/out" "$LEAN" "$WORK/AuditAxioms.lean" > outputs/lean_axioms.txt
# where Classical.choice enters under 4.35.0 (a meta-program of ours; `import Lean` and `#eval` in our own file)
LEAN_PATH="$WORK/B/out" "$LEAN" lean_audit/AxiomTrace.lean > outputs/lean_axiom_trace.txt
# our own checks against the source's definitions
LEAN_PATH="$WORK/B/out" "$LEAN" lean_audit/ResolveChecks.lean > outputs/lean_resolve_checks.txt
python3 -I lean_resolve.py "$SRC" --axioms outputs/lean_axioms.txt --checks outputs/lean_resolve_checks.txt \
    --build-record outputs/lean_build_B.tsv --build-record-a outputs/lean_build_A.tsv \
    --tsv outputs/lean_resolve.tsv > outputs/lean_resolve.txt
