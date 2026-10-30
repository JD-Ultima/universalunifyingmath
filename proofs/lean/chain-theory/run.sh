#!/usr/bin/env bash
# Usage: LEAN=/path/to/lean-4.35.0/bin/lean UPSTREAM=/path/to/chain-theory ./run.sh
set -euo pipefail
up="$UPSTREAM/papers/alpha/lean/ChainAlpha.lean"
echo "522844e911cf5f6091b603db0d1d042f50544d2ec8d10ed9336b722945c331da  $up" | sha256sum -c -
work=$(mktemp -d)
cat "$up" "$(dirname "$0")/ConventionCheck.lean" > "$work/Combined.lean"
"$LEAN" "$work/Combined.lean"
