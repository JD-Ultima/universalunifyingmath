# Lean checks — Chain Theory

`ConventionCheck.lean` is appended to the author's `ChainAlpha.lean`, which is pinned by sha256 and not modified. The combined file is checked with Lean 4.35.0, core library only.

Run:

```bash
LEAN=/path/to/lean-4.35.0-linux/bin/lean UPSTREAM=/path/to/chain-theory ./run.sh
```

- `UPSTREAM` is a clone of `github.com/alexwg/chain-theory` at commit `a28db89`.
- Lean binary: `lean-4.35.0-linux.tar.zst` from the leanprover/lean4 GitHub release, sha256 `7289597a44543246c9227b8cae501767e095ce17f6d448b95fc4a8dc97268f6d`.
- Wall time is about 3 minutes on 4 cores. The output is in `run.out`.

Theorems added (see `formalization/chain-theory/alpha/audit.md`, F-CT-3):

- `actA_ne_zero_iff` — the upstream `actingCharge` is exactly "the primitive central element acts nontrivially under the unconjugated convention".
- `uniqueness_fails_under_conjugate_convention` — under the conjugate convention, a second configuration (2R 2R 2U at a shared qubit) satisfies all hypotheses of upstream `configuration_unique`.
- `upstream_excludes_qubitConfig` — the upstream convention excludes it.
