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
- `weight_is_class_function` — Prop. 2 of the audit (F-CT-2): on the object, the upstream `kOf` equals a function of (commutator?, `parity` under `exchAssign`, `sameEnd`) only. Proved by `decide +kernel` over all 172 records.
- `ReviewCheck.third_config_under_conjugate_convention` (referee addition, review CT-4) — the upstream `readable` is $d_1,d_2\le3$, not the read cap. So under the conjugate convention the (3,3) configuration 3R 3R 3U also satisfies the hypotheses of `configuration_unique`. The `#eval`s give the readouts of the all-qubit and all-qutrit worlds: $(30+30\sqrt2)/2304$, and $(184+200\sqrt3)/419904$.
