# Howard — computations

**Environment:** Python 3.13.17 with sympy 1.14.0, mpmath 1.3.0, numpy 2.5.3, z3-solver 5.3.3.0 and mendeleev 1.3.0. Install:

```bash
pip install sympy==1.14.0 mpmath==1.3.0 numpy==2.5.3 z3-solver==5.3.3.0 mendeleev==1.3.0
```

The scripts use no randomness, except one Monte-Carlo check in `lynchpin_geometry.py`, which uses numpy `default_rng(0)`. They read no downloaded files. Every number they use is typed in from the quotes listed in `claims/howard.md`.

| Script | Purpose (claim IDs) | Runtime | Output |
|---|---|---|---|
| `arith_candidates.py` | Eight candidate operations $\otimes$ tested over $\mathbb R$ and $\mathbb Z$ against H1, H2, LD, RD, AS, CM, identity, INC and BAL. z3 decides each property, giving a proof or a counterexample. Also gives $\otimes$-square roots (TH-001–004, 007) | ~2 s | `arith_candidates.out` |
| `axiom_lattice.py` | z3 ground certificates (UNSAT of a subset's own ground instances) for nine cores; all 128 subsets of the seven rules H1, H2, LD, RD, AS, CM, INC classified, with counts per certificate type (46 z3, 12 Theorem A6 only, 70 witnesses); BAL as an optional eighth rule (Part 2b); models on $\mathbb Z/n$; brute force over all tables on $\mathbb Z/2$ and $\mathbb Z/3$ | <1 s | `axiom_lattice.out` |
| `semirings.py` | Log₂ and tropical semirings that realise H1 and H2 with distributivity (`arithmetic.md` §5) | <1 s | `semirings.out` |
| `numerics.py` | √2 loop, "most numbers", speed of light, element "tones" and element data, the 64-circle count, 3-smoothness (TH-005, 006, 008–010; P-TH-1) | ~1 s | `numerics.out` |
| `lynchpin_geometry.py` | T-cone geometry, pentagon mismatch, an explicit pentagonal Lynchpin, ball interstices, cluster (V,E,F) table (TH-020, 022, 024) | ~12 s | `lynchpin_geometry.out` |

Run all of them from this directory:

```bash
for s in arith_candidates axiom_lattice semirings numerics lynchpin_geometry; do python3 -I $s.py > $s.out 2>&1; done
```

## Lean companion

`proofs/lean/howard/HowardArith.lean` uses Lean 4 core only (no Mathlib, no imports). The current file was checked with **Lean 4.35.0** (commit `eb6cafb7`) on linux x86_64: exit 0, about 1.5 s. The theorems that predate the 2026-10-30 revision were also checked with **Lean 4.24.0** (commit `797c613e`). The additions (generic one-sided no-go, ℚ instances via core `Rat`, A6) were checked with 4.35.0 only.

```bash
lean proofs/lean/howard/HowardArith.lean   # from the repo root
```

The expected output is the `#print axioms` lines recorded in `proofs/lean/howard/HowardArith.lean.out`, one for each of the file's 80 theorems. It contains no errors and no `sorry`. The generic theorems use no axioms. The ℤ, ℤ/4 and Euler theorems use at most `propext` and `Quot.sound`. The ℚ instances also use `Classical.choice`, which core's own `Rat` lemmas (e.g. `Rat.le_trans`, `Rat.natCast_add`) depend on.

Mathlib's cache hosts were blocked in this sandbox, and so was `release.lean-lang.org`. The toolchain was therefore installed from the GitHub release tarball `lean-4.24.0-linux.tar.zst` and registered with `elan toolchain link`. The 4.35.0 check used the core `lean` binary from the `lean-4.35.0-linux` release archive.

The repository README lists `proofs/lean/` as the home for Lean files. This agent's write scope was `formalization/howard/`, so the coordinator may move the file there.
