# Geometric Unity — computations

Environment: Python 3.13.17, numpy 2.5.3, sympy 1.14.0 (`pip install numpy==2.5.3 sympy==1.14.0`).
Randomness: fixed seeds only (`random.seed(20210401)`, `numpy.random.default_rng(20210401)` / `(1)`).
No script reads the source PDF; inputs are transcribed in the code with section/equation references.

| Script | Purpose | Claims | Exact? | Runtime | Output |
|---|---|---|---|---|---|
| `fiber_signature.py` | Signature of the DeWitt family $G_\lambda$ on $\mathrm{Sym}^2$; chimeric signature on $V\oplus H^*$; classification of neutral cases; GL-invariance; Stiefel–Whitney classes of $C$ | GU-010…015, GU-022, GU-023 | exact (Fractions, sympy) | ~3 s | `fiber_signature.out` |
| `spinor_branching.py` | Weight-multiset branching $D_7\to D_2\times D_5\to D_2\times D_3\times D_2\to$ SM; $\nu\oplus\zeta$ content; comparison with the draft's tables (pp. 52–53) | GU-025…027, GU-029, GU-058…065 | exact (integers) | <1 s | `spinor_branching.out` |
| `clifford77.py` | Explicit real $\mathrm{Cl}(7,7)$ on $\mathbb R^{128}$; Majorana–Weyl; invariant forms; grade rules for $\mathfrak{so}(64,64)$, $\mathfrak u(64,64)$; complex structure and commutant on $S^+_\mathbb R$ | GU-016…020 | exact except one SVD (tolerance printed) | ~25 s | `clifford77.out` |
| `gauge_group_identities.py` | Inhomogeneous gauge group $\mathcal G=\mathcal H\ltimes\mathcal N$: group law, actions, $\tau$, stabiliser, $\pi$, $\mu$, augmented torsion (draft Secs. 5–7); sign repair | GU-034…044 | exact (sympy polynomial identities) | ~65 s | `gauge_group_identities.out` |
| `xcheck_dewitt.py` | Independent coordinator cross-check of the fiber/chimeric signature tables | GU-010…015 | numerical (eigvalsh, tol $10^{-9}$) | <1 s | `xcheck_dewitt.out` |
| `register_stats.py` | Counts of claim types/statuses and gaps quoted in the report | — | exact | <1 s | `register_stats.out` |
| `shiab_checks.py` | Shiab operator family in $d=14$, $\mathfrak u(64,64)\subset \mathrm{Cl}(7,7)\otimes\mathbb C$: degrees, Einstein contraction, Weyl annihilation, tilted equivariance, formal self-adjointness of $\circledast\circ d_B$, cubic term | GU-045, 047, 048, 050, 051, 053 | floating point (complex128), relative errors printed | ~35 s | `shiab_checks.out` |

Run all (from this directory; `shiab_checks.py` imports `clifford77.py`):

```bash
for s in fiber_signature spinor_branching clifford77 gauge_group_identities shiab_checks register_stats; do
  python3 $s.py > $s.out 2>&1
done
```

Notes on reading the outputs

- Lines marked "(expected false)" in `gauge_group_identities.out` are controls; `FAILS` there is the intended outcome. The two unexpected failures are G8 and G9′, the draft's Lemmas 6.5 and 7.2 for the printed sign of $\tau$; see `formalization/geometric-unity/C-gauge-shiab.md`, Prop. C3.
- In `shiab_checks.out` S3/S4, the random-sample values change with the RNG stream. The conclusions do not: the defects are $O(1)$ and there is no constant ratio. Prop. C6 gives an exact counterexample that needs no randomness. It is printed as "explicit counterexample".
- The weight code in `spinor_branching.py` labels Weyl modules by the parity of minus signs in the half-integer weights, so $S^+$ has an even number. The SM embedding is $Y=(B-L)/2+T_{3R}$ with $B-L=-\tfrac23(w_1+w_2+w_3)$, $T_{3L}=(w_4+w_5)/2$ and $T_{3R}=(w_4-w_5)/2$. With these conventions $16_+$ has the SM content of one left-handed generation (check B5).
