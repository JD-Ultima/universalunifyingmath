# Geometric Unity — computations

Environment: Python 3.13.17, numpy 2.5.3, sympy 1.14.0 (`pip install numpy==2.5.3 sympy==1.14.0`); `G_symbols.py`/`G_fermion.py` also use scipy 1.18.1 (sparse connected components only).
Randomness: fixed seeds only (`random.seed(20210401)`, `numpy.random.default_rng(20210401)` / `(1)`).
No script reads the source PDF; inputs are transcribed in the code with section/equation references.

| Script | Purpose | Claims | Exact? | Runtime | Output |
|---|---|---|---|---|---|
| `fiber_signature.py` | Signature of the DeWitt family $G_\lambda$ on $\mathrm{Sym}^2$; chimeric signature on $V\oplus H^*$; classification of neutral cases; GL-invariance; Stiefel–Whitney classes of $C$ | GU-010…015, GU-022, GU-023 | exact (Fractions, sympy) | ~3 s | `fiber_signature.out` |
| `spinor_branching.py` | Weight-multiset branching $D_7\to D_2\times D_5\to D_2\times D_3\times D_2\to$ SM; $\nu\oplus\zeta$ content; comparison with the draft's tables (pp. 52–53) | GU-025…027, GU-029, GU-058…065 | exact (integers) | <1 s | `spinor_branching.out` |
| `clifford77.py` | Explicit real $\mathrm{Cl}(7,7)$ on $\mathbb R^{128}$; Majorana–Weyl; invariant forms; grade rules for $\mathfrak{so}(64,64)$, $\mathfrak u(64,64)$; complex structure and commutant on $S^+_\mathbb R$ | GU-016…020 | exact except one SVD (tolerance printed) | ~25 s | `clifford77.out` |
| `gauge_group_identities.py` | Inhomogeneous gauge group $\mathcal G=\mathcal H\ltimes\mathcal N$: group law, actions, $\tau$, stabiliser, $\pi$, $\mu$, augmented torsion (draft Secs. 5–7); sign repair | GU-034…044 | exact (sympy polynomial identities) | ~65 s | `gauge_group_identities.out` |
| `xcheck_dewitt.py` | Independent coordinator cross-check of the fiber/chimeric signature tables | GU-010…015 | numerical (eigvalsh, tol $10^{-9}$) | <1 s | `xcheck_dewitt.out` |
| `xcheck_144.py` | Independent coordinator cross-check: SM content of $144_\pm$ via explicit Cl(10) and Casimirs | GU-058…065 | numerical (eigh, tol $10^{-9}$) | ~2 s | `xcheck_144.out` |
| `D_rep_counts.py` | Exact multiplicities for $\mathrm{Hom}_{\mathrm{Spin}(7,7)}(\Lambda^2\otimes\mathfrak g,\Lambda^{13}\otimes\mathfrak g)$ via Weyl character formula (D1) | GU-045, GU-P09 | exact | ~3 min | `D_rep_counts.out` |
| `D_selfadjoint_shiab.py` | Symbol criterion, self-adjoint subspace, Einstein/Weyl, Helmholtz no-go, true EL equation (D2–D5) | GU-047, 053, GU-P09 | exact (rational) | ~2 min | `D_selfadjoint_shiab.out` |
| `register_stats.py` | Counts of claim types/statuses and gaps quoted in the report | — | exact | <1 s | `register_stats.out` |
| `shiab_checks.py` | Shiab operator family in $d=14$, $\mathfrak u(64,64)\subset \mathrm{Cl}(7,7)\otimes\mathbb C$: degrees, Einstein contraction, Weyl annihilation, tilted equivariance, formal self-adjointness of $\circledast\circ d_B$, cubic term | GU-045, 047, 048, 050, 051, 053 | floating point (complex128), relative errors printed | ~35 s | `shiab_checks.out` |
| `F_jets.py` | Library (no output): exact truncated-Taylor jets over $\mathbb Q(\lambda,c_H,c_V)$; curvature of a metric at a point from its 2-jet; DeWitt metric; Zorro metrics Z1/Z2 on $\mathrm{Met}(X)$ | — | exact | — | — |
| `F_fiber_curvature.py` | Fiber $(\mathrm{Met}(\mathbb R^n),c_VG_\lambda)$: Ricci $=-\tfrac n4G_{1/n}$, $s_F=-n(n-1)(n+2)/(8c_V)$, full curvature $-\tfrac{c_V}4\mathrm{tr}([[H,K],L]M)$ ($\lambda$-free), geodesic symmetry, $\mathbb R\times SL(n)/SO(p,q)$ splitting | GU-010…013 (geometry), GU-G05 | exact (sympy fields; FC3/FC4 for $n\ge3$ at random rational points) | ~60 s | `F_fiber_curvature.out` |
| `F_zorro_reduction.py` | Zorro metric (Z1 tautological / Z2 base) along a section: pull-back, second fundamental form, Ricci blocks and scalar curvature on $s(X)$ for $n=2$ (general metric, symbolic) and $n=3,4,5$ (random metrics, 7 signatures); O'Neill off-section control; 4D cosmological terms $\Lambda_P,\Lambda_A$ vs $\lambda$ | GU-006, 021, 024, 054; GU-G06, G14 | exact | ~27 min | `F_zorro_reduction.out` |
| `F_zorro_einstein.py` | $n=1$ check; Einstein condition on $s(X)$ ($\lambda=1/(n+1)$, Ricci-flat base) incl. off-section failure; Kaluza–Klein isometry with $SL(n+1)/SO$; fiber scale for observed $\Lambda$ | GU-054, GU-G14 | exact except E3 (floating point, inputs printed) | ~40 s | `F_zorro_einstein.out` |
| `G_el_noether.py` | Complete EL system of (9.4) under R-C1′: Stueckelberg factorisation $I=J\circ\pi'$, Noether identity $\Xi=d_A\Upsilon$, true vs printed $\Upsilon$ (exact Euler operator in a non-abelian $d=3$ model, 66 jet directions) | GU-051, 052, 053 | exact (integers / rationals) | ~45 s | `G_el_noether.out` |
| `G_fermion.py` | Fermion operator (9.16): spinor Shiab = Rarita–Schwinger, chirality flip, Grassmann symbol criterion, symmetry of $C\gamma_I$ on $S^\pm$, cross block, Hermitian reading, exact symbol ranks in (7,7) and (14,0) | GU-056, 063, 072 | exact (integers, Gaussian integers) | ~3 min | `G_fermion.out` |
| `G_symbols.py` | Deformation complex (10.1)–(10.9): exact ranks of $\sigma(\circledast\circ d)$ and $\sigma(\delta_2)$ on $\Lambda^1\otimes\mathfrak u(64,64)$ (and $\mathfrak{so}(64,64)$) for $W_A$, $W'_A$, combinations, R-S, at non-null and null covectors, (7,7) and (14,0); monomial Clifford algebra cross-checked against the $128\times128$ matrices | GU-057, 072 | exact (integers) | ~17 min | `G_symbols.out` |

Run all (from this directory; `shiab_checks.py` imports `clifford77.py`):

```bash
for s in fiber_signature spinor_branching clifford77 gauge_group_identities shiab_checks register_stats \
         F_fiber_curvature F_zorro_reduction F_zorro_einstein G_el_noether G_fermion G_symbols; do
  python3 $s.py > $s.out 2>&1
done
```

Notes on reading the outputs

- Lines marked "(expected false)" in `gauge_group_identities.out` are controls; `FAILS` there is the intended outcome. The two unexpected failures are G8 and G9′, the draft's Lemmas 6.5 and 7.2 for the printed sign of $\tau$; see `formalization/geometric-unity/C-gauge-shiab.md`, Prop. C3.
- In `shiab_checks.out` S3/S4, the random-sample values change with the RNG stream. The conclusions do not: the defects are $O(1)$ and there is no constant ratio. Prop. C6 gives an exact counterexample that needs no randomness. It is printed as "explicit counterexample".
- The weight code in `spinor_branching.py` labels Weyl modules by the parity of minus signs in the half-integer weights, so $S^+$ has an even number. The SM embedding is $Y=(B-L)/2+T_{3R}$ with $B-L=-\tfrac23(w_1+w_2+w_3)$, $T_{3L}=(w_4+w_5)/2$ and $T_{3R}=(w_4-w_5)/2$. With these conventions $16_+$ has the SM content of one left-handed generation (check B5).
- `F_*` outputs print `True` for each exact identity checked. In `F_zorro_reduction.out` R3, each line tests the closed forms of Theorems F3 (Z2) and F4 (Z1) of `formalization/geometric-unity/F-zorro-reduction.md` on one random metric; `nabla Ric != 0: True` confirms that the vanishing of the mixed Ricci block is not an artefact of a parallel test metric. `F_*.py` import `F_jets.py`. In `F_zorro_einstein.out` E2 (ii)–(v), `False` is the predicted outcome (Theorem F6: the Einstein condition holds only for Z1, $\lambda=1/(n+1)$, Ricci-flat base, on the section).
- `G_*.out` print `True` for each exact identity checked. In `G_el_noether.out`, the line "Xi = d_A Upsilon^pr" is predicted `False` except for the Chern–Simons control (`M_id`), and "without the D-term" is predicted `False` for `M_sa` (its $Q$ is not $A_0$-parallel; in GU the $\Phi$'s are). In `G_fermion.out` F2, a non-zero "Grassmann" defect for every Shiab is the predicted outcome (Theorem G7). In `G_symbols.out`, `exact: False` everywhere is the result (Theorem G11); RS is implemented as 2·RS (integer coefficients, same rank).
