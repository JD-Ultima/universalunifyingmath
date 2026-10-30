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

## `GaugeSelection.lean`: gauge-group selection (compendium 64.2, 65.1)

This file formalizes G1-G3 of `formalization/chain-theory/compendium/03-gauge-group.md`. It is standalone: it does not read the upstream file and has no imports. It is checked with Lean 4.35.0, core library only.

Run:

```bash
LEAN=/path/to/lean-4.35.0-linux/bin/lean ./run_gauge.sh > GaugeSelection.out
```

- Wall time is about 10 s.
- The output is in `GaugeSelection.out`: exit 0, no warnings.
- 66 theorems, each followed by `#print axioms`:
  - 33 use no axioms;
  - 7 use `propext` only;
  - 26 use `propext` and `Quot.sound`.
- No theorem uses `Classical.choice` or `sorryAx`.
- The two `#eval` tables in the output (the selection for caps 0-16, and the six kernel triples) are for reading only; the theorems are the certificates.

**What is certified.**

- **G3, the selection.** The inputs are the read cap $d_1d_2\le c$, the acting centre $d_1\ne d_2$, and the alphabet $\{2,3,4\}$, which is 9.1's "admitted alphabet". It is not 2.7 "The act alphabet", which is about the two move letters.
  - `selection`: the three inputs together select $(3,2)$. `readable_types`: the cap alone selects $(2,2),(3,2)$; this is also the result of dropping the acting centre. `drop_cap`: the acting centre alone selects $(3,2),(4,2),(4,3)$.
  - `cap_profile` holds for every cap $c$: $(3,2)$ enters at 6, $(4,2)$ at 8 and $(4,3)$ at 12. `unique_iff_cap`: the selection is exactly $\{(3,2)\}$ iff $6\le c\le7$. `unique_iff_cap_unbounded` proves the same for every alphabet of integers $\ge2$. `cap_4`, `cap_8`, `cap_9` and `cap_12` give the sensitivity values in the note.
  - The upper end of the alphabet does not matter. `alphabet_2_to_8` gives 28 genres, with the same selection. `selection_any_alphabet` holds for any alphabet of integers $\ge2$ that contains 2 and 3. The lower end does matter: admitting $d=1$ adds $(2,1),(3,1),(4,1)$, and, with no upper end, also $(5,1),(6,1)$ (`alphabet_with_one`, `cap6_acting_with_one`).
  - Under the conjugate convention of F-CT-3, the acting criterion holds at every genre (`actingConj_always`). The selection is then $(2,2),(3,2)$ (`conj_selection`).
  - Under the upstream readability predicate (both ends $\le3$), $(3,2)$ is still the only selection (`upstream_selection`). Under the conjugate convention with that predicate, $(2,2),(3,2),(3,3)$ survive (`upstream_conj_selection`; cf. `ReviewCheck`).
- **G2 and G3, the formulas.** The formulas are $\dim f=d_1^2+d_2^2-1$, derived dimension $=\dim f-1$, rank $=d_1+d_2-1$, and central scalar $d_2-d_1$.
  - `genre_matrix` reproduces the printed matrix $(7;0),(12;1),(19;2),(17;0),(24;1),(31;0)$.
  - In the "dim 12, rank 4" cell, `dim12_rank4` shows that the two conditions select $(3,2)$. `dim12_alone` and `rank4_alone` show that each condition alone also does. `dim12_unbounded` holds over all pairs of naturals ($13=3^2+2^2$ only). `rank4_unbounded` holds for $d_2\ge2$.
  - `targets_are_SM`: 12 and 4 are the dimension and rank of $su(3)\oplus su(2)\oplus u(1)$.
  - `primitive_scalar`: the printed scalar $|d_2-d_1|$ uses the generator $(d_2,-d_1)$, which is not primitive when $\gcd>1$. The primitive generator gives 1 at $(4,2)$, not 2. Only zero versus nonzero enters the selection.
- **G1 and G2 at the matrix level.** Matrices are `Fin n → Fin n → α` over any `α` with a core `Lean.Grind.CommRing` structure. Determinants are the Leibniz formulas for $n=2,3$.
  - `face_kernel`, for every genre: $X\otimes1+1\otimes Y=0$ iff $(X,Y)=(c1,-c1)$.
  - `centre_acts`: the generator $(a1,-b1)$ acts on the link as $(a-b)\,1$. `centre_traceless`: it is traceless on $C^{d_1}\oplus C^{d_2}$.
  - `kron_kernel`, for every genre: $A\otimes B=1$ iff $A=\lambda1$, $B=\mu1$ with $\lambda\mu=1$.
  - `tensor_kernel_trivial_23`: on $\det A\det B=1$, that kernel is trivial for $(2,3)$.
  - `rescale_into_slice`: $(cA,c^{-1}B)$ with $c=\det A\det B$ lies in the slice and has the same image.
  - The covering $(A,B,z)\mapsto(z^3A,z^{-2}B)$ maps into the slice (`cover_into_slice`). `cover_kernel_23`: $(A,B,z)$ with $\det A=\det B=1$ lies in its kernel iff $z^6=1$ and $(A,B)=(z^31,z^21)$.
  - For every genre at scalar level: `tensor_kernel_scalar` says the slice kernel is $\lambda^{d_1-d_2}=1$. `cover_kernel_scalar` says the covering with exponents $(e_1,e_2)$, $e_2d_1=e_1d_2=L$, has kernel $z^L=1$.
  - `slice_representative`: 65.1's $(3,3,-2,-2,-2)$ is the generator at $(2,3)$, and it is traceless.
- **Finite exponent models.** $\mu_M$ is modelled as $\mathbb Z/M$.
  - `six_triples`: the kernel of the $(2,3)$ covering consists of the six triples $k\cdot(3,2,1)$ mod 6, as 65.1 states. `six_triples_in_mu12`: it still has six elements inside $\mu_{12}$.
  - `cover_kernel_orders`: the kernel orders are $\mathrm{lcm}(d_1,d_2)=2,6,4,3,12,4$ for the six primitive coverings. `cover_kernel_normalisation`: the unreduced covering at $(4,2)$ has kernel order 8.
  - `tensor_kernel_orders`: the tensor kernels have $|d_1-d_2|$ elements, and the whole $\mu_{12}$ when $d_1=d_2$.

**What is not certified.**

- **ℂ and unitarity.** ℂ is not in core, so the instantiation of §3 at ℂ is not checked. Unitarity is never used. The matrix theorems hold over every commutative ring. The file instantiates them only at `Int`, as a non-vacuity check (`kron_kernel_int_instance`, `minus_one_not_in_slice`).
- **Lie-group statements.** None of these is expressed:
  - $U(n)$, $SU(n)$ and $S(U(2)\times U(3))$ as Lie groups;
  - the isomorphism $S(U(2)\times U(3))\cong(SU(2)\times SU(3)\times U(1))/\mathbb Z_6$;
  - surjectivity of the covering, which needs sixth roots in $U(1)$;
  - that $\mu_N\subset U(1)$ has $N$ elements; §4 uses exponent models.
- **Paper steps.** These are the standard steps between what is proved and the Lie-algebra numbers:
  - the mixed-product property $(A\otimes B)(C\otimes D)=AC\otimes BD$, which turns "trivial kernel" into "injective";
  - rank-nullity and $\dim_\mathbb R u(d)=d^2$, which turn `face_kernel` into $\dim f=d_1^2+d_2^2-1$;
  - that the derived algebra is $su(d_1)\oplus su(d_2)$, and that the centre is one-dimensional;
  - $\det(\lambda1_d)=\lambda^d$ for $d\ne2,3$.

  For $(2,3)$ the upstream alpha Lean checks the ranks 12 and 11 with its own `rank` function: `ChainAlpha.generators_independent` and `ChainAlpha.center_one_dimensional`, built by `run.sh`.
- **Inputs and physics.** The read cap itself (9.1, CONDITIONAL) is a parameter here; the file shows what each value selects. That the readability and acting-centre criteria are the right ones is not a mathematical statement. The physics identifications are not addressed (G5).

`ChainAlpha.selection_under_read_cap` in the upstream alpha Lean already proves the core arithmetic: under the cap, as a hypothesis, $d_1,d_2\ge2$ with $d_1\ne d_2$ gives $(2,3)$ or $(3,2)$. This file adds the dependence on each input and the G1/G2 algebra.
