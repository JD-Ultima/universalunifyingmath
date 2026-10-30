# CT-C013: "the gauge group SU(3) x SU(2) x U(1) ... as theorems" (statements 64.2, 65.1)

Source: compendium 0.1, commit `a28db899`. Abstract `tex/front.tex:25`. Theorem 64.2 "The gauge group as
the face group of the tether", `tex/parts/05-gauge-structure.tex:30-31`; Theorem 65.1 "The gauge group as a
corollary of one construction", lines 54-55, derivation 63-79; genre matrix `tex/parts/01-framework.tex:52-56`;
9.1 (CONDITIONAL, readability) line 344; 9.2 (THEOREM, acting centre) line 350. Computations:
`computation/chain-theory/compendium/gauge_audit.py` (output `outputs/gauge_audit.txt`).

## 1. What the source asserts (a)

64.2: "the face group of the qutrit-qubit tether, the image of the local unitaries U(2) x U(3) on the
link, is the connected determinant slice S(U(2) x U(3)) with Lie algebra of dimension twelve, derived
algebra of dimension eleven and a one-dimensional centre ..., isomorphic to (SU(2) x SU(3) x U(1))/Z6 by
an explicit covering whose kernel is the cyclic group of order six; the uniqueness cell at the dimension
twelve and the rank four picks the qutrit-qubit type".
65.1: applying "the stabiliser of the Hermitian unit inside the isometries of the split" at three end kinds
gives su(2) (dim 3, rank 1, centre Z2) at "hub-leg qubit copies", su(3) (dim 8, rank 2, centre Z3) at
"qutrit ends open at the hub leg", and u(1) at "the dock qubit ends" with slice representative
(3,3,-2,-2,-2); "the kernel of the covering is re-met ... the cyclic group of order six".
Both graded THEOREM. Direct premises of 64.2 include OPEN 2.4, 2.5, 2.7, 2.8, 9.8, 18.5, 18.6 and
CONDITIONAL 2.2, 9.1, 9.7, 18.1, 18.3, 18.4; of 65.1, OPEN 9.8, D.1, D.2, D.5-D.7 and CONDITIONAL 9.1, 9.7.

## 2. Definitions (b)

For a genre (d1, d2): the face algebra f(d1,d2) = {X (x) 1 + 1 (x) Y : X in u(d1), Y in u(d2)} inside
u(d1 d2); the determinant slice S(U(d1) x U(d2)) = {(A,B) : det A det B = 1}; its centre generator
z = (i d2 1_{d1}, -i d1 1_{d2}). "Readable" (9.1): d1 d2 <= 6 ("the complete-read bound").
"Acting centre" (9.2): z acts on C^{d1} (x) C^{d2} as a nonzero scalar.

## 3. Results

**Proposition G1 (the group; proven, c; standard).** The map S(U(2) x U(3)) -> U(6), (A,B) -> A (x) B,
is injective, and its image equals the image of U(2) x U(3). Proof: (lambda 1, mu 1) acts trivially iff
lambda mu = 1; in the slice lambda^2 mu^3 = 1, so lambda = mu = 1. Given (A,B) in U(2) x U(3), put
c = det A det B; then (c A, c^-1 B) lies in the slice and has the same image. The covering
SU(2) x SU(3) x U(1) -> S(U(2) x U(3)), (A,B,z) -> (z^3 A, z^-2 B), has kernel {z : z^6 = 1} = Z6 (all six
elements verified exactly). dim f(3,2) = 12, derived algebra dim 11, centre dim 1 (computed by exact rank).
This is the textbook isomorphism S(U(2) x U(3)) = (SU(2) x SU(3) x U(1))/Z6 (e.g. Baez-Huerta 2010,
which the source cites). Novelty: none.

**Proposition G2 (all six genres; proven, c).** dim f = d1^2 + d2^2 - 1, derived algebra dim
d1^2 + d2^2 - 2, and z acts on the tensor product as i(d2 - d1) 1. Hence the kernel of
S(U(d1) x U(d2)) -> U(d1 d2) is U(1) for d1 = d2, Z_{|d1-d2|} otherwise (trivial for |d1-d2| = 1), and
|ker(SU x SU x U(1) -> S(U x U))| = lcm(d1,d2). The printed genre matrix (dimension; central scalar) =
(7;0), (12;1), (19;2), (17;0), (24;1), (31;0) is reproduced exactly as (dim f; |d1 - d2|).

**Proposition G3 (the selection of (3,2); c, with its dependence).** Among the six genres:
- dim 12 and rank 4 select (3,2) only. These targets are the Standard Model's; the criterion is a match,
  not a derivation.
- "readable and acting centre" (d1 d2 <= 6 and d1 != d2) selects (3,2) only. Sensitivity: cap 4 selects
  none; cap 8 or 9 selects (3,2) and (4,2); cap 12 selects (3,2), (4,2), (4,3). Enlarging the alphabet
  to {2..8} at cap 6 does not change the selection. Uniqueness therefore rests on the value 6 of the read
  cap, which enters through 9.1, graded CONDITIONAL.

**Proposition G4 (the u(1) and the global form; c, e).** diag(3,3,-2,-2,-2) is traceless, commutes with
su(2) + su(3), and equals 6Y on the 5 = (1,2)_{1/2} + (3,1)_{-1/3}. Every SM field (Q_L, u_R, d_R, L, e_R,
H, gauge bosons) satisfies t/3 + d/2 + Y in Z (t triality, d duality), so the Z6 quotient acts
faithfully. The quotients by Z1, Z2, Z3 are equally consistent with all observed fields (Tong 2017):
the Z6 global form is a consistent but presently untested output.
*Failure condition:* observation of any state with t/3 + d/2 + Y not in Z (for instance a colour-singlet,
weak-singlet particle of charge e/2) falsifies the Z6 form.

**Proposition G5 (what is not derived; a).** The theorems give a group of local unitaries on one link.
They do not give: gauge fields or locality in spacetime (the source names gauge invariance and the group
"by correspondence"); the matter representations (66.1-66.3, 66.5 OPEN); the hypercharge values
(163.x CONDITIONAL/OPEN); the assignment of the three factors to three end kinds (a choice in 65.1).
Tension: the abstract's correspondence assigns (3,2) to electromagnetism and (3,3), (4,2) to the strong
and weak interactions, while the full SU(3) x SU(2) x U(1) is read at (3,2).

## 4. Grade consistency

The front matter states that "a composed claim carr[ies] the lowest status in its chain". 64.2 and 65.1
list OPEN statements among their premises and are graded THEOREM. Either the premise lists are not
logical premises (then the dependency data overstate dependence) or the grade violates the stated rule.

## 5. Choices vs structure matched

Discrete choices: (1) alphabet {2,3,4} (inherited from the 64 listed webs, 9.1); (2) read cap 6 (9.1);
(3) acting-centre criterion d1 != d2 (9.2); (4)-(6) assignment of SU(2), SU(3), U(1) to hub-leg qubit
copies, open qutrit ends and dock qubit slice (65.1). Six choices. Output matched: the Lie algebra
su(3) + su(2) + u(1) and the Z6 centre, which are forced once (3,2) is chosen; no further number.

## 6. Lean (task 4)

64.1 and 64.2 name `Foundations.Ch02.row_521259_GIKC` and `StandardModel.Ch10.row_510060_GGRL`, both
`opaque ... : Prop` (no content). 65.1 names `Foundations.Ch05.row_530016_GSTB`, which restates
`D530016_GSTB.checks` (an 8-row key/value `String` table; the theorem proves that `cell.C-EQGS` has value
`"YES"`, `eqgs.letter` has value `"EQUAL"`, etc.), and `row_530426_UB4H : cond_530426_UB4H -> ...`, the same
for another 8-row string table under an opaque hypothesis. No group, matrix, Lie algebra or kernel is
represented in Lean. G1-G3 are, however, readily formalizable (Mathlib has `Matrix.unitaryGroup`,
`Matrix.specialUnitaryGroup`, Kronecker products and `LieAlgebra`).

## 7. Repo-side Lean for G1-G3 (2026-10-30)

The proofs are in `proofs/lean/chain-theory/GaugeSelection.lean`. It is checked with Lean 4.35.0, core library only, by `run_gauge.sh`; the output is in `GaugeSelection.out`. There are 71 theorems: 36 use no axioms, 7 `propext` only, 28 `propext` and `Quot.sound`. None uses `Classical.choice` or `sorryAx`. README.md lists what is and is not certified. Results beyond §3:

- **The cap enters as an interval.** "Readable and acting centre" selects exactly $\{(3,2)\}$ iff the cap $c$ satisfies $6\le c\le7$ (`unique_iff_cap`). This holds over $\{2,3,4\}$ and over the alphabet of all integers $\ge2$ (`unique_iff_cap_unbounded`), hence over every alphabet of integers $\ge2$ that contains 2, 3 and 4. Alphabets without 4 give other intervals: every cap $\ge6$ over $\{2,3\}$, caps 6-9 over $\{2,3,5\}$ (`unique_iff_cap_23`, `unique_iff_cap_235`). Caps $\le5$ select nothing, caps 8-11 select $(3,2),(4,2)$, and caps $\ge12$ add $(4,3)$ (`cap_profile`). G3's sensitivity values are special cases.
- **At cap 6 the alphabet enters only through its lower end and through $2,3\in A$.** No upper end is needed (`selection_any_alphabet`); the cap interval above does depend on $4\in A$. Admitting $d=1$ adds $(2,1),(3,1),(4,1)$, and with no upper end also $(5,1),(6,1)$.
  - The alphabet input is 9.1's "admitted alphabet". The box 2.7 "The act alphabet" (OPEN) is a different object, the two move letters attach and join.
- **The acting centre is the only input that separates $(3,2)$ from $(2,2)$.** Under the conjugate convention of F-CT-3 it holds at every genre, and the selection becomes $(2,2),(3,2)$ (`conj_selection`).
- **"Dim 12 and rank 4": each condition alone suffices.** Each selects $(3,2)$ among the six genres. Dim 12 alone does so among all pairs of naturals, since $13=3^2+2^2$ is the only representation (`dim12_unbounded`). Rank 4 alone does so for $d_2\ge2$. The two targets are $\dim$ and rank of $su(3)\oplus su(2)\oplus u(1)$ (`targets_are_SM`). G3's "a match, not a derivation" is unchanged.
- **Normalisation.** §2's generator $z=(i d_2 1,-i d_1 1)$ is the normalisation of the printed scalars $|d_2-d_1|$. It is not primitive when $\gcd(d_1,d_2)>1$: the primitive generator gives 1 at $(4,2)$, not 2 (`primitive_scalar`).

  G2's $|\ker|=\mathrm{lcm}(d_1,d_2)$ holds for the primitive covering $(A,B,z)\mapsto(z^{d_2/g}A,z^{-d_1/g}B)$. Exponentiating §2's $z$ instead gives kernel order $d_1d_2$: 8 rather than 4 at $(4,2)$ (`cover_kernel_normalisation`, `cover_kernel_scalar`). At $(3,2)$ the two coincide (6).
- **G1 algebra over any commutative ring:**
  - kernel of $(X,Y)\mapsto X\otimes1+1\otimes Y$ equal to $\{(c1,-c1)\}$, and kernel of $\otimes$ equal to the scalar pairs with $\lambda\mu=1$, for every genre;
  - trivial kernel on the $(2,3)$ slice;
  - the rescaling $(cA,c^{-1}B)$;
  - the covering kernel $\{z:z^6=1\}$, with $(A,B)=(z^31,z^21)$;
  - the six triples $k\cdot(3,2,1)$ of 65.1 in the exponent model.

**Not certified.** These remain paper steps or are outside core: ℂ, unitarity, the Lie groups, the Lie-group isomorphism with $(SU(2)\times SU(3)\times U(1))/\mathbb Z_6$, surjectivity of the covering, the mixed-product property, rank-nullity, the derived algebra and the centre dimension. §6's Mathlib route remains the way to state these.
