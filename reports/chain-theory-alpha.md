# Chain Theory — α(0) Letter: Phase 1–2 summary

**Source:** Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, v1.0.0, 2026-10-23, with Supplement and Lean file (`sources/chain-theory.md`, commit `a28db89`).
**Details:** `formalization/chain-theory/alpha/audit.md`.
**Registers:** `claims/chain-theory-alpha.md` (CT-001…023) and `gaps/chain-theory-alpha.md` (CT-G01…G10).

## What the Letter claims

1. One axiom (pure maximally entangled pairs of qudits, known only by projection onto the pair state) plus five definitions select a unique configuration: a qutrit shared by three qubit links, two of them "read". This is Theorem 1.
2. A weighted count over the 172 "measurement records" of that configuration gives $(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$, with $1/x = 137.035\,999\,189\,566\ldots$. This is Theorem 2.
3. This number is conjectured to equal $\alpha(0)$ (the Photon Conjecture). It lies $0.60\sigma$ from CODATA 2022.

## What we established

| Item | Result | Grade |
|---|---|---|
| Counting and value | Reproduced exactly by an independent $36\times36$-operator implementation. Agrees with the author's `records.csv` class by class. 60 digits agree. | **Verified** |
| Author's Lean file | Builds under Lean 4.35.0 in 170 s, axioms `propext`/`Quot.sound`/`Classical.choice` only. `main` is an equality of coefficient lists, i.e. arithmetic of the definitions. | **Verified** |
| Two exchanges | $X$ and $Y$ are the same involution $E$ on all product records. The weight law reduces to a class function: odd → 1, even mixed → $\sqrt2$, even same-end → $\sqrt3$, commutator → $\sqrt6$. | **Proven** (Prop. 2) |
| Selection of the (2,3) link | "Acting central charge" holds only for $d_1\ne d_2$ under the Letter's convention (A). Under the conjugate far-end convention (A*), $(2,2)$ also qualifies. Under a matched-direction reading (B), no type qualifies. Lean: `uniqueness_fails_under_conjugate_convention`. | **Proven**; uniqueness is convention-dependent |
| Hypercharge of the link | $Y=1/6$ ($Q_L$) under A; $Y=5/6$ (SU(5) leptoquark block) under A*. | **Proven** |
| Alternatives | All 13 closed-form alternatives in the author's table are reproduced. Four normalization alternatives the table omits give $1/x \in\{205.6, 60.9, 129.5, 188.9\}$. | **Computed** |
| Trial factor | In a declared 46 080-rule space containing every published alternative, the expected number of chance hits within $1\sigma$ is $4.4\times10^{-6}$; a chance hit needs ~$10^{10}$ rules. In $\tfrac1{31104}\mathbb Z_{\ge0}[\sqrt2,\sqrt3]$ with coefficients ≤200 there are $1.6\times10^9$ elements; exactly one is within $1\sigma$ (the author's), against 0.024 expected. | **Computed** |
| Covariance | The Letter's plain exchange does not commute with the symmetry group of its own joint state, and its record span is not invariant. The covariant exchange (with transpose) restores invariance and gives identical counts and value. | **Proven** (F-CT-10) |
| Freedom left by symmetry | Under $SU(3)\times SU(2)\times\langle E'\rangle$ the record space has 10 isotypic components (commutant 22). An invariant dressing has a 10-parameter trace, and the Letter's weights are one point of that family. | **Computed** (F-CT-10) |
| Configurations | 76 / 42 / 23 reproduced. Selected: 1 configuration under A, 2 under A* (the competitor has $1/\text{value}=31.81$), none under B. | **Computed** |
| Read-cap lemmas | Choi's map is positive and non-decomposable by SDP, including after embedding. | **Cross-checked** |
| Algebraic data | $\alpha_\text{cand}$ has degree 4. The minimal polynomial and Galois conjugates are given. | **Proven** |

## Assessment, by layer

- **(c) Proven.** The arithmetic and the selection theorem follow from the stated definitions under convention A.
- **(b) Reconstruction gaps.** The number rests on choices the axiom does not fix:
  - a conjugation convention and a direct-sum trace (CT-G03/G04);
  - three wording clauses assigning $\sqrt3$ vs $\sqrt2$ (CT-G02);
  - amplitude rather than Born weights (CT-G08);
  - the definition of $N$ (CT-G05);
  - an exception admitting single-slot commutator records (CT-G09).

  None of these is shown to be forced. Their published alternatives spread $1/x$ over $[97.5,\,180.8]$.
- **(d) Conjectured.** The Photon Conjecture. Physically it lacks:
  - a map from the hypercharge readout to $\alpha_\text{em}$, given that $e=g'\cos\theta_W$ (CT-G06);
  - a propagator or momentum dependence for "on-shell" (CT-G07).
- **(e) Testable.** The candidate fixes every digit.
  - It is refuted if consistent determinations of $\alpha(0)$ with $u\lesssim2.5\times10^{-9}$ (8× better than CODATA 2022) exclude it at $5\sigma$.
  - Resolving the Rb ($-1.49\sigma$) vs electron $g-2$ ($+1.57\sigma$) tension is the nearer test.

**Evidential weight of the 10-digit match.** If the rules were fixed blind to $\alpha$, the agreement has $p\sim10^{-6}$ against the declared rule space. If an effective search of ~$10^{9}$–$10^{10}$ variants preceded publication, the agreement is expected by chance. The paper does not document which case holds (CT-G10). Neither conclusion is supported over the other by the material in hand.

## Proposed next steps (formalization order)

1. **CT-G01/G02.** Attempt to construct the dressed insertion $M$ as an operator built from projections and twirls, with no reference to the record list. Test whether its diagonal is $\sqrt{k_r}$. This would settle whether the weight law is derived or chosen.
2. **CT-G03/G04.** Determine whether any axiom-level statement fixes the far-end conjugation. If none does, the uniqueness theorem should carry the convention as an explicit hypothesis.
3. Extend the Lean companion to prove the class-function form (Prop. 2) inside the upstream definitions.
4. Rerun the trial-factor analysis on a rule space generated mechanically from the Lean definitions (mutating each definition), rather than our hand-declared space.
5. Compendium: see `reports/chain-theory-compendium.md`.
