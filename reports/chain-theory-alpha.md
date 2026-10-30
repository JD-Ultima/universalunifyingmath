# Chain Theory — α(0) Letter: Phase 1–2 summary

**Source:** Wissner-Gross, *A Parameter-Free Candidate for the Fine-Structure Constant from Entangled Qudits*, v1.0.0, 2026-10-23, with Supplement and Lean file (`sources/chain-theory.md`, commit `a28db89`).
**Details:** `formalization/chain-theory/alpha/audit.md`.
**Registers:** `claims/chain-theory-alpha.md` (CT-001…023) and `gaps/chain-theory-alpha.md` (CT-G01…G12).

## What the Letter claims

1. One axiom (pure maximally entangled pairs of qudits, known only by projection onto the pair state) plus five definitions select a unique configuration: a qutrit shared by three qubit links, two of them "read". This is Theorem 1.
2. A weighted count over the 172 "measurement records" of that configuration gives $(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$, with $1/x = 137.035\,999\,189\,566\ldots$. This is Theorem 2.
3. This number is conjectured to equal $\alpha(0)$ (the Photon Conjecture). It lies $0.60\sigma$ from CODATA 2022.

## What we established

| Item | Result | Grade |
|---|---|---|
| Counting and value | Reproduced exactly by an independent $36\times36$-operator implementation. Agrees with the author's `records.csv` class by class. 51 printed digits (60 ancillary) agree. | **Verified** |
| Author's Lean file | Builds under Lean 4.35.0 in 170 s, axioms `propext`/`Quot.sound`/`Classical.choice` only. `main` is an equality of coefficient lists, i.e. arithmetic of the definitions. | **Verified** |
| Two exchanges | $X$ and $Y$ are the same involution $E$ on all product records. The weight law reduces to a class function: odd → 1, even mixed → $\sqrt2$, even same-end → $\sqrt3$, commutator → $\sqrt6$. | **Proven** (Prop. 2) |
| Selection of the (2,3) link | "Acting central charge" holds only for $d_1\ne d_2$ under the Letter's convention (A). Under the conjugate far-end convention (A*), $(2,2)$ also qualifies. Under a matched-direction reading (B), no type qualifies. Lean: `uniqueness_fails_under_conjugate_convention`. | **Proven**; uniqueness is convention-dependent |
| Hypercharge of the link | $Y=1/6$ ($Q_L$) under A; $Y=5/6$ (SU(5) leptoquark block) under A*. | **Proven** |
| Alternatives | The 12 weight-law alternatives of the author's Table S5 are reproduced; its 3 other rows lie outside the declared space. Four normalization alternatives the table omits give $1/x \in\{205.6, 60.9, 129.5, 188.9\}$. | **Computed** |
| Trial factor | In a declared 46 080-rule space containing the 12 weight-law alternatives of Table S5 (its 3 other rows change the counting unit, $N$ or the configuration and lie outside it), the expected number of chance hits within $1\sigma$ is $4.4\times10^{-6}$, even allowing selection of the best rule; a chance hit needs ~$10^{10}$ rules. In $\tfrac1{31104}\mathbb Z_{\ge0}[\sqrt2,\sqrt3]$ with coefficients ≤200 there are $1.6\times10^9$ elements; exactly one is within $1\sigma$ (the author's), against 0.024 expected. | **Computed** |
| Covariance | The Letter's plain exchange does not commute with the symmetry group of its own joint state, and its record span is not invariant. The covariant exchange (with transpose) restores invariance and gives identical counts and value. | **Proven** (F-CT-10) |
| Freedom left by symmetry | Under $SU(3)\times SU(2)\times\langle E'\rangle$ the record space has 10 isotypic components over ℂ (real commutant 22). A Hermitian-preserving invariant dressing has a 9-parameter real trace, and the Letter's weights are one point of that family. | **Computed** (F-CT-10) |
| Configurations | 76 / 42 / 23 reproduced. Selected: 1 configuration under A, 2 under A* (the competitor has $1/\text{value}=31.81$), none under B. | **Computed** |
| Read-cap lemmas | Choi's map is positive and non-decomposable by SDP, including after embedding. | **Cross-checked** |
| Algebraic data | $\alpha_\text{cand}$ has degree 4. The minimal polynomial and Galois conjugates are given. | **Proven** |
| Operator origin of the weights (CT-G01) | The dressed insertion has an exact operator form in the model's primitives, but its qutrit projector is the end rule (CT-G02). No word of length ≤ 3 in 38 primitives and no dressing local in the original frame reproduces the weights. A Galois lemma extends the $\sqrt{k_r}$ and trace negatives to all lengths, a factorization theorem the $k_r$ negative to all products of factorizing letters; for words with $\Pi_\pm$ or a read link the $k_r$ negative holds up to length 3. Transporting qubit letters through the read links reproduces the weights with a frame-local dressing, but only with two further choices: which exchange is assessed in which frame, and a rank-3 join credit for letters on a rank-2 support. | **Proven** (negative within declared families; F-CT-11, `operator-origin.md`, `transport.md`) |

## Assessment, by layer

- **(c) Proven.** The arithmetic and the selection theorem follow from the stated definitions under convention A.
- **(b) Reconstruction gaps.** The number rests on choices the axiom does not fix:
  - a conjugation convention and a direct-sum trace (CT-G03/G04);
  - three wording clauses assigning $\sqrt3$ vs $\sqrt2$ (CT-G02);
  - amplitude rather than Born weights (CT-G08);
  - the definition of $N$ (CT-G05);
  - an exception admitting single-slot commutator records (CT-G09).

  None of these is shown to be forced. Their published weight-law alternatives spread $1/x$ over $[97.5,\,180.8]$ (Table S5 as a whole: $[0.758,\,180.8]$).
- **(d) Conjectured.** The Photon Conjecture. Physically it lacks:
  - a map from the hypercharge readout to $\alpha_\text{em}$, given that $e=g'\cos\theta_W$ (CT-G06);
  - a propagator or momentum dependence for "on-shell" (CT-G07).
- **(e) Testable.** The candidate fixes every digit.
  - It is refuted if consistent determinations of $\alpha(0)$ with $u\lesssim2.5\times10^{-9}$ (8× better than CODATA 2022) exclude it at $5\sigma$.
  - Resolving the Rb ($-1.49\sigma$) vs electron $g-2$ ($+1.57\sigma$) tension is the nearer test.

**Evidential weight of the 10-digit match.** If the rules were fixed blind to $\alpha$, the agreement has $p\approx4\times10^{-6}$ against the declared rule space even allowing selection of the best rule in it, and of order $10^{-10}$ for a single rule fixed in advance. If an effective search of ~$10^{9}$–$10^{10}$ variants preceded publication, the agreement is expected by chance. The paper does not document which case holds (CT-G10). Neither conclusion is supported over the other by the material in hand.

## Proposed next steps (formalization order)

1. **CT-G02.** Done for CT-G01 (F-CT-11): the operator construction exists but carries the end rule as an input. What remains is a principle that fixes the end attribution (equivalently, the frame choice and join credit of `transport.md`), or a proof that none of a declared kind exists.
2. **CT-G03/G04.** Determine whether any axiom-level statement fixes the far-end conjugation. If none does, the uniqueness theorem should carry the convention as an explicit hypothesis.
3. Done: `weight_is_class_function` in `proofs/lean/chain-theory/ConventionCheck.lean` proves the class-function form inside the upstream definitions.
4. Partly done: the author's 19 Lean mutations are classified against the declared space (`mutations_classified.py`; 13 lie in it). A rule space generated mechanically from all Lean definitions remains to be built.
5. Compendium: see `reports/chain-theory-compendium.md`.
