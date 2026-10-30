# Synthesis: three unification frameworks, graded

Status as of 2026-10-30. Phase 1 is complete. Phase 2 has finished four rounds on GU and three on the Chain Theory compendium. The detailed reports, which this page only summarises, are:

- `reports/chain-theory-alpha.md`
- `reports/chain-theory-compendium.md`
- `reports/geometric-unity-phase1.md`
- `reports/howard-phase1.md`

"Proven" means that the result follows from stated definitions, with code and/or Lean in the repo. It never means physically confirmed.

## 1. Ledger

| Framework | What is precise | Proven by us | Found false / inconsistent | Open (well-posed) | Empirically testable now |
|---|---|---|---|---|---|
| **Chain Theory — α Letter** | All of it, as definitions. Lean checks the arithmetic. | Counts and value reproduced independently (F-CT-1). The two exchanges are one involution (F-CT-2, Lean). Uniqueness depends on a convention (F-CT-3, Lean). The record set is not covariant but has a covariant repair with the same value (F-CT-10). | Theorem 1 does not hold outside convention A. Under A* a second configuration qualifies, with $1/\text{value}=31.81$. | Whether the weight law has an operator origin (CT-G01; symmetry leaves 9 free real parameters). Hypercharge → $\alpha_\text{em}$ (CT-G06). | **Yes.** $1/\alpha=137.035\,999\,189\,566\ldots$ is falsifiable at $u\lesssim2.5\times10^{-9}$. |
| **Chain Theory — compendium 0.1** | The census and the premise graph. Some finite counts. | The census reproduces. 66/100 theorems rest directly on open premises, and the α theorem sits on premise cycles. | Weak angle (1/2, 1067σ); $a_e$ third order = 0, against $C_3=1.18$; hydrogen scale off by ~$10^9$. Scorecard over all 53 named quantities: at the readings the source adopts, 1 within 3σ, 1 more within 1%, 3 within a factor 10, 11 worse. | Statement-level premises. Identification maps (record order → $\alpha/\pi$, units). | Only where identification maps exist. Where they do, they fail. |
| **Geometric Unity — 2021 draft** | Finite-dimensional algebra: observerse signature, Clifford/spinor content, gauge-group algebra. | (7,7) iff Lorentzian base, $\lambda>\tfrac14$ and a fixed relative sign (cross-checked). One Majorana–Weyl spinor = one SM generation with $\nu_R$; the SM tables are reproduced; the 144 carries a further generation (cross-checked). The space of Shiab operators is computed exactly ($\dim 200$, self-adjoint part 48). Stueckelberg factorisation and the off-shell Noether identity hold for the true Υ. The reduction along a section gives 4D Einstein plus Λ. A generic self-adjoint Shiab gives an elliptic complex in (14,0). | Lemmas 6.5 and 7.2 false as printed (sign; repairable). Field equation (9.7) is not the variational equation of any functional, for *every* nonzero admissible Shiab operator (round 2, D4). Self-adjoint Einstein Shiab operators do exist (D3). "Spin(14,ℂ)→U(128)" is ill-typed. The claimed Euclidean-elliptic deformation complex fails for the draft's form (8.1), but holds for a generic self-adjoint Shiab operator in (14,0) with $\mathfrak u(128)$ (round 4). A single Majorana–Weyl reading admits no fermion kinetic term in (7,7), so the content doubles to complex Weyl ($6\times16\oplus2\times144$). | Which self-adjoint, elliptic Shiab (if any) the theory should use; the ℷ-variation. Choice between Z1 and Z2, and $\lambda$, $c_V/c_H$. Under Z2 the reduction gives 4D Einstein + $\Lambda=9c_H/2c_V$ (AdS sign). | No. Exotic states (charges $\pm\tfrac43,\pm\tfrac53,\pm2$) are predicted without masses. |
| **Howard** | The arithmetic, once rules are fixed. The numerical statements. One geometric reconstruction. | No ring, semiring or near-ring on ℕ, ℤ, ℚ, ℝ (more generally, in characteristic not dividing 4) realises $1\otimes1=2$ with $2\otimes2=4$ (Lean; ℤ/4 is a ring model). Bi-additive readings (on ℝ under regularity) are relabelled standard arithmetic; other consistent witnesses such as $2b$ are not. The axiom lattice is fully decided. | "Most numbers satisfy $2\sqrt n>n$"; the speed-of-light figures; the 6-pentagon closed solid. | Requires the primary book (not obtained). | Chemistry/physics claims are out of the mathematical scope. |

## 2. Cross-framework findings

1. **Sign and convention choices are load-bearing in all three.**
   - Chain Theory selects its unique link only under one far-end conjugation convention (F-CT-3).
   - GU's (7,7) signature needs a fixed relative block sign and $\lambda>\tfrac14$, and its gauge-equivariance lemmas hold only after one sign repair (C3/C3′).
   - Howard's no-go depends on reading $2\otimes2=4$ from one sentence.

   In each case the repo states the result separately for each choice, so the dependence is explicit rather than hidden.
2. **Two routes to the same group theory.** Chain Theory gets $S(U(2)\times U(3))\subset SU(5)$ from the $2\otimes3$ link. GU gets Pati–Salam from the maximal compact subgroup $\mathrm{Spin}(6)\times\mathrm{Spin}(4)\subset\mathrm{Spin}(6,4)$, and a 16 of Spin(10) type from one Majorana–Weyl spinor. These are the two classical maximal-subgroup routes into Spin(10) (Georgi–Glashow vs Pati–Salam; Baez–Huerta 2010). Neither framework's group-theoretic output goes beyond what those routes already give.
   - **Cross-constraint.** Under Chain Theory's convention A the link carries $(3,2)_{1/6}$, a matter multiplet that sits in GU's 16. Under A* it carries $(3,2)_{-5/6}$, which sits in the adjoint 45 (gauge sector). So whether CT's photon link is "matter" or "gauge" in Spin(10) language is exactly the unresolved convention CT-G03.
3. **Where numbers meet experiment, there is one striking agreement and several clear disagreements.**
   - The α Letter agrees with CODATA at $0.6\sigma$. Relative to every rule variant listed in the paper, a chance agreement this close has $p\sim4\times10^{-6}$ (F-CT-4). The construction history needed to interpret that figure is not documented.
   - The compendium's own adopted identifications mostly fail: the weak angle, $a_e$, the hydrogen scale (scorecard: 11 of 16 adopted readings are off by more than a factor 10). Unselected printed readings come close for $m_p/m_e=(2871/67)^2$ ($+1.8\times10^{-5}$), $m_\mu/m_e$ ($-0.09\%$) and on-shell $\sin^2\theta_W=(17/36)^2$ ($-0.8\sigma$). The expected chance hits within each declared space are $0.0025$, $0.007$ and $0.02$–$0.14$. Each needs a reading the source does not adopt. A round-3 derivation audit found that $2871$ is *entered* as a candidate, the integer nearest $67\sqrt{m_p/m_e}$, and that the source's own printed mass laws spread the proton prediction over a $1.2\times10^{-4}$ band containing CODATA. The lepton agreements need a face choice the compendium's own definitions exclude. Graded: selection or conditional, not derivation.
   - Howard's numerical claims fail or are not specified.
4. **Formal verification checks the definitions, not the physics.** The Lean coverage in both Chain Theory sources certifies finite arithmetic of definitions, and in the compendium often only label counts. Our own Lean files (Chain Theory convention check and class-function theorem; GU tilt-sign and Stueckelberg lemmas; Howard arithmetic) are stated with that limitation in mind.

## 3. Candidates for peer-reviewed output

| Draft | Content | Location |
|---|---|---|
| Comment on the α Letter | F-CT-1–4 and F-CT-10: reproduction, single involution, convention dependence (Lean), trial factor, covariance repair | `proofs/latex/ct-alpha-comment/comment.tex` |
| Note on GU's algebraic core | Signature classification; spinor content (fourth generation-equivalent; complex-Weyl doubling); sign error with repair; non-variational field equation for every admissible Shiab; Stueckelberg/Noether; ellipticity (generic yes, draft form no); 14→4 reduction | `proofs/latex/gu-note/note.tex` (draft) |
| Note on "$1\times1=2$" | Complete axiom lattice of seven rules, with z3 and Lean certificates. Low novelty; suitable as an expository note. | `proofs/latex/howard-note/note.tex` (draft) |

## 4. Next steps, in priority order

1. GU: the ℷ-variation; the choice between Z1 and Z2; the reduction and fermion content under the elliptic operator $W_A+W_C\circ\mathrm{Ad}_\Gamma$ of H5. (The ε-variation and the fermion operator are done: round 3.)
2. CT-G01/G02: the 9-parameter (real) family from F-CT-10. Is there any principle beyond the wording clauses that selects the Letter's point in it?
3. Howard: obtain *One Times One Equals Two* and rerun the axiom lattice on its stated products.
4. Compendium: statement-level premise extraction and re-grading.
