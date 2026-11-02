# Research synthesis

Current account as of November 2, 2026. This project seeks precise rules that could explain observed structure; a common underlying mechanism remains an open research question. The case studies below separate what follows mathematically from what has a physical interpretation.

Results have undergone internal, AI-assisted review. They have not received external peer review. Detailed evidence and correction history remain in the linked reports, formalizations and [review records](reviews/).

## 1. Current findings

### Chain Theory: the α(0) Letter

The 172-record construction, sector counts and candidate value are independently reproduced:

$$\alpha_{\rm cand}^{-1}=137.035\,999\,189\,566\ldots$$

The candidate lies 0.60σ from the source's CODATA 2022 reference. Its identification with $\alpha_{\rm em}(0)$ remains a conjecture: the construction supplies a hypercharge direction, while the physical coupling and electroweak identification require further input.

The mathematical audit narrows the defining choices:

- **Selection.** The conjugate far-end convention admits a second configuration. Requiring matching closure, with matchings trivial off the common support, fixes the convention and trace condition. Within the Letter's stated proper-algebra definition, that closure requirement is equivalent to the desired selection; it is not derived from the axiom.
- **Records.** Covariant repairs retain the sector counts and value. A uniform definition replaces the commutator exception under the additional order and identification rules OD and ID. In a covariant frame, reproducing the value also chooses the qubit exchange sign R.
- **Weights.** An exact operator form exists, but its attribution projector encodes the end rule. The declared operator and transport searches do not derive an independent principle fixing the weights or the amplitude-level reading.
- **Rule sensitivity.** The candidate is the only distinct value within 1σ in both bounded rule families. In the mechanically generated family it is also the only value within relative $10^{-4}$; the hand-declared family has four such values. The density reference counts use illustrative random-target models and are not calibrated p-values or evidence of a construction search history.

The fixed numerical conjecture is testable. If a consistent determination remains near the CODATA 2022 central value, uncertainty in $1/\alpha$ below about $2.5\times10^{-9}$ would separate it from the candidate at 5σ. This is a conditional sensitivity estimate, not a forecast of a measurement.

[Report](chain-theory-alpha.md) · [Audit](../formalization/chain-theory/alpha/audit.md) · [Conjugation](../formalization/chain-theory/alpha/conjugation.md) · [Uniform records](../formalization/chain-theory/alpha/records-uniform.md) · [Mechanical family](../formalization/chain-theory/alpha/mechanical-trial-factor.md)

### Chain Theory: compendium 0.1

The audit distinguishes the source's index from statement-level derivations. Of 100 THEOREM boxes, 41 explicitly declare an OPEN premise in their derivation and at least 76 violate the source's chain rule under the documented reading. Several index cycles are chapter-citation artefacts; the separate cycle through the α theorem and gain lemma depends on the stated interpretation of chapter designators. Shared derivations and ambiguous any-of citations remain unresolved.

The Lean audit also distinguishes a successful build from coverage of the written claim. The published tree does not build in the audited 4.35.0 environment; five recorded proof edits yield a build of all 493 modules with no `sorry`. The source's pinned 4.32.1 environment was not rebuilt here. Of 64 theorem rows tied to Lean, the principal clause-counting rule gives **1 full match, 7 partial matches and 56 non-matches**; alternative rules give 1/8/55 or 0/7/57. Of 618 table-check lemmas, 560 check the occurrence of an entered string.

Some finite results do have substantive formal coverage, including gauge-group selection from the declared alphabet and read cap. Physical comparisons are assessed separately: 11 of the 16 adopted numerical readings differ from their reference quantities by more than a factor of ten, while other printed readings give closer values. The scorecard records these alternatives and the identification choices required to compare them.

[Report](chain-theory-compendium.md) · [Statement premises](../formalization/chain-theory/compendium/07-statement-premises.md) · [Lean content](../formalization/chain-theory/compendium/10-lean-table-rows.md) · [Numeric scorecard](../formalization/chain-theory/compendium/05-numeric-scorecard.md) · [Mass cells](../formalization/chain-theory/compendium/06-mass-cells.md)

### Geometric Unity

The algebraic work reconstructs bundle signatures, spinor branching, gauge identities and the Shiab-operator space. A global sign repair corrects the printed gauge-equivariance identities. When the connection is varied at fixed $\varepsilon$, the operator is independent of that connection, and the pairing is constant and invariant, the printed field equation is not variational for any nonzero admissible Shiab operator (D4). The reconstructed variational equation is therefore treated separately. Generic self-adjoint operators give the claimed elliptic complex in signature (14,0), with compact gauge algebra $\mathfrak u(128)$ and $\kappa\ne0$; the draft's own operator form does not.

The reduction results concern the **four-parameter self-adjoint family of J.0**, an invariant odd field $T_o=t\Phi^1+u\Phi^1\Gamma$, $\kappa\ne0$, and an even field slaved in the vector sector with $\det M\ne0$. They are not a classification of arbitrary fields in the theory.

| Metric and real form | Current result for the full equation $\Upsilon=0$ in that sector |
|---|---|
| Z2, (14,0) or (7,7) | No solution on any open subset for members with nonzero Einstein coefficients $(x,y)$. |
| Z1, (7,7) | No solution on any open subset for $(x,y)\ne0$. |
| Z1, (14,0) | A solution requires $\lambda=1/5$. Except for the remaining twisted-cone case, solutions are the classified flat-base, constant-odd-field solutions with $T_e=0$, without an analyticity assumption. |

For $(x,y)=0$, the only solution in this sector is $T=0$ under either metric. Some steps in these results are computer-assisted; the exact certificates and independent checks are identified in the notes.

Three boundaries matter:

- **Existence threshold.** For a flat (4,0) base, Z1, $\lambda=1/5$ and $c_V>0$, the constant-field equations have real solutions iff $c_V\kappa^2\ge9360(x'^2+y'^2)$. Requiring $\det M\ne0$ makes this inequality **strict** for the rotated member with $a^2=125/576$: the equality case is a singular double root. Q15 defines the primed coefficients and rotation parameter.
- **Extra even fields.** R computes the full even operator, including non-vector kernels. The slaved-sector conclusions extend to arbitrary even fields if $\det M\ne0$ and the odd coefficients meet the other blocks' singular set only on a set with empty interior. T1–T3 extend R10 to a fixed **possibly nonconstant solved slaved baseline** in (14,0), compact $\mathfrak u(128)$, real fields, $\kappa\ne0$ and $\det M\ne0$: no additional even field solves the equation with that same odd field. For an unsolved baseline, the residual identity gives only a necessary sign condition. Other singular curves and indefinite real forms remain open.
- **Specified singular lines.** For $W_A$, T4 excludes $t=-\kappa/176$ and $-\kappa/208$ when $\operatorname{tr}_{g_Y}G^Y\ge0$. T5 forces $u$ to be constant on each of six hook lines with $t/\kappa>0$. V makes the exception set explicit: each fixed hook line has nine possible constant odd pairs, with additional necessary norm constraints. Its grade-five counterexample passes the two invariant trace tests but fails the full equation. Constant-field hook solutions remain open.

The remaining Z1 twisted cone is in (14,0), at $\lambda=1/5$, away from the parallel set, with base Ricci rank at most one. U and W extend S15's single excluded jet to metrics with U1's quadratic part, zero cubic part, **arbitrary symmetric quartic terms**, and arbitrary smooth higher terms. They cover every nonzero Ricci amplitude $\rho$, unit direction $p$ and positive metric scales. At $0$, $\mathrm{Ric}_g=\rho p\otimes p$ and Weyl curvature vanishes; no global rank-one condition is asserted. W's three exact equations remain inconsistent because a gradient cancellation removes every quartic coefficient from those rows. Covariance and fibre rationality exclude a fibre-open consistent patch. Continuity excludes a solution domain meeting the fibre over $0$ within S.0's invariant-odd/slaved-vector sector. General cubic jets, arbitrary Weyl curvature and the Ricci-flat nonflat case remain open.

Projection to a section and variation of a weighted action remain separate questions. The section-variation calculation gives a 4D vacuum Einstein equation for specified fibre weights in the treated sector, but this equation alone does not establish a joint solution of the full system. No new quantified physical prediction is established by the representation or reduction results.

[Report](geometric-unity-phase1.md) · [Variational equation](../formalization/geometric-unity/D-selfadjoint-shiab.md) · [Ellipticity](../formalization/geometric-unity/H-ellipticity.md) · [Sector and hypotheses](../formalization/geometric-unity/J-slaved-sector.md) · [Z2 classification](../formalization/geometric-unity/P-z2-last-exception.md) · [Threshold](../formalization/geometric-unity/Q-z1-remaining-members.md) · [Full even operator](../formalization/geometric-unity/R-even-operator.md) · [Z1 off-section analysis](../formalization/geometric-unity/S-z1-off-section.md) · [Singular residual](../formalization/geometric-unity/T-singular-residual.md) · [Quadratic family](../formalization/geometric-unity/U-rank-one-frontier.md) · [Hook catalogue](../formalization/geometric-unity/V-hook-intersections.md) · [Arbitrary quartic terms](../formalization/geometric-unity/W-quartic-rank-one.md) · [Section variation](../formalization/geometric-unity/O-gimel-variation.md)

### Howard: specified arithmetic and geometry

With ordinary addition fixed, no ring, semiring or near-ring on ℕ, ℤ, ℚ or ℝ satisfies both $1\otimes1=2$ and $2\otimes2=4$. The seven-axiom consistency lattice is decided with z3 and Lean certificates. The conclusion is carrier-dependent: ℤ/4 supplies a ring model, and weaker axiom sets admit other operations. Several numerical statements and the specified six-pentagon closed-solid reconstruction fail their stated tests.

The primary book and full relevant patent texts were not obtained. In particular, the reading of $2\otimes2=4$ relies on an interpretation of the available wording. Physical and chemistry claims remain outside this mathematical analysis.

[Report](howard-phase1.md) · [Arithmetic](../formalization/howard/arithmetic.md) · [Geometry](../formalization/howard/geometry.md) · [Sources](../sources/howard.md)

## 2. What the comparisons support

The recurring question is which defining choices determine a result. Convention changes, alternate operations, field sectors and observable-identification maps are tested explicitly. A counterexample applies to the statement and assumptions tested; it does not classify every possible reconstruction.

The group-theoretic connections are also precise but limited. Chain Theory's convention-A link lies in the $\mathbf{10}$ of SU(5), while its conjugate-convention link lies in the adjoint $\mathbf{24}$. GU's branching uses Spin(10)-type representations and the Pati–Salam subgroup. These embeddings provide a mathematical comparison; they do not establish a common dynamics or physical mechanism. [Pattern catalogue](../patterns/chain-theory.md)

## 3. Manuscript drafts

| Draft | Current scope | Source |
|---|---|---|
| Comment on the α Letter | Reproduction, convention dependence, covariance, operator weights and bounded rule-family sensitivity | [TeX](../proofs/latex/ct-alpha-comment/comment.tex) |
| Geometric Unity's algebraic core | Signatures, spinors, gauge identities, variation, ellipticity and reduction through the A–Q analyses. The newer R–W results remain in research notes. | [TeX](../proofs/latex/gu-note/note.tex) |
| Arithmetic reconstruction of “1 × 1 = 2” | Seven-axiom consistency lattice, alternative models and formal certificates; no novelty claim | [TeX](../proofs/latex/howard-note/note.tex) |

All three sources compiled with the desktop editor's compiler on November 1, 2026. Compilation checks document construction, not mathematical validity or external acceptance. The [verification guide](../README.md#reproduce-selected-checks) distinguishes the selected executable suite from long searches, fresh Lean builds and manuscript work.

<a id="4-next-steps-in-priority-order"></a>

## 4. Open questions and next work

1. **GU non-vector singular branches.** Use T's residual identity to solve or obstruct the odd equation for $T_e=\Sigma(dt,du)+k$ on the remaining singular curves. The constant $W_A$ hook cases and their coincident kernels also remain open. Keep $\det M=0$, indefinite pairings and non-invariant odd fields as separate cases.
2. **GU remaining Z1 metrics.** Extend W's four-jet obstruction to nonzero cubic data or exhibit consistent data in S15's twisted cone. Arbitrary Weyl components and Ricci-flat nonflat metrics remain untreated. Choices of Z1/Z2, the off-section field rule and the fibre weight remain part of the reconstruction.
3. **Chain Theory inputs.** Seek an independent rule fixing the weight attribution, matching closure, and the uniform-record assumptions. Their conditional mathematics is reproduced; an axiom-level derivation remains open.
4. **Chain Theory scope.** Extend the mechanical family only with a declared mutation depth and resource bound. Resolve shared derivations and ambiguous premise citations before strengthening the compendium dependency conclusions. Physical comparisons need explicit identification maps.
5. **Howard definitions.** Obtain the missing primary sources and apply the consistency analysis to the operations actually defined there.

The [gap registers](../gaps/) retain the individual questions and closure criteria. Current independent reviews include [GU R](reviews/REVIEW12-GU-R-2026-11-01.md), [GU S](reviews/REVIEW12-GU-S-2026-11-01.md), [GU T](reviews/REVIEW13-GU-T-2026-11-02.md), [GU U](reviews/REVIEW14-GU-U-2026-11-02.md), [GU V](reviews/REVIEW16-GU-V-2026-11-02.md), [GU W](reviews/REVIEW15-GU-W-2026-11-02.md) and [Chain Theory](reviews/REVIEW6-CT-2026-11-01.md).
