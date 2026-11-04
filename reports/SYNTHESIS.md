# Research synthesis

Current account as of November 3, 2026. This project seeks precise rules that could explain observed structure; a common underlying mechanism remains an open research question. The case studies below separate what follows mathematically from what has a physical interpretation.

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

If the measured central value remains near CODATA 2022, uncertainty in $1/\alpha$ below about $2.5\times10^{-9}$ would separate it from the candidate at 5σ. This is a conditional sensitivity estimate.

[Report](chain-theory-alpha.md) · [Audit](../formalization/chain-theory/alpha/audit.md) · [Conjugation](../formalization/chain-theory/alpha/conjugation.md) · [Uniform records](../formalization/chain-theory/alpha/records-uniform.md) · [Mechanical family](../formalization/chain-theory/alpha/mechanical-trial-factor.md)

### Chain Theory: compendium 0.1

The audit distinguishes the source's index from statement-level derivations. Of 100 THEOREM boxes, 41 explicitly declare an OPEN premise in their derivation and at least 76 violate the source's chain rule under the documented reading. Several index cycles are chapter-citation artefacts; the separate cycle through the α theorem and gain lemma depends on the stated interpretation of chapter designators. Shared derivations and ambiguous any-of citations remain unresolved.

The Lean audit also distinguishes a successful build from coverage of the written claim. The published tree does not build in the audited 4.35.0 environment; five recorded proof edits yield a build of all 493 modules with no `sorry`. The source's pinned 4.32.1 environment was not rebuilt here. Of 64 theorem rows tied to Lean, the principal clause-counting rule gives **1 full match, 7 partial matches and 56 non-matches**; alternative counting rules are recorded in the coverage audit. Of 618 table-check lemmas, 560 check the occurrence of an entered string.

Some finite results do have substantive formal coverage, including gauge-group selection from the declared alphabet and read cap. The separate [read-cap reuse check](chain-theory-readcap-reuse.md) freshly compiles the universal positive-map theorem and its adapter to the author's selection theorem. It retains the assumptions on dimensions and unequal ends; operational read semantics and physical identifications remain separate. Physical comparisons are assessed separately: 11 of the 16 adopted numerical readings differ from their reference quantities by more than a factor of ten, while other printed readings give closer values. The scorecard records these alternatives and the identification choices required to compare them.

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

The following qualifications apply:

- **Existence threshold.** For a flat (4,0) base, Z1, $\lambda=1/5$ and $c_V>0$, the constant-field equations have real solutions iff $c_V\kappa^2\ge9360(x'^2+y'^2)$. Requiring $\det M\ne0$ makes this inequality **strict** for the rotated member with $a^2=125/576$: the equality case is a singular double root. Q15 defines the primed coefficients and rotation parameter.
- **Extra even fields.** R computes the full even operator, including non-vector kernels. The slaved-sector conclusions extend to arbitrary even fields if $\det M\ne0$ and the odd coefficients meet the other blocks' singular set only on a set with empty interior. T1–T3 extend R10 to a fixed **possibly nonconstant solved slaved baseline** in (14,0), compact $\mathfrak u(128)$, real fields, $\kappa\ne0$ and $\det M\ne0$: no additional even field solves the equation with that same odd field. For an unsolved baseline, the residual identity gives only a necessary sign condition. Other singular curves and indefinite real forms remain open.
- **Hook and other singular lines.** For $W_A$, T4 excludes $t=-\kappa/176$ and $-\kappa/208$ when $\operatorname{tr}_{g_Y}G^Y\ge0$. T5 forces $u$ to be constant on six positive hook lines. V lists nine necessary constant odd pairs per line before subsequent exclusions; its counterexample shows that scalar traces do not imply the full equation. [AZ](../formalization/geometric-unity/AZ-first-hook-classification.md) excludes the entire first pure hook ($t=3\kappa/136$, $u=0$) over every smooth positive four-base and every positive Z1 metric parameter $\lambda<1/4$. Exact polynomial minors prove the flat comparison map has rank 896 throughout that range; full metric rescaling and fibre density then force the field to vanish, contradicting the Euler Ricci identity. The other five pure hooks and all coincident factors remain open.
- **Lambda3 singular fields.** [AK](../formalization/geometric-unity/AK-lambda3-classification.md) classifies both single $W_A$ branches over every smooth positive four-base and positive Z1 metric with $\lambda<1/4$. The first has exactly $T_o=-i\kappa\Phi^1/400$, $T_e=0$, requiring a flat projected base, $\lambda=1/5$ and $c_V\kappa^2=150000/11$; the opposite branch has no solution. No curvature-trace sign or analyticity is imposed. [AF](../formalization/geometric-unity/AF-symmetric-gamma.md)/[AG](../formalization/geometric-unity/AG-other-lambda3.md) supply the endpoint reduction; AK combines a radial multiplicity criterion and the negative Euler Ricci eigenvalue to exclude three endpoints and fix $\lambda$ at the fourth. [AJ](../formalization/geometric-unity/AJ-curved-endpoint.md) supplies the remaining flatness result, and [AI](../formalization/geometric-unity/AI-flat-lambda3.md) the full-equation converse. Coincident kernels and other sectors are outside this classification.
- **Single Lambda5 fields.** [AW](../formalization/geometric-unity/AW-lambda5-conic-classification.md) completes the opposite-conic exclusion and combines it with [AP](../formalization/geometric-unity/AP-first-lambda5-obstruction.md): neither single $W_A$ Lambda5 branch admits a smooth full solution on an open positive-Z1 domain, under the compact-real, true-equation, invariant-odd, $\kappa\ne0$, regular-vector and noncoincident-factor assumptions. A normal Clifford moment forces a rational form-norm profile; Euler alignment and conservation require an incompatible rational identity. Smooth persistence also excludes any nonzero extra Lambda5 component of a full solution away from factor intersections, while zero-component conic contacts are not excluded by this result. The intermediate Gram, jet, crossing and spectral results are retained. Separate zero-order algebraic feasibility and AV's cubic implication remain undecided; domain nonexistence does not decide them.
- **Complete global examples.** [AM](../formalization/geometric-unity/AM-global-flat-models.md) realizes F7's positive observer over a flat base as the complete space $SL(5,\mathbb R)/SO(5)$ and constructs complete flat-torus quotients. [AO](../formalization/geometric-unity/AO-global-invariant-family.md) extends the known Q15 constant invariant family globally for every fixed compact J.0 member with nonzero Einstein coefficients, including axial fields and all sixteen torus spin choices. Its hand factorization gives the exact exceptional scale and regular configuration count, preserving Q15's strict-threshold exception. The singular vector block can occur in a smooth field solving the underlying constant-ansatz equation; it is excluded when regular slaving is assumed. Joint base variation, gauge inequivalence and physical interpretation remain separate.

The remaining slaved-vector Z1 twisted-cone case is in (14,0), at $\lambda=1/5$, away from the parallel set, with base Ricci rank at most one. At a nonzero rank-one Ricci point, write
$$\operatorname{Ric}_g=\rho\,p^\flat\otimes p^\flat,\qquad \rho\ne0,\quad |p|=1.$$
Three curvature conditions are excluded; electric and magnetic Weyl matrices are taken on $p^\perp$, with the magnetic matrix encoding the mixed part:

- **Vanishing mixed Weyl part (CC):** $W_g(p,u,v,w)=0$ for $u,v,w\perp p$, with arbitrary electric part.
- **Vanishing electric Weyl part (EE):** $W_g(p,u,p,v)=0$ for $u,v\perp p$, with arbitrary mixed part.
- **Joint Weyl part (FF):** the electric and magnetic matrices commute and have a common nonzero kernel vector. Both may be nonzero; this condition is narrower than arbitrary commuting matrices.

All three results allow arbitrary higher metric jets. Exact finite compatibility covers and coefficient-rank certificates at those same positive fibres prove the obstruction. Under S.0's invariant-odd/slaved-vector assumptions, a solution domain cannot meet the fibre over a point satisfying any of these conditions. Joint Weyl cases beyond FF's common-kernel family remain open at nonzero Ricci.

At **every Ricci-zero point with nonzero Riemann curvature**, [HH](../formalization/geometric-unity/HH-ricci-zero.md) excludes a fibre-open relaxed-system patch and, under S.0, a full-field solution domain meeting that fibre. All higher metric jets are allowed; Ricci zero is required only at the initial point. [GG](../formalization/geometric-unity/GG-ricci-zero-weyl.md) supplies the simpler case where both chiral blocks are singular. Nonzero Riemann curvature, rather than nonzero Ricci, supplies local persistence. These results do not cover zero-curvature points.

At a zero-curvature base point, Y excludes a nonzero conformal quartic direction from the relaxed system on the origin fibre. Its coefficient/augmented ranks are 28/29. This gives a fibre obstruction, without W's full-rank base-neighbourhood persistence or a local full-field exclusion.

Projection to a section and variation of a weighted action remain separate questions. The section-variation calculation gives a 4D vacuum Einstein equation for specified fibre weights in the treated sector, but this equation alone does not establish a joint solution of the full system. No new quantified physical prediction is established by the representation or reduction results.

The [GU report](geometric-unity-phase1.md) preserves the audit history. The principal evidence is grouped below.

| Subject | Proofs and certificates |
|---|---|
| Field equation and sector | [Variational equation](../formalization/geometric-unity/D-selfadjoint-shiab.md) · [Ellipticity](../formalization/geometric-unity/H-ellipticity.md) · [Slaved sector](../formalization/geometric-unity/J-slaved-sector.md) |
| Metric reductions | [Z2](../formalization/geometric-unity/P-z2-last-exception.md) · [Z1 threshold](../formalization/geometric-unity/Q-z1-remaining-members.md) · [Z1 off-section](../formalization/geometric-unity/S-z1-off-section.md) |
| Even operator and hooks | [Full operator](../formalization/geometric-unity/R-even-operator.md) · [Residual identity](../formalization/geometric-unity/T-singular-residual.md) · [Hook catalogue](../formalization/geometric-unity/V-hook-intersections.md) · [First pure-hook classification](../formalization/geometric-unity/AZ-first-hook-classification.md) |
| Lambda3 kernels | [Positive-Z1 classification and its inputs](../formalization/geometric-unity/AK-lambda3-classification.md) |
| Lambda5 kernels | [Completed single-conic classification](../formalization/geometric-unity/AW-lambda5-conic-classification.md) · [Gram constraints](../formalization/geometric-unity/AH-lambda5-gram.md) · [Full opposite system](../formalization/geometric-unity/AS-opposite-lambda5-system.md) · [Cubic transport and its algebraic gap](../formalization/geometric-unity/AV-opposite-gram-transport.md) |
| Global examples | [Complete metric and torus quotients](../formalization/geometric-unity/AM-global-flat-models.md) · [Constant invariant family and exact singular scale](../formalization/geometric-unity/AO-global-invariant-family.md) |
| Nonzero rank-one Ricci | [Vanishing mixed part](../formalization/geometric-unity/CC-diagonal-weyl.md) · [Vanishing electric part](../formalization/geometric-unity/EE-zero-electric-weyl.md) · [Joint Weyl family](../formalization/geometric-unity/FF-joint-weyl.md) |
| Earlier curvature certificates | [Zero Weyl](../formalization/geometric-unity/Z-cubic-rank-one.md) · [One electric direction](../formalization/geometric-unity/BB-weyl-direction.md) · [One mixed direction](../formalization/geometric-unity/DD-mixed-weyl.md) |
| Ricci-zero and zero-curvature points | [All nonzero-curvature Ricci-zero points](../formalization/geometric-unity/HH-ricci-zero.md) · [Singular-block subcase](../formalization/geometric-unity/GG-ricci-zero-weyl.md) · [Zero-curvature quartic control](../formalization/geometric-unity/Y-zero-curvature-quartic.md) |
| Base-metric variation | [Weighted section variation](../formalization/geometric-unity/O-gimel-variation.md) |

## 2. What the comparisons support

The recurring question is which defining choices determine a result. Convention changes, field sectors and observable-identification maps are tested explicitly. A counterexample applies to the statement and assumptions tested; it does not classify every possible reconstruction.

The group-theoretic connections are also precise but limited. Chain Theory's convention-A link lies in the $\mathbf{10}$ of SU(5), while its conjugate-convention link lies in the adjoint $\mathbf{24}$. GU's branching uses Spin(10)-type representations and Pati–Salam factors, but its internal real form is noncompact Spin(6,4). That distinction is explicit in the [branching analysis](../formalization/geometric-unity/B-spinor-branching.md). These representation comparisons do not establish common dynamics or a physical mechanism. [Pattern catalogue](../patterns/chain-theory.md)

## 3. Manuscript drafts

| Draft | Current scope | Source |
|---|---|---|
| Comment on the α Letter | Reproduction, convention dependence, covariance, operator weights and bounded rule-family sensitivity | [TeX](../proofs/latex/ct-alpha-comment/comment.tex) |
| Geometric Unity's algebraic core | Signatures, spinors, gauge identities, variation, ellipticity and reduction through the A–Q analyses. The newer results from R onward remain in research notes. | [TeX](../proofs/latex/gu-note/note.tex) |

The GU manuscript source last compiled on November 2, 2026; the Chain Theory manuscript source on November 1. Compilation checks document construction, not mathematical validity. The [verification guide](../computation/README.md) separates selected checks from long searches, fresh Lean builds and manuscript work.

<a id="4-next-steps-in-priority-order"></a>

## 4. Open questions and next work

The [prior-work survey](RELATED-WORK.md) identifies overlapping results and verification targets. In particular, the universal operator-theoretic read cap already has Lean proofs; the [fresh pinned rebuild and adapter](chain-theory-readcap-reuse.md) now verify that mathematical premise. Further work concerns operational read semantics and independent selection/weight inputs.

1. **GU remaining singular fields.** Investigate a declared remaining pure hook or intersection case using the full odd tensor equations. AZ closes the first pure hook over every smooth positive base for all $\lambda<1/4$; its former parameter problem is no longer a search target. The other pure cases, all coincident kernels, $\det M=0$, operator members, indefinite pairings and non-invariant odd fields retain separate scope. AW closes both single Lambda5 branches, so their former realization problem is no longer a search target. The AS/AV zero-order cubic implication remains a secondary algebraic question.
2. **GU remaining Z1 metrics.** Investigate nonzero rank-one Ricci points with both Weyl parts nonzero outside FF's commuting common-kernel family, and the zero-curvature higher-jet boundary. HH settles the nonzero-curvature Ricci-zero case within S.0; the general S15 twisted cone remains open. Choices of Z1/Z2, the off-section field rule and the fibre weight remain part of the reconstruction.
3. **Chain Theory inputs.** Seek an independent rule fixing the weight attribution, matching closure, and the uniform-record assumptions. Their conditional mathematics is reproduced; an axiom-level derivation remains open.
4. **Chain Theory scope.** Extend the mechanical family only with a declared mutation depth and resource bound. Resolve shared derivations and ambiguous premise citations before strengthening the compendium dependency conclusions. Physical comparisons need explicit identification maps.

The [gap registers](../gaps/) retain individual questions and closure criteria. Detailed internal cross-reviews are in the [review records](reviews/).
