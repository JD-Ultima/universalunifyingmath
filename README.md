# Universal Unifying Math

Research into the explicit mathematical rules that might explain physical structure. Finding an underlying “source code” remains an open goal.

The project examines **Chain Theory** and **Geometric Unity** through sourced claims, mathematical reconstructions, proofs and reproducible computations. The [research charter](RESEARCH_PROMPT.md) combines structural reasoning, thought experiments, useful approximations and empirical testing with explicit evidence standards.

Start with the [plain-language overview](reports/PLAIN-LANGUAGE-SYNOPSIS.md), the [current findings and open questions](reports/SYNTHESIS.md), or the [verification guide](computation/README.md). The [prior-work survey](reports/RELATED-WORK.md) identifies existing proofs and implementations to verify and reuse.

## Current findings

| Investigation | Established within the stated assumptions | Unresolved |
|---|---|---|
| **Chain Theory: α Letter** | Reproduces the 172 records and candidate inverse 137.035999189566…; verifies the imported [operator read cap](reports/chain-theory-readcap-reuse.md); identifies convention, record and weight choices. | An independent principle fixing those choices and the identification with electromagnetic coupling. [Report](reports/chain-theory-alpha.md) |
| **Chain Theory: compendium 0.1** | Audits premises, numerical comparisons and formal coverage. Of 64 Lean-tied theorem rows, the documented rule gives 1 full and 7 partial matches. | Ambiguous premises and the gap between the formal statements and attached physical claims. [Report](reports/chain-theory-compendium.md) |
| **Geometric Unity: 2021 draft** | Reconstructs equations, repairs gauge signs and gives scoped [field classifications](reports/SYNTHESIS.md#geometric-unity), [global models](formalization/geometric-unity/AO-global-invariant-family.md), and [singular-field and curvature results](reports/SYNTHESIS.md#geometric-unity). | Other singular fields, kernel intersections, remaining Z1 curvature cases and physical interpretation. [Current scope](reports/SYNTHESIS.md#geometric-unity) |

The [positive-Z1 hook results](reports/SYNTHESIS.md#geometric-unity) give the [full 6,834-dimensional H6 curvature kernel](formalization/geometric-unity/BO-h6-full-curvature-kernel.md) for the flat observer at $\lambda=1/5$. Exact [invariant quadratic constraints](formalization/geometric-unity/BQ-h6-invariant-moments.md) exclude each individual isotypic type and several sums, one of dimension 5,290. They force at least 41.4% of the squared norm of any full solution outside the original contained module. The [Cartan-square obstruction](formalization/geometric-unity/BP-h6-cartan-grade-nine.md) excludes a separate 184-dimensional sum using grade nine. A [rational certificate](formalization/geometric-unity/BL-mixed-grade-one-existence.md) gives a mixed point value satisfying grade one. General pointwise feasibility, full fields and physical identification remain open.

Internal reviews include AI-assisted work. External peer review remains outstanding. Source versions and acquisition limits are recorded in [sources](sources/), including the byte-for-byte comparison of the GU draft mirror with the recorded author-host download.

The [upper triangular closure construction](formalization/geometric-unity/BR-h6-closed-field-lifts.md) gives a 952-dimensional family of invariant closed real fields of constant norm. Finite sums of differently rotated fields provide analytic closed extensions of every value in a 6,534-dimensional submodule. An [exact positive receiver](formalization/geometric-unity/BV-h6-triangular-positive-obstruction.md) and a [diagonal limit](formalization/geometric-unity/BW-h6-triangular-orbit-averages.md) exclude every such sum and every fixed finite signed average of these fields from the full nonlinear equation on any open set.

An exact [second-derivative obstruction](formalization/geometric-unity/BU-h6-second-jet-obstruction.md) rules out the 14-dimensional type $(2,0)$ as a component of any smooth closed hook on this observer. This leaves a 6,820-dimensional upper bound on possible closed values.

An [exact partial point value](formalization/geometric-unity/BS-h6-certified-partial-value.md) now satisfies all grade-one and grade-thirteen equations and both invariant higher-grade projections. A required grade-nine coefficient is positive at that value, so full pointwise feasibility remains open.

## Reproduce selected checks

After preparing the pinned environment described in the [verification guide](computation/README.md#selected-suite), run:

```bash
python -X utf8 scripts/verify.py --recent
```

The guide covers exact upstream sources, completed validation and checks outside this suite. [CI](.github/workflows/verify.yml) runs the selected checks on Linux and Windows with Python 3.12 and 3.13.

## Evidence

- [Claims](claims/), [gaps](gaps/) and [patterns](patterns/): attributed statements, unresolved inputs and exploratory comparisons.
- [Formalizations](formalization/), [computations](computation/) and [Lean proofs](proofs/lean/): definitions, derivations, code and recorded results.
- [Internal reviews](reports/reviews/) and [corrections](reports/CORRECTIONS.md): objections, checks and revisions.
- [Manuscript drafts](proofs/latex/): coverage and status are listed in the [synthesis](reports/SYNTHESIS.md#3-manuscript-drafts).

## Contributing

State the question, source version, assumptions and reproducible check. Explain limitations and connect the result to its claim or gap record. Confirmed errors update the active account and enter the correction log.

No repository license has been selected. Upstream materials retain their respective terms.
