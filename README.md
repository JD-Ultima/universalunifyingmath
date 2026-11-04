# Universal Unifying Math

Research into the explicit mathematical rules that might explain physical structure. Finding an underlying “source code” remains an open goal.

The project examines **Chain Theory** and **Geometric Unity** through sourced claims, mathematical reconstructions, proofs and reproducible computations. The [research charter](RESEARCH_PROMPT.md) combines structural reasoning, thought experiments, useful approximations and empirical testing with explicit evidence standards.

Start with the [plain-language overview](reports/PLAIN-LANGUAGE-SYNOPSIS.md), the [current findings and open questions](reports/SYNTHESIS.md), or the [verification guide](computation/README.md). The [prior-work survey](reports/RELATED-WORK.md) identifies existing proofs and implementations to verify and reuse.

## Current findings

| Investigation | Established within the stated assumptions | Unresolved |
|---|---|---|
| **Chain Theory: α Letter** | Reproduces the 172 records and candidate inverse 137.035999189566…; verifies the imported [operator read cap](reports/chain-theory-readcap-reuse.md); identifies convention, record and weight choices. | An independent principle fixing those choices and the identification with electromagnetic coupling. [Report](reports/chain-theory-alpha.md) |
| **Chain Theory: compendium 0.1** | Audits premises, numerical comparisons and formal coverage. Of 64 Lean-tied theorem rows, the documented rule gives 1 full and 7 partial matches. | Ambiguous premises and the gap between the formal statements and attached physical claims. [Report](reports/chain-theory-compendium.md) |
| **Geometric Unity: 2021 draft** | Reconstructs equations, repairs gauge signs and gives scoped [field classifications](reports/SYNTHESIS.md#geometric-unity), [global models](formalization/geometric-unity/AO-global-invariant-family.md), and [Lambda5](formalization/geometric-unity/AW-lambda5-conic-classification.md) and [first-hook](formalization/geometric-unity/AZ-first-hook-classification.md) exclusions. | Other singular fields, kernel intersections, remaining Z1 curvature cases and physical interpretation. [Current scope](reports/SYNTHESIS.md#geometric-unity) |

The first pure-hook exclusion covers every smooth positive four-base throughout the positive Z1 range $\lambda<1/4$. Further [hook exclusions](formalization/geometric-unity/BD-affine-h4.md) cover affine-invariant pure $H_4$ and $H_6\Gamma$ fields over a flat base at $\lambda=1/5$; arbitrary fields in those branches remain open. GU conclusions depend on explicit metric, operator, real-form and field-sector assumptions. Numerical agreement, a finite check and a formal proof establish different things; physical identification requires additional evidence.

Internal reviews include AI-assisted work. External peer review remains outstanding. Source versions and acquisition limits are recorded in [sources](sources/), including the byte-for-byte comparison of the GU draft mirror with the recorded author-host download.

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
