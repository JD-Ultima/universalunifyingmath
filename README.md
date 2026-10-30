# universalunifyingmath

**A reproducible mathematical audit of three proposed unification frameworks**

| Framework | Author | Sources audited |
|---|---|---|
| Chain Theory | A. D. Wissner-Gross | the α(0) Letter (v1.0.0) with Supplement and Lean file; the compendium (release 0.1.0) |
| Geometric Unity | E. Weinstein | the 2021 author's working draft |
| "Terryology" / *One times one equals two* | T. Howard | interviews, a podcast transcript, patent filings known from summaries (the book itself was not obtainable) |

For each framework the repository reconstructs the strongest precise mathematics the sources support, proves or refutes what can be decided, and states plainly what remains open. Mathematical validity and physical validity are graded separately. Every number in the write-ups is printed by a committed script, or checked by Lean, or cited.

> **Status (2026-10-30).** Working research repository. The results have been through three internal referee rounds ([`reports/reviews/`](reports/reviews/)) but **not external peer review**. Corrections are logged in [`reports/CORRECTIONS.md`](reports/CORRECTIONS.md).

---

## Contents

- [Start here](#start-here)
- [Summary of findings](#summary-of-findings)
- [Reproducing the results](#reproducing-the-results)
- [Repository layout](#repository-layout)
- [Grading conventions](#grading-conventions)
- [Manuscript drafts](#manuscript-drafts)
- [Scope and limitations](#scope-and-limitations)
- [Reporting an error](#reporting-an-error)

## Start here

| If you want… | Read |
|---|---|
| The graded overview of all three frameworks | [`reports/SYNTHESIS.md`](reports/SYNTHESIS.md) |
| The Chain Theory α(0) audit | [`reports/chain-theory-alpha.md`](reports/chain-theory-alpha.md) → [`formalization/chain-theory/alpha/audit.md`](formalization/chain-theory/alpha/audit.md) |
| The Chain Theory compendium audit | [`reports/chain-theory-compendium.md`](reports/chain-theory-compendium.md) → [`formalization/chain-theory/compendium/`](formalization/chain-theory/compendium/) |
| The Geometric Unity audit | [`reports/geometric-unity-phase1.md`](reports/geometric-unity-phase1.md) → [`formalization/geometric-unity/`](formalization/geometric-unity/) |
| The Howard audit | [`reports/howard-phase1.md`](reports/howard-phase1.md) → [`formalization/howard/`](formalization/howard/) |
| Every individual claim and its status | [`claims/`](claims/) (IDs `CT-…`, `GU-…`, `TH-…`) |
| Where the sources are underspecified | [`gaps/`](gaps/) |
| What was found wrong in our own earlier statements | [`reports/CORRECTIONS.md`](reports/CORRECTIONS.md) |

## Summary of findings

The statements below are summaries. Each has a precise version, with hypotheses, in the linked file.

### Chain Theory — the α(0) Letter

The Letter proposes $\alpha_{\rm cand}=(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$, i.e. $1/\alpha_{\rm cand}=137.035\,999\,189\,566\ldots$, from an axiom about maximally entangled qudit pairs.

- **Reproduced exactly.** An independent operator-level implementation reproduces the 172 records, the sector sizes 79/48/42/3 and the value. It lies 0.6σ from CODATA 2022. ([F-CT-1](formalization/chain-theory/alpha/audit.md))
- **The uniqueness theorem depends on a convention.** The selection of the qutrit–qubit link holds under one far-end conjugation convention. Under the conjugate convention a second configuration satisfies every hypothesis. This is checked in Lean against the author's own definitions. ([F-CT-3](formalization/chain-theory/alpha/audit.md), [`ConventionCheck.lean`](proofs/lean/chain-theory/ConventionCheck.lean))
- **Covariance.** The records as defined are not covariant under the symmetry group of their own joint state. A covariant repair leaves the value unchanged. Symmetry alone leaves a 9-parameter (real) family of readouts. ([F-CT-10](formalization/chain-theory/alpha/audit.md))
- **The weight law is an input, not a derivation.** The weights have an exact operator form in the model's primitives, but its qutrit projector is the "end rule" written as an operator. No local or short multiplicative construction reproduces them. ([F-CT-11](formalization/chain-theory/alpha/audit.md); [`operator-origin.md`](formalization/chain-theory/alpha/operator-origin.md), [`transport.md`](formalization/chain-theory/alpha/transport.md))
- **Trial factor.** The declared rule space has 46 080 variants and contains all twelve weight-law alternatives the author lists. In it the expected number of chance values within 1σ is $4.4\times10^{-6}$. A space generated mechanically from the author's Lean file (332 308 single and paired mutants, checked against Lean) gives $1.3\times10^{-6}$ ([mechanical space](formalization/chain-theory/alpha/mechanical-trial-factor.md)). Whether the rules were fixed before the value was computed is not documented. ([F-CT-4](formalization/chain-theory/alpha/audit.md))
- **Testability.** The candidate fixes every digit. It is falsifiable by determinations of $\alpha(0)$ with absolute uncertainty $\lesssim2.5\times10^{-9}$ in $1/\alpha$.

### Chain Theory — the compendium (release 0.1)

- **Grades against the source's own chain rule.** At statement level, 41 of 100 THEOREMs declare an OPEN premise in their own derivation, and at least 76 violate the rule. ([`07-statement-premises.md`](formalization/chain-theory/compendium/07-statement-premises.md))
- **Premise cycle.** Under the better-supported reading of the source's chapter designators, a premise cycle runs through the α theorem: 124.1 → the gain lemma (chapter 18) → 124.1. The cycles shown in the source's own index are an artefact of chapter-level citation.
- **Numbers against experiment.** Where the source's adopted identifications give numbers, they mostly disagree with experiment: the weak angle, the third-order $a_e$ coefficient, and the hydrogen scale (off by about $10^9$). ([`05-numeric-scorecard.md`](formalization/chain-theory/compendium/05-numeric-scorecard.md), [`08-hydrogen-registers.md`](formalization/chain-theory/compendium/08-hydrogen-registers.md))
- **The compendium's Lean.** The source pins Lean 4.32.1, and the tree does not build under Lean 4.35.0 as published. With five proof-only edits all 493 modules build, with no `sorry`. Of the 64 theorem rows the compendium ties to Lean, 1 is fully stated by its Lean, 7 partly, and 56 not at all: most table checks only confirm that an entered string occurs. ([`09-lean-audit-at-scale.md`](formalization/chain-theory/compendium/09-lean-audit-at-scale.md), [`10-lean-table-rows.md`](formalization/chain-theory/compendium/10-lean-table-rows.md))
- **Gauge-group selection.** Certified in core Lean (71 theorems). With alphabet {2,3,4}, the link (3,2) is selected uniquely exactly for read caps 6 and 7. ([`GaugeSelection.lean`](proofs/lean/chain-theory/GaugeSelection.lean))

### Geometric Unity — 2021 draft

- **Signature.** The chimeric bundle has signature (7,7) for bases (1,3), (3,1) and (4,0) iff the base is Lorentzian, $\lambda>\tfrac14$, and a fixed relative block sign holds; a (2,2) base also gives (7,7). It has no Spin(7,7) structure over $\mathrm{Met}(X)$. Over the time-oriented double cover it has one iff $w_2(X)=0$. ([A](formalization/geometric-unity/A-chimeric-signature.md))
- **Spinor content.** One Majorana–Weyl spinor gives one Standard Model generation with $\nu_R$, and the 144 carries a further generation's worth of states (independently cross-checked). In (7,7) a single Majorana–Weyl field has no fermion kinetic term, so the content doubles to complex Weyl. ([B](formalization/geometric-unity/B-spinor-branching.md), [G](formalization/geometric-unity/G-full-el-fermions.md))
- **Sign error.** Lemmas 6.5 and 7.2 are false as printed. One global sign choice repairs them, and it changes at least two printed formulas. Lean-checked. ([C](formalization/geometric-unity/C-gauge-shiab.md), [`GUTilt.lean`](proofs/lean/geometric-unity/GUTilt.lean))
- **Field equation.** The space of Shiab operators is computed exactly: dimension 200, with a self-adjoint part of dimension 48. Under stated hypotheses, the printed field equation (9.7) is the variational equation of no functional, for every nonzero admissible Shiab operator. ([D](formalization/geometric-unity/D-selfadjoint-shiab.md))
- **Ellipticity.** The claimed Euclidean-elliptic deformation complex fails for the draft's own Shiab form. It holds for a generic self-adjoint Shiab operator in (14,0). ([H](formalization/geometric-unity/H-ellipticity.md))
- **14→4 reduction.** Along a section, 4D Einstein with a cosmological term appears in the projected sense. With the displacement field included, the field equation has no solution near the section in which the displacement vanishes. In the invariant sector in (7,7) it has no solution near the section at all. ([F](formalization/geometric-unity/F-zorro-reduction.md), [I](formalization/geometric-unity/I-elliptic-reduction.md))

### Howard

- **"1 × 1 = 2".** No ring, semiring or near-ring on ℕ, ℤ, ℚ or ℝ (more generally, in characteristic not dividing 4) has $1\otimes1=2$ together with $2\otimes2=4$. This is Lean-checked; ℤ/4 is a ring model of both. The lattice of seven candidate axioms is fully decided, with z3 and Lean certificates. ([`arithmetic.md`](formalization/howard/arithmetic.md), [`HowardArith.lean`](proofs/lean/howard/HowardArith.lean))
- **Numerical statements.** Several are false as stated, e.g. "most numbers satisfy $2\sqrt n>n$" and the speed-of-light figures. ([`numerics.md`](formalization/howard/numerics.md))
- **Primary source missing.** The primary book was not obtainable, so the audit is limited to the sources listed in [`sources/howard.md`](sources/howard.md).

## Reproducing the results

**Software**

- **Python 3.13** with the versions used:

```bash
pip install numpy==2.5.3 sympy==1.14.0 mpmath==1.3.0 scipy==1.18.1 cvxpy==1.6.0 z3-solver==5.3.3.0 mendeleev==1.3.0
```

  Most scripts need only numpy, sympy and mpmath. scipy, cvxpy, z3 and mendeleev are each used by a few scripts, and each folder's README says which.
- **Lean 4.35.0**, core library only (no Mathlib).

**Computations.** Each folder under [`computation/`](computation/) has a README with the exact command, runtime and output file for every script. The outputs are committed beside the scripts, so a result is reproduced by rerunning and diffing:

```bash
cd computation/chain-theory/alpha
python3 records_independent.py | diff - records_independent.out
```

Randomness, where used, has fixed seeds. Most scripts run in seconds. The longest searches take about 15–20 minutes on 4 cores, and longer on a shared machine.

**Lean.** Each check prints `#print axioms` for every theorem.

```bash
lean proofs/lean/howard/HowardArith.lean
```

```bash
lean proofs/lean/geometric-unity/GUTilt.lean
```

```bash
lean proofs/lean/geometric-unity/GStueckelberg.lean
```

```bash
LEAN=/path/to/lean ./proofs/lean/chain-theory/run_gauge.sh
```

```bash
LEAN=/path/to/lean UPSTREAM=/path/to/chain-theory ./proofs/lean/chain-theory/run.sh
```

`run.sh` checks the author's `ChainAlpha.lean` against its pinned sha256 before appending our file.

**Upstream sources** are not redistributed. They are pinned by commit and sha256 in [`sources/`](sources/). Chain Theory is public at `github.com/alexwg/chain-theory` (commit `a28db89`), and the scripts that read it take its path as an argument.

## Repository layout

```
sources/         annotated, version-pinned bibliographies (sha256) and the external theorems used
claims/          claim registers: every source claim, its type, and its current status
gaps/            points where a source is underspecified, with what would close each one
patterns/        conjectures and cross-framework observations (exploration, not results)
formalization/   statements and proofs, one file per topic
proofs/lean/     Lean 4 checks (core library), with recorded outputs
proofs/latex/    manuscript drafts
computation/     scripts and their recorded outputs; a README in each folder
metrics/         one YAML record per result (status, reproducibility, confidence)
reports/         phase summaries, SYNTHESIS.md, CORRECTIONS.md, reviews/
```

## Grading conventions

- **Proven**: follows from stated definitions, with a proof in `formalization/` and code or Lean in the repo. It never means "physically confirmed".
- **Computed / Verified**: an exact or numerical computation, committed with its output. Verified means the source's own statement was reproduced.
- **Conditional**: holds under a stated reading or convention that the source does not fix. Each reading is reported separately.
- **False / Inconsistent**: contradicted by a proof or an exact counterexample.
- **Open**: well-posed but unresolved, recorded in `gaps/`.

Negative results from searches are stated relative to a declared family, whose size is given.

## Manuscript drafts

These are LaTeX sources. They were not compiled in the environment used to build this repository.

| Draft | Path |
|---|---|
| Comment on the α(0) Letter: reproduction, convention dependence, covariance, weight law, trial factor | [`proofs/latex/ct-alpha-comment/comment.tex`](proofs/latex/ct-alpha-comment/comment.tex) |
| The algebraic core of Geometric Unity: signatures, spinor content, gauge equivariance, field equation, ellipticity, 14→4 reduction | [`proofs/latex/gu-note/note.tex`](proofs/latex/gu-note/note.tex) |
| Which arithmetics satisfy "1 × 1 = 2"? | [`proofs/latex/howard-note/note.tex`](proofs/latex/howard-note/note.tex) |

## Scope and limitations

- **Mathematics, not physics.** Physical claims are assessed only where a source gives a number that can be compared with measurement.
- **Source quality varies.** The Geometric Unity draft was obtained from a third-party mirror, and the Howard material is mostly informal. [`sources/`](sources/) records exactly what was and was not obtained.
- **Reconstruction choices are explicit.** Where a source admits more than one reading, the choice is named and the alternatives are reported.
- **How this was produced.** The analysis and code were produced with an AI coding assistant (Claude Code), directed by the repository owner. Every claim is tied to code, a proof or a citation so that it can be checked independently of how it was produced.

## Reporting an error

Please open an issue that names the file, the claim ID or line, and, if possible, a failing command. Confirmed errors are fixed and logged in [`reports/CORRECTIONS.md`](reports/CORRECTIONS.md); earlier statements are never silently rewritten.

The research charter is in [`RESEARCH_PROMPT.md`](RESEARCH_PROMPT.md). No license file has been added yet.
