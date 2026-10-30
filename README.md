# universalunifyingmath

This repository formally investigates three proposed unification frameworks:

1. **Chain Theory** (A. D. Wissner-Gross)
2. **Geometric Unity** (E. Weinstein)
3. **Terrence Howard's "Terryology"**

The investigation is hypothesis-neutral. For each framework we reconstruct the strongest precise mathematics it supports, then grade what is proven, conjectured, inconsistent or open. Mathematical validity and physical validity are reported separately. The charter is in [`RESEARCH_PROMPT.md`](RESEARCH_PROMPT.md).

## Where to start

**Overview and graded ledger:** [`reports/SYNTHESIS.md`](reports/SYNTHESIS.md).

| Framework | Phase summary | Main formal results |
|---|---|---|
| Chain Theory — α(0) Letter | [`reports/chain-theory-alpha.md`](reports/chain-theory-alpha.md) | [`formalization/chain-theory/alpha/audit.md`](formalization/chain-theory/alpha/audit.md) |
| Chain Theory — compendium | [`reports/chain-theory-compendium.md`](reports/chain-theory-compendium.md) | `formalization/chain-theory/compendium/` |
| Geometric Unity | [`reports/geometric-unity-phase1.md`](reports/geometric-unity-phase1.md) | `formalization/geometric-unity/` |
| Howard | [`reports/howard-phase1.md`](reports/howard-phase1.md) | `formalization/howard/` |

## Manuscript drafts (LaTeX, uncompiled in this environment)

| Draft | Path |
|---|---|
| Comment on the α(0) Letter: selection conventions and trial factor | `proofs/latex/ct-alpha-comment/comment.tex` |
| The algebraic core of Geometric Unity | `proofs/latex/gu-note/note.tex` |
| Which arithmetics satisfy "1×1=2"? | `proofs/latex/howard-note/note.tex` |

## Layout

```
sources/         annotated bibliographies, version-pinned (sha256)
claims/          Claim Registers  (IDs CT-…, GU-…, TH-…)
gaps/            Open Specification Gaps
patterns/        conjecture catalogue (exploration, not results)
formalization/   theorem/proof write-ups per framework
proofs/lean/     Lean 4 files
computation/     scripts + recorded outputs (each folder has a README with run commands)
metrics/         per-result metric records (YAML)
reports/         phase summaries; CORRECTIONS.md (log of corrections and retractions)
```

## Conventions

- Every number stated in `formalization/` or `reports/` is either printed by a committed script, whose output is committed beside it, or cited.
- Imported theorems are listed in the relevant `sources/` file.
- Corrections are appended to `reports/CORRECTIONS.md`, never silently edited.
