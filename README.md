# Universal Unifying Math

An open investigation into the rules behind mathematical and physical structure: the definitions, relationships, algorithms and mechanisms that might serve as its underlying **“source code.”** Finding such a structure remains an open research aim.

The first case studies are Chain Theory, Geometric Unity, and arithmetic and geometry associated with Terrence Howard. We reconstruct explicit claims, test their consequences, and record both supporting results and failures. The [research charter](RESEARCH_PROMPT.md) sets the method; the [synthesis](reports/SYNTHESIS.md) gives the current findings and research frontier.

**Status — November 2, 2026:** active research with reproducible computations, Lean proofs and manuscript drafts. Internal reviews include AI-assisted work. External peer review remains outstanding. Revisions are tracked in the [correction log](reports/CORRECTIONS.md).

## Current findings

| Investigation | What the work establishes | What remains open |
|---|---|---|
| **Chain Theory: α(0) Letter** | Independently reproduces the 172 records and candidate inverse $137.035\,999\,189\,566\ldots$. Classifies convention changes and covariant record constructions. | A principle fixing the weight attribution and other defining inputs; the physical identification with electromagnetic coupling. [Report](reports/chain-theory-alpha.md) |
| **Chain Theory: compendium 0.1** | Audits statement dependencies, numerical identifications and Lean coverage. Of 64 Lean-tied theorem rows, 1 matches fully and 7 partly under the documented counting rule. | Ambiguous premise references, missing identification maps, and the gap between finite checks and the physical claims attached to them. [Report](reports/chain-theory-compendium.md) |
| **Geometric Unity: 2021 draft** | Reconstructs signatures, representations and field equations. Establishes a gauge-sign repair, variational limits, and scoped reduction results, including singular-field rigidity and an obstruction for declared quadratic metrics with rank-one Ricci at the reference point. | Remaining Z1 and singular-field branches, more general odd fields, and choices needed to define the full physical system. [Report](reports/geometric-unity-phase1.md) |
| **Howard: specified arithmetic and geometry** | Classifies explicit arithmetic reconstructions and checks numerical and geometric statements. With ordinary addition, no ring, semiring or near-ring on ℕ, ℤ, ℚ or ℝ satisfies both $1\otimes1=2$ and $2\otimes2=4$. | The primary definitions: the book and full relevant patent texts were not obtained. Other consistent operations exist under weaker axioms. [Report](reports/howard-phase1.md) |

The α candidate is the **only distinct value within 1σ** in each of the hand-declared and mechanically generated rule families. These bounded comparisons describe sensitivity to the chosen rules; their illustrative density reference counts are not calibrated significance. The candidate itself fixes every digit and can be compared with future measurements. See the [α audit](formalization/chain-theory/alpha/audit.md).

For Geometric Unity, projection to a section, the full field equation on an open set, and variation of a weighted action are distinct tests. The strongest reduction conclusions require a specified metric, operator family and field sector. The [synthesis](reports/SYNTHESIS.md#geometric-unity) keeps these hypotheses and the singular-threshold exception together.

## Read the evidence

- [Synthesis and open questions](reports/SYNTHESIS.md): a current account across the case studies.
- [Formalizations](formalization/): precise definitions, hypotheses and proofs.
- [Claim registers](claims/) and [gaps](gaps/): attributed statements and unresolved inputs, with stable IDs.
- [Computations](computation/) and [Lean proofs](proofs/lean/): implementations, recorded outputs and axiom reports.
- [Source records](sources/): versions, hashes, citations and acquisition limits.
- [Internal reviews](reports/reviews/) and [corrections](reports/CORRECTIONS.md): checks, objections and subsequent changes.
- [Manuscript drafts](proofs/latex/): proposed written outputs; their scope and status are listed in the [synthesis](reports/SYNTHESIS.md#3-manuscript-drafts).

A mathematical result applies to its stated definitions and hypotheses. A finite search applies to its declared family. A Lean proof establishes its formal statement. Numerical agreement requires an explicit identification with an observable before it supports a physical claim.

**Source limits matter.** Chain Theory uses pinned public sources. The Geometric Unity draft came from a recorded third-party mirror and was not verified against the author's host. The Howard investigation relies partly on a machine transcript not checked against audio, plus patent metadata and summaries. These limits are recorded alongside the relevant claims.

## Reproduce selected checks

The selected suite supports Python 3.12 and 3.13 with pinned dependencies. Create an environment from the repository root:

```bash
python -m venv .venv
```

Activate it with `source .venv/bin/activate` on Linux/macOS or `.\.venv\Scripts\Activate.ps1` in PowerShell, then run:

```bash
python -m pip install -r requirements.txt
python -X utf8 scripts/verify.py --recent
```

The verifier runs assertions and regression checks without overwriting recorded outputs. `--recent` adds the newer Chain Theory record checks and independent GU checks. The [CI workflow](.github/workflows/verify.yml) runs the selected suite on Linux and Windows with both Python versions.

For source-dependent checks, prepare a pristine checkout of the exact Chain Theory snapshot:

```bash
git clone --config core.autocrlf=false https://github.com/alexwg/chain-theory.git ../chain-theory
git -C ../chain-theory checkout --detach a28db8998bce5bd5ef11195c05f13c86cd9dac5c
python -X utf8 scripts/verify.py --recent --upstream ../chain-theory
```

Keep that checkout free of modified, untracked and ignored files, including build artifacts. The alpha checks require canonical LF bytes and retain strict SHA-256 pins. These checks read the upstream source as data; the Lean-coverage classifier uses recorded build and axiom logs.

### Full reproduction

Each area documents its commands, dependencies, outputs and runtimes:

| Area | Instructions |
|---|---|
| Chain Theory α(0) | [Computation guide](computation/chain-theory/alpha/README.md) · [Lean guide](proofs/lean/chain-theory/README.md) |
| Chain Theory compendium | [Computation guide](computation/chain-theory/compendium/README.md) |
| Geometric Unity | [Computation guide](computation/geometric-unity/README.md) |
| Howard | [Computation guide](computation/howard/README.md) |

The recorded formal checks use **Lean 4.35.0**. The upstream compendium pins 4.32.1; its documented 4.35.0 compatibility build uses five recorded proof edits in scratch copies. The selected Python suite does not rebuild Lean, compile manuscripts, run long searches, or cover every research result. Optional dependencies are in [requirements-optional.txt](requirements-optional.txt).

The full transport search requires the `fork` multiprocessing method; recorded runs used Linux. The compendium build driver and shell pipeline require Linux or WSL. The selected Python checks run natively on Windows with UTF-8 output.

## Contributing

Name the question, source version, definitions and assumptions. Provide a proof, counterexample or reproducible command, and state the limits of the result. For physical comparisons, specify the observable, units and uncertainty; for significance claims, specify the probability and selection models.

Error reports should identify the claim ID or file location and the statement that fails. Confirmed corrections update the active account and enter the [correction log](reports/CORRECTIONS.md); historical reviews remain available.

No license has been selected for this repository. Upstream materials retain their respective terms, recorded in [sources](sources/).
