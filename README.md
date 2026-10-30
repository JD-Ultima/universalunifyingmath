# Universal Unifying Math

**Evidence-led exploration of mathematical structures and physical models.**

This project investigates which explicit rules can explain observed structure, remain mathematically consistent, and produce testable predictions. Our goal is to find the underlying **“source code”**: definitions, relationships, algorithms, and mechanisms that can be stated precisely and examined through proof, computation, and measurement. The form of that underlying structure remains an open research question.

The research begins with questions. A framework's author, reputation, or agreement with an existing view does not determine its result here. Conclusions apply to the specified source, definitions, assumptions, and tests. Both supporting and conflicting evidence are recorded, and conclusions are revised when the evidence changes.

**Status — October 30, 2026:** active exploratory research. The repository contains source reconstructions, computations, formal proofs, internal review records, and manuscript drafts. Internal reviews include AI-assisted analysis; external peer review remains outstanding. The [research charter](RESEARCH_PROMPT.md) defines the method, and the [correction log](reports/CORRECTIONS.md) records revisions.

## Research method

Every investigation separates five things:

1. **Source statement:** what an identified, versioned source actually says.
2. **Reconstruction:** the definitions and interpretive choices used to make the statement precise.
3. **Verification:** the proof, computation, counterexample, or comparison supporting a result.
4. **Open question:** what the available definitions or evidence do not settle.
5. **Physical test:** a prediction with a specified observable, units, uncertainty, and failure condition.

Numerical agreement is recorded at its stated precision. Statistical significance requires a justified probability model and an account of selection choices. A finite search establishes a result within its declared family. A Lean proof establishes its formal statement under its definitions and hypotheses. Physical interpretations require their own evidence.

Exploration is welcome when it produces explicit questions and a route to testing them. Unverified interpretations are kept out of factual conclusions. Conclusions require traceable mathematical or empirical evidence.

## Current investigations

These are the initial case studies; the research scope can expand as suitable sources and questions become available.

| Investigation | Material examined | Results and open questions |
|---|---|---|
| **Chain Theory — fine-structure constant** | A. D. Wissner-Gross's α(0) Letter v1.0.0, supplement, Lean definitions, and record table | Reproduces the 172-record construction and its numerical value. Examines conjugation conventions, covariance, weight assignments, and sensitivity to alternative rules. The physical identification and a calibrated statistical significance remain unresolved. [Report](reports/chain-theory-alpha.md) · [Detailed analysis](formalization/chain-theory/alpha/audit.md) |
| **Chain Theory — compendium** | Release 0.1.0 at the recorded source commit | Examines statement dependencies, numerical identifications, and the content covered by Lean declarations. Coverage and dependency counts include documented interpretation rules. [Report](reports/chain-theory-compendium.md) · [Analysis](formalization/chain-theory/compendium/) |
| **Geometric Unity** | E. Weinstein's 2021 working draft, obtained from a recorded mirror | Reconstructs signatures, spinor representations, gauge identities, field equations, ellipticity, and dimensional reduction under stated hypotheses. The projected and full field equations have different constraints; matching reduced algebraic equations does not establish identical solution sets. [Report](reports/geometric-unity-phase1.md) · [Analysis](formalization/geometric-unity/) |
| **Arithmetic and geometry associated with Terrence Howard** | Identified interview material, a third-party transcript, and available patent metadata or summaries | Tests explicit arithmetic reconstructions, numerical statements, and geometric models. Results depend on the chosen operations and carrier. The primary book and full patent texts were unavailable in the recorded investigation. [Report](reports/howard-phase1.md) · [Analysis](formalization/howard/) |

See the [research synthesis](reports/SYNTHESIS.md) for the result ledger, [claim registers](claims/) for individual statements, [source records](sources/) for provenance, and [specification gaps](gaps/) for unresolved inputs.

## Reproduce selected checks

The recorded research environment used Python 3.13 and Lean 4.35.0. The selected Python verification suite targets Python 3.12 and 3.13, with pinned dependencies in [requirements.txt](requirements.txt).

From the repository root:

~~~bash
python -m venv .venv
~~~

Activate the environment:

~~~bash
# Linux / macOS
source .venv/bin/activate
~~~

~~~powershell
# Windows PowerShell
.\.venv\Scripts\Activate.ps1
~~~

Then install dependencies and run the selected checks:

~~~bash
python -m pip install -r requirements.txt
python -X utf8 scripts/verify.py
~~~

The verifier runs script assertions and focused regression tests, including the projected-equation counterexample. It prints failures and returns a nonzero exit status if a check fails. It does not overwrite recorded outputs. The [GitHub Actions workflow](.github/workflows/verify.yml) runs this suite on Linux and Windows with Python 3.12 and 3.13.

For optional checks against the exact Chain Theory source snapshot:

~~~bash
git clone https://github.com/alexwg/chain-theory.git ../chain-theory
git -C ../chain-theory checkout a28db8998bce5bd5ef11195c05f13c86cd9dac5c
python -X utf8 scripts/verify.py --upstream ../chain-theory
~~~

Use a pristine upstream checkout with no modified, untracked, or ignored files, including build artifacts. These additional checks inspect the source as data. The Lean-coverage classifier reads committed build and axiom logs; running it is separate from rebuilding Lean.

### Full computations and formal proofs

Each computation directory documents commands, dependencies, expected outputs, and runtimes:

- [Chain Theory α(0)](computation/chain-theory/alpha/README.md)
- [Chain Theory compendium](computation/chain-theory/compendium/README.md)
- [Geometric Unity](computation/geometric-unity/README.md)
- [Howard arithmetic and geometry](computation/howard/README.md)

[requirements-optional.txt](requirements-optional.txt) adds the recorded CVXPY and element-data dependencies for scripts outside the selected suite. Those optional checks and the long searches are not run by the default verifier.

**Platform requirements:** the full transport search requires the fork multiprocessing method; its recorded runs used Linux. The compendium build driver and shell pipeline require Linux, including WSL, for process accounting and shell tools. Both Python tools reject unsupported platforms before beginning expensive work. Other selected Python checks run natively on Windows; the UTF-8 flag handles mathematical output consistently.

The [Lean instructions](proofs/lean/chain-theory/README.md) describe the Chain Theory checks. Standalone checks include:

~~~bash
lean proofs/lean/howard/HowardArith.lean
lean proofs/lean/geometric-unity/GUTilt.lean
lean proofs/lean/geometric-unity/GStueckelberg.lean
~~~

Use Lean **4.35.0** for the recorded checks. The upstream compendium pins **4.32.1**; the documented 4.35.0 compatibility build applies five recorded proof edits to scratch copies. Its build result applies to that configuration. The selected Python suite does not certify a fresh Lean build or reproduce every manuscript result. All three manuscript sources compiled successfully with the desktop editor's compiler on October 30, 2026; compilation is separate from mathematical verification.

## Evidence labels

| Label | Meaning |
|---|---|
| **Source statement** | An attributed assertion, with its source and acquisition status recorded. |
| **Definition / reconstruction** | A specified object or interpretation used in an investigation. |
| **Proven** | A statement established under its explicit mathematical hypotheses. |
| **Computed / reproduced** | A calculation or source result checked at the stated precision and scope. |
| **Conditional** | A result that depends on an identified convention, assumption, or interpretation. |
| **Counterexample / inconsistent** | A specified statement fails under a demonstrated example or incompatible set of assumptions. |
| **Open** | A question that the available definitions, proof, or evidence do not settle. |

Established mathematical counterexamples remain part of the evidence. They concern the precise statements tested and do not decide every possible formulation of a research program.

## Repository guide

| Directory | Contents |
|---|---|
| [sources/](sources/) | Versioned bibliographies, acquisition limits, and imported results. |
| [claims/](claims/) | Attributed statements with stable claim IDs. |
| [gaps/](gaps/) | Missing definitions, inputs, and conditions needed to resolve questions. |
| [patterns/](patterns/) | Documented observations and explicitly identified research questions. |
| [formalization/](formalization/) | Definitions, hypotheses, derivations, and proofs. |
| [computation/](computation/) | Python implementations and recorded outputs. |
| [proofs/lean/](proofs/lean/) | Formal proofs and recorded axiom reports. |
| [proofs/latex/](proofs/latex/) | Manuscript sources; compilation and review status are recorded separately. |
| [metrics/](metrics/) | Result status, dependencies, and reproducibility records. |
| [reports/](reports/) | Summaries, internal reviews, and correction history. |
| [scripts/](scripts/) · [tests/](tests/) | Selected verification runner and regression checks. |

## Contributing

Contributions should identify the question, source version, definitions, assumptions, method, and evidence. Include a reproducible command or proof where applicable. Distinguish a source quotation from an interpretation, and a mathematical result from a physical prediction. Report numerical uncertainty and any choices made after inspecting a target value.

For an error report, provide the claim ID or file location, the failing statement, and a counterexample or reproduction command when available. Corrections update the active account and are appended to [CORRECTIONS.md](reports/CORRECTIONS.md). Historical reviews remain available with their context. See the [research charter](RESEARCH_PROMPT.md) for the full protocol.

No license has been selected for this repository. Upstream materials retain their respective terms, as recorded in [sources/](sources/).
