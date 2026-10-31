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

### Main results for review

Each result is stated under the hypotheses given in the linked note, with code or Lean in the repository. Internal reviews are AI-assisted referee reports in [reports/reviews/](reports/reviews/); none of these results has had external peer review.

- **Chain Theory, α Letter.**
  - The 172-record construction and its value $1/\alpha=137.035\,999\,189\,566\ldots$ are reproduced independently (F-CT-1).
  - The uniqueness theorem holds only under one far-end conjugation convention. Under the conjugate convention a second configuration qualifies, with $1/\text{value}=31.81$ (F-CT-3, Lean). Requiring the identification algebra to contain the matchings $U\otimes\bar U$ (acting trivially off the common support) forces the Letter's convention and its direct-sum trace. Within the Letter's own definition of that algebra, this one requirement is equivalent to the uniqueness theorem's selection; the Letter does not state it ([conjugation](formalization/chain-theory/alpha/conjugation.md)).
  - The record set is not covariant under the symmetry group of the joint state. A covariant repair gives the same value (F-CT-10).
  - The weights have an exact operator form, but its qutrit projector restates the end-attribution rule. No construction found in the declared search families avoids that rule (F-CT-11).
  - The three commutator records enter the 172 by an exception clause. One rule replaces it, conditional on two readings the Letter does not state. With a covariant exchange it gives the Letter's sectors and value (commutator-type records $S_x,D_y,S_z$) only if the read qubits are identified by $\sigma_y$, not by the transpose ([records-uniform](formalization/chain-theory/alpha/records-uniform.md)).
  - In a declared family of 46 080 rules and a mechanically generated family of 332 308 mutants of the author's Lean definitions, the Letter's value is the only one within relative $10^{-4}$. The accompanying reference counts are not calibrated significance (F-CT-4).
  - Details: [audit](formalization/chain-theory/alpha/audit.md).
- **Chain Theory, compendium.** Of 64 Lean-tied theorem rows, 1 is matched clause by clause, 7 partly and 56 not at all (1/8/55 and 0/7/57 under the alternative counting rules). 560 of the 618 table-check lemmas only confirm that an entered string occurs. Details: [report](reports/chain-theory-compendium.md).
- **Geometric Unity, algebra.**
  - The observerse signature is (7,7) exactly under the conditions of A4 and A5.
  - The spinor content matches one Standard Model generation with $\nu_R$, plus a further generation in the 144.
  - Lemmas 6.5 and 7.2 are false as printed, and repairable by one global sign choice.
  - The field equation (9.7) is the variational equation of no functional, for every admissible Shiab operator (D4).
  - A generic self-adjoint Shiab operator gives an elliptic complex in (14,0) for $\mathfrak u(128)$ with $\kappa\neq0$; the draft's own form (8.1) does not (H5).
  - Details: [report](reports/geometric-unity-phase1.md).
- **Geometric Unity, 14→4 reduction.**
  - Theorem F4(c), the Ricci blocks of the Z1 metric along a section, holds for every $n$ (K).
  - In the invariant odd sector, with the even displacement slaved to it in the vector sector ($\kappa\neq0$, $\det M\neq0$), $\Upsilon=0$ in Z2 has no solution on any open subset of $Y$ meeting the section, for every member with $(x,y)\neq0$ in (14,0) and (7,7), and none on any open subset at all; for $(x,y)=0$ the only solution is $T=0$, for every metric (round 6: J5; J5d; round 8: [L](formalization/geometric-unity/L-nonradial-z2.md); round 9: [M](formalization/geometric-unity/M-z2-remaining-members.md); round 12: [P](formalization/geometric-unity/P-z2-last-exception.md), which closes the last exception: off the section the fibre equation has no local solution for any base curvature). Round 8 was refereed in [REVIEW6-GU](reports/reviews/REVIEW6-GU-2026-10-30.md): correct as stated, with two summary overstatements corrected. Round 9 was refereed in [REVIEW7-GU](reports/reviews/REVIEW7-GU-2026-10-30.md): correct, with one side claim corrected (the member $c^2=\tfrac{11}6$ in (14,0) is settled on every open set) and summary wording corrected. Round 12 has not yet been internally reviewed.
  - In Z1 (round 10, [N](formalization/geometric-unity/N-z1-neighbourhood.md)), for $W_A+W_C\circ\mathrm{Ad}_\Gamma$, $W_A$ and the other members treated in round 8: no solution near the section in (7,7); in (14,0) only $\lambda=\tfrac15$ over a Ricci-flat base is possible, and for analytic fields exactly the flat-base solutions with constant odd part and zero even part. For $W_A+W_C\circ\mathrm{Ad}_\Gamma$ and $W_A$ these exist when $c_V\kappa^2\ge9360$, and they solve the field equation on all of $Y$, with $\Lambda_4=0$. One step (N6(f)) is computer-assisted. Round 10 was refereed in [REVIEW8-GU](reports/reviews/REVIEW8-GU-2026-10-30.md): no errors; N6(f) completed and summary wording corrected. Round 13 ([Q](formalization/geometric-unity/Q-z1-remaining-members.md)) extends this to every member of the family with $(x,y)\neq0$ (the cone $xz\pm yw=0$ and $x=\pm y$ in (7,7)); for $(x,y)=0$ only $T=0$, for every metric; the flat-base solutions exist for every (14,0) member iff $c_V\kappa^2\ge9360(x'^2+y'^2)$. Besides N6(f), round 13 uses exact computations (a resultant certificate and one pairing constant). Round 13 was refereed in [REVIEW10-GU](reports/reviews/REVIEW10-GU-2026-10-31.md): one edge case corrected (at $a^2=\tfrac{125}{576}$ the existence threshold is strict).
  - The variation in the section ℷ (round 11, [O](formalization/geometric-unity/O-gimel-variation.md)) needs a weight on the non-compact fibre. In every signature that gives (7,7), no natural weight function converges; the convergent natural weights are measures on the scaling ray through the section. With them, and with $\Upsilon=0$ kept as the draft's equation on $Y$, the ℷ-equation in the invariant sector is exactly 4D vacuum GR, with $\Lambda\to\Lambda_A$ (the reduced-action value of F5) as the displacement goes to zero; in Z1 this is $r_\mu\Lambda_A$, with a factor of the weight equivalent to a change of the free fibre scale, and equal to $\Lambda_A$ for the weight localised at the section. This is the ℷ-equation alone: in this sector $\Upsilon=0$ itself holds only for the flat-base solutions above. Using the weight also for the $Y$-fields removes the Einstein term. Round 11 was refereed in [REVIEW9-GU](reports/reviews/REVIEW9-GU-2026-10-30.md): one minor error in a side claim, off-shell qualifiers added.
- **Howard arithmetic.** No ring, semiring or near-ring on ℕ, ℤ, ℚ or ℝ realises $1\otimes1=2$ together with $2\otimes2=4$ (Lean). The lattice of the seven axioms is fully decided. Details: [report](reports/howard-phase1.md).

The open, well-posed problems are listed in order of priority in [SYNTHESIS §4](reports/SYNTHESIS.md#4-next-steps-in-priority-order) and, item by item, in [gaps/](gaps/).

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

Use Lean **4.35.0** for the recorded checks. The upstream compendium pins **4.32.1**; the documented 4.35.0 compatibility build applies five recorded proof edits to scratch copies. Its build result applies to that configuration. The selected Python suite does not certify a fresh Lean build or reproduce every manuscript result. All three manuscript sources compiled with the desktop editor's compiler on October 30, 2026, before the pre-submission audit edits of the same day; the edited sources have not been recompiled. Compilation is separate from mathematical verification.

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
