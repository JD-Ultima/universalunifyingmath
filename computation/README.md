# Verification and reproduction

The selected suite checks specified Chain Theory and Geometric Unity results. It does not cover every statement in the repository.

## Selected suite

Use Python 3.12 or 3.13. From the repository root, create an environment:

```bash
python -m venv .venv
```

Activate it with the appropriate command:

```bash
source .venv/bin/activate
```

```powershell
.\.venv\Scripts\Activate.ps1
```

Then install the pinned dependencies and run the suite:

```bash
python -m pip install -r requirements.txt
python -X utf8 scripts/verify.py --recent
```

The verifier runs assertions and regressions without overwriting recorded outputs. It checks installed dependency versions. Individual commands have a 180-second default limit.

| Group | Commands |
|---|---:|
| Focused regressions, containing 12 tests | 1 |
| Chain Theory computations | 9 |
| Geometric Unity computations | 46 |
| **Self-contained total** | **56** |
| Optional pinned-source checks below | 6 |
| **Total with source checks** | **62** |

Inspect the actual command list without executing it:

```bash
python -X utf8 scripts/verify.py --recent --list
```

The [CI workflow](../.github/workflows/verify.yml) runs the 56 self-contained commands on Linux and Windows with Python 3.12 and 3.13. Source-dependent checks additionally require the exact checkout below.

## Pinned upstream sources

Prepare a pristine LF checkout of the recorded Chain Theory snapshot:

```bash
git clone --config core.autocrlf=false https://github.com/alexwg/chain-theory.git ../chain-theory
git -C ../chain-theory checkout --detach a28db8998bce5bd5ef11195c05f13c86cd9dac5c
python -X utf8 scripts/verify.py --recent --upstream ../chain-theory
```

Keep that checkout free of modified, untracked and ignored files, including build artifacts. The verifier checks the commit and raw source hashes; CRLF conversion is rejected. These commands read upstream material as data. The Lean-coverage classifier uses recorded build and axiom logs.

## Completed validation

The [AZ classification](../formalization/geometric-unity/AZ-first-hook-classification.md) adds a standalone exact affine-geometry and polynomial-parameter certificate. Ten symbolic maximal minors cover five character-weight orbits; their common roots reduce the remaining positive parameters to three exact cases, all with full rational rank. The full parameter range follows from that exhaustive cover, not a rank sample. [Review59](../reports/reviews/REVIEW59-GU-AZ-2026-11-03.md) and [Review60](../reports/reviews/REVIEW60-GU-AZ-2026-11-03.md) assess the model, certificate and general-base field argument. The current manifest has 56 self-contained commands, 62 with optional source checks.

The [AX proof](../formalization/geometric-unity/AX-hook-curvature.md) has three standalone certificates: all **12,544 derivative columns** on the 896-dimensional hook, sixteen exact rational full-rank curvature blocks with nonzero determinant witnesses, and a separate maximal-rank certificate using a different matrix basis and modular elimination. The rational helper also checks the full constant-curvature counterexample to general injectivity. [Review55](../reports/reviews/REVIEW55-GU-AX-2026-11-03.md) and [Review56](../reports/reviews/REVIEW56-GU-AX-2026-11-03.md) record the new reproductions and scope audits. The curved-base extension AX4–AX5 is a hand proof using those unchanged certificates; [Review57](../reports/reviews/REVIEW57-GU-AX-2026-11-03.md) and [Review58](../reports/reviews/REVIEW58-GU-AX-2026-11-03.md) assess its full-metric rescaling, fibre density and Euler contradiction. No helper or recorded output changed, and no old local computation was rerun for the extension. That milestone retained 55 self-contained commands, 61 with optional source checks.

The AW helper passes **34** new exact controls. [Review53](../reports/reviews/REVIEW53-GU-AW-2026-11-03.md) accepts the normal-moment/conservation proof and independently reproduces the public helper; [Review54](../reports/reviews/REVIEW54-GU-AW-2026-11-03.md) records the phase/scope audit and discloses its contributor's role in the convention controls. The standalone helper checks noncommuting Clifford pairings and exact rational certificates; the universal trace and actual-interval arguments are hand proofs. That milestone contained 58 commands with optional source checks. Its focused runs did not repeat an unchanged full local suite.

The preceding [AT/AU/AV milestone](https://github.com/JD-Ultima/universalunifyingmath/pull/15) passed 13, 24 and 13 new exact controls and all eight CI jobs, each running its 51-command self-contained suite. [Review50](../reports/reviews/REVIEW50-GU-AT-2026-11-03.md), [Review51](../reports/reviews/REVIEW51-GU-AU-2026-11-03.md) and [Review52](../reports/reviews/REVIEW52-GU-AV-2026-11-03.md) retain the analytic, metric-jet and independent Clifford convention assessments.

The preceding [AR/AS milestone](https://github.com/JD-Ultima/universalunifyingmath/pull/14) passed 34 and 52 exact controls and all eight CI jobs, each running its 48-command self-contained suite. [Review48](../reports/reviews/REVIEW48-GU-AR-2026-11-03.md) and [Review49](../reports/reviews/REVIEW49-GU-AS-2026-11-03.md) retain the geometric continuation and full field-system assessments.

The preceding [AP/AQ milestone](https://github.com/JD-Ultima/universalunifyingmath/pull/13) passed 40 and 23 exact controls and all eight final-head CI jobs, each running its 46-command self-contained suite. [Review46](../reports/reviews/REVIEW46-GU-AP-2026-11-03.md) and [Review47](../reports/reviews/REVIEW47-GU-AQ-2026-11-03.md) retain the corresponding proof assessments and shared-engine limits.

The November 3 AL/AM milestone passed all **50/50** selected commands with pinned Windows Python 3.12 and the pristine LF upstream snapshot. It includes **59** exact endpoint-system controls and **12** global-model controls. [Review42](../reports/reviews/REVIEW42-GU-AM-2026-11-03.md) and [Review43](../reports/reviews/REVIEW43-GU-AL-2026-11-03.md) record their proof assessments and limits.

The preceding [AH9 milestone](https://github.com/JD-Ultima/universalunifyingmath/pull/9) passed 48 local commands and all eight final-head CI jobs before merging; its AH helper retains 144 controls. The separate [read-cap reuse verification](../reports/chain-theory-readcap-reuse.md) records a fresh Lean target rebuild and adapter check. That build is outside the Python suite.

Per-result proof arguments, imported lemmas and additional controls are documented in the [formalizations](../formalization/) and [internal reviews](../reports/reviews/). For example, [Review36](../reports/reviews/REVIEW36-GU-HH-2026-11-02.md) reproduces HH's exact rank/orbit certificates and separately checks nonunit scaling and higher-jet perturbations.

Command counts are not theorem counts. Passing tests confirms the expected result of each check, including counterexamples. Some checks share an underlying calculation engine. These records do not constitute external peer review or physical confirmation.

## Full reproduction and exclusions

| Area | Instructions |
|---|---|
| Chain Theory α Letter | [Computations](chain-theory/alpha/README.md) · [Lean](../proofs/lean/chain-theory/README.md) |
| Chain Theory compendium | [Computations and build audit](chain-theory/compendium/README.md) |
| Geometric Unity | [Computations](geometric-unity/README.md) |

Long searches, optional SDP checks, manuscript compilation and fresh Lean builds run separately. Additional dependencies are in [requirements-optional.txt](../requirements-optional.txt).

Earlier recorded formal checks use Lean 4.35.0; the read-cap reuse check pins Lean 4.35.0-rc4 and its own Mathlib/dependency revisions. The upstream compendium pins 4.32.1; its documented 4.35.0 compatibility build uses five proof edits in scratch copies. The selected suite does not repeat that build. See the [coverage audit](../formalization/chain-theory/compendium/10-lean-table-rows.md) for what the formal statements establish.

The transport search requires the fork multiprocessing method; recorded runs used Linux. The compendium build driver and shell pipeline require Linux or WSL. Selected Python checks run natively on Windows with UTF-8 output. Manuscript scope and recorded build dates are listed in the [synthesis](../reports/SYNTHESIS.md#3-manuscript-drafts).
