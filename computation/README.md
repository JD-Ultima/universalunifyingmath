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
| Geometric Unity computations | 55 |
| **Self-contained total** | **65** |
| Optional pinned-source checks below | 6 |
| **Total with source checks** | **71** |

Inspect the actual command list without executing it:

```bash
python -X utf8 scripts/verify.py --recent --list
```

The [CI workflow](../.github/workflows/verify.yml) runs the 65 self-contained commands on Linux and Windows with Python 3.12 and 3.13. Source-dependent checks additionally require the exact checkout below.

## Pinned upstream sources

Prepare a pristine LF checkout of the recorded Chain Theory snapshot:

```bash
git clone --config core.autocrlf=false https://github.com/alexwg/chain-theory.git ../chain-theory
git -C ../chain-theory checkout --detach a28db8998bce5bd5ef11195c05f13c86cd9dac5c
python -X utf8 scripts/verify.py --recent --upstream ../chain-theory
```

Keep that checkout free of modified, untracked and ignored files, including build artifacts. The verifier checks the commit and raw source hashes; CRLF conversion is rejected. These commands read upstream material as data. The Lean-coverage classifier uses recorded build and axiom logs.

## Completed validation

[BG](../formalization/geometric-unity/BG-h6-fixed-values.md) reconstructs the complete stabilizer-fixed H6 curvature kernel and checks its nonlinear positive-square obstruction. [BH](../formalization/geometric-unity/BH-h6-degree2.md) covers every value in its degree-two module using a complete diagonal section and a hand spectrum lemma. Three new exact helpers retain full curvature channels, all quadratic cross terms and targeted comparisons with a separate Clifford implementation. Reviews70–71 in the [review records](../reports/reviews/) assess the computation/proof boundary and disclose shared contributions. These are restricted mathematical exclusions, not a classification of arbitrary fields.

The preceding [BE](../formalization/geometric-unity/BE-affine-h2.md) all-positive-parameter affine H2 exclusion and [BF](../formalization/geometric-unity/BF-h6-curvature-module.md) contained-module calculation retain their recorded validation in Reviews66–69. BG extends BF's 121-dimensional bound to 135; neither computes the whole ambient curvature kernel.

Earlier validation remains in the per-result notes and reviews:

| Work | Evidence |
|---|---|
| BA closure; BC/BD affine hook exclusions | [BA](../formalization/geometric-unity/BA-hook-closure.md), [BC](../formalization/geometric-unity/BC-affine-hooks.md), [BD](../formalization/geometric-unity/BD-affine-h4.md) and Reviews61–65 in the [review records](../reports/reviews/) |
| AZ positive-parameter first-hook classification | [AZ proof](../formalization/geometric-unity/AZ-first-hook-classification.md), exact polynomial cover and [Review59](../reports/reviews/REVIEW59-GU-AZ-2026-11-03.md)/[Review60](../reports/reviews/REVIEW60-GU-AZ-2026-11-03.md) |
| AX derivative and curvature certificates; curved-base extension | [AX proof](../formalization/geometric-unity/AX-hook-curvature.md), Reviews55–58 in [review records](../reports/reviews/) |
| AW single-Lambda5 classification; preceding AP–AV reductions | [AW proof](../formalization/geometric-unity/AW-lambda5-conic-classification.md), linked controls and Reviews46–54 |
| AL endpoint equations and AM global models | [AL](../formalization/geometric-unity/AL-opposite-lambda5-system.md), [AM](../formalization/geometric-unity/AM-global-flat-models.md), Reviews42–43; recorded 50/50 local checks with pinned sources |
| AH algebraic constraints | [AH proof](../formalization/geometric-unity/AH-lambda5-gram.md), 144 controls and [PR #9](https://github.com/JD-Ultima/universalunifyingmath/pull/9) |
| HH Ricci-zero obstruction | [HH proof](../formalization/geometric-unity/HH-ricci-zero.md) and [Review36](../reports/reviews/REVIEW36-GU-HH-2026-11-02.md) |
| Reused read-cap theorem | [Fresh Lean rebuild and adapter](../reports/chain-theory-readcap-reuse.md); outside the Python suite |

Command counts are not theorem counts. Some checks share an underlying engine. Passing confirms each implemented statement, including expected counterexamples; it does not supply external peer review or physical confirmation. Recorded Lean builds and manuscript checks have separate scope.

## Full reproduction and exclusions

| Area | Instructions |
|---|---|
| Chain Theory α Letter | [Computations](chain-theory/alpha/README.md) · [Lean](../proofs/lean/chain-theory/README.md) |
| Chain Theory compendium | [Computations and build audit](chain-theory/compendium/README.md) |
| Geometric Unity | [Computations](geometric-unity/README.md) |

Long searches, optional SDP checks, manuscript compilation and fresh Lean builds run separately. Additional dependencies are in [requirements-optional.txt](../requirements-optional.txt).

Earlier recorded formal checks use Lean 4.35.0; the read-cap reuse check pins Lean 4.35.0-rc4 and its own Mathlib/dependency revisions. The upstream compendium pins 4.32.1; its documented 4.35.0 compatibility build uses five proof edits in scratch copies. The selected suite does not repeat that build. See the [coverage audit](../formalization/chain-theory/compendium/10-lean-table-rows.md) for what the formal statements establish.

The transport search requires the fork multiprocessing method; recorded runs used Linux. The compendium build driver and shell pipeline require Linux or WSL. Selected Python checks run natively on Windows with UTF-8 output. Manuscript scope and recorded build dates are listed in the [synthesis](../reports/SYNTHESIS.md#3-manuscript-drafts).
