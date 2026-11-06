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
| Geometric Unity computations | 58 |
| **Self-contained total** | **68** |
| Optional pinned-source checks below | 6 |
| **Total with source checks** | **74** |

Inspect the actual command list without executing it:

```bash
python -X utf8 scripts/verify.py --recent --list
```

The [CI workflow](../.github/workflows/verify.yml) runs the 68 self-contained commands on Linux and Windows with Python 3.12 and 3.13. Source-dependent checks additionally require the exact checkout below.

## Pinned upstream sources

Prepare a pristine LF checkout of the recorded Chain Theory snapshot:

```bash
git clone --config core.autocrlf=false https://github.com/alexwg/chain-theory.git ../chain-theory
git -C ../chain-theory checkout --detach a28db8998bce5bd5ef11195c05f13c86cd9dac5c
python -X utf8 scripts/verify.py --recent --upstream ../chain-theory
```

Keep that checkout free of modified, untracked and ignored files, including build artifacts. The verifier checks the commit and raw source hashes; CRLF conversion is rejected. These commands read upstream material as data. The Lean-coverage classifier uses recorded build and axiom logs.

## Completed validation

[BI](../formalization/geometric-unity/BI-h6-cubic-and-mixed.md) combines a universal four-contraction proof with exact coefficient controls and a Gram-rank obstruction on the whole cubic module. Its separate mixed-family helper checks complete diagonal quadratics and a positive-square identity. [BJ](../formalization/geometric-unity/BJ-h6-orbit-averages.md) is a geometric/analytic proof for finite-measure orbit averages; its standalone helper checks only the sphere constants, not the global theorem. Reviews72–73 in the [review records](../reports/reviews/) distinguish hand arguments, finite calculations and shared contributions.

The preceding [BG](../formalization/geometric-unity/BG-h6-fixed-values.md)/[BH](../formalization/geometric-unity/BH-h6-degree2.md) kernel and degree-two results retain their validation in Reviews70–71. [BE](../formalization/geometric-unity/BE-affine-h2.md)/[BF](../formalization/geometric-unity/BF-h6-curvature-module.md) retain Reviews66–69. The 135-dimensional contained kernel is not a full ambient-rank calculation, and none of these checks confirms the proposed physics.

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
