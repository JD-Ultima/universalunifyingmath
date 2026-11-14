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
| Geometric Unity computations | 73 |
| **Self-contained total** | **83** |
| Optional pinned-source checks below | 6 |
| **Total with source checks** | **89** |

Inspect the actual command list without executing it:

```bash
python -X utf8 scripts/verify.py --recent --list
```

The [CI workflow](../.github/workflows/verify.yml) runs the 83 self-contained commands on Linux and Windows with Python 3.12 and 3.13. Source-dependent checks additionally require the exact checkout below.

## Pinned upstream sources

Prepare a pristine LF checkout of the recorded Chain Theory snapshot:

```bash
git clone --config core.autocrlf=false https://github.com/alexwg/chain-theory.git ../chain-theory
git -C ../chain-theory checkout --detach a28db8998bce5bd5ef11195c05f13c86cd9dac5c
python -X utf8 scripts/verify.py --recent --upstream ../chain-theory
```

Keep that checkout free of modified, untracked and ignored files, including build artifacts. The verifier checks the commit and raw source hashes; CRLF conversion is rejected. These commands read upstream material as data. The Lean-coverage classifier uses recorded build and axiom logs.

## Completed validation

[BV](../formalization/geometric-unity/BV-h6-triangular-positive-obstruction.md) excludes the full 952-dimensional triangular closed family. The selected check verifies the receiver, positivity of all stored matrices, and exact rational and modular reconstruction on one block. Use `--full` on `BV_h6_an_positive.py` to reconstruct all sixteen blocks.

[BW](../formalization/geometric-unity/BW-h6-triangular-orbit-averages.md) excludes fixed finite signed averages of rotated triangular closed fields. The diagonal limit is analytic and uses the BR and BV certificates.

[BX](../formalization/geometric-unity/BX-h4-curvature-and-closed-fields.md) verifies the exact dimensions and curvature representation types for the pure H4 branch. Its selected check runs the complete rational upper and modular lower rank bounds and every triangular null vector.

[BY](../formalization/geometric-unity/BY-h4-low-isotype-obstruction.md) verifies the 24-dimensional grade-five exclusion. The selected check verifies both modular receiver and principal-minor ranks at two primes.

[CA](../formalization/geometric-unity/CA-h4-full-analytic-lifts.md) checks the full H4 analytic closed lift span using the BX character. The selected check verifies all 28 real triangular witnesses, 57 extremals and the 22 nonzero pairing determinants at two primes.

[CF](../formalization/geometric-unity/CF-h6-fourth-linear-jet.md) uses an exact rational coefficient comparison to verify a fourth-order H6 closure jet at a nonzero value in the remaining type $(5,5)$.

[CB](../formalization/geometric-unity/CB-h4-certified-partial-value.md) verifies the real H4 partial point value, the complete grade-one and grade-nine detectors, and its separating grade-five coefficient. The selected check reconstructs all invariant grade-five terms.

[CD](../formalization/geometric-unity/CD-hook-clifford-identities.md) compares two exact Clifford implementations and verifies the exterior identities and operator sum of squares on rational controls.

[CG](../formalization/geometric-unity/CG-hook-quadratic-separation-limits.md) verifies the exact invariant receiver counts and moment witnesses that limit quadratic separation on the full H4 curvature kernel and the 6,534-dimensional H6 liftable submodule.

[CH](../formalization/geometric-unity/CH-h4-second-order-stress.md) verifies the exact second-order stress bound on the CB box and computes the inertia of the curvature source on all H4 curvature types.

[CI](../formalization/geometric-unity/CI-h4-grade-five-cartan-support.md) verifies the exact highest-weight exclusions extending BY to a 512-dimensional H4 support. The selected check compares all 136 necessary coefficients in two exact Clifford engines.

[CJ](../formalization/geometric-unity/CJ-h6-universal-stress-gradient.md) verifies the H6 gradient estimate using two positive invariant curvature pencils across all 15 types.

[BU](../formalization/geometric-unity/BU-h6-second-jet-obstruction.md) removes the full 14-dimensional type $(2,0)$ from every closed H6 field. The verifier reconstructs the complete seven-dimensional equivariant jet test and proves the cokernel obstruction in exact rational arithmetic.

[BT](../formalization/geometric-unity/BT-h6-triangular-character-obstruction.md) excludes each of two 73-dimensional closed character spaces. The selected check rebuilds the closed basis, tests independent Clifford coordinates and ranks 16 selected quadratic columns. The full 2,701-column certificate requires `--full` and is outside the selected suite.

[BS](../formalization/geometric-unity/BS-h6-certified-partial-value.md) certifies a real unit point value satisfying complete grade one and grade thirteen, together with the invariant grade-five and grade-nine projections. The default check verifies the rational contraction, detector scope, and separation from a full equation. The helper's `--full` mode reconstructs all affine tensor and polynomial data and takes several minutes.

[BR](../formalization/geometric-unity/BR-h6-closed-field-lifts.md) constructs a 952-dimensional family of closed fields. The default helper checks the seven real closure witnesses, the orbit minors for the 6,534-dimensional lifted submodule and the exact 73-dimensional sign block. The selected suite also checks the affine character identity. The helper's `--full-exact` mode reconstructs every sign block.

[BQ](../formalization/geometric-unity/BQ-h6-invariant-moments.md) evaluates two invariant pairings on all 15 isotypic types. The selected suite reconstructs the receivers and three representative blocks, and checks the exact positivity and energy certificates on every listed block. The helper's `--full` option recalibrates every matrix.

[BO](../formalization/geometric-unity/BO-h6-full-curvature-kernel.md) computes the full 6,834-dimensional kernel from 185 modular ranks and 18 rational extremal vectors. It also verifies a negative grade-nine witness in the physical frame. [BP](../formalization/geometric-unity/BP-h6-cartan-grade-nine.md) transports four of the extremal vectors and checks their nonzero Cartan-square outputs. [Review79](../reports/reviews/REVIEW79-GU-BO-2026-11-14.md) and [Review80](../reports/reviews/REVIEW80-GU-BP-2026-11-14.md) record independent checks.

[BN](../formalization/geometric-unity/BN-h6-invariant-grade-nine.md) excludes the whole 135-dimensional BG module by a positive invariant grade-nine projection. Its helper constructs the receiver and calibrates all three summands from public inputs. The optional `--full-nested` flag on `computation/geometric-unity/BN_h6_invariant_grade9.py` checks all 3,429 ordered coefficients with the separate nested Clifford evaluator. [Review78](../reports/reviews/REVIEW78-GU-BN-2026-11-14.md) records the shared conventions and the complete comparison.

[BL](../formalization/geometric-unity/BL-mixed-grade-one-existence.md) supplies the complete cubic/quintic cross formula and an exact real mixed solution of the grade-one equations. Its three helpers calibrate ordered Clifford coefficients, reduce all 196 finite-family scalarity equations to four quadrics, and verify a rational contraction certificate. The pure-cubic full-symmetry completion and contraction theorem are hand arguments. [BM](../formalization/geometric-unity/BM-mixed-grade-nine-obstruction.md) reconstructs the actual six H6 tensors and proves that candidate fails grade nine throughout its certified box. [Review76](../reports/reviews/REVIEW76-GU-BL-2026-11-07.md) and [Review77](../reports/reviews/REVIEW77-GU-BM-2026-11-07.md) record fresh reproduction and its limits. Other scalarity zeros and full field solutions remain open.

[BK](../formalization/geometric-unity/BK-quintic-gram-and-symmetry.md) gives a universal four-contraction quintic identity, exact coefficient controls and the complete grade-one Gram criterion. Its second helper reconstructs two full symmetry classes, checks every Gram coefficient and expands four rational polynomial certificates to one. The constrained-isometry equivalence and scope transfer are hand proofs. Reviews74–75 in the [review records](../reports/reviews/) distinguish these arguments from finite controls and disclose shared contributions.

The [BI](../formalization/geometric-unity/BI-h6-cubic-and-mixed.md)/[BJ](../formalization/geometric-unity/BJ-h6-orbit-averages.md) cubic, mixed and global orbit results retain Reviews72–73. Their scalar sphere helper does not verify the global analytic theorem. [BG](../formalization/geometric-unity/BG-h6-fixed-values.md)/[BH](../formalization/geometric-unity/BH-h6-degree2.md) retain Reviews70–71, and [BE](../formalization/geometric-unity/BE-affine-h2.md)/[BF](../formalization/geometric-unity/BF-h6-curvature-module.md) Reviews66–69. None computes the full ambient curvature kernel or confirms the proposed physics.

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
