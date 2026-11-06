# Review72: BI cubic/mixed values and BJ orbit averages

2026-11-06. **Accepted under the stated hypotheses and boundaries.** Read
[BI](../../formalization/geometric-unity/BI-h6-cubic-and-mixed.md),
[BJ](../../formalization/geometric-unity/BJ-h6-orbit-averages.md), both metric
records, the three new helpers/outputs and their shared inputs. No unresolved
mathematical correction was found in the reviewed package.

This is an internal AI-assisted, contribution-aware review. I did not author
the cubic or mixed-family calculations, but authored BJ and its private
precursors, collaborating with root on the profile, collapse and measure
extension. I also contributed earlier X/BC/BF convention audits and secondary
Clifford checks. BJ is consequently a self-audit here;
[Review73](REVIEW73-GU-BI-BJ-2026-11-06.md) provides the other reviewer's BJ
assessment. Shared library reproduction is not an independent symbolic engine
or external peer review.

## Mathematical checkpoints

- **BI1 completeness:** the trace-free projections extend the coefficient
  tensor to ten standard-vector slots. In odd dimension, central `-I` acts
  trivially in degree ten, extending SO5 invariance to O5. I checked the cited
  [Goodman Theorem7.4/Corollary7.5, printed p.26](https://sites.math.rutgers.edu/~goodman/pub/nuslect.pdf#page=26):
  its metric-pairing span supplies the imported invariant-theory input.
  Valence counting gives exactly the four surviving contractions. Thus the
  nonsingular calibration determines a universal identity, not a sample fit.
- **BI1 operator argument:** cross traces `14I/5` and `27H/20` give the
  scalar/H2 subtraction. Full symmetrization of the trace-free remainder and
  its complement exhaust all105 self-adjoint operators. The projected
  reshuffle factors `3/10,1/10,1,-1/2` follow from these explicit formulas.
  The three nonscalar stress multipliers are nonzero. Scalar stress forces
  `L_T*L_T` scalar; rank at most5 on dimension14 forces zero, then real
  positive norm forces `T=0`. No ambient-kernel completeness is assumed.
- **BI conventions:** `D_i=2*A_i` agrees with the calibrated Phi action and
  the covariant polynomial rotation sign. The actual compact phase is
  `iB6Gamma`; the unnormalized-volume quadratic requires `detGram`, with no
  additional sign. The metric isometry, weighted hooks, full ten-term Casimir,
  norm ratio and separate ordered off-diagonal controls fix these conventions.
- **BI2 finite coverage:** the quarter-turn polynomials interpolate the odd
  frequencies1/3 and1/3/5 on their checked cyclic annihilators. They reconstruct
  the five rotated zonals with the declared normalization. Distinct modules
  and the two exact Lie-action ranks10 imply zero infinitesimal stabilizer
  for every nonzero amplitude pair. All14 diagonal moments include the full
  mixed term. The displayed positive-square difference excludes every real
  nonzero pair, including either zero-amplitude chart. Independently rotating
  the two components, or assigning an irreducible ten/four splitting, is not
  part of the argument.
- **BJ profiles and measures:** L fixes the whole seed, so `b_n` is a genuine
  smooth map on K/L. The positive parabolic and inverse-congruence action give
  exactly `B(s^T n/||s^T n||)`, without a character or changed integration
  variable. Finite total variation permits every local derivative integral.
  Tonelli with `|mu|` supplies an uncharged equator even for atomic/singular
  measures; signed cancellation is not substituted for that condition.
- **BJ limit and domains:** dominated convergence occurs in an isometrically
  identified tensor fiber and uses only pointwise quadratic continuity.
  A nonzero limiting seed multiple fails its rotated anisotropy; a zero
  multiple fails the positive baseline. The uniform positive real part of
  the complexified denominator justifies analytic integration for every
  finite-TV measure. The identity theorem applies to this globally defined
  analytic residual, not to arbitrary smooth solutions or their domains.
- **BJ4 construction:** compact orthogonality gives the121-dimensional
  evaluation right inverse and the two exact hemisphere constants. During
  review, root requested the now-included non-affine proof: P-fixity would
  force `w=alpha(w)B*`, while `alpha(cB*)=33c/128` forces `w=0`. This completes
  the opening assertion without importing the separate infinite-orbit proof.
- **Fields and scales:** both records retain true WA, compact positive real
  conventions, the pure pair and constant nonzero kappa. Derivative grade7
  cannot cancel grade1, including with variable amplitudes. Dilation and
  homothety give the common positive moment factor `cV^-7` and actual positive
  baseline `15/(2cV)+15*kappa^2/52`. BI is a conditional one-point exclusion;
  BJ concerns open restrictions of its specified fixed-measure global fields.

## Reproduction and exact pins

Pinned Python3.12.10/SymPy1.14.0, `-X utf8 -I`. Each new helper was executed
once; the earlier unchanged cubic run was reused. Stdout was captured in
separate private `Math/<stem>_review72.out` files, preserving author outputs.
Every non-timing line matched exactly; no stderr or failed assertion occurred.

| Helper | Fresh seconds | Author seconds | Matching non-timing lines |
|---|---:|---:|---:|
| BI_h3_stress | 4.059 | 4.076 | 15 |
| BI_mixed_h6 | 16.593 | 17.284 | 36 |
| BJ_sphere_boundary | 0.101 | 0.103 | 7 |

The cubic checks include all625 reshuffle entries,28 diagonal controls and
two secondary off-diagonal targets. The mixed helper checks all14 full
amplitude polynomials, their trace and three secondary targets at `(2,-3)`.
The standalone BJ helper checks only scalar zonal/eigenvalue/hemisphere
normalizations; its output is not a computational proof of BJ1–BJ3.

SHA-256 after CRLF-to-LF normalization; paths are relative to
`computation/geometric-unity/`. Proof/metadata review-status edits are not
pinned by this table.

| Source, author output or shared dependency | SHA-256 |
|---|---|
| BI_h3_stress.py | a05a625ee582a619bffd15132978a9aa6bf95e3705d43f18944a6aa1ab3fbd2e |
| BI_h3_stress.out | 2b0d87c4829ef922ea6edde364864f677066f7664af69a32c0670c7a525c69b1 |
| BI_mixed_h6.py | 16caa50a3d4e7a0faa8ba08654d471d335a3b30a22cc9d41e9a4c19b11325b55 |
| BI_mixed_h6.out | 18eaf5049db9861b0ffdf914509530e0bdfa26371d7b0dcedf9c4d8958539693 |
| BJ_sphere_boundary.py | b37f277407440d509de892354240bf1dcb9b18d8addf50682dcff39b052bff71 |
| BJ_sphere_boundary.out | c60b9f8ba03c3498e5654e0a72ce2c0f8adcf32536c8a7d742c215d8a41748f1 |
| BC_affine_hooks.py | 9c8871a93e3da243057248a3757f51a39e1dd4710f51108694e50eceb017caba |
| BF_h6_casimir.py | 72f509389a2292e55a7bb0b7769991ea7b6abbccf2bcf01818b747486a5514c5 |
| BG_h6_curvature.py | 6d8c25a290757afd4576896bb44155ef152bc9957ed7cdd26d35415828380736 |
| BC_quadratic_xcheck.py | 8b2849350623974b2a5a3900e3a7491372dd38ba4c93939c06f9b297cf4c5f11 |
| BC_h6_closed_generator.json | 2cac3aa23dac4cdedfcdcf030c8b866f4c6dc8bed583a4006377b3734b2f0eb1 |

## Limits of acceptance

Metadata and proofs agree on real positivity, fixed lambda1/5, all positive
constant scales and the restricted classes. Neither excludes the entire
degree-five or mixed135 module, all curvature-kernel values, arbitrary
closed fields or other extensions with the same value. BJ excludes neither
point-dependent measures nor positive-order distributions under its current
hypotheses, and makes no quotient-descent claim. Curved bases, other
parameters, coincident kernels and other real forms/operators remain outside
the theorem. Classical invariant theory and compact averaging are credited.
No old main, full suite or independent integrated validation was run here;
no physical or literature-priority conclusion is accepted.
