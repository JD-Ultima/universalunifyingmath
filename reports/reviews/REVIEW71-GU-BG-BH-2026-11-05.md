# Review71: H6 fixed values and the degree-two module

2026-11-05, America/New_York. **Accepted: BG1–BG3 and BH1–BH3 within their
stated positive real setting.** I read the final proofs, metadata, sources,
certificate and outputs, and ran only the three new helpers once each. No
mathematical correction is requested. This is internal AI-assisted review;
external peer review remains outstanding.

I previously reviewed BC/BF, audited the cubic phase/metric conventions and
authored the targeted Clifford implementation reused as BC's cross-check.
Those contributions are disclosed. I did not author these new curvature,
moment or diagonal-spectrum proofs. The separate nested implementation is
independent of the primary coefficient algorithm, not of this reviewer or
of the shared cubic formula. This reproduction is not another full geometry
or original-equation engine.

## Complete fixed-value curvature map

The domain is the full23-dimensional SO4-fixed hook space, not BC's smaller
closed-field kernel. The SO4 invariants are reconstructed completely and
both hook constraints are retained. Hodge sends the30 raw degree6 invariants
to a complete degree8 invariant target. Its nonzero coordinate minor makes
the selected30 outputs a coordinate isomorphism, not sampled receivers.

The curvature is full ambient `-[[Phi_i,Phi_j],Phi_a]/4`. With Q=5Phi,
projection against Qd indeed gives denominator2500*g_d. Metric isometry,
all six stabilizer calibrations, skewness, Bianchi and Ric=-5g/4 fix the
model and sign. The coefficient-covector action is -R transpose; it is not
the sum of curvature actions on all tensor slots. Both output routines have
the correct wedge shuffle sign. Horizontal, mixed and vertical inputs all
occur.

The20x20 nonzero minor proves rank at least20. Three independent reconstructed
nullvectors pass the entire sparse curvature output, proving rank at most20
on the complete23-dimensional source. Thus completeness also has a direct
lower-rank/full-nullvector argument; it is not inferred from selected output
vanishing alone. Source-coordinate reconstruction and the entire JSON are
compared exactly. The tensors have738,1358,878 nonzero coefficients.

The positive SO5 Casimir is applied only after verifying SO4 fixation, as
required by BF's shortened formula. Full tensor reconstruction checks every
matrix column and eigenvector, giving distinct eigenvalues10,18,40. BF's
injective matrix-coefficient/sphere argument applies to each nonzero fixed
eigenvector; it does not depend on a closed-field lift of the new one.
The orbit dimensions are14,30,91, with disjoint Casimir eigenspaces. Their
sum135 is the exact generated orbit span and only a lower bound for the
full curvature kernel. Rank at most36473 is consequently an upper bound,
not a computed full ambient rank.

## Fixed three-space and the single-point consequence

The moment calculation uses the actual compact field iB6 Gamma. Right
coordinate volume has square-1/detGram; multiplying its quadratic by the
positive determinant gives normalized Gamma with no extra sign. Complete
polarization retains all six amplitude coefficients. The orthogonal norm,
trace `-(52/3)N`, old BC moments and12 nested-algorithm targets agree.

I checked the completed square for qI-qH directly: its mixed coefficient is
960 and the remaining B squared coefficient is26540/53. All square weights,
including28575/32 for the new direction, are strictly positive. Thus real
scalarity forces all three amplitudes to zero. The assertion concerns the
complete fixed curvature three-space; it is not a conclusion about arbitrary
sums throughout the135-dimensional orbit span.

For an arbitrary smooth full field on the pure branch, BA supplies closure
and coefficient-curvature compatibility, with derivative grade7 only. Varying
amplitudes cannot introduce grade1. At any point with an SO4-fixed value,
the complete curvature restriction therefore places it in this three-space.
Scalar grade1 forces its value to zero, contradicting the positive baseline
`ell=15/(2cV)+(15/52)kappa^2`. SO5 equivariance gives the same result for any
conjugate stabilizer. Neither global affine invariance, constancy nor domain
access to a distinguished origin is assumed for this point-value exclusion.

## Entire degree-two module

BF's spherical proof identifies the Casimir-ten cyclic module equivariantly
with real Sym0(R5). The unique fixed line permits normalization by S4.
The identity `rho(J_i4)^2 S4=2(Si-S4)` gives the actual quarter rotations
`vi=seed+2Ai^2 seed`. The checked cubic action identity justifies the finite
rotation formula. Crucially, the helper restores all six stabilizer terms
when checking Casimir on rotated values; BF's four-term shortcut alone would
not apply there.

The five zonals sum to zero and the first four have positive Gram determinant.
Their eigenvalue coordinates `si=zi-mean(z), s4=-mean(z)` cover every tracefree
diagonal matrix, with inverse zi=si-s4. All four pure and six mixed terms of
each of the ten stress quadratics are computed exactly. Their common formula
has `a=-2375/32,b=1125/16,c=-3625/2`; all coefficient identities are checked.
The ten nested-Clifford targets supply separate finite convention checks.
These are necessary scalarity tests, not an assertion that ten moments
determine the whole stress tensor.

The spectrum lemma covers every multiplicity case. From one unequal pair,
b!=0 forces the other three entries equal. Reapplying the difference identity
shows at most one exceptional entry. Trace zero then gives(-4r,r,r,r,r), and
pair equality gives `-5(b-3a)r^2=0`. Here b-3a=9375/32!=0, so r=0. No division
by an unknown amplitude or repeated-eigenvalue difference loses a chart.
Real symmetric diagonalization by SO5 is legitimate: an eigenvector sign
repairs orientation without changing the diagonal. Therefore the lemma
covers every value of the14-dimensional module, not just its fixed line
or four diagonal amplitudes. A closed-field lift is unnecessary for this
negative pointwise test.

Affine horizontal dilation and positive homothety preserve the kernel,
module identification and scalarity. Fixed seven-covariant tensors give a
common cV^-7 factor on metric-relative quadratic moments; the actual
baseline stays positive. Both exclusions hold at any one point for all
positive scales and both signs of nonzero kappa, at this fixed lambda1/5.
They do not classify other parameters, curved bases, coincident sectors,
the full ambient kernel or mixtures with the other harmonic summands.
No global field-Casimir polynomial or nonlinear solution is inferred.

## Fresh reproduction

Pinned interpreter: `C:\Users\User-2\Documents\ChatGPT\Math\.review-env\Scripts\python.exe`,
with `-X utf8 -I`; each new helper ran once. No imported main or old suite ran.
All author outputs were preserved; separate stdout files are
`Math/BG_h6_curvature_review71.out`, `Math/BG_h6_moments_review71.out` and
`Math/BH_h6_degree2_review71.out`.

| New helper | Fresh time | Author time | Identical non-timing lines |
|---|---:|---:|---:|
| BG_h6_curvature.py | 7.617 s | 7.742 s | 32 |
| BG_h6_moments.py | 2.846 s | 2.806 s | 34 |
| BH_h6_degree2.py | 6.451 s | 6.423 s | 27 |

The curvature helper reconstructs the certificate; the other two read it.
Exact arithmetic checks support the stated finite identities. Completeness,
spherical representation theory, diagonalization, spectrum cases and field
scope are separately audited hand inputs. No literature-priority or physical
claim follows from these results or their internal reproduction.

SHA-256 after CRLF-to-LF normalization:

| Computational file | SHA-256 |
|---|---|
| BG_h6_curvature.py | `6d8c25a290757afd4576896bb44155ef152bc9957ed7cdd26d35415828380736` |
| BG_h6_curvature.out | `945e737feb3e6c10e2291a5bdba1a6d74cd5b094b3828e4cd08a1705a5bd1533` |
| BG_h6_moments.py | `6fa3e73dc5bc4eadf33de10467360ff6ee0133ba5597ca7a5fae92dacc858399` |
| BG_h6_moments.out | `f208ef43ad6c1e5549ddd9fedf3f40311eca76443daf0b6de41f9148a99fa2a5` |
| BH_h6_degree2.py | `8fb3dc6da9d37377145016e291fb8da073a9bc0e8119131ea19f7ee2abaae098` |
| BH_h6_degree2.out | `c2c1dc889e2fc31d91d5314bfdcc2ace03ef6b351969562254d9e1a310138c55` |
| BG_h6_curvature_kernel.json | `eac97d439360d4e63468e16ad356f2e92d1e198e46a9f309d7a99157f44482e6` |
| BC_affine_hooks.py dependency | `9c8871a93e3da243057248a3757f51a39e1dd4710f51108694e50eceb017caba` |
| BF_h6_casimir.py dependency | `72f509389a2292e55a7bb0b7769991ea7b6abbccf2bcf01818b747486a5514c5` |
| BC_quadratic_xcheck.py dependency | `8b2849350623974b2a5a3900e3a7491372dd38ba4c93939c06f9b297cf4c5f11` |
| BC_h6_closed_generator.json input | `2cac3aa23dac4cdedfcdcf030c8b866f4c6dc8bed583a4006377b3734b2f0eb1` |
