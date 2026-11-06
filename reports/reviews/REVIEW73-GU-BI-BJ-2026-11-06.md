# Internal review of the cubic, mixed-value and orbit-average obstructions

2026-11-06. **BI1–BI3 and BJ1–BJ4 are accepted under their stated hypotheses.** I read both final proofs and metadata, all three helpers and outputs, and the classical invariant-theory reference. Three fresh helper runs passed. Every non-timing line matches the author output:

| Helper | Fresh runtime | Matching non-timing lines |
|---|---:|---:|
| `BI_h3_stress.py` | 4.328 seconds | 15 |
| `BI_mixed_h6.py` | 17.697 seconds | 36 |
| `BJ_sphere_boundary.py` | 0.105 seconds | 7 |

No mathematical correction or theorem narrowing is needed. This is internal AI-assisted review, not external peer review. The reviewer authored the mixed helper and its private precursor, contributed earlier BC tensor and spherical-module work, and independently hand-audited the finite-measure argument before publication. The contributions and shared libraries make this a collaborative assessment rather than blind independent discovery. The cubic and mixed helpers use selected checks from a separately implemented Clifford algorithm, but both algorithms implement the same mathematical cubic formula.

## BI: the whole cubic module

The universality step is valid. After trace-free projection and polarization there are ten standard-vector slots. Since dimension five is odd and the tensor power is even, adjoining $-I$ extends $SO(5)$ invariance to $O(5)$. Complexification of the real invariant coefficient fits the pairing-spanning theorem in [Goodman's Theorem 7.4 and Corollary 7.5, printed page 26](https://sites.math.rutgers.edu/~goodman/pub/nuslect.pdf#page=26), which I checked directly. No independence or stable-dimension assertion about all pairings is required.

Self-pairings vanish by trace-freeness. Counting edges between the two cubic and two quadratic factors gives $k=j+1$. The cases $(3,2),(2,1),(1,0)$ leave respectively one, one and two contraction patterns. Thus the four displayed invariants span the actual stress quadratic. The nonzero four-control determinant determines their coefficients universally; it is not a sample fit replacing completeness. The zonal and $x_0x_1x_2$ tensor reconstructions have the correct generator normalization and compact-real volume phase.

I checked the operator decomposition directly. The identity projection has cross trace $(14/5)I$; the displayed embedding of a trace-free symmetric $H$ has cross trace $(27/20)H$. Subtracting these leaves every trace zero, and full symmetrization splits it into the symmetric-four and zero-symmetrization components. Projected reshuffling acts by $3/10,1/10,1,-1/2$ on the four components, respectively. These follow from the formulas for arbitrary tensors; the 625-entry controls corroborate their implementation.

The three nonscalar multipliers are exactly the stated nonzero values. Scalar symmetric stress therefore forces $K_T=L_T^*L_T$ scalar. This positive Gram operator has rank at most five on a 14-dimensional space, so it must be zero. Its trace is the positive real tensor norm, forcing $T=0$. Only the symmetric stress was needed; no untested antisymmetric part is assumed absent. Realness and the positive metric are essential.

## BI: the separate mixed family

Both harmonic power sums are reconstructed by calibrated quarter rotations, with complete tensor annihilators and polynomial-action checks. The polynomial Lie-action ranks establish zero infinitesimal stabilizers; distinct Casimir components prevent cancellation of that conclusion in a nonzero mixed value. The family consequently goes beyond conjugate-$SO(4)$-fixed values.

All fourteen diagonal quadratic polynomials retain the cross terms. The exact norm, trace and alternative Clifford targets agree. The strictly positive-square difference excludes every nonzero real amplitude pair, including both projective charts. Necessary diagonal tests suffice; no full off-diagonal stress or irreducible ten/four splitting is inferred. This family does not permit arbitrary or independently rotated components.

## BJ: finite measures and open restrictions

The construction is well-defined on $K/L$ because $L$ fixes the entire global seed. Compact-parameter bounds on all local derivatives justify integration and closure for any fixed real measure of finite total variation, including singular and atomic measures. Different measures need not define different fields.

The positive hyperplane parabolic and $S=KP$ give the exact inverse-transpose sphere profile. No character or Jacobian is missing: the whole seed is $P$-fixed and the integration variable is unchanged. The spherical degrees three/five belong to the value map; they are not a false global Casimir polynomial for the original field.

Tonelli applied to total variation supplies an uncharged equator. The collapse then uses a uniformly bounded integrand in one tensor fiber after isometric transport. Dominated convergence gives a real multiple of a rotated seed. Comparing the rotated seed directions forces that coefficient to zero, and the positive baseline then contradicts the limiting equation. Zero collapse coefficients and measures representing the zero field are included. No hemisphere-transform injectivity or derivative-jet limit is assumed.

The analyticity step also holds for singular measures. Positive Cholesky coordinates give an analytic section; a complex-neighborhood shrink keeps every real-sphere denominator in one uniformly bounded right-half-plane branch. Finite-total-variation integration therefore preserves the local holomorphic extensions. The analytic identity theorem applies to the residual of an already globally defined fixed-measure profile. It does not continue an arbitrary smooth field, assume its domain meets the collapse curve, or add the degenerate endpoint to the metric manifold.

BJ4's compact orthogonality proof gives a genuine injective 121-dimensional closed-field lift. The final non-affine paragraph is valid: $P$-invariance would force $w=\alpha(w)B_*$, while $\alpha(B_*)=33/128$ then forces $w=0$. The scalar helper verifies its optional normalization and hemisphere constants only; it does not verify BJ1–BJ3's global theorem. Positive-scale transfer and the common $c_V^{-7}$ moment weight preserve all stated contradictions and the positive baseline for either nonzero $\kappa$ sign.

## Reproduction, receipts and limits

Each invocation was `python -X utf8 -I computation/geometric-unity/<helper>.py`, using pinned Python 3.12/SymPy 1.14. Separate transcripts preserve the fresh runs; author outputs were untouched. Only these three new helpers ran once each in this review. No older main, prototype, full suite or source-pin change occurred.

Canonical LF SHA-256 receipts, without BOM; paths are in `computation/geometric-unity/`:

| Artifact | SHA-256 |
|---|---|
| BI cubic source | a05a625ee582a619bffd15132978a9aa6bf95e3705d43f18944a6aa1ab3fbd2e |
| BI cubic output | 2b0d87c4829ef922ea6edde364864f677066f7664af69a32c0670c7a525c69b1 |
| BI mixed source | 16caa50a3d4e7a0faa8ba08654d471d335a3b30a22cc9d41e9a4c19b11325b55 |
| BI mixed output | 18eaf5049db9861b0ffdf914509530e0bdfa26371d7b0dcedf9c4d8958539693 |
| BJ scalar source | b37f277407440d509de892354240bf1dcb9b18d8addf50682dcff39b052bff71 |
| BJ scalar output | c60b9f8ba03c3498e5654e0a72ce2c0f8adcf32536c8a7d742c215d8a41748f1 |
| BC library | 9c8871a93e3da243057248a3757f51a39e1dd4710f51108694e50eceb017caba |
| BF library | 72f509389a2292e55a7bb0b7769991ea7b6abbccf2bcf01818b747486a5514c5 |
| BG curvature library | 6d8c25a290757afd4576896bb44155ef152bc9957ed7cdd26d35415828380736 |
| Separate Clifford library | 8b2849350623974b2a5a3900e3a7491372dd38ba4c93939c06f9b297cf4c5f11 |
| BC seed certificate | 2cac3aa23dac4cdedfcdcf030c8b866f4c6dc8bed583a4006377b3734b2f0eb1 |

The pointwise results retain the positive flat $\lambda=1/5$ model and declared pure compact-real $W_A$ branch. BJ additionally requires a fixed global finite-real-measure representation. Arbitrary quintics, general mixed curvature values, other closed extensions, point-dependent measures, positive-order distributions, curved bases, other parameters, coincident sectors and quotient descent are not classified. No physical conclusion, unrestricted field exclusion or literature priority follows. External peer review remains outstanding.
