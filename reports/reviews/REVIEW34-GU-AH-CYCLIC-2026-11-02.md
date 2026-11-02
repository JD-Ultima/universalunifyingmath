# Internal review of AH5: the explicit cyclic five-form orbit

Scope: AH5 in `formalization/geometric-unity/AH-lambda5-gram.md`, its updated metric record and `AH_cyclic_obstruction.py/.out`. This is an independent internal proof and executable-certificate review. Review31 retains its AH1–AH4 scope; neither report is external peer review.

## Assessment and scope

AH5 is valid as stated. At a nonzero-w point on the first single Lambda5 c/wGamma branch, a full solution cannot have its five-form in the nonzero real scaling/SO(14) orbit of AH4's exact cyclic form. No correction or narrowing was required.

The result keeps the real compact (14,0), u(128), fixed parallel WA, regular slaving and single-kernel assumptions of AH. It excludes one specified pointwise orbit, not arbitrary five-forms, the opposite pair, zero-w endpoints or coincident kernels. The original Gram-channel control remains correct: AH5 shows why admission by that channel was insufficient for a full solution.

## Reduction without discarding derivatives

The cyclic form has N=14 and Gram 5I. A nonzero real scaling and orthogonal rotation leave its Gram isotropic. Substitution in AH2 makes a scalar identity equal to `3328*w*C tensor C`. At nonzero w, rank fourteen versus rank at most one forces C=0 and hence v=0 at the point.

The proof correctly does not differentiate this pointwise conclusion. The vector derivative may be nonzero. Its input type is `V tensor V`, whose scalar, alternating-two and symmetric-traceless constituents have no overlap with `V tensor Lambda5`, or the right-Gamma equivalent grade-nine target. Therefore it supplies neither selected grade even for arbitrary derivative jets.

For the first pair, the kernel fields have grades four and eight with compact phase i. The WA triple anticommutator permits grades three/seven and seven/eleven, respectively. The parallel connection and volume element preserve this exclusion for form, scale, rotation and chart-weight derivatives. The invariant terms have only grades one/thirteen and the WA curvature term grade one, while quadratic vector and mixed terms vanish at v=0. Thus the full grade-five/nine equations at the point reduce to the corresponding components of Q(k). No other derivative equation is claimed solved.

## Exact coefficient and phase certificate

The raw integer polarizations have the correct compact conversion:

```text
Q(i*(r*c+s*d))
= -i*(r^2*M(c,c)+2*r*s*M(c,d)+s^2*M(d,d))/6,
d = w_lift(Omega)*Gamma.
```

This follows directly from T's WA coefficient and quadratic normalization. The cross coefficient contains the necessary factor of two. The direct phased control checks this reconstruction in every grade-five and grade-nine slot, not only in the three displayed entries.

At external index zero, the grade-nine coefficient with Clifford indices `(0,1,3,4,6,8,9,10,11)` is

```text
-(16i/3)*(4*r^2+17*s^2).
```

It cannot vanish for a real chart with weights not both zero. Multiplying the form by a nonzero real scale multiplies the coefficient tensor by its square. SO(14) equivariance preserves nonvanishing of that tensor, although the particular nonzero coordinate component may change. This proves the entire declared orbit exclusion, including fields with varying pointwise rotations and scales, since the derivative exclusions allow arbitrary jets.

The additional two entries and the exact ideal `(r^2,rs,s^2)` are consistent checks; the real positive quadratic already suffices. No all-five-form contraction ansatz is used. In particular the sixteen-slot grade-five target cannot inherit AA6's twelve-slot uniqueness argument.

## Fresh reproduction and limits

The ordinary portable command

```text
python -X utf8 -I computation/geometric-unity/AH_cyclic_obstruction.py
```

passed all fourteen reported checks in 5.5 seconds with Python 3.12 and SymPy 1.14.0. The exact output reproduces the three coefficient polynomials, their ideal, compact phase and direct control. At AH4's regular weights the entries are `80/3`, `-320*sqrt(14)` and `-1984/3`. The same direct computation confirms its symmetric Gamma channel still balances, while the remaining grade-nine equation fails. Nonzero kernel-derivative and vector-derivative controls verify the stated grade exclusions.

The certificate shares T's accepted Clifford implementation; it is not a second field-equation engine or new physical evidence. Universality in the two real weights follows from exact homogeneous polarization, and extension along this orbit from scaling and equivariance. The calculation does not classify the remaining five-form space or prove a full Lambda5 exclusion. All of AH's other nonlinear and differential questions outside this one orbit retain their stated boundaries.
