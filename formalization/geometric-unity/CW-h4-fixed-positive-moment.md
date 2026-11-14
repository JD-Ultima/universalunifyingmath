# CW. Positive all-grade moments on the H4 permutation-fixed space

The full permutation-fixed subspace of the H4 curvature kernel
has dimension 68. It admits a rational positive definite covariance
annihilating every coefficient of the homogeneous grade-one,
grade-five and grade-nine residual. Consequently every positive
semidefinite linear combination of those quadratics vanishes
identically on this space.

## CW1. Exact input and output spaces

Use the ordinary permutation action by conjugation on
\(V=\operatorname{Sym}^2_0\mathbb R^5\), as in
[CN](CN-h4-permutation-cubic-exclusions.md). The rational data
give 68 independent physical-frame tensors. The verifier checks
both hook traces, the complete coefficient curvature on their
Hodge duals and all four adjacent \(S_5\) generators. The full
[BX](BX-h4-curvature-and-closed-fields.md) character and exact
Weyl restriction give dimension 68 for the fixed space.

Use the homogeneous residual
\(q(B)_a=Q(B)_a+(6/7)\|B\|^2\gamma_a\). The exterior identity
of [CD](CD-hook-clifford-identities.md) and the completely
alternating top output give these bounds for the fixed quadratic
output space:

| Output | Dimension bound |
| --- | ---: |
| Trace-free symmetric grade one | \(5-1=4\) |
| Grade five with zero exterior trace | \(223-28=195\) |
| Grade nine from the eight-form trace | 28 |
| Total | 227 |

The dimensions follow by exact permutation character averaging on
\(\operatorname{Sym}^2V\), \(V\otimes\Lambda^5V\),
\(\Lambda^6V\), and \(\Lambda^8V\). The exterior map splits
equivariantly, including on the fixed space. Choose four actual
grade-one coefficients, 195 grade-five coefficients and 28
grade-nine coefficients. The square minor certified below proves
that these 227 coordinates detect the complete homogeneous output.

## CW2. Exact covariance

In the exact basis \(B\), the dyadic matrix \(L\) gives a
first covariance \(M_0=LL^T\) with \(M_0>I/4\) by
fraction-free Sylvester elimination. Write \(F_j=(BL)_j\)
for its 68 rank-one factors. The selected 227 covariance
coordinates define a square matrix \(A\), whose diagonal
columns are \(q(B_i)\) and off-diagonal columns are
\(2q(B_i,B_j)\) by polarization.

The verifier rounds the unit-frame coordinates of \(B_j\)
and \(F_j\) to points \(D_j,T_j\) with dyadic denominator
\(2^{36}\). The physical coframe metric factors are bounded
by 80-bit rational square-root intervals. Independently for
each column the approximation error is less than
\(\delta=1/(2\cdot10^9)\). The same bounds give
\(\|D_j\|<1.09\), \(\|B_j\|<1.1\), and
\(\sum_j\|T_j\|^2<80\).

Evaluate the selected dyadic Clifford columns \(A_D\) and
the dyadic residual \(v_D=\sum_jq(T_j)\) exactly.
On ambient tensors the polynomial from the five reduced
Clifford terms obeys
\(\|q(U)\|\le200\|U\|^2\) and
\(\|q(U,W)\|\le200\|U\|\|W\|\).
Indeed, normalized Hilbert–Schmidt norms give a product factor
\(\sqrt{128}\). Cauchy–Schwarz over the 14 external indices
and the scalar correction bound the quadratic by
\((10/3)\sqrt{1792}+(6/7)\sqrt{14}<200\).
Polarization and rescaling give the bilinear bound.
Thus an off-diagonal entry of \(A-A_D\) has absolute
value at most \(880\delta\), and
\(\|A-A_D\|_F\le227\cdot880\delta\).
The residual error is at most
\(200\delta(148+68\delta)\).

A dyadic approximate right inverse \(Y\) satisfies
\(\|Y\|_F<450\). Exact integer row and column norms for
\(I-A_DY\), followed by the preceding perturbation bounds,
give the rational inequalities
$$
 \|I-AY\|<1/20,\qquad
 \frac{2\|Y\|\|\sum_jq(F_j)\|}{1-\|I-AY\|}<1/64.
$$
Neumann inversion proves \(A\) invertible. Put the vector
\(-A^{-1}\sum_j q(F_j)\) into the selected covariance
coordinates and zero elsewhere. Its symmetric matrix
\(\Delta\) has operator norm less than \(1/64\), so
\(M=M_0+\Delta>(15/64)I\). Its output vanishes in the
selected coordinates and hence completely. The fixed unit
coframe scalings cancel row by row in the linear solve, so
\(\Delta\) and \(M\) are rational in the exact basis.

If a receiver gives a positive semidefinite quadratic
matrix \(H\), then \(\operatorname{tr}(HM)=0\) and \(M>0\)
give \(H=0\). The covariance represents a linear moment
relaxation. Existence of a single nonzero real value satisfying
the full equations remains open.

## Reproduction

    python -X utf8 computation/geometric-unity/CW_h4_s5_moment.py --full

The integer Clifford helper uses the original 128-dimensional
matrix algebra in two parity blocks. For each matrix entry,
two signed 21-bit limbs give products of dot length at most
896, with \(2\cdot896(2^{21}-1)^2<2^{53}\).
Every binary64 integer product and partial sum is exact.
The helper separately checks all 227 receiver coefficients on
a generic rational hook using the original Clifford convention.
The full option recomputes all 363 integer evaluations.
The default verifies the geometry, rational bounds and six
representative evaluations against the archive.

[Review108](../../reports/reviews/REVIEW108-GU-CW-CX-2026-11-14.md)
records the checks. External peer review remains outstanding.
