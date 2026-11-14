# DC. Four-value H4 moment and a quartic separation

There are four linearly independent real algebraic H4 values
\(F_1,\ldots,F_4\) in the 68-dimensional permutation-fixed
curvature kernel for which \(\sum_r q(F_r)=0\), where
\(q(B)=Q(B)+(6/7)\|B\|^2\gamma\) is the full homogeneous
Clifford residual. Their covariance has rank four. On their
four-dimensional span the same equations exclude every nonzero
point value even over the complex numbers. The squared residual
has a strictly positive quartic lower bound on that real span.

## DC1. An exact four-value covariance

Use the exact real rational basis \(B=(B_i)_{i=1}^{68}\) and
complete 227 actual receiver coordinates of
[CW](CW-h4-fixed-positive-moment.md). For a real
\(68\times4\) matrix \(X\), write \(F_r=BX_r\).
Take the 228 equations consisting of the selected coefficients
of \(\sum_r q(BX_r)\) and \(\|X\|_F^2-1=0\).
The supplied affine slice varies 228 of the 272 entries
near an exact dyadic matrix \(X_0\) and fixes the rest.
Fraction-free elimination gives \(X_0^T X_0>I_4/5\).

The verifier reconstructs the exact basis \(B\) and four
anchor tensors, then uses interval square roots to round
their unit-frame coordinates to dyadics. Every column is
within \(10^{-15}\) of its exact tensor. The dyadic basis
Gram is bounded above by \(3I/2\), so \(\|B\|<\sqrt2\).
Exact Clifford evaluations give the center residual and the
complete square Jacobian of the affine system by ordinary
polarization. The final norm row is rational and exact.
The dyadic polynomial uses the ambient five-term reduction
of the Clifford formula. Its restrictions to exact hooks
agree with the original formula by the hook identities.

As in CW, the homogeneous Clifford map and its polarization
are bounded by \(200\|U\|^2\) and \(200\|U\|\|W\|\).
The basis bound gives a Jacobian Lipschitz bound \(802\)
including the scalar norm row. The interval arithmetic
bounds the difference between the exact and dyadic Jacobians
by \(16\cdot840\cdot10^{-15}\), and the center residual
difference by \(1600\cdot10^{-15}\).
Use the dyadic inverse approximation \(Y\) with
\(\|Y\|_F<1800\). Its exact inverse defect at the center
is below \(10^{-3}\), and the Newton displacement there
is at most half of the radius \(\rho=10^{-8}\).
Throughout the closed radius-\(\rho\) ball the associated
map \(t\mapsto t-Yf(X_0+St)\) has derivative norm below
\(1/50\). Here \(S\) inserts the 228 selected affine entries.
The contraction gives an exact zero and an invertible
derivative on the slice.

The factor at the root has smallest singular value greater
than \(2/5\), so its four tensor columns are linearly
independent. The unit receiver rows have fixed nonzero
scalings from the rational coframe. Removing these scalings
gives a rational polynomial system with an isolated nonsingular
root on the fixed rational slice. Thus the actual factor
entries are real algebraic. The exact covariance
\(\sum_r F_r\otimes F_r\) annihilates the full residual by
the receiver completeness of CW.

## DC2. The equations on the four-dimensional span

Let \(L:\operatorname{Sym}_4\to\mathrm{output}\) be the
moment map in the actual four root tensors, so \(L(I_4)=0\).
On the traceless subspace use the three diagonal differences
\(E_{ii}-E_{44}\) and six off-diagonal symmetric generators.
Six additional exact anchor-pair calculations recover the
corresponding nine columns at the dyadic anchors.
Each true root anchor is within \(2\cdot10^{-8}\) of the
dyadic one. The verifier bounds the corresponding perturbation
of a selected square minor, giving an exact inverse norm below
266 throughout the root ball.

Hence \(\ker L=\mathbb R I_4\), and the complexified kernel
is \(\mathbb C I_4\). If \(q(\sum_i x_iF_i)=0\) over
\(\mathbb C\), then \(xx^T\) is a scalar matrix. Its rank
is at most one, so that scalar and \(x\) both vanish.
Equivalently the equations span the mixed monomials and
the differences of coordinate squares. Their cubic multiples
give every degree-three monomial. A nonnegative quadratic
receiver has trace zero against \(I_4\), so it is identically
zero on the span.

For real \(x\) the traceless matrix
\(M=xx^T-\|x\|^2I_4/4\) has squared Frobenius norm
\(3\|x\|^4/4\). In the stated basis its coordinate vector
\(z\) satisfies \(\|M\|_F^2\le4\|z\|^2\).
The inverse norm gives \(\|q(\sum_i x_i F_i)\|^2\ge
3\|x\|^4/1132096\). Since \(\|B\|<\sqrt2\) and the
factor \(X\) has Frobenius norm one, the physical bound is
$$
 \|q(U)\|^2\ge\frac3{4528384}\|U\|^4
 \qquad\bigl(U\in\operatorname{span}_{\mathbb R}\{F_1,\ldots,F_4\}\bigr).
$$

## Reproduction

    python -X utf8 computation/geometric-unity/DC_h4_rank4_moment.py --full

The full check reconstructs 306 integer Clifford evaluations
and all rational bounds. The default checks the same geometry
and inequalities, six representative derivative or center
evaluations, and all six quartic pair evaluations. The
dyadic calculation uses three signed 21-bit limbs and
the 128-dimensional Clifford algebra in parity blocks.

The covariance concerns four values. The quartic obstruction
concerns their span. Existence or exclusion of a full H4 value
on the complete curvature kernel remains open.
[Review110](../../reports/reviews/REVIEW110-GU-DC-2026-11-14.md).
