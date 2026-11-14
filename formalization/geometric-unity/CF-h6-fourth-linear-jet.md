# CF. An exact fourth-order H6 closure jet at type (5,5)

The unresolved 286-dimensional type $(5,5)$ in
[BU](BU-h6-second-jet-obstruction.md) contains a nonzero real value
with a verified fourth-order Taylor jet satisfying the linear
covariant closure equations through degree three. The certificate
checks all coefficients over $\mathbb Q$ in the split model and
takes the compact real part. This is finite-order compatibility.
The existence of a genuine closed field with that value remains
open.

## CF1. Radial frame and the closure identities

Use the trace metric of [BO](BO-h6-full-curvature-kernel.md) on
$V=\operatorname{Sym}^2_0(\mathbb R^5)$ and
$Y=SL(5,\mathbb R)/SO(5)$. Let $X\in V$ be the radial coordinate
in the section $s(X)=\exp(X/2)$, and put $A=\operatorname{ad}_X$.
This section transports the frame parallelly along radial
geodesics. In the normalization of [BR](BR-h6-closed-field-lifts.md),
the tangent projection is twice the symmetric projection of the
left Maurer–Cartan form. Applying its even and odd series on the
two symmetric-space summands gives the inverse soldering frame
and the connection on frame directions
$$
 E=\frac{A/2}{\sinh(A/2)}=I-\frac{A^2}{24}+O(X^4),\qquad
 \omega=-\frac12\tanh(A/4)=-\frac{A}{8}+\frac{A^3}{384}+O(X^5).
$$
The connection series follows by multiplying the skew part of
$s^{-1}ds$ by $E$. Its action $\rho$ is on all seven covariant
slots. The curvature on the tangent space is
$R_{X,i}=-\frac14\operatorname{ad}_{[X,P_i]}$, where $P_i$
is a frame vector at the identity.

Write the proposed radial frame coefficients as
$$
 B(X)=b_0+b_2(X)+b_4(X),\qquad
 b_d(X)=\sum_{i_1\le\cdots\le i_d}X^{i_1}\cdots X^{i_d}b_{i_1\cdots i_d}.
$$
These are ordinary monomial coefficients. Define
$\delta=\sum_i e^i\wedge\partial_i$ and $T_i=R_{X,i}X$.
Use $A^2 P_i=4T_i$ and the full tensor curvature on the
coefficient value itself. The degree-one and degree-three
closure equations from $\sum_i e^i\wedge(\partial_{EP_i}+\rho(\omega_i))B$
are
$$
 \delta b_2=-\frac12\sum_i e^i\wedge R^{\rm tot}_{X,i}b_0,
$$
$$
 \delta b_4=\sum_i e^i\wedge\left(
 \frac16\partial_{T_i}b_2-\frac12R^{\rm tot}_{X,i}b_2
 +\frac1{24}R^{\rm tot}_{X,T_i}b_0\right).
$$
The inverse frame is even in $X$, and the connection is odd.
An even coefficient jet gives zero closure in degrees zero and
two. The two displayed identities exhaust closure through degree
three for this jet.

## CF2. Exact rational certificate

The [verifier](../../computation/geometric-unity/CF_h6_fourth_jet.py)
uses the exact split model of BO. The adjacent data give a
nonzero type $(5,5)$ highest value and 2,442 quadratic and
116,510 quartic nonzero coefficients. The verifier checks every
coefficient for both hook traces and the full coefficient
curvature wedge. It checks the $(5,5)$ total weight and both
raising equations on the derivative and value slots.

The unique type $(2,0)$ in the curvature kernel is the tangent
module. Weyl-denominator extraction gives zero multiplicity of
$(5,5)$ in $\operatorname{Sym}^d(V)\otimes V$ for $d=0,2,4$.
The checked highest polynomial jets therefore have no coefficient
in this excluded type. They lie in the 6,820-dimensional
subbundle allowed by BU.

Direct polynomial assembly checks every degree-one and
degree-three closure coefficient. The latter equality includes
23,744 nonzero coordinates and all other coordinates equal to
zero. Compact complex conjugation is induced by swapping the
two conjugate pairs of isotropic fundamental axes. The verifier
forms the real part of $b_0$ and obtains its exact squared norm
$3/131072$ in the positive model. Taking compact real parts of
the full jet preserves all the real linear identities. By
$SO(5)$ transport and linear combination, the real evaluation
span of such finite jets contains the entire irreducible type.

## Reproduction and scope

    python -X utf8 computation/geometric-unity/CF_h6_fourth_jet.py

[Review93](../../reports/reviews/REVIEW93-GU-CF-2026-11-14.md)
records the independent commutator and normalization checks.
The exact verifier uses rational identities throughout. Prime
arithmetic was used to propose the certificate. The calculation
gives finite Taylor data for linear closure. Higher-order
compatibility in characteristic zero, convergent field lifts and
the nonlinear condition remain open. External peer review
remains outstanding.
