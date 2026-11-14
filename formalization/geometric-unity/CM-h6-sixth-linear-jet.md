# CM. Exact sixth-order H6 closure at type (5,5)

The open 286-dimensional type \((5,5)\) in the H6 linear
extension problem has a nonzero real sixth-order Taylor jet
satisfying covariant closure through degree five. The certificate
extends the [CF](CF-h6-fourth-linear-jet.md) fourth jet, preserves
both hook traces and the full curvature condition, and excludes
the 14-dimensional type eliminated by [BU](BU-h6-second-jet-obstruction.md)
from every coefficient.

## CM1. Degree-five closure

In CF's radial parallel frame \(s(X)=\exp(X/2)\) write the
ordinary polynomial jet as \(B(X)=b_0+b_2(X)+b_4(X)+b_6(X)\).
With \(A=\operatorname{ad}_X\), the inverse soldering and connection
have expansions
$$
 E=\frac{A/2}{\sinh(A/2)}
   =1-\frac{A^2}{24}+\frac{7A^4}{5760}+O(X^6),\qquad
 \omega=-\frac12\tanh(A/4)
   =-\frac{A}{8}+\frac{A^3}{384}-\frac{A^5}{15360}+O(X^7).
$$
The connection acts on all seven covariant slots. Put
\(T_i=R_{X,i}X\) and \(U_i=R_{X,T_i}X\) for the tangent curvature.
The same convention gives \(A^2P_i=4T_i\), \(A^4P_i=16U_i\), and
\(R_{X,i}^{\mathrm{tot}}=-\rho([X,P_i])/4\).
Extracting degree five from the closure operator gives
$$
 \delta b_6=\sum_i e^i\wedge\left(
 \frac16\partial_{T_i}b_4-\frac12R_{X,i}^{\mathrm{tot}} b_4
 -\frac7{360}\partial_{U_i}b_2+\frac1{24}R_{X,T_i}^{\mathrm{tot}}b_2
 -\frac1{240}R_{X,U_i}^{\mathrm{tot}}b_0
 \right),
$$
where \(\delta=\sum_i e^i\wedge\partial_i\) is the coordinate
exterior derivative. The frame coefficients and the inverse
soldering are even in \(X\), and the connection is odd. Closure
therefore has only odd degrees. The CF equations and this last
equation exhaust the constraints through degree five.

## CM2. Exhaustive rational certificate

The sixth coefficient has 1,936,021 nonzero tensor entries in
13,817 monomials. The full verifier checks every coefficient for
both hook traces, the complete coefficient-curvature wedge
equation and projection onto the excluded \((2,0)\) module.
It checks the two highest-weight raising equations on the
polynomial and tensor slots together. The direct projection uses
an independently verified basis of the 14-dimensional \((2,0)\)
module and its nondegenerate invariant Gram matrix.

The assembled left and right sides of the quintic equation have
604,246 nonzero coordinates and coincide as rational polynomials.
A separate backwards-coordinate calculation agrees on 986
quintic coordinates. The compact real part has base squared norm
\(3/131072\) by CF. Taking the real part preserves all the
verified linear identities. Its compact translates span the
whole \((5,5)\) module at the base point. Together with BR's
6,534-dimensional span of genuine field values, this gives
sixth-order Taylor data over the full H2-free 6,820-dimensional
value module.

## Reproduction

    python -X utf8 computation/geometric-unity/CM_h6_sixth_jet.py --full

The [verifier](../../computation/geometric-unity/CM_h6_sixth_jet.py)
reads lossless rational data from three checksum-verified gzip
streams. It evaluates 296,775 scalar factors against 3,963
coefficient tensors and then tests the expanded polynomial.
The binary streams encode base-128 integer lengths and signed
little-endian integer magnitudes. The adjacent
[decoder](../../computation/geometric-unity/CM_h6_sixth_codec.py)
defines the format. The proof only uses the resulting rational
tensors.

The default selected regression checks all of CF, the weights
throughout the factored certificate, and selected sixth-order
coefficients and identities over the rationals. For the exhaustive
check use \(\texttt{--full}\). Optional \(\texttt{python-flint}\)
accelerates the exact rational calculation. The
\(\texttt{--fraction}\) option runs the same check using the
standard-library rational backend.

[Review104](../../reports/reviews/REVIEW104-GU-CM-2026-11-14.md)
records the two exhaustive runs and the independent geometry
controls. The result supplies finite-order data for linear
closure. [CU](CU-h6-global-analytic-closure.md) subsequently
constructs analytic fields by a different method. General
formal integrability and the nonlinear equations remain open.
External peer review remains outstanding.
