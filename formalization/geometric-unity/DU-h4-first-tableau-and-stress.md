# DU. The H4 first tableau and sharp derivative stress

Let \(K\) be the 10,388-dimensional [BX](BX-h4-curvature-and-closed-fields.md)
H4 hook and dual coefficient-curvature kernel. The H4 derivative
equation is covariant coclosure. The symbol of contraction on
\(V^*\otimes K\) has rank 4,810 and kernel \(g_1\) of dimension
140,622. Exterior multiplication on \(g_1\) has rank 24,024.
The joint closed and coclosed symbol has dimension 116,598.

## Exact first-derivative ranks

The 145,432-dimensional source has 577 highest copies across 36
output types. The verifier reconstructs the multiplicities by the
B2 Racah–Speiser rule from the full BX character. For every supplied
highest tensor it checks the weight, both raising equations, both
hooks and the Hodge-dual coefficient-curvature equation on each
derivative slice. Full column rank modulo 65,537 proves
independence over the rationals. The exact character count then
proves completeness.

For each highest block, contraction is computed using the inverse
metric. Exact rational ranks and independently checked rational
kernel relations identify its full kernel. Applying exterior
multiplication to that kernel and multiplying the resulting ranks
by the irreducible dimensions gives the stated totals.

## Sharp derivative spectrum

For \(J=\nabla C\), let \(\tau\) exchange the derivative and coefficient
slots. The [CH](CH-h4-second-order-stress.md) derivative quadratic is
\[
 F(J)=2\langle J,\tau J\rangle-\|J\|^2+\|dJ\|^2 .
\]
On \(g_1\), the sharp bounds are
\[
 -3\|J\|^2\le F(J)\le\tfrac{17}{5}\|J\|^2 .
\]
The two equality spaces have dimensions 42,537 and 19,019.
The verifier uses the compact Hermitian metric and exterior Gram
minors to rebuild the restricted norm and stress matrices in all
36 blocks. Exact positive semidefinite congruences give the
endpoint bounds and the kernel dimensions. The identity
\(F+3\|J\|^2=\|J+\tau J\|^2+\|dJ\|^2\) also characterizes
lower equality as skewness in the two slots and zero exterior symbol.

At a point preserving scalarity and norm through second order,
\(H(C)+F(J)=0\) by CH. In particular \(\|J\|^2\ge H(C)/3\)
when \(H(C)\ge0\), and \(\|J\|^2\ge-5H(C)/17\) when
\(H(C)\le0\). [CK](CK-h4-invariant-grade-five.md) proves that the
curvature source has both signs even on the invariant value equation.

## Low derivative rank

For a real rank-one derivative \(J=\xi\otimes C\) with unit
\(\xi\), contraction forces \(\iota_\xi C_a=0\). In the alternating
hook, the component along \(\xi\) then forces \(C_\xi=0\).
Hence the exchange contribution is zero and \(\|dJ\|^2=\|J\|^2\).
Thus \(F=0\) at derivative rank zero or one. At the certified CB
partial box the positive CH source rules out both ranks under the
second-order scalarity conditions.

At real derivative rank at most two, the sharp bounds are
\[
 -\tfrac75\|J\|^2\le F(J)\le\tfrac75\|J\|^2 .
\]
Before imposing curvature, rotate a Euclidean orthonormal frame so
that both derivatives lie in a coordinate plane. The differentiated
hooks, contraction and stress form are orthogonal invariants. Six
nonempty sign-block types suffice by permutations within the plane
and its complement. The verifier computes each exact restricted
characteristic polynomial and checks both endpoints.

Two explicit rational ten-entry derivatives in the physical
BC plane \((12,13)\) attain both endpoints after imposing dual
coefficient-curvature. The verifier checks both hooks,
contraction, all curvature coefficients and the exact stress ratios.
Their two derivative slices have disjoint support, giving rank two.
At a rank-two point preserving scalarity and constant norm through
second order, the bound gives \(\|J\|^2\ge(5/7)|H(C)|\).
The CB-box lower bound is therefore
\(\|J\|^2>(5/7)\,1.46150093\).

    python -X utf8 computation/geometric-unity/DU_h4_first_tableau.py --full
    python -X utf8 computation/geometric-unity/DU_h4_stress_spectrum.py
    python -X utf8 computation/geometric-unity/DU_h4_rank_two.py

The spectra concern the linear first tableau. A full field also
requires tangency to the nonlinear point equations and higher-order
compatibility.
