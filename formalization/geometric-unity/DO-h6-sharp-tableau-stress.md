# DO. The sharp H6 first-tableau stress bound

Put \(\beta=2+1/\sqrt5\). On the complete first-derivative tableau
\(g_1\) of [DM](DM-h6-first-tableau.md), the derivative term in the
[CJ](CJ-h6-universal-stress-gradient.md) stress-divergence identity
satisfies the sharp algebraic inequality
\[
 2X-Y\ge-\beta\|\nabla B\|^2 .
\]
Equality in this linear tableau holds on precisely one copy each of
types \((8,1)\) and \((8,3)\), of total dimension 2,565.
For every full field,
\[
 \|\nabla B\|^2>\frac{5193}{25000\beta}\|B\|^2,
 \qquad
 H(*B)+\frac{1}{\sqrt5}\|\nabla\Pi_{55}B\|^2
 \le\beta\|\nabla B\|^2 .
\]
The first coefficient is approximately \(0.0848802084\).
A nonzero full-field point also requires real derivative rank at
least two.

## DO1. Spectrum of the derivative quadratic

Write \(J=\nabla B\) and let \(\tau\) exchange its derivative and
coefficient slots. The actual Clifford eight-form is the form-slot
star \(C=*B\) up to sign. The star is parallel, isometric and
commutes with \(\tau\). Thus CJ's term \(X\) equals
\(\langle J,\tau J\rangle\). The exterior norm identity gives
\(Y=\|\nabla C\|^2-\|d^\nabla C\|^2\), and
\(\|d^\nabla C\|=\|\delta^\nabla B\|\). On all derivative tensors,
\[
 F(J):=2X-Y=2\langle J,\tau J\rangle-\|J\|^2+
 \|\delta J\|^2.
\]

The 34 nonzero tableau multiplicity blocks of DM determine its
spectrum. The verifier computes the compact Hermitian Gram,
the slot-exchange form and the contraction Gram directly from
raw witnesses and exterior metric minors. Restrict to the exact
closure kernel in each block, and call the resulting norm and
quadratic matrices \(G,Q\). All \(G\) are positive definite.
The two exceptional blocks have normalized determinants
\[
 \begin{array}{ll}
 (8,1):&\det(Q-tG)/\det G=(t-1)^2(5t^2+20t+19)/5,\\
 (8,3):&\det(Q-tG)/\det G=-(t-1)(5t^2+20t+19)/5.
 \end{array}
\]
Their least eigenvalue is \(-2-1/\sqrt5\), simple in both
multiplicity spaces. On every other block \(Q+(12/5)G\)
is positive by exact rational elimination. Thus the least
eigenvalue is \(-\beta\) and the equality dimensions are
\(836+1729=2565\).

## DO2. Stress compatibility and the mixed gradient penalty

CJ gives \(H(*B)+F(\nabla B)=0\) for a closed field preserving
scalarity and norm through second order. Invariants \(L_5=L_9=0\)
give \(H>(5193/25000)\|B\|^2\). With just scalarity and constant norm,
CJ gives the lower bound \(19\|B\|^2/(100\beta)\) for the
squared gradient. At every value in the certified BS box the
squared gradient must exceed \(1.67097264/\beta\) if those
second-order conditions hold.

By [DK](DK-h6-mixed-gradients.md), the \((5,5)\) fiber
contribution to \(J\) takes place in the three output blocks
\((7,5),(6,5),(6,4)\). On each entire block the matrix
\(Q+2G\) is positive semidefinite by the same exact calculation.
Let \(P\) project onto the sum of those output blocks.
On their orthogonal complement use the sharp bound of DO1.
It follows that
\[
 F(J)+\beta\|J\|^2\ge(\beta-2)\|PJ\|^2
 \ge\tfrac1{\sqrt5}\|\nabla\Pi_{55}B\|^2 .
\]
Together with stress compatibility this gives the stated penalty.

## DO3. Rank-one derivatives have zero stress contribution

If \(J=\xi\otimes b\) is a real rank-one element of \(g_1\),
normalize \(\xi\) to have unit norm. Closure gives
\(b_a=\xi\wedge c_a\) with \(\iota_\xi c_a=0\).
The contraction hook is
\(c_\xi-\xi\wedge\sum_a\iota_a c_a=0\).
Its two terms lie in orthogonal form summands, so \(b_\xi=0\).
Therefore \(X=\|b_\xi\|^2=0\) and
\(\|\delta J\|^2=\|J\|^2\), giving \(F(J)=0\).
The same is true for zero derivative. CJ has \(H(*B)>0\) on
every nonzero H6 admissible value, so neither derivative rank can
satisfy the double-divergence identity at a full-field point.

Real rank two can give a negative derivative term in the linear
system. In [CU](CU-h6-global-analytic-closure.md)'s complex field,
multiply by a single transverse coordinate and take the real
part. At the origin the field is zero and its derivative is
\(\operatorname{Re}(e^j\otimes b_{55})\). In the first five
transverse covectors, the choices \(j=3,9\) give quotient
\(F/\|J\|^2=-1\), and \(j=0,1,8\) give \(-1/3\).
These jets satisfy linear closure. Their compatibility with
nonlinear tangent equations at a nonzero value remains open.

    python -X utf8 computation/geometric-unity/DO_h6_stress_spectrum.py
    python -X utf8 computation/geometric-unity/DO_h6_stress_spectrum.py --verify-tableau

The first command rebuilds all 34 rational Gram pencils.
The second first rechecks the complete DM tableau membership and ranks.
