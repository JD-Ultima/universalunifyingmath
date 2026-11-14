# BQ. Invariant quadratic constraints on the full H6 curvature kernel

Two exact scalar projections of grades five and nine give necessary
conditions on the full 6,834-dimensional curvature kernel of
[BO](BO-h6-full-curvature-kernel.md). These conditions exclude each
of the 15 individual isotypic components and at least 94 of the
105 possible two-type supports. One excluded sum has dimension
5,290. Every potential solution must carry at least $532/1285$
of its squared norm outside the contained module of
[BN](BN-h6-invariant-grade-nine.md).

The [verifier](../../computation/geometric-unity/BQ_h6_invariant_forms.py)
reconstructs both forms on each multiplicity space with `--full`, checks the [exact
matrices and pencils](../../computation/geometric-unity/BQ_h6_invariant_forms.json),
and reproduces the energy bound. [Review81](../../reports/reviews/REVIEW81-GU-BQ-2026-11-14.md)
records the independent checks.

## BQ1. The invariant forms

Use the observer, real branch and physical quadratic of BN. As
in BN, set $Q_a=5\Phi_a$ and $F=\sum_a Q_a e^a$. Form
$\Omega_r(e_a,\cdot)=\operatorname{tr}(Q_a F^r)/c_r$ for $r=5,9$,
with increasing-index exterior coefficients and positive greatest
common divisors $c_5=50000$ and $c_9=10^{10}$. Exact reconstruction
gives 1,143 nonzero coefficients in each receiver. The trace of
ordered matrix products gives invariance under $SO(5)$.
Set $L_r(B)=\langle\Omega_r,Q_r(iB\Gamma)/i\rangle$ with the
positive metric on every slot. Both scalars must vanish on a full
solution, because curvature and baseline have grade one, and the
variable derivative has grade seven.

The real kernel is a direct sum of isotypic components. The irreducible
types in BO have real forms and irreducible complexifications. In the
split model, the specified highest vectors have rational coefficients.
Compact conjugation is implemented on those vectors as the
$SO(5)$ action of $s=\operatorname{diag}(1,-1,1,-1,1)$ in the
original five-axis basis. Indeed, this matrix interchanges
$e_0+i e_1$ with $e_0-i e_1$ and similarly for axes two and three.
Its action on the split coframe is conjugation on the transported
rational coefficients. Accordingly, taking the real part of each
highest vector applies the same $(1+s)/2$ on the common
representation factor of all copies. The highest and opposite
weights differ for these nonzero types, so this factor is nonzero.
Primitive rescaling of each real seed chooses a multiplicity basis
with a fixed common real direction in the representation factor.

Distinct types are orthogonal for the norm and all invariant bilinear
forms. On a given type with multiplicity $k$, each invariant symmetric
bilinear form is a symmetric $k$ by $k$ matrix on its multiplicity
space times the same invariant inner product on the irreducible
factor. This follows by applying complex Schur's lemma to each pair
of copies and using the corresponding real forms. Let $G_w,A_w,D_w$
be the norm, $L_5$ and $384L_9$ matrices obtained by evaluating on a
common unit direction as above. Diagonal entries come from individual
seeds, and off-diagonal entries come from polarization on their sums.
All formulas use exact inverse-Gram weights, the physical unit volume,
and real rational tensor arithmetic.

## BQ2. Support exclusions

For a real value in one type, write its representation components
in an orthonormal basis of the irreducible factor. The coefficient
vectors in the $k$ multiplicity slots have a positive semidefinite
covariance $M_w$. Thus the total squared norm and the two moments
are sums of $\operatorname{tr}(G_w M_w)$,
$\operatorname{tr}(A_w M_w)$ and $\operatorname{tr}(D_w M_w)$.
The matrices have size one or two. In particular, a combination
$\alpha A_w+\beta D_w$ that is strictly positive on every included
type excludes their entire sum from $L_5=L_9=0$ away from zero.

The certificate gives nine such pencils. Exact leading principal
minors verify positivity on every displayed sum. These sums cover
each of the 15 individual types and 94 of the 105 pairs. The pencil
$2L_5-384L_9$ is strictly positive on the sum of types
$$
(2,0), (3,2), (5,0), (5,2), (6,1), (6,3), (7,1), (7,2), (7,3),
$$
including both copies when the multiplicity is two. That sum has
dimension 5,290. The positive baseline rules out the zero value
in the field equation. These calculations certify the displayed
exclusions without requiring a search to prove they exhaust the
possible support exclusions.

## BQ3. A sharp bound for the two-moment relaxation

Let $J_w$ be one on weights $(2,0),(3,0),(5,0)$ and zero elsewhere,
and set $t=753/1285$. For every type the following exact matrix is
positive semidefinite:
$$
(t-J_w)G_w-\frac{72}{6425}A_w+\frac{51}{257}D_w.
$$
Pair with $M_w$ and sum over all types. At a zero of both invariant
moments this yields
$$
\|B_{(2,0)}\|^2+\|B_{(3,0)}\|^2+\|B_{(5,0)}\|^2
\le \frac{753}{1285}\|B\|^2.
$$
The verifier also finds positive squared norm shares on types
$(3,0),(5,0),(6,3)$ that attain equality and cancel both moments
exactly. Thus the bound is sharp for these two moment equations.
Full grade-five and grade-nine tensor equations impose further
conditions. This does not certify that a full solution attains
the bound, or that one exists.
