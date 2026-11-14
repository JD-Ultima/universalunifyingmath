# DH. A strict bundle spectral gap for H6 type \((5,5)\)

On \(SL(5,\mathbb R)/SO(5)\) with the trace metric, the
self-adjoint rough Laplacian on \(L^2\) sections of the pure
\((5,5)\) bundle satisfies
$$
 \nabla^*\nabla\ge29/4.
$$
The linear closure equation on this type implies rough
eigenvalue \(7\) by [DE](DE-h6-elliptic-periodic.md).
Thus the gap to the relevant \(L^2\) bundle spectrum is
at least \(1/4\).

## DH1. Tempered support of the bundle

Use the full group Plancherel decomposition and trace Casimir
of [DF](DF-h6-l2-vanishing.md). Standard real parabolics of
split \(SL(5,\mathbb R)\) have ordered block sizes \(n_j\)
summing to five. Once the real split center is removed,
the connected semisimple Levi factors are \(SL(n_j,\mathbb R)\).
A block with \(n_j\ge3\) has complex rank \(n_j-1\) strictly
greater than the compact rank \(\lfloor n_j/2\rfloor\).
The compact-Cartan criterion thus leaves only sizes one
and two for discrete inducing data. In particular there
can be at most two \(SL(2,\mathbb R)\) factors.
[Harish-Chandra, Theorem 13](https://doi.org/10.1007/BF02392813).

For an \(SL(2,\mathbb R)\) discrete series of lowest absolute
rotation weight \(k\ge2\), the absolute Harish-Chandra
parameter is \(\ell=k-1\). Its compact weights have
absolute value at least \(k\). In the compact picture
of normalized induction, restriction to \(K\) is induction
from \(K\cap M\). Therefore a contributing representation
that contains the \((5,5)\) type must meet its allowed
weights on each of the \(SO(2)\) factors of \(K\cap M\).
These factors rotate disjoint defining planes. All weights
of \((5,5)\) have individual absolute frequencies at most
five, from its Weyl weight polytope. Hence \(\ell\le4\)
on each two-block. Finite Levi components can exchange signs
and contribute finite characters. They do not change the
absolute frequency condition or trace Casimir.

The squared trace norm of \(\rho_G=(2,1,0,-1,-2)\) is
10. In an embedded two-block the diagonal coroot
\(\operatorname{diag}(1,-1)\) has trace norm squared two.
The trace Casimir on this tempered family is therefore
$$
 \Omega_{\mathrm{tr}}
    =\tfrac12\sum_{\text{two-blocks}}\ell_j^2-10-\|\nu\|_{\mathrm{tr}}^2
    \le 6.
$$
Here \(\nu\) is the real unitary split parameter, orthogonal
to the fixed Levi part. The full Plancherel theorem transfers
this upper bound to the regular \(L^2\) bundle. Since
\(\Omega_{\mathrm{tr}}=-4\nabla^*\nabla+35\) there, it gives
\(\nabla^*\nabla\ge29/4\).
[Herb and Wolf, Theorem 6.17](https://www.numdam.org/item/CM_1986__57_3_271_0.pdf).

## DH2. Weighted obstruction for global closed fields

For a closed pure field \(B\) and compactly supported scalar
cutoff \(\eta\), the eigenvalue equation gives
$$
 \|\nabla(\eta B)\|_2^2
  =7\|\eta B\|_2^2+\||d\eta|B\|_2^2 .
$$
The gap implies
\(\tfrac14\|\eta B\|_2^2\le\||d\eta|B\|_2^2\).
Suppose \(e^{-sr}B\in L^2\) for \(0\le s<1/2\) and
\(r=d(o,\cdot)\). Test with \(e^{-sr}\) times expanding
cutoffs whose derivative bounds tend to zero, and use
\(|dr|\le1\) almost everywhere. The inequality tends to
\(\tfrac14\|e^{-sr}B\|_2^2\le s^2\|e^{-sr}B\|_2^2\).
Thus every such weighted square-integrable field is zero.
For a nonzero closed field, the \(L^2\) mass of the
radius-\(R\) ball consequently has
\(\limsup_{R\to\infty}R^{-1}\log\int_{B_R}\|B\|^2\ge1\).
Otherwise shell summation would give weighted integrability
for some \(s<1/2\).

    python -X utf8 computation/geometric-unity/DH_h6_bundle_gap.py

The verifier checks the exact weights and trace factors on
all relevant ordered Levi block patterns. The spectral
proof uses the full regular representation and addresses this
pure type. It leaves open optimality of the bound and general
mixed fields.
[Review111](../../reports/reviews/REVIEW111-GU-DS-DH-DI-2026-11-14.md).
