# DF. Square-integrable pure \((5,5)\) closure on \(SL(5,\mathbb R)/SO(5)\)

Every globally square-integrable pure H6 \((5,5)\) field on
the flat positive observer with \(d^\nabla B=0\) in
distributions vanishes. In particular the nonzero decaying
lifts of [DB](DB-h6-bounded-analytic-lifts.md) are outside
\(L^2\) for the invariant metric volume.

## DF1. The group Casimir on the bundle

[DE](DE-h6-elliptic-periodic.md) gives analytic regularity
and the identity \(\nabla^*\nabla B=7B\) in the trace metric.
Write \(G=SL(5,\mathbb R)\), \(K=SO(5)\), and let
\(\Omega_{\rm tr}\) be the Casimir for the invariant form
\(\operatorname{tr}(XY)\) on \(\mathfrak{sl}_5\).
In the trace-orthonormal tangent frame a group generator
\(P\in\mathfrak p\) gives twice the horizontal derivative,
since the unit tangent \(P\) is represented by \(\exp(tP/2)\).
On the right \(K\)-equivariant functions for the bundle, the
vertical contribution is one half the compact Casimir for
\(-\operatorname{tr}(KL)/2\). Its eigenvalue is \(70/2=35\).
Thus
$$
 \Omega_{\rm tr}=-4\nabla^*\nabla+35,\qquad
 \Omega_{\rm tr}B=7B .
$$
The associated-bundle identification gives an \(E_o\)-valued
function \(f(g)\) with \(f(gk)=k^{-1}f(g)\).
On normalizing Haar measure of \(K\) to one, \(f\) is in
\(L^2(G,E_o)\) with the same norm as the section. Each scalar
component satisfies the Casimir eigenvalue equation in the
regular \(L^2(G)\) representation. The two invariant
differential realizations of this central element agree.

## DF2. The Plancherel parameter has no atom at this level

The full group Plancherel theorem decomposes the regular
representation into tempered representations induced from
cuspidal parabolics \(MAN\), with discrete inducing data
on \(M\) and a unitary parameter \(\nu\in\mathfrak a^*\).
On each such family the measure has a density with respect
to Lebesgue measure in \(\nu\).
[Herb and Wolf, Theorem 6.17](https://www.numdam.org/item/CM_1986__57_3_271_0.pdf).

There is no group discrete series because the complex rank
of \(\mathfrak{sl}_5\) is four and the compact rank is two.
Thus every contributing parabolic is proper, with
\(\dim\mathfrak a>0\).
[Harish-Chandra, Theorem 13](https://doi.org/10.1007/BF02392813).

The infinitesimal character of normalized induction has the
fixed Levi component plus \(i\nu\). The trace form makes
the split center orthogonal to the Levi Cartan and is positive
on the split center. The quadratic Casimir scalar therefore
has the form \(c_\sigma-\|\nu\|^2_{\rm tr}\) for each discrete
inducing datum. Every fixed level has Lebesgue measure zero
on the positive-dimensional split parameter space. The
Plancherel transform intertwines the distributional Casimir
equation with multiplication by this scalar. Equivalently,
the central differential operator acts by the scalar on the
smooth vectors of each family, so distributional testing
localizes the transformed vector to the same level set.
Absolute continuity makes that vector zero. Applying this
to each scalar component of \(f\) proves the result.

    python -X utf8 computation/geometric-unity/DF_h6_l2_vanishing.py

The exact check verifies the Casimir constant and the positive
split-center Gram matrices for all fifteen proper standard
block parabolics.

The argument uses the full regular representation and the
calibrated bundle Casimir. It concerns the entire symmetric
space. General mixed hooks and fields on its quotients remain
open. External peer review remains outstanding.
