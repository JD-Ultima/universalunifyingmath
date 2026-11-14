# DE. Ellipticity and periodic analytic H6 lifts at type \((5,5)\)

On the flat positive observer, exterior covariant closure restricted
to the pure \((5,5)\) subbundle has a uniformly injective real symbol.
Every closed section of this subbundle is covariantly coclosed.
For a sufficiently long cyclic geodesic translation, every value
in the type has a bounded analytic closed lift invariant under the
translation. A lift through a nonzero value stays nonzero along
the corresponding escaping orbit.

## DE1. Exact symbols

Let \(V\) be the 14-dimensional tangent representation and \(E\)
the 286-dimensional real irreducible \((5,5)\) subbundle. For
\(\sigma:V^*\otimes E\to V^*\otimes\Lambda^7 V^*\) given by
exterior multiplication in the form slots, the exact \(SO(5)\)
decomposition and metric eigenvalues are as follows.

| Type | Dimension | \(\sigma^*\sigma\) | Cross Casimir |
| --- | ---: | ---: | ---: |
| \((7,5)\) | 1309 | 0 | 20 |
| \((6,5)\) | 715 | 0 | 4 |
| \((6,4)\) | 810 | \(26/27\) | \(-6\) |
| \((5,5)\) | 286 | \(7/3\) | \(-10\) |
| \((5,4)\) | 429 | \(5/6\) | \(-20\) |
| \((5,3)\) | 455 | \(223/210\) | \(-28\) |

The [verifier](../../computation/geometric-unity/DE_h6_elliptic_periodic.py)
checks one rational highest vector for each type. It checks every
derivative slice for hook and curvature, the pure-value Casimir
70 and both total simple-root equations. The dimensions sum
to \(14\cdot286\). Thus the decomposition is complete and
multiplicity free. Exact compact Hermitian norms give the
exterior-symbol eigenvalues. The two kernel types also give
zero contraction by an exact check on their highest vectors.

Normalize the compact Casimir by \(C(l,m)=l(l+3)+m(m+1)\),
using the Lie algebra metric \(-\operatorname{tr}(KL)/2\).
The cross Casimir \(C(V^*\otimes E)-C(V^*)-C(E)\)
is \(-2\sum_j\rho_V(K_j)\otimes\rho_E(K_j)\) in a compact
orthonormal basis. Its expectation on \(\xi\otimes b\) is
zero for every real \(\xi\) and complex \(b\), since each
compact generator on \(V\) is real skew. Let \(a,d\) be
the squared norms of its projections to the symbol kernel
and its orthogonal complement. Then \(0\ge4a-28d\), so
\(d\ge(a+d)/8\). The smallest nonzero symbol eigenvalue is
\(5/6\). Hence
$$
 \|\xi\wedge b\|^2\ge\frac5{48}\|\xi\|^2\|b\|^2
 \qquad(\xi\ \text{real}).
$$

## DE2. Covariant harmonicity

The subbundle \(E\) is parallel. Closure puts \(\nabla B\)
in the first two types of DE1, whose contraction is zero.
Thus \(\delta^\nabla B=0\) and
\((d^\nabla\delta^\nabla+\delta^\nabla d^\nabla)B=0\).
In the trace metric this operator is
\(\nabla^*\nabla+q(R)\), where in a normal frame
\(q(R)=-\sum_{ij}e^i\wedge\iota_j R_{ij}^{\rm tot}\).
The verifier applies the full curvature to the exact BO
highest tensor and obtains \(q(R)b=-7b\). Equivariance
extends that identity to all of \(E\). Therefore
\(\nabla^*\nabla B=7B\). In particular, for a square-integrable
closed pure field, the cutoff identity
\(\|\nabla(\eta B)\|_{L^2}^2=7\|\eta B\|_{L^2}^2+
\|d\eta\otimes B\|_{L^2}^2\) first gives finite derivative
norm on exhaustion. Letting the cutoffs tend to one gives
\(\|\nabla B\|_{L^2}^2=7\|B\|_{L^2}^2\).
[DF](DF-h6-l2-vanishing.md) excludes nonzero \(L^2\)
fields of this pure type by the group Plancherel theorem.

## DE3. Periodic bounded lifts

[DB](DB-h6-bounded-analytic-lifts.md) gives a linear family
\(F_v\) with \(F_v(o)=v\) and
$$
 \|F_v(Y)\|\le C\|v\|\exp(-\kappa d(o,Y)),\qquad
 C=286(1+\sqrt{3/2}),\qquad \kappa=\frac5{24\sqrt{20}}.
$$
Indeed in the trace metric \(d(o,Y)=\|\log Y\|_F\).
If the smallest log eigenvalue is \(-a\), the sum-zero
condition bounds the squared sum of all five by \(20a^2\).
This converts the spectral decay estimate of DB to the
displayed distance bound.

Choose a geodesic through \(o\) and a transvection \(\gamma\)
of length \(L\) along it. Put \(u=\exp(-\kappa L)\) and choose
\(L\) so \(u=1/(8C+1)\). The bilateral sum
$$
 P_w=\sum_{n\in\mathbb Z}\gamma^n\cdot F_w
$$
converges locally uniformly, since for \(d(o,x)\le R\)
each term is bounded by \(C\|w\|\exp(\kappa R)u^{|n|}\).
Nearest-point projection to the geodesic is nonexpansive
on this Hadamard symmetric space. The distance from \(x\)
to \(\gamma^n o\) is consequently at least \(|s(x)-nL|\)
for the projected real geodesic coordinate \(s(x)\).
Summing the geometric series gives
\(\|P_w\|_\infty\le C\|w\|(1+u)/(1-u)\).

The sum is \(\gamma\)-invariant by shifting the index.
Local uniform convergence gives distributional closure and
coclosure, as well as all algebraic hook conditions.
The covariant Laplace operator in DE2 is analytic and
elliptic, so \(P_w\) is real analytic by interior analytic
elliptic regularity. [Morrey and Nirenberg](https://doi.org/10.1002/cpa.3160100204),
[distributional formulation](https://ocw.mit.edu/courses/18-757-representations-of-lie-groups-fall-2023/mit18_757_f23_lec06.pdf).

At the base point the evaluation \(T_L:w\mapsto P_w(o)\)
satisfies \(\|T_L-I\|\le2Cu/(1-u)=1/4\). Hence \(T_L\)
is invertible and \(P_{T_L^{-1}v}(o)=v\). For \(v\ne0\)
the field has norm \(\|v\|\) at every \(\gamma^n o\), and
therefore does not vanish at infinity. This lies beyond the
decaying averages in DB. The [BQ](BQ-h6-invariant-moments.md)
invariant grade-five coefficient still excludes the pure
type by itself from the nonlinear equations.

    python -X utf8 computation/geometric-unity/DE_h6_elliptic_periodic.py

The verifier checks the exact norm ratios, curvature and
periodization normalization. The uniform convergence,
analytic regularity and interpolation arguments are given
above. External peer review remains outstanding.

[Review109](../../reports/reviews/REVIEW109-GU-DE-2026-11-14.md)
records the independent checks.
