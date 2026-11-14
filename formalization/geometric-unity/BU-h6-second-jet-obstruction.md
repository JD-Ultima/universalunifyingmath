# BU. A second-derivative obstruction to the H2 component of a closed H6 field

Every smooth closed real $H_6$ hook on the flat positive observer has
zero component in the 14-dimensional type $(2,0)$ of the curvature
kernel of [BO](BO-h6-full-curvature-kernel.md). Its values lie in the
remaining 6,820-dimensional subbundle. In contrast to an ansatz
restriction, this applies to any covariantly closed hook, and therefore
to any full field on the branch if one exists.

The [verifier](../../computation/geometric-unity/BU_h6_h2_prolongation.py)
reconstructs the exact seven-dimensional space of equivariant test
jets and verifies a nonzero rational obstruction.

## BU1. Differentiated closure

Let $V=\operatorname{Sym}^2_0(\mathbb R^5)$ in the normalized observer
frame. By [BA](BA-hook-closure.md), a closed hook $B$ takes values in
the full linear curvature kernel
$\mathcal K=\ker(\operatorname{Alt},\operatorname{Ctr},\mathcal R)$.
The metric and curvature are parallel on the symmetric space, so
$\mathcal K$ is a parallel subbundle. Let $b$ be the value at a
normal-frame point, and put
$$
h_{ij}=\tfrac12(\nabla_i\nabla_jB+\nabla_j\nabla_iB).
$$
The symmetric tensor $h$ belongs to $\operatorname{Sym}^2V^*\otimes
\mathcal K$. The antisymmetric second derivative is fixed by the
full tensor curvature $R^{\mathrm{tot}}_{ij}$ acting on the external
covector and all six exterior slots. Differentiating $d^\nabla B=0$
gives the linear equation
$$
(A h)_{i,a} = \sum_j e^j\wedge (h_{ij})_a
= -\tfrac12\sum_j e^j\wedge(R^{\mathrm{tot}}_{ij}b)_a
=F(b)_{i,a}.
$$
Here the curvature in $F$ is that of the seven-slot value tensor.
The earlier curvature-wedge condition on $b$ uses the coefficient
covector connection. Both conventions follow from the same
symmetric-space curvature.

The maps $A,F$ commute with $SO(5)$. Therefore the kernel of the
induced map $b\mapsto[F(b)]$ into the cokernel of $A$ is an invariant
subspace of $\mathcal K$. To test its unique type $V_{(2,0)}=V$, it
suffices to ask whether $F$ restricted to that type has an
equivariant lift through $A$. Indeed, if the restriction had image
in $\operatorname{im} A$, an invariant complement of $\ker A$ would
give such a lift by complete reducibility. The obstruction has no
cross-cancellation with inequivalent input types.

## BU2. Seven equivariant lifts to test

The tangent module is self-dual, so
$$
\operatorname{Hom}_{SO(5)}\big(V,\operatorname{Sym}^2V^*\otimes\mathcal K\big)
\simeq
\operatorname{Hom}_{SO(5)}\big(\mathcal K,\operatorname{Sym}^2V\otimes V\big).
$$
The character of the latter target, of dimension 1,470, is
$$
(6,0)+(5,1)+2(4,2)+2(4,0)+(3,2)+3(3,1)+2(2,2)
+(2,1)+4(2,0)+(1,1)+(0,0).
$$
The multiplicities are obtained by summing tangent weights
on the symmetric pair and the third factor and subtracting exact
Weyl characters. The only types shared with BO are $(5,1),(3,2),
(2,1),(2,0)$, each with one kernel copy. The corresponding target
multiplicities are 1, 1, 1, and 4. Thus every possible equivariant
lift is accounted for by seven parameters in these four kernel types.

The verifier generates their 414 dimensions by negative root
actions on BO's four exact highest vectors. In the resulting
symmetric-jet space there are 862 columns of weight $(2,0)$.
The two positive simple root equations have rank 855. Their
seven exact independent nullvectors are the highest jets, and
therefore give all intertwiner images of the chosen highest $b$.

## BU3. Cokernel certificate

In the split rational coframe the following eight output rows
of $A$ suffice:
$$
(0,0,254),(0,0,1468),(0,0,2001),(0,0,2020),
(0,0,3892),(0,0,4572),(0,0,5976),(0,1,253).
$$
Each label records the derivative slot, coefficient covector slot,
and seven-form mask. The verifier computes the seven column images
and the curvature source on the same rows. After exact nonzero
column rescalings, the augmented determinant is 241629972480.
Equivalently, the row functional
$$
(915357947,-4246099842,0,-228528579,-717216798,-1258489440,
1006791552,0)
$$
annihilates all seven closure columns and pairs with $F(b)$ to
give 1006791552. Consequently $F$ has no equivariant lift on this
type. Its cokernel map is injective on the unique irreducible
14-dimensional copy and every closed field has zero projection
to that copy.

Positive constant scales transfer by the established dilation and
homothety, preserving covariant closure and the kernel types.
BR gives analytic closed lifts on a 6,534-dimensional submodule,
leaving the 286-dimensional type $(5,5)$ between the explicit
lower bound and the 6,820-dimensional upper bound.
[CF](CF-h6-fourth-linear-jet.md) proves finite compatibility
through the fourth Taylor order for that type. The later
[CU construction](CU-h6-global-analytic-closure.md) gives
global analytic lifts and attains this exact upper bound.
The nonlinear equations remain open. External peer review
remains outstanding.
