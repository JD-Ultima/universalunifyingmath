# BC. Affine-invariant pure hook fields

The pure $H_6\Gamma$ branch has no affine-invariant full solution on the flat-base positive Z1 observer at $\lambda=1/5$, for any positive observer scales. Its derivative equation leaves one real line; the nonlinear grade-one equation excludes that line. The corresponding derivative test on the pure $H_4$ branch leaves six dimensions and does not decide its full equations.

Evidence: [construction and nonlinear certificate](../../computation/geometric-unity/BC_affine_hooks.py), [output](../../computation/geometric-unity/BC_affine_hooks.out), [explicit closed generator](../../computation/geometric-unity/BC_h6_closed_generator.json), and a [separate Clifford implementation](../../computation/geometric-unity/BC_quadratic_xcheck.py) with [output](../../computation/geometric-unity/BC_quadratic_xcheck.out). [BA](BA-hook-closure.md) supplies the derivative reduction.

## Scope and homogeneous model

Retain BA's true dualized bosonic equation, fixed parallel real $W_A$, compact real $(14,0)/\mathfrak u(128)$, real smooth fields and constant $\kappa\ne0$. The base is flat and positive; $\lambda=1/5$ and $c_H,c_V>0$ are constant. Restrict the even field to the restriction of an invariant tensor under
$$
G=GL^+(4,\mathbb R)\ltimes\mathbb R^4,
\qquad Y=G/SO(4).
$$
This assumption concerns the tensor field under the affine action, including its spin/Clifford identification. It is stronger than smoothness. Invariance does not mean that its Levi-Civita derivative is zero.

At the homogeneous identity, [AZ1](AZ-first-hook-classification.md) identifies
$$
\mathfrak p=\mathbb R^4\oplus\operatorname{Sym}_4,
\qquad (v,S)\longmapsto(v,-2S)
$$
in actual horizontal/vertical tangent coordinates. Normalize $c_H=2,c_V=1$. The metric and Nomizu maps are
$$
g(v,w)=2v\cdot w,\qquad
g(S,T)=4\left(\operatorname{tr}(ST)-\frac15\operatorname{tr}S\operatorname{tr}T\right),
$$
$$
\Lambda_S=0,\qquad \Lambda_vS=-Sv,\qquad
\Lambda_vw=\frac12\operatorname{sym}(vw)+\frac12(v\cdot w)I,
\quad \operatorname{sym}(vw)=\frac{vw^T+wv^T}{2}.
$$
A $G$-invariant covariant tensor is determined by an $SO(4)$-fixed value at the identity. Its covariant derivative is
$$
(\nabla_X B)(a;X_1,\ldots,X_q)
=-B(\Lambda_Xa;X_1,\ldots,X_q)
-\sum_jB(a;X_1,\ldots,\Lambda_XX_j,\ldots,X_q).
\tag{BC.1}
$$
Every slot is present. The ten vertical derivative directions vanish in this homogeneous description because $\Lambda_S=0$; the four horizontal ones generally do not.

## BC1. Complete linear reduction within the invariant ansatz

Use horizontal axes $H_i$, six off-diagonal symmetric matrices $O_{ij}=E_{ij}+E_{ji}$, then
$$
D_1=\operatorname{diag}(1,-1,0,0),\quad
D_2=\operatorname{diag}(1,1,-2,0),\quad
D_3=\operatorname{diag}(1,1,1,-3),\quad I.
$$
The Gram diagonal is
$$
(2,2,2,2,8,8,8,8,8,8,8,24,48,16/5).
$$
All tensor components are covariant. Contraction uses its inverse, not the identity matrix.

The exact construction gives:

| Form degree of the actual hook | Finite-subgroup candidates | Full $SO(4)$ tensor invariants | Hook invariants | Necessary derivative rank | Kernel dimension |
|---|---:|---:|---:|---:|---:|
| $C_4$ | 80 | 20 | 13 | 7 | 6 |
| $B_6$ | 216 | 30 | 23 | 22 | 1 |

The first derivative map is $C_4\mapsto d_\nabla(*C_4)$, equivalently $\delta_\nabla C_4$. The second is $B_6\mapsto d_\nabla B_6$. These are the necessary equations for $t=\kappa/24,u=0$ and $t=3\kappa/104,u=0$, respectively. R/V give a regular vector block and no other even kernel at either pair.

The calculation covers the complete affine-invariant spaces, rather than selected invariant examples:

1. Begin in the monomial horizontal/off-diagonal/$E_{ii}$ basis. The even coordinate-sign subgroup of $SO(4)$ fixes precisely the combined characters `0000` and `1111`. Retaining both allows orientation-dependent invariants. There are 1,736 degree-four and 5,268 degree-six raw components in these sectors.
2. Average each coordinate orbit over the 24 orientation-preserving signed representatives of coordinate permutations. On the even-sign fixed space these represent the quotient by the sign subgroup. Nonzero averages give independent bases of its finite-group invariants: 80 and 216. The helper explicitly checks each surviving average under every representative and exhausts the raw coordinates.
3. Change every covariant slot to the orthogonal $D/I$ basis. Impose the six infinitesimal generators of $\mathfrak{so}(4)$ on the whole tensor. Since $SO(4)$ is connected, their joint kernel is exactly its invariant subspace. Exact rational elimination leaves 20 and 30 independent tensors.
4. Impose both alternating and inverse-Gram contraction traces. This leaves 13 and 23 hooks. Apply (BC.1), then the full derivative map. Exact rational row reduction computes the complete nullspaces; nonzero square minors certify the matching lower ranks.

The recorded degree-four $7\times7$ minor is
$$-\frac{284375}{2918332558536081408},$$
and the degree-six $22\times22$ minor is
$$-132966527198300209152000.$$
The output supplies row/column labels and the complete nullspace amplitude vectors in the deterministically constructed bases. No floating-point or sampled-rank inference is used.

For degree four, the helper separately verifies the metric Hodge square, commutation with every nonzero Nomizu map and
$$d_\nabla(*C_4)=*\delta_\nabla C_4.$$
It uses star divided by the nonzero constant $\sqrt{\det g}$ to keep coefficients rational. The same factor occurs on both sides. All six reconstructed degree-four kernel tensors and the degree-six kernel tensor pass their entire derivative maps.

## BC2. The nonlinear obstruction on the $H_6\Gamma$ line

Let $B_*$ be the primitive generator reconstructed by BC1. The JSON certificate records its 1,358 nonzero covariant coefficients in the basis above. The primary helper constructs $B_*$ from the invariant equations and then compares every entry with that certificate; the certificate is not used to choose the kernel.

Every affine-invariant field on the pure branch has
$$T_o=i\frac{3\kappa}{104}\Phi^1,\qquad K=izB_*\Gamma,\quad z\in\mathbb R,$$
after its derivative equation has been imposed. In the normalized metric,
$$N_*:=\|B_*\|^2=\frac{14650875}{4096}>0.$$
Write $Q_1(K)/i$ for the grade-one quadratic tensor. Its diagonal values relative to the metric on the horizontal four-space, trace-free symmetric nine-space and scalar vertical line are, respectively,
$$
z^2N_*\left(
-\frac{13589}{13023},\quad
-\frac{162854}{117207},\quad
-\frac{8522}{13023}\right).
\tag{BC2.1}
$$
These three distinct irreducible $SO(4)$ summands have multiplicity one, so the invariant tensor is scalar on each. The primary helper checks all fourteen diagonal entries. The separate implementation computes one representative of each summand using different Clifford multiplication and nested-gradient code.

For clarity, the full cubic-gradient formula used to obtain the coefficients is as follows. For any real even Clifford field $C$, in an orthonormal basis,
$$
\frac{Q(iC)_a}{i}
=-\frac13\sum_{b<c}\left(
\{A_{abc},[C_b,C_c]\}
+[\{A_{abc},C_b\},C_c]
-[\{A_{abc},C_c\},C_b]
\right).
\tag{BC2.2}
$$
The three terms are the derivatives of the three field positions in the cubic potential. The alternating trivector vanishes for repeated indices. In the displayed non-unit covector basis, Clifford squares are $1/g_a$ and the right side acquires the factor $g_a$. Orthonormal conversion contributes $\sqrt{g_ag_bg_c}$ from the trivector and $1/\sqrt{g_bg_c}$ from the two field slots; converting the output supplies the remaining $\sqrt{g_a}$.

For $C=B_*\Gamma$, let $E$ be the product of the fourteen coordinate covectors. Then
$$E^2=-1/\det g,\qquad \Gamma=\sqrt{\det g}\,E.$$
Evaluate (BC2.2) on $B_*E$ and multiply the quadratic by $\det g$. No additional sign is inserted. As a normalization check, the moments obey
$$4q_H+9q_S+q_I=-\frac{52}{3}z^2N_*,$$
which independently agrees with T2's kernel-pairing identity. Both implementations also check the $312\Phi^1$ baseline normalization or compare selected coefficients with the complete cubic gradient. These tests fix conventions; the obstruction uses the exact values in (BC2.1).

The flat observer has $G^Y=15g_Y/(2c_V)$. On this pure branch the grade-one baseline is therefore the positive scalar
$$
L=\frac{15}{2c_V}+\kappa t+312t^2
=\frac{15}{2c_V}+\frac{15}{52}\kappa^2>0.
$$
In the normalized metric, subtracting the horizontal and trace-free vertical grade-one equations gives
$$
0=z^2N_*\left(-\frac{13589}{13023}+\frac{162854}{117207}\right)
=z^2N_*\frac{40553}{117207}.
$$
Thus $z=0$. The remaining grade-one equation would require $L=0$, a contradiction. Higher quadratic grades cannot cancel grade one; none is assumed solved. This proves nonexistence of a full field within the declared invariant ansatz.

## BC3. Positive scales and remaining gap

Arbitrary positive $c_H,c_V$ are recovered from the normalized model by
$$g(c_H,c_V)=c_VF_s^*g(2,1),\qquad F_s(x,y)=(sx,y),\quad s=\sqrt{c_H/(2c_V)}.$$
The map conjugates the affine action by $(A,b)\mapsto(A,sb)$, preserving the invariant ansatz. Constant homothety leaves the Levi-Civita connection and the hook/closure kernels unchanged up to the natural metric identifications. For a fixed seven-covariant $B_6$, each metric-relative quadratic moment scales by $c_V^{-7}$; the common positive factor preserves their unequal values. The positive scalar $L$ above remains the actual baseline. Hence the exclusion holds for every positive pair of scales at this fixed $\lambda$.

This is a classification of the affine-invariant pure $H_6\Gamma$ ansatz at $\lambda=1/5$ over a flat base. It is not an exclusion of arbitrary smooth $H_6\Gamma$ fields, other metric parameters, curved bases or coincident kernels. The six-dimensional invariant $H_4$ derivative kernel still requires its nonlinear and metric equations. The other hook copies, singular vector blocks, other real forms/operators and non-invariant odd fields retain their separate scope. No physical conclusion or literature priority is asserted.
