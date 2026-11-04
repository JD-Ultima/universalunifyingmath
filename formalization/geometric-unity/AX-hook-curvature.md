# AX. Curvature excludes the first pure hook over flat-base Z1

The first pure hook candidate in [V3](V-hook-intersections.md#v3-the-grade-five-equation-is-not-automatic) admits no smooth full solution over the positive flat-base Z1 model at $\lambda=1/5$. The argument covers its entire 896-dimensional even kernel, extending [X](X-symmetric-hook.md)'s exclusion of a specified common-axis family.

The derivative equation forces covariant closure of a hook tensor. Full ambient observer curvature then forces that tensor to vanish, while the grade-one equation requires a nonzero quadratic term. The curvature step is an exact, computer-assisted linear-algebra proof, with a separate implementation and basis.

## Setting

Use the reconstructed true dualized equation, signature $(14,0)$, compact gauge real form $\mathfrak u(128)$, fixed parallel real $W_A$, positive constant $c_H,c_V$, and constant real $\kappa\ne0$. Fix
$$
T_o=it\Phi^1,\qquad t=3\kappa/136,\qquad u=0.
$$
[R2](R-even-operator.md) and V3 give exactly
$$
T_e=k=iB\Gamma,\qquad B\in H_2,\qquad
\det M=(1313/17)\kappa^4>0.
$$
Here $B_a$ is a real two-form, $\Gamma=\Gamma_{14}$, $\Gamma^2=-1$, and
$$
\sum_a e^a\wedge B_a=0,\qquad \sum_a\iota_a B_a=0. \tag{AX0}
$$
These are the alternating and contraction conditions in an orthonormal frame. The grade-12 matrices $B_a\Gamma$ are Hermitian, so $k$ is anti-Hermitian. Every other even block is invertible at this odd pair; no additional even component is omitted.

The derivative reduction below holds for any positive observer metric. The curvature calculation and exclusion require a flat positive four-base, Z1 and $\lambda=1/5$. No analyticity, boundary condition, compactness or parallel-field assumption is imposed on $B$.

## AX1. The full derivative equation is covariant closure

Treat $B$ as a $T^*Y$-valued two-form. Define
$$
D=\sum_b\nabla_bB_b,\qquad
(d_\nabla B)_a=\sum_b e^b\wedge\nabla_bB_a.
$$
The derivative acts on all three tensor indices, including the coefficient covector $a$. With $\mathcal Dk={*}^{-1}\circledast(d_0k)$ in V3's convention,
$$
(\mathcal Dk)_a\Gamma
=2\big[(d_\nabla B)_a-e^a\wedge D\big]. \tag{AX1a}
$$
Consequently
$$
\mathcal Dk=0\quad\Longleftrightarrow\quad d_\nabla B=0,
\qquad\text{and then }D=0. \tag{AX1b}
$$

*Proof.* Work at a normal-frame point. Put $H_{bc}=\nabla_bB_c-\nabla_cB_b$ and $\gamma_{abc}=e_a\wedge e_b\wedge e_c$ in the Clifford algebra. Repeated indices give zero. The real operator and field phases give
$$
(\mathcal Dk)_a\Gamma=\sum_{b<c}[\gamma_{abc},H_{bc}].
$$
For a bivector $H$,
$$
[\gamma_{abc},H]=2\big[e^a\wedge e^b\wedge\iota_cH
-e^a\wedge e^c\wedge\iota_bH
+e^b\wedge e^c\wedge\iota_aH\big].
$$
Differentiate both parallel identities AX0. The first two summed terms reduce to $-2e^a\wedge D$; the last reduces to $2(d_\nabla B)_a$. This proves AX1a without dropping derivatives of any component.

The same identities give
$$
\sum_a\iota_a(d_\nabla B)_a=D,\qquad
\sum_a\iota_a(e^a\wedge D)=12D.
$$
Contracting AX1a yields $-22D$. Thus its vanishing forces $D=0$ and then $d_\nabla B=0$. Conversely, contraction of $d_\nabla B=0$ gives $D=0$, so AX1a vanishes. ∎

For a nonzero covector $\xi$, the principal-symbol kernel is $\operatorname{Sym}^2_0(\xi^\perp)$, of dimension 90: $\xi\wedge B_a=0$ gives $B_a=\xi\wedge v_a$, and AX0 makes the tangent matrix symmetric and tracefree. This describes derivative symbols, not the pointwise shape of an actual field. X's common-axis obstruction alone therefore does not prove the result below.

## AX2. The ambient curvature map is injective on the hook

For any solution of AX1,
$$
0=(d_\nabla)^2B=F^{T^*Y}\wedge B. \tag{AX2a}
$$
The curvature acts on the coefficient covector index. This is the usual exterior-covariant-derivative identity; the form-slot contributions cancel by the torsion-free first Bianchi identity.

For the flat-base positive Z1 model at $\lambda=1/5$, the linear map
$$
H_2\longrightarrow T^*Y\otimes\Lambda^4T^*Y,
\qquad B\longmapsto F^{T^*Y}\wedge B \tag{AX2b}
$$
has rank **896**, hence zero kernel, at every observer point.

*Exact certificate.* [F7](F-zorro-reduction.md) and [AM1](AM-global-flat-models.md#am1-a-global-complete-realization) identify this full observer metric, for any $c_H,c_V>0$, with $SL(5,\mathbb R)/SO(5)$. At the identity, take
$$
\mathfrak p=\operatorname{Sym}^2_0(\mathbb R^5),\qquad g(X,Y)=c_V\operatorname{tr}(XY).
$$
In AM's positive-matrix tangent coordinates,
$$
R(X,Y)Z=-\tfrac14[[X,Y],Z].
$$
This also follows directly by differentiating AM1's connection. The helper uses four times this curvature, $\widehat R=-[[X,Y],Z]$; multiplying the map by a nonzero scalar leaves its kernel unchanged. The factor $c_V$ likewise does not change the connection or the hook's zero contraction condition. No vertical-fibre curvature is substituted.

Use ten off-diagonal matrices $O_{ij}=E_{ij}+E_{ji}$ and four diagonal matrices
$$
D_k=\operatorname{diag}(\underbrace{1,\ldots,1}_{k},-k,0,\ldots,0),\qquad 1\le k\le4.
$$
Their trace-metric diagonal is $(2,\ldots,2,2,6,12,20)$. Construct $\widehat R_{ij}$ by matrix commutators and rational coordinate extraction. On coefficient covectors the curvature matrix is $-\widehat R_{ij}^{\mathsf T}$.

There are $14\binom{14}{2}=1274$ raw covariant components $B_{a,bc}$, $b<c$. Impose the 364 alternating equations and 14 metric-weighted contraction equations. For each $i<j<k<l$, impose all 14 components of
$$
F_{ij}B_{kl}-F_{ik}B_{jl}+F_{il}B_{jk}
+F_{jk}B_{il}-F_{jl}B_{ik}+F_{kl}B_{ij}=0. \tag{AX2c}
$$

Diagonal sign changes in $SO(5)$ partition the raw variables into sixteen character blocks. Each row has only one character, checked during construction. Exact rational elimination and nonzero square determinant witnesses give:

| Character type | Blocks | Raw columns per block | Hook-constraint rank | Hook dimension | Combined rank |
|---|---:|---:|---:|---:|---:|
| Empty subset | 1 | 94 | 18 | 76 | 94 |
| Two-element subset | 10 | 85 | 25 | 60 | 85 |
| Four-element subset | 5 | 66 | 22 | 44 | 66 |
| Total | 16 | 1274 | 378 | 896 | 1274 |

All entries are exact. The [recorded output](../../computation/geometric-unity/AX_hook_curvature.out) retains each selected square row list and its nonzero determinant. Row denominators are cleared by nonzero scalars. The combined matrix has full column rank, so its real kernel is zero. Homogeneity transfers this to every observer point. ∎

A separate implementation uses the nonorthogonal diagonal basis $E_{ii}-E_{44}$, direct matrix-entry coordinates and its full inverse metric. Its combined integer matrix has maximal column rank 1274 modulo the prime 104729. This also proves real injectivity: a maximal integer minor nonzero modulo a prime is nonzero over the integers. No conclusion from a deficient modular rank is used.

Both implementations check metric skewness, first Bianchi and $\widehat{\operatorname{Ric}}=-5\operatorname{tr}(XY)$. The rational implementation also verifies that constant-sectional-curvature tensors annihilate the entire 896-dimensional hook. Thus injectivity is a property of this ambient curvature, not a general consequence of the hook identities.

## AX3. No full solution on this pure branch

There is no smooth solution of the reconstructed full equation on a nonempty open subset of this flat-base positive Z1 model with $t=3\kappa/136$, $u=0$ and arbitrary real even field.

*Proof.* The even equation forces $T_e=iB\Gamma$ with $B\in H_2$. V3's grade-eleven equation gives AX1, hence AX2a. Injectivity gives $B=0$ throughout the domain. But F7 and V3's grade-one equation require
$$
\frac{[Q(iB\Gamma)]_1}{i}
=-\left(\frac{15}{2c_V}+\frac{201}{1156}\kappa^2\right)g_Y.
$$
The left side is zero; the coefficient on the right is strictly positive. This is a contradiction for either sign of $\kappa$. ∎

The grade-five equation has not been replaced by its scalar traces or assumed automatic. The necessary derivative and grade-one equations already contradict each other; imposing grade five cannot repair that contradiction.

## Evidence and boundary

- [Derivative certificate](../../computation/geometric-unity/AX_hook_derivative.py) · [output](../../computation/geometric-unity/AX_hook_derivative.out): all 896 hook basis tensors and 14 derivative directions, **12,544 exact columns**. The basis contains two columns per distinct triple and twelve tracefree repeated-index columns per remaining index. It spans the whole hook, rather than sampled supports.
- [Rational curvature certificate](../../computation/geometric-unity/AX_hook_curvature.py) · [output](../../computation/geometric-unity/AX_hook_curvature.out): all sixteen exact rank and determinant witnesses, with the full constant-curvature control.
- [Independent curvature implementation](../../computation/geometric-unity/AX_hook_curvature_xcheck.py) · [output](../../computation/geometric-unity/AX_hook_curvature_xcheck.out): different matrix basis and elimination method, proving maximal rank by an integer minor modulo a prime.
- [Result record](../../metrics/geometric-unity/AX-hook-curvature.yaml). Internal reviews are linked there; external peer review remains outstanding.

AX3 closes only the first pure hook over this flat-base metric. The other five pure hook cases, all hook–Lambda intersections, curved bases, other $\lambda$, operators and real forms remain outside the exclusion. AX1 is a useful necessary equation in the more general positive first-hook problem; neither its symbol kernel nor curvature compatibility constructs a field. No literature-first or physical claim is made.
