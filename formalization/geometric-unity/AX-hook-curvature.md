# AX. Curvature excludes the first pure hook on positive Z1 at $\lambda=1/5$

[AZ](AZ-first-hook-classification.md) subsequently extends the exclusion to every positive Z1 parameter $\lambda<1/4$. This note preserves the earlier $\lambda=1/5$ proofs and their exact inputs.

The first pure hook candidate in [V3](V-hook-intersections.md#v3-the-grade-five-equation-is-not-automatic) admits no smooth full solution over any smooth positive four-base with reconstructed positive Z1 at $\lambda=1/5$. The argument covers its entire 896-dimensional even kernel, extending [X](X-symmetric-hook.md)'s exclusion of a specified common-axis family.

The derivative equation forces covariant closure of a hook tensor. An exact, computer-assisted calculation proves that the flat-base observer's ambient curvature annihilates no nonzero such tensor. An exact radial rescaling then gives curvature injectivity on a dense open subset of each positive fibre for an arbitrary smooth base. Continuity forces the actual field to vanish on its own domain, contradicting the grade-one equation.

AX1–AX3 retain the derivative reduction, flat curvature certificate and flat-base exclusion. AX4–AX5 give the rescaling, fibre-density argument and general-base consequence. No field or base analyticity, continuation of an unknown field, or access of its domain to a radial point is assumed.

## Setting

Let $g$ be a smooth positive four-metric, and let $g_Y$ be its reconstructed Z1 observer metric with $\lambda=1/5$ and constant $c_H,c_V>0$. Use the true dualized equation, observer signature $(14,0)$, compact gauge real form $\mathfrak u(128)$, fixed parallel real $W_A$, real smooth fields, and constant real $\kappa\ne0$. On a nonempty open observer domain fix
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

The derivative reduction holds for any positive observer metric. AX2–AX3 concern the flat-base Z1 model; AX4–AX5 allow arbitrary smooth positive bases at the same $\lambda$. No boundary condition, compactness or parallel-field assumption is imposed on $B$.

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

## AX4. Radial rescaling and dense fibre injectivity

For every smooth positive four-base and every base point $p$, the curvature map AX2b is injective on an open dense subset of the positive fibre over $p$. Injectivity at every curved-base observer point is not asserted.

*Full ambient curvature limit.* Choose $g$-normal base coordinates $x$ centred at $p$, with $g(0)=I$. Write the Z1 metric as
$$
g_Y=c_Hy_{ij}\,dx^i dx^j+c_VG_{\lambda,y}(D_gy,D_gy),\qquad
D_gy=dy-(\Gamma_g^{\mathsf T}y+y\Gamma_g),
$$
where $\Gamma_g$ is a matrix-valued connection one-form. For $\varepsilon>0$ put
$$
d_\varepsilon(z)=\varepsilon z,\qquad
g_\varepsilon=\varepsilon^{-2}d_\varepsilon^*g,\qquad
\Psi_\varepsilon(z,h)=(\varepsilon z,\varepsilon^{-2}h). \tag{AX4a}
$$
The component functions of $g_\varepsilon$ are $g_{ij}(\varepsilon z)$. Constant metric scaling preserves the Levi-Civita connection, so
$$
(\Gamma_\varepsilon)_\mu(z)=\varepsilon(\Gamma_g)_\mu(\varepsilon z).
$$
The fibre factor in $\Psi_\varepsilon$ is the inverse-pullback factor for covariant metrics under $d_\varepsilon$. Direct substitution gives
$$
\Psi_\varepsilon^*(D_gy)=\varepsilon^{-2}D_{g_\varepsilon}h,\qquad
\boxed{\Psi_\varepsilon^*g_Y(g;c_H,c_V)=g_Y(g_\varepsilon;c_H,c_V).} \tag{AX4b}
$$
Indeed the horizontal factors $\varepsilon^{-2}$ and $\varepsilon^2$ cancel, while
$G_{\lambda,ah}(aU,aV)=G_{\lambda,h}(U,V)$ for $a>0$. This identity includes every mixed metric term and keeps $c_H,c_V,\lambda$ fixed.

At $(z,h)=(0,I)$, $\Gamma_\varepsilon=0$, so all the metrics on the right of AX4b have exactly the same positive tangent inner product
$$
c_H I_4\oplus c_VG_{1/5,I}.
$$
Also
$$
\partial_z^\beta\Gamma_\varepsilon(0)
=\varepsilon^{|\beta|+1}\partial_x^\beta\Gamma_g(0). \tag{AX4c}
$$
The full metric two-jet therefore converges to the constant-base model's two-jet. To include the vertical and mixed derivatives explicitly, set
$K_\mu=\Gamma_{\varepsilon,\mu}^{\mathsf T}h+h\Gamma_{\varepsilon,\mu}$.
Its three coordinate blocks are
$$
\begin{aligned}
(g_Y)_{\mu\nu}&=c_Hh_{\mu\nu}+c_VG_{\lambda,h}(K_\mu,K_\nu),\\
(g_Y)_{\mu A}&=-c_VG_{\lambda,h}(K_\mu,E_A),\\
(g_Y)_{AB}&=c_VG_{\lambda,h}(E_A,E_B),
\end{aligned} \tag{AX4d}
$$
where $E_A$ are coordinate symmetric matrices. On a fixed relatively compact positive neighbourhood of $h=I$, differentiate the difference from the constant-base metric up to total order two in $(z,h)$. Every resulting term is a smooth function of $h,h^{-1}$ multiplying $\Gamma_\varepsilon$, its first two base derivatives, or their products. Those connection jets tend to zero. The tautological horizontal term $c_Hh$ and the vertical block are already independent of $\varepsilon$. Smoothness of the base metric supplies the finite derivative bounds.

Curvature depends continuously on the metric two-jet. In one fixed orthonormal frame for the displayed common tangent inner product, its entire ambient tensor consequently converges to the flat-base tensor, with no omitted horizontal, mixed or vertical channel. Under $d\Psi_\varepsilon$, that frame becomes an orthonormal frame at
$(p,qg_p)$ for $q=\varepsilon^{-2}$: horizontal vectors scale by $q^{-1/2}$ and vertical vectors by $q$. Thus there is no changing-frame normalization left in the limit. The limiting curvature is AM1's full $SL(5,\mathbb R)/SO(5)$ curvature with the $-1/4$ normalization used in AX2.

The hook subspace at the comparison point is fixed because its metric is fixed. In fixed complete tensor bases, the flat-limit combined alternating, contraction and curvature map is injective by AX2, so choose **one nonzero global 1274-column maximal minor**. Keep that row selection; continuity preserves its nonvanishing for sufficiently small positive $\varepsilon$. The curved matrix need not retain any character-block symmetry; no separate continuation of the sixteen flat character blocks is used. Hence, for each fixed $p$, some sufficiently large **finite** $q$ gives injectivity at $(p,qg_p)$. No uniform choice of $q$ over base points is prescribed.

*Metric-only fibre density.* Fix $p$ and use all 1274 raw **covariant** components of $B$ in the coordinate tangent/coframe basis. The combined linear map is
$$
\mathcal L_p(y)B
=\big(\operatorname{Alt}B,\operatorname{Ctr}_{g_Y}B,
F^{T^*Y}\wedge B\big). \tag{AX4e}
$$
Its alternating row has no metric weight, its contraction uses the full inverse observer metric, and curvature acts on the coefficient covector by the negative transpose of its tangent action. These are the tensor conventions of AX2, also in nonorthonormal coordinates.

Every entry of $\mathcal L_p(y)$ is rational, hence real analytic, in the positive fibre matrix $y$. The coefficients involve only a finite smooth base jet at $p$: the Z1 metric, inverse metric, connection and curvature are obtained from the displayed metric formula by finite coordinate derivatives and matrix inversions. All denominators are nonzero on the positive cone under the standing assumptions. No analyticity in the base variables is required.

At the finite radial witness, changing from the comparison frame to the coordinate tensor bases is invertible. Thus select one global maximal coordinate minor of $\mathcal L_p(y)$ that is nonzero there. It is not the zero analytic function on that connected cone. Its nonzero locus is open dense, and there the combined map has full column rank. Equivalently, the curvature map restricted to the hook has zero kernel there. The chosen minor, like the radial witness, may depend on $p$. This proves the assertion. Only a matrix built from the metric has been used outside an actual solution domain; no field or field equation is continued. ∎

## AX5. No full solution over an arbitrary smooth positive base

Under the standing assumptions, there is no smooth full solution on a nonempty open positive-Z1 observer domain at $\lambda=1/5$ with
$t=3\kappa/136$, $u=0$, and arbitrary real even field.

*Proof.* The even equation again forces $T_e=iB\Gamma$ with $B\in H_2$. V3 and AX1 give $d_\nabla B=0$, hence $F^{T^*Y}\wedge B=0$. If the solution domain is $U$, then for each $p\in\pi(U)$ the intersection $U_p$ is open in its positive fibre. That slice may be disconnected: density on the whole connected positive cone still makes the injective locus dense in every open slice. Thus $B=0$ on a dense subset of $U_p$, and continuity gives $B=0$ on all of $U_p$. Doing this for every $p$ gives $B=0$ throughout $U$. The radial witnesses need not belong to $U$.

V3's full grade-one equation now becomes
$$
G^Y=-b_Hg_Y,\qquad b_H=\frac{201}{1156}\kappa^2>0.
$$
In dimension fourteen, trace reversal gives
$\operatorname{Ric}^Y=(b_H/6)g_Y$, so Ricci is positive definite.
But [AU1](AU-opposite-ricci-spectrum.md#au1-the-vertical-ricci-tensor) is a general positive-Z1 metric identity, independent of its later opposite-conic field hypotheses. Applied to the nonzero vertical Euler vector $E=y$, its curvature-action correction vanishes and gives
$$
\operatorname{Ric}^Y(E,E)
=-\frac54\operatorname{tr}I+\frac14(\operatorname{tr}I)^2
=-1,\qquad |E|^2=\frac{4c_V}{5}>0.
$$
This contradicts $\operatorname{Ric}^Y(E,E)=(b_H/6)|E|^2>0$. Both signs of $\kappa$ are included. ∎

As in AX3, the full grade-five equation remains a necessary equation; derivative compatibility and grade one already contradict each other. No grade-five classification or pointwise shape assumption for $B$ is used.

## Evidence and boundary

- [Derivative certificate](../../computation/geometric-unity/AX_hook_derivative.py) · [output](../../computation/geometric-unity/AX_hook_derivative.out): all 896 hook basis tensors and 14 derivative directions, **12,544 exact columns**. The basis contains two columns per distinct triple and twelve tracefree repeated-index columns per remaining index. It spans the whole hook, rather than sampled supports.
- [Rational curvature certificate](../../computation/geometric-unity/AX_hook_curvature.py) · [output](../../computation/geometric-unity/AX_hook_curvature.out): all sixteen exact rank and determinant witnesses, with the full constant-curvature control.
- [Independent curvature implementation](../../computation/geometric-unity/AX_hook_curvature_xcheck.py) · [output](../../computation/geometric-unity/AX_hook_curvature_xcheck.out): different matrix basis and elimination method, proving maximal rank by an integer minor modulo a prime.
- [Result record](../../metrics/geometric-unity/AX-hook-curvature.yaml). Internal reviews are linked there; external peer review remains outstanding.

AX4–AX5 are hand scaling, continuity and fibre-analytic arguments using the unchanged AX1–AX2 certificates and AU1's metric identity. They add no executable check. Smooth base metrics and fields suffice; the analytic object is the metric-only matrix as a function of a fixed fibre. Curvature injectivity is asserted on an open dense subset of each fibre, not at every observer point.

AX5 closes the first pure hook at $\lambda=1/5$ over arbitrary smooth positive four-bases. AZ supplies the subsequent full positive-parameter result. The other five pure hook cases, all hook–Lambda intersections, other operators and real forms remain outside that classification. In particular AW does not remove the coincident hook–Lambda5 candidates. Neither a symbol-kernel description nor a formal curvature constraint constructs a field. No literature-first or physical claim is made.
