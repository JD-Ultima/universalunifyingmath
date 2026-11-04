# AZ. The first pure hook on every positive Z1 metric

The first pure $W_A$ hook, $t=3\kappa/136$ and $u=0$, admits no smooth full bosonic solution over any smooth positive four-base with reconstructed Z1 and any $\lambda<1/4$. The argument treats the entire 896-dimensional hook, without a support, parallel-field or analyticity assumption.

[AX](AX-hook-curvature.md) supplies the derivative reduction and the radial-fibre argument. The additional input here is an exact curvature certificate for the flat observer at every positive metric parameter, not only the $SL(5,\mathbb R)/SO(5)$ value. Metric-only fibre density then forces the actual hook field to vanish on its own domain. Its remaining grade-one equation contradicts the negative Euler Ricci direction.

Code: [AZ helper](../../computation/geometric-unity/AZ_flat_lambda.py) and [output](../../computation/geometric-unity/AZ_flat_lambda.out). Record: [AZ metadata](../../metrics/geometric-unity/AZ-first-hook-classification.yaml). Internal [Review59](../../reports/reviews/REVIEW59-GU-AZ-2026-11-03.md) and [Review60](../../reports/reviews/REVIEW60-GU-AZ-2026-11-03.md) accept the complete argument and reproduce the public helper, disclosing their earlier private contributions. External peer review remains outstanding.

## Setting

Let $g$ be a smooth positive four-metric, and let $g_Y$ be reconstructed Z1 with constant $c_H,c_V>0$ and $\lambda<1/4$. Keep the true dualized equation, observer signature $(14,0)$, compact gauge real form $\mathfrak u(128)$, fixed parallel real $W_A$, real smooth fields and constant $\kappa\ne0$. On an open observer domain fix
$$
T_o=it\Phi^1,\qquad t=3\kappa/136,\qquad u=0.
$$
[R2](R-even-operator.md) and [V3](V-hook-intersections.md#v3-the-grade-five-equation-is-not-automatic) give
$$
T_e=iB\Gamma,\qquad B\in H_2,\qquad
\det M=(1313/17)\kappa^4>0.
$$
Every other even block is invertible at this odd pair, so the even field is unrestricted before that algebraic reduction. Here $B$ is a real coefficient-covector-valued two-form, $\Gamma=\Gamma_{14}$, $\Gamma^2=-1$, and in an orthonormal frame
$$
\sum_a e^a\wedge B_a=0,\qquad \sum_a\iota_aB_a=0.
$$
The grade-12 matrices $B_a\Gamma$ are Hermitian; multiplication by $i$ makes the field anti-Hermitian. Put
$$
\mu=\frac1{1-4\lambda}>0.
$$
This parametrizes the whole positive-Z1 range, including negative $\lambda$.

## AZ1. The flat observer's affine model

For a flat base, the full observer is the affine homogeneous space
$$
(\operatorname{GL}^+(4,\mathbb R)\ltimes\mathbb R^4)/SO(4),
$$
with reductive tangent space $\mathfrak p=\mathbb R^4\oplus\operatorname{Sym}^2(\mathbb R^4)$. Normalize $c_H=2,c_V=1$. If $v,w$ are vectors and $S,T$ symmetric matrices, its tangent metric is
$$
\langle v,w\rangle=2v\cdot w,\qquad
\langle S,T\rangle=4\big[\operatorname{tr}(ST)-\lambda\operatorname{tr}S\operatorname{tr}T\big]. \tag{AZ1a}
$$
Write
$$
\operatorname{sym}(vw)=\frac{vw^{\mathsf T}+wv^{\mathsf T}}2,\qquad
K(v,w)=\frac12\left[\operatorname{sym}(vw)+\frac{\mu-1}{4}(v\cdot w)I\right].
$$
The Levi-Civita Nomizu maps and curvature are
$$
\Lambda_S=0,\qquad \Lambda_vS=-Sv,\qquad \Lambda_vw=K(v,w), \tag{AZ1b}
$$
$$
R(S,T)=-\operatorname{ad}_{[S,T]}|_{\mathfrak p},\qquad
R(S,v)=-\Lambda_{Sv},\qquad R(v,w)=[\Lambda_v,\Lambda_w]. \tag{AZ1c}
$$
In the middle expression $Sv$ indexes a connection map; it does not mean $\Lambda_Sv$.

*Proof and normalization.* The affine action is
$$
(A,b):(x,y)\longmapsto(Ax+b,A^{-\mathsf T}yA^{-1}).
$$
It preserves both terms of the flat Z1 metric and is transitive, with stabilizer $SO(4)$ at $(0,I)$. A symmetric generator $S$ has actual vertical tangent $-2S$, which explains the factor four in AZ1a. The reductive brackets are
$$
[S,T]\in\mathfrak{so}(4),\qquad [S,v]_{\mathfrak p}=Sv,\qquad [v,w]=0.
$$
Use the conventions
$$
\operatorname{Tor}(X,Y)=\Lambda_XY-\Lambda_YX-[X,Y]_{\mathfrak p},\qquad
R(X,Y)=[\Lambda_X,\Lambda_Y]-\Lambda_{[X,Y]_{\mathfrak p}}
-\operatorname{ad}_{[X,Y]_{\mathfrak{so}(4)}}.
$$
This is the classical homogeneous-space curvature identity; see Lorenzo-Naveiro's author-hosted [*Some useful formulae for homogeneous spaces*](https://jnaveiro.github.io/Documents/Lorenzo-Naveiro%20-%20Some%20useful%20formulae%20for%20homogeneous%20spaces.pdf), Theorem 6, §4, equation (25).
AZ1b has zero torsion. Its metric-skew property follows from
$\langle K(v,w),S\rangle=2v\cdot Sw$: indeed
$K-\lambda\operatorname{tr}(K)I=\operatorname{sym}(vw)/2$.
The invariant torsion-free metric connection is therefore the Levi-Civita connection. The curvature convention gives AZ1c.

For clarity, all ambient input/output types are
$$
\begin{aligned}
R(S,T)U&=-[[S,T],U], & R(S,T)w&=-[S,T]w,\\
R(S,v)T&=TSv, & R(S,v)w&=-K(Sv,w),\\
R(v,w)T&=-\tfrac14[P(v,w),T],&
R(v,w)z&=\tfrac{3-\mu}{8}P(v,w)z,
\end{aligned} \tag{AZ1d}
$$
where $P(v,w)=vw^{\mathsf T}-wv^{\mathsf T}$. No mixed or horizontal curvature is replaced by intrinsic fibre curvature.

Arbitrary positive scales satisfy
$$
g_Y(c_H,c_V)=c_V F_s^*g_Y(2,1),\qquad
F_s(x,y)=(sx,y),\quad s=\sqrt{c_H/(2c_V)}.
$$
This is a constant positive homothety and a coordinate change. It preserves the hook's zero contraction condition and the curvature-wedge kernel, so the normalization loses no parameters. ∎

The helper checks metric compatibility and first Bianchi symbolically, and obtains
$$
\operatorname{Ric}_{HH}(v,w)=-\frac\mu2v\cdot w,\qquad
\operatorname{Ric}_{VV}(S,T)=-5\operatorname{tr}(ST)+\operatorname{tr}S\operatorname{tr}T.
$$
At $\mu=5$, the map
$$
\Phi(v,S)=
\begin{pmatrix}
-2S+\frac25\operatorname{tr}(S)I_4&-v\\
-v^{\mathsf T}&\frac25\operatorname{tr}(S)
\end{pmatrix}
$$
matches the metric and every one of the 1,274 curvature columns with
$-\tfrac14[[\Phi X,\Phi Y],\Phi Z]$. This calibrates the full tensor against [AM1](AM-global-flat-models.md#am1-a-global-complete-realization). The $SL(5)/SO(5)$ identification is used only at $\mu=5$; the general model is the affine model above.

## AZ2. Exact curvature injectivity for every $\mu>0$

For every positive parameter in the flat model, the full ambient map
$$
H_2\longrightarrow T^*Y\otimes\Lambda^4T^*Y,\qquad
B\longmapsto F^{T^*Y}\wedge B \tag{AZ2a}
$$
has rank **896** and zero kernel.

*Complete parameter certificate.* Use four horizontal coordinate vectors; six vertical $O_{ij}=E_{ij}+E_{ji}$; three $D_k=\operatorname{diag}(1,\ldots,1,-k,0,\ldots,0)$ for $k=1,2,3$; and the vertical identity. The Gram diagonal is
$$
(2,2,2,2,\;8,8,8,8,8,8,\;8,24,48,16/\mu).
$$
Write the 1,274 raw components as
$B=\sum_a e^a\otimes\sum_{b<c}B_{a,bc}e^b\wedge e^c$.
All slots are covariant. Alternation has no metric weight; contraction uses the inverse Gram, including $\mu/16$ on the identity direction. Coefficient-covector curvature has matrix $-R_{ij}^{\mathsf T}$, and its wedge includes all six alternating shuffles, as in AX2c.

The combined alternating, contraction and curvature matrix is affine in $\mu$. To see the curvature degree before any rank calculation, the $\mu$-dependent part of $\Lambda_v$ maps horizontal vectors into the scalar vertical line and annihilates vertical vectors. Two such parts compose to zero, so AZ1c has degree at most one. The contraction is also affine. Clearing rational coefficient denominators and integer coefficient gcds uses only nonzero **constant** row factors, never a parameter factor.

The four coordinate-sign changes preserve the metric and AZ1d. They partition the matrix into sixteen characters. Coordinate permutations intertwine blocks of the same sign weight; they may mix the chosen diagonal basis, which does not affect ranks. Thus five representatives suffice:

| Sign weight | Multiplicity | Raw columns | Hook-constraint rank | Hook dimension |
|---|---:|---:|---:|---:|
| 0 | 1 | 94 | 18 | 76 |
| 1 | 4 | 85 | 25 | 60 |
| 2 | 6 | 85 | 25 | 60 |
| 3 | 4 | 66 | 22 | 44 |
| 4 | 1 | 66 | 22 | 44 |

The total hook dimension is 896. The metric decomposition
$V\otimes\Lambda^2V=\Lambda^3V\oplus V\oplus H_2$
gives alternating and contraction summands of dimensions 364 and 14. Thus their constraints have rank 378 for any positive metric, independently of a sampled value of $\mu$.

For each representative, select square row witnesses at $\mu=5$ and $\mu=2$, then compute their determinants **as polynomials in $\mathbb Q[\mu]$**. The [output](../../computation/geometric-unity/AZ_flat_lambda.out) records both full row lists and factorizations. The two selection parameters need not produce different witnesses; they only select rows. Up to nonzero numerical constants, the five determinant gcds are
$$
\begin{array}{c|l}
0&\mu^3(\mu-7)(\mu-1)^3\\
1&\mu^5(\mu-7)(\mu-3)^6(\mu-1)^6\\
2&\mu^2(\mu-7)^3(\mu-1)^2\\
3&\mu(\mu-7)(\mu-3)^5(\mu-1)^3(\mu+1)\\
4&1.
\end{array} \tag{AZ2b}
$$
The helper also checks that each monic gcd divides a power of
$\mu(\mu-1)(\mu-3)(\mu-7)(\mu+1)$.
Hence, outside $\{1,3,7\}$ in the positive parameter range, at least one maximal minor in each representative block is nonzero.

At each of $\mu=1,3,7$, exact elimination over $\mathbb Q$ gives the five combined ranks
$(94,85,85,66,66)$, all maximal. These fifteen cases exhaust the remaining possible positive failures. Their role is not sampled-rank extrapolation: the polynomial certificate first isolates this finite set.
The combined matrix therefore has rank 1,274 for every $\mu>0$. Restricting to the 896-dimensional hook proves AZ2a. No modular or floating-point rank is used. ∎

This is an exact computational linear-algebra proof. Its geometric model, coordinate symmetries and parameter-degree argument are the hand inputs. The helper is not an independent derivation of every geometric statement.

## AZ3. No first pure-hook full solution

Under the standing assumptions, no smooth full solution exists on a nonempty open observer domain for **any** smooth positive four-base and **any** $\lambda<1/4$.

*Radial limit at fixed parameter.* Fix $\lambda$, a base point $p$, and normal coordinates with $p=0$, $g(0)=I$. The exact rescaling in AX4 does not require $\lambda=1/5$:
$$
g_\varepsilon=\varepsilon^{-2}d_\varepsilon^*g,\qquad
d_\varepsilon(z)=\varepsilon z,\qquad
\Psi_\varepsilon(z,h)=(\varepsilon z,\varepsilon^{-2}h),
$$
$$
\Psi_\varepsilon^*g_Y(g;c_H,c_V,\lambda)
=g_Y(g_\varepsilon;c_H,c_V,\lambda),\qquad
(\Gamma_\varepsilon)_\nu(z)=\varepsilon(\Gamma_g)_\nu(\varepsilon z). \tag{AZ3a}
$$
At $(0,I)$ all these metrics have exactly the same positive inner product
$c_HI_4\oplus c_VG_{\lambda,I}$.
The connection and its first two base derivatives tend to zero. AX4d's explicit three metric blocks consequently give convergence of the entire observer two-jet, including mixed and vertical derivatives, in a fixed tangent frame. The limiting full ambient curvature is the flat affine model AZ1 with the **same** $\lambda$, not the $\mu=5$ symmetric-space tensor.

In fixed complete tensor bases, AZ2 makes the flat-limit combined map injective. Choose one nonzero global 1,274-column maximal minor and retain its row selection under the limit. Continuity gives an injective curvature map at some finite radial point $(p,qg_p)$, $q=\varepsilon^{-2}$. Curved jets need not preserve any sign-character symmetry; separate flat blocks are not continued.

*Density on the actual fibre slice.* At fixed $p$, the combined alternating, full inverse-metric contraction and coefficient-covector curvature matrix has rational, hence real-analytic, entries in the positive fibre matrix $y$. Its coefficients use a finite smooth base jet. Positivity makes all metric denominators nonzero. The finite radial witness supplies a global coordinate maximal minor not identically zero on the connected positive cone, so its nonzero set is open dense. The witness and minor may depend on $p$; injectivity at every curved observer point is neither required nor asserted.

By V3 and AX1 an actual field satisfies
$$
d_\nabla B=0,\qquad F^{T^*Y}\wedge B=0.
$$
For every $p$ in its projection, its open fibre slice $U_p$ meets the dense injective set, and $B$ vanishes on a dense subset of $U_p$. Continuity gives $B=0$ on the whole slice, even if it is disconnected, and hence throughout the observer domain. No domain access to the radial point, field continuation, or analyticity of the base or field has been used. Only the metric matrix is analytic along a fixed fibre.

*Euler contradiction.* With $B=0$, V3's grade-one equation is
$$
G^Y=-b_Hg_Y,\qquad b_H=\frac{201}{1156}\kappa^2>0.
$$
In fourteen dimensions trace reversal gives
$\operatorname{Ric}^Y=(b_H/6)g_Y$.
But [AU1](AU-opposite-ricci-spectrum.md#au1-the-vertical-ricci-tensor) is a general positive-Z1 metric identity, independent of its later opposite-conic field assumptions. For the nonzero vertical Euler vector $E=y$, its curvature correction vanishes and gives
$$
\operatorname{Ric}^Y(E,E)
=-\frac54\operatorname{tr}I+\frac14(\operatorname{tr}I)^2=-1,\qquad
|E|^2=4c_V(1-4\lambda)>0.
$$
This contradicts positive Einstein Ricci for either sign of $\kappa$. ∎

The full grade-five equation remains necessary; it is not replaced by invariant traces or assumed automatic. Derivative compatibility and grade one already contradict each other. No pointwise classification of its quadratic zero locus is needed.

## Evidence and boundary

Run the standalone helper from the repository root:

    python -X utf8 -I computation/geometric-unity/AZ_flat_lambda.py

The recorded author run completes in 43.867 seconds. It checks zero torsion, symbolic metric and curvature identities, the full $\mu=5$ calibration, affine matrix degree, all sixteen reference blocks, ten polynomial minors, their gcd cover and fifteen exact exceptional ranks. Public-review receipts belong to the metadata. The universal derivative input remains AX1's complete basis certificate; the model, symmetry, rescaling, density and continuity arguments are hand proofs. No independent symbolic or observer metric engine is claimed.

The conclusion covers only the first pure $W_A$ hook $t=3\kappa/136,u=0$ in the declared positive, compact, true-equation setting. The other five pure hooks, all hook–Lambda intersections, other operators, non-invariant odd fields, singular vector blocks and indefinite real forms remain outside it. In particular AW does not remove coincident hook–Lambda5 candidates. No construction, unrestricted hook classification, literature-first claim or physical confirmation follows.
