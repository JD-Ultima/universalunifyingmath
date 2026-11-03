# AU. Ricci spectrum on the opposite Lambda5 conic

The vertical Ricci tensor of positive Z1 is its flat-base value plus a positive semidefinite correction of bounded rank. In the opposite-conic field sector, this excludes open Euler-collinear scalar patches and imposes a repeated Ricci eigenvalue and finite critical-value restrictions. [AW](AW-lambda5-conic-classification.md) subsequently excludes actual fields in that sector. AU1 remains a general metric identity; AU2–AU3 record the intermediate necessary field consequences.

Code: [AU helper](../../computation/geometric-unity/AU_opposite_ricci_spectrum.py) and [output](../../computation/geometric-unity/AU_opposite_ricci_spectrum.out). Record: [AU metadata](../../metrics/geometric-unity/AU-opposite-ricci-spectrum.yaml). [Independent internal Review51](../../reports/reviews/REVIEW51-GU-AU-2026-11-03.md) accepts AU1–AU3 and independently reproduces the helper. External peer review remains outstanding.

## Setting

AU1 concerns Z1 over a smooth positive $n$-metric, with $c_H,c_V>0$ and $\lambda<1/n$. Write $g_Y|_H=c_Hy$ and $g_Y|_V=c_VG_\lambda$, where $y$ is the positive fibre matrix and vertical vectors are symmetric matrices. The connection defining the splitting is the base Levi-Civita connection.

AU2–AU3 retain all [AR](AR-opposite-lambda5-geometry.md) assumptions: positive four-base and positive Z1; compact $(14,0)$ and $\mathfrak u(128)$; the true bosonic equation, fixed parallel real $W_A$, real smooth fields, constant nonzero $\kappa$, $\det M\ne0$, and an entire connected open single opposite $\Lambda^5$ conic domain away from coincident non-vector factors. AR supplies $w\ne0$, $\lambda=1/5$, and
$$c_0=4c_V/5,\qquad C=d\chi=24\,da/\delta,\qquad
\operatorname{Ric}=24\nabla C+Lg_Y-WC\otimes C,$$
$$L=2\operatorname{div}C+(\kappa^2U+V|C|^2)/6+W|C|^2/12,
\qquad \operatorname{Ric}(C^\sharp)=-C^\sharp/c_0.$$
The Euler field satisfies $E=y$, $|E|^2=c_0$, $\nabla E=P_H/2$ and $\operatorname{Ric}(E)=-E/c_0$. Curvature without a base subscript is that of $g_Y$.

## AU1. The vertical Ricci tensor

Put $\mathcal A(X,Y)=[X^H,Y^H]^V/2$ and define the skew horizontal operator $K_U$ by
$$\langle K_UX,Y\rangle=-\langle\mathcal A(X,Y),U\rangle.$$
For vertical $U,V$,
$$\boxed{\operatorname{Ric}(U,V)
=-\frac{n+1}{4}\operatorname{tr}(y^{-1}Uy^{-1}V)
+\frac14\operatorname{tr}(y^{-1}U)\operatorname{tr}(y^{-1}V)
+B(U,V),}\tag{AU1}$$
where $B(U,V)=\langle K_U,K_V\rangle_{\rm HS}$. Thus $B\ge0$, $B(E,\cdot)=0$, and $\operatorname{rank}B\le n(n-1)/2$. At $n=4$, $\lambda=1/5$,
$$\operatorname{Ric}|_{V\times V}=-g_Y|_V/c_0+B,\qquad\operatorname{rank}B\le6.$$

*Proof.* Use normal base coordinates at $x_0$. Trivialize the symmetric-matrix fibres in those coordinates and extend a local fibre-geodesic vector field $U(y)$ independently of $x$. Every fibre has the same coordinate metric $c_VG_\lambda$, and the fibres are totally geodesic by [N1](N-z1-neighbourhood.md). Hence $\nabla_UU=0$ nearby. For a projectable horizontal field $X$, $[U,X]=0$ along the entire fibre over $x_0$, where the base connection coefficients vanish.

Koszul, along this fibre, gives
$$\nabla_XU=\nabla_UX=L_UX,\qquad L_U=S_U+K_U,\qquad S_U=\tfrac12y^{-1}U.$$
There is no vertical component: the fibre metric and $U$ have zero horizontal coordinate derivative there. Consequently
$$R(X,U)U|_H=-(\partial_UL_U+L_U^2)X.$$
Fibre geodesics satisfy $y''=y'y^{-1}y'$ ([F1](F-zorro-reduction.md)), so $\operatorname{tr}(y^{-1}U)$ is constant. Since $S_U$ is self-adjoint and $K_U$ skew-adjoint,
$$\operatorname{tr}_H R(\cdot,U)U
=-\partial_U\operatorname{tr}L_U-\operatorname{tr}L_U^2
=-\tfrac14\operatorname{tr}((y^{-1}U)^2)+\|K_U\|^2.$$
Add the intrinsic fibre Ricci tensor $-nG_{1/n}/4$ from F1 and polarize. This proves (AU1) without a Riemannian-submersion assumption on Z1.

The map $U\mapsto K_U$ takes values in $\mathfrak{so}(H)$, proving the rank bound and positivity. The curvature-action bracket formula gives $\operatorname{tr}(y^{-1}\mathcal A(X,Y))=0$, since a base curvature endomorphism has trace zero. Thus $K_E=0$. ∎

## AU2. No open Euler-collinear patch

The set where $E$ and $C^\sharp$ are independent is open dense. Everywhere on the domain,
$$\dim\ker(\operatorname{Ric}^\sharp+I/c_0)\ge2.\tag{AU2}$$

*Proof.* Suppose $C^\sharp=qE$ on a nonempty open patch. Closedness of $C$ and $E^\flat$ gives $dq\wedge E^\flat=0$. Hence
$$\nabla C=\tfrac q2g_Y|_H+dq\otimes E^\flat.$$
On the nine-dimensional tracefree vertical space $V_0$, AR's metric identity gives $\operatorname{Ric}|_{V_0}=Lg_Y|_{V_0}$. AU1 implies
$$B|_{V_0}=(L+1/c_0)g_Y|_{V_0}.$$
A nonzero scalar multiple has rank nine, whereas $B$ has rank at most six. Therefore $L=-1/c_0$ and $B|_{V_0}=0$. Since $B(E,\cdot)=0$, $B=0$, so every $K_U$ and $\mathcal A$ vanish on this actual patch.

At each base point in its projection, each base curvature endomorphism is consequently skew for a fibre-open set of positive matrices $y$. This is a polynomial identity in $y$, so it holds for every positive $y$. An endomorphism skew for every positive scalar product must vanish. The projected base is flat on an open neighbourhood.

[N2](N-z1-neighbourhood.md) now gives $\operatorname{Ric}=-g_Y/c_0$. Comparing horizontal blocks in AR's equation yields $12q+L=-1/c_0$, so $q=0$. The resulting constant-$a$ open patch contradicts AR1. Thus the closed collinearity set has empty interior. On its open dense complement, N1 and AR5a give two independent eigenvectors with eigenvalue $-1/c_0$. Continuity of the thirteen-by-thirteen minors extends (AU2) everywhere. ∎

## AU3. Critical Hessian rank and finitely many critical values

At any critical point, write $\mathsf H=(\nabla C)^\sharp$. Then
$$\mathsf H=\eta P,\qquad
\eta=-\frac{\kappa^2U(a)/6+1/c_0}{24+2m}<0,\qquad
4\le m=\operatorname{rank}\mathsf H\le14,\tag{AU3a}$$
where $P$ is an orthogonal projection of rank $m$.

*Proof.* The [AT2](AT-opposite-profile-crossing.md) spectral argument applies at any $C=0$ point: differentiating AR5a and using AR0 gives $24\mathsf H^2+(L+1/c_0)\mathsf H=0$. Positivity and $U>0$ give the displayed nonzero eigenvalue, with $E$ in its eigenspace. At such a point,
$$\mathcal S:=\operatorname{Ric}^\sharp+I/c_0=24\eta(P-I)\ge0.$$
Its vertical compression is $B$ by AU1. If $U\in\ker B\subset V$, then $\langle\mathcal SU,U\rangle=0$. Positive semidefiniteness implies $\mathcal SU=0$, hence $U\in\operatorname{im}P$. Since $\dim\ker B\ge10-6=4$, $m\ge4$. ∎

For a fixed scale $t=\kappa^2c_0>0$, every critical scalar value lies among the roots of the eleven cubics
$$\boxed{t(208a+1)H_m(a)-360N_m(a)=0,\qquad m=4,\ldots,14,}\tag{AU3b}$$
where
$$\begin{aligned}
H_m={}&(4629031680m+48048349440)a^2\\
&+(48524736m+404087008)a+128681m+822893,\\
N_m={}&(43264m+3266432)a-2537m+20620.
\end{aligned}$$
Thus there are at most **33 distinct critical scalar values** before the conic and regularity exclusions. This bounds values, not the number of critical points.

*Proof.* AT excludes the sole zero of $D$. At any other critical value $a_0$, its actual Euler curve has the Taylor expansions used in AT2, with nonzero quadratic change in $a$. On a punctured half-curve this gives a bounded one-sided regular profile $y=c_0|C|^2$ with $y(a_0)=0$. The ordinary linear profile equation, whose leading coefficient is now nonzero, extends this scalar solution uniquely and gives
$$y'(a_0)=48c_0\eta/\delta(a_0).$$
Substitute into AR5b and eliminate $c_0\eta$ using (AU3a). With $q_m=24+2m$, the denominator-free identity is
$$t\,[q_m(\delta U'-2WU)-8DU]
=48D+12q_m(W-48).$$
Its two sides factor respectively as $-16t(208a+1)H_m/4275$ and $-128N_m/95$, giving (AU3b). Each polynomial has nonzero cubic coefficient $208t(4629031680m+48048349440)$, so each contributes at most three real values. No zero of $\delta$, $D$ or a vector factor is divided away. ∎

## Evidence and remaining question

Run `python -X utf8 -I computation/geometric-unity/AU_opposite_ricci_spectrum.py`. The recorded run passes 24 exact controls in 2.161 seconds with Python 3.12 and SymPy 1.14.0; Review51 reproduces them in 2.148 seconds. It checks fresh off-section Ricci tensors in dimensions $n=2,3,4$, using the accepted `F_jets` library without running its old suites: all 245 vertical tensor entries for four declared positive curved-base jets, the correction ranks and Euler kernel, and critical-curve factorizations. These finite controls support the hand proof; they do not classify all curvature tensors by sampling. The left-arc scalar candidates retained by the coefficient checks supply no actual field realization.

[AW](AW-lambda5-conic-classification.md) now excludes actual full fields in this single-conic setting. A permitted algebraic Ricci spectrum or cubic root never supplied a solution. Separate zero-order form feasibility, other singular kernels and the remaining slaved-vector Z1 curvature problem retain their own scope.