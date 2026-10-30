# K. The Ricci blocks of the Zorro metric Z1 along the section, for every $n$ (a proof of Theorem F4(c))

Claims: GU-054 (context: GU-006, GU-021). Gaps: GU-G14. Code: `computation/geometric-unity/K_z1_ricci.py` → `.out` (sections K0–K5; exact: Python Fractions and sympy rational-function fields; fixed seed `random.Random(20261030)`; 480 s (8 min) on one core; a second run gave identical output apart from the timing annotations, 433 s). The script imports no repository code. Its curvature routine, jets of the Zorro metrics and diagram calculus are written from the definitions, independently of `F_jets.py`. Record: `metrics/geometric-unity/K-z1-ricci.yaml`. Builds on F (R-Z, F2, F3, F4(a)(b), L-F1). Layers: (b) reconstruction (that of F, unchanged), (c) proven.

## K.0 Question and setting (layer b)

Theorem F4(c) gives the Ricci blocks of the draft-literal reading Z1 along the section. Until now it was proven only for $n=1,2$ (E1, R2) and verified exactly for $n=3,4,5$ (R3). The general-$n$ argument was a sketch. Every Z1 statement downstream ($\Lambda_P$ of F5, the Z1 half of F6, Cor. I2.1 and I7 in Z1, J3–J4 in Z1) carried that status. This note proves F4(c) for every $n\ge1$, every signature of $g$, every $\lambda\neq1/n$ and every $c_H,c_V\neq0$.

**Setting (R-Z of F).** $X^n$ carries a metric $g$ of signature $(p,q)$; $Y=\mathrm{Met}_{(p,q)}(X)$ has coordinates $(x^\mu,y_A)$, $A=(i\le j)$, and $Dy=dy-\Psi\,dx$ with $\Psi_{(ij)\mu}=\Gamma^k_{\mu i}y_{kj}+\Gamma^k_{\mu j}y_{ik}$ (the Levi-Civita connection of $g$ induced on $\mathrm{Sym}^2T^*X$). Then
$$h:=\mathrm{Z2}=c_H\,g(\pi_*,\pi_*)+c_VG_\lambda|_y(Dy,Dy),\qquad \hat h:=\mathrm{Z1}=c_H\,y(\pi_*,\pi_*)+c_VG_\lambda|_y(Dy,Dy),$$
$$\delta:=\hat h-h=c_H\,(y-g)(\pi_*,\pi_*),\qquad G_\lambda|_y(U,W)=\mathrm{tr}(y^{-1}Uy^{-1}W)-\lambda\,\mathrm{tr}(y^{-1}U)\,\mathrm{tr}(y^{-1}W).$$
Curvature conventions are those of F: $\mathrm{Ric}_{bd}=R^a{}_{bad}$, and the unit sphere has $\mathrm{Ric}=(n-1)g$.

**Notation at a point.** Fix $p=s(x)\in s(X)$ and write $g=g_x$. The connection map $K:=Dy:T_pY\to\mathrm{Sym}^2T^*_xX$ identifies $V_p=\ker\pi_*$ with $\mathrm{Sym}^2T^*_xX$, and $\pi_*$ identifies $H_p=\ker K$ with $T_xX$. Letters $X,Y$ denote horizontal vectors and $U,W$ vertical ones. In the adapted frame ($e_\mu=\partial_\mu+\Psi_{A\mu}\partial_A$, $e_A=\partial_A\leftrightarrow E_A$) one has $h_p=\hat h_p=c_Hg\oplus c_VG_\lambda|_g$ (F2.3). $G_0,G_{1/n},G_{1/(n+1)}$ below are evaluated at $y=g$, and $\mu:=\lambda/(1-n\lambda)$.

**Statement to prove (F4(c)).** On $s(X)$,
$$\mathrm{Ric}^{\rm Z1}_{HH}=\mathrm{Ric}_g+\frac{c_H}{4c_V(n\lambda-1)}\,g,\qquad \mathrm{Ric}^{\rm Z1}_{HV}=0,\qquad \mathrm{Ric}^{\rm Z1}_{VV}=-\tfrac n4G_{1/n}-\tfrac14G_0=-\tfrac{n+1}4G_{1/(n+1)} .$$

**Inputs from F (all proven there).**
- (F3-i) $\pi:(Y,h)\to(X,c_Hg)$ is a pseudo-Riemannian submersion with horizontal distribution $\ker D$, nondegenerate fibres $(\mathrm{Met}_{(p,q)}(T_xX),c_VG_\lambda)$ and $T\equiv0$.
- (F3-ii) Its $A$-tensor vanishes on $s(X)$.
- (F3-iii) On $s(X)$: $\mathrm{Ric}^{\rm Z2}_{HH}=\mathrm{Ric}_g$, $\mathrm{Ric}^{\rm Z2}_{HV}=0$, $\mathrm{Ric}^{\rm Z2}_{VV}=-\tfrac n4G_{1/n}$.
- O'Neill's identities for a pseudo-Riemannian submersion with nondegenerate fibres (L-F1), for basic $X,Y$ and vertical $U,W$: $\mathcal H\nabla_XY=(\nabla^g_XY)^H$, $\mathcal V\nabla_XY=A_XY$, $\mathcal H\nabla_XU=A_XU$, $\mathcal H\nabla_UX=A_XU$, $\mathcal V\nabla_UX=T_UX$, $\mathcal H\nabla_UW=T_UW$, $\mathcal V\nabla_UW=\hat\nabla_UW$ (the fibre's Levi-Civita connection); $A_X$ is $h$-skew.

## K.1 Results (layer c)

**Lemma K1 (Ricci at a zero of the metric difference).** Let $h$ and $\hat h=h+\delta$ be pseudo-Riemannian metrics on a manifold, $\nabla=\nabla^h$, and let $\delta$ vanish at the point $p$. Put $T_{dbc}:=\nabla_b\delta_{dc}+\nabla_c\delta_{db}-\nabla_d\delta_{bc}$. Then at $p$, with indices raised by $h_p$,
$$\begin{aligned}
\mathrm{Ric}(\hat h)_{bc}-\mathrm{Ric}(h)_{bc}
={}&\tfrac12h^{ad}\big(\nabla_a\nabla_b\delta_{dc}+\nabla_a\nabla_c\delta_{bd}-\nabla_a\nabla_d\delta_{bc}-\nabla_b\nabla_c\delta_{ad}\big)&&\text{(terms 1–4)}\\
&-\tfrac12h^{ap}h^{dq}\,\nabla_a\delta_{pq}\,T_{dbc}+\tfrac12h^{ap}h^{dq}\,\nabla_b\delta_{pq}\,\nabla_c\delta_{ad}&&\text{(terms 5, 6)}\\
&+\tfrac14h^{ap}h^{dq}\,T_{pad}\,T_{qbc}-\tfrac14h^{ap}h^{dq}\,T_{pbd}\,T_{qac}&&\text{(terms 7, 8)}.
\end{aligned}$$
At $p$, $\nabla^2\delta$ is symmetric in its two derivative slots, so the order of the derivatives in terms 1–4 does not matter.

*Proof.*
- $C:=\nabla^{\hat h}-\nabla$ is a tensor. Subtracting the Koszul formulas of the two connections gives, everywhere, $C^a{}_{bc}=\tfrac12\hat h^{ad}T_{dbc}$.
- For a connection $\nabla+C$ one has $\hat R(X,Y)Z=R(X,Y)Z+(\nabla_XC)(Y,Z)-(\nabla_YC)(X,Z)+C(X,C(Y,Z))-C(Y,C(X,Z))$ ($\nabla$ torsion-free). Contracting with $\mathrm{Ric}_{bd}=R^a{}_{bad}$ gives
$$\mathrm{Ric}(\hat h)_{bc}=\mathrm{Ric}(h)_{bc}+\nabla_aC^a{}_{bc}-\nabla_bC^a{}_{ac}+C^a{}_{ad}C^d{}_{bc}-C^a{}_{bd}C^d{}_{ac}.$$
- $C^a{}_{ac}=\tfrac12\hat h^{ad}\nabla_c\delta_{ad}$: the other two terms of $T_{dac}$ cancel against each other under the symmetric $\hat h^{ad}$.
- At $p$, $\hat h^{-1}=h^{-1}$ and $\nabla\hat h^{-1}=-h^{-1}(\nabla\delta)h^{-1}$. Differentiating $C$ by the product rule and substituting gives the formula. Terms 5 and 6 are the $\nabla\hat h^{-1}$ parts of $\nabla_aC^a{}_{bc}$ and $-\nabla_bC^a{}_{ac}$; terms 7 and 8 are the $CC$ terms.
- Symmetry: $\nabla^2_{a,b}\delta-\nabla^2_{b,a}\delta$ is $R(\partial_a,\partial_b)$ acting on δ, which vanishes where δ does. ∎

Terms 1–4 are the classical first variation of the Ricci tensor (Besse, *Einstein Manifolds*, 1.174). The formula is exact at $p$, not an expansion: at a zero of δ, $\mathrm{Ric}(\hat h)_p$ depends only on $h$ near $p$ and on $(\nabla\delta)_p,(\nabla^2\delta)_p$.

**Lemma K2 (the jets of δ along $s(X)$).** Let $\nabla=\nabla^{\rm Z2}$ and $p\in s(X)$. With respect to $T_pY=H\oplus V$, the only non-zero components of $(\nabla\delta)_p$ and $(\nabla^2\delta)_p$ are
$$(\nabla_U\delta)(X,Y)=c_H\,U(X,Y),\qquad (\nabla^2_{U,W}\delta)(X,Y)=\tfrac{c_H}2\,\big(Ug^{-1}W+Wg^{-1}U\big)(X,Y).$$
In particular all components with a horizontal derivative slot or a vertical δ-slot vanish: $\nabla_X\delta$, $(\nabla_U\delta)(U',\cdot)$, $\nabla^2_{X,Y}\delta$, $\nabla^2_{X,U}\delta$ and $\nabla^2_{U,X}\delta$. This holds for every base metric $g$, not only in the flat model.

*Proof.* Let $E:=\pi^*\mathrm{Sym}^2T^*X$ carry the pull-back $\nabla^E$ of the Levi-Civita connection of $g$. Let $\tau\in\Gamma(E)$ be the tautological section $\tau(y)=y$, and put $\varphi:=\tau-\pi^*g$. Let $\iota(k)(Z,W):=k(\pi_*Z,\pi_*W)$. Then $\delta=c_H\,\iota(\varphi)$, and $\varphi$ vanishes exactly on $s(X)$.

1. *$\nabla^E\varphi=K$.* $\nabla^E\pi^*g=\pi^*(\nabla g)=0$. In coordinates, $\nabla^E_{\partial_\mu}\tau=-(\Gamma_\mu^Ty+y\Gamma_\mu)=-\Psi_\mu$ and $\nabla^E_{\partial_A}\tau=E_A$, so $\nabla^E\tau=dy-\Psi\,dx=K$.
2. *$\nabla\pi_*=0$ at $p$.* View $\pi_*$ as a section of $T^*Y\otimes\pi^*TX$ with the connection induced by $\nabla$ and $\nabla^{\pi^*}$. By O'Neill's identities, for basic $X,Y$ and vertical $U,W$:
   - $(\nabla_X\pi_*)Y=\nabla^g_XY-\pi_*\mathcal H\nabla_XY=0$;
   - $(\nabla_U\pi_*)Y=-\pi_*\mathcal H\nabla_UY=-\pi_*A_YU$, since $\pi_*Y=Y\circ\pi$ is constant along the fibre, so $\nabla^{\pi^*}_U(\pi_*Y)=0$;
   - $(\nabla_X\pi_*)U=-\pi_*A_XU$;
   - $(\nabla_U\pi_*)W=-\pi_*T_UW$.

   $T\equiv0$ (F3-i). $A_XY=0$ at $p$ (F3-ii), and $h(A_XU,Y)=-h(U,A_XY)=0$ gives $A_XU=0$ at $p$.
3. *$\nabla K$ at $p$.* Here $(\nabla_ZK)(W)=\nabla^E_Z(K(W))-K(\nabla_ZW)$. For basic $X,Y$, vertical $U$ and the vertical lift $k^V$ of a section $k$ of $\mathrm{Sym}^2T^*X$ (these lifts span $V$):
   - $(\nabla_XK)(Y)=-K(\mathcal V\nabla_XY)=-A_XY$, which vanishes at $p$.
   - $(\nabla_UK)(Y)=-K(T_UY)=0$.
   - $(\nabla_XK)(k^V)=\pi^*(\nabla_Xk)-K(\nabla_Xk^V)=0$. Indeed $\nabla_Xk^V=\nabla_{k^V}X+[X,k^V]$, $\mathcal V\nabla_{k^V}X=T_{k^V}X=0$, and $[X^H,k^V]=(\nabla_Xk)^V$. For the last identity: in coordinates $X^H=X^\mu(\partial_\mu+\Psi_{A\mu}\partial_A)$ and $\Psi$ is linear in $y$, so $[X^H,k^V]_A=X(k_A)-X^\mu\Psi_{A\mu}(k)=(\nabla_Xk)_A$.
   - $(\nabla_UK)(k^V)=-K(\hat\nabla_Uk^V)=-\Gamma^F_y(U,k)$. On the fibre, an open subset of the vector space $\mathrm{Sym}^2T^*_xX$, $k^V$ is a constant field, and $\Gamma^F$ is the Christoffel symbol of $c_VG_\lambda$ in linear coordinates. Koszul's formula for constant fields gives
   $$\Gamma^F_y(U,W)=-\tfrac12\big(Uy^{-1}W+Wy^{-1}U\big)\qquad\text{for every }\lambda:$$
   with $u=y^{-1}U$ etc., $2G_\lambda(\Gamma^F(U,W),Z)=-\mathrm{tr}(uwz)-\mathrm{tr}(wuz)+2\lambda\,\mathrm{tr}(uw)\,\mathrm{tr}\,z=UG_\lambda(W,Z)+WG_\lambda(U,Z)-ZG_\lambda(U,W)$.

   So at $p$ the only component of $\nabla K$ is $(\nabla_UK)(W)=\tfrac12(Ug^{-1}W+Wg^{-1}U)$.
4. *Assembly.* By the Leibniz rule, $\nabla_Z(\iota\varphi)=\iota(\nabla^E_Z\varphi)+\iota'_Z(\varphi)$ with $\iota'_Z(\varphi)(W_1,W_2)=\varphi((\nabla_Z\pi_*)W_1,\pi_*W_2)+\varphi(\pi_*W_1,(\nabla_Z\pi_*)W_2)$. Differentiating once more, $\nabla^2(\iota\varphi)=\iota(\nabla^2\varphi)$ plus terms containing a factor $\nabla\pi_*$ or a factor $\varphi$. At $p$ both factors vanish (step 2, and $\varphi_p=0$). Hence $(\nabla\delta)_p=c_H\,\iota(K)$ and $(\nabla^2\delta)_p=c_H\,\iota((\nabla K)_p)$, with $\nabla^2\varphi=\nabla K$ by step 1. Since $\iota$ kills vertical arguments and $K$ kills horizontal ones, steps 1 and 3 give the statement. ∎

Off $s(X)$ the argument breaks at step 2: $\nabla\pi_*\ne0$ and $(\nabla_XK)(Y)=-A_XY\neq0$ there. Besides $\delta_p=0$ (K1) and $\varphi_p=0$ (step 4), this is the only place where the restriction to the section enters. (Corrected after review 2026-10-30, R7-4.)

**Lemma K3 (the algebra, for every $n$).** Let $(H,g)$ be a nondegenerate quadratic space of dimension $n\ge1$ and $V=\mathrm{Sym}^2H^*$. Put $h=c_Hg\oplus c_VG_\lambda|_g$ on $H\oplus V$ ($\lambda\ne1/n$, $c_H,c_V\neq0$), and let $\nabla\delta,\nabla^2\delta$ be the tensors of Lemma K2. Then the right-hand side of Lemma K1 equals
$$\Delta\mathrm{Ric}_{HH}=\frac{c_H}{4c_V(n\lambda-1)}\,g,\qquad \Delta\mathrm{Ric}_{HV}=0,\qquad \Delta\mathrm{Ric}_{VV}=-\tfrac14G_0 .$$

*Proof.* Write $E_A$ for a basis of $V$ and $G_\lambda^{AB}$ for the inverse of $G_\lambda(E_A,E_B)$. Four facts are used.
- (i) Sherman–Morrison, since $G_\lambda=G_0-\lambda\,\mathrm{tr}\otimes\mathrm{tr}$, $\mathrm{tr}=G_0(g,\cdot)$ and $G_0(g,g)=n$:
$$\textstyle\sum_{AB}G_\lambda^{AB}(E_A)_{ij}(E_B)_{kl}=\tfrac12(g_{ik}g_{jl}+g_{il}g_{jk})+\mu\,g_{ij}g_{kl}.$$
- (ii) Contracting (i): $\sum G_\lambda^{AB}\,E_Ag^{-1}E_B=\big(\tfrac{n+1}2+\mu\big)g$.
- (iii) Contracting (i) the other way: $\sum G_\lambda^{AB}\,\mathrm{tr}(g^{-1}E_A)\,E_B=(1+n\mu)\,g=g/(1-n\lambda)$.
- (iv) By Lemma K2, $T$ has exactly one vertical slot, with the components $T_{\mu\nu A}=T_{\mu A\nu}=c_H(E_A)_{\mu\nu}$ and $T_{A\mu\nu}=-c_H(E_A)_{\mu\nu}$. $\nabla\delta$ is non-zero only with a vertical derivative slot and horizontal δ-slots, and $\nabla^2\delta$ only with vertical derivative slots and horizontal δ-slots. Since $h^{-1}$ is block-diagonal, each contraction joins slots of the same type.

Going through the eight terms of Lemma K1:

| term | $\Delta\mathrm{Ric}_{HH}$ (coefficient of $g$) | $\Delta\mathrm{Ric}_{VV}$ | reason |
|---|---|---|---|
| 1, 2 | 0 | 0 | $h^{ad}$ pairs a derivative slot of $\nabla^2\delta$ (non-zero only if vertical) with a δ-slot (non-zero only if horizontal) |
| 3: $-\tfrac12h^{ad}\nabla_a\nabla_d\delta_{bc}$ | $-\tfrac{c_H}{2c_V}\big(\tfrac{n+1}2+\mu\big)$ | 0 | (ii) |
| 4: $-\tfrac12h^{ad}\nabla_b\nabla_c\delta_{ad}$ | 0 | $-\tfrac12G_0$ | $\tfrac1{c_H}g^{\rho\sigma}\tfrac{c_H}2(Ug^{-1}W+Wg^{-1}U)_{\rho\sigma}$ |
| 5 | 0 | 0 | $h^{ap}\nabla_a\delta_{pq}$ pairs the derivative slot (vertical) with a δ-slot (horizontal) |
| 6: $\tfrac12h^{ap}h^{dq}\nabla_b\delta_{pq}\nabla_c\delta_{ad}$ | 0 | $+\tfrac12G_0$ | $\tfrac12\mathrm{tr}(g^{-1}Ug^{-1}W)$ |
| 7: $\tfrac14h^{ap}h^{dq}T_{pad}T_{qbc}$ | $-\tfrac{c_H}{4c_V(1-n\lambda)}$ | 0 | (iii); in VV, $T_{qbc}$ would have two vertical slots |
| 8: $-\tfrac14h^{ap}h^{dq}T_{pbd}T_{qac}$ | $+\tfrac{c_H}{2c_V}\big(\tfrac{n+1}2+\mu\big)$ | $-\tfrac14G_0$ | HH: two equal contributions by (ii); VV: $-\tfrac14\mathrm{tr}(g^{-1}Ug^{-1}W)$ |

The HH entries sum to $-c_H/(4c_V(1-n\lambda))$ and the VV entries to $-\tfrac14G_0$.

HV block: $-1\in O(g)$ acts on $H\oplus V$ by $(-1,+1)$ and preserves $h$, $\nabla\delta$ and $\nabla^2\delta$. Since $\Delta\mathrm{Ric}$ is a contraction of these, it is invariant, so $\Delta\mathrm{Ric}(X,U)=\Delta\mathrm{Ric}(-X,U)=0$. Equivalently, $P:=\mathrm{id}_H\oplus(-\mathrm{id}_V)$ preserves $h$ and $\nabla^2\delta$ and reverses $\nabla\delta$. Every term is linear in $\nabla^2\delta$ or quadratic in $\nabla\delta$, so $\Delta\mathrm{Ric}$ is $P$-invariant, and its HV block vanishes. ∎

This is the precise form of the sketch's "∇δ is odd under Φ". At $p$, $d\Phi_p=P$, since Φ fixes $s(X)$ pointwise and is the geodesic symmetry of the fibre. But Φ is an isometry of Z2, not of Z1, so it cannot be applied to $\mathrm{Ric}(\mathrm{Z1})$ directly. What is used is that $P$ acts on the data $(h_p,\nabla\delta_p,\nabla^2\delta_p)$ by $(+,-,+)$.

*Computation (K3).* The script evaluates the eight terms of Lemma K1 on the data of Lemma K2 by the $O(g)$-invariant strand-diagram (Brauer) calculus. An $H$-slot has one endpoint and a $V$-slot two; every strand is $g^{\pm1}$; a closed loop counts $\mathrm{tr}\,\mathrm{id}=n$, with $n$ a symbolic loop weight. Coefficients lie in $\mathbb Q(n,\lambda,c_H,c_V)$. Only the evaluation map from diagrams to tensors is used, which is valid for every $n\ge1$; linear independence of diagrams is never assumed. The output reproduces every entry of the table and the three block values. An HV form has three endpoints, so no diagram exists. The same engine at loop weights $1,2,3,4$ agrees with a direct evaluation on explicit $(n+N)$-dimensional tensors.

**Theorem K (F4(c) for every $n$).** For every $n\ge1$, every signature of $g$, every $\lambda\ne1/n$ and every $c_H,c_V\neq0$, the Ricci tensor of Z1 along $s(X)$ is as stated in K.0:
$$\mathrm{Ric}^{\rm Z1}_{HH}=\mathrm{Ric}_g+\frac{c_H}{4c_V(n\lambda-1)}\,g,\qquad \mathrm{Ric}^{\rm Z1}_{HV}=0,\qquad \mathrm{Ric}^{\rm Z1}_{VV}=-\tfrac{n+1}4\,G_{1/(n+1)} .$$

*Proof.* Let $p\in s(X)$. $\delta_p=0$, so Lemma K1 applies with $h=$ Z2, $\hat h=$ Z1. Lemma K2 gives $(\nabla\delta)_p,(\nabla^2\delta)_p$ in terms of $g_x$ alone, and Lemma K3 evaluates the formula. Adding $\mathrm{Ric}(\mathrm{Z2})_p$ from F3-iii gives the blocks. For the VV block, $-\tfrac n4G_{1/n}-\tfrac14G_0=-\tfrac{n+1}4G_0+\tfrac14\mathrm{tr}\otimes\mathrm{tr}=-\tfrac{n+1}4G_{1/(n+1)}$. ∎

**Corollary K4 (now proven for every $n$).**
- (a) $\Lambda_P$ for Z1 in Cor. F5: $\Lambda_P=\Lambda_A+c_H/(4c_V(n\lambda-1))$, from the HH block and F4(b).
- (b) The Z1 half of Theorem F6 ($n\ge2$). On $s(X)$, $\mathrm{Ric}(\mathrm{Z1})=\kappa\,\mathrm{Z1}$ iff $\lambda=1/(n+1)$ and $\mathrm{Ric}(g)=0$, with $\kappa=-(n+1)/(4c_V)$. The VV block is proportional to $c_VG_\lambda$ iff the $G_0$ and $\mathrm{tr}\otimes\mathrm{tr}$ coefficients match (these two forms are linearly independent for $n\ge2$; for $n=1$ they coincide), which forces $\kappa c_V=-(n+1)/4$ and $\lambda=1/(n+1)$. Then the HH block is $\kappa c_Hg+\mathrm{Ric}_g$. Output K3 solves this symbolically.
- (c) Consistency with F4(b), which F proves independently by Gauss and Walczak: the $h^{-1}$-trace of the blocks is $R/c_H+\sigma_1/c_V$. Output K3 checks this in $\mathbb Q(n,\lambda,c_H,c_V)$; the fibre part contributes $s_F=-n(n-1)(n+2)/(8c_V)$.
- (d) The vertical part of the Christoffel difference is $C(X,Y)^V=-\tfrac{c_H}{2c_V}[X^\flat\odot Y^\flat+\mu\,g(X,Y)\,g]$, which is II of F4(a), since $s(X)$ is totally geodesic for Z2 (output K3).

*Remark (why the flat model of the sketch gives the general answer).* By Lemma K2, the data $(h_p,\nabla\delta_p,\nabla^2\delta_p)$ in the adapted frame are the same functions of $g_x$ for every base metric. The curvature of $g$ could enter $(\nabla^2\delta)_p$ only through $A$ (steps 2 and 3), and $A$ vanishes on $s(X)$. In the flat model (base $\mathbb R^n$ with the constant metric $g_x$) Z2 is the product $c_Hg_x\,dx^2+c_VG_\lambda(dy,dy)$, so Lemma K2 holds there trivially and $\Delta\mathrm{Ric}$ is the same. The flat model has $\mathrm{Ric}(\mathrm{Z2})=0\oplus0\oplus(-\tfrac n4G_{1/n})$, so its $\mathrm{Ric}(\mathrm{Z1})-\mathrm{Ric}(\mathrm{Z2})$ is exactly $\Delta\mathrm{Ric}$. That is what the sketch evaluated with O'Neill for the projection $(x,y)\mapsto y$. Both conditions the sketch left implicit are used: $A=0$ on $s(X)$, and the fact that the curvature of the curved base enters the jets of δ at $p$ through $A$ alone. The sketch's flat-model formulas hold at every $y$, by the $GL(n)$-homogeneity of the flat model; K5 checks them at a point $y_0\ne\eta$.

## K.2 What is verified, precisely

All checks are exact and asserted. Signatures are written (positive, negative).

| Statement | How (section of `K_z1_ricci.out`) | Scope |
|---|---|---|
| curvature routine | K0: $4|dx|^2/(1\pm|x|^2)^2$ has $\mathrm{Ric}=\pm(m-1)g$; K4: two independent routes to $\mathrm{Ric}_g(0)$ agree | $m=2,3,4$; every K4 base |
| Lemma K1 | K1: random exact 2-jets of $h$ and δ with $\delta(0)=0$; both Ricci tensors computed directly and compared with the formula. The quadratic terms are non-zero in every case | $m=2,3,4,5$; signatures (2,0), (1,1), (3,0), (1,2), (2,2), (3,2) |
| Lemma K2 | K2: Z2 on random curved bases (3-jets, non-diagonal $g_0$, $\Gamma_g(0)\ne0$, $\mathrm{Ric}_g(0)\ne0$), coefficients in $\mathbb Q(\lambda,c_H,c_V)$. Checked: (a) $s_*$ is the horizontal lift; (b) $\nabla\pi_*=0$ at $p$; (c) $\nabla K$ at $p$; (d) all components of $\nabla\delta,\nabla^2\delta$ at $p$ in the adapted frame; (e) fibre Christoffel symbols, λ-free | $n=2,3,4$; signatures (1,1), (2,1), (1,3); every λ |
| Lemma K3 | K3: strand-diagram evaluation of all eight terms with symbolic loop weight $n$, in $\mathbb Q(n,\lambda,c_H,c_V)$; term-by-term table printed; engine cross-checked against explicit tensors | all $n$ (symbolic); cross-check at $n=1,2,3,4$ |
| Theorem K | proof. Independently, K4: brute-force Ricci tensor of Z1 at $p\in s(X)$ from the exact coordinate 2-jet of the metric on $Y$ (dimension 5 to 27); random curved base 3-jets in non-normal coordinates, non-diagonal $g_0$, $\nabla\mathrm{Ric}_g(0)\ne0$, $\mathrm{Ric}_g(0)$ not Einstein for $n\ge3$ | $n=2,3$: symbolic in $(\lambda,c_H,c_V)$, signatures (2,0), (1,1), (3,0), (2,1). $n=4$: (1,3), (2,2), (4,0). $n=5$: (5,0), (3,2). $n=6$: (5,1), (3,3). For $n\ge4$, two parameter sets per base from $(\lambda,c_H,c_V)\in\{(\tfrac12,1,1),(\tfrac27,-\tfrac53,\tfrac34),(1,\tfrac32,-\tfrac23),(-\tfrac35,2,\tfrac52),(\tfrac15,\tfrac74,-\tfrac13),(\tfrac38,-\tfrac12,-\tfrac45)\}$ |
| F3 (input) | K4: the same brute force for Z2 | same cases |
| chain K1+K2 | K4: $\mathrm{Ric}(\mathrm{Z1})-\mathrm{Ric}(\mathrm{Z2})$ (brute force) equals the K1 formula evaluated on the brute-force $\nabla^{\rm Z2}\delta$, $\nabla^{\rm Z2}\nabla^{\rm Z2}\delta$ in coordinates | every K4 case |
| Cor. K4, F4(a)(b) consistency | K3: diagram calculus and sympy | all $n$ |
| flat model of the sketch | K5: brute force at $y_0\ne\eta$ | $n=3$ (symbolic), $n=4$ |

For $n=6$ these are the first exact checks of F4(c) in the repository, and for $n=2,\dots,5$ they are independent of R2/R3. They agree with F.

## K.3 Limits (layers b, c)

- **Inputs.** The proof uses F3 (Z2) and O'Neill's identities (L-F1), both proven in F. K4 re-verifies F3's blocks independently for $n\le6$, but the general-$n$ statement rests on F3's proof.
- **Scope: on the section.** The argument gives $\mathrm{Ric}(\mathrm{Z1})$ only at points of $s(X)$. Off $s(X)$, $\delta\ne0$ and $\nabla\pi_*\neq0$, and $\nabla^2_{X,Y}\delta$ contains, among other terms, $-c_H\,\iota(A_XY)$. Nothing here bears on the off-section behaviour (F6, E2(v), Conjecture F9), nor on derivatives of $\mathrm{Ric}$ along $s(X)$.
- **Nondegeneracy.** $\lambda\neq1/n$ and $c_Hc_V\neq0$ are needed for $h$ to be a metric. No other restriction on λ, $c_H$, $c_V$ or the signature is used.
- **Computational checks.** These are exact but cover finitely many bases and, for $n\ge4$, finitely many parameter values. The statement for all $n$ and all parameters rests on the proof; the symbolic-$n$ evaluation of K3 covers only its algebraic step.
- **Reconstruction.** The reading R-Z of F (Levi-Civita splitting, DeWitt fibre metric, Z1 as the draft-literal reading) is assumed as in F. This note changes no layer-(b) choice and no conclusion about the draft's intent.

## K.4 Consequences for gaps and claims

- **Theorem F4(c):** proven for every $n$, in every signature (previously: proven for $n=1,2$, exact verification for $n=3,4,5$, general-$n$ sketch). The closed form is unchanged; no discrepancy with F was found.
- **Downstream.** Every statement whose Z1 part "inherits the status of F4(c)" now has that part proven, under its own other hypotheses: $\Lambda_P$ (Z1) of Cor. F5; the Z1 half of Theorem F6; Cor. I2.1 and I7 in Z1; J3–J4 in Z1.
- **GU-054:** status C unchanged. Its Z1 statements (route dependence $\Lambda_P\ne\Lambda_A$; Einstein along $s(X)$ only for $\lambda=1/(n+1)$ and $\mathrm{Ric}(g)=0$) are now proven for every $n$, not only verified for $n\le5$.
- **GU-G14:** content unchanged.
