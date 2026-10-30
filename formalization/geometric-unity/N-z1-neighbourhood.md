# N. Υ = 0 on a neighbourhood of a point of $s(X)$ in the Zorro metric Z1

Claims: GU-054 (context). Gaps: GU-G14 (and F9). Code: `computation/geometric-unity/N_z1_neighbourhood.py` → `.out` (sections N0–N9; exact: Python Fractions, also as coefficients of truncated power series, and sympy rationals; fixed seeds; 450 s (7.5 min) on one core; a second run gave identical output apart from the timing annotations, 473 s). It imports `F_jets.py` and `L_nonradial.py` (and through it `J_slaved_sector.py`, …), all unmodified. Output sections: §N0 the line equation, §N1 Lemma N1, §N2 Lemma N2, §N3 Lemma N5, §N4 Lemma N6, §N5 Theorem N7, §N6 Proposition N8, §N7 Propositions N9–N10, §N8 Proposition N11, §N9 Theorem N12(c). Record: `metrics/geometric-unity/N-z1-neighbourhood.yaml`. Builds on F (F1, F4, F6, F7, R-Z), K (K1, K2, Theorem K), J (J1, J2, J4), J5d and L (L1, L3, L5, L8). Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## N.0 Question and setting (layer b)

J.3 and L.3 list as open $\Upsilon=0$ on a neighbourhood in the draft-literal metric Z1, "where the Euler field is not Killing, so $s$ is not parallel; L6 and the radial structure are unavailable". J4 shows that (E) at the points of $s(X)$ is solvable in Z1 with free second normal derivatives. This note asks at which order of the normal jets the neighbourhood equations constrain the base, and answers the question for the members treated by L8.

**Hypotheses (H1)** (those of L8, with Z1 for Z2):
- Z1 with $n=4$, $\lambda\neq\tfrac14$, $c_Hc_V\neq0$; the real forms of J.0: (14,0) (base (4,0), $\lambda<\tfrac14$) and (7,7) (base (1,3), (3,1) or (2,2), $\lambda>\tfrac14$; A4, A5).
- $\kappa\neq0$; $T_o=t\Phi^1+u\Phi^1\Gamma$ with $t,u$ smooth on $U$; $T_e=\Sigma\,(dt,du)^\top$ the slaved field of J1 in $\mathcal V$, $\det M(t,u)\neq0$ on $U$.
- $U\subset Y$ connected and open, $\Upsilon=0$ on $U$, $U\cap s(X)\neq\emptyset$.

**Notation.** As in L.0: $E^{(i)}$, $H$, $\mathcal Q,\mathcal Q'$, $\varpi$, $x_1=x$, $x_{13}=y$, $N=x^2\pm y^2$, $\ell=xz\pm yw$; the classes (A), (B) of L8. Further:
- $E=y\,\partial_y$, $s=\tfrac14\log|\det(g^{-1}y)|$, $P_H$ the projection onto $H$ along $V$, $V_0=\{U\in V:\mathrm{tr}(g^{-1}U)=0\}$; for $X\in H$, $U\in V$ write $\tilde U=g^{-1}U$.
- $c:=4c_V(1-4\lambda)$, $\ b:=\dfrac{1-5\lambda}{c_V(1-4\lambda)}$, $\ \beta_S:=-\dfrac5{4c_V}$, $\ \mu_1:=-\tfrac7{18}$, $\ \sigma_g:=\mathrm{scal}_g/c_H$.
- $Z:=\{q\in U:dt=du=0\}$, $U^*:=U\setminus Z$.
- $\Sigma:=\{(x,e^\theta g(x))\}\subset Y$, the radial set. Since $e^\theta g$ has the Levi-Civita connection of $g$, and Z1 depends on $g$ only through ℵ (F2.4), each point of Σ lies on the section of the parallel metric $e^\theta g$. Every statement proven at the points of $s(X)$ therefore holds at the points of Σ with $g$ replaced by $e^\theta g$ (and $\mathrm{Ric}_g$ unchanged, $\mathrm{scal}_g$ multiplied by $e^{-\theta}$).

**Real forms.** As in L.0, all algebraic statements are polynomial identities, valid in both real forms.

**Inputs (proven earlier).** F4 (Z1 on $s(X)$: II, Ricci blocks, scalar curvature $\mathrm{scal}_g/c_H+\sigma_1/c_V$), F6, F7, Theorem K and K1; J1–J2 and the tables; L1 and L5 (metric-free), L3 (fibre lemma), the steps of L8.

## N.1 Results (layer c)

**Lemma N1 (the Euler field of Z1).** *[proven]* At every point of $Y$, for every base metric:

(a) The fibres are totally geodesic ($T\equiv0$). The symmetric part of $\mathcal V\nabla_XY$ ($X,Y$ basic) is $-\tfrac{c_H}{2c_V}\big[X^\flat\odot Y^\flat+\mu\,y(X,Y)\,y\big]$ with $\flat$ taken with $y$; its horizontal trace is $-E/(2c_V(1-4\lambda))$.

(b) $E^\flat=c\,ds$, $|E|^2=c$, and $\nabla E=\tfrac12P_H$. Equivalently $\mathcal L_Eg_\aleph=g_\aleph|_H$ and $\nabla^2s=g_\aleph|_H/(2c)$. The fibre scaling $\phi_\theta(y)=e^\theta y$ satisfies $\phi_\theta^*\mathrm{Z1}(c_H)=\mathrm{Z1}(e^\theta c_H)$.

(c) $E$ is neither Killing nor homothetic, and Z1 has no nonzero conformal Killing field tangent to the fibres.

(d) $\mathrm{Ric}(E,\cdot)=-\tfrac1c\,E^\flat=-ds$.

(e) $R(X,Y)E=\tfrac12[(\nabla_XP_H)Y-(\nabla_YP_H)X]$. At $p\in s(X)$ (and at the points of Σ), for $X,Y\in H$ and $U,W\in V$: $R(X,Y)E=0$, $R(X,U)E=-\tfrac14\tilde UX$, $R(U,W)E=0$, $R(E,X)E=\tfrac14X$.

*Proof.*
- (a) For vertical $U,W$ and basic $X$, Koszul gives $2g(\nabla_UW,X)=-Xg(U,W)+g([U,W],X)-g([U,X],W)-g([W,X],U)$, which involves only the vertical metric and brackets; these are those of Z2, where $T=0$ (F3). For basic $X,Y$ and vertical $U$, $2g(\nabla_XY,U)=g([X,Y],U)-U\,g(X,Y)$ and $U\,g(X,Y)=c_HU(X,Y)$; Sherman–Morrison inverts $G_\lambda$ as in F4(a), now at the point $y$.
- (b) $E^\flat=c_VG_\lambda(y,Dy)=c_V(1-4\lambda)\,\mathrm{tr}(y^{-1}Dy)=4c_V(1-4\lambda)\,ds$, as in J5(a). $\phi_\theta^*$ multiplies $c_Hy(dx,dx)$ by $e^\theta$ and leaves $G_\lambda(Dy,Dy)$ unchanged, so $\mathcal L_Eg_\aleph=g_\aleph|_H$. As $E^\flat$ is closed, $\nabla E^\flat=\tfrac12\mathcal L_Eg_\aleph$.
- (c) For vertical $K$ and basic $X,Y$, $(\mathcal L_Kg_\aleph)(X,Y)=K\,g_\aleph(X,Y)=c_HK(X,Y)$ since $[K,X]$ is vertical. If $\mathcal L_Kg_\aleph=fg_\aleph$, then $K=fy=fE$ on each fibre. On $V$, $\mathcal L_{fE}g_\aleph=f\mathcal L_Eg_\aleph+2\,df\odot E^\flat=2\,df\odot E^\flat$ has rank $\le2$, against $fg_\aleph|_V$ of rank 10, unless $f=0$.
- (d) For a gradient field with $\nabla E=B$ ($B$ symmetric), $\mathrm{Ric}(E,\cdot)=\mathrm{div}B-d\,\mathrm{tr}B$. Here $B=\tfrac12P_H$, $\mathrm{tr}B=2$, and $\mathrm{div}\,P_H$ is the horizontal mean curvature vector (the fibres are totally geodesic), $-E/(2c_V(1-4\lambda))$ by (a). So $\mathrm{Ric}(E,\cdot)=-E^\flat/(4c_V(1-4\lambda))=-E^\flat/c$.
- (e) $R(X,Y)E=(\nabla_XB)Y-(\nabla_YB)X$ with $B=\tfrac12P_H$. Use $(\nabla_ZP_H)W=(\nabla_ZW)^V$ for $W$ horizontal and $-(\nabla_ZW)^H$ for $W$ vertical, $T=0$, $[X^H,Y^H]^V=-(R_g(X,Y)\cdot y)^V=0$ at $y=g$ (L-F5), and $(\nabla_XU)^H=\tfrac12\tilde UX$ on $s(X)$ (from $g(\nabla_XU,Y)=-g(U,\mathrm{II}(X,Y))=\tfrac{c_H}2U(X,Y)$). ∎

Output §N1 checks (a)–(d) at points **off** the section ($y_0\neq g(x_0)$, $\Gamma(x_0)\neq0$, cubic random bases) in base signatures (4,0), (1,3), (2,2); §N4 checks (e) at $p$.

So what replaces the parallel Killing field of Z2 is a gradient field with $\nabla E=\tfrac12P_H$ that is still a Ricci eigenvector everywhere, now with eigenvalue $-1/c\neq0$ (in Z2: 0).

**Lemma N2 (flat model; F9 for flat ℵ).** *[proven]* If ℵ is flat near $\pi(q)$, then near $q$
$$\mathrm{Ric}(\mathrm{Z1})=\beta_S\,g_\aleph+b\,\big(g_\aleph|_H+c\,ds\otimes ds\big)=\beta_S\,g_\aleph+4bc\,e^{-s/2}\nabla^2e^{s/2}.$$
Hence Z1 over a flat base is Einstein on some open set iff it is Einstein everywhere iff $\lambda=\tfrac15$; then $\mathrm{Ric}=-\tfrac54c_V^{-1}g_\aleph$ and $G=\tfrac{15}{2c_V}g_\aleph$ (F7).

*Proof.* The flat model is homogeneous under $GL(4)\ltimes\mathbb R^4$, so F4(c) at $y=g$ transports to every point: $\mathrm{Ric}_{HH}=\tfrac{c_H}{4c_V(4\lambda-1)}y=-\tfrac1cg_\aleph|_H$ and $\mathrm{Ric}_{VV}=-\tfrac54G_{1/5}=\beta_Sg_\aleph|_V+4(1-5\lambda)\,ds^2$, with $4(1-5\lambda)=bc$ and $\beta_S+b=-1/c$. The second form uses N1(b). ∎ (Output §N2 at $y_0\neq\eta$, five values of λ.)

**Lemma N3 (what transfers from L).** *[proven]*

(a) *Metric-free steps.* L1 ($dt\wedge du=0$ on open solution sets for $\circledast_\sharp$, $xW_A$, $xW_A+zW_C\mathrm{Ad}_\Gamma$ and every member with $\ell\neq0$), L5 (the curvature-free channel, $R'$, the line and its invariance), and the line reduction of L8 step 5 are identities in $(t,u)$ and the J1–J2 tables. They do not involve the metric and hold verbatim in Z1.

(b) *Concircular functions on Z1.* Let $W\subset Y$ be open and $r$ a function on $W$ with $\nabla^2r=\psi\,g_\aleph$ and $\mathrm{Ric}\in\mathrm{span}(g_\aleph,dr\otimes dr)$ pointwise. Then $dr=0$.

*Proof of (b).* The fibres are totally geodesic (N1(a)) with the metric $c_VG_\lambda$ of Z2, so L3(a) gives $\psi=0$: $Z:=\nabla r$ is parallel. Write $\mathrm{Ric}=\rho_1g+\rho_2\,dr\otimes dr$. Then $\mathrm{Ric}(Z)=0$ gives $\rho_1+\rho_2|Z|^2=0$, and N1(d) gives $\rho_1E+\rho_2\langle Z,E\rangle Z=-E/c$.
- If $Z,E$ are independent on an open set: $\rho_1=-1/c$, so $\rho_2|Z|^2=1/c\neq0$ and $\langle Z,E\rangle=0$. With $f:=\langle Z,E\rangle$, N1(b) gives $\nabla f=\tfrac12Z^H$, so $Z$ is vertical. For a vertical parallel $Z$ and basic $X,Y$, $0=g(\nabla_XZ,Y)_{\rm sym}=\tfrac{c_H}2Z(X,Y)$ by N1(a), so $Z=0$.
- If $Z=\phi E$ on an open set: $0=\nabla Z=d\phi\otimes E+\tfrac\phi2P_H$ forces $\phi=0$. ∎

(c) *Consequences for (H1)* (as in L8 steps 1–5, with (b) for L4).
- Class (A): $U^*=\emptyset$, so $(t,u)$ is constant on $U$.
- Class (B), including $\circledast_\sharp$: on $U^*$, $(t,u)=\sigma\,(x,y)$.
- $xW_A$ ($x\neq0$): on $U^*$, $u=0$. (L5(c) gives $u\,dt=0$. On $\{u\neq0\}\cap U^*$, $dt=0$ and $T_e=0$, so $\mathrm{Ric}=-\tfrac{c_1}{6x}g_\aleph$, hence $c_1=6x/c$ by N1(d) and Schur, $c_{13}=0$ gives $t=-\kappa/(208x)$, so $u^2$ and then $u$ are locally constant, against $(dt,du)\neq0$.) The line $u=0$ carries the same coefficients as for $\circledast_\sharp$ (L5, output §L5).

(d) *The line equation.* On the line, with $r:=-1/(224N(\kappa+224N\sigma))$ (smooth and nonzero where $\det M\neq0$), $\Upsilon=0$ is equivalent to
$$(\star)\quad \nabla^2r=\mu_1\,r\,\mathrm{Ric}+\nu\,g_\aleph,\qquad (\mathrm{Tr})\quad \nu=\tfrac7{468}\,r\,s_Y-\tfrac1{13}\Phi(r)+\tfrac3{28}\,\frac{|dr|^2}r ,$$
with $s_Y$ the scalar curvature of Z1 and $\Phi(r)=\tfrac7{18}\,r\,(\kappa\sigma+312N\sigma^2)$. Here $r\Phi(r)$ is a polynomial in $r$ with constant term $A=13/(269746176\,N^3)\neq0$ (output §N0, from L's own `reduced`).

*Proof of (d).* The coefficients of L5 on the line (Hessian $576N/(\kappa+224N\sigma)$, gradient terms $-258048N^2/(\ldots)^2$ and $78336N^2/(\ldots)^2$, zeroth order $\kappa\sigma+312N\sigma^2$) give $\tilde E=-G+\tilde c\,g-\tfrac{18}{7r}(\nabla^2r-\Box r\,g)-\tfrac{351}{98r^2}|dr|^2g$. Its traceless part is $(\star)$ with $\mu_1=-\tfrac7{18}$ (L8 step 5), and its trace, with $\Box r=\mu_1rs_Y+14\nu$, is (Tr). ∎

**Proposition N4 (trichotomy at the points of $s(X)$).** *[proven]* Assume (H1), member $\circledast_\sharp$, $xW_A$, or in class (B). Every $q\in U\cap s(X)$ is of one of three types:
- (i) $q\in\mathrm{int}\,Z$;
- (ii) $q$ has a neighbourhood $N$ with $N\cap\mathrm{int}\,Z=\emptyset$; then $(\star)$, (Tr) hold on $N$;
- (iii) $q\in\overline{\mathrm{int}\,Z}\setminus\mathrm{int}\,Z$.

In cases (i) and (iii), Z1 is Einstein at $q$ with $\mathrm{Ric}=-\tfrac1cg_\aleph$, hence (F6) $\lambda=\tfrac15$ and $\mathrm{Ric}_g(\pi(q))=0$.

*Proof.* On each component of $\mathrm{int}\,Z$, $T_e=0$ and the odd equation is I6's, so Z1 is Einstein there with constant $-1/c$ (N1(d)). By continuity it is Einstein at every point of $\overline{\mathrm{int}\,Z}$. In case (ii), $U^*$ is dense in $N$. N3(c) puts $(t,u)$ on the line on $U^*\cap N$, hence on $N$ by continuity, so $\sigma$ and $r$ are smooth on $N$. The equation $\tilde E=0$ holds on a dense subset, hence on $N$. ∎

**Lemma N5 (fibre jets of Ric(Z1); side result).** *[computed]* Let $p=g(x_0)$ and $h=y-g(x_0)$ along the fibre over $x_0$. In the adapted frame,
$$\mathrm{Ric}(\mathrm{Z1})-\mathrm{Ric}(\mathrm{Z2})=\Big[-\tfrac{c_H}c\,y\Big]_{HH}\oplus\Big[-\tfrac54G_{1/5}+G_{1/4}\Big]_{VV}+O(|h|^3)\ \text{(HH, VV)},\qquad O(|h|^2)\ \text{(HV)} .$$
So, to second order in the fibre, the curvature corrections of Z1 off the section are exactly the O'Neill A-terms of Z2 ($\mathcal A$ in VV, $-2(A_X,A_Y)$ in HH). At order 3 (HH, VV) and order 2 (HV) they differ.

*Status: computed* (output §N3: random $R$ and cubic terms, so $\nabla R\neq0$, in base signatures (4,0), (1,3), (2,2); exact to order 3 for (4,0) and order 2 for the others). No hand proof. N5 is not used below.

**Lemma N6 (data at $p\in s(X)$).** *[(a)–(e) proven; (f) computed]* At $p$, for $X,Y,Z\in H$, $U,W\in V$:

(a) $\mathrm{Ric}=\big(\mathrm{Ric}_g-\tfrac{c_H}cg\big)\oplus(-\tfrac1c\text{ on }E)\oplus(\beta_S\text{ on }V_0)$, $\mathrm{Ric}_{HV}=0$ (F4(c)).

(b) $(\nabla_X\mathrm{Ric})(Y,Z)=(\nabla_X\mathrm{Ric}_g)(Y,Z)$; $(\nabla_X\mathrm{Ric})(Y,U)=-\mathrm{Ric}(\mathrm{II}(X,Y),U)-\mathrm{Ric}(Y,\tfrac12\tilde UX)$; $(\nabla_X\mathrm{Ric})(U,W)=0$.

(c) $(\nabla_E\mathrm{Ric})(X,Y)=-\mathrm{Ric}_g(X,Y)$, and $\nabla_E\mathrm{Ric}$ vanishes on $H\times V$ and $V\times V$.

(d) $(\nabla_\zeta\mathrm{Ric})(E,\cdot)=-\tfrac1c\,g(\nabla_\zeta E,\cdot)-\mathrm{Ric}(\nabla_\zeta E,\cdot)$ for every ζ.

(e) $R(U,W)U'=R^S(U,W)U'$ for $U,W,U'\in V_0$ ($S$ the unimodular slice, $R^S(U,W)U'=-\tfrac14[[\tilde U,\tilde W],\tilde U']$).

(f) For $U\in V_0$: $(\nabla_U\mathrm{Ric})|_{V\times V}=0$ and $(\nabla_U\mathrm{Ric})(X,Z)=-\tfrac12\big[\mathrm{Ric}_g(X,\tilde UZ)+\mathrm{Ric}_g(Z,\tilde UX)\big]$.

*Proof.*
- (a) is F4(c).
- (b) Differentiate the identities $\mathrm{Ric}(Y,Z)=(\mathrm{Ric}_g-\tfrac{c_H}cg)(Y,Z)$, $\mathrm{Ric}(U,Z)=0$ and $\mathrm{Ric}(U,W)=-\tfrac54G_{1/5}(U,W)$ along $s(X)$ (with $U,W$ sections of the normal bundle, parallel for the ℵ-induced normal connection, which preserves $G_{1/5}|_g$). Use $\nabla_XZ=\nabla^s_XZ+\mathrm{II}(X,Z)$ and $(\nabla_XU)^H=\tfrac12\tilde UX$ (N1(e)); $\nabla^s$ is the Levi-Civita connection of $c_Hg$.
- (c) $\mathcal L_E\mathrm{Ric}=\partial_\theta\mathrm{Ric}(\mathrm{Z1}(e^\theta c_H))$ by N1(b) and naturality. At $p$ it equals $c_H\partial_{c_H}$ of F4(c): $-\tfrac1cg$ on HH, 0 on HV and VV. Then subtract $\mathrm{Ric}(\nabla_\cdot E,\cdot)+\mathrm{Ric}(\cdot,\nabla_\cdot E)$ with $\nabla E=\tfrac12P_H$.
- (d) Differentiate N1(d).
- (e) The fibres are totally geodesic and $F=\mathbb R\times S$ (F1).
- (f) *Computed* (output §N4: exact on random curved bases in signatures (4,0), (1,3), (2,2), with and without cubic terms; all other items (a)–(e) are also re-checked there). It is reduced to finitely many evaluations as follows.
  - The isotropy element $-1\in O(g_x)$ acts as $-1$ on $H$ and $+1$ on $V$, and sends $\nabla R\mapsto-\nabla R$. So the components $(\nabla_U\mathrm{Ric})_{VV}$ and $(\nabla_U\mathrm{Ric})_{HH}$ (an even number of $H$-slots) are even in $\nabla R$.
  - At $p$, $\nabla\mathrm{Ric}$ is affine in the base data. Ric(Z2) is affine by O'Neill, with $A=0$ at $p$; the K1-type exact formula for $\mathrm{Ric}(\mathrm{Z2}+\delta)$ contains at most one $\nabla^3\delta$ factor, $\nabla^3\delta_p$ is affine in $R$, and $\nabla\delta_p$, $\nabla^2\delta_p$ are $R$-free (K2).
  - The $R$-free part is the flat value, which is 0 by N2.
  - The $R$-linear part is an $SO(g)$-equivariant linear map on algebraic curvature tensors. Over ℂ its four constituents ($W^\pm$, traceless Ricci, scalar) are pairwise non-isomorphic. By Schur, its value at one $R$ with all four components nonzero determines it; §N4 checks this genericity.
  - The affinity step is argued here, not written out in full; hence the status. ∎

**Theorem N7 (degree-1 obstruction).** *[proven, given N6(f)]* Assume $(\star)$ on a neighbourhood of $p\in s(X)$ with $r(p)=r_0\neq0$ (ν arbitrary). Put $v=\nabla r(p)$, $n=\nabla\nu(p)$. Then:

(a) $v=v_E\,E$ and $n=n_E\,E$, with $v_E=\mu_1r_0\big(2b-\tfrac{\sigma_g}2\big)$ and $n_E=-\mu_1\beta_Sv_E$;

(b) $g$ is Einstein at $x_0=\pi(p)$, and $d^\nabla\mathrm{Ric}_g(x_0)=0$;

(c) $\sigma_g^2=\dfrac{4b(4\mu_1bc-1)}{\mu_1c}=\dfrac{10(1-5\lambda)(13-56\lambda)}{7c_V^2(1-4\lambda)^2}$.

(d) Consequently, if $(\star)$ holds on a neighbourhood of $p$ in $Y$: $\mathrm{Ric}_g(x_0)=0$, $\lambda\in\{\tfrac15,\tfrac{13}{56}\}$, and $\nabla r=2\mu_1b\,r\,E$ at the points of Σ near $p$.

*Proof.* The Ricci identity applied to $(\star)$ gives, at every point where $(\star)$ holds nearby,
$$(\ast)\quad R(X,Y)v-\mu_1\big[\langle v,X\rangle\widehat{\mathrm{Ric}}\,Y-\langle v,Y\rangle\widehat{\mathrm{Ric}}\,X\big]-\big[\langle n,X\rangle Y-\langle n,Y\rangle X\big]=\mu_1r_0\big[(\nabla_X\widehat{\mathrm{Ric}})Y-(\nabla_Y\widehat{\mathrm{Ric}})X\big].$$
Write $v=v_H+v_EE+v_0$ ($v_0\in V_0$), likewise $n$. Components of $(\ast)$ at $p$, using N1 and N6:
1. *$E$-components.* $\langle R(X,Y)v,E\rangle=-\langle R(X,Y)E,v\rangle$ is $R$-free by N1(e), and the right side vanishes by N6(d) and $\mathrm{Ric}_{HV}(p)=0$.
   - The pair $(X,U)$, $U\in V_0$, gives $U(X,v_H)=0$ for all $X$ and $U$, so $v_H=0$.
   - The pair $(X,E)$ gives $n_H=0$.
   - The pair $(U,E)$ gives $n_0=\tfrac{\mu_1}cv_0$.
2. *$(U,W)$, $U,W\in V_0$, output in $V_0$.* By N6(e), (f): $R^S(U,W)v_0+\mu_1b\big[\langle v_0,U\rangle W-\langle v_0,W\rangle U\big]=0$, since $\beta_S+1/c=-b$. Take commuting independent $U,W$ (as in L7). If $b\neq0$ this gives $v_0\perp U$ for generic $U$; if $b=0$, the trivial commutant of L3 gives the same. So $v_0=0$ and $n_0=0$.
3. *$(X,U)$, $U\in V_0$, output in $H$.* With $v=v_EE$:
$$-\tfrac{c_H}4v_E\,U(X,Z)=\mu_1r_0\Big[-\tfrac{c_Hb}2U(X,Z)+\tfrac12\mathrm{Ric}_g(X,\tilde UZ)\Big]\qquad(\text{N1(e), N6(b), (f)}).$$
   The left side is symmetric in $X,Z$. So $[\mathrm{Ric}_g^\sharp,\tilde U]=0$ for every traceless $\tilde U$, i.e. $\mathrm{Ric}_g=\tfrac{\mathrm{scal}_g}4g$, and then $v_E=\mu_1r_0(2b-\sigma_g/2)$.
4. *$(E,U)$, output in $V_0$*: $-\mu_1v_Ec\beta_S-n_Ec=0$ (N6(c), (d)).
5. *$(X,E)$, output in $H$*: $-\tfrac14v_E-\mu_1v_E(1+c\beta_S)+\tfrac{\mu_1c\sigma_g}4v_E=\tfrac{\mu_1r_0\sigma_g}8$ (N1(e), N6(b), (c)). With $1+c\beta_S=-bc$ and step 3 this is $2kb^2-\tfrac b2-\tfrac k8\sigma_g^2=0$, $k=\mu_1c$, which is (c).
6. *$(X,Y)$, output in $H$*: $d^\nabla\mathrm{Ric}_g=0$ (N6(b); the left side vanishes because $v_H=n_H=0$).

(d): apply (a)–(c) at the points $(x_0,e^\theta g(x_0))\in\Sigma$ near $p$. There $\mathrm{Ric}_{e^\theta g}=\mathrm{Ric}_g$ and $\sigma_{e^\theta g}=e^{-\theta}\sigma_g$, while $b,c$ do not depend on θ. So $e^{-2\theta}\sigma_g^2$ is constant on an interval, hence $\sigma_g=0$ and $b(4\mu_1bc-1)=0$. Now $b=0$ iff $\lambda=\tfrac15$, and $4\mu_1bc=-\tfrac{56}9(1-5\lambda)=1$ iff $\lambda=\tfrac{13}{56}$. ∎

Output §N5 confirms Theorem N7 on the full linear system $(\ast)$ in $(v,n)$ (1274 equations, 28 unknowns, rank 28):
- *flat base:* solvable exactly for $\lambda\in\{\tfrac{13}{56},\tfrac15\}$ among seven values;
- *random Weyl tensor, $\nabla R(x_0)=0$:* in (4,0), (1,3), (2,2), solvable at $\tfrac{13}{56}$ and $\tfrac15$ with the predicted unique solution, and not at $\tfrac17$, $\tfrac12$;
- *constant curvature at $\lambda=\tfrac17$ (and at $\lambda=\tfrac{13}{30}$, on the (7,7) side), $c_H=2$, $c_V=3$:* solvable for the two values $\sigma_g=\pm\tfrac{10}9$ (resp. $\pm\tfrac{65}{33}$) of (c). It is not solvable for $2\sigma_g$, and not at the point $(x_0,2g(x_0))$ of Σ;
- *random $R$, an Einstein $R$ with generic $\sigma_g$, a traceless-Ricci $R$, and a Weyl $R$ with non-Codazzi $\nabla\mathrm{Ric}$:* not solvable.

So at a single point of $s(X)$ the first-order (3-jet) conditions **can** be met by special Einstein bases, also for $\lambda>\tfrac14$. It is their persistence along Σ that forces $\mathrm{Ric}_g=0$ and quantises λ.

**Proposition N8 ($\lambda=\tfrac{13}{56}$ is excluded).** *[proven, given N6(f)]* Under (H1), $(\star)$ and (Tr) have no solution on a neighbourhood of a point of $s(X)$ when $\lambda=\tfrac{13}{56}$.

*Proof.*
- By N7(d), at the points $q_\theta=(x_0,e^\theta g(x_0))$ of Σ near $p$: $\mathrm{Ric}_g(x_0)=0$ and $\nabla r=\tfrac r{2c}E$ (here $4\mu_1bc=1$). So $Er=\tfrac r2$, $r(q_\theta)=r_0e^{\theta/2}$ and $|dr|^2=r^2/(4c)$.
- By N1(b) and $(\star)$, $\nabla^2r(E,E)=E(Er)$ and $\mathrm{Ric}(E,E)=-1$, so $c\,\nu=E(Er)+\mu_1r=r(\tfrac14+\mu_1)$ on Σ.
- F4(b) for $e^\theta g$ gives $s_Y=\sigma_1/c_V$ on Σ.
- Inserted in (Tr), $r\,(\nu_{\rm Tr}-\nu)$ is a polynomial in $r$ with constant term $-A/13\neq0$ (output §N6). It cannot vanish for all $r$ in the interval $r_0e^{\theta/2}$. ∎

**Proposition N9 ($\lambda=\tfrac15$, line type, definite base).** *[proven, given N6(f); α computed exactly]* Under (H1) in (14,0), let $(\star)$ hold on a neighbourhood of $p\in s(X)$ with $\lambda=\tfrac15$. Then the Weyl tensor of $g$ vanishes at $x_0$, so $R_g(x_0)=0$.

*Proof.*
1. By N7(d), $\mathrm{Ric}_g=0$ near $x_0$ and $\nabla r=0$ on Σ near $p$. So $Er\equiv0$ on the $E$-orbit through $p$, and $c\nu(p)=\mu_1r_0$ (as in N8). At $\lambda=\tfrac15$, Z1 is Einstein at $p$ with constant $\kappa_E=-\tfrac1c$ (F6), so $\nabla^2r(p)=0$.
2. Restrict $(\star)$ to the unimodular slice $S$ through $p$, which is totally geodesic in $Y$. With $w=r|_S$: $\mathrm{Hess}^Sw=\psi\,g_S+\mu_1w\,\hat{\mathcal A}$, where $\psi=\mu_1\kappa_Ew+\nu$ and $\hat{\mathcal A}$ is the $TS$-part of $g_F^{-1}(\mathrm{Ric}_{VV}-\kappa_Eg_F)$. Then $\hat{\mathcal A}(p)=0$, $\nabla\hat{\mathcal A}(p)=0$ (N6(f)), $\nabla w(p)=0$ and $\psi(p)=0$.
3. Differentiating the Ricci identity of $S$ once at $p$, exactly as in L7, all $c_0$- and $c_1$-terms drop (they carry $\nabla^2w(p)=0$). What is left is
$$\mathrm{Hess}\,\psi(X,Y)Z-\mathrm{Hess}\,\psi(X,Z)Y+\mu_1r_0\,q(X,Y,Z)=0,\qquad q(X,Y,Z)=(\nabla^2_{X,Y}\hat{\mathcal A})Z-(\nabla^2_{X,Z}\hat{\mathcal A})Y,$$
   for $X,Y,Z\in T_pS$. $\mathrm{Hess}\,\psi(p)$ is an unknown symmetric form.
4. Pair with $T':=c_0-\tfrac{\langle c_0,c_1\rangle}{\langle c_1,c_1\rangle}c_1$ (L7's $c_0,c_1$ and pairing, now with the metric $g_\aleph$). By invariance, $\langle N\wedge,T'\rangle=0$ for every symmetric $N$: the trace part is a multiple of $c_1$, and an invariant functional on traceless symmetric forms on the irreducible $T_pS$ vanishes.
5. Hence $\langle q,T'\rangle=0$. Now $\langle q,T'\rangle$ is a quadratic form in the Weyl tensor $W=R_g(x_0)$; $\mathrm{Ric}_{VV}$ over $x_0$ is a polynomial of degree ≤2 in $R_g(x_0)$. It is $O(4)$-invariant, by the naturality of Z1 under base isometries, including orientation-reversing ones. So it equals $\alpha|W|^2$.
6. Output §N7 computes $\alpha$ exactly: $\alpha=\dfrac{45}{8c_H^2c_V}$ for $(c_H,c_V)=(2,3),(1,1),(3,\tfrac12)$, with the same value for two random Weyl tensors at $(2,3)$. (By the scaling of Z1 under $(c_H,c_V)\mapsto(tc_H,tc_V)$ and the base dilation of F7, α is a monomial in $c_H,c_V$, so this formula holds for all $c_H,c_V$; only $\alpha\neq0$ is used.) For the definite base $|W|^2=0$ gives $W=0$. ∎

So the Weyl tensor first enters at order $\varepsilon^2$ in the fibre (degree 2 of the equation), as in Z2 (L7). The pairing differs from L7's: there $T=4c_0-3c_1$ killed the A-terms, here $T'$ kills the concircular part.

**Proposition N10 (Conjecture F9 for Z1 near $s(X)$, definite base).** *[proven; computer-assisted: α is one exact computation, Schur does the rest]* Let the base be definite and let Z1 be Einstein on a neighbourhood of a point $p\in s(X)$. Then $\lambda=\tfrac15$ and $R_g(\pi(p))=0$. Together with N2: near $s(X)$, Z1 over a definite base is Einstein on an open set iff $\lambda=\tfrac15$ and the base is flat there.

*Proof.* F6 gives $\lambda=\tfrac15$ and $\mathrm{Ric}_g=0$ on an open set. Einstein near $p$ means $\hat{\mathcal A}\equiv0$ on $S$ near $p$, so $q=0$ and $\alpha|W|^2=0$ (N9, step 5). ∎

**Proposition N11 (flat base: no line solution on any open set).** *[proven]* If ℵ is flat on $\pi(W)$, $W\subset Y$ open, then $(\star)$ and (Tr) have no solution on $W$ with $dr\neq0$.

*Proof.*
- By N2, $\mathrm{Ric}|_{TS}=\beta_Sg_S$, so on each slice $S$ the restriction $r|_S$ is concircular. A concircular function on $S=SL(4)/SO(p,q)$ is constant: $\mathrm{Ric}_S=\rho_Sg_S$ with $\rho_S\neq0$ gives $d\psi\propto dr$; then $R^S(X,Y)\nabla r\propto\langle\nabla r,X\rangle Y-\langle\nabla r,Y\rangle X$, and commuting $X,Y$ give $\nabla r|_S=0$.
- In affine coordinates the HV part of $(\star)$ is $-\tfrac12(\tilde U)^\nu{}_\mu\partial_\nu r=0$ for all $U\in V_0$, so $\partial_xr=0$ and $r=f(s)$.
- Then $\nabla^2r=f''ds^2+\tfrac{f'}{2c}g_\aleph|_H$, and $(\star)$ requires $f'=2c\mu_1bf$ and $f''=\mu_1bcf$. So either $b=0$ and $f$ is constant ($dr=0$), or $4\mu_1bc=1$, i.e. $\lambda=\tfrac{13}{56}$, and $f=Ce^{s/2}$.
- For the latter, (Tr) multiplied by $r$ is a polynomial in $r$ with constant term $14A\neq0$, and $r$ ranges over an interval (output §N8). ∎

*Remark (the radial ansatz in Z1).* $t=f(s)$, $u=h(s)$ is a special case of the above. On a flat base, the traceless part of the equation has the exact radial solution $r=Ce^{s/2}$ at $\lambda=\tfrac{13}{56}$. It is Z1's analogue of the Z2 leaf structure, and the trace equation removes it. N7(a) also shows that the isotropy-invariant form of the normal gradient (J4's (J)) is **forced** at degree 1: $\nabla r(p)\in\mathbb RE$.

**Theorem N12 (Z1 neighbourhood solutions).** *[(a), (b) proven, given N6(f); (c) proven]* Assume (H1) with member $\circledast_\sharp$, $xW_A$ ($x\neq0$), or in class (A) or (B) of L8.

(a) **(7,7): no solution.** $\Upsilon=0$ has no solution on any connected open $U$ meeting $s(X)$.

(b) **(14,0):**
- $\lambda=\tfrac15$, and $g$ is Ricci-flat on $\pi(U\cap s(X))$.
- No point of $U\cap s(X)$ is of type (ii), i.e. every point lies in $\overline{\mathrm{int}\,Z}$.
- For class (A), $(t,u)$ is constant on $U$, $T_e\equiv0$, and $g$ is flat on $\pi(U\cap s(X))$.
- For $\circledast_\sharp$, $xW_A$ and class (B): if $(t,u)$ is constant near $p\in U\cap s(X)$ — in particular whenever $t,u$ are real-analytic — then $T_e\equiv0$ near $p$ and $g$ is flat near $\pi(p)$.

(c) **Existence.** Let $g$ be flat (base (4,0)), $\lambda=\tfrac15$ and $c_V\kappa^2\ge9360$ ($c_V>0$). Then the constant $(t,u)$ of I7(b) ($u=0$ and a root of $312s^2+\kappa s+\tfrac{15}{2c_V}=0$, in I7(b)'s variables), with $T_e\equiv0$, solve $\Upsilon=0$ on all of $Y$. $\det M\neq0$ at both roots, except at one root for the single value $c_V\kappa^2=277440/29$.

*Proof.*
- (a) By N3(c), class (A) has $U^*=\emptyset$, so every point of $U\cap s(X)$ is of type (i). For the other members use N4.
  - Types (i) and (iii) need $\lambda=\tfrac15<\tfrac14$.
  - Type (ii) needs $\lambda\in\{\tfrac15,\tfrac{13}{56}\}$ (N7(d)). Both are below $\tfrac14$.
- (b) Type (ii) is excluded:
  - $\lambda=\tfrac{13}{56}$ by N8.
  - For $\lambda=\tfrac15$, N9 at every point of the open type-(ii) set gives $R_g=0$ on an open set of $X$. On its preimage $(\star)$ has no solution with $dr\neq0$ (N11), and $dr\neq0$ on the dense set $U^*\cap N$.
  - Types (i), (iii) give $\lambda=\tfrac15$ and $\mathrm{Ric}_g=0$ (N4).
  - At a type-(i) point Z1 is Einstein on a neighbourhood, so $g$ is flat there (N10).
  - Real-analytic $(t,u)$: either $Z=U$, or $\mathrm{int}\,Z=\emptyset$ and every point of $U\cap s(X)$ is of type (ii), which is excluded.
- (c) Flat ℵ and $\lambda=\tfrac15$ make Z1 Einstein on $Y$ with $G=\tfrac{15}{2c_V}g_\aleph$ (N2, F7). With $T_e=0$, $\Upsilon=0$ is I6's $c_i=\tfrac{15}{2c_V}x_i$ at every point, and I7(b) solves it in (14,0) iff $c_V\kappa^2\ge9360$. On $u=0$, $\det M=\kappa^2(\kappa+224s)(\kappa+544s)$ (J4(d) proof). Output §N9 locates its zeros on the root set. ∎

**Answers.**
1. *Geometry of Z1.*
   - $E$ is a gradient with $\nabla E=\tfrac12P_H$; it is neither Killing nor homothetic, and no conformal Killing field is tangent to the fibres. It is still a Ricci eigenvector everywhere, with eigenvalue $-1/c$ (N1).
   - Over a flat base, $\mathrm{Ric}=\beta_Sg+4bc\,e^{-s/2}\nabla^2e^{s/2}$, which is Einstein iff $\lambda=\tfrac15$ (N2).
   - Off the section, the curvature corrections agree with Z2's O'Neill A-terms to second order in the fibre and differ at third order (N5, computed).
   - Near $s(X)$ over a definite base, Z1 is Einstein only when the base is flat and $\lambda=\tfrac15$ (N10).
2. *$\circledast_\sharp$ and $W_A$.*
   - In (7,7), no solution on any neighbourhood of a point of $s(X)$.
   - In (14,0), only $\lambda=\tfrac15$ over a Ricci-flat base is possible. Under analyticity, exactly the flat-base, $T_e\equiv0$ solutions of I7(b) occur, and they exist (N12).
   - These are the first neighbourhood solutions of $\Upsilon=0$ recorded in the audit's invariant/slaved sector; they follow from F7 and I7(b), but had not been stated (I7(b) left (E) on a neighbourhood "not established"). They have $\Lambda_4=0$ and a flat base, and they live in the Euclidean real form (14,0), not in (7,7).
3. *Other members.* N12 covers classes (A) and (B) of L8 in Z1. The members left open by L.3 remain open.
4. *The order.*
   - Degree 0 (2-jets) constrains nothing. Without J4's ansatz (J), $(\star)$ at $p$ just fixes the Hessian.
   - Degree 1 (3-jets) forces an Einstein base at each point of Σ. Hence $\mathrm{Ric}_g=0$ and $\lambda\in\{\tfrac15,\tfrac{13}{56}\}$ (N7).
   - The trace equation along Σ then removes $\tfrac{13}{56}$ (N8).
   - Degree 2 (order $\varepsilon^2$ in the fibre) brings in the Weyl tensor through $\alpha|W|^2$ (N9).

## N.2 What is verified, precisely

| Statement | How (section of `N_z1_neighbourhood.out`) | Scope |
|---|---|---|
| N3(d), inputs | §N0: L's `reduced` on the line ($\circledast_\sharp$, $W_A$, both real forms) equals the line formulas at $N=1$; reparametrisation identities; $\mu_1=-\tfrac7{18}$; the ν formula; constant term $A$ of $r\Phi$ (sympy, κ, $N$ symbolic) | both real forms |
| N1 | proof; §N1: brute-force metric 2-jet and Ricci of Z1 at $y_0\neq g(x_0)$, cubic random bases in general coordinates ($\Gamma(x_0)\neq0$): $T=0$, $E^\flat=c\,ds$, $\nabla E=\tfrac12P_H$, $\lvert E\rvert^2=c$, $\mathrm{Ric}(E,\cdot)=-E^\flat/c$ | (4,0) λ=0; (1,3) λ=½; (2,2) λ=⅔ |
| N2 | proof; §N2: brute force at $y_0\neq\eta$, five λ; $G=\tfrac{15}{2c_V}g$ at λ=⅕ | (4,0) |
| N5 | §N3: Taylor series of Ric(Z1) and Ric(Z2) in the 10 fibre variables (truncated-series coefficient ring in `F_jets`), random $R$ and cubic terms | (4,0) to order 3; (1,3), (2,2) to order 2 |
| N6 | proof of (a)–(e); §N4: Riemann and the full 1-jet of Ric at $p$ (shifted base point, 14 series variables), all blocks against the closed forms; genericity of $R$ for the Schur step | (4,0), (1,3), (2,2); with and without $\nabla R$ |
| N7 | proof; §N5: exact ranks of the 1274×28 system $(\ast)$; predicted solutions; the special $\sigma_g$; the point $(x_0,2g)$ of Σ; failing bases | flat (7 λ); Weyl in three signatures; constant curvature at λ=1/7, 13/30 |
| N8 | proof; §N6: sympy | — |
| N9, N10 | proof; §N7: second fibre derivatives of $\hat{\mathcal A}$ from the series, the pairings $\langle c_0,c_0\rangle,\langle c_0,c_1\rangle,\langle c_1,c_1\rangle$, $\langle c_1,T'\rangle=0$, $\langle N\wedge,T'\rangle=0$ for random traceless $N$, $\langle q,T'\rangle=\alpha\lvert W\rvert^2$ with $\alpha\neq0$ | (4,0); $(c_H,c_V)=(2,3),(1,1),(3,\tfrac12)$ |
| N11 | proof; §N8: radial ODE and trace equation (sympy) | — |
| N12(c) | §N9: discriminant and $\det M$ (sympy); Einstein flat model in §N2 | (14,0) |

Run: `cd computation/geometric-unity && python3 N_z1_neighbourhood.py > N_z1_neighbourhood.out 2>&1`. Result: 44 checks, 44 True, each asserted; 450 s on one core, a second run 473 s, with identical output apart from the timing annotations.

## N.3 Limits (layers b, d)

- **One computed input.** N6(f), the vertical derivative of Ric at $p$ in the $V_0$-directions, is computed exactly on generic bases. It is reduced to those evaluations by Schur, but the affine dependence on the base data that this needs is argued, not written out. Theorem N7 and everything resting on it (N8, N9, N12(a), (b)) depend on it. N1–N4, N10, N11 and N12(c) do not.
- **N5** (agreement with Z2's A-terms to second order in the fibre) is computed only. It is a side result, used nowhere.
- **Type-(iii) points.** In (14,0), without analyticity, $U\cap s(X)$ could consist of limits of Einstein pieces that do not reach it. Only $\mathrm{Ric}_g=0$ is proven there, not flatness.
- **Indefinite bases in F9.** N10 uses $|W|^2$, which is definite only for a definite base. For a Lorentzian or (2,2) base at $\lambda=\tfrac15$ (signatures (8,6), (10,4), …, outside the setting) it gives $|W|^2=0$ only. F9 for those bases, and F9 on open sets away from $s(X)$, remain open.
- **Members.** As L.3: the non-elliptic members other than the multiples of $W_A$, $x=\pm y$ in (7,7). The Z2 treatment of these members in round 9 (`M-z2-remaining-members.md`) was not transferred to Z1. For $(x,y)=0$, M8 states $T=0$ for every metric.
- **Hypotheses kept:** $n=4$; the real forms of J1–J2; $\kappa\neq0$; $T_e\in\mathcal V$ slaved with $\det M\neq0$; $T_o$ in the invariant sector.
- **Reconstruction:** reading R-Z of F with Z1. This note changes no layer-(b) choice.

## N.4 Consequences for gaps and claims

- **GU-G14.**
  - Under the hypotheses of J5/L8 with Z1, for $\circledast_\sharp$, $W_A$ and classes (A), (B):
    - no neighbourhood solution exists in (7,7);
    - in (14,0) a neighbourhood solution forces $\lambda=\tfrac15$ and a Ricci-flat base, and (for analytic $t,u$, or for class (A)) it is exactly a flat-base solution with $T_e\equiv0$ (N12).
  - These solutions exist (F7 with I7(b), $c_V\kappa^2\ge9360$). They reproduce 4D vacuum GR with $\Lambda_4=0$ only in the trivial sense of a flat base.
  - So in Z1 as in Z2, the neighbourhood equations do not reduce to a 4D theory with curved solutions in this sector. The obstruction now sits at the 3-jet (Ricci, along the radial set Σ) and at the 4-jet (Weyl, along the fibre slice).
- **GU-054.** Status unchanged (C). Round 5's (E) solutions (I7(b)) extend to all of $Y$ when the base is flat (N12(c)). Over a curved definite base, (E) holds on no neighbourhood of a point of $s(X)$ (N10).
- **F9.**
  - For Z1: proven for flat ℵ (Einstein iff $\lambda=\tfrac15$, N2) and near $s(X)$ for definite bases (N10).
  - Open for indefinite bases and on open sets away from $s(X)$.
  - The statement in F after Theorem F6, "So (E) holds only along $s(X)$, not on a neighbourhood in $Y$", holds for the Weyl base tested there but not in general. For flat ℵ at $\lambda=\tfrac1{n+1}$, (E) holds on all of $Y$ (F7).
- **J.3, L.3.** The item "Z1" is settled for the members of L8 except for type-(iii) points (N.3). No statement made before this note is contradicted, apart from the scope of the F sentence above.
