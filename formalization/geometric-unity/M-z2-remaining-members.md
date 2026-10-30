# M. Υ = 0 on a neighbourhood in Z2: the members round 8 leaves open

Claims: GU-054 (context). Gaps: GU-G14, GU-G16. Code: `computation/geometric-unity/M_remaining.py` → `.out` (sections M0–M8; exact: Python Fractions, sympy rationals and Gaussian rationals; fixed seeds; 270 s on one core, a second run 260 s, identical apart from the timing annotations). It imports `L_nonradial.py` (and through it `J_slaved_sector.py`, `F_jets.py`, …) unmodified. Output sections: §M0 inputs, §M1 Lemma M1, §M2 Proposition M2, §M3 Proposition M3, §M4 Proposition M4, §M5 Lemma M5, §M6 Theorem M7, §M7 Proposition M8, §M8 Proposition M9. Record: `metrics/geometric-unity/M-z2-remaining-members.yaml`. Builds on L (L1–L9, as amended after REVIEW6-GU), J (J1, J2, Cor. J2.1), J5d (L0, L1), F (F1, F3, L-F1). Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## M.0 Question and setting (layer b)

L8 and L9 decide $\Upsilon=0$ on open sets under the hypotheses (H) of L.0 for the members with $\ell\neq0$, $N\neq0$ and for the multiples of $W_A$. L.3 lists three classes as open: $\ell=0$ other than the multiples of $W_A$; $N=0$ in (7,7); $(x,y)=0$. This note treats all three.

**Hypotheses.** (H) of L.0: Z2, $n=4$, $\lambda\neq\tfrac14$, $c_Hc_V\neq0$, the real forms (14,0) and (7,7) (base (1,3), (3,1) or (2,2) in (7,7)), $\kappa\neq0$, $T_o=t\Phi^1+u\Phi^1\Gamma$ with $t,u$ arbitrary, $T_e=\Sigma(dt,du)^\top$ the slaved field in $\mathcal V$ with $\det M\neq0$, $U$ connected and open, $\Upsilon=0$ on $U$. (H) includes $(x,y)\neq0$ (as in J5; added to L.0 after REVIEW6-GU, R8-5(c)); M8 drops only that. "Elliptic" is used only in (14,0) (R8-3).

**Notation.** As in L.0: $\ell=xz+\epsilon yw$, $N=x^2+\epsilon y^2$, $\delta=ty-ux$, $U^*$, $E'=yE^{(1)}-xE^{(13)}$, its Hessian row $h=(y,-x)H$, kernel field $v$ and residual $R'$ (L5). In addition:
- $\sigma:=tx+\epsilon uy$ and $\theta:=tw-uz$ (the entries of $M$, J1); $J(x,y)=(-y,x)$, $K(x,y)=(y,x)$.
- $G_\epsilon$ is $SO(2)$ for $\epsilon=1$ and $SO^+(1,1)$ for $\epsilon=-1$, acting by $R=\begin{pmatrix}c&-\epsilon s\\ s&c\end{pmatrix}$, $c^2+\epsilon s^2=1$, i.e. $R^\top\mathrm{diag}(1,\epsilon)R=\mathrm{diag}(1,\epsilon)$.
- $K_{\rm tr}:=v^\top(\mathcal Q''_{\rm sym}+\mathcal Q''')v$, with $\mathcal Q'''=y\mathcal Q'^{(1)}-x\mathcal Q'^{(13)}$: the trace coefficient of $E'$ along an integral curve of $v$.
- (C1) is L7's order-1 condition, $c_0+\mu'c_1+q(\mathfrak B)=0$; (C1)$_0$ is (C1) with $\mu'=0$.
- Real forms as in L.0. Where a sign matters (M4), the variables are $\tau:=xt/\kappa$ and $c:=w/x$, which are real in both real forms ($xt=x't'$ and $w/x=w'/x'$ in (14,0)).

**The cone $\ell=0$.** In (14,0) it is $x'z'+y'w'=0$, in (7,7) $xz-yw=0$: a quadric cone of signature (2,2) in $\mathbb R^4$. Its points with $(x,y)\neq0$ are $(z,w)=c\,J(x,y)$ in (14,0) and $(z,w)=c\,K(x,y)$ in (7,7), $c\in\mathbb R$; there $xw-yz=cN$. $c=0$ is $\mathrm{span}(W_A,W'_A)$. The remaining points have $(x,y)=0$.

## M.1 Results (layer c)

**Lemma M1 (symmetries of the reduced system).** *Proven* (identities computed, §M1).

(a) For $R\in G_\epsilon$ let $Rp:=(R(x,y),R(z,w))$. Then $M(Rp;R(t,u))=M(p;(t,u))$ and $S(Rp)R=S(p)$. Moreover the columns $(\ell^{(1)}_k,\ell^{(13)}_k)$, $(\ell^{(3)}_k,\ell^{(11)}_k)$, $(b^{(1)},b^{(13)})$, $(b^{(3)},b^{(11)})$, $(c_1,c_{13})$ and $(x_1,x_{13})$ all transform by $R$. Hence $T_e$ and $\det M$ are invariant, and $(E^{(1)},E^{(13)})$, $(\omega^{(3)},\omega^{(11)})$ transform by $R$ ($dt\wedge du$ is invariant). So $(t,u)$ solves $\Upsilon=0$ under (H) for $p$ iff $R(t,u)$ does for $Rp$, on the same $(U,g_\aleph)$.

(b) In (7,7) let $p':=-Kp=(-y,-x,-w,-z)$ and $(t',u'):=(u,t)$. Then $\det M$ is invariant. Also $H'_{1\phi'}=-H_{13,K\phi}$, $H'_{13\phi'}=-H_{1,K\phi}$, and the same holds for $\mathcal Q,\mathcal Q'$; $c'_1=c_{13}$, $c'_{13}=c_1$, $x'_1=-x_{13}$, and the $\varpi'$ are those of $p$ with $1\leftrightarrow13$, $3\leftrightarrow11$. Replacing $g_\aleph$ by $-g_\aleph$ flips only the terms $c_ig$; the tensors $G$, $\nabla^2\phi$, $\Box\phi\,g$, $d\phi\otimes d\psi$ and $\langle d\phi,d\psi\rangle g$ are invariant. Hence $E'^{(1)}[-g_\aleph]=-E^{(13)}[g_\aleph]$ and $E'^{(13)}[-g_\aleph]=-E^{(1)}[g_\aleph]$. So $(t,u)$ solves $\Upsilon=0$ for $p$ on $(U,g_\aleph)$ iff $(u,t)$ does for $-Kp$ on $(U,-g_\aleph)$. Here $-g_\aleph$ is the Z2 metric with constants $(-c_H,-c_V)$, the same $\lambda$ and signature (7,7).

*Proof.* The reduced system (J1, J2 and "Substituting the slaving") is assembled from $S$, $M$, the two tables, $c_i$ and $x_i$ alone. The listed identities are polynomial (rational in the group parameter, symbolic one-parameter group in (a)) and checked exactly in §M1. ∎

(a) is presumably right multiplication by $e^{\varphi\Gamma}$; only its action on the reduced system is used.

**Proposition M2 (the cone: orbits and $\Lambda^2$ coefficients).** *Proven* (§M2).

(a) *Orbits.*
- (14,0): every member with $\ell=0$ and $(x,y)\neq0$ is $G_+$-equivalent to $(x,0,0,w)=x(W_A+c\,W'_C\mathrm{Ad}_\Gamma)$, $x\neq0$, $c=w'/x'$. This includes $W'_A\sim W_A$, $\mathrm{span}(W_A,W'_A)\sim$ multiples of $W_A$, and $(1,0,0,1)$ ($c=1$).
- (7,7), $N>0$: $G_-$-equivalent to $(x,0,0,w)$.
- (7,7), $N<0$: equivalent to $(0,y,z,0)$, which M1(b) maps to $(-y,0,0,-z)$ on $(U,-g_\aleph)$. For $z=0$ (e.g. $W'_A$) this is a multiple of $W_A$.
- (7,7), $N=0$: the members $(x,\pm x,\pm cx,cx)$. Here $(z,w)\parallel(x,y)$, so they are L8(B)'s family at $N=0$.

(b) *$\Lambda^2$ coefficients.* For $(x,0,0,w)$, $x\neq0$: $\varpi^{(13)}=0$, and $\varpi^{(1)}=\mp\mathcal Q^{(13)}_{tt}$ (upper sign (14,0)) is a nonzero polynomial times $\kappa ux^2/(\det M)^2$. So L1 applies on the whole cone with $N\neq0$ (by M1 the pairs of $\varpi$'s transform as vectors). For $N=0$, $\ell=0$ in (7,7), and for $(x,y)=0$, all four $\varpi$ vanish identically, so L1 does not apply there. It is not needed: M7(a) and M8 do not use it. This agrees with the referee's observation table (REVIEW6-GU R8-1, `xcheck_round8.out` §X1: $\varpi\not\equiv0$ for $W'_A$, $W_A+W'_C\mathrm{Ad}_\Gamma$, $W'_A+W_C\mathrm{Ad}_\Gamma$; $\varpi\equiv0$ for $W_C\mathrm{Ad}_\Gamma$, $W'_C\mathrm{Ad}_\Gamma$ and, in (7,7), $W_A+W'_A$, $(1,1,1,1)$), and it settles the members listed there without L1 (M3, M7(a), M8).

**Proposition M3 (the first-order channel of $(x,0,0,w)$).** *Proven* (§M3). Let $p=(x,0,0,w)$, $x\neq0$, in either real form. Exactly:
- $\Sigma$ has no $du$ column, and $H_{1u}=H_{13t}=H_{13u}=0$;
- all $E^{(13)}$ coefficients except the following vanish, and
$$E^{(13)}=u\,(\kappa+208xt)\,g+\mathcal Q^{(13)}_{tt}\big(dt\otimes dt-|dt|^2g\big).$$
Under (H), $\Upsilon=0$ on $U$ forces $u\equiv0$. Then $\Upsilon=0$ is the single equation
$$E^{(1)}=-xG+(\kappa t+312xt^2)\,g+H(t)\big(\nabla^2t-\Box t\,g\big)+q(t)\,dt\otimes dt+q'(t)\,|dt|^2g=0,$$
$$H(t)=\frac{576\big(\kappa(x^2+\epsilon w^2)+544\,tx(x^2-\epsilon w^2)\big)}{\kappa^2+768\kappa tx+121856x^2t^2-147456\epsilon w^2t^2},\qquad \det M|_{u=0}=\kappa^2\big(\kappa^2+768\kappa tx+121856x^2t^2-147456\epsilon w^2t^2\big),$$
with $q,q'$ printed in §M4.

*Proof.*
1. If $dt=0$ on an open $V\subset U$, then $T_e=\Sigma_t\,dt=0$ and $E^{(1)}=-xG+c_1g$ on $V$. So $\mathrm{Ric}^Y\propto g$ there, which L2(a) excludes. Hence $\{dt\neq0\}$ is dense.
2. On $\{dt\neq0\}$ the traceless part of $E^{(13)}$ is $\mathcal Q^{(13)}_{tt}(dt\otimes dt)_0$. Since $(dt\otimes dt)_0\neq0$ in dimension 14, also for null $dt$, this forces $\mathcal Q^{(13)}_{tt}=0$.
3. Then $E^{(13)}=u(\kappa+208xt)g$. On the open set $\{dt\neq0,\ u\neq0\}$ this forces $t=-\kappa/(208x)$, against $dt\neq0$. So $u=0$ on a dense set, hence on $U$.
4. With $du=0$ every $du$ term and $dt\wedge du$ drop. ∎

For $w=0$ this is L8 step 3 ($xW_A$). The mirror statement for $(0,y,z,0)$ in (7,7) is checked directly in §M3: $E^{(1)}=t(\kappa-208yu)g+\mathcal Q^{(1)}_{uu}(du\otimes du-|du|^2g)$, which forces $t\equiv0$. It is also the M1(b) image of M3.

**Proposition M4 (the remaining equation).** *Computed* (§M4). Where $H(t)\neq0$, the reparametrisation $r=F(t)$, $F''/F'=q/H$, turns $E^{(1)}=0$ into $(\star)$ $\nabla^2r=\mu(r)\mathrm{Ric}^Y+\nu g$ with $\mu=xF'/H$ (as in L8 step 5), and
$$\mu'(r)=\frac{x\,(q-dH/dt)}{H^2}=-\frac{(7-3\epsilon c^2)+(7616-3008\epsilon c^2-768c^4)\,\tau+(2071552-3226624\epsilon c^2+1155072c^4)\,\tau^2}{18\,\big((1+\epsilon c^2)+544(1-\epsilon c^2)\,\tau\big)^2}.$$
- $c=0$ (multiples of $W_A$): $\mu'=-\tfrac7{18}$, L8 step 5.
- $c\neq0$: $\mu'$ is a nonconstant rational function of $\tau$ and never identically zero, so $\mu$ is not linear in $r$. The discriminant of the numerator in $\tau$ is $16384c^2(c^2+35\epsilon)(6c^2-11\epsilon)^2$. In (14,0) it is $\ge0$, so the numerator has real roots for every $c\neq0$. In (7,7) it is negative for $0<c^2<35$, and the leading coefficient $1024(c^2+1)(1128c^2+2023)$ is positive, so there $\mu'<0$ everywhere.
- $H$ vanishes exactly at $\tau_H=(1+\epsilon c^2)/(544(\epsilon c^2-1))$, and nowhere if $\epsilon c^2=1$.

So for $c\neq0$ L9's continuation, which needs a linear $\mu$, is unavailable, and L7, which needs $\mu'(w(p))\neq0$, fails at the zeros of $\mu'$ and of $H$. M5 and M6 deal with both.

**Lemma M5 (fibre obstruction at $s(X)$ for any $\mu$).** *Proven* ((c) by an exact certificate, §M5). Let $p=g(x_0)\in s(X)$, let $\Omega\subset F_{x_0}$ be a neighbourhood of $p$, let $\mu$ be smooth near $w(p)$, and let $w$ on $\Omega$ satisfy $\mathrm{Hess}^Fw=\mu(w)(\mathrm{Ric}_F+\mathcal A_{x_0})$ with $\mu(w(p))\neq0$. Then:

(a) $\nabla w(p)\in\mathbb R\,\partial_s$, $\mu'(w(p))=0$, and $R=R_g(x_0)$ satisfies (C1)$_0$: $c_0+q(\mathfrak B)=0$.

(b) (C1)$_0$ implies $|R|^2=-2c_H^2/c_V^2$. In particular no such $w$ exists if $g(x_0)$ is definite.

(c) (C1)$_0$ has no solution with $\mathrm{Ric}(R)=0$ in bases (1,3), (3,1), (2,2): its Gröbner basis over $\mathbb Q$ is $\{1\}$, so there is none even over $\mathbb C$. In base (4,0) this follows from (b).

*Proof.*
- *Order 0.* Use L7's Ricci identity at $p$, where $\mathcal A(p)=0$ and $\nabla\mathcal A(p)=0$. If $\mu'(w(p))\neq0$, L7 gives $\nabla w(p)=0$. If $\mu'(w(p))=0$, then $R^F(Y,Z)\nabla w(p)=0$ for all $Y,Z$, and L3's commutant argument gives $\nabla w(p)\in\mathbb R\partial_s$. In both cases $Xw(p)=0$ for $X\in T_pS$.
- *Order 1.* Differentiate along $X\in T_pS$. The terms with $\mu''$ or $\mu'\nabla\mathcal A$ carry a factor $Xw(p)$, $Yw(p)$, $Zw(p)$ ($Y,Z\in T_pS$) or $\nabla\mathcal A(p)$, so they vanish. What remains is L7's (C1) with $\mu'=\mu'(w(p))$. Pairing it with $4c_0-3c_1$ gives $-360\mu'=0$ (L7). Hence $\mu'(w(p))=0$ and (C1)$_0$ holds.
- *(b).* Pair (C1)$_0$ with $c_1$. With $\langle c_0,c_1\rangle=18$, $\langle q(K\otimes K'),c_1\rangle=-9\,\mathrm{tr}(KK')$ and $\sum\epsilon_i\epsilon_j\mathrm{tr}R_{ij}^2=-|R|^2$ (L7; re-checked in all four base signatures, §M5) this gives $18+9(c_V/c_H)^2|R|^2=0$.
- *(c).* For the generic Ricci-flat tensor (10 Weyl parameters), the 2796 scalar equations of (C1)$_0$ have rank 12. The Gröbner basis of 12 independent ones is $\{1\}$. The factor $(c_V/c_H)^2$ is absorbed by rescaling $R$.
- *Independent check.* A brute-force Taylor solution, built as in L §L6 but with $\mu'=0$, is inconsistent for a random rational Weyl tensor with $|W|^2=-204$ and $(c_V/c_H)^2=1/102$, so that (b) holds. The control $R=R_1$ with the A-term multiplied by $-\tfrac1{12}$ is consistent. ∎

So L7 extends to every $\mu$: a point of $s(X)$ escapes the obstruction only if $\mu'(w(p))=0$ and $R_g(x_0)$ is indefinite, not Ricci-flat, and has $|R|^2=-2c_H^2/c_V^2$.

**Theorem M6 (the cone with $N\neq0$).** *Proven.* Assume (H), $\ell=0$, $(x,y)\neq0$ and $N\neq0$ (automatic in (14,0)).

(a) If $U$ meets $s(X)$, $\Upsilon=0$ has **no solution on $U$**.

(b) If $(z,w)=0$, i.e. on $\mathrm{span}(W_A,W'_A)$, there is **no solution on any nonempty open $U$**.

In (14,0), (a) covers every non-elliptic member with $(x,y)\neq0$, for example $W'_A$, $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$ and $(1,0,0,1)$.

*Proof.* By M2(a) and M1 it suffices to treat $p=(x,0,0,w)$, $x\neq0$, on $(U,\pm g_\aleph)$. The sign is needed only for $N<0$ in (7,7). Every ingredient used holds for all constants $c_H,c_V$ of the setting.
1. By M3, $u\equiv0$ and one equation $E^{(1)}[t]=0$ remains, with $\{dt\neq0\}$ dense.
2. $W:=\{dt\neq0,\ H(t)\neq0,\ \mu'(t)\neq0\}$ is open and dense ($H$ and $\mu'$ are nonzero rational functions of $t$, and $t$ is a submersion on $\{dt\neq0\}$). On each component of $W$, $r=F(t)$ satisfies $(\star)$ with $\mu'(r)\neq0$, so L6 gives $Er=0$ and $\nu=0$. Hence $Et=0$ on $U$.
3. (b): $c=0$ gives $\mu=-\tfrac7{18}r$ (M4), which is L8(B) and L9 for $xW_A$.
4. (a): let $p\in U\cap s(X)$, $x_0=\pi(p)$. Suppose first $H(t(p))=0$. Restrict $E^{(1)}(p)$ to $V_p\times V_p$, using $\mathrm{Ric}^Y|_{VV}(p)=\mathrm{Ric}_F(p)$ ($\mathcal A(p)=0$), $\mathrm{Ric}_F(E,\cdot)=0$ and $\mathrm{Ric}_F|_{V_0}=-c_V^{-1}g_\aleph$. This gives $-x\mathrm{Ric}_F+\alpha g_{VV}+q\,\vartheta\otimes\vartheta=0$, with $\vartheta=dt|_V$ and $\alpha=\tfrac x2\mathrm{scal}+c_1+q'|dt|^2$. The $V_0$ block (dimension 9) forces $\alpha=-x/c_V$. The $E$ component then gives $q\,\vartheta(E)^2=(x/c_V)|E|^2\neq0$, so $Et(p)\neq0$, against step 2.
5. Suppose instead $H(t(p))\neq0$. Near $p$, $r=F(t)$ is defined and $\nu=0$ (by continuity from $W$). So $w=r|_{F_{x_0}}$ solves $\mathrm{Hess}^Fw=\mu(w)(\mathrm{Ric}_F+\mathcal A_{x_0})$, since the fibres are totally geodesic (O'Neill). M5 gives $\mu'(t(p))=0$ and (C1)$_0$ at $x_0$. In (14,0) this is impossible by M5(b).
6. Hence every $p\in U\cap s(X)$ has $t(p)$ in the finite zero set of $\mu'$, and $t$ is constant on each component $O$ of $U\cap s(X)$. Then $dr=0$ on $O$:
   - tangentially, because $T_pO=H_p$ ($s(X)$ is horizontal) and $r|_O$ is constant;
   - vertically, by step 5 ($\nabla w\in\mathbb R\partial_s$) and $Et=0$.

   For $X,X'$ tangent to $O$, $\nabla^2r(X,X')=X(X'r)-dr(\nabla_XX')=0$. So $\mu\,\mathrm{Ric}^Y_{HH}=0$ on $O$, i.e. $\mathrm{Ric}_g=0$ on $\pi(O)$, since $A=0$ on $s(X)$. Together with (C1)$_0$ this contradicts M5(c). ∎

**Theorem M7 ($N=0$ in (7,7)).** *Proven.* Assume (H) in (7,7) with $y=\pm x\neq0$, and put $\mathrm{sgn}:=y/x$.

(a) If $\ell=0$ (then $(z,w)=c(x,y)$), $\Upsilon=0$ has **no solution on any nonempty open $U$**.

(b) If $\ell\neq0$, $\Upsilon=0$ has **no solution on any nonempty open $U$**.

*Proof of (a)* (§M6).
- $T_e=\Sigma_t(dt-\mathrm{sgn}\,du)=x^{-1}\Sigma_t\,d\sigma$.
- The derivative coefficients of $E^{(1)}$ and $\mathrm{sgn}\,E^{(13)}$ coincide. So $E^{(1)}-\mathrm{sgn}\,E^{(13)}=(c_1-\mathrm{sgn}\,c_{13})g=x^{-1}\sigma(\kappa+208\sigma)g$.
- Hence $\sigma\in\{0,-\kappa/208\}$, so $\sigma$ is constant on $U$, $d\sigma=0$ and $T_e=0$.
- Then $E^{(1)}=-xG+c_1g$ gives $\mathrm{Ric}^Y\propto g$, which L2(a) excludes. ∎

*Proof of (b)* (§M6).
1. L1 ($\ell\neq0$) gives $dt\wedge du=0$, and $U^*$ is dense (L8 step 2).
2. L8 step 3 (the concircular branch, L4) gives $\hat h=0$ on $U^*$. So $\gamma'\parallel v$, where $v$ is polynomial and nonzero wherever $\det M\neq0$ (L5(a)).
3. Along $v$, $R'=0$ (L5(a) with $N=0$), and in addition $K_{\rm tr}=0$ identically. So $E'=(yc_1-xc_{13})g=\delta(\kappa+208\sigma)g$. Since $\delta=\mathrm{sgn}\,\sigma$, $\sigma\in\{0,-\kappa/208\}$ is constant on $U$.
4. $v\cdot\nabla\delta=-\kappa\ell(\kappa+16\sigma)(\kappa+208\sigma)$. At $\sigma=0$ this is $-\kappa^3\ell\neq0$: an image curve tangent to $v$ and lying in that line would be constant. So $\sigma\equiv-\kappa/208$, and this line is $v$-invariant.
5. On the line ($t=\mathrm{sgn}\,u-\kappa/(208x)$, $dt=\mathrm{sgn}\,du$): $\det M=576\kappa^2\zeta^2/(169x^2)$ with $\zeta:=\kappa w+208\ell u$. Also $E^{(13)}=\mathrm{sgn}\,E^{(1)}$, and
$$E^{(1)}=-xG+c_1g+\frac{312x\ell}{\zeta}\big(\nabla^2u-\Box u\,g\big)-\frac{105456x\ell^2}{\zeta^2}du\otimes du+\frac{52728x\ell^2}{\zeta^2}|du|^2g .$$
   With $F'=|\zeta|^{-13/8}$ the function $r=F(u)\propto|\zeta|^{-5/8}$ is defined on all of $U$, is nonzero, and satisfies $\nabla^2r=-\tfrac5{12}r\,\mathrm{Ric}^Y+\nu g$.
6. L6 gives $Er=0$ and $\nu=0$. L9's argument (linear, $\mu_1=-\tfrac5{12}\neq0$, L7's linear case) gives $r\equiv0$ on $U\cap F_{x_0}$, a contradiction. ∎

*Remark.* L.3's "on the line the Hessian coefficient vanishes" refers to L8's line $\delta=0$, where the coefficient is $576N/(\kappa+224N\sigma)=0$. For $\ell\neq0$ that line is not $v$-invariant (step 4). The candidate line is $\sigma=-\kappa/208$ instead, where the Hessian coefficient $312x\ell/\zeta$ does not vanish.

**Proposition M8 ($(x,y)=0$).** *Proven* (§M7). Assume (H) except that $(x,y)=0$. Then $\Upsilon=0$ on a connected open $U$ iff $t=u=0$ on $U$, i.e. $T=0$. Conversely, $T=0$ solves $\Upsilon=0$ on every $U$ for every $g_\aleph$ of the setting: with $x_1=x_{13}=0$ the curvature drops from both channels (J2, I6).

*Proof.*
- $(z,w)=0$: then $\circledast=0$, $S=0$, $T_e=0$ and $E^{(i)}=c_ig=\kappa(t,u)_ig$.
- $(z,w)\neq0$: everything depends on $\theta$.
  - $\Sigma(dt,du)$ is a multiple of $d\theta$, and $\det M=\kappa^2(\kappa^2-147456\epsilon\theta^2)=:D$.
  - $E^{(1)}=\kappa t\,g+\epsilon w\,\Xi_0$ and $E^{(13)}=\kappa u\,g-z\,\Xi_0$, with
$$\Xi_0=\frac{576\kappa^3}{D}\big(\nabla^2\theta-\Box\theta\,g\big)+\epsilon\frac{184025088\kappa^5\theta}{D^2}\big(d\theta\otimes d\theta-\tfrac12|d\theta|^2g\big).$$
  - So $zE^{(1)}+\epsilon wE^{(13)}=\kappa(zt+\epsilon wu)g$, and $zt+\epsilon wu=0$.
  - If $z^2+\epsilon w^2=0$ (only in (7,7), $w=\pm z$), this already gives $\theta=0$.
  - Otherwise $(t,u)=\theta(\epsilon w,-z)/(z^2+\epsilon w^2)$, and one equation $\Xi:=\kappa\theta(z^2+\epsilon w^2)^{-1}g+\Xi_0=0$ remains.
  - Where $d\theta\neq0$, $F'=|\kappa^2-147456\epsilon\theta^2|^{-13/12}$ makes $r=F(\theta)$ concircular. L3(b) gives $\nabla dr=0$ and $|dr|^2=C$ constant.
  - The trace of $\Xi$ is then $\theta\big(\kappa/(z^2+\epsilon w^2)+92012544\epsilon\kappa^5|d\theta|^2/D^2\big)$, with $|d\theta|^2=C|\kappa^2-147456\epsilon\theta^2|^{13/6}$. For $\theta\neq0$ this forces $\theta^2$ constant if $C\neq0$, and is impossible if $C=0$.
  - So $d\theta=0$ wherever $\theta\neq0$, and then the trace gives $\theta=0$. Hence $\theta\equiv0$.
- $\theta\equiv0$ gives $T_e=0$, and then $E^{(1)}=\kappa t\,g$, $E^{(13)}=\kappa u\,g$ give $t=u=0$. ∎

So for $(x,y)=0$, $\Upsilon=0$ imposes no condition on $g_\aleph$ and forces the trivial connection in the sector.

**Proposition M9 (the first integrability condition off the section).** *Computed* (§M8); the consequence is proven from the computation.

Setting:
- Let $y_0=Pg(x_0)P^\top$. The fibre isometry $y\mapsto PyP^\top$ pulls the fibre equation $\mathrm{Hess}^Fw=\mu(w)(\mathrm{Ric}_F+\mathcal A_{x_0})$ at $y_0$ back to the same equation at $g(x_0)$, with $R_{ij}$ replaced by $P^\top R_{ij}P^{-\top}$.
- Normalise $\mu(w(y_0))=1$ (via $w\mapsto w/\mu(w(y_0))$, which leaves $\mu'$ unchanged).
- The 3-jet of a solution must satisfy a linear system in which $\mu'$ enters only through $\mu'\,dw$.
- Replace $\mu'\,dw$ by a free covector. Call the resulting system the *relaxed system*; it has 295 unknowns.

Findings:
- The relaxed system is **inconsistent**, exactly (ranks 294 and 295), at every sampled off-section point: two random curvature tensors and $R=R_1$, each with a random $P$, in bases (1,3), (3,1), (4,0), (2,2).
- It is consistent at $P=1$ (the section), as L7 requires.
- The Taylor coefficients of $s$ span a kernel of its coefficient matrix ($\mathrm{Hess}^Fs=0$), so the rank is at most 294 identically.

Consequences:
- Inconsistency holds on a nonempty Zariski-open set of $(R,P)$.
- For every curvature tensor outside a proper real-algebraic subset $\mathcal Z$ (in each base signature; $\mathcal Z\not\ni$ the sampled tensors), the fibre equation has, for every $\mu$ with $\mu\neq0$, no solution on any open subset of the fibre. Here $R$ stands for $(c_V/c_H)R_g(x_0)$.
- $0\in\mathcal Z$. The flat case is excluded separately: $\mathcal A\equiv0$, and the Ricci identity on the symmetric fibre then forces $\nabla w=0$ at every point (as in L7, order 0), so $\mu=0$.

**Corollary M10 (the whole family under (H)).** *Proven.* Assume (H). For every member with $(x,y)\neq0$, in (14,0) and in (7,7):

(a) $\Upsilon=0$ has no solution on any open $U$ that meets $s(X)$.

(b) $\Upsilon=0$ has no solution on any nonempty open $U$ at all, with one possible exception. The exception is the members with $\ell=0$, $N\neq0$ and $(z,w)\neq0$, i.e. up to M1 the members $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$ and, in (7,7) with $N<0$, $y(W'_A+cW_C\mathrm{Ad}_\Gamma)$, with $c\neq0$. For them a solution on a $U$ with $U\cap s(X)=\emptyset$ would need $(c_V/c_H)R_g(x)\in\mathcal Z\setminus\{0\}$ for every $x\in\pi(U)$ (M9).

For $(x,y)=0$ the solutions are exactly $T=0$, for every metric (M8).

*Proof.*
- Each member with $(x,y)\neq0$ falls into one of five classes:
  - $\ell\neq0$, $N\neq0$, $xw\neq yz$: L8(A).
  - $\ell\neq0$, $N\neq0$, $xw=yz$: then $(z,w)=(\ell/N)(x,y)$, which is L8(B) and L9.
  - $\ell=0$, $N\neq0$: M6.
  - $N=0$ (only in (7,7)): M7.
- In (14,0), $N\neq0$ whenever $(x,y)\neq0$. ∎

**Answers.**
1. *Members with $\ell=0$ (task (a)).*
   - The cone $\ell=0$ consists, up to the rotation/boost of M1(a) and the (7,7) swap of M1(b), of the one-parameter family $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, the null members $(x,\pm x,\pm cx,cx)$ of (7,7), and the boundary $(x,y)=0$.
   - On it the curvature-free row vanishes, so the analysis runs through the first-order channel. That channel kills $u$ (M3), not $dt\wedge du$; L1 holds on the cone with $N\neq0$ but is not needed.
   - $\mathrm{span}(W_A,W'_A)$ and the null members: no solution on any open set (M6(b), M7(a)).
   - $c\neq0$: no solution on any open set meeting $s(X)$ (M6(a)). Here the remaining scalar equation has a nonlinear $\mu$, and the fibre obstruction had to be extended to $\mu'=0$ (M5).
2. *$N=0$ in (7,7) (task (b)).* The curvature-free channel degenerates further: both $R'$ and its trace coefficient vanish along $v$. It still forces $\sigma(\kappa+208\sigma)=0$. The only invariant line is $\sigma=-\kappa/208$, where the remaining equation is linear, $\mu=-\tfrac5{12}r$, so L9's continuation applies (M7(b)).
3. *$(x,y)=0$ (task (c)).* $\Upsilon=0$ says nothing about the metric and forces $T=0$ (M8).
4. *No solution was found* on any member with $(x,y)\neq0$; the only solutions found are the trivial ones $T=0$ of the members with $(x,y)=0$ (M8).

## M.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| Inputs (§M0) | L's §L0 run silently (8 asserted checks: J's `j0`/`j1` for $S$, $M$, $c_1$, $c_{13}$; the J2 tables against J's evaluator) | both signatures |
| M1(a) (§M1) | $M$, $S$, the $\ell$ and $b$ tables, $c_i$, $x_i$ under the symbolic one-parameter group $R(m)$, exact rational-function identities | (14,0) $SO(2)$, (7,7) $SO^+(1,1)$ |
| M1(b) (§M1) | the reduced numerators ($\det M$, $H\det M$, $\mathcal Q(\det M)^2$, $\mathcal Q'(\det M)^2$, $\varpi(\det M)^2$) of $p$ and $-Kp$, symbolic in $x,y,z,w,t,u,\kappa$ | (7,7) |
| M2 (§M2) | exact rational orbit maps; cone algebra; $\varpi$ for $(x,0,0,w)$ symbolic; all $\varpi\equiv0$ for $N=0=\ell$ and $(x,y)=0$ | both signatures |
| M3 (§M3) | exact symbolic ($x,w,\kappa,t,u$) vanishing pattern and the form of $E^{(13)}$; mirror $(0,y,z,0)$ | both; mirror (7,7) |
| M4 (§M4) | $\det M$, $H$, $c_1$ on $u=0$; $\mu'(\tau,c)$ in closed form; $\mu'=-\tfrac7{18}$ at $c=0$; discriminant; zero of $H$ | both |
| M5 (§M5) | pairing constants on all 36 basis elements of $\mathfrak{so}\otimes\mathfrak{so}$; Gröbner certificate for (C1)$_0$ with $\mathrm{Ric}=0$ (2796 equations, rank 12, basis $\{1\}$); brute-force Taylor system (1000 unknowns) with a random Weyl tensor, inconsistent; control $R_1$, consistent | pairings (4,0), (1,3), (3,1), (2,2); certificate (1,3), (3,1), (2,2); Taylor (1,3) |
| M6 | proof from M1–M5, L2, L6, L8(B), L9 | (H) |
| M7 (§M6) | exact symbolic ($x,z,w,c,\kappa,t,u$): the (a) identities; $v$ polynomial, $h\cdot v=0$, $R'=0$, $K_{\rm tr}=0$, $yc_1-xc_{13}=\delta(\kappa+208\sigma)$, $v\cdot\nabla\delta$; coefficients on $\sigma=-\kappa/208$, $\mu'=-\tfrac5{12}$, $F$ and $\mu=-\tfrac5{12}r$ | (7,7), $y=\pm x$ |
| M8 (§M7) | exact symbolic ($z,w,\kappa,t,u$): $\det M$, the $\theta$-structure of $\Sigma,H,\mathcal Q,\mathcal Q'$, the derivative-free combination, $F''/F'$, the trace identity | both |
| M9 (§M8) | exact ranks of the relaxed order-$\varepsilon$ system (295 unknowns) at 3 off-section samples per base signature and at the section; the $s$-kernel vector | (1,3), (3,1), (4,0), (2,2) |
| M10 | proof from L8, L9, M6, M7 | (H) |

Run: `cd computation/geometric-unity && python3 M_remaining.py > M_remaining.out 2>&1`. Result: 35 checks, 35 True, each asserted; 270 s on one core, a second run 260 s, with identical output apart from the timing annotations.

## M.3 Limits (layers b, d)

- **Open.**
  - The members $\ell=0$, $N\neq0$, $(z,w)\neq0$ on open sets $U$ with $U\cap s(X)=\emptyset$. By M3 and L6 this is rigorously the scalar problem $\nabla^2r=\mu(r)\mathrm{Ric}^Y$, $Er=0$, with nonlinear $\mu$ (M4), and fibrewise $\mathrm{Hess}^Fw=\mu(w)(\mathrm{Ric}_F+\mathcal A)$.
  - L9's continuation along the fibre uses linearity. For nonlinear $\mu$ the continued solution can reach the singular values of $\mu$ ($H=0$ or $\det M=0$).
  - M9 reduces the open case to base curvatures in the proper algebraic set $\mathcal Z$. $\mathcal Z$ is not determined: that would need the relaxed system with symbolic $R$.
- **Hypotheses kept from (H).** Z2; the real forms (14,0), (7,7) of J1–J2; $\kappa\neq0$; $T_e$ the slaved field in $\mathcal V$ with $\det M\neq0$ on $U$; $T_o$ in the invariant sector. M8 uses L3, hence Z2.
- **M1** is a symmetry of the reduced system of J1–J2. It is not claimed for the full $\Upsilon$.
- **M5(c)** is a certificate for the Ricci-flat case only. Indefinite, non-Ricci-flat solutions of (C1)$_0$ are not classified; M6 does not need them.
- **M9** is a computation at sample points, together with a genericity argument. It is not a proof for every base metric.
- **Untouched:** Z1; reading (A); the (P) system of J3(a); J3(b), (c) and J4(d′) for the members of this note.

## M.4 Consequences for gaps and claims

- **GU-G14.**
  - Under (H), the non-radial Z2 question is now settled for every member with $(x,y)\neq0$ on every open set meeting $s(X)$ (M10(a)).
  - It is settled on every open set except for the twisted members $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, $c\neq0$, with $\ell=0$ (up to M1), away from $s(X)$ (M10(b)). Even there M9 confines a solution to special base curvature.
  - This includes $W'_A$ and all of $\mathrm{span}(W_A,W'_A)$ (M6(b)), $(1,0,0,1)$ in (14,0) near $s(X)$ (M6(a)), and every member with $x=\pm y$ in (7,7) (M7).
  - For $(x,y)=0$, $\Upsilon=0$ forces $T=0$ and leaves the metric free (M8). These members carry no gravitational equation at all in this sector.
  - Open: the exception of M10(b); Z1; other sectors.
- **GU-G16.**
  - $\ell$ is again what separates the cases. For $\ell=0$ the curvature-free row vanishes, the first-order channel removes one scalar (M3), and the remaining equation is nonlinear ($\mu'$ is a nonconstant function of $c=w/x$ and $\tau$) unless $c=0$.
  - The only weaker result in the family sits there.
  - $N=0$ is not special in the end: the line $\sigma=-\kappa/208$ carries a linear equation (M7(b)).
- **GU-054.** Status unchanged (C). No statement made before this note is contradicted. L.3's list of open members and J.3's "Still open there" are answered as above.
