# L. Υ = 0 on a neighbourhood beyond the radial Z2 ansatz: non-radial $t,u$

Claims: GU-054 (context). Gaps: GU-G14 (and GU-G16 for the role of $xz\pm yw$). Code: `computation/geometric-unity/L_nonradial.py` → `.out` (sections L0–L6; exact: Python Fractions, sympy rationals and Gaussian rationals; fixed seeds; 380 s on one core, a second run 380 s, identical apart from the timing annotations). It imports `J_slaved_sector.py` (and through it `I_elliptic_reduction.py`, `G_symbols.py`, …) and `F_jets.py`, both unmodified. Output sections: §L0 inputs, §L1 Lemma L1, §L2 Lemma L2, §L3 Lemma L3, §L4 Proposition L5, §L5 Theorem L8 step 5, §L6 Lemma L7. Record: `metrics/geometric-unity/L-nonradial-z2.yaml`. Builds on J (J1, J2, Cor. J2.1, J5), J5d (L0, L1), F (F1, F3, L-F1, L-F5) and A (A4, A5). Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## L.0 Question and setting (layer b)

J5 and J5d exclude $\Upsilon=0$ on a neighbourhood of a point of $s(X)$ only for the radial ansatz $t=f(s)$, $u=h(s)$ in Z2. J.3 lists as open the non-radial case, where one only knew $dt\wedge du=0$ at points with some $\varpi\neq0$ (Cor. J2.1(e)). This note removes the radial ansatz.

**Hypotheses (H)** (those of J5, minus the radial ansatz):
- Z2 with $n=4$, $\lambda\neq\tfrac14$, $c_Hc_V\neq0$; the two real forms of J.0: (14,0) (base (4,0), $\lambda<\tfrac14$) and (7,7) (base (1,3), (3,1) or (2,2), $\lambda>\tfrac14$; A4, A5).
- $\kappa\neq0$; $T_o=t\Phi^1+u\Phi^1\Gamma$ with $t,u$ arbitrary smooth functions on $U$; $T_e=\Sigma\,(dt,du)^\top$ the slaved field of J1 in $\mathcal V$, with $\det M(t,u)\neq0$ on $U$.
- $U\subset Y$ connected and open; $\Upsilon=0$ on $U$. Where stated, $U$ meets $s(X)$.
- $(x,y)\neq0$ (as in J5); every member of L8 satisfies it.

**Notation.** $E^{(1)},E^{(13)},\omega^{(3)},\omega^{(11)}$, $H$, $\mathcal Q,\mathcal Q'$, $\varpi$ are those of J2 ("Substituting the slaving"); $x_1=x$, $x_{13}=y$; $\epsilon=\pm1$ in (14,0)/(7,7). Put
$$\ell:=xz+\epsilon yw\ (\text{the ellipticity factor of H7}),\qquad N:=x^2+\epsilon y^2,\qquad \delta:=ty-ux,\qquad U^*:=\{q\in U:(dt,du)(q)\neq0\}.$$
$s=\tfrac14\log|\det(g^{-1}y)|$ and the Euler field $E=y\partial_y$ are those of J5(a): $E$ is Killing and parallel in Z2, $E^\flat=4c_V(1-4\lambda)\,ds$, $|E|^2=4c_V(1-4\lambda)\neq0$, $\nabla ds=0$. $V_0\subset V$ is the $g_\aleph$-orthogonal complement of $E$ (traceless directions).

**Real forms.** In (7,7) everything is real. In (14,0), $(x,y,z,w)=-i(x',y',z',w')$ and $t,u\in i\mathbb R$ (I.0). Every reduced equation is a real tensor equation times a constant phase, and $\delta$, $N\,\sigma$, … are real multiples of real quantities. So below $t,u$ are treated as real functions and all algebraic statements are polynomial identities, hence hold in both real forms. Coefficient functions such as $\alpha,\beta$ in "$\mathrm{Ric}=\alpha g+\beta\,\theta\otimes\theta$" are then real: their imaginary parts multiply real tensors of different rank.

**Inputs (all proven earlier).**
- F1: the fibre $F_x=(\mathrm{Met}_{(p,q)}(T_xX),c_VG_\lambda)$ is $\mathbb R\times S$, $S=SL(4)/SO(p,q)$, with $\partial_s$ parallel; $\mathrm{Ric}_F=-G_{1/4}$, i.e. $0$ on $\partial_s$ and $-\tfrac1{c_V}g$ on $V_0$; $R^S(X,Y)Z=-\tfrac14[[X,Y],Z]$ in the endomorphism picture $H=y^{-1}h$.
- F3: $\pi:Y\to(X,c_Hg)$ is a pseudo-Riemannian submersion with totally geodesic fibres ($T=0$) and $A_XY=-\tfrac12(R_g(X,Y)\cdot y)^{\rm vert}$ (L-F5), which vanishes on $s(X)$.
- O'Neill (L-F1, $T=0$): at every point of $Y$,
$$\mathrm{Ric}^Y(U,W)=\mathrm{Ric}_F(U,W)+\mathcal A(U,W),\qquad \mathrm{Ric}^Y(X,X')=\mathrm{Ric}_g(X,X')-2(A_X,A_{X'}),$$
  with $\mathcal A(U,W)=\sum_i\epsilon_ig_\aleph(A_{E_i}U,A_{E_i}W)=\tfrac{c_V^2}{4c_H^2}\sum_{ij}\epsilon_i\epsilon_j\,G_y(U,R_{ij}\cdot y)\,G_y(W,R_{ij}\cdot y)$ and $(A_X,A_{X'})=\tfrac{c_V}{4c_H}\Phi_y(X,X')$ (J5d L0(a)); $R_{ij}=R_g(e_i,e_j)$. Both blocks are checked by brute force **off** the section in output §L2.
- J5d: $\det(y)^2\Phi_y$ is a polynomial in $y$, $\Phi_g=0$, and $\Phi\equiv0$ on an open set of $y$ forces $F\equiv0$ (L0(b), (c)), hence $|\mathrm{Ric}_g(x_0)|^2=0$ (L1(a), $n=4$).
- Cor. J2.1(e): the antisymmetric parts of $E^{(1)},E^{(13)}$ are $\varpi^{(1)},\varpi^{(13)}$ times $dt\wedge du$, and $\omega^{(j)}=\varpi^{(j)}dt\wedge du$.

## L.1 Results (layer c)

**Lemma L1 (the Λ² equations hold everywhere).** Assume (H). If some $\varpi^{(j)}(t,u)$ is not the zero rational function, then $dt\wedge du=0$ on all of $U$. This holds:
- for $\circledast_\sharp$ ($\varpi^{(13)}_\sharp$, the exact formula of Cor. J2.1(e));
- for $xW_A$, $x\neq0$ ($\varpi^{(1)}=\pm1916928\,\kappa u x^4/(\kappa^2+224\kappa tx\mp3328u^2x^2)^2$);
- for $xW_A+zW_C\mathrm{Ad}_\Gamma$ with $xz\neq0$ ($\varpi^{(13)}=xz^2\cdot(\text{nonzero})/(\det M)^2$);
- for every member with $\ell\neq0$ (at leading order $\varpi^{(3)}=\mp88320\,i\,y\ell/\kappa^2$, $\varpi^{(11)}=88320\,i\,x\ell/\kappa^2$, and $(x,y)\neq0$ since $\ell\neq0$).

*Proof.* On the open set $O=\{dt\wedge du\neq0\}$ the map $(t,u)$ is a submersion, so its image contains an open set $\Omega$ of the $(t,u)$-plane. At each point of $O$, $\Upsilon=0$ forces $\varpi^{(j)}(t,u)=0$. So $\varpi^{(j)}$ vanishes on $\Omega$, hence identically, a contradiction. So $O=\emptyset$. ∎ (Output §L1.)

So the restriction "at points where some $\varpi\neq0$" of J.3 is not needed on open solution sets.

**Lemma L2 (Z2 has no Ricci-special open sets).** Let $W\subset Y$ be nonempty and open, and $\alpha$ a function on $W$ with

(i) $\mathrm{Ric}^Y(U,U')=\alpha\,g_\aleph(U,U')$ for $U,U'\in V_0$, and (ii) $\mathrm{tr}_H\mathrm{Ric}^Y=4\alpha$.

Then for every $x_0\in\pi(W)$: $|A|^2\equiv0$ on the whole fibre over $x_0$, $\alpha=-1/c_V$ on $W\cap\pi^{-1}(x_0)$, and the scalar curvature is $\mathrm{scal}_g(x_0)=-4c_H/c_V$. If moreover

(iii) $\mathrm{Ric}^Y(X,X')=\alpha\,g_\aleph(X,X')$ for horizontal $X,X'$,

there is no such $W$; in (14,0) already (i)+(ii) are impossible. In particular there is no nonempty open $W$ on which
- (a) $\mathrm{Ric}^Y=\alpha g_\aleph+\beta\,ds\otimes ds$ for functions $\alpha,\beta$ (this includes Einstein and Ricci-flat), or
- (b) $\mathrm{Ric}^Y=\beta\,\theta\otimes\theta$ with $\theta$ horizontal and null.

*Proof.*
1. By O'Neill, the trace of (i) over $V_0$ (dimension 9) is $9\alpha=-\tfrac9{c_V}+|A|^2$. Here $\mathcal A(E,\cdot)=0$, because $g(A_XE,X')=-g(E,A_XX')=0$ (the vector $R\cdot y$ is $y$-traceless). So $\sum_{V_0}\epsilon_U\mathcal A(U,U)=\sum_{ij}\epsilon_i\epsilon_j|A_{E_i}E_j|^2=|A|^2$.
2. The trace of (ii) is $4\alpha=\mathrm{scal}_g/c_H-2|A|^2$. Together, $\alpha=\tfrac1{11}\big(\tfrac{\mathrm{scal}_g}{2c_H}-\tfrac9{c_V}\big)$, which depends on $x_0$ only. Hence $|A|^2=9\alpha+9/c_V$ is constant on the open set $W\cap\pi^{-1}(x_0)$ of the fibre.
3. $|A|^2(x_0,y)=\tfrac{c_V}{4c_H^2}\sum_i\epsilon_i\Phi_y(e_i,e_i)$, and $\det(y)^2(|A|^2-C)$ is a polynomial in $y$ that vanishes on an open subset of $\mathrm{Sym}^2$. So it vanishes identically, and at $y=g(x_0)$ ($A=0$, F3) one gets $C=0$. This gives $\alpha=-1/c_V$ and $\mathrm{scal}_g(x_0)=c_H(22\alpha+18/c_V)=-4c_H/c_V$.
4. In (14,0), $g_\aleph$ is definite, so $|A|^2\equiv0$ gives $A\equiv0$ on the fibre, i.e. $R_g(x_0)\cdot y=0$ for all $y$, so $R_g(x_0)=0$ (as in J5(c)). This contradicts $\mathrm{scal}_g(x_0)=-4c_H/c_V$.
5. With (iii): $(A_X,A_{X'})=\tfrac12(\mathrm{Ric}_g-\alpha c_Hg)(X,X')$ is independent of $y$. By the same polynomial argument it equals its value at $y=g$, which is 0. So $\Phi\equiv0$ on the fibre and $\mathrm{Ric}_g(x_0)=\alpha c_Hg=-\tfrac{c_H}{c_V}g$. J5d L0(c) and L1(a) give $|\mathrm{Ric}_g(x_0)|^2=0$, against $4c_H^2/c_V^2$.
6. (a): $ds$ vanishes on $H$ and on $V_0$, so (i)–(iii) hold. (b): (i) and (ii) hold with $\alpha=0$, against $\alpha=-1/c_V$. ∎

Inputs checked in output §L2 (brute-force Ricci of Z2 at points off $s(X)$, bases (4,0), (1,3), (2,2)).

*Remark (side results).*
- (a) with $\beta=0$: for $n=4$ and every $\lambda\neq\tfrac14$, **Z2 is Einstein on no open subset of $Y$**. F6 had this only on $s(X)$. This proves Conjecture F9 for Z2 at $n=4$, in a stronger form (flat bases included); its open content is Z1, where the flat model F7 is Einstein (R8-4).
- (a) with $\beta\neq0$ is the conclusion J5(b) draws from the radial ansatz. So **J5(c) holds without the hypothesis that $U$ meets $s(X)$**: no radial solution exists on any open set.

**Lemma L3 (fibre lemma: concircular functions on Z2 are parallel).**
- (a) Let $\Omega$ be a connected open subset of a fibre $F_x$, and let $f$ on $\Omega$ satisfy $\mathrm{Hess}^F f=\psi\,g_F$. Then $\psi=0$ and $f=as+b$ with constants $a,b$.
- (b) If $r$ on an open $W\subset Y$ is concircular, $\mathrm{Hess}^Yr=\psi\,g_\aleph$, then $\psi\equiv0$, i.e. $\nabla dr=0$. Moreover $a:=g(\nabla r,\nabla s)/|ds|^2$ is constant on each component, and $Z:=\nabla r-a\nabla s$ is parallel and horizontal.

*Proof.*
- (a) For $\nabla V=\psi\,\mathrm{Id}$ ($V=\nabla f$) one has $R(X,Y)V=(X\psi)Y-(Y\psi)X$. Since $\partial_s$ is parallel, $R(\partial_s,Y)=0$. Taking $Y\perp\partial_s$ gives $d\psi=0$, so $R(X,Y)V=0$ for all $X,Y$.
- With $V=V_s\partial_s+V_S$: $R^S(X,Y)V_S=-\tfrac14[[X,Y],V_S]=0$ for all $X,Y\in T S$. Since $[\mathfrak p,\mathfrak p]=\mathfrak{so}(\eta)$, $V_S$ commutes with $\mathfrak{so}(\eta)$, whose commutant in $\mathfrak{gl}(4)$ is $\mathbb R1$ ($\mathbb R^4$ is absolutely irreducible: its complexification is the vector representation of $\mathfrak{so}(4,\mathbb C)$; output §L3). So $V_S$ is a multiple of $1$, and it is traceless, so $V_S=0$.
- Then $\nabla_Y(V_s\partial_s)=(YV_s)\partial_s=\psi Y$; for $Y\in TS$ this forces $\psi=0$ and $V_s$ constant.
- (b) The fibres are totally geodesic ($T=0$), so $\mathrm{Hess}^Yr|_{VV}=\mathrm{Hess}^F(r|_F)$ and (a) applies on every fibre piece: $\psi=0$, and the vertical part of $\nabla r$ is $a\nabla s$ ($\nabla s$ is vertical). $da=0$ because $\nabla dr=\nabla ds=0$. ∎

Output §L3: the commutant and $[\mathfrak p,\mathfrak p]$ in bases (4,0), (1,3), (2,2), and a Taylor-jet check of $\mathrm{Hess}^Fw=\psi G_\lambda$ at $y=\eta$: $\psi(\eta)=0$, $d\psi(\eta)=0$, $dw(\eta)\in\mathbb R\,ds$.

So the warped-product structure suggested by the curvature-free channel (J.3) cannot occur on Z2: the concircular function is forced to have parallel gradient.

**Proposition L4 (the concircular branch is empty).** Assume (H) and let $W\subset U$ be open. There is no function $r$ on $W$ with $dr\neq0$, $\nabla^2r=\psi g_\aleph$ and $\mathrm{Ric}^Y\in\mathrm{span}(g_\aleph,dr\otimes dr)$ pointwise.

*Proof.*
- By L3(b), $\nabla dr=0$, $a$ is constant and $Z$ is parallel and horizontal. Write $\mathrm{Ric}=\alpha g+\beta\,dr\otimes dr$.
- Parallel fields are in the kernel of Ric. So $0=\mathrm{Ric}(\nabla s)=\alpha\,ds+\beta a|ds|^2dr$ and $0=\mathrm{Ric}(Z)=\alpha Z^\flat+\beta|Z|^2dr$.
- *$Z=0$:* $a\neq0$ since $dr\neq0$, and $\mathrm{Ric}=\alpha g+\beta a^2ds\otimes ds$. This is L2(a): the radial case, i.e. J5.
- *$Z\neq0$, $a\neq0$:* $ds$ and $Z^\flat$ are independent, so $\beta=0=\alpha$ and $\mathrm{Ric}=0$, which is L2(a).
- *$a=0$:* $\alpha=0$ and $\beta|Z|^2=0$. Either $\mathrm{Ric}=0$ (L2(a)) or $Z$ is null and $\mathrm{Ric}=\beta Z^\flat\otimes Z^\flat$ (L2(b)). This covers the null-gradient case (only possible in (7,7)) and the "basic" case $r=r(x)$.
- In all cases L2 gives a contradiction. ∎

**Proposition L5 (the curvature-free channel).** For every member, $E':=yE^{(1)}-xE^{(13)}$ contains no curvature term. Its Hessian row is $h:=(y,-x)H$, and its gradient coefficients are $\mathcal Q'':=y\mathcal Q^{(1)}-x\mathcal Q^{(13)}$.

(a) *Kernel field and residual* (output §L4, exact, symbolic in $x,y,z,w,t,u$ at $\kappa=1$). For $\ell\neq0$,
$$v:=\frac{\det M}{576\,\ell}\,(h_u,-h_t)$$
is polynomial, and nonzero wherever $\det M\neq0$ (because $\det H=576^2\ell^2\kappa^2/\det M\neq0$, Cor. J2.1(b′)). Moreover, exactly,
$$R'(t,u):=h\cdot(Dv)v+v^\top\mathcal Q''_{\rm sym}v=-1916928\,N\,\ell^2\,\delta .$$
Since Σ is homogeneous of degree $-1$ in $(\kappa,t,u)$, $R'$ is homogeneous of degree 4, so for general κ it equals $-1916928\,\kappa^3N\ell^2\delta$; its zero set is the line $\delta=0$, i.e. $(t,u)\parallel(x,y)$.

(b) *Invariance of the line.* At $\kappa=1$, on $(t,u)=\sigma(x,y)$: $\;v\cdot\nabla\delta=-(xw-yz)(160N\sigma-1)$. So the line is $v$-invariant iff $xw=yz$, i.e. $(z,w)=c(x,y)$. Otherwise $v$ is tangent to it only at one point ($\sigma_*=\kappa/(160N)$ for general κ, by homogeneity).

(c) *The family $xW_A+zW_C\mathrm{Ad}_\Gamma$* (κ symbolic). $H_{13}=z^2\cdot(\ldots)/\det M$. For $z\neq0$ the kernel field can be taken as $v=(\kappa^2+224\kappa tx\mp3328u^2x^2,\ 384\kappa ux)$ with residual $-1916928\,\kappa ux^4$ of the $E^{(13)}$ row (in the normalisation of (a), $v=-\kappa z\,v_c$ and $R'=-x(\kappa z)^2$ times it); $v$ is tangent to $u=0$, and $H_{13}=0$ only where $\det M=0$. For $z=0$ (the members $xW_A$), $H_{13}\equiv0$, Σ has no $du$ column, and
$$E^{(13)}=c_{13}\,g-16x\,\Sigma_C\Sigma_B\,(dt\otimes dt-|dt|^2g),\qquad \Sigma_C=\frac{24\kappa x}{\kappa^2+224\kappa tx\mp3328u^2x^2},\quad \Sigma_B=\frac{4992\,ux^2}{\kappa^2+224\kappa tx\mp3328u^2x^2}.$$

(d) *Meaning.* On $U^*$, by L1, $(t,u)=\gamma(\rho)$ locally with $\gamma'\neq0$, $d\rho\neq0$. Then the traceless part of $E'$ is $\hat h\,(\nabla^2\rho)_0+\hat q\,(d\rho\otimes d\rho)_0$, with $\hat h=h\cdot\gamma'$ and $\hat q=h\cdot\gamma''+\gamma'^\top\mathcal Q''\gamma'$. If $\hat h\equiv0$ along γ, then $\gamma'=\mu v$ with $\mu\neq0$ and $\hat q=\mu^2R'(\gamma)$. (Checked at random exact points and curve jets, member $(2,-1,3,1)$, both signatures, output §L4.)

*Proof.* (a)–(c) are exact polynomial identities (output §L4). (d): $\nabla^2\phi=\gamma'_\phi\nabla^2\rho+\gamma''_\phi d\rho\otimes d\rho$. If $h\cdot\gamma'=0$ then $\gamma'\parallel v$. Differentiating $h(\gamma)\cdot v(\gamma)=0$ gives $h\cdot\gamma''=\mu^2h\cdot(Dv)v$. ∎

**Lemma L6 (the residual equation).** Let $W\subset Y$ be open, $r$ a function on $W$, and μ a smooth function of one variable with $\mu'(r)\neq0$ on $W$. Suppose
$$\nabla^2r=\mu(r)\,\mathrm{Ric}^Y+\nu\,g_\aleph \qquad (\star)$$
for some function ν. Then $Er=0$ and $\nu=0$ on $W$.

*Proof.*
- $E$ is parallel, so $\nabla^2r(E,\cdot)=d(Er)$ and $\mathrm{Ric}(E,\cdot)=0$. Hence $dk=\nu E^\flat=c\,\nu\,ds$, with $k:=Er$ and $c=4c_V(1-4\lambda)\neq0$.
- So $d\nu\wedge ds=0$ and $\nabla^2k=c\,d\nu\otimes ds=c\nu_s\,ds\otimes ds$.
- $E$ is Killing, so $\mathcal L_E$ commutes with $\nabla^2$ and kills Ric and $g$. Applying $\mathcal L_E$ to $(\star)$: $c\nu_s\,ds\otimes ds=\mu'(r)k\,\mathrm{Ric}+(E\nu)g$.
- On the open set $\{\mu'(r)k\neq0\}$ this puts Ric in $\mathrm{span}(g,ds\otimes ds)$, which L2(a) excludes. So $k\equiv0$, and then $\nu=0$ from $dk=c\nu\,ds$. ∎

**Lemma L7 (fibre obstruction at a point of $s(X)$).** Let $p=g(x_0)\in s(X)$ and Ω a neighbourhood of $p$ in the fibre $F_{x_0}$. Let μ be smooth with $\mu(w(p))\neq0\neq\mu'(w(p))$. There is no $w$ on Ω with
$$\mathrm{Hess}^Fw=\mu(w)\,(\mathrm{Ric}_F+\mathcal A_{x_0}) .$$
If $\mu(w)=\mu_1w$ with $\mu_1\neq0$ (linear case), every solution near $p$ vanishes identically.

*Proof.* Write $\hat{\mathcal R}=g_F^{-1}(\mathrm{Ric}_F+\mathcal A)$ and $\hat{\mathcal A}=g_F^{-1}\mathcal A$, so that $\nabla_Y\nabla w=\mu(w)\hat{\mathcal R}Y$. The Ricci identity gives
$$R^F(Y,Z)\nabla w=\mu'(w)\big[(Yw)\hat{\mathcal R}Z-(Zw)\hat{\mathcal R}Y\big]+\mu(w)\big[(\nabla_Y\hat{\mathcal A})Z-(\nabla_Z\hat{\mathcal A})Y\big],$$
using $\nabla\mathrm{Ric}_F=0$ (F is symmetric). At $p$: $\mathcal A=O(\varepsilon^2)$ along the fibre ($A=0$ on $s(X)$ and $A$ is linear in $R\cdot y$), so $\mathcal A(p)=0$ and $\nabla\mathcal A(p)=0$.
- *Order 0.* Take $Y=\partial_s$ (with $R(\partial_s,\cdot)=0$ and $\hat{\mathcal R}\partial_s=0$) and $Z\in TS$, where $\hat{\mathcal R}Z=-Z/c_V$: this gives $\partial_sw(p)=0$. For $Y,Z\in TS$: $-\tfrac14[[Y,Z],v]=-\mu'[\langle v,Y\rangle Z-\langle v,Z\rangle Y]$ with $\langle\cdot,\cdot\rangle=G_0$. Choosing commuting independent $Y,Z$ (e.g. $Y$ generic, $Z=Y^2-\tfrac14\mathrm{tr}Y^2$) gives $\langle v,Y\rangle=0$ for generic $Y$, so $v=0$. Hence $\nabla w(p)=0$.
- *Order 1.* Differentiate the identity at $p$ along $X\in T_pS$, using $\nabla w(p)=0$ and $\nabla_X\nabla w=-\tfrac{\mu_0}{c_V}X$ ($\mu_0=\mu(w(p))$, $\mu'=\mu'(w(p))$). After division by $-\mu_0/c_V$, for all $X,Y,Z\in\mathfrak p$:
$$\mathrm{(C1)}\qquad c_0(X,Y,Z)+\mu'\,c_1(X,Y,Z)+q(\mathfrak B)(X,Y,Z)=0,$$
  where
  - $c_0=-\tfrac14[[Y,Z],X]$ and $c_1=\mathrm{tr}(XY)Z-\mathrm{tr}(XZ)Y$;
  - $q(\mathfrak B)(X,Y,Z)=d^2\!\mathcal A(X,Y)Z-d^2\!\mathcal A(X,Z)Y$, with $d^2\!\mathcal A(X,Y)Z=\tfrac14\sum_s[\mathrm{tr}(Z[K_s,X])[K'_s,Y]+\mathrm{tr}(Z[K_s,Y])[K'_s,X]]$;
  - $\mathfrak B=\sum_sK_s\otimes K'_s=\big(\tfrac{c_V}{c_H}\big)^2\sum_{ij}\epsilon_i\epsilon_j R_{ij}\otimes R_{ij}\in\mathfrak{so}(\eta)\otimes\mathfrak{so}(\eta)$.
- *Pairing.* With $\langle E,F\rangle=\sum_{abc}\mathrm{tr}\,E(X_a,X_b,X_c)F(X^a,X^b,X^c)$ over dual bases of $\mathfrak p$, exactly (output §L6, all four base signatures):
  - $\langle c_0,c_0\rangle=\tfrac{27}2$, $\langle c_0,c_1\rangle=18$, $\langle c_1,c_1\rangle=144$;
  - for every $K\otimes K'$: $\langle q(K\otimes K'),c_1\rangle=-9\,\mathrm{tr}(KK')$ and $\langle q(K\otimes K'),c_0\rangle=-\tfrac{27}4\mathrm{tr}(KK')$.
- Hence with $T:=4c_0-3c_1$: $\langle q(\mathfrak B),T\rangle=0$ for **every** $\mathfrak B$, $\langle c_0,T\rangle=0$ and $\langle c_1,T\rangle=-360$. Pairing (C1) with $T$ gives $-360\,\mu'=0$, a contradiction.
- *Linear case.* Order 0 gives $\nabla w(p)=0$ without using $w(p)$. If $w(p)\neq0$, then $\mu_0\neq0$ and order 1 contradicts as above. So $w(p)=0$ and $(w,\nabla w)(p)=0$, and a solution of this finite-type linear equation vanishes on a neighbourhood. ∎

Two independent confirmations (output §L6):
- the span certificate: for the generic curvature tensor (20 components, 2988 scalar equations) the quadratic parts of (C1) have rank 21, $c_0$ lies in their span and $c_1$ does not, so (C1) is solvable only for $\mu'=0$;
- a brute-force Taylor solution of the fibre equation at $y=\eta$ to order $\varepsilon^2$ (1000 unknowns, Christoffel symbols from $G_\lambda$ itself):
  - it is inconsistent for $\mu'=-\tfrac7{18}$ with random $R$ in bases (1,3), (4,0), (2,2);
  - it is consistent exactly in the control case that (C1) predicts to be solvable ($\mu'=0$, $R$ the unit constant-curvature tensor, A-term multiplied by $-\tfrac1{12}$), and inconsistent again in the same control with $\mu'=-\tfrac7{18}$;
  - at order ε it forces $dw(\eta)=0$.

**Theorem L8 (non-radial Z2).** Assume (H).

(A) *Members with $\ell\neq0$ (in (14,0): the elliptic members), $N\neq0$ and $xw\neq yz$.* If $\ell\neq0$, $N\neq0$ and $xw\neq yz$, then $\Upsilon=0$ has **no solution on $U$**, for every nonempty connected open $U$ (it need not meet $s(X)$).

(B) *The members $x(W_A+cW_C\mathrm{Ad}_\Gamma)+y(W'_A+cW'_C\mathrm{Ad}_\Gamma)$ with $cN\neq0$, and $xW_A$ with $x\neq0$.* This includes $\circledast_\sharp$ ($x=c=1$, $y=0$) and $W_A$. If $U$ meets $s(X)$, $\Upsilon=0$ has **no solution on $U$**.

In (14,0), $N=-(x'^2+y'^2)\neq0$ for every member with $(x,y)\neq0$. So (A) and (B) together cover **every elliptic member** ($x'z'+y'w'\neq0$, H7) and the multiples of $W_A$. In (7,7) they cover every member with $xz-yw\neq0$ and $x\neq\pm y$, and the multiples of $W_A$.

*Proof.*
1. *Λ² step.* $dt\wedge du=0$ on $U$ (L1).
2. *No locally constant pieces.* On an open set where $dt=du=0$, the field $T_e=\Sigma(dt,du)$ vanishes. The odd equation is then I6's, $-x_iG+c_ig=0$, so $\mathrm{Ric}\propto g$, which L2(a) excludes. Hence $U^*$ is dense in $U$.
3. *Concircular branch.* Let $q\in U^*$ with $\hat h=h\cdot\gamma'\neq0$ (well defined up to a factor). Near $q$, take $r=F(\rho)$ with $F''/F'=\hat q/\hat h$ (real; see L.0). By L5(d), $(\nabla^2r)_0=0$: $r$ is concircular. Then every term of $E^{(i)}$ lies in $\mathrm{span}(g,dr\otimes dr)$, and some $x_i\neq0$, so $\mathrm{Ric}\in\mathrm{span}(g,dr\otimes dr)$. This is excluded by L4. So $h\cdot(dt,du)=0$ on $U^*$.
   - For $xW_A$, $h\equiv0$ and this step is empty. There, L5(c) gives $u\,dt=0$ pointwise, because $(dt\otimes dt)_0=0$ forces $dt=0$ in dimension 14.
   - On $\{u\neq0\}$ one then has $dt=0$, hence $T_e=0$ (no $du$ column) and $\mathrm{Ric}\propto g$, excluded by L2(a). So $u\equiv0$.
4. *Residual.* For $\ell\neq0$: on $U^*$, $\gamma'=\mu v$ and the traceless part of $E'$ is $\mu^2R'(t,u)\,(d\rho\otimes d\rho)_0=0$, so $\delta(t,u)=0$ on $U^*$ (L5(a); $\kappa N\ell\neq0$).
   - (A): along γ the curve lies in the line $\delta=0$ and is tangent to $v$, so $v$ is tangent to the line along a whole interval. By L5(b) this happens only at the single point $\sigma_*$, so γ is constant, against $\gamma'\neq0$. Hence $U^*=\emptyset$, and step 2 gives a contradiction.
   - (B): $\delta=0$ on $U^*$, hence on $U$ by density. So $(t,u)=\sigma(x,y)$ with $\sigma=(xt+\epsilon yu)/N$ smooth.
5. *The single remaining equation* (B; output §L5). On the line, $\det M=\kappa^2(\kappa+224N\sigma)(\kappa+544N\sigma)\neq0$, and each $E^{(i)}$ equals $x_i$ times one tensor equation. Its Hessian coefficient is $576N/(\kappa+224N\sigma)$, its gradient coefficients are $-258048N^2/(\ldots)^2$ and $78336N^2/(\ldots)^2$, and its zeroth-order term is $(\kappa\sigma+312N\sigma^2)g$, all independent of $c$. For $xW_A$ the same holds on $u=0$ with $N=x^2$; there $\circledast_\sharp$ and $W_A$ have identical coefficients.
   - With $r=F(\sigma)=-1/(224N(\kappa+224N\sigma))$, the traceless part becomes $(\nabla^2r)_0=\mu(r)(\mathrm{Ric})_0$ with $\mu(r)=-\tfrac7{18}r$ exactly. So $(\star)$ holds with $\mu'=-\tfrac7{18}$, and $r\neq0$ everywhere.
6. By L6, $Er=0$ and $\nu=0$. At $p\in U\cap s(X)$ the fibre is totally geodesic, so $\mathrm{Hess}^F(r|_F)=-\tfrac7{18}\,r\,(\mathrm{Ric}_F+\mathcal A_{x_0})$ on a neighbourhood of $p$ in $F_{x_0}$ (O'Neill). This contradicts L7. ∎

**Corollary L9 (no open set at all).** In case (B) the hypothesis "$U$ meets $s(X)$" can be dropped: $\Upsilon=0$ has no solution on any nonempty open $U$.

*Proof.*
- Steps 1–6 up to L6 do not use $s(X)$. Fix $x_0\in\pi(U)$. On the connected open set $\Omega\subset F_{x_0}$, $w=r|_F$ solves the **linear** equation $\mathrm{Hess}^Fw=-\tfrac7{18}w(\mathrm{Ric}_F+\mathcal A_{x_0})$.
- Its coefficients are rational in $y$, with poles only on $\det y=0$, hence analytic on the connected fibre $\mathrm{Met}_{(p,q)}(T_{x_0}X)$ ($=GL^+(4)/SO(p,q)$, since $O(p,q)$ meets both components of $GL(4)$). Solutions are the parallel sections of an analytic connection $D$ on $\mathbb R\oplus TF$.
- Transport $(w,\nabla w)$ along a polygonal path from Ω to $p=g(x_0)$. Each $(\nabla^kR^D)(\cdot)\,(w,\nabla w)$ is analytic along each segment and vanishes near its start, so it vanishes along the segment.
- At each vertex the infinitesimal holonomy algebra annihilates $(w,\nabla w)$. For analytic connections it equals the local holonomy algebra (Kobayashi–Nomizu I, Thm II.10.8: for an analytic connection the infinitesimal holonomy algebra equals the local, indeed the restricted, holonomy algebra), so the transported section is a local solution near the next vertex, and finally near $p$.
- By L7 (linear case) it vanishes there. Transporting back, $w\equiv0$ on Ω, against $r\neq0$. ∎

**Answers.**
1. *Λ² equations.* For the members of L1 (in particular $\circledast_\sharp$, $xW_A$ and every member with $xz\pm yw\neq0$), $dt\wedge du=0$ holds on every open solution set, not only where some $\varpi\neq0$ (L1). (Corrected after review 2026-10-30, R8-1; the text had no condition on the member.) Locally $(t,u)=\gamma(\rho)$ on the dense set $U^*$. The complement of $U^*$ has empty interior, because constant $(t,u)$ on an open set means round 5 with $T_e=0$, i.e. an Einstein open set, which Z2 does not have (L2).
2. *The two Einstein-type channels.*
   - The curvature-free channel ($E^{(13)}$ for $y=0$; $yE^{(1)}-xE^{(13)}$ in general) does make the scalar concircular after reparametrisation, as the referee noted, wherever its reduced Hessian coefficient $\hat h\neq0$.
   - But on Z2 a concircular function has parallel gradient (L3), because the fibre is $\mathbb R\times SL(4)/SO$ with no warped factor. With channel $E^{(1)}$, $\mathrm{Ric}\in\mathrm{span}(g,dr\otimes dr)$ then reduces to Ricci-flat, radial (J5) or null-horizontal cases, all excluded by L2 on every open set.
   - Where $\hat h=0$ the channel forces $R'=0$, i.e. $(t,u)\parallel(x,y)$ (L5). That ends the analysis unless the line is invariant ($xw=yz$).
   - On the invariant line one equation is left, $\nabla^2r=-\tfrac7{18}r\,\mathrm{Ric}^Y$ after L6. L6 uses the Killing, parallel Euler field. The equation fails at order $\varepsilon^2$ in the fibre at any point of $s(X)$ (L7, using $A=0$ on $s(X)$ and the A-tensor formula L-F5).
   - Degenerate cases: a null gradient is covered (L3 and L2(b)); $\nabla t$ parallel is L4's case analysis (it contains J5); constant $(t,u)$ on an open set is excluded by L2(a).
3. *Decided:* $\circledast_\sharp$, $W_A$, and every member in the classes (A), (B) of L8: no solution on a neighbourhood of a point of $s(X)$, and in fact on no open set (L8(A), L9). No further ansatz was needed. What remains open is listed in L.3.
4. *Other members.* L8(A) needs no point of $s(X)$ and no fibre obstruction: for generic elliptic members the curvature-free channel, together with L2–L4 for the concircular branch, already excludes non-constant scalars.

## L.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| Inputs (§L0) | J's own `j0`, `j1` give $S$, $M$, $c_1$, $c_{13}$ equal to the transcription; J's evaluator reproduces the transcribed $\ell$ (basis tensors $(0,1),(2,2),(3,9)$, all $k$) and $b$ (pairs $(e^0,e^1),(e^0,e^0)$, all ten $(k,l)$), four operators. The full tables are J2's (all bases) and xcheck_round6's | both signatures |
| L1 (§L1) | exact $\varpi$ for $\circledast_\sharp$, $W_A$ (equal to J's formulas), $xW_A$, the family $xW_A+zW_C\mathrm{Ad}_\Gamma$; leading $\varpi^{(3)},\varpi^{(11)}$ for the symbolic member | both signatures |
| L2 (§L2) | proof; O'Neill blocks VV $=\mathrm{Ric}_F+\mathcal A$ and HH $=\mathrm{Ric}_g-\tfrac{c_V}{2c_H}\Phi_y$ by brute-force Ricci of Z2 at a point **off** $s(X)$ (`F_jets`, 2-jets in 14 variables, polynomial base metric of degree 3 in general coordinates, $c_H=2$, $c_V=3$); $\mathrm{Ric}(E,\cdot)=0$, $|E|^2=4c_V(1-4\lambda)$; polynomial continuation and L1 of J5d | bases (4,0) $\lambda=0$, (1,3) and (2,2) $\lambda=\tfrac12$ |
| L3 (§L3) | proof; $[\mathfrak p,\mathfrak p]=\mathfrak{so}(\eta)$ and trivial commutant (exact linear algebra); Taylor-jet solution space of $\mathrm{Hess}^Fw=\psi G_\lambda$ at $y=\eta$ to order 2 (1066 unknowns): $\psi,d\psi$ vanish, $dw\in\mathbb R\,ds$ | (4,0), (1,3), (2,2); jets: (4,0), (1,3) |
| L4 | proof (case analysis on L2, L3) | both signatures |
| L5 (§L4) | exact polynomial identities in $x,y,z,w,t,u$ ($\kappa=1$; homogeneity): $v$ polynomial, $R'=-1916928N\ell^2\delta$, $v\cdot\nabla\delta$ on the line; the family with κ symbolic; $xW_A$'s $E^{(13)}$ and Σ; $\hat q=\mu^2R'$ at random exact curve jets | both signatures; every member |
| L6 | proof (Killing and parallel $E$, L2) | both |
| L5/L8 line reduction (§L5) | exact, symbolic in $x,y,c,\kappa,\sigma$: coefficients on the line, $c_i=x_i(\kappa\sigma+312N\sigma^2)$, $\det M$ on the line, $F$, $\mu=-\tfrac7{18}r$, $\mu'=-\tfrac7{18}$; $\circledast_\sharp=W_A$ on $u=0$ | both signatures |
| L7 (§L6) | proof; the pairing identities on all 36 basis elements of $\mathfrak{so}\otimes\mathfrak{so}$ and the three constants; $\sum\epsilon_i\epsilon_j\mathrm{tr}R_{ij}^2=-\lvert R\rvert^2$; span ranks 21/21/22 for the generic curvature tensor; brute-force Taylor solutions (inconsistent for $\mu'=-\tfrac7{18}$, three random $R$; control consistent; order ε forces $dw=0$) | pairings: (4,0), (1,3), (3,1), (2,2); span: (4,0), (1,3), (2,2); Taylor: (1,3), (4,0), (2,2) |
| L8, L9 | proof from L1–L7 (L9 also uses Kobayashi–Nomizu I, II.10.8) | (H) |

Run: `cd computation/geometric-unity && python3 L_nonradial.py > L_nonradial.out 2>&1`. Result: 73 checks, 73 True, each asserted; 380 s on one core, a second run 380 s, with identical output apart from the timing annotations.

## L.3 Limits (layers b, d)

- **Members not covered.**
  - For $W'_A$, $W_A+W'_C\mathrm{Ad}_\Gamma$, $W'_A+W_C\mathrm{Ad}_\Gamma$ (and $W_A+W'_A$ in (14,0)) some ϖ is a nonzero rational function, so L1 applies. For $W_C\mathrm{Ad}_\Gamma$, $W'_C\mathrm{Ad}_\Gamma$ and, in (7,7), $W_A+W'_A$ and $(1,1,1,1)$, all four ϖ vanish identically, and even $dt\wedge du=0$ is open (REVIEW6-GU R8-1, `xcheck_round8.out` §X1).
  - Non-elliptic members ($\ell=0$) other than the multiples of $W_A$: for instance $W'_A$, $\mathrm{span}(W_A,W'_A)$ with $y\neq0$, $(1,0,0,1)$ in (14,0). There the curvature-free row $h$ vanishes identically (exact image condition, Cor. J2.1(b′)), so that channel has no second-order term; its first-order constraints were analysed only for the multiples of $W_A$ (L5(c)).
  - In (7,7), elliptic members with $x=\pm y$ ($N=0$): $R'\equiv0$, and on the line the Hessian coefficient vanishes.
  - $(x,y)=0$ (no curvature in either channel; outside J5's hypotheses as well).
  - (Round 9, `M-z2-remaining-members.md`.) All three classes are now treated. Up to a rotation/boost symmetry of the reduced system and a (7,7) swap with $g_\aleph\to-g_\aleph$ (M1, M2), the members with $\ell=0$ are $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, the (7,7) null members and $(x,y)=0$. There is no solution on an open set meeting $s(X)$ for any member with $(x,y)\neq0$ (M6(a), M7; with L8, L9: M10(a)), and none on any open set for $\mathrm{span}(W_A,W'_A)$ (M6(b)) and for $x=\pm y$ in (7,7) (M7); for $(x,y)=0$ only $T=0$ (M8). The members above where L1 does not apply are settled without it (M7(a), M8). Still open: $c\neq0$ ($c^2\neq\tfrac{11}6$ in (14,0)) away from $s(X)$ (M10(b)).
- **Hypotheses kept from J5.**
  - Z2: in Z1 the Euler field is not Killing, so L6 and the radial structure are unavailable. Z1 is treated in `N-z1-neighbourhood.md` (round 10): there $E$ is a gradient with $\nabla E=\tfrac12P_H$ and $\mathrm{Ric}(E,\cdot)=-E^\flat/c$; for the members of L8 there is no solution near $s(X)$ in (7,7), and in (14,0) only $\lambda=\tfrac15$ over a Ricci-flat base, flat with $T_e\equiv0$ for analytic fields (N12).
  - The real forms (14,0) and (7,7) of J1–J2; the signatures (8,6), (10,4), … are not set up.
  - $\kappa\neq0$; $T_e$ the slaved field in $\mathcal V$ with $\det M\neq0$ on $U$. Solutions crossing $\{\det M=0\}$, and $T_e$ outside $\mathcal V$ for large $|T_o|$, are not covered.
  - $T_o$ in the invariant sector.
- L2 is stated and proven for $n=4$ (traces 9 and 4; J5d L1 at $n=4$).
- L7 uses only the order-$\varepsilon^2$ jet of the fibre equation at one point of $s(X)$. L9 replaces "meets $s(X)$" by analytic continuation along one fibre.
- Reading (A), the ℷ-variation and the (P) system on $X$ of J3(a) (whose classification J.3 leaves open) are untouched. L concerns $\Upsilon=0$ on open subsets of $Y$.

## L.4 Consequences for gaps and claims

- **GU-G14.**
  - Under the hypotheses of J5 without the radial ansatz (Z2, $\kappa\neq0$, slaved $T_e\in\mathcal V$, $\det M\neq0$), $\Upsilon=0$ has no solution on any open subset of $Y$ for $\circledast_\sharp$ and $W_A$. The same holds in (14,0) for every elliptic member and in (7,7) for every member with $xz-yw\neq0$, $x\neq\pm y$, and for the multiples of $W_A$ (L8, L9).
  - So the J5/J5d no-go is not an artefact of the radial ansatz. Within this sector Z2 admits no neighbourhood solution at all. Solutions at the points of $s(X)$ (J4) remain pointwise statements with free second normal derivatives.
  - Open: Z1 beyond N12 and Q16 (limit points of Einstein pieces; open sets away from $s(X)$) and other sectors; of the remaining members only $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$ up to the symmetries M1, $c\neq0$, $c^2\neq\tfrac{11}6$ in (14,0), away from $s(X)$ (round 9, M10).
- **GU-G16.** The factor $\ell^2=(xz\pm yw)^2$ reappears in the residual $R'=-1916928\kappa^3(x^2\pm y^2)\ell^2(ty-ux)$. The cases left open are exactly those with $\ell=0$ (apart from $W_A$) or $N=0$.
- **J5.** J5(c) is recovered as the case $Z=0$ of L4 and holds without "$U$ meets $s(X)$" (L2(a)).
- **F.** Side result: Z2 is Einstein on no open subset of $Y$ ($n=4$, every $\lambda\neq\tfrac14$; L2). F6 had this on $s(X)$; Conjecture F9 holds for Z2 at $n=4$ (in the stronger form: for no base metric); for Z1 ($n=4$) it is proven near $s(X)$ over definite bases (`N-z1-neighbourhood.md`, N10); over flat bases its exception occurs exactly at $\lambda=\tfrac15$ (N2).
- **GU-054.** Status unchanged (C). No statement made before this note is contradicted. J.3 and J5d.3 listed the non-radial case as open; it is now settled in the classes above.
