# P. The last Z2 exception: M9's exceptional set is the null-curvature set, and M10(b) has no exception

Claims: GU-054 (context). Gaps: GU-G14, GU-G16. Code: `computation/geometric-unity/P_last_exception.py` → `.out` (sections P0, P1, P3, P4; exact: Python Fractions, sympy rationals, $\mathbb Q[k]$; fixed seeds; 715 s on one core, a second run 728 s, identical apart from the timing annotations). It imports `L_nonradial.py` (M9's own pipeline: `rand_curv`, `fibre_data`, `hessian_system`) and the referee scripts `xcheck_round9.py`, `xcheck_round8.py` (`relaxed_system`, `c0t`, `qt`, `curv_op`, `BR_coeff`, `pq_bases`), all unmodified. Output sections: §P0 the construction, §P1 Proposition P1, §P3 Lemma P3, §P4 Lemma P4 (Lemma P2 is checked inside §P0 and §P3). Record: `metrics/geometric-unity/P-z2-last-exception.yaml`. Builds on M (M1–M4, M6 steps 1–2, M9, M10, as corrected after REVIEW7-GU), L (L.0, L2, L6, L7's Ricci identity), J5d (L0(a)), F (F1, F3, L-F1) and REVIEW7-GU §2 ($\sum\epsilon\epsilon R_{cd}\otimes R_{cd}=-4\mathcal R^2\,\mathrm{Cas}$). Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## P.0 Setting (layer b)

**The exception.** M10(b) leaves open, under (H) of L.0, the members with $\ell=0$, $N\neq0$, $(z,w)\neq0$, up to M1 $p=(x,0,0,cx)$ with $c\neq0$ ($c^2\neq\tfrac{11}6$ in (14,0)), on open $U$ with $U\cap s(X)=\emptyset$. M3 and M6 steps 1–2 use no point of $s(X)$. They give a dense open $W\subset U$ on which $r=F(t)$ satisfies $(\star)$ with $\mu'(r)\neq0$, hence (L6) $Er=0$, $\nu=0$, and $\mu=xF'/H\neq0$. The fibres are totally geodesic, so on each $W\cap F_x$
$$\mathrm{(F)}\qquad \mathrm{Hess}^Fw=\mu(w)\,(\mathrm{Ric}_F+\mathcal A_x),\qquad w=r|_{F_x},\quad \mu(w)\neq0 .$$
M9 showed that the first integrability condition of (F) fails off the section unless $(c_V/c_H)R_g(x)$ lies in a proper algebraic set $\mathcal Z$, which it did not determine.

**Notation.**
- Fix $x$ and a $g_x$-orthonormal frame, so $g_x=\eta$, base signature (4,0), (1,3), (3,1) or (2,2). Then $F=F_x=\mathrm{Met}_{(p,q)}(\mathbb R^4)\ni y$.
- $G:=G_\lambda$; the fibre metric is $c_VG$, and the Hessian and the Ricci identity do not see $c_V$. $\Gamma(h,k)=-\tfrac12(hy^{-1}k+ky^{-1}h)$ is λ-free, and $\mathrm{Ric}_F=-G_{1/4}$ is parallel (F1).
- For $K\in\mathfrak{gl}(4)$, $\xi_K(y):=-(K^\top y+yK)$ is the vector $K\cdot y$ of L-F5. It is a Killing field, since $GL(4)$ acts by isometries.
- $\mathfrak k:=\mathfrak{so}(\eta)$, so $\xi_K(\eta)=0$ for $K\in\mathfrak k$. $K_{(ab)}$ ($a<b$) is the basis with entries $[a][b]=\eta_a$, $[b][a]=-\eta_b$; $kg$ is the trace form and $\mathrm{Cas}:=\sum(kg^{-1})^{ab}K_a\otimes K_b$.
- For $\mathfrak b\in\mathrm{Sym}^2\mathfrak{gl}(4)$ put $\mathcal A_{\mathfrak b}(U,W):=\sum\mathfrak b^{ab}G(U,\xi_{K_a})G(W,\xi_{K_b})$. By L.0 (J5d L0(a)), $\mathcal A_x=\mathcal A_{\mathfrak b_R}$ with
$$\mathfrak b_R:=c_A\sum_{i,j}\epsilon_i\epsilon_j\,R_{ij}\otimes R_{ij}\in\mathrm{Sym}^2\mathfrak k,\qquad c_A=\frac{c_V^2}{4c_H^2},\quad R_{ij}=R_g(x)(e_i,e_j).$$
- **The relaxed system** $\Sigma_y(\mathfrak b)$ at $y\in F$ has unknowns $(V,\xi)\in T_yF\times T_y^*F$ (20 in all) and the 450 scalar rows
$$G(R^F(X,Y)V,Z)-\xi(X)B(Y,Z)+\xi(Y)B(X,Z)=(\nabla_X\mathcal A_{\mathfrak b})(Y,Z)-(\nabla_Y\mathcal A_{\mathfrak b})(X,Z),\qquad B:=-G_{1/4}+\mathcal A_{\mathfrak b}.$$
  Write it as $A_y(\mathfrak b)\zeta=b_y(\mathfrak b)$. If $w$ solves (F) near $y$ with $\mu(w(y))\neq0$, then $(\nabla w/\mu,\ \mu'dw/\mu)(y)$ solves $\Sigma_y(\mathfrak b_R)$; this is the Ricci identity of L7, with $\nabla\mathrm{Ric}_F=0$.
- M9's relaxed system (295 unknowns: the 1-, 2- and 3-jets of $w$ and a covector) is equivalent to $\Sigma_y$. Its order-0 rows fix the 2-jet. Its order-1 rows can be solved for a symmetric 3-jet iff their antisymmetric part holds, and that part is $\Sigma_y$. So its ranks are those of $\Sigma_y$ plus 275 (§P0).
- **$\mathcal Z$ (now defined).** $R\in\mathcal Z$ iff $\Sigma_y(\mathfrak b_R)$ is consistent for all $y$ in some nonempty open subset of $F$. $\mathcal Z_{\mathfrak k}$ is the corresponding set of $\mathfrak b\in\mathrm{Sym}^2\mathfrak k$.

## P.1 Results (layer c)

**Proposition P1 (the A-term sees only $\mathcal R^2$).** *Proven* (checks §P1).

(a) $\mathcal A_x$ depends on $R_g(x)$ only through $\mathfrak b_R$, and linearly. Moreover $\mathfrak b_R=-4c_A\,\mathcal R^2\cdot\mathrm{Cas}$, with $\mathcal R$ the curvature operator on $\Lambda^2\cong\mathfrak k$ (REVIEW7-GU §2). So $\mathfrak b_R=0$ iff $\mathcal R^2=0$.

(b) In base (4,0), $\mathfrak b_R=0$ iff $R=0$. In bases (1,3), (3,1), (2,2) there are $R\neq0$ with $\mathfrak b_R=0$, for example
$$R=\sum H_{ij}\,(\ell\wedge e_i)\otimes(\ell\wedge e_j),$$
with $\ell$ null, $e_i\perp\ell$ and $H$ symmetric (pp-wave type). Their Ricci tensor is a multiple of $\ell\otimes\ell$, nonzero in general. For them $\mathcal A_x\equiv0$ on the whole fibre, although $A\neq0$ off $s(X)$.

(c) For $P\in GL(4)$, the isometry $y\mapsto PyP^\top$ identifies $\Sigma_{P\eta P^\top}(\mathfrak b)$ with $\Sigma_\eta(\mathfrak b^P)$. Here $\mathfrak b^P$ is $\mathfrak b$ with every $K$ replaced by $P^\top KP^{-\top}$ (M9's pull-back).

*Proof.*
- (a) $\mathcal A_x=c_A\sum\epsilon_i\epsilon_j\,\xi_{R_{ij}}^\flat\otimes\xi_{R_{ij}}^\flat$ is the image of $\mathfrak b_R$ under the linear map $K\otimes K'\mapsto\xi_K^\flat\otimes\xi_{K'}^\flat$. The identity is REVIEW7-GU's, re-checked in four signatures.
- (b) $\mathcal R$ is self-adjoint for the trace form, which is definite on $\mathfrak{so}(4)$. So $\mathcal R$ is diagonalisable over $\mathbb R$, and $\mathcal R^2=0$ forces $\mathcal R=0$.
- (b), examples. The 2-forms $\alpha_i=\ell\wedge e_i$ satisfy $\langle\alpha_i,\alpha_j\rangle=0$, so $\mathcal R^2=0$; first Bianchi holds because $\alpha_i\wedge\alpha_j=0$.
- (c) $\phi_*\xi_K=\xi_{P^{-\top}KP^\top}$ for $\phi(y)=PyP^\top$. ∎

**Lemma P2 (algebraic structure of $\Sigma$).** *Proven* (rank checks §P0, §P3).

(a) The entries of $A_y(\mathfrak b)$ and $b_y(\mathfrak b)$ are polynomials in $(y,(\det y)^{-1})$. $A$ is affine in $\mathfrak b$, and $b$ is linear in $\mathfrak b$.

(b) The dilation field $y$ (a constant multiple of $\nabla s$), with $\xi=0$, lies in $\ker A_y(\mathfrak b)$ for every $(y,\mathfrak b)$. Hence $\mathrm{rank}\,A\le19$.

(c) At the section $y=\eta$: $\mathcal A_{\mathfrak b}(\eta)=0$ and $\nabla\mathcal A_{\mathfrak b}(\eta)=0$ for $\mathfrak b\in\mathrm{Sym}^2\mathfrak k$. So $A_\eta=:A_0$ does not depend on $\mathfrak b$, it has rank 19, and $b_\eta=0$.

(d) *Dichotomy.* For $\mathfrak b\in\mathrm{Sym}^2\mathfrak k$ exactly one of the following holds.
- (i) Some 20-minor of $[A|b]$ is not identically zero on $F$. Then $\Sigma_y(\mathfrak b)$ is inconsistent on a dense open subset of $F$.
- (ii) All 20-minors vanish identically. Then $\Sigma_y(\mathfrak b)$ is consistent wherever $\mathrm{rank}\,A_y=19$, in particular on a neighbourhood of $\eta$.

  $\mathfrak b\in\mathcal Z_{\mathfrak k}$ iff (ii) holds.

*Proof.*
- (a) $G$, Γ, the Killing fields and their covariant derivatives are polynomial in $y$ and $y^{-1}=\mathrm{adj}(y)/\det y$.
- (b) $\nabla s$ is parallel, so $R^F(\cdot,\cdot)\nabla s=0$.
- (c) $\xi_K(\eta)=0$ for $K\in\mathfrak k$. The rank is computed.
- (d) A 20-minor is a polynomial times a power of $(\det y)^{-1}$. Its numerator is a polynomial on $\mathrm{Sym}^2$ and $F$ is open in $\mathrm{Sym}^2$, so a minor that is not identically zero on $F$ is nonzero on a dense open subset of $F$, and one that vanishes on a nonempty open subset vanishes identically (connectedness is not needed; R12-5(a)). There $\mathrm{rank}[A|b]=20>\mathrm{rank}\,A$.
- In case (ii), $\mathrm{rank}[A|b]\le19$ everywhere, so the system is consistent wherever $\mathrm{rank}\,A=19$. Near η, $\mathrm{rank}\,A\ge19$ by semicontinuity.
- Consistency on an open set gives $\mathrm{rank}[A|b]\le19$ there, so every 20-minor vanishes there, hence identically. ∎

**Lemma P3 (first order at the section).** *Computed* (exact, two constructions, §P3); the consequence is proven. If $\mathfrak b\in\mathcal Z_{\mathfrak k}$, then $\mathfrak b\in\mathbb R\,\mathrm{Cas}$.

*Proof.*
1. Let $X\in\mathfrak{gl}(4)$ and $P_\varepsilon=1+\varepsilon X$. By P1(c), $\Sigma_{P_\varepsilon\eta P_\varepsilon^\top}(\mathfrak b)$ is $\Sigma_\eta(\mathfrak b^{P_\varepsilon})$; write it $A(\varepsilon)\zeta=b(\varepsilon)$.
2. $K'_a(\varepsilon)=K_a+\varepsilon[X^\top,K_a]+O(\varepsilon^2)$, so $\xi_{K'_a}(\eta)=O(\varepsilon)$ and $\mathcal A=O(\varepsilon^2)$. Hence $A(\varepsilon)=A_0+O(\varepsilon^2)$ and $b(\varepsilon)=\varepsilon\,b_1(\mathfrak b,X)+O(\varepsilon^2)$, with $b_1$ linear in $\mathfrak b$.
3. By P2(c), (d)(ii) the system is consistent and $\mathrm{rank}\,A(\varepsilon)=19$ for small ε. Fix a nonsingular $19\times19$ minor $A_{0,IJ}$ of $A_0$. For small ε, $A_{IJ}(\varepsilon)$ is invertible and the rows $I$ span the row space of $[A(\varepsilon)|b(\varepsilon)]$ (both have rank 19). So $\zeta_J=A_{IJ}(\varepsilon)^{-1}b_I(\varepsilon)$, $\zeta_{J^c}=0$ is an analytic solution, with $\zeta(0)=0$ since $b(0)=0$ (R12-5(b)). At order ε: $A_0\zeta_1(X)=b_1(\mathfrak b,X)$ for every $X$.
4. Computed, in each of the four base signatures: the linear system "$b_1(\mathfrak b,X)=A_0\zeta_1(X)$ for the 16 matrix units $X$" has 7200 rows and 341 unknowns ($\mathfrak b\in\mathrm{Sym}^2\mathfrak k$ and the $\zeta_1(X)$). Its solution space has dimension $17=16+1$, and its $\mathfrak b$-projection is $\mathbb R\,\mathrm{Cas}$. ∎

*Second construction.* At first order $b_1(\mathfrak b,X)$ is L7's $q(\mathfrak b)(X,\cdot,\cdot)$. The image of $A_0$ on $T S$ is spanned by $c_0(V\,\cdot,\cdot,\cdot)$ and $\Xi(Y)Z-\Xi(Z)Y$. So P3 says that the **relaxed (C1)**
$$q(\mathfrak b)(X,Y,Z)=c_0(VX,Y,Z)+\Xi(X)(Y)\,Z-\Xi(X)(Z)\,Y,\qquad V\in\mathrm{End}\,\mathfrak p,\ \Xi\in\mathrm{Hom}(\mathfrak p,\mathfrak p^*)\text{ arbitrary},$$
holds only for $\mathfrak b\in\mathbb R\,\mathrm{Cas}$. Recomputed with the referee's (C1) tensors (`xcheck_round9` `c0t`, `qt`; 9044 equations, 183 unknowns): nullity 1, $\mathfrak b\in\mathbb R\,\mathrm{Cas}$, four signatures. L7's (C1) is the special case $V\propto\mathrm{Id}$, $\Xi\propto\mu'\,(\cdot)^\flat$, and $q(\mathrm{Cas})=-3c_0$ (REVIEW7-GU).

**Lemma P4 (the Casimir line).** *Computed* (exact, §P4); the consequence is proven. If $k\,\mathrm{Cas}\in\mathcal Z_{\mathfrak k}$, then $k=0$.

*Proof.*
1. Continue the expansion of P3 for $\mathfrak b=k\,\mathrm{Cas}$. Then $b_m=k\,b_m^{\rm Cas}$, and $A(\varepsilon)=A_0+\varepsilon^2kA_2^{\rm Cas}+O(\varepsilon^3)$, where $A_2$ (the A-term at order $\varepsilon^2$) acts only on the ξ-columns. Also $\zeta_1=k\,\zeta_1^{\rm Cas}+(\ker A_0)$, and $A_2$ annihilates $\ker A_0$.
2. At order $\varepsilon^3$: $A_0\zeta_3+A_2\zeta_1=b_3$, i.e. $A_0\zeta_3=k\,b_3^{\rm Cas}-k^2A_2^{\rm Cas}\zeta_1^{\rm Cas}$.
3. Computed at a seeded rational $X$ in each signature: $A_2^{\rm Cas}\zeta_1^{\rm Cas}\in\mathrm{Im}\,A_0$, but $b_3^{\rm Cas}\notin\mathrm{Im}\,A_0$ (rank 20).
4. So $k\,b_3^{\rm Cas}\in\mathrm{Im}\,A_0$, and $k=0$. ∎

Independent confirmations (§P4):
- *Order 2.* $b_2^{\rm Cas}\in\mathrm{Im}\,A_0$: Cas passes order 2, so order 3 is needed.
- *Global.* At two seeded off-section points, a 20-minor of $[A(k)|b(k)]$ (with the right-hand column) is a polynomial in $k$. The two polynomials have gcd $k^2$ in each signature, so by P2(d)(i) $k\,\mathrm{Cas}\notin\mathcal Z_{\mathfrak k}$ for $k\neq0$. This uses no expansion.
- *At the other roots.* At the nonzero roots of the first minor (rational, two per signature), three constructions are inconsistent at the first point: the own system, the referee's direct Taylor system (`relaxed_system`, $R=R_1$, $c_A=-k/4$) and M9's pull-back Taylor system (`L.fibre_data`).
- *Example.* $k=\tfrac1{12}$ is M5(d)'s value: the referee's $R_L$ at $(c_V/c_H)^2=\tfrac1{12}$ has $\mathfrak b_{R_L}=\mathrm{Cas}/12$. It passes (C1)$_0$ at the section point, but fails on every open subset of the fibre.

**Theorem P5 ($\mathcal Z$ determined).** *Proven* from P1–P4.
- $\mathcal Z_{\mathfrak k}=\{0\}$. Hence
$$\mathcal Z=\{R:\ \mathfrak b_R=0\}=\{R:\ \mathcal R^2=0\}=\{R:\ \mathcal A_x\equiv0\text{ on }F_x\}.$$
- In base (4,0), $\mathcal Z=\{0\}$. In bases (1,3), (3,1), (2,2), $\mathcal Z\supsetneq\{0\}$: it contains the null curvatures of P1(b).
- For $R\notin\mathcal Z$, $\Sigma_y(\mathfrak b_R)$ is inconsistent on a dense open subset of $F_x$.

*Proof.*
- $\mathfrak b\in\mathcal Z_{\mathfrak k}$ gives P2(d)(ii), then $\mathfrak b\in\mathbb R\,\mathrm{Cas}$ by P3 and $\mathfrak b=0$ by P4. Conversely, $\mathfrak b=0$ gives $b\equiv0$.
- If $\mathcal A_x\equiv0$ on $F_x$, then $b\equiv0$, so $\mathfrak b_R\in\mathcal Z_{\mathfrak k}=\{0\}$. Hence $\mathcal A_x\equiv0\iff\mathfrak b_R=0\iff\mathcal R^2=0$ (P1(a)). ∎

**Theorem P6 (the fibre equation has no local solution with $\mu\neq0$).** *Proven.* Data:
- $x\in X$, with base signature (4,0), (1,3), (3,1) or (2,2); any $\lambda\neq\tfrac14$, any $c_Hc_V\neq0$ and any curvature $R_g(x)$;
- $\Omega\subset F_x$ nonempty and open;
- $\mu$ a $C^1$ function on an interval containing $w(\Omega)$, and $w\in C^3(\Omega)$ with $\mathrm{Hess}^Fw=\mu(w)(\mathrm{Ric}_F+\mathcal A_x)$ on Ω.

Then $\mu(w)\equiv0$ on Ω.

*Proof.* Suppose $\mu(w(q_0))\neq0$, and shrink Ω so that $\mu(w)\neq0$ on Ω.
- If $R_g(x)\notin\mathcal Z$, P5 gives a $q\in\Omega$ where $\Sigma_q$ is inconsistent. But $(\nabla w/\mu,\mu'dw/\mu)(q)$ solves it.
- If $R_g(x)\in\mathcal Z$, then $\mathcal A_x\equiv0$ (P5). M9's flat-case argument (as corrected after R9-2) uses only $\mathcal A_x\equiv0$, and it gives $\mu(w)=0$. ∎

P6 needs no point of $s(X)$. It implies L7 (the case $\Omega\ni g(x)$). It also strengthens M5, whose order-$\varepsilon^2$ test at the section admits Lorentzian $R$ with $\mathcal R^2=-\tfrac1{12}(c_H/c_V)^2\mathrm{Id}$ (M5(d)): no Ω, with or without $g(x)$, carries a solution with $\mu(w)\neq0$.

**Corollary P7 (M10(b) has no exception).** *Proven.* Assume (H). For every member with $(x,y)\neq0$, in (14,0) and in (7,7), $\Upsilon=0$ has **no solution on any nonempty open $U\subset Y$**. With M8: under (H) without $(x,y)\neq0$, $\Upsilon=0$ holds on a nonempty open set iff $(x,y)=0$ and $T=0$.

*Proof.*
1. M10's proof covers every member except the exception of M10(b): L8(A), L8(B) with L9, M6(b), M6(c), M7.
2. Exception. By M1 and M2(a) it suffices to take $p=(x,0,0,cx)$ on $(U,\pm g_\aleph)$. The constants $(\pm c_H,\pm c_V)$ belong to the setting, and P6 holds for all of them.
3. M3 ($u\equiv0$, $\{dt\neq0\}$ dense) and M6 step 2 use no point of $s(X)$. M6 step 2 gives the dense open $W=\{dt\neq0,\ H(t)\neq0,\ \mu'(t)\neq0\}$. On $W$, $(\star)$ holds with $\mu'(r)\neq0$, and L6 gives $Er=0$, $\nu=0$.
4. On $W\cap F_x$ the function $w=r|_{F_x}$ solves (F) (O'Neill, $T=0$), with $\mu(w)=xF'/H\neq0$: $x\neq0$, $F'\neq0$, and $H$ is finite and nonzero on $W$.
5. $W\neq\emptyset$. For $x_0\in\pi(W)$, $W\cap F_{x_0}$ is a nonempty open subset of $F_{x_0}$, and P6 gives a contradiction. ∎

**Remark P8 (one argument for every line equation).** *Proven.* Wherever L and M reduce $\Upsilon=0$ on an open set to $(\star)$ with $\mu'\neq0$ and $\mu(r)\neq0$ on a nonempty open subset, L6 and P6 give the contradiction directly. The cases are:
- L8(B) step 5: $\mu=-\tfrac7{18}r$;
- M6(b), (c): $\mu=-\tfrac7{18}r$ and $-\tfrac3{289}r$;
- M7(b): $\mu=-\tfrac5{12}r$;
- M6(a), on $W$.

So L7, M5, L9's analytic continuation and M6 steps 4–6 are not needed for L9, M6, M7(b) and M10. They remain correct; they are the order-$\varepsilon^2$ statements at the section that P6 strengthens. For the members with linear μ the route through L7 and L9 needs neither P3 nor P4 and remains their computation-light proof (R12-5(f)). The hypothesis "$U$ meets $s(X)$" appears nowhere.

**Answers.**
1. *$\mathcal Z$.* It is $\{\mathcal R^2=0\}$: $\{0\}$ in base (4,0), and in indefinite bases also the null curvatures ($R\neq0$, $\mathcal R^2=0$). They include the pp-wave type tensors of P1(b) and, in Lorentzian bases, tensors of other type, e.g. $F_1\otimes F_1-F_2\otimes F_2$ with $F_1=e^{01}+e^{23}$, $F_2=e^{02}+e^{31}$ (REVIEW11-GU R12-3). On $\mathcal Z$ the A-term vanishes identically on the fibre, and the flat-case argument applies verbatim (P5, P6). In M5(d)'s terms: of the conditions $\mathcal R^2=k\,\mathrm{Id}$ (the Casimir line, which contains M5(d)'s value), only $k=0$ survives. The line passes orders 1 and 2 at the section but fails at order 3 (P3, P4).
2. *Horizontal equations, nonlinearity of μ.* Neither is needed. P6 holds for every μ, because the relaxation already discards $\mu'$.
3. *Flat base.* M9's flat case gives $\mu(w)=0$, against $\mu=xF'/H\neq0$. So flat bases are excluded on $U$, and by the same argument every $R\in\mathcal Z$.
4. *No solution exists.* The exception of M10(b) is empty, and under (H) the Z2 family has no neighbourhood solution for $(x,y)\neq0$ (P7).

## P.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| construction (§P0) | at a random fibre point: Γ metric and $G_{1/4}$ parallel (exact derivatives by dual numbers), $\mathrm{Ric}(R^F)=-G_{1/4}$ (F1), $\xi_K$ Killing ($K\in\mathfrak{gl}(4)$ random), the Killing-field formula for $\nabla\mathcal A$ against the exact derivative; a manufactured solution ($w$ cubic in $y$, $B:=\mathrm{Hess}\,w$) satisfies all 450 rows; the rows are identical for λ = 0 or ½, −1, 2 | (4,0), (1,3), (3,1), (2,2) |
| equivalence with M9 (§P0) | M9's own samples (same seeds, same pull-back): $(19,20)$ three times and $(19,19)$ at the section, against M9's $(294,295)$, $(294,294)$; dilation kernel; the referee's direct construction at a random $(R,y_0)$: $(294,295)$ against own $(19,20)$ | M9: four signatures; referee: (4,0), (3,1) (xcheck_round9 §Z8 has (1,3), (2,2)) |
| P1 (§P1) | $\sum\epsilon\epsilon R_{cd}\otimes R_{cd}=-4\mathcal R^2kg^{-1}$ (referee's `curv_op`, `BR_coeff`) and $\mathcal A,\nabla\mathcal A$ from $R_{cd}$ = from $\mathfrak b_R$, random $R$; pp-wave $R$: Bianchi, $R\neq0$, $\mathrm{Ric}\neq0$, $\mathcal R^2=0$, $\mathfrak b_R=0$, $\mathcal A\equiv0$ and $\Sigma$ consistent at two fibre points; trace form definite in (4,0) | four signatures; pp-wave (1,3), (3,1), (2,2) |
| P2 | proof; rank 19 at the section, dilation kernel (§P0, §P3) | four signatures |
| P3 (§P3) | first-order system, 7200 × 341, nullity 17, $\mathfrak b$-projection $\mathbb R\,\mathrm{Cas}$; relaxed (C1) with the referee's tensors, 9044 × 183, nullity 1, $\mathbb R\,\mathrm{Cas}$ | four signatures each |
| P4 (§P4) | order 1, 2 pass; $A_2\zeta_1\in\mathrm{Im}A_0$, $b_3\notin\mathrm{Im}A_0$; two 20-minors in $\mathbb Q[k]$ with gcd $k^2$; at the nonzero rational roots own, referee's and M9's constructions inconsistent; $\mathfrak b_{R_L}=\mathrm{Cas}/12$ at $(c_V/c_H)^2=\tfrac1{12}$ | four signatures; $R_L$ (1,3) |
| P5–P8 | proof from P1–P4, M3, M6 steps 1–2, L6, M9's flat case | (H) |

Run: `cd computation/geometric-unity && python3 P_last_exception.py > P_last_exception.out 2>&1`. Result: 39 checks, 39 True, each asserted; 715 s on one core, a second run 728 s, with identical output apart from the timing annotations.

## P.3 Limits (layers b, d)

- **P3, P4 are computations.** They are finite exact linear algebra over $\mathbb Q$ and $\mathbb Q[k]$, done in each base signature. There is no hand proof.
  - Representation theory gives a frame: over $\mathbb C$, $\mathrm{Sym}^2\mathfrak k=(5,1)\oplus(1,5)\oplus(3,3)\oplus2(1,1)$, and the solution space of P3 is a $K$-submodule.
  - The computation shows that only one trivial summand survives, and P4 removes it. Neither step is derived by hand.
- **P6 uses only the first integrability condition** (plus the flat case). It says nothing about solutions with $\mu\equiv0$ on the range of $w$, e.g. $w=as+b$.
- **Hypotheses kept from (H).** Z2; $n=4$; the real forms (14,0), (7,7) of J1–J2 (bases (4,0); (1,3), (3,1), (2,2)); $\kappa\neq0$; $T_e$ the slaved field in $\mathcal V$ with $\det M\neq0$ on $U$; $T_o$ in the invariant sector. Not covered:
  - solutions crossing $\det M=0$;
  - $T_e$ outside $\mathcal V$ (round 14, `R-even-operator.md` R5, not yet internally reviewed: excluded under (R), a condition on $(t,u)$; still open on open sets where $(t,u)\in Z_\perp$; in (7,7) excluded for $(x,y)=0$ and, up to M1, $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$ with $c^2>\tfrac{25}{144}$, R7);
  - other sectors;
  - the Z2 signatures (8,6), (10,4), ….
- **M1** is a symmetry of the reduced system only (as in M.3).
- **Z1 is not covered.** P6 is a statement about Z2's fibre and its O'Neill block $\mathrm{Ric}_F+\mathcal A$ (F1, F3). Whether the same pointwise method closes the Z1 cases of N.3 was not examined.
- **Untouched:** reading (A); the (P) system of J3(a); J3(b), (c) and J4(d′) for these members.

## P.4 Consequences for gaps and claims

- **GU-G14.**
  - Under (H), the non-radial Z2 neighbourhood question is settled for the whole family. For $(x,y)\neq0$, $\Upsilon=0$ has no solution on any nonempty open subset of $Y$ (P7). For $(x,y)=0$ the only solution is $T=0$, for every metric (M8).
  - The exception of M10(b) is empty, so the slaved vector sector of Z2 carries no neighbourhood solution except $T=0$ for the members without curvature.
  - The fibre equation $\mathrm{Hess}^Fw=\mu(w)(\mathrm{Ric}_F+\mathcal A)$ has no local solution with $\mu\neq0$, for any base curvature and any μ (P6). This replaces the section-based arguments of L7, L9, M5 and M6 steps 4–6 (P8).
  - Open: Z1 beyond N; other sectors; the rule for $T$ off $s(X)$.
- **GU-G16.** $\ell=xz\pm yw$ decides which channel does the work: the curvature-free channel for $\ell\neq0$ (L5, L8(A)), and for $\ell=0$ the first-order channel and then the fibre theorem (M3, P6). It no longer separates solvable from unsolved cases: in Z2 under (H) no member with $(x,y)\neq0$ has a neighbourhood solution.
- **GU-054.** Status unchanged (C).
  - No committed statement is contradicted. M9's $\mathcal Z$ is now determined. It is not $\{0\}$ in indefinite bases, so the "$\mathcal Z\setminus\{0\}$" of M10(b) is nonempty in (7,7). But those curvatures have $\mathcal A\equiv0$ and fail by M9's own flat-case argument.
  - M.3's open item and L.3/L.4's "still open" are answered.
