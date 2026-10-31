# J5d. The leaf condition of Theorem J5 in indefinite signature: no radial Z2 solution in (7,7)

Claims: GU-054 (context). Gaps: GU-G14. Code: `computation/geometric-unity/J5d_leaf.py` → `.out` (exact, sympy rationals; fixed seed; 181 s on one core, a second run 170 s, identical apart from the timing annotations). The output sections carry the labels of the results: §L0 Lemma L0, §L1 Theorem L1, §L2 Corollary L2, §L3 Proposition L3, §L4 Proposition L4. It imports no repository code. Record: `metrics/geometric-unity/J5d-indefinite-leaf.yaml`. Builds on J (Theorem J5 and its proof), F (Theorem F3, L-F1, L-F5) and A (A4, A5). Refereed in `reports/reviews/REVIEW5-GU-2026-10-30.md` (R7-1, R7-2, R7-5 applied here). Layers: (b) reconstruction, (c) proven.

## J5d.0 Question and setting (layer b)

Theorem J5 (`J-slaved-sector.md`) concerns the radial ansatz $t=f(s)$, $u=h(s)$ in Z2.
- Hypotheses: $\lambda\neq\tfrac14$, $(x,y)\neq0$, $\kappa\neq0$, $T_e$ the slaved field of J1 with $\det M\neq0$, and $\Upsilon=0$ on a connected open $U$ meeting $s(X)$.
- The setting is that of I.0/J.0: J1–J2, and hence J5(b), are set up and computed for the two real forms (14,0) and (7,7).
- J5(b): the leaves $\{s={\rm const}\}$ are Einstein, and at $x_0=\pi(p)$, $p\in U\cap s(X)$, $\mathrm{Ric}(g)=\Lambda_g\,g$ with $\Lambda_g=-c_H/c_V\neq0$.
- J5(c) proves "no solution" when $g_\aleph$ is definite, i.e. in (14,0).
- J5(d) leaves the indefinite case open. In the setting this is (7,7), which comes from a base (1,3) or (3,1) with $\lambda>\tfrac14$ (A4), or from a base (2,2) with $\lambda>\tfrac14$ (A5).
- In the indefinite case a solution needs, besides J5(b), the order-$\varepsilon^2$ leaf condition $F\equiv0$ below. That condition fails for every base of nonzero constant curvature.

This note decides the open case at the level of the algebraic curvature tensor $R=R_g(x_0)$.

**Notation.**
- $(V,g)=(T_{x_0}X,g(x_0))$, $n=\dim V$, any signature. $e_a$ is an orthonormal frame with $g(e_a,e_b)=\eta_a\delta_{ab}$, and $\theta^a$ its dual coframe.
- Conventions of F: $R_{abcd}=g(R(e_c,e_d)e_b,e_a)$, $\mathrm{Ric}_{bd}=g^{ac}R_{abcd}$, and $R(X,Y,Z,W):=g(R(X,Y)Z,W)$. So $R(X,Y)\in\mathfrak{so}(g)$.
- An endomorphism $B$ acts on symmetric bilinear forms by $(B\cdot y)(Z,W)=-y(BZ,W)-y(Z,BW)$, i.e. $B\cdot y=-(B^{\mathsf T}y+yB)$ in matrices. This is the action in L-F5's formula $A_XY=-\tfrac12(R_g(X,Y)\cdot y)^{\rm vert}$ (F3).
- For a nondegenerate symmetric $y$, $X,Y\in V$ and $h,h'\in\mathrm{Sym}^2V^*$:
$$\Phi_y(X,Y)=\sum_{ij}g^{ij}\,\mathrm{tr}\big(y^{-1}(R(X,e_i)\cdot y)\,y^{-1}(R(Y,e_j)\cdot y)\big),$$
$$F(X,Y;h,h')=\sum_{ij}g^{ij}\,G_0\big(R(X,e_i)\cdot h,\,R(Y,e_j)\cdot h'\big),\qquad G_0(h,k)=\mathrm{tr}(g^{-1}hg^{-1}k),$$
  and $F(X,Y;h):=F(X,Y;h,h)$. This is J5(d)'s form: $\sum_i\epsilon_i(\cdot)(E_i,E_i)$ over a $g_\aleph$-orthonormal horizontal frame equals $c_H^{-1}\sum_{ij}g^{ij}(\cdot)(e_i,e_j)$, and the factor does not matter.
- **"$F\equiv0$" means $F(X,Y;h,h)=0$ for all $X,Y,h$.** This is what the leaf condition gives (L0(c)).
  - Polarising in $h$ gives only $F(X,Y;h,h')+F(X,Y;h',h)=0$.
  - The bilinear form itself is **not** symmetric under $h\leftrightarrow h'$. Since $F(X,Y;h',h)=F(Y,X;h,h')$, this is the same as asymmetry under $X\leftrightarrow Y$ (Theorem L1(c)).
  - So only contractions built from the quadratic values $F(X,Y;h,h)$ are forced to vanish.
- Two such contractions, with $\theta^p\odot\theta^q=\tfrac12(\theta^p\otimes\theta^q+\theta^q\otimes\theta^p)$, $h_{aq}:=\theta^a\odot\theta^q$ and $Q_{ab}(h):=F(e_a,e_b;h,h)$:
$$C_1=\sum_{a,p,q}\eta_a\eta_p\eta_q\,Q_{aa}(\theta^p\odot\theta^q),\qquad C_2^{\rm sym}=\sum_{a,b,q}\eta_q\,\tfrac12\big[Q_{ab}(h_{aq}+h_{bq})-Q_{ab}(h_{aq})-Q_{ab}(h_{bq})\big].$$
  - $C_1$ is the trace over $X$ together with the $G_0$-trace over $h$.
  - $C_2^{\rm sym}$ contracts the $X$-slot with an index of $h$.
  - Both are linear combinations of values of $Q$, so $F\equiv0$ gives $C_1=C_2^{\rm sym}=0$.
  - In terms of the bilinear form, $C_2^{\rm sym}=\tfrac12(C_2+C_2')$ with $C_2=\sum_{a,b,q}\eta_qF(e_a,e_b;h_{aq},h_{bq})$ and $C_2'$ the same sum with $h_{aq},h_{bq}$ swapped.
  - $C_2$ and $C_2'$ separately are **not** forced to vanish by the leaf condition.
  - (Corrected after review 2026-10-30, R7-1; see `reports/CORRECTIONS.md`. The first version used $C_2$ and stated "Both are linear in $F$, so $F\equiv0$ gives $C_1=C_2=0$".)
- $|R|^2=R_{abcd}R^{abcd}$, $|\mathrm{Ric}|^2=\mathrm{Ric}_{ab}\mathrm{Ric}^{ab}$, and $R_1$ is the unit constant-curvature tensor, $(R_1)_{abcd}=g_{ac}g_{bd}-g_{ad}g_{bc}$. An Einstein tensor with $\mathrm{Ric}=\Lambda g$ in $n=4$ is $R=W+\tfrac\Lambda3R_1$, with $W$ an arbitrary Weyl tensor (10 parameters in every signature), and $|R_1|^2=24$.

## J5d.1 Results (layer c)

**Lemma L0 (J5(d)'s reduction, re-derived).** Assume the hypotheses of J5 and let $p\in U\cap s(X)$, $x_0=\pi(p)$. Let $\mathcal W$ be a connected open neighbourhood of $g(x_0)$ in the unimodular slice $\{\det(g(x_0)^{-1}y)=1\}$ such that $\{x_0\}\times\mathcal W$ lies in the connected component of $(\text{leaf through }p)\cap U$ that contains $p$. On that component the Einstein constant of J5(b) is constant (Schur). (Corrected after review 2026-10-30, R7-5; the text had "the open neighbourhood".) Put $R=R_g(x_0)$.

(a) For horizontal $X,Y$ at $(x_0,y)$, $(A_X,A_Y)=\dfrac{c_V}{4c_H}\,\Phi_y(X,Y)$.
- The λ-term of $G_\lambda$ drops at every $y$, not only at $y=g$, because $\mathrm{tr}(y^{-1}(B\cdot y))=-2\,\mathrm{tr}\,B=0$ for $B\in\mathfrak{so}(g)$. For the same reason $A_XY$ is tangent to the unimodular slice.
- With $M=g^{-1}y$, $B_i=R(X,e_i)$ and $C_j=R(Y,e_j)$:
$$\Phi_y(X,Y)=2\sum_{ij}g^{ij}\big[\mathrm{tr}(B_iC_j)-\mathrm{tr}(M^{-1}B_iMC_j)\big],\qquad \Phi_{cy}=\Phi_y\ (c\neq0).$$

(b) The leaf condition is $\Phi_y(X,Y)=0$ for all $X,Y\in V$ and all $y\in\mathcal W$. $\det(y)^2\Phi_y$ is a polynomial in $y$ and $\Phi_y$ is homogeneous of degree 0. So the condition is equivalent to $\Phi_y\equiv0$ for every nondegenerate symmetric $y$, a finite system of polynomial identities in $R$ (the "full condition").

(c) At $y=g+\varepsilon h$ with $\mathrm{tr}(g^{-1}h)=0$:
- $\Phi_{g+\varepsilon h}(X,Y)=\varepsilon^2F(X,Y;h,h)+O(\varepsilon^3)$; the $\varepsilon^0$ and $\varepsilon^1$ terms vanish identically.
- Hence $F(X,Y;h,h)=0$ for all $X,Y$ and all $g$-traceless $h$.
- $R\cdot g=0$ gives $F(X,Y;h+cg,h+cg)=F(X,Y;h,h)$, so $F\equiv0$ for all $X,Y$ and all $h\in\mathrm{Sym}^2V^*$.
- $F(X,Y;h,h)$ is symmetric in $X,Y$ and homogeneous of degree 2 in $R$.

So the reduction stated in J5(d) is correct: the action $R\cdot h$, the signs $\epsilon_i$, the restriction to traceless $h$ (immaterial) and the second-order reduction all check out.

*Proof.*
- (a) $A_XE_i=-\tfrac12(R(X,E_i)\cdot y)^{\rm vert}$ (L-F5, F3), and $g_\aleph|_V=c_VG_\lambda|_y$. Also $y^{-1}(B\cdot y)=-(B^*+B)$ with $B^*:=y^{-1}B^{\mathsf T}y$, whose trace is $-2\,\mathrm{tr}B$. Expanding $\mathrm{tr}((B^*+B)(C^*+C))$ with $\mathrm{tr}(B^*C^*)=\mathrm{tr}(BC)$ and $\mathrm{tr}(BC^*)=\mathrm{tr}(B^*C)$ gives $2\mathrm{tr}(BC)+2\mathrm{tr}(B^*C)$. For $B\in\mathfrak{so}(g)$, $B^{\mathsf T}=-gBg^{-1}$, so $B^*=-M^{-1}BM$.
- (b) As in the proof of J5(c): $Y'\to X$ is a pseudo-Riemannian submersion with totally geodesic, nondegenerate fibres ($c_VG_0$ on the traceless part) and the A-tensor of $Y$. O'Neill (L-F1, $T=0$) gives $\mathrm{Ric}^{Y'}(X,Y)=\mathrm{Ric}_g(X,Y)-2(A_X,A_Y)$. At the points $(x_0,y)$, $y\in\mathcal W$, both Ricci terms equal $-\tfrac{c_H}{c_V}g(X,Y)$ (J5(b), on the component containing $p$).
- (c) $y^{-1}=g^{-1}+O(\varepsilon)$, $R\cdot y=\varepsilon R\cdot h$, and $T_g\mathcal W=\{h:\mathrm{tr}(g^{-1}h)=0\}$. Symmetry of $F(X,Y;h,h)$ in $X,Y$ follows from the cyclicity of the trace and $g^{ij}=g^{ji}$. ∎

Output §L0 checks (a)–(c) exactly at random rational data in bases (4,0), (1,3), (3,1), (2,2). It also reproduces the three constant-curvature polynomials $F(e_0,e_0;h,h)$ printed in `J_slaved_sector.out` §J5 in the same parametrisation, so the conventions agree with J5(d)'s computation.

**Theorem L1 (two traces of the quadratic form $F$).** Let $R$ be any algebraic curvature tensor on $(V,g)$, any $n\ge2$, any signature.

(a) The traces satisfy
$$C_1=(n+2)\,|R|^2,\qquad C_2^{\rm sym}=\tfrac{n+8}8\,|R|^2+\tfrac{n+4}4\,|\mathrm{Ric}|^2,$$
and hence
$$8(n+2)\,C_2^{\rm sym}-(n+8)\,C_1=2(n+2)(n+4)\,|\mathrm{Ric}|^2 .$$
For $n=4$: $4C_2^{\rm sym}-C_1=8\,|\mathrm{Ric}|^2$.

(b) If $\mathrm{Ric}=\Lambda g$, then $|\mathrm{Ric}|^2=n\Lambda^2$. So for $n=4$ and $R=W+\tfrac\Lambda3R_1$:
$$4\,C_2^{\rm sym}-C_1=32\,\Lambda^2,\qquad C_1=6|W|^2+16\Lambda^2,\qquad C_2^{\rm sym}=\tfrac32|W|^2+12\Lambda^2 .$$
In particular, **$F\equiv0$ forces $\Lambda=0$, for every Einstein tensor in every signature.**

(c) *Side identities (not consequences of the leaf condition).* The bilinear traces are
$$C_2=\tfrac34|R|^2+\tfrac{n+3}2|\mathrm{Ric}|^2,\qquad C_2'=\tfrac{n+5}4|R|^2+\tfrac12|\mathrm{Ric}|^2 .$$
- $C_2-C_2'=\tfrac{n+2}2|\mathrm{Ric}|^2-\tfrac{n+2}4|R|^2\neq0$ in general.
- The antisymmetric part of the bilinear form is, with $H=g^{-1}h$ and $H'=g^{-1}h'$,
$$F(X,Y;h,h')-F(X,Y;h',h)=-\sum_ig^{ii}\,\mathrm{tr}\big([B_i,C_i]\,[H',H]\big).$$

(Corrected after review 2026-10-30, R7-1; see `reports/CORRECTIONS.md`. The first version stated L1 for $C_2$, with $8C_2-C_1=28|\mathrm{Ric}|^2$ and $112\Lambda^2$ for Einstein tensors, and drew the conclusion from $C_2=0$, which the leaf condition does not give. The $C_2$ identity is correct but is now only a side identity.)

Remarks on (b):
- The normalisation of $R$ plays no role ($F$ is homogeneous of degree 2), and neither does the sign of $\Lambda$.
- $\Lambda^2$ is a $\mathbb Q$-linear combination of values of the quadratic form, hence of the coefficients of $F(X,X;h,h)$. Output §L4 also certifies this directly from those coefficients. So the ideal they generate contains $\Lambda^2$, and with $\Lambda$ normalised to 1 its Gröbner basis is $\{1\}$; no elimination is needed.
- The argument works because $C_1$ and $C_2^{\rm sym}$ weight $|W|^2$ and $\Lambda^2$ in different ratios, $6:16$ and $\tfrac32:12$. Neither trace suffices alone, since in indefinite signature $|W|^2$ takes both signs.

*Proof.*
- For $B\in\mathfrak{so}(g)$ and $h$ symmetric, put $H=g^{-1}h$, which is $g$-self-adjoint. Then $g^{-1}(B\cdot h)=[B,H]$, so
$$F(X,Y;h,h')=\sum_{ij}g^{ij}\,\mathrm{tr}\big([B_i,H][C_j,H']\big).$$
- *$C_1$.* It uses only the diagonal values $Q_{aa}(\theta^p\odot\theta^q)=F(e_a,e_a;\theta^p\odot\theta^q,\theta^p\odot\theta^q)$.
  - $\mathrm{tr}([B,H][C,H'])=-\mathrm{tr}(H\,\mathrm{ad}_B\mathrm{ad}_CH')$. Summing over $G_0$-dual bases of the $g$-self-adjoint endomorphisms $\mathrm{Sym}_g$ gives $-\mathrm{Tr}_{\mathrm{Sym}_g}(\mathrm{ad}_B\mathrm{ad}_C)$.
  - $\mathrm{Tr}_{\mathfrak{gl}}(\mathrm{ad}_B\mathrm{ad}_C)=2n\,\mathrm{tr}(BC)-2\,\mathrm{tr}B\,\mathrm{tr}C$ and $\mathrm{Tr}_{\mathfrak{so}(g)}(\mathrm{ad}_B\mathrm{ad}_C)=(n-2)\,\mathrm{tr}(BC)$. So the sum is $-(n+2)\,\mathrm{tr}(BC)$.
  - Finally $\sum_{a,i}g^{aa}g^{ii}\,\mathrm{tr}(R(e_a,e_i)R(e_a,e_i))=R_{eqai}R^{qeai}=-|R|^2$.
- *$C_2$ (bilinear).*
  - For vectors $u,v$ let $u\odot v$ be the endomorphism $Z\mapsto\tfrac12(g(v,Z)u+g(u,Z)v)$; it is $g^{-1}(u^\flat\odot v^\flat)$. For $B\in\mathfrak{so}(g)$, $[B,u\odot v]=(Bu)\odot v+u\odot(Bv)$.
  - Since $g^{-1}h_{aq}=\eta_a\eta_q\,e_a\odot e_q$, the sums over $a$ and $b$ factorise:
$$C_2=\sum_{i,q}\eta_i\eta_q\,\mathrm{tr}(U_{iq}^2),\qquad U_{iq}:=\sum_a\eta_a\,[R(e_a,e_i),\,e_a\odot e_q].$$
  - $\sum_a\eta_aR(e_a,e_i)e_a=-\mathrm{Ric}(e_i)^\sharp$, and $\sum_a\eta_a\,e_a\odot(R(e_a,e_i)e_q)$ is $g^{-1}$ of $(Y,Z)\mapsto\tfrac12[R(Y,e_i,e_q,Z)+R(Z,e_i,e_q,Y)]$. So $U_{iq}=g^{-1}u_{iq}$ with
$$u_{iq}(Y,Z)=-\underbrace{\tfrac12\big[\mathrm{Ric}(e_i,Y)g(e_q,Z)+\mathrm{Ric}(e_i,Z)g(e_q,Y)\big]}_{P_{iq}(Y,Z)}+\underbrace{\tfrac12\big[R(Y,e_i,e_q,Z)+R(Z,e_i,e_q,Y)\big]}_{S_{iq}(Y,Z)},$$
    and $C_2=|u|^2=|P|^2-2\langle P,S\rangle+|S|^2$, all norms with $g$.
  - $|P|^2=\tfrac{n+1}2|\mathrm{Ric}|^2$.
  - $\langle P,S\rangle=\sum\mathrm{Ric}^{iy}S_{iq,y}{}^{q}$, and $\sum_qS_{iq}(e_y,e^q)=\tfrac12\sum_qR(e^q,e_i,e_q,e_y)=-\tfrac12\mathrm{Ric}(e_i,e_y)$ (the other term vanishes by skew-symmetry in the last pair). So $\langle P,S\rangle=-\tfrac12|\mathrm{Ric}|^2$.
  - $|S|^2=\tfrac14\big(2|R|^2+2R_{yiqz}R^{ziqy}\big)=\tfrac34|R|^2$, by the standard consequence $R_{abcd}R^{acbd}=\tfrac12|R|^2$ of the first Bianchi identity (pair symmetry and the two skew-symmetries turn $R^{ziqy}$ into $R^{yqiz}$).
  - Adding up: $C_2=\tfrac{n+1}2|\mathrm{Ric}|^2+|\mathrm{Ric}|^2+\tfrac34|R|^2$.
- *The antisymmetric part.*
  - $\mathrm{tr}([B,H][C,H'])=-\mathrm{tr}(H[B,[C,H']])$ and $\mathrm{tr}([B,H'][C,H])=-\mathrm{tr}(H[C,[B,H']])$. By Jacobi their difference is $-\mathrm{tr}(H[[B,C],H'])=-\mathrm{tr}([B,C][H',H])$. This gives (c)'s formula.
  - For $C_2-C_2'$ one needs $\sum_q\eta_q[g^{-1}h_{bq},g^{-1}h_{aq}]=\eta_a\eta_b\,\tfrac{n+2}4\,(e_b\wedge e_a)$, with $(u\wedge v)(Z)=g(v,Z)u-g(u,Z)v$. This follows from $\sum_q\eta_q(e_b\odot e_q)(e_a\odot e_q)=\tfrac14[(n+2)\,e_b\,g(e_a,\cdot)+g_{ab}\,1]$.
  - With $D_{ab}=\sum_i\eta_i[R(e_a,e_i),R(e_b,e_i)]$ and $\mathrm{tr}(K(u\wedge v))=-2g(u,Kv)$ for $K\in\mathfrak{so}(g)$:
$$C_2-C_2'=-\tfrac{n+2}2\sum_{a,b}\eta_a\eta_b\,g(e_a,D_{ab}e_b)=-\tfrac{n+2}2\Big(\tfrac12|R|^2-|\mathrm{Ric}|^2\Big).$$
    Here $\sum_{ab}\eta_a\eta_b\,g(e_a,[R(e_a,e_i),R(e_b,e_i)]e_b)$ summed with $\eta_i$ splits as $-\sum_i\eta_i|\mathrm{Ric}(e_i)|^2$ (from $\sum_a\eta_aR(e_a,e_i)e_a=-\mathrm{Ric}(e_i)^\sharp$) plus $\sum\eta_a\eta_b\eta_i\,g(R(e_b,e_i)e_a,R(e_a,e_i)e_b)=\tfrac12|R|^2$ (Bianchi, as above).
- *$C_2^{\rm sym}$.*
  - $Q_{ab}(h+h')-Q_{ab}(h)-Q_{ab}(h')=F(e_a,e_b;h,h')+F(e_a,e_b;h',h)$, so $C_2^{\rm sym}=\tfrac12(C_2+C_2')=C_2-\tfrac12(C_2-C_2')$.
  - This is $\tfrac34|R|^2+\tfrac{n+3}2|\mathrm{Ric}|^2+\tfrac{n+2}8|R|^2-\tfrac{n+2}4|\mathrm{Ric}|^2$, which is (a).
  - $C_2'=2C_2^{\rm sym}-C_2$ gives the value in (c).
- *Einstein, $n=4$.* $\langle W,R_1\rangle=0$ (Ricci-flatness of $W$), so $|R|^2=|W|^2+\tfrac{24}9\Lambda^2$, and $|\mathrm{Ric}|^2=4\Lambda^2$. ∎

Output §L1 checks (a) and (c) exactly as polynomial identities for the generic algebraic curvature tensor:
- $n=3$ (6 components), signatures (3,0) and (1,2);
- $n=4$ (20 components), signatures (4,0), (1,3), (3,1) and (2,2);
- $n=5$ (50 components), signatures (5,0), (1,4) and (2,3).

The code computes $C_1$ and $C_2^{\rm sym}$ from values $Q_{ab}(h)$ of the quadratic form only. It checks the asymmetry $C_2\neq C_2'$, as an expected result, together with the antisymmetry formula at random rational data in the four $n=4$ signatures. Output §L2 checks (b) for an arbitrary Weyl tensor in all four $n=4$ signatures. For constant curvature $+1$, $C_1=144$ and $C_2^{\rm sym}=108$, consistent with J5(d)'s "fails".

**Corollary L2 (Theorem J5(c) holds in both signatures of the setting).** Assume the hypotheses of J5: Z2, $\lambda\neq\tfrac14$, $(x,y)\neq0$, $\kappa\neq0$, $U$ connected, open and meeting $s(X)$, $t=f(s)$ and $u=h(s)$ on $U$, $\det M(t,u)\neq0$, $T_e$ the slaved field of J1, and $\Upsilon=0$ on $U$. Then there is no such solution, for any $(z,w)$, in either signature of the setting:
- **(14,0)**, from a definite base;
- **(7,7)**, from a base (1,3) or (3,1) (A4) or from a base (2,2) (A5).

The leaf argument (L0, L1) does not depend on the signature, so it applies wherever J5(b) holds. The other Z2 signatures, such as (8,6), (10,4), (9,5), (13,1) and their mirrors (A4), are not covered, because J1–J2, and with them $M$, $\det M$ and J5(b), are not set up there. (Corrected after review 2026-10-30, R7-2; see `reports/CORRECTIONS.md`. The heading had "holds in every signature", and the statement "whatever the signature of $g_\aleph$".)

*Proof.*
1. By Lemma L0, $F\equiv0$ for $R=R_g(x_0)$.
2. By J5(b), $\mathrm{Ric}_g(x_0)=\Lambda_gg$ with $\Lambda_g=-c_H/c_V\neq0$.
3. By Theorem L1, $0=4C_2^{\rm sym}-C_1=32\Lambda_g^2$. This is a contradiction. ∎ (Corrected after review 2026-10-30, R7-1; see `reports/CORRECTIONS.md`. Step 3 had "$0=8C_2-C_1=112\Lambda_g^2$", but $C_2=0$ does not follow from $F\equiv0$.)

The proof uses only the order-$\varepsilon^2$ part of the horizontal–horizontal O'Neill block over the single point $x_0$. It does not use the VV and HV blocks of the leaf's Einstein condition, $\nabla R$, the other leaves, or the higher orders. The sign of $\Lambda_g$, i.e. of $c_Hc_V$ (fixed by A4 in (7,7)), does not matter.

Two independent routes reach the same conclusion. For a Lorentzian or definite base, Proposition L3 needs a single value of the quadratic form. For every base, the span certificate of output §L4 shows that $\Lambda^2$ lies in the span of the coefficients of $F(X,X;h,h)$.

**Proposition L3 (positivity: one value of $F$ for a Lorentzian or definite base).** Let $R$ be any algebraic curvature tensor, $T$ non-null with $T^\perp$ definite, and $e_i$ an orthonormal frame of $T^\perp$. Then
$$F(T,T;\,T^\flat\otimes T^\flat)=2\,g(T,T)\sum_{i,j}\eta_i\eta_j\,R(T,e_i,T,e_j)^2,$$
which is a definite function of $R(T,\cdot,T,\cdot)|_{T^\perp}$, since all $\eta_i\eta_j=1$.
- This is the case for $g$ definite (any $T$) and for $g$ Lorentzian with $T$ timelike.
- So for a base of signature (4,0), (1,3) or (3,1), $F\equiv0$ forces $R=0$, with no Einstein hypothesis.

*Proof.*
- $R(T,e_i)\cdot(T^\flat\otimes T^\flat)=v_i^\flat\otimes T^\flat+T^\flat\otimes v_i^\flat$ with $v_i=R(T,e_i)T\perp T$. Its $G_0$-square is $2g(v_i,v_i)g(T,T)$, and $g(v_i,v_i)=\sum_j\eta_jR(T,e_i,T,e_j)^2$. Summing with $\eta_i$ gives the formula; the term $i=T$ vanishes.
- If $F\equiv0$, then $R(T,X,T,Y)=0$ for all timelike $T$ and all $X,Y$; the components along $T$ drop by the symmetries of $R$.
- This is a polynomial in $T$ vanishing on an open cone, so it vanishes identically.
- The form $R(T,X,T,Y)$ determines $R$ (output §L3: rank 20). ∎

For a Lorentzian base this is a second, one-line proof of Corollary L2. It also sharpens J5(c)'s definite-signature argument, which used $A=0$ on all of $\mathcal W$; here only the $\varepsilon^2$ term is used. In base (2,2), $T^\perp$ has signature (1,2) or (2,1) for every non-null $T$, and no such value is definite. Proposition L4 shows that this is not a defect of the method.

**Proposition L4 (the Einstein solutions of $F\equiv0$, $n=4$).** Let $R$ be Einstein. Then $F\equiv0$ if and only if:
- **base (4,0), (1,3) or (3,1):** $R=0$;
- **base (2,2):** $R_{abcd}=c\,n_{ab}n_{cd}$ with $c\in\mathbb R$ and $n$ a null self-dual or anti-self-dual 2-form. Equivalently: $\Lambda=0$, one of $W^\pm$ vanishes, and the other, as a self-adjoint endomorphism of $\Lambda^\mp$ (induced form of signature (1,2)), has square 0.

These (2,2) tensors also satisfy the full condition: $\det(y)^2\Phi_y(X,Y)=0$ identically in $y$, $X$, $Y$.

*Proof.*
- The (4,0), (1,3) and (3,1) cases are Proposition L3.
- (2,2), "if": $R(X,e_i)=c\,n(X,e_i)\,N$, with $N$ the endomorphism of $n$. So every term of $F$ (quadratic or bilinear) and of $\Phi_y$ carries the factor $\sum_{ij}g^{ij}n(X,e_i)n(Y,e_j)=(n\circ n)(X,Y)$.
  - For an (anti-)self-dual 2-form in dimension 4, $(n\circ n)_{ab}=n_{ac}n_b{}^c=\tfrac14n_{cd}n^{cd}\,g_{ab}$, which is 0 for $n$ null.
  - The same identity gives $\mathrm{Ric}=0$. First Bianchi holds since $n\wedge n=\pm\langle n,n\rangle\,\mathrm{vol}=0$.
- (2,2), "only if", by exact certificates (output §L4). They use only the coefficients of the quadratic form $F(X,X;h,h)$. Here $\mathcal R^\pm=W^\pm+\tfrac\Lambda3$ are the blocks of the curvature operator on $\Lambda^\pm$; $[\mathcal R,*]=0$ and $\mathrm{tr}\,\mathcal R^\pm=\Lambda$.
  1. The 81 products $(\mathcal R^+)_{ij}(\mathcal R^-)_{ab}$ are $\mathbb Q$-linear combinations of those coefficients. So $F\equiv0$ forces $\mathcal R^+=0$ or $\mathcal R^-=0$, and then $\Lambda=0$ by the trace.
  2. On $\mathcal R^+=0$ the coefficients span exactly the span of the entries of $(W^-)^2$ (dimension 6 each), and symmetrically.
  3. A self-adjoint square-zero endomorphism of a nondegenerate 3-space has rank $\le1$, so it is $c\,n\langle n,\cdot\rangle$ with $\langle n,n\rangle=0$. ∎

*Structural explanation (not used in the proofs).*
- Over $\mathbb C$, $V=S^+\otimes S^-$, with $\mathfrak{so}(g)=\mathfrak{sl}(S^+)\oplus\mathfrak{sl}(S^-)$ and $\mathrm{Sym}_{g,0}=\mathfrak{sl}(S^+)\otimes\mathfrak{sl}(S^-)$ as endomorphisms.
- By Schur, $\sum_{ij}g^{ij}(X\wedge e_i)\otimes(X\wedge e_j)$ has $\Lambda^\pm\otimes\Lambda^\pm$ parts proportional to $g(X,X)$ and a $\Lambda^+\otimes\Lambda^-$ part proportional to $(X\otimes X)_0$.
- The part of $F(X,X;S,S)$ that is traceless in $X$ is, for $S=\sum\sigma_a\otimes\tau_a$, the bilinear form $(P,N)\mapsto2\sum_{a,c}\mathrm{tr}(P[\sigma_a,\sigma_c])\,\mathrm{tr}(N[\tau_c,\tau_a])$ evaluated on $(\mathcal R^+\otimes\mathcal R^-)(X\otimes X)_0$. It vanishes for all $S$ only if $\mathcal R^+=0$ or $\mathcal R^-=0$.
- The $g(X,X)$ part then reads $\langle\tau,(\mathrm{tr}\,Q-Q)\tau\rangle=0$ with $Q=(W^\mp)^2$, i.e. $Q=0$.
- In (1,3), $\Lambda^\pm$ are complex conjugate, so $\mathcal R^+=0$ forces $\mathcal R=0$, as in L3.
- The coefficients of $F(X,X;h,h)$ span a 47-dimensional space of quadratic forms in $(W,\Lambda)$, out of 66. These are the 36 products $\mathcal R^+\mathcal R^-$ (with $\Lambda$ shared) and 11 conditions on $(W^+)^2$, $(W^-)^2$.
  - The rank was computed in (1,3), (2,2) and (4,0) (output §L4).
  - (3,1) and (0,4) follow by $g\to-g$, which leaves the endomorphisms $R(X,Y)$ and $G_0$ unchanged and sends $F\to-F$ (and $\Lambda\to-\Lambda$).
  - (Corrected after review 2026-10-30, R7-5; the text had "in every signature".)

**Remark (higher orders and the full condition).** Order $\varepsilon^2$ already fails for every Einstein base with $\Lambda\neq0$, so the orders $\varepsilon^3,\varepsilon^4$ and the full condition on $\mathcal W$ are not needed for J5. The full condition is exactly computable (L0(b)). The Ricci-flat (2,2) tensors of L4 satisfy it to all orders (output §L4: $\det(y)^2\Phi_y\equiv0$ for the example $R=n\otimes n$). So the horizontal leaf condition alone does not force $R=0$ in (2,2). What excludes those tensors in J5 is $\Lambda_g=-c_H/c_V\neq0$ from J5(b).

## J5d.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| L0(a), (c) (output §L0) | exact at 3 random rational $y$ (trace identity, closed form of $\Phi_y$, scale invariance) and at random rational $R$, $h$, $X$, $Y$ (Taylor coefficients $\varepsilon^0,\varepsilon^1,\varepsilon^2$ of $\det(y)^2\Phi_y$); $F(h+cg,h+cg)=F(h,h)$ and $X\leftrightarrow Y$ symmetry of the quadratic form, symbolically | bases (4,0), (1,3), (3,1), (2,2) |
| L0, conventions (§L0) | the three constant-curvature polynomials of `J_slaved_sector.out` §J5 reproduced exactly | (4,0), (1,3), (2,2) |
| L1(a) (§L1) | exact polynomial identities for the generic algebraic curvature tensor; $C_1$, $C_2^{\rm sym}$ computed from values of the quadratic form only | $n=3$: 2 signatures; $n=4$: 4; $n=5$: 3 |
| L1(c) (§L1) | side identities for $C_2$, $C_2'$; $C_2\neq C_2'$ (expected); the antisymmetric part of the bilinear form at random rational data, nonzero | same; antisymmetry: the four $n=4$ signatures |
| L1(b), L2 (§L2) | Einstein forms for an arbitrary Weyl tensor; $4C_2^{\rm sym}-C_1=32\Lambda^2$; the contradiction with $\Lambda_g=-c_H/c_V$ | four $n=4$ signatures; J5(b) in (14,0) and (7,7) |
| L3 (§L3) | the identity for generic $R$ (20 components); rank 20 of $R\mapsto R(T,X,T,Y)$ | (4,0), (1,3), (3,1); (2,2) shown indefinite |
| L4 (§L4) | proof ("if"); exact span certificates from the coefficients of $F(X,X;h,h)$ ("only if"): rank 47 of 66 and $\Lambda^2$ in the span ((1,3), (2,2), (4,0); (3,1), (0,4) by $g\to-g$); the 81 block products in the span (they span 36); the restricted spans equal those of the entries of $(W^\mp)^2$; the null example, including the full condition symbolically in $y$ | $n=4$; (2,2) for the block certificates |
| Review | independent cross-check `xcheck_round7.py` §X4–X6 (REVIEW5-GU), which found R7-1 | L0–L4 |

Run: `cd computation/geometric-unity && python3 J5d_leaf.py > J5d_leaf.out 2>&1`. Result: 57 checks, 57 True, each asserted; 181 s and 170 s on one core (the first version: 44 checks, 102 s).

## J5d.3 Limits (layers b, d)

- **Scope.** Corollary L2 removes the hypothesis "$g_\aleph$ definite" from J5(c) within the setting of J1–J2, i.e. it adds (7,7) to (14,0). All other hypotheses of J5 stay: Z2, the radial ansatz, the slaved $T_e\in\mathcal V$ with $\det M\neq0$, $\kappa\neq0$, $(x,y)\neq0$ and $\lambda\neq\tfrac14$.
- **Other signatures.** The Z2 signatures (8,6), (10,4), (9,5), (13,1) and their mirrors are not covered. The leaf argument would apply there, but J5(b) is not established, since J1–J2 are not set up. (Added after review 2026-10-30, R7-2.)
- **Still not covered** (unchanged from J.3):
  - $\Upsilon=0$ on a neighbourhood beyond the radial Z2 ansatz: Z1, where the Euler field is not Killing and $s$ is not parallel; non-radial $t,u$ (settled for the members of L8 in `L-nonradial-z2.md`);
  - components of $T_e$ outside $\mathcal V$ (round 14, `R-even-operator.md` R5, not yet internally reviewed: excluded under (R), a condition on $(t,u)$; open only on open sets where $(t,u)\in Z_\perp$);
  - $T_o$ outside the invariant sector.
- O'Neill's formulas are used in indefinite signature, as in J5. They are algebraic identities for nondegenerate fibres (L-F1; O'Neill 1983). Here the fibres of $Y'\to X$ carry $c_VG_0$ on traceless forms, which is nondegenerate, and the horizontal block $c_Hg$ is nondegenerate.
- Theorem L1 and Proposition L3 hold for every algebraic curvature tensor; Proposition L4 is stated for Einstein $R$ and $n=4$. J5 itself is stated for $n=4$; L1 would serve any $n$.

## J5d.4 Consequences for gaps and claims

- **J5(c), (d).** J5(c) holds in both signatures of the setting, (14,0) and (7,7). J5(d)'s open question has a negative answer: no Einstein base with $\Lambda_g\neq0$, of any signature and of constant curvature or not, satisfies the order-$\varepsilon^2$ leaf condition. (Corrected after review 2026-10-30, R7-2; the text had "J5(c) holds in every signature".)
  - J5(d)'s sentence "a radial solution near $s(X)$ needs a non-constant-curvature Einstein base with $\Lambda_g=-c_H/c_V$ satisfying this quadratic condition" is true but vacuous. J5(d) asserted openness, not existence, so no statement made before J5d is false.
  - The proof step of Corollary L2 as first committed was invalid (R7-1); its conclusion stands. That correction has its own `reports/CORRECTIONS.md` entry.
  - The summaries that call the (7,7) case open are now outdated.
- **GU-G14.** For the radial Z2 ansatz ($\kappa\neq0$, slaved $T_e\in\mathcal V$, $\det M\neq0$, $(x,y)\neq0$), $\Upsilon=0$ has no solution on a neighbourhood of a point of $s(X)$ in (14,0) or in (7,7), the latter from a base (1,3), (3,1) or (2,2). The open neighbourhood questions are those beyond the radial Z2 ansatz (Z1; non-radial $t,u$, settled for the members of L8 in `L-nonradial-z2.md`).
- **GU-054.** Status unchanged (C).
