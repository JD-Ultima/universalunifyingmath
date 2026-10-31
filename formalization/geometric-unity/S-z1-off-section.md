# S. Υ = 0 in the Zorro metric Z1 on open sets that need not meet $s(X)$

Claims: GU-054 (context). Gaps: GU-G14 (and F9, GU-G16). Code: `computation/geometric-unity/S_z1_off_section.py` → `.out` (sections S0–S5; exact: Python Fractions, also as coefficients of N's truncated power series, and sympy rationals; fixed seeds; 260 s on one core, a second run 262 s, identical apart from the timing annotations). It imports `F_jets.py` and `N_z1_neighbourhood.py` (and through it `L_nonradial.py`, `J_slaved_sector.py`, …), all unmodified. Output sections: §S0 Lemma S1, §S1 Lemma S2 and Corollary S3, §S2 Lemma S4, §S3 Theorem S5 and Lemma S6, §S4 Propositions S7–S9, §S5 Lemma S13 and Theorem S14. Record: `metrics/geometric-unity/S-z1-off-section.yaml`. Builds on N (N1–N4, N6, N7, N9–N11), Q (Q1–Q8, Q12–Q15), P (P1(a), P2(d)), F (F1, F2.4, F4, F6), L (L9) and M (M3, M4). Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## S.0 Question and setting (layer b)

N.3 and Q.3 list as open in Z1 "open sets away from $s(X)$" and "type-(iii) points without analyticity". There L6 is unavailable: the Euler field is not Killing. This note drops the hypothesis $U\cap s(X)\neq\emptyset$.

**Hypotheses (H1′).** (H1) of N.0, for every member with $(x,y)\neq0$ (Q.0), without $U\cap s(X)\neq\emptyset$:
- Z1, $n=4$, $\lambda\neq\tfrac14$, $c_Hc_V\neq0$; (14,0) (base (4,0), $\lambda<\tfrac14$) and (7,7) (base (1,3), (3,1) or (2,2), $\lambda>\tfrac14$);
- $\kappa\neq0$; $T_o=t\Phi^1+u\Phi^1\Gamma$; $T_e=\Sigma(dt,du)^\top$ slaved in $\mathcal V$ with $\det M\neq0$;
- $U\subset Y$ connected and open, $\Upsilon=0$ on $U$.

**Notation.** As in N.0 and Q.0 ($E$, $s$, $P_H$, $V_0$, $c$, $b$, $\beta_S$, $B=4bc$, $\sigma_g$, $Z$, $U^*$, Σ, the types (i)–(iii) of N4, now for every point of $U$). Further:
- $F_x=\pi^{-1}(x)$; $q_\theta:=(x,e^\theta g(x))$, so $\Sigma_x:=\Sigma\cap F_x$ is the $E$-orbit through $g(x)$; $\rho_\theta:=e^{-\theta}\sigma_g/(4b)$.
- $\Pi\subset Y$, the parallel set (S2); $F_x^{\mathbb C}:=\{y\in\mathrm{Sym}^2(\mathbb C^4):\det y\neq0\}$.
- $\mathfrak L_k(q)$, for $q\in Y$ and $k\in\mathbb C$: the linear system $R(X,Y)v-k\big[\langle v,X\rangle\widehat{\mathrm{Ric}}Y-\langle v,Y\rangle\widehat{\mathrm{Ric}}X\big]-\big[\langle n,X\rangle Y-\langle n,Y\rangle X\big]=(\nabla_X\widehat{\mathrm{Ric}})Y-(\nabla_Y\widehat{\mathrm{Ric}})X$ at $q$, unknowns $(v,n)\in T_qY^2$ (28). This is N's 1274×28 system with $\mu_1\to k$, right side divided by $\mu_0$ (Q6). If $(\star)$ $\nabla^2r=\mu(r)\mathrm{Ric}+\nu g_\aleph$ holds near $q$ with $\mu(r(q))\neq0$, then $(\nabla r,\nabla\nu)/\mu(r(q))$ solves $\mathfrak L_{\mu'(r(q))}(q)$.
- $\mathfrak L^{\rm rel}(q)$: the same with $k\langle v,\cdot\rangle$ replaced by a free covector ξ; unknowns $(v,\xi,n)$ (42). A solution of $\mathfrak L_k$ gives the solution $(v,kv^\flat,n)$ of $\mathfrak L^{\rm rel}$.

**Real forms.** As in L.0. Statements about $\mathfrak L$ and the fibre jets are polynomial identities, valid in both real forms.

**Inputs (proven earlier).** N1, N2, N3, N4's argument, N6 ((f) computer-assisted, REVIEW8-GU R10-1), N7's proof, N9 (α), N10, N11; Q1–Q8, Q12, Q13, Q15; P1(a), P2(d); F1, F2.4, F4, F6; L9 (the fibre is connected); M3, M4.

## S.1 Results (layer c)

**Lemma S1 (the Euler field off the section; what replaces L6).** *[proven]* At every point of $Y$:

(a) $[E,X^H]=0$ for basic $X$, $\nabla_EX^H=\tfrac12X^H$, and $\nabla_EU=[E,U]$ for vertical $U$ ($=-\partial_A$ for the coordinate fields).

(b) If $(\star)$ holds on an open $W$ (any μ), then on $W$
$$(\mathrm{S1H})\quad E(X^Hr)=\tfrac12\,X^Hr,\qquad (\mathrm{S1V})\quad d(Er)|_V=\big(c\nu-\mu(r)\big)\,ds|_V .$$
So $X^Hr$ is homogeneous of degree ½ along the $E$-orbits; on each fibre $\nabla^2r(E,U)=0$ for $U\in V_0$, $d(Er)\wedge ds=0$, and $c\nu=E(Er)+\mu(r)$.

(c) $\mathcal L_E\mathrm{Ric}=c_H\partial_{c_H}\mathrm{Ric}$ (N1(b)) is not zero (already on $s(X)$, N6(c)). So $(\star)$ is not $E$-invariant, and L6's conclusion ($Er=0$, $\nu=0$) has no analogue.

(d) The fibre equation alone carries no obstruction. Over a flat base, for every λ, $w=f(s)$ with $f''=\mu_1bc\,f$ solves $\mathrm{Hess}^Fw=\mu_1w\,\mathrm{Ric}_{VV}+\nu g_F$ with $\nu=-\mu_1\beta_Sw$. So round 12's fibre-only method (P6) does not transfer to Z1. What transfers is P2(d)'s dichotomy, applied to the full degree-1 system on $Y$ (S4, S7).

*Proof.*
- (a) The flow of $E$ is $\phi_\theta(y)=e^\theta y$. It commutes with π and preserves $H=\ker D$, since $D$ is linear in $y$. So $[E,X^H]=0$ and $\nabla_EX^H=\nabla_{X^H}E=\tfrac12X^H$ (N1(b)). For vertical $U$, $\nabla_UE=\tfrac12P_HU=0$.
- (b) $(\star)(E,X^H)$: the left side is $E(X^Hr)-\tfrac12X^Hr$; the right side is $-\mu\,ds(X^H)=0$ (N1(d); $ds|_H=0$). $(\star)(U,E)$: the left side is $U(Er)$; the right side is $(c\nu-\mu)\,ds(U)$.
- (d) By N2, $\mathrm{Ric}_{VV}=\beta_Sg_F+bc\,ds^2$ and $g_F=g_{V_0}+c\,ds^2$. Also $\mathrm{Hess}^Ff(s)=f''ds^2$, since $ds$ is parallel on $F$ (F1). ∎

Output §S0 checks (a) at points off the section, over curved bases in general coordinates, and checks (d) symbolically.

**Lemma S2 (the parallel set).** *[proven]* For $q=(x,y)\in Y$ the following are equivalent:
- (i) $q$ lies on the section of a metric $g'$, defined near $x$, with $\nabla^{g'}=\aleph$; then $\mathrm{Z1}(g')=\mathrm{Z1}(g)$ near $q$ (F2.4);
- (ii) $y$ extends to an ℵ-parallel symmetric 2-tensor near $x$;
- (iii) $y$ is annihilated by the local holonomy algebra $\mathfrak{hol}^*_x$ of ℵ.

Let Π be the set of such $q$ and $\Pi_x:=\Pi\cap F_x$. Then $\Pi_x$ is $\mathrm{Met}_{(p,q)}\cap(\mathrm{Sym}^2T_x^*X)^{\mathfrak{hol}^*_x}$, an open subset of a linear subspace, and $\Pi_x\supseteq\Sigma_x$.
- $\Pi_x=\Sigma_x$ (for base (2,2): $\Sigma_x\cup(-\Sigma_x)$) iff the $\mathfrak{hol}^*_x$-invariant symmetric forms are $\mathbb Rg(x)$. This holds when the curvature endomorphisms $R(X,Y)$ at $x$ have only $\mathbb Rg(x)$ as common invariant, in particular for generic curvature (§S1).
- $\Pi_x\supsetneq\Sigma_x$ for local products ($g_1\oplus ag_2$) and for bases with a parallel null 1-form ℓ ($g+a\,\ell\otimes\ell$; pp-waves). For flat ℵ, $\Pi=Y$.

*Proof.*
- (i)⇔(ii): ℵ is torsion-free, so $\nabla^{g'}=\aleph$ iff $\aleph g'=0$.
- (ii)⇔(iii): on a small simply connected neighbourhood $U$ of $x$, parallel tensors correspond to the tensors at $x$ invariant under $\mathrm{Hol}^0(U)$, and $\mathrm{Hol}^0(U)=\mathrm{Hol}^*(x)$ for $U$ small (Kobayashi–Nomizu I, II.10).
- $\mathfrak{hol}^*_x$ contains the $R(X,Y)$ at $x$. The invariant symmetric forms are $g(A\cdot,\cdot)$ with $A$ $g$-self-adjoint in the commutant. In base (2,2), $-g(x)$ has signature (2,2) and lies in $F_x$. ∎

**Corollary S3 (transfer from $s(X)$ to Π).** *[proven]* N12, Q14 and Q16, and every statement of N and Q at points of $s(X)$ or Σ, hold with "$U\cap s(X)\neq\emptyset$" replaced by "$U\cap\Pi\neq\emptyset$". In the conclusions $g$ is replaced by the parallel metric $g'$ through a point of $U\cap\Pi$, with $R_{g'}=R_g$ and $\mathrm{Ric}_{g'}=\mathrm{Ric}_g$.

*Proof.* The fields and operators of the setting are built from $g_\aleph$ and the chimeric bundle with its soldering form $\Theta_\aleph$ ($\Phi^1$ and Γ are $\nabla^0$-parallel, J1). These depend on $g$ only through ℵ (F2.4). So (H1) for $(g,U)$ with $U\cap s_{g'}(X)\neq\emptyset$ is (H1) for $(g',U)$. ∎

Output §S1:
- generic curvature in three signatures has only $\mathbb R\,g$ as invariant form;
- at points of Π over a product base and over a plane-wave base, $\mathrm{Ric}(\mathrm{Z1})$ is F4(c) for $g'$, and at the plane-wave points Z1 is Einstein exactly for $\lambda=\tfrac15$ (F6 transferred);
- the bases $g$ and $3g$ give the same 2-jet of Z1.

**Lemma S4 (fibre polynomiality; the dichotomy).** *[proven; (a) checked]*

(a) Fix $x\in X$ and the 3-jet of ℵ at $x$. In the coordinates $(x,y)$ the components at $(x,y)$ of $g_\aleph$, $g_\aleph^{-1}$, the Christoffel symbols, Riemann, Ric, $s_Y$ and $\nabla\mathrm{Ric}$ of Z1 are polynomials in $(y,(\det y)^{-1})$. (Riemann, Ric and $s_Y$ need only the 2-jet of ℵ.) So the coefficients of $\mathfrak L_k$, $\mathfrak L^{\rm rel}$ and of the Einstein condition on $F_x$ are polynomials in $(y,(\det y)^{-1},k)$.

(b) *Dichotomy (P2(d) for Z1).* Let $A(y)\zeta=b(y)$ be such a linear system on $F_x$, and $r$ the maximal rank of $A$ on $F_x$.
- If the system is consistent on a nonempty open $\Omega\subset F_x$, it is consistent at every $y\in F_x$ with $\mathrm{rank}\,A(y)=r$.
- If $r$ is the number of unknowns, the solution $\zeta(y)$ is a rational function on $F_x$. Every polynomial identity in $y$, ζ and the fibre derivatives of ζ that holds on Ω holds on $F_x$, and on $F_x^{\mathbb C}$ off the poles.

*Proof.*
- (a) $g_\aleph=c_Hy(dx,dx)+c_VG_\lambda|_y(Dy,Dy)$. Here $G_\lambda|_y$ is polynomial in $y^{-1}=\mathrm{adj}(y)/\det y$, and $Dy=dy-\Psi\,dx$ with $\Psi_\mu=\Gamma_\mu^\top y+y\Gamma_\mu$ linear in $y$, the jets of Γ at $x$ being fixed. In the frame of horizontal lifts $\det g_\aleph$ is a nonzero constant times $(\det y)^{-4}$, since $\det G_\lambda|_y\propto(1-4\lambda)(\det y)^{-5}$. Derivatives of $y^{-1}$ are polynomial in $y^{-1}$.
- (b) On Ω, $\mathrm{rank}[A|b]\le r$, so every $(r+1)$-minor of $[A|b]$ vanishes there. These minors are polynomials in $(y,(\det y)^{-1})$, so they vanish on $F_x$, which is connected (L9). Where $\mathrm{rank}\,A=r$, $\mathrm{rank}[A|b]\le\mathrm{rank}\,A$. The rational solution follows from Cramer's rule. ∎

Output §S2: along a line in the fibre over a curved (4,0) base, every component of $(\det y)^3\,\mathrm{Ric}(\mathrm{Z1})$ is a polynomial of degree ≤ 13 (exact at 20 nodes).

**Theorem S5 (Einstein pieces; Conjecture F9 for Z1 at $n=4$).** Let $W\subset Y$ be open and Z1 Einstein on $W$.

(a) *[proven]* For every $x\in\pi(W)$, Z1 is Einstein at every point of $F_x$. Hence $\lambda=\tfrac15$ and $\mathrm{Ric}_g=0$ on $\pi(W)$. In particular, in (7,7) Z1 is Einstein on no open set.

(b) *[proven; α computed (N9)]* Over a definite base, $g$ is flat on $\pi(W)$.

(c) *[proven, computer-assisted: S6 and two exact evaluations]* Over a base (1,3) or (3,1), $g$ is flat on $\pi(W)$.

(d) *[proven, computer-assisted: S6 and sixteen exact evaluations]* Over a base (2,2), $g$ is flat on $\pi(W)$.

So F9 holds for Z1 at $n=4$, on every open set and over every base: Z1 is Einstein on an open set $W$ iff $\lambda=\tfrac15$ and ℵ is flat on $\pi(W)$ (the converse is N2).

*Proof.*
- (a) The Einstein constant is $-1/c$ (N1(d)). On $F_x$ the components of $\mathrm{Ric}+g_\aleph/c$ are polynomials in $(y,(\det y)^{-1})$ (S4(a)) that vanish on the open set $W\cap F_x$, hence on $F_x$. At $g(x)$, F6 applies.
- (b) Z1 is Einstein along $F_x$ near $p=g(x)$, so $\hat{\mathcal A}\equiv0$ on the slice through $p$. N10's proof uses only this and $\mathrm{Ric}_g(x)=0$: $\alpha|W|^2=0$ with $\alpha=\tfrac{45}{8c_H^2c_V}\neq0$.
- (c), (d), step 1.
  - By S6 the order-2 fibre coefficient at $p$ of $\mathrm{Ric}_{VV}+\tfrac54G_{1/5}$ is that of $\mathcal A_x=\mathcal A_{\mathfrak b_W}$, so it vanishes.
  - That coefficient is $\mathfrak b^{ab}G_\eta(\cdot,\xi_{K_a}(h))\,G_\eta(\cdot,\xi_{K_b}(h))$. For generic $h$ the six covectors $G_\eta(\cdot,\xi_K(h))$, $K\in\mathfrak{so}(\eta)$, are independent (§S3, four signatures).
  - So $\mathfrak b_W=0$, i.e. $\mathcal R_W^2=0$ (P1(a)).
- (c), step 2. On $\Lambda^2\cong\mathbb C^3$ (for Ricci-flat $R$), $\mathcal R=W$ is complex symmetric and traceless.
  - $W^2=0$ makes $\mathrm{Im}\,W$ totally isotropic, so $W=\psi\,o^4$ for a spinor $o$: Petrov type N.
  - $SL(2,\mathbb C)$ is transitive on nonzero spinors, and $o\mapsto ao$ multiplies ψ by $a^4$. So the type-N tensors form one $SO_0(1,3)$-orbit, that of the plane-wave tensor $R_N$.
  - The condition "Einstein along all of $F_x$" is invariant under:
    - base isometries (naturality);
    - $(c_H,c_V)\mapsto(tc_H,tc_V)$;
    - the fibre scaling (N1(b): it maps $F_x$ to itself and multiplies $c_H$ by $e^\theta$);
    - $(c_H,c_V)\mapsto(-c_H,-c_V)$ (Q1(b)).
  - So the signs $c_H/c_V>0$ and $<0$ decide it. Output §S3: for $R_N$ the VV block of $\mathrm{Ric}+g_\aleph/c$ is nonzero at a fibre point in both cases. Hence $W=0$.
  - (3,1): $y\mapsto-y$ is an isometry from $\mathrm{Met}_{(3,1)}$ with $(c_H,c_V)$ to $\mathrm{Met}_{(1,3)}$ with $(-c_H,c_V)$, with the same ℵ.
- (d), step 2. $\Lambda^\pm$ have signature (1,2), and $W=Q_+\oplus Q_-$ with $Q_\pm^2=0$.
  - $Q_\pm$ is self-adjoint, so $\mathrm{Im}\,Q_\pm=(\ker Q_\pm)^\perp\supseteq\mathrm{Im}\,Q_\pm$ is totally null, of dimension ≤ 1.
  - So $Q_\pm=s_\pm\,\alpha_\pm\otimes\alpha_\pm^\flat$ with $\alpha_\pm=\ell\wedge m_\pm$ null and decomposable, spanning totally null planes of the two families.
  - $SO_0(2,2)$ acts through $SO_0(1,2)\times SO_0(1,2)$, transitively on each null cone up to positive scaling. So the nonzero null $W$ form eight orbits, $s_\pm\in\{-1,0,1\}$, not both 0.
  - Output §S3 evaluates a representative of each, for $c_V=\pm3$: none is Einstein at the fibre point. Hence $W=0$. ∎

Over a plane-wave base Z1 is therefore Einstein (at $\lambda=\tfrac15$) exactly on the 6-dimensional set Π (S2, §S1), and on no open set. In Z2 the A-term vanishes identically there (P5); in Z1 the fibre correction does not.

**Lemma S6 (N5 at order 2 in the VV block, for Weyl tensors).** *[proven, computer-assisted]* Let the base be Ricci-flat at $x_0$, with Weyl tensor $W$, in any signature, and $c_Hc_V\neq0$. Then the order-2 Taylor coefficient in $h=y-g(x_0)$ of $\mathrm{Ric}_{VV}(\mathrm{Z1})$ along $F_{x_0}$ equals that of $-\tfrac54G_{1/5}|_y+\mathcal A_{x_0}$, at $\lambda=\tfrac15$.

*Proof.*
- Let Δ be the order-≤2 part of $\mathrm{Ric}_{VV}(\mathrm{Z1})-\mathrm{Ric}_{VV}(\mathrm{Z2})-(-\tfrac54G_{1/5}+G_{1/4})$, where $\mathrm{Ric}_{VV}(\mathrm{Z2})=-G_{1/4}+\mathcal A$ (O'Neill; $\mathcal A$ given by L.0's formula). The argument below covers orders 0, 1 and 2 together, so N6(f) is not used.
- Along $F_{x_0}$, Ric has degree ≤ 2 in $R_g(x_0)$, and its VV block does not depend on $\nabla R_g(x_0)$ (N9 step 5). So $\Delta=\Delta_1+\Delta_2$, linear plus quadratic in $W$.
- Δ is $O(\eta)$-equivariant (naturality) and polynomial, hence $SO(4,\mathbb C)$-equivariant.
- Over ℂ, $\mathrm{Weyl}=(5,1)\oplus(1,5)$, so $\Delta_1=0$ as soon as $\Delta_1(W_1)=0$ for one $W_1$ with $W_1^\pm\neq0$ (Schur, as in N6(f)).
- $\mathrm{Sym}^2\mathrm{Weyl}=[(9,1)+(5,1)+(1,1)]\oplus[(1,9)+(1,5)+(1,1)]\oplus(5,5)$. Only $(1,1)$ repeats, twice.
- For $W^\pm\neq0$ and $W$ not null, the projections of $W\otimes W$ onto $(9,1)$, $(5,1)$, $(5,5)$ and their mirrors are nonzero:
  - $(9,1)$: a product of nonzero binary quartics;
  - $(5,1)$: $(Q_+^2)_0$, which vanishes for a traceless $3\times3$ matrix only if $Q_+^2=0$;
  - $(5,5)$: $W^+\otimes W^-$.
  Each of these generates its summand. The projection onto $2(1,1)$ is $(\langle W^+,W^+\rangle,\langle W^-,W^-\rangle)$.
- So $\Delta_2=0$ as soon as $\Delta_2(W_1)=\Delta_2(W_2)=0$ for two such tensors with independent invariant pairs. The values $\Delta(W_1)$, $\Delta(2W_1)$, $\Delta(W_2)$ separate $\Delta_1$ from $\Delta_2$.
- Output §S3 evaluates them in bases (1,3) and (2,2) at $(c_H,c_V)=(2,3)$, with $\nabla R=0$, and checks the Schur data. (4,0) is not needed: S5(b) uses α.
- The fibre scaling and the transfer to Σ carry the identity (for all $W$) to every $c_H/c_V>0$. Since the coefficients are rational in $c_H/c_V$, it then holds for every value. ∎

**Proposition S7 (linear μ: the degree-1 conditions continue to Σ).** *[proven; uses N6(f)]* Let $W\subset Y$ be open, $\mu_1\neq0$, and $r$ smooth and nonzero on $W$ with $\nabla^2r=\mu_1r\,\mathrm{Ric}+\nu g_\aleph$ on $W$. Then for every $x\in\pi(W)$:
- $\mathrm{Ric}_g(x)=0$ and $b(4\mu_1bc-1)=0$;
- the unique solution $\zeta=(v',n')$ of $\mathfrak L_{\mu_1}$ is rational on $F_x$;
- on a dense open subset of $W\cap F_x$ it equals $(\nabla r,\nabla\nu)/(\mu_1r)$;
- at every $q_\theta$ it equals $(2bE,-2\mu_1\beta_Sb\,E)$.

*Proof.*
1. $\mathfrak L_{\mu_1}(q)$ is consistent at every $q\in W\cap F_x$.
2. At $q_\theta$, a point of $s(X)$ for $e^\theta g$, the coefficient matrix has rank 28 for every base. Indeed N7's steps 1–4 for the homogeneous system give $v_H=n_H=0$, $v_0=n_0=0$, then $v_E=0$ (step 3 with right side 0) and $n_E=0$ (output §S2, §S4).
3. So the maximal rank on $F_x$ is 28, and S4(b) gives consistency at every $q_\theta$, $\theta\in\mathbb R$.
4. N7(a)–(c) (Q6 with $k=\mu_1$) are consequences of the consistency of $\mathfrak L_{\mu_1}(q_\theta)$. They give an Einstein base at $x$ and $\mu_1c\,e^{-2\theta}\sigma_g^2=4b(4\mu_1bc-1)$ for all θ. Hence $\sigma_g=0$, $b(4\mu_1bc-1)=0$, and the solution is Q6(a)'s with $\sigma=0$. ∎

**Proposition S8 ($\lambda=\tfrac{13}{56}$ is excluded on every open set).** *[proven]* Under S7 with $\mu_1=-\tfrac7{18}$ and the trace equation (Tr) of N3(d) on $W$, $\lambda\neq\tfrac{13}{56}$.

*Proof.* Suppose $4\mu_1bc=1$ and fix $x\in\pi(W)$. Then $\mathrm{Ric}_g(x)=0$ (S7).
1. *A rational ν.* Put $\nu'(y):=\big(\langle\nabla_Ev',E\rangle+1+\mu_1\langle v',E\rangle^2\big)/c$, a rational function on $F_x$.
   - On the solution $\nabla v'=\widehat{\mathrm{Ric}}+\tfrac\nu{\mu_1r}\mathrm{I}-\mu_1v'\otimes v'^\flat$, so $\nu'=\nu/(\mu_1r)$ there.
   - Hence, for vertical $U$, $(\mathrm H)\ \nabla_Uv'=\widehat{\mathrm{Ric}}\,U+\nu'U-\mu_1\langle v',U\rangle v'$ holds on $W\cap F_x$, so on $F_x$ (S4(b)).
2. *(Tr) as a relation between $r$ and $y$.* With $\rho(y):=13\big[\tfrac7{468}s_Y+\tfrac3{28}\mu_1^2|v'|^2-\mu_1\nu'\big]$, (Tr) reads $P_0(r)=\rho$. Here $P_0(r):=\Phi(r)/r=Ar^{-2}+a_0r^{-1}+a_1$ with $A\neq0$ (N3(d)). On a ball $B\subset W\cap F_x$ where $v'$ is regular: $d\log r=\mu_1v'^\flat$ and $P_0(r)=\rho$.
3. *Continuation.* $F_x^{\mathbb C}$ is connected, and $v'$, ρ are holomorphic off a proper analytic subset $Z$, whose complement is connected.
   - Continue along a path γ in $F_x^{\mathbb C}\setminus Z$ from $B$ to $q_{\theta_0}$, which is regular (rank 28): $r_\gamma=r(y_1)\exp(\mu_1\int_\gamma v'^\flat)$ is holomorphic and nonzero.
   - $P_0(r_\gamma)-\rho$ vanishes near the start, so along γ.
4. *Along the orbit.*
   - $v'=2bE$, so $\tfrac d{d\theta}\log r_\gamma(q_\theta)=2\mu_1bc=\tfrac12$.
   - $\nabla_Ev'=2b\nabla_EE=0$, so $\nu'(q_\theta)=(1+4\mu_1b^2c^2)/c$, which is N8's value.
   - $s_Y(q_\theta)=\sigma_1/c_V$ (F4(b)) and $|v'|^2=4b^2c$.
   - So ρ is a constant $\rho_\Sigma$ along the orbit, and the Laurent polynomial $P_0-\rho_\Sigma$ vanishes on the arc $r_\gamma(q_{\theta_0})e^{(\theta-\theta_0)/2}$. Hence $A=0$, a contradiction. ∎

Output §S4 gives the constants: $d\log r/d\theta=\tfrac12$, $\nu'$ equal to N8's value, and $r^2(P_0-\rho_\Sigma)$ with constant term $A=13/(269746176N^3)$.

**Proposition S9 ($\lambda=\tfrac15$, linear μ, definite base).** *[proven; uses N6(f); α (N9)]* Under S7 in (14,0) with $\lambda=\tfrac15$, $R_g=0$ on $\pi(W)$.

*Proof.* Fix $x$ and put $p=g(x)$.
1. By S7, $\mathrm{Ric}_g(x)=0$ and $v'=0$ on $\Sigma_x$. (H) holds on $F_x$ near $p$ (S8 step 1), and $v'^\flat|_{F_x}$ is closed, since the right side of (H) is symmetric.
2. Let $w$ solve $d\log w=\mu_1v'^\flat$ on $F_x$ near $p$ with $w(p)=1$. As $F_x$ is totally geodesic, (H) gives $\mathrm{Hess}^Fw=\mu_1w(\mathrm{Ric}_{VV}+\nu'g_F)$.
3. At $p$: $\nabla w(p)=0$, $\nu'(p)=1/c$ (since $\nabla_Ev'=0$ on Σ) and $\mathrm{Ric}(p)=-g/c$ (F6). So $\mathrm{Hess}^Fw(p)=0$, and $\psi:=\mu_1w(\kappa_E+\nu')$ vanishes at $p$.
4. On the slice $S$ through $p$: $\mathrm{Hess}^Sw=\psi g_S+\mu_1w\hat{\mathcal A}$ with $\hat{\mathcal A}(p)=0$ and $\nabla\hat{\mathcal A}(p)=0$ (N6(f)).
5. N9 steps 2–6 apply to $w$: $\langle q,T'\rangle=\alpha|W|^2=0$ with $\alpha\neq0$. So $W(x)=0$. ∎

**Theorem S10 (the members with a linear line equation, and the Einstein classes, on every open set).** *[proven; uses N6(f) and α]* Assume (H1′) with member $\circledast_\sharp$, $xW_A$ ($x\neq0$), class (A) or (B) of L8, (C2) or (C3) (Q.0).

(a) (7,7): **no solution** on $U$.

(b) (14,0):
- $\lambda=\tfrac15$, $g$ is flat on $\pi(U)$, $(t,u)$ is constant on $U$, and $T_e\equiv0$.
- So the solutions are exactly those of Q15. For $c_V>0$ they exist iff $c_V\kappa^2\ge9360(x'^2+y'^2)$, with strict inequality when $a^2=\tfrac{125}{576}$.
- No analyticity is assumed.

*Proof.*
1. *Einstein pieces.* Every $q\in U$ is of type (i), (ii) or (iii) (N4's argument at points of $U$). The components of $\mathrm{int}\,Z$ are Einstein (N4; Q4; Q13 step 1). By S5: in (7,7), $\mathrm{int}\,Z=\emptyset$; in (14,0), $\lambda=\tfrac15$ and $g$ is flat on $\pi(\mathrm{int}\,Z)$.
2. *No type-(ii) points.* Class (A) and (C2) have $U^*=\emptyset$ (N3(c), Q4). At a type-(ii) point of the other members, on a neighbourhood $N$, $(\star)$ holds with linear μ and $r\neq0$, and $dr\neq0$ on the dense set $U^*\cap N$. The cases are:
   - $\circledast_\sharp$, $xW_A$, (B): N4, with $\mu_1=-\tfrac7{18}$ and (Tr);
   - (C3): Q13 steps 2–3, with $\mu_1=-\tfrac5{12}$.

   S7 then gives $\mathrm{Ric}_g=0$ on $\pi(N)$ and $\lambda\in\{\tfrac15,\tfrac{13}{56}\}$ (resp. $\{\tfrac15,\tfrac{23}{100}\}$), all below $\tfrac14$. So there are none in (7,7). In (14,0):
   - $\tfrac{13}{56}$ is excluded by S8;
   - at $\tfrac15$, S9 gives $R_g=0$ on $\pi(N)$, and N11 excludes $(\star)$, (Tr) with $dr\neq0$ over a flat base.
3. *(7,7).* Every point would be of type (ii), so $U=\emptyset$.
4. *(14,0).*
   - $\mathrm{int}\,Z$ is dense in $U$, so $\pi(\mathrm{int}\,Z)$ is dense in $\pi(U)$, and $g$ is flat on $\pi(U)$.
   - $U^*$ is open and disjoint from $\mathrm{int}\,Z$, so it is empty. So $(t,u)$ is locally constant, hence constant, and $T_e=\Sigma(dt,du)=0$.
   - Over the flat base Z1 is Einstein (N2), and $\Upsilon=0$ becomes the algebra of Q15. ∎

**Corollary S11 (type-(iii) points; analyticity).** *[proven]*
- In N12(b) and Q14(b), every type-(iii) point of $U\cap s(X)$ lies over a point where $R_g=0$: it is a limit of Einstein pieces, which lie over flat regions (S5).
- For the members of S10 the analyticity hypothesis of N12(b) and Q16(b) is not needed.
- For (C1) near $s(X)$, see S12.

**Proposition S12 ((C1) near $s(X)$, without analyticity).** *[proven; uses the inputs of Q14(b)]* Assume (H1) (so $U\cap s(X)\neq\emptyset$), (C1) and (14,0). Then:
- $g$ is flat on the open set $O:=\pi(U\cap s(X))$;
- on $U_1:=U\cap\pi^{-1}(O)$, a neighbourhood of $U\cap s(X)$, $(t,u)$ is locally constant and $T_e\equiv0$.

*Proof.*
1. Q14(b): $\lambda=\tfrac15$, and the points of $U\cap s(X)$ are of types (i) and (iii).
2. Type-(i) points have a neighbourhood over a flat base (S5(b)). Type-(iii) points lie over $\overline{\pi(\mathrm{int}\,Z)}$, where $R_g=0$. So $g$ is flat on $O$, and Z1 is Einstein on $\pi^{-1}(O)$ (N2).
3. A type-(ii) point of $U_1$ would give, on a neighbourhood $N\subset U_1$, $(\star)$ with $\mu\neq0$ on $N\cap\{H(t)\neq0\}$ and $dr\neq0$ on a nonempty open subset. Q12 excludes this, since $b=0$.
4. So $\mathrm{int}\,Z$ is dense in $U_1$, and $U^*\cap U_1=\emptyset$. ∎

**Lemma S13 (the relaxed system at Σ).** *[proven; ranks checked]* At $q_\theta$, for every base:

(a) *Kernel.* For $b\neq0$, the homogeneous $\mathfrak L^{\rm rel}(q_\theta)$ has only the trivial solution unless $\mathrm{Ric}_g(x)=-c_Hb\,e^\theta g(x)$. For $b=0$, it has only the trivial solution if $\mathrm{Ric}_g(x)$ has rank ≥ 2.

(b) *Consistency.* If $\mathfrak L^{\rm rel}(q_\theta)$ is consistent, then $g$ is Einstein at $x$, and every solution has $v=v_EE$ with $v_E=2b-e^{-\theta}\sigma_g/2$ ($=2b(1-\rho_\theta)$ for $b\neq0$).

*Proof.* Follow N7's steps with ξ free. Put $m:=\xi/c-n$ and $\rho:=\mathrm{Ric}_g^\sharp/c_H$ (for $e^\theta g$), and use the data N1(e), N6 at $q_\theta$.
1. *$E$-components.* The right side vanishes (N6(d); $\mathrm{Ric}_{HV}=0$).
   - The pair $(X,U)$ gives $v_H=0$.
   - The pair $(X,E)$ gives $m_H=0$.
   - The pair $(U,E)$ gives $m_0=0$.
2. *$(U,W)\subset V_0$, output $V_0$* (right side 0, N6(f)). $R^S(U,W)v_0=m'(U)W-m'(W)U$ with $m':=\beta_S\xi_0+n_0$. Commuting $U,W$ give $m'=0$, and the trivial commutant gives $v_0=0$. With $m_0=0$: $b\,\xi_0=0$.
3. *$(X,U)$, output $H$.* $-\tfrac{c_H}4v_EU(X,Z)+\xi_0(U)\mathrm{Ric}_g(X,Z)=[\text{right side}]=-\tfrac{c_Hb}2U(X,Z)+\tfrac12\mathrm{Ric}_g(X,\tilde UZ)$.
   - The antisymmetric part in $(X,Z)$ forces an Einstein base. The $U$-part then gives $v_E=2b-\sigma/2$, and the trace part gives $\xi_0\sigma=0$.
   - Homogeneous case: $v_E=0$, and $\xi_0=0$ if $\mathrm{Ric}_g\neq0$.
4. *$(X,U)$, output $V_0$, homogeneous:* $b\,\xi_H=0$.
5. *$(X,Y)$, output $H$, homogeneous:* $\xi_H(X)\rho Y=\xi_H(Y)\rho X$, which forces $\xi_H=0$ if $\mathrm{rank}\,\rho\ge2$.
6. *$(X,E)$, output $H$, with $(E,U)$, output $V_0$:* $\xi(E)(\rho+b)=0$.
7. $n$ is then determined. ∎

Output §S5 gives:
- rank 42 at section and radial points for random bases ((1,3) at $\lambda=\tfrac12$, (4,0) at $\tfrac17$ and $\tfrac15$), where the system is inconsistent: these bases are not Einstein;
- rank 41 at $\mathrm{Ric}_g=-c_Hbg$, and at $\lambda=\tfrac15$ with $\mathrm{Ric}_g$ of rank 1;
- rank 28 at $\lambda=\tfrac15$ with $\mathrm{Ric}_g=0$;
- for constant-curvature bases, the unique solution, which is Q6's: $\xi=k\,v'^\flat$ with $kc\sigma^2=4b(4kbc-1)$;
- generic rank 42, and inconsistency, off the section.

**Theorem S14 ((C1) away from $s(X)$, $\lambda\neq\tfrac15$).** *[proven; uses N6(f), S13 and the certificate of Q8]* Assume (H1′), (C1) and $\lambda\neq\tfrac15$. Then $U$ has no type-(ii) point. Hence in (7,7) there is **no solution**, and in (14,0) there is no solution unless $\lambda=\tfrac15$.

*Proof.* Let $q$ be of type (ii), with $N$ as in Q5(c): $u=0$ and $(E_t)$ hold on $N$. Let Ω be a connected open subset of $N\cap U^*\cap\{H(t)\neq0\}$ on which $\mathfrak L^{\rm rel}$ has rank 42. There:
- $r=F(t)$, $\mu=xF'/H\neq0$ and $\mu'(r)=\kappa(\tau)=-N_2/(18D_1^2)$ (M4, Q.0);
- $(\nabla r,\mu'\nabla r,\nabla\nu)/\mu$ solves $\mathfrak L^{\rm rel}$.

Fix $x\in\pi(\Omega)$; $b\neq0$.
1. *Continuation to Σ.*
   - By S13(a) the rank is 42 at all $q_\theta$ but at most one. So by S4(b) the solution $\zeta=(v',\xi,n')$ is rational on $F_x$, consistent at those $q_\theta$, and equal on $\Omega\cap F_x$ to the solution above.
   - By S13(b), $g$ is Einstein at $x$ and $v'(q_\theta)=v_E(\theta)E$ with $v_E=2b(1-\rho_\theta)$.
2. *The coefficient $k$.*
   - On Ω, $\xi=\mu'(r)v'^\flat$. So the 2×2 minors of $(\xi,v'^\flat)$ vanish on $F_x$, and $k:=\xi(Z)/v'^\flat(Z)$, for any fixed $Z$ with $v'^\flat(Z)\not\equiv0$, is a rational function on $F_x$ equal to $\kappa(\tau)$ on Ω.
   - Where $v_E(\theta)\neq0$, $\xi=k\,v'^\flat$ at $q_\theta$, so $(v',n')$ solves $\mathfrak L_{k(\theta)}(q_\theta)$. Q6(c) then gives $(\mathrm C')\ Bk(1-\rho_\theta^2)=1$.
3. *An algebraic $t$.*
   - On $\Omega\cap F_x$, τ is a root of the quadratic $18k(y)D_1(\tau)^2+N_2(\tau)=0$ (§S5): a branch of an algebraic function on $F_x$.
   - $(H(t)/x)\,d_Ft=v'^\flat|_F$ there, since $dr/\mu=(H/x)dt$. The form $v'^\flat|_F$ is closed (it is exact on Ω).
   - Let $m$ be its continued primitive and $\mathcal G:=\int H\,dt/x=576\int D_1/Q_2\,d\tau$. $\mathcal G$ has logarithmic singularities at the zeros of $Q_2$ (where $\det M=0$) and, for $\epsilon A\neq1$, at $\tau=\infty$; the residue sum is $306(\epsilon A-1)/(144\epsilon A-119)$ (§S5).
   - The relations $k=\kappa(\tau)$, $(H/x)dt=v'^\flat$ and $\mathcal G(t)-m=\mathrm{const}$ continue along every path in $F_x^{\mathbb C}$ that avoids the branch and pole locus, a proper algebraic subset with connected complement.
4. *Along the orbit*, for the continued branch:
   $$(\mathrm B)\quad \frac{H(t)}x\,\frac{dt}{d\theta}=v'^\flat(E)=c\,v_E=\tfrac B2(1-\rho_\theta),\qquad k(\theta)=\kappa(\tau(\theta)).$$
   - *$\sigma_g\neq0$.* $k(\theta)=1/(B(1-\rho_\theta^2))$ is nonconstant. So, off finitely many θ, the branch is analytic and nonconstant.
     - Differentiating (C′) with (B) gives Q7's (D), $\Phi(1-\rho)+4\rho^2=0$, with $\Phi=(x\kappa'/H)/\kappa^2=\mu\mu''/\mu'^2$ (§S5, for any complex branch).
     - So $\tilde R(\tau(\theta))=0$ at infinitely many distinct τ, and $\tilde R\equiv0$ in τ.
     - Q8's certificate is an algebraic statement in $(A,B)$. It leaves only $(A,B)\in\{(0,-\tfrac{18}7),(\tfrac{11}6\epsilon,-\tfrac{289}3)\}$, where $\rho\equiv0$, against $\sigma_g\neq0$ ($A=a^2>0$).
   - *$\sigma_g=0$.* $k(\theta)\equiv1/B$, so $t(q_\theta)$ is a continuous function with values in the finite set $\kappa^{-1}(1/B)\cup\{\infty\}$, hence constant, $=t_*$.
     - $t_*$ is not a logarithmic singularity of $\mathcal G$: there $|\mathcal G|\to\infty$, while $m$ is finite near $q_\theta$, where $v'$ is regular.
     - So $\mathcal G(t(q_\theta))$, and with it $m(q_\theta)$, is constant along the orbit. But $dm/d\theta=v'^\flat(E)=B/2\neq0$.
5. *Conclusion.* In (7,7), $b\neq0$ always, and $\mathrm{int}\,Z=\emptyset$ (S5(a)), so $U=\emptyset$. In (14,0) with $\lambda\neq\tfrac15$, likewise. ∎

**Proposition S15 ((C1), $\lambda=\tfrac15$, away from $s(X)$).** *[proven; partial]* Assume (H1′), (C1), (14,0) and $\lambda=\tfrac15$. Then:
- $g$ is flat on $\pi(\mathrm{int}\,Z)$;
- no type-(ii) point lies in Π (Q14(b) through S3);
- at every type-(ii) point $q$, $\mathrm{Ric}_g(\pi(q))$ has rank ≤ 1.

*Proof of the last item.* If $\mathrm{rank}\,\mathrm{Ric}_g(x)\ge2$, S13(a) gives rank 42 at every $q_\theta$. As in S14:
- $g$ is Einstein at $x$ with $\sigma_g\neq0$;
- (C) reads $k(\theta)c\,e^{-2\theta}\sigma_g^2=0$, so $k\equiv0$ along the orbit;
- $t(q_\theta)$ is then constant, with values in $\kappa^{-1}(0)\cup\{\infty\}$, while $dm/d\theta=c\,v_E=-ce^{-\theta}\sigma_g/2\neq0$. ∎

**Corollary S16 (Z1 under (H1′): the family on every open set).** *[proven; uses N6(f), α, α″, the certificate of Q8]*

(a) **(7,7):** for every member with $(x,y)\neq0$, $\Upsilon=0$ has no solution on any nonempty open $U\subset Y$. With Q3: $\Upsilon=0$ holds on a nonempty open set iff $(x,y)=0$ and $T=0$. This is the same as in Z2 (P7).

(b) **(14,0):**
- For every member with $(x,y)\neq0$, $\lambda=\tfrac15$.
- For every member other than the twisted cone (C1), the solutions on a connected open $U$ are exactly those of Q15: flat base on $\pi(U)$, constant $(t,u)$, $T_e\equiv0$. They exist iff $c_V\kappa^2\ge9360(x'^2+y'^2)$ ($c_V>0$), with strict inequality when $a^2=\tfrac{125}{576}$.
- For (C1): the same holds near $U\cap\Pi$ (S12 with S3), and away from Π S15 applies.

(c) **$(x,y)=0$:** only $T=0$ (Q3).

*Proof.* M10's classification as in Q16: S10 for L8(A), L8(B), the multiples of $W_A$, (C2) and (C3); S14 and S15 for (C1); Q3. ∎

**Answers.**
1. *What replaces L6.*
   - Not $Er=0$, $\nu=0$: $(\star)(E,\cdot)$ only gives $E(X^Hr)=\tfrac12X^Hr$ and $d(Er)|_V=(c\nu-\mu)ds|_V$ (S1).
   - The working replacement is the continuation of the pointwise degree-1 conditions along the fibre to the radial points, where N7/Q6 apply (S4, S7, S13). The Euler field enters there through $\mathrm{Ric}(E,\cdot)=-ds$ and $\nabla E=\tfrac12P_H$.
2. *Which points lie on sections with the same connection.*
   - Exactly the points $(x,y)$ with $y$ invariant under the local holonomy algebra (S2).
   - For generic curvature these are Σ (and $-\Sigma$ in base (2,2)). Products and pp-waves have more.
   - All of N and Q transfer to open sets meeting Π (S3).
3. *Round 12's method in Z1.*
   - The fibre equation alone has solutions over the flat model for every λ (S1(d)). So no fibre-only obstruction exists.
   - P2(d)'s dichotomy does apply to the full degree-1 system on $Y$. Its coefficients are polynomial in $(y,(\det y)^{-1})$ along a whole fibre (S4); the Z1 Ricci tensor along a fibre is algebraic, with $(\det y)^3\mathrm{Ric}$ polynomial along lines.
   - It carries consistency from any open subset of a fibre to Σ (S7, S13). For linear μ this closes everything (S8–S10). For the twisted cone it closes $\lambda\neq\tfrac15$ (S14).
   - N5's order-2 VV statement is proven for Weyl tensors (S6).
4. *Einstein pieces (F9).* Z1 at $n=4$ is Einstein on an open set iff $\lambda=\tfrac15$ and ℵ is flat over it, for every base signature (S5). In (7,7) never.
5. *Solutions.* No new solution: the flat-base, constant-$(t,u)$ solutions of Q15 are the only ones, on every open set, except possibly for (C1) at type-(ii) points away from Π at $\lambda=\tfrac15$ (S16).

## S.2 What is verified, precisely

| Statement | How (section of `S_z1_off_section.out`) | Scope |
|---|---|---|
| S1 | proof; §S0: $\nabla_EX^H=\tfrac12X^H$, $\nabla_E\partial_A=-\partial_A$ by brute-force 2-jets off the section, cubic random bases in general coordinates ($\Gamma(x_0)\neq0$); (d) symbolic in λ, $c_V$ | (4,0) λ=1/7; (1,3) λ=½ |
| S2, S3 | proof; §S1: invariant forms of the curvature endomorphisms (dimension 1 for random $R$; 2 for a product base and a plane-wave base); Ric(Z1) at points of Π equal to F4(c) for $g'$ (control off Π fails); Einstein at those plane-wave points exactly at λ=1/5; base $3g$ gives the same Z1 | (4,0), (1,3), (2,2); product (4,0); plane wave (1,3) |
| S4 | proof; §S2: $(\det y)^3\mathrm{Ric}(\mathrm{Z1})$ polynomial of degree ≤ 13 along a fibre line (20 exact nodes, all 105 components); ranks of $\mathfrak L_{\mu_1}$ (28/29 at section, Σ and off-section points over a curved base); flat base: consistency exactly as predicted (radial solutions of N11) | (4,0) cubic base; (1,3) λ=½; flat (4,0) |
| S5 | proof; §S3: the injectivity of $K\mapsto G_\eta(\cdot,\xi_K(h))$; type N (1,3) and the eight (2,2) null classes: properties (Bianchi, Ricci-flat, $\mathcal R^2=0$, duality), and nonzero VV Einstein defect at a fibre point for $c_V=\pm3$ | four signatures; (1,3), (2,2) |
| S6 | §S3: order-≤2 fibre series of Ric(Z1) and Ric(Z2) (N's truncated series) for $W_1$, $2W_1$, $W_2$; Schur data ($W^\pm\neq0$, not null, independent invariant pairs) | (1,3), (2,2); λ=1/5, $(c_H,c_V)=(2,3)$ |
| S7 | proof (N7's homogeneous and inhomogeneous steps, S4(b)); §S2, §S4: rank 28 at section and radial points (random, Weyl, Ricci-flat at 13/56 and 1/5; $\mu_1=-\tfrac7{18}$, $-\tfrac5{12}$) | (4,0), (1,3), (2,2) |
| S8 | proof; §S4: $d\log r/d\theta=\tfrac12$, ν′ equal to N8's value, constant term $A$ (sympy) | — |
| S9 | proof; §S4: at $p$, λ=1/5, Weyl base: $(v',n')=0$ the unique solution, $\mathrm{Ric}(p)=-g/c$ | (4,0) |
| S10–S12, S16 | proof | (H1′) |
| S13 | proof; §S5: ranks and consistency of the 42-unknown system at section and radial points (random, $\mathrm{Ric}_g=-c_Hbg$, rank-1 Ricci and Ricci-flat at λ=1/5), Q6's solution with $\xi=kv'^\flat$ for constant curvature, generic rank and inconsistency off the section | (1,3) λ=½; (4,0) λ=1/7, 1/5 |
| S14, S15 | proof; §S5: residues of $576D_1/Q_2$; the quadratic $18kD_1^2+N_2$; the derivation of (D) from (C′), (B) for an arbitrary branch (sympy) | both real forms |

Run: `cd computation/geometric-unity && python3 S_z1_off_section.py > S_z1_off_section.out 2>&1`. Result: 23 checks, 23 True, each asserted; 260 s on one core, a second run 262 s, with identical output apart from the timing annotations.

## S.3 Limits (layers b, d)

- **Open.** (C1), the twisted cone $x(W_A+aW'_C\mathrm{Ad}_\Gamma)$ with $a\neq0$, in (14,0) at $\lambda=\tfrac15$, at type-(ii) points away from Π. There the base has $\mathrm{rank}\,\mathrm{Ric}_g\le1$ (S15). The relaxed system degenerates at Σ when $\mathrm{rank}\,\mathrm{Ric}_g\le1$ (S13(a)), and Q10/Q11 need $(\star)$ near a point of $s(X)$.
- **Computer-assisted inputs.**
  - N6(f), as in N.3. S7–S10 and S12–S16 rest on it; S1–S6 do not. S5(b) uses the computed constant α (N9) but not N6(f), since $\hat{\mathcal A}\equiv0$ there.
  - S6 is reduced by Schur to three exact evaluations per signature, and S5(c), (d) to 2 + 16 exact evaluations. They are used only for F9 over indefinite bases, which lie outside the GU setting at $\lambda=\tfrac15$.
  - S14 uses Q8's certificate (exact resultants).
- **Continuation arguments.**
  - S8 and S14 continue $r$, resp. an algebraic branch $t$, along paths in the complexified fibre $F_x^{\mathbb C}$.
  - The identities continued are polynomial (S8: $P_0(r)=\rho$; S14: $k=\kappa(\tau)$, $(H/x)dt=v'^\flat$, $\mathcal G(t)-m$ constant). The paths avoid proper analytic subsets, whose complements are connected.
  - This is the step a referee should check most carefully.
- **The dichotomy needs a full-rank point on Σ.** S7 has it for every base. S13 has it except for $\mathrm{Ric}_g=-c_Hbe^\theta g$ (at most one θ) and, at $\lambda=\tfrac15$, for $\mathrm{rank}\,\mathrm{Ric}_g\le1$.
- **Hypotheses kept:** $n=4$; the real forms of J1–J2; $\kappa\neq0$; $T_e\in\mathcal V$ slaved with $\det M\neq0$ on $U$ (solutions crossing $\det M=0$ are not covered; the continued branches in S14 may pass $\det M=0$ in $F_x^{\mathbb C}$, which is harmless); $T_o$ in the invariant sector; for existence, $c_V>0$.
- **M1** is a symmetry of the reduced system only (as M.3).
- **Reconstruction:** reading R-Z of F with Z1. No layer-(b) choice is changed.

## S.4 Consequences for gaps and claims

- **GU-G14.**
  - Z1 under (H1′), with no assumption that $U$ meets $s(X)$:
    - (7,7): no member with $(x,y)\neq0$ has a solution on any open set (S16(a)), as in Z2 (P7).
    - (14,0): $\lambda=\tfrac15$, and for every member except (C1) the only solutions are Q15's flat-base solutions with constant $(t,u)$ and $T_e\equiv0$, with no analyticity assumption (S10).
  - So the limitations "open sets away from $s(X)$" and "type-(iii) points without analyticity" of N.3 and Q.3 are removed, except for (C1) in (14,0) at $\lambda=\tfrac15$ away from Π (S15).
  - The two metrics now differ only there and in the existence of Q15's flat solutions in (14,0).
- **GU-054.** Status unchanged (C). In the slaved vector sector of Z1 no curved 4D solution of $\Upsilon=0$ exists on any open subset of $Y$, apart from the open (C1) case.
- **F9.** Proven for Z1 at $n=4$, over every base signature and on every open set: Einstein iff $\lambda=\tfrac15$ and ℵ is flat (S5). With L2 (Z2) the conjecture is proven at $n=4$ for both metrics. For indefinite bases the proof is computer-assisted (S6 and 18 exact evaluations).
- **GU-G16.** The twisted cone remains the one member where the proof is incomplete, now only in (14,0) at $\lambda=\tfrac15$ away from Π.
- **N.3, N.4, Q.3, Q.4, F (F9), J.3.** The open items "open sets away from $s(X)$", "type-(iii) points" and "F9 for indefinite bases and off $s(X)$" are answered as above.
  - Q4's remark ("Without $U\cap s(X)\neq\emptyset$ this would need F9 for indefinite bases") is too pessimistic: F6 at the section of the same fibre suffices (S5(a)). It is not wrong as a statement about Q's proof.
  - No statement made before this note is contradicted.
