# Q. Υ = 0 near $s(X)$ in the Zorro metric Z1: the members round 10 leaves open

Claims: GU-054 (context). Gaps: GU-G14 (and GU-G16). Code: `computation/geometric-unity/Q_z1_remaining.py` → `.out` (sections Q0–Q9; exact: Python Fractions, also as coefficients of N's truncated power series, and sympy rationals; fixed seeds; 258 s on one core, a second run 259 s, identical apart from the timing annotations). It imports `N_z1_neighbourhood.py` and `M_remaining.py` (and through them `L_nonradial.py`, `J_slaved_sector.py`, `F_jets.py`, …), all unmodified. Output sections: §Q0 inputs (M's metric-free sections), §Q1 Lemma Q1(b), §Q2 Lemmas Q2, Q5, §Q3 Lemma Q6, §Q4 Lemma Q7, §Q5 Proposition Q8, §Q6 Lemma Q9, §Q7 Proposition Q11, §Q8 Propositions Q12, Q13, §Q9 Proposition Q15. Record: `metrics/geometric-unity/Q-z1-remaining-members.yaml`. Builds on N (N1, N2, N4, N6, N7, N9–N12), M (M1–M4, M7, M8), L (L1, L3, L5, L8), F (F4, F6, F7) and I (I6, I7(b)). Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## Q.0 Question and setting (layer b)

N.3 lists as open in Z1 "the members left open by L.3": those M treats in Z2, namely the cone $\ell=0$ other than the multiples of $W_A$, $N=0$ in (7,7), and $(x,y)=0$. This note treats them in Z1, near $s(X)$.

**Hypotheses.** (H1) of N.0: Z1, $n=4$, $\lambda\neq\tfrac14$, $c_Hc_V\neq0$; (14,0) (base (4,0), $\lambda<\tfrac14$) and (7,7) (base (1,3), (3,1) or (2,2), $\lambda>\tfrac14$); $\kappa\neq0$; $T_o=t\Phi^1+u\Phi^1\Gamma$; $T_e=\Sigma\,(dt,du)^\top$ slaved in $\mathcal V$ with $\det M\neq0$; $U$ connected and open, $\Upsilon=0$ on $U$, $U\cap s(X)\neq\emptyset$. (H1) includes $(x,y)\neq0$; Q3 drops it and $U\cap s(X)\neq\emptyset$.

**Notation.** As in N.0 ($E$, $s$, $P_H$, $V_0$, $c=4c_V(1-4\lambda)$, $b$, $\beta_S$, $\sigma_g$, $Z$, $U^*$, Σ) and M.0 ($\ell$, $N$, $\sigma$, $\delta$, $\theta$, $G_\epsilon$). Further:
- $a:=w/x$ (M's $c$, renamed because $c$ is N's constant), $A:=a^2$, $\tau:=xt/\kappa$. In (14,0) $a=w'/x'$ and τ are real (M.0).
- $B:=4bc=16(1-5\lambda)$: $B=0$ iff $\lambda=\tfrac15$, $B>-4$ in (14,0), $B<-4$ in (7,7); $B$ does not depend on $c_H,c_V$.
- $D_1:=(1+\epsilon A)+544(1-\epsilon A)\tau$, $L_1:=(11+\epsilon A)+2464(1-\epsilon A)\tau$, $Q_2:=1+768\tau+(121856-147456\epsilon A)\tau^2$,
  $N_2:=(7-3\epsilon A)+(7616-3008\epsilon A-768A^2)\tau+(2071552-3226624\epsilon A+1155072A^2)\tau^2$,
  $N_q:=(7+19\epsilon A)+(7616-2304\epsilon A-4992A^2)\tau+(2071552-3902464\epsilon A+1830912A^2)\tau^2$.
- $q_\theta:=(x_0,e^\theta g(x_0))\in\Sigma$. For a solution of $(\star)$ $\nabla^2r=\mu(r)\mathrm{Ric}+\nu g_\aleph$: $\mu_0:=\mu(r_0)$, $k:=\mu'(r_0)$, $\Phi:=\mu\mu''/\mu'^2$.
- $P_{ii}\subset U\cap s(X)$: the points of type (ii) of N4; $P_H:=\{p\in P_{ii}:H(t(p))\neq0\}$; $P_0:=\mathrm{int}(P_{ii}\setminus P_H)$ (interior in $s(X)$).

**Classes** (M2, M7, M8; up to the symmetries M1):
- (C1) *the twisted cone* $p=(x,0,0,ax)$, $x\neq0$, $a\neq0$, i.e. $x(W_A+aW'_C\mathrm{Ad}_\Gamma)$, in both real forms. In (7,7) this covers $N>0$ directly and $N<0$ through M1(b): $(0,y,z,0)\mapsto(-y,0,0,-z)$ on $(U,-g_\aleph)$, $a=z/y$. For $z=0$ (multiples of $W'_A$) the image is a multiple of $W_A$ (N12).
- (C2) the (7,7) null cone members $(x,\pm x,\pm ax,ax)$, $x\neq0$ ($\ell=N=0$).
- (C3) (7,7), $y=\pm x\neq0$, $\ell\neq0$.
- (C4) $(x,y)=0$.

The cone members with $a=0$, $\mathrm{span}(W_A,W'_A)$, are multiples of $W_A$ up to M1 and are covered by N12.

**Real forms.** As in L.0, the algebraic statements are polynomial identities, valid in both real forms. Where a sign matters the variables τ, $a$ and (Q15) $x't'$, $x'u'$ are real.

**Inputs (proven earlier).** N1, N2, N4's argument, N6 ((f) computer-assisted, REVIEW8-GU R10-1), N7's proof, N9, N10, the first two steps of N11; M1, M2, M3 (steps 2–4), M4, the reduced identities of M7 and M8; L1, L3(a), L5; F4, F6; I6, I7(b).

## Q.1 Results (layer c)

**Lemma Q1 (what transfers; the swap in Z1).** *[proven; (b) also checked]*

(a) The reduced system of J1–J2 ($S$, $M$, Σ, $H$, $\mathcal Q$, $\mathcal Q'$, ϖ, $c_i$, $x_i$) does not involve the metric. So M1, M2, M3's statements about coefficients (Σ has no $du$ column; $E^{(13)}=u(\kappa+208xt)g+\mathcal Q^{(13)}_{tt}(dt\otimes dt-|dt|^2g)$), M4 ($H$, $q$, $q'$, μ′), M7's identities and M8's θ-structure hold verbatim in Z1.

(b) $g_\aleph=c_H\,y(dx,dx)+c_V\,G_\lambda(Dy,Dy)$ is linear in $(c_H,c_V)$, so $-g_\aleph(c_H,c_V)=g_\aleph(-c_H,-c_V)$: the Z1 metric of the same ℵ and λ with constants $(-c_H,-c_V)$. It has the same $c_Hc_V\neq0$, signature (7,7) again, and the same Levi-Civita connection and Ricci tensor. $s(X)$, Σ, $E$ and $s$ are unchanged; $c$, $b$, $\beta_S$, $\sigma_g$ change sign; $B$ and N1(d) do not. So (H1) is preserved, and every statement below that is proven for all $c_H,c_V$ with $c_Hc_V\neq0$ applies to the image of M1(b).

*Proof.* (a) N3(a); output §Q0 re-runs M's sections M0–M4, M6, M7 (22 asserted checks). (b) $Dy$ does not depend on $(c_H,c_V)$. Output §Q1 checks the 2-jets off the section (curved (1,3) base, λ=½), Ric, N1(d) and the signature at a point of $s(X)$. ∎

**Lemma Q2 (Z1 has no nonconstant concircular functions).** *[proven]* Let $W\subset Y$ be open and $r$ a function on $W$ with $\nabla^2r=\psi\,g_\aleph$. Then $dr=0$.

*Proof.*
- The fibres are totally geodesic with the metric $c_VG_\lambda$ of Z2 (N1(a)), so L3(a) on each fibre piece gives $\psi=0$. Hence $Z:=\nabla r$ is parallel.
- Then $R(\cdot,\cdot)Z=0$, so $\mathrm{Ric}(Z,\cdot)=0$, and N1(d) gives $\langle E,Z\rangle=-c\,\mathrm{Ric}(E,Z)=0$.
- Differentiating: $0=\langle\nabla_YE,Z\rangle=\tfrac12\langle P_HY,Z\rangle=\tfrac12\langle Y,Z^H\rangle$ (N1(b)). So $Z$ is vertical.
- For basic $X,Y$: $0=\langle\nabla_XZ,Y\rangle=-\langle Z,(\nabla_XY)^V\rangle$. By N1(a) and the identity $G_\lambda(Z,X^\flat\odot Y^\flat+\mu\,y(X,Y)y)=Z(X,Y)$ (output §Q2), the symmetric part of the right side is $\tfrac{c_H}2Z(X,Y)$. So $Z=0$. ∎

This is N3(b) without its Ricci hypothesis. In Z2 a concircular function can have a parallel horizontal gradient (L3(b)); in Z1 the Euler field rules that out.

**Proposition Q3 ($(x,y)=0$; class (C4)).** *[proven]* Assume (H1) except $(x,y)\neq0$ and $U\cap s(X)\neq\emptyset$. Then $\Upsilon=0$ on a connected open $U$ iff $t=u=0$ on $U$, i.e. $T=0$, for every Z1 metric of the setting. Conversely $T=0$ solves $\Upsilon=0$ on every $U$.

*Proof.* As M8, with Q2 in place of L3(b).
- $(z,w)=0$: $E^{(i)}=\kappa(t,u)_ig$.
- $(z,w)\neq0$: $zE^{(1)}+\epsilon wE^{(13)}=\kappa(zt+\epsilon wu)g$ gives $zt+\epsilon wu=0$. If $z^2+\epsilon w^2=0$, then $\theta=0$. Otherwise $(t,u)=\theta(\epsilon w,-z)/(z^2+\epsilon w^2)$ and $\Xi=\kappa\theta(z^2+\epsilon w^2)^{-1}g+\Xi_0=0$ remains; it has no curvature term.
- On the open set $\{d\theta\neq0\}$, $r=F(\theta)$ with $F'=|\kappa^2-147456\epsilon\theta^2|^{-13/12}$ removes the traceless part of Ξ, so $\nabla^2r=\psi g_\aleph$, and Q2 gives $dr=0$, against $F'\neq0$. So $d\theta=0$ on $U$, θ is constant, and the trace $\kappa\theta/(z^2+\epsilon w^2)=0$ gives $\theta=0$.
- Then $T_e=0$, $E^{(1)}=\kappa t\,g$, $E^{(13)}=\kappa u\,g$, so $t=u=0$. Conversely, with $x_1=x_{13}=0$ and $t=u=0$ every term vanishes. ∎

So M8 holds in Z1 on every open set; its trace step is not needed.

**Proposition Q4 (null cone members; class (C2)).** *[proven]* Assume (H1) in (7,7) with $p=(x,sx,sax,ax)$, $s=\pm1$, $x\neq0$. Then $\Upsilon=0$ has no solution on $U$.

*Proof.* By M7(a)'s identities (Q1(a)), $E^{(1)}-sE^{(13)}=x^{-1}\sigma(\kappa+208\sigma)g$. So $\sigma\in\{0,-\kappa/208\}$ is constant on $U$, $d\sigma=0$ and $T_e=x^{-1}\Sigma_t\,d\sigma=0$. Then $E^{(1)}=-xG+c_1g=0$: Z1 is Einstein on $U$. At a point of $U\cap s(X)$ this forces $\lambda=\tfrac15$ (F6), against $\lambda>\tfrac14$. ∎

Without $U\cap s(X)\neq\emptyset$ this would need F9 for indefinite bases, which is open.

**Lemma Q5 (the γΓ-channel of (C1) in Z1; trichotomy).** *[proven]* Assume (H1), $p=(x,0,0,w)$, $x\neq0$.

(a) $U^*=\{dt\neq0\}$, and $u\equiv0$ on $\overline{U^*}$.

(b) On $\mathrm{int}\,Z$ ($Z=\{dt=du=0\}=\{dt=0\}$), $(t,u)$ is locally constant, $T_e=0$, and Z1 is Einstein with $\mathrm{Ric}=-g_\aleph/c$. So $c_1(t,u)=6x/c$ and $u(\kappa+208xt)=0$ there.

(c) Every point $p'\in U\cap s(X)$ is of type (i), (ii) or (iii) of N4. In cases (i), (iii), $\lambda=\tfrac15$ and $\mathrm{Ric}_g(\pi(p'))=0$. In case (ii), with $N$ a neighbourhood of $p'$ such that $N\cap\mathrm{int}\,Z=\emptyset$: $u\equiv0$ on $N$, and on $N$
$$(E_t)\qquad -xG+(\kappa t+312xt^2)g+H(t)\big(\nabla^2t-\Box t\,g\big)+q(t)\,dt\otimes dt+q'(t)\,|dt|^2g=0 .$$
Where $H(t)\neq0$, $r=F(t)$ with $F''/F'=q/H$ satisfies $(\star)$ with $\mu=xF'/H\neq0$ and $\mu'=-N_2/(18D_1^2)$ (M4).

*Proof.*
- On $W:=\{dt\neq0\}$, M3 steps 2–3 (metric-free) give $\mathcal Q^{(13)}_{tt}=0$ and $u(\kappa+208xt)=0$, so $u=0$ on $W$ (else $t=-\kappa/(208x)$ on an open set with $dt\neq0$), and on $\overline W$.
- On $V:=\mathrm{int}\{dt=0\}$, $T_e=\Sigma_tdt=0$ (Σ has no $du$ column). So $E^{(13)}=u(\kappa+208xt)g$ and $E^{(1)}=-xG+c_1g$, and Z1 is Einstein on $V$. By N1(d) its constant is $-1/c$, so $G=\tfrac6cg_\aleph$ and $c_1=6x/c$ (output §Q2). With $dt=0$, the constancy of $c_1=\kappa t+x(312t^2+104\epsilon u^2)$ makes $u^2$, hence $u$, locally constant: $du=0$ on $V$.
- If $dt(p')=0\neq du(p')$, then $p'\notin\overline W$ ($\overline W\subset\{u=0\}$, a hypersurface near $p'$), so $p'\in V$, against $du=0$ on $V$. Hence $U^*=W$ and $\mathrm{int}\,Z=V$: (a), (b).
- (c): types (i), (iii) as in N4 (Einstein at $p'$ with constant $-1/c$, then F6). In case (ii), $U^*\cap N$ is dense in $N$, so $u=0$ on $N$; M3's $E^{(1)}$ on $u=0$ is $(E_t)$, and the reparametrisation is M4's. ∎

Compared with M3, step 1 (L2(a)) is replaced by (b): in Z1 Einstein pieces exist (N2, F7).

**Lemma Q6 (degree 1 at a point, any μ).** *[proven; uses N6(f), computer-assisted]* Let $(\star)$ hold on a neighbourhood of $p\in s(X)$, with μ smooth near $r_0=r(p)$, $\mu_0=\mu(r_0)\neq0$ and $k=\mu'(r_0)$. Put $v=\nabla r(p)$, $n=\nabla\nu(p)$. Then:

(a) $v=v_EE$ and $n=n_EE$, with $v_E=\mu_0\big(2b-\tfrac{\sigma_g}2\big)$ and $n_E=-k\beta_Sv_E$;

(b) $g$ is Einstein at $x_0$, and $d^\nabla\mathrm{Ric}_g(x_0)=0$;

(c) $k\,c\,\sigma_g^2=4b\,(4kbc-1)$.

The same holds at the points of Σ, with $g$ replaced by $e^\theta g$ and $\sigma_g$ by $e^{-\theta}\sigma_g$.

*Proof.* The Ricci identity for $\nabla_X\nabla r=\mu(r)\widehat{\mathrm{Ric}}X+\nu X$ is N7's $(\ast)$ with $\mu_1$ replaced by $k$ in the gradient term and $\mu_1r_0$ by $\mu_0$ on the right. N7's steps 1–6 treat the two coefficients separately, so they go through with these replacements:
1. $v_H=n_H=0$ and $n_0=\tfrac kcv_0$;
2. $R^S(U,W)v_0+kb\big[\langle v_0,U\rangle W-\langle v_0,W\rangle U\big]=0$, so $v_0=0$ (commuting $U,W$ if $kb\neq0$; the commutant of L3 if $kb=0$);
3. $\mathrm{Ric}_g=\tfrac{\mathrm{scal}_g}4g$ and $v_E=\mu_0(2b-\sigma_g/2)$;
4. $n_E=-k\beta_Sv_E$;
5. $v_E\big(-\tfrac14+kbc+\tfrac{kc\sigma_g}4\big)=\tfrac{\mu_0\sigma_g}8$, which with step 3 is (c);
6. $d^\nabla\mathrm{Ric}_g=0$. ∎

For $k=\mu_1=-\tfrac7{18}$, (c) is N7(c). Output §Q3 solves the full system (N's 1274×28 `star_system` with $\mu_1\to k$, right side divided by $\mu_0$):
- (4,0), $\lambda=\tfrac17$, constant curvature $\sigma_g=1$: solvable exactly at the predicted $k=-\tfrac{14}{17}$, with the predicted solution; not at $k=0,-\tfrac7{18}$;
- the same at $\sigma_g=4b$, where (c) has no solution: not solvable for three $k$;
- flat base: only at $k=1/B=\tfrac7{32}$; at $\lambda=\tfrac15$ for every sampled $k$, and with $\sigma_g=1$ only at $k=0$;
- (1,3), $\lambda=\tfrac12$, $\sigma_g=\tfrac12$: exactly at $k=-\tfrac2{45}$.

**Lemma Q7 (persistence along Σ).** *[proven; uses N6(f)]* Let Θ be an open interval with $0\in\overline\Theta$ such that $(\star)$ holds near $q_\theta$ for every $\theta\in\Theta$, with $\mu(r)\neq0$ there, and put $r(\theta):=r(q_\theta)$. Then $g$ is Einstein at $x_0$ and, for $\theta\in\Theta$,
$$(B)\quad \frac{dr}{d\theta}=c\,\mu(r)\big(2b-\tfrac12e^{-\theta}\sigma_g\big),\qquad (C)\quad c\,\mu'(r)\big(e^{-2\theta}\sigma_g^2-16b^2\big)=-4b .$$
Consequently:

(a) $b=0$ ($\lambda=\tfrac15$): if μ′ vanishes on no interval, then $\sigma_g=0$, i.e. $\mathrm{Ric}_g(x_0)=0$, and $\nabla r=0$ at every $q_\theta$.

(b) $b\neq0$, $\sigma_g=0$: $r$ is strictly monotone on Θ, and $\mu'=1/B$ on the interval $r(\Theta)$.

(c) $b\neq0$, $\sigma_g\neq0$: with $\rho:=e^{-\theta}\sigma_g/(4b)$, $\rho\notin\{0,\pm1\}$, $r$ is strictly monotone on Θ, and on Θ
$$(C')\quad B\mu'(1-\rho^2)=1,\qquad (D)\quad \Phi\,(1-\rho)+4\rho^2=0 .$$

(d) For linear $\mu=\mu_1r$: $\sigma_g=0$ and $\lambda\in\{\tfrac15\}\cup\{4\mu_1bc=1\}$. For $\mu_1=-\tfrac7{18}$ this is N7(d) ($\tfrac{13}{56}$).

*Proof.*
- Q6 at $q_\theta$, which lies on the section of $e^\theta g$ (N.0): $\mathrm{Ric}_{e^\theta g}=\mathrm{Ric}_g$, $\sigma\mapsto e^{-\theta}\sigma_g$, $b,c$ unchanged. $\theta\mapsto q_\theta$ is an integral curve of $E$, so $dr/d\theta=Er=c\,v_E$: (B). Q6(c) at $q_\theta$ is (C).
- (a): (C) reads $\mu'(r)\,\sigma_g^2e^{-2\theta}=0$. If $\sigma_g\neq0$, then $\mu'(r(\theta))=0$ on Θ while $dr/d\theta=-\tfrac12c\mu e^{-\theta}\sigma_g\neq0$, so μ′ vanishes on the interval $r(\Theta)$. Hence $\sigma_g=0$, and $v_E=0$.
- (b): (C) gives $\mu'=4b/(16b^2c)=1/B$, and $dr/d\theta=2bc\mu\neq0$.
- (c): (C) is (C′), so $\rho^2\neq1$, and (B) is $dr/d\theta=\tfrac12B\mu(1-\rho)\neq0$. Differentiating (C′) along Θ, with $\rho'=-\rho$, (B) and (C′), gives $\tfrac12B\mu'\big[\Phi(1-\rho)+4\rho^2\big]=0$ (output §Q4).
- (d): $\Phi=0$, so (D) gives $\rho=0$; (b) then gives $4\mu_1bc=1$ unless $b=0$. ∎

Output §Q4 also gives the linear values $\lambda=\tfrac{13}{56},\tfrac{23}{100},\tfrac{337}{240}$ for $\mu_1=-\tfrac7{18},-\tfrac5{12},-\tfrac3{289}$. So for a nonlinear μ the quantisation of λ in N7(d) becomes the ODE condition (B), (C) on μ along $r(\Sigma)$.

**Proposition Q8 ((C1): no persistence with $b\neq0$).** *[computed: exact certificate; the rest proven]* For (C1) on $u=0$, as functions of $(\tau,A)$:
$$\mu'=-\frac{N_2}{18D_1^2}\ (\text{M4}),\qquad \mu\mu''=\frac{A(6A-11\epsilon)\,L_1Q_2}{81D_1^4},\qquad \Phi=\frac{4A(6A-11\epsilon)\,L_1Q_2}{N_2^2}.$$

(a) μ′ is constant on an interval iff $A(6A-11\epsilon)=0$: $a=0$ ($\mu'=-\tfrac7{18}$) or, in (14,0) only, $a^2=\tfrac{11}6$ ($\mu'=-\tfrac3{289}$), as in M4.

(b) For real $a\neq0$ and admissible $\lambda\neq\tfrac15$, neither Q7(b) nor Q7(c) holds on any interval.

*Proof.*
- $\mu\mu''=\mu\,(d\mu'/dt)/F'=x\,(d\mu'/dt)/H$; with M's $H$ and μ′ this is the formula (output §Q5, from M's reduced tables). (a) follows, since $L_1,Q_2\not\equiv0$.
- Q7(b) needs $\mu'\equiv1/B$, i.e. by (a) $a^2=\tfrac{11}6$ in (14,0) and $B=-\tfrac{289}3$, i.e. $\lambda=\tfrac{337}{240}>\tfrac14$: not on the (14,0) side.
- Q7(c): $\Phi\not\equiv0$, since otherwise (D) gives $\rho=0$. On a subinterval where $\Phi\neq0$, (D) gives $\rho=1+4X/\Phi$ with $X:=\rho^2=1-1/(B\mu')$. So $\tilde R:=(\Phi+4X)^2-X\Phi^2$ vanishes on an interval of τ (τ is strictly monotone along Θ, since $r$ is and $F'\neq0$), hence identically in τ.
- Certificate (output §Q5): $\tilde R$ is rational in $(\tau,A,B)$, with numerator of degree 12 in τ. If $\tilde R\equiv0$ in τ at $(A_0,B_0)$, then $B_0$ is a common root of the specialisations $\tilde R(\tau_0)$, $\tau_0\in\{0,1,-1,2\}$, so $A_0$ is a root of $\mathrm{Res}_B(\tilde R(0),\tilde R(1))$ and of $\mathrm{Res}_B(\tilde R(0),\tilde R(-1))$. Their gcd is a constant times $A^4(3A-7\epsilon)^2(6A-11\epsilon)^4$. At $A=\tfrac73\epsilon$ (where $N_2(0)=0$), $\tilde R(0)$ vanishes identically in $B$ but $\tilde R(1),\tilde R(-1),\tilde R(2)$ have no common root. At $A=0$ and $A=\tfrac{11}6\epsilon$ the only common root is $B=-\tfrac{18}7$, resp. $-\tfrac{289}3$, and there $\tilde R\equiv0$. At both points $X\equiv0$, i.e. $\rho\equiv0$, against $\sigma_g\neq0$. ∎

At $(0,-\tfrac{18}7)$, i.e. $W_A$ at $\lambda=\tfrac{13}{56}$, the certificate recovers N7(d)'s second value.

**Lemma Q9 (points with $H(t)=0$).** *[proven; resultants computed]* Assume (H1) and (C1), and let $p'\in N$ (a neighbourhood as in Q5(c)) lie on $s(X)$ or on Σ, with $H(t(p'))=0$. Write $\vartheta:=dt|_V$, $\vartheta_0$ for its $V_0$-part, $dt_H:=dt|_H$, and $g'$ for the base metric ($g$ or $e^\theta g$) whose section passes through $p'$.

(a) $\tfrac x2s_Y+c_1+q'|dt|^2=x\beta_S$, $q\,\vartheta_0\otimes\vartheta_0=0$, $q\,(Et)^2=xbc$, $q\,(Et)\vartheta_0=0$, $q\,dt_H\otimes\vartheta=0$, and $x\,\mathrm{Ric}_{g'}+xc_Hb\,g'=q\,dt_H\otimes dt_H$.

(b) $H$ and $q$ have no common zero where $\det M\neq0$, for real $a\neq0$.

So: if $b\neq0$, then $Et(p')\neq0$. If $b=0$ and $\mathrm{Ric}_{g'}(\pi(p'))=0$, then $dt(p')=0$. If $p'\in s(X)$ and $t$ is constant along $s(X)$ near $p'$, then $\mathrm{Ric}_g(\pi(p'))=-c_Hb\,g$.

*Proof.*
- (a) At $p'$ the Hessian term of $(E_t)$ drops: $-x\,\mathrm{Ric}+\alpha g_\aleph+q\,dt\otimes dt=0$ with $\alpha:=\tfrac x2s_Y+c_1+q'|dt|^2$. Use N6(a) for $\mathrm{Ric}(p')$: $\mathrm{Ric}_{HV}=0$, $\mathrm{Ric}(E,E)=-1$, $\mathrm{Ric}|_{V_0}=\beta_Sg_\aleph$, $\mathrm{Ric}_{HH}=\mathrm{Ric}_{g'}-\tfrac{c_H}cg'$; and $1+c\beta_S=-bc$. The $V_0$ block (dimension 9) forces $\alpha=x\beta_S$ and $q\,\vartheta_0\otimes\vartheta_0=0$. The $(E,E)$ component is $x+\alpha c+q(Et)^2=0$. The other blocks are read off.
- (b) On $u=0$: $H=576x^2D_1/(\kappa Q_2)$, $q=-36864x^3N_q/(\kappa^2Q_2^2)$, $\det M=\kappa^4Q_2$ (output §Q6). $\mathrm{Res}_\tau(D_1,N_q)\propto A(A-\epsilon)(6A-11\epsilon)^2$ and $\mathrm{Res}_\tau(D_1,Q_2)\propto A(6A-11\epsilon)^2$. For real $a\neq0$ they vanish only in (14,0), at $A=1$, where $D_1=2$ has no zero, and at $A=\tfrac{11}6$, where the zero of $D_1$ is a zero of $\det M$. ∎

**Proposition Q10 ($\lambda=\tfrac15$, $H\neq0$: N9 for any μ).** *[proven; uses N6(f); α computed (N9)]* Under (H1) in (14,0) with $\lambda=\tfrac15$, let $(\star)$ hold near $p\in s(X)$ with $\mu(r(p))\neq0$, $\mathrm{Ric}_g=0$ near $x_0$ and $\nabla r=0$ on Σ near $p$. Then $R_g(x_0)=0$.

*Proof.* N9 with $\mu_1r_0$ replaced by $\mu_0$:
1. $\nabla_EE=0$ at $p$ and $Er\equiv0$ on the $E$-orbit, so $\nabla^2r(E,E)=0$, i.e. $c\,\nu(p)=\mu_0$. With $\mathrm{Ric}(p)=-g_\aleph/c$ (F6), $\nabla^2r(p)=0$.
2. On the slice $S$: $\mathrm{Hess}^Sw=\psi g_S+\mu(w)\hat{\mathcal A}$, $\psi=\mu(w)\kappa_E+\nu$, with $\psi(p)=0$, $\hat{\mathcal A}(p)=0$, $\nabla\hat{\mathcal A}(p)=0$ (N6(f)) and $\nabla w(p)=0$.
3. The Ricci identity of $S$ is $R^S(Y,Z)\nabla w=(Y\psi)Z-(Z\psi)Y+\mu'(w)\big[(Yw)\hat{\mathcal A}Z-(Zw)\hat{\mathcal A}Y\big]+\mu(w)\big[(\nabla_Y\hat{\mathcal A})Z-(\nabla_Z\hat{\mathcal A})Y\big]$. Differentiated once at $p$, every term with μ″ or μ′ carries a factor $Xw(p)$, $Yw(p)$, $Zw(p)$, $\hat{\mathcal A}(p)$ or $\nabla^2w(p)$, all zero. What is left is N9's identity with $\mu_0$ for $\mu_1r_0$.
4. N9 steps 4–6: $\langle q,T'\rangle=\alpha|W|^2=0$ with $\alpha=\tfrac{45}{8c_H^2c_V}$, so $W=0$. ∎

**Proposition Q11 ($\lambda=\tfrac15$ at the zeros of $H$).** *[proven; uses N6(f); α″ computed exactly]* Under (H1) in (14,0), class (C1), $\lambda=\tfrac15$, let $P'\subset U\cap s(X)$ be open and contained in a neighbourhood $N$ as in Q5(c) (so of type (ii)), with $H(t)=0$ on $P'$. Then $R_g=0$ on $\pi(P')$.

*Proof.*
1. $t$ is constant on each component of $P'$ (a zero $t_H$ of $H$). As $s(X)$ is horizontal, $dt_H=0$ on $P'$, and Q9 gives $\mathrm{Ric}_g=0$ on $\pi(P')$.
2. Let $\Sigma':=\{(x,e^\theta g(x)):x\in\pi(P'),\ |\theta|\text{ small}\}\subset N$: a 5-dimensional submanifold through $P'$ with $T_p\Sigma'=H_p\oplus\mathbb RE$ at $p\in P'$. At its points with $H(t)\neq0$, Q6(a) gives $\nabla r=0$ ($v_E=\mu_0(2b-\sigma/2)=0$), so $dt=0$. At its points with $H(t)=0$, Q9 gives $dt=0$ ($q(t)\neq0$ by Q9(b)). So $dt\equiv0$ on Σ′.
3. Fix $p\in P'$ and put $K:=\nabla^2t(p)=\partial^2t(p)$ (as $dt(p)=0$). By step 2, $K(A,\cdot)=0$ for $A\in H_p\oplus\mathbb RE$.
4. *Second-order expansion of $(E_t)$ at $p$*, in any coordinates. Put $\tilde T:=\mathrm{Ric}+g_\aleph/c$. At $p$: $\tilde T=0$ (F6); $ds_Y|_V=0$ (N6(c), (f)); $G-\tfrac6cg_\aleph=\tilde T-\tfrac12(\mathrm{tr}\tilde T)g_\aleph$; $H(t_H)=0$, $dt(p)=0$, and $c_1(t_H)=6x/c$ (order 0). For $A,B\in V_0$ and all $i,j$:
$$-x\big[\partial_A\partial_B\tilde T_{ij}-\tfrac12\partial_A\partial_B(\mathrm{tr}\tilde T)\,g_{ij}\big]+\beta_{AB}\,g_{ij}+H'(t_H)K_{AB}K_{ij}+q(t_H)\big(K_{Ai}K_{Bj}+K_{Bi}K_{Aj}\big)=0,$$
   with $\beta_{AB}:=c_1'K_{AB}-H'K_{AB}\,\mathrm{tr}K+2q'(K^2)_{AB}$. For $i,j\in H$, and for $i=j=E$, all $K_{\cdot i}$ terms vanish. Eliminating $\beta_{AB}$ between the trace of the $HH$ block and the $(E,E)$ component gives the $K$-free condition
$$(O2)\qquad \frac1{4c_H}\,\mathrm{tr}_g\,\partial_A\partial_B\tilde T\big|_{HH}=\frac1c\,\partial_A\partial_B\tilde T(E,E)\qquad(A,B\in V_0).$$
5. Contract (O2) over a basis of $V_0$ and its dual. The result is an $O(4)$-invariant function of $W=R_g(x_0)$, and it is quadratic: along the fibre over $x_0$, $g_\aleph$ is $R$-free and Ric has degree ≤2 in $R_g(x_0)$ and is at most linear in $\nabla R_g(x_0)$; the isometry induced by $x\mapsto-x$ (N9 step 5) fixes that fibre, acts by $-1$ on $H$ and $+1$ on $V$, so the $HH$ and $(E,E)$ components and their fibre derivatives do not depend on $\nabla R_g(x_0)$. There is no invariant linear form on Weyl tensors, so it equals $\alpha''|W|^2$. Output §Q7 gives $\alpha''=-\dfrac{3c_V}{2c_H^2}\neq0$. For the definite base $W=0$. ∎

Output §Q7: $\tilde T$ vanishes at $p$ and its $HH$, $VV$ blocks to first order along the fibre; over the flat base it vanishes to second order (N2); α″ is exact for $(c_H,c_V)=(2,3),(1,1),(3,\tfrac12)$, the same for two Weyl tensors at $(2,3)$, four times as large for $2W$, and unchanged with $\nabla R(x_0)\neq0$. The traceless $HH$ part of the second-order condition, $(\partial_A\partial_B\tilde T|_{HH})_0=0$, is a further necessary condition; its only invariant pairing, with $(U\eta^{-1}W+W\eta^{-1}U)_0$, gives 0 (§Q7), so it is not used.

**Proposition Q12 (flat base, any μ).** *[proven]* If ℵ is flat on $\pi(W)$, $W\subset Y$ open, and $(\star)$ holds on $W$ with $\mu(r)\neq0$, then $dr=0$, unless $b\neq0$ and $\mu'=1/B$ on an interval. In particular $dr=0$ at $\lambda=\tfrac15$.

*Proof.* N11's first two steps do not use μ and give $r=f(s)$. Then $\nabla^2r=f''ds^2+\tfrac{f'}{2c}g_\aleph|_H$, and N2's $\mathrm{Ric}=\beta_Sg_\aleph+b(g_\aleph|_H+c\,ds^2)$ give, on $V_0$, $H$ and $E$: $\nu=-\mu\beta_S$, $f'=2bc\,\mu(f)$, $f''=bc\,\mu(f)$. With $f''=2bc\,\mu'(f)f'$ this is $bc\,\mu(f)\,(4bc\,\mu'(f)-1)=0$ (output §Q8). ∎

**Proposition Q13 ((C3): $x=\pm y$, $\ell\neq0$).** *[proven; uses N6(f)]* Under (H1) in (7,7) with $y=\pm x\neq0$ and $\ell\neq0$, $\Upsilon=0$ has no solution on $U$.

*Proof.*
1. On $\mathrm{int}\,Z$, $T_e=0$ and Z1 is Einstein ($x\neq0$); at $s(X)$ this needs $\lambda=\tfrac15$ (F6). So every point of $U\cap s(X)$ is of type (ii) (N4), with a neighbourhood $N$ on which $U^*$ is dense.
2. L1 ($\ell\neq0$) gives $dt\wedge du=0$, so $(t,u)=\gamma(\rho)$ locally on $U^*$ (L5(d)). Where $\hat h\neq0$, $r=F(\rho)$ with $F''/F'=\hat q/\hat h$ is concircular (L8 step 3), so $dr=0$ by Q2, against $d\rho\neq0$. Hence $\hat h=0$ on $U^*\cap N$.
3. M7(b) steps 2–5 (metric-free; Q1(a)): $\sigma=-\kappa/208$ on $U^*\cap N$, hence on $N$, and $r\propto|\zeta|^{-5/8}$ ($\zeta\neq0$, as $\det M=576\kappa^2\zeta^2/(169x^2)$) satisfies $(\star)$ with $\mu=-\tfrac5{12}r$ on $N$.
4. Q7(d) at $p\in U\cap s(X)$: $\lambda\in\{\tfrac15,\tfrac{23}{100}\}$, both $<\tfrac14$ (output §Q8). ∎

**Theorem Q14 (Z1: the remaining members near $s(X)$).** Assume (H1).

(a) *(C1) in (7,7)*, including the $N<0$ family through M1(b): **no solution** on $U$. *[proven; uses N6(f) and the certificate of Q8]*

(b) *(C1) in (14,0)*: *[proven; uses N6(f), the certificate of Q8, α (N9) and α″ (Q11)]*
- $\lambda=\tfrac15$, and $g$ is Ricci-flat on $\pi(U\cap s(X))$.
- No point of $U\cap s(X)$ is of type (ii): every point lies in $\overline{\mathrm{int}\,Z}$.
- If $(t,u)$ is constant near $p\in U\cap s(X)$, in particular whenever $t,u$ are real-analytic, then near $p$: $T_e\equiv0$, $c_1=6x/c$, $u(\kappa+208xt)=0$, and $g$ is flat near $\pi(p)$.

(c) *(C2) and (C3)*: **no solution** on $U$. *[(C2) proven; (C3) proven, uses N6(f)]*

(d) *(C4)*: exactly $T=0$, on every open set, for every Z1 metric (Q3). *[proven]*

*Proof.*
- (a) Types (i), (iii) need $\lambda=\tfrac15$ (Q5(c)). Let $p$ be of type (ii); as $\lambda>\tfrac14$, $b\neq0$. There is an interval Θ, $0\in\overline\Theta$, with $H(t(q_\theta))\neq0$ on Θ: if $H(t(p))\neq0$, a neighbourhood of 0; if $H(t(p))=0$, Q9 gives $Et(p)\neq0$, so $t(q_\theta)\neq t(p)$ for small $\theta\neq0$, and since the zeros of $H$ are isolated, $\Theta=(0,\delta)$ works. On Θ, $(\star)$ holds near $q_\theta$ (Q5(c)), so Q7 applies. $\sigma_g=0$ leads to Q7(b) and $\sigma_g\neq0$ to Q7(c), both excluded by Q8(b).
- (b) If $\lambda\neq\tfrac15$, the argument of (a) excludes type (ii), and types (i), (iii) need $\lambda=\tfrac15$. So $\lambda=\tfrac15$. Suppose $P_{ii}\neq\emptyset$; take $p\in P_{ii}$ with a neighbourhood $N$ as in N4, and $P_N:=N\cap s(X)\subset P_{ii}$.
  - On $P_H\cap P_N$: Q7(a) (μ′ vanishes on no interval, M4) gives $\mathrm{Ric}_g=0$ and $\nabla r=0$ on Σ near each point, and Q10 gives $R_g=0$.
  - On $P_0\cap P_N$: Q11 (with $P_0\cap P_N$ for $P'$) gives $R_g=0$.
  - $P_H\cup P_0$ is dense in $P_{ii}$ ($P_{ii}\setminus P_H$ is closed in $P_{ii}$), so $R_g=0$ on $\pi(P_N)$.
  - On $N':=N\cap\pi^{-1}(\pi(P_N))$, a neighbourhood of $p$ over a flat base, $(\star)$ holds with $dr\neq0$ on the nonempty open set $N'\cap U^*\cap\{H(t)\neq0\}$, against Q12 ($b=0$).

  So $P_{ii}=\emptyset$, and Ricci-flatness follows from types (i), (iii). If $(t,u)$ is constant near $p$, then $p\in\mathrm{int}\,Z$: Q5(b) holds there, Z1 is Einstein near $p$, and $g$ is flat near $\pi(p)$ (N10). Real-analytic $(t,u)$: either $Z=U$, or $\mathrm{int}\,Z=\emptyset$ and every point of $U\cap s(X)$ would be of type (ii).
- (c) Q4, Q13. (d) Q3. ∎

**Proposition Q15 (flat-base constant solutions for every (14,0) member).** *[proven; the $\det M$ statement computed]* Let $g$ be flat (base (4,0)), $\lambda=\tfrac15$, $c_V>0$, and let the member have $(x',y')\neq0$. Up to M1(a), $y=0$ and $x'^2$ is replaced by $x'^2+y'^2$; let $a:=w'/x'$ for the rotated member. With $T_e\equiv0$ and constant $(t,u)$, $\Upsilon=0$ holds on all of $Y$ iff
$$r\,(\kappa+208T)=0,\qquad 312T^2+\kappa T+104R^2+\frac{15x'^2}{2c_V}=0\qquad(T=x's,\ R=x'r,\ t=is,\ u=ir).$$
Real solutions exist iff $c_V\kappa^2\ge9360(x'^2+y'^2)$, and all have $u=0$. $\det M\neq0$ at one of the two roots for every $a$, and at both except for at most two values of $c_V\kappa^2/(x'^2+y'^2)$, which depend on $a^2$ (at $a=0$ the single value $277440/29$ of N12(c)).

*Proof.* N2 at $\lambda=\tfrac15$: Z1 is Einstein on $Y$ with $G=\tfrac{15}{2c_V}g_\aleph$. By I6, with $T_e=0$ and constant $(t,u)$, $\Upsilon=0$ iff $c_1=\tfrac{15}{2c_V}x$ and $c_{13}=\tfrac{15}{2c_V}y$; $(z,w)$ do not occur. M1(a) rotates $(x',y')$ to $(\sqrt{x'^2+y'^2},0)$ and leaves $\det M$ invariant. Output §Q9: the real form of the two equations; the discriminant $\kappa^2-9360x'^2/c_V$ on $r=0$; on $\kappa+208T=0$ the left side is $\kappa^2/416+104R^2+15x'^2/(2c_V)>0$; $\det M$ on $u=0$ is $\kappa^2(\kappa^2+768\kappa T+(121856-147456a^2)T^2)$, independent of $z$; its resultant with the quadratic is a polynomial of degree 2 in $\omega=9360x'^2/(c_V\kappa^2)$ with negative constant term, and both roots on $\det M=0$ would need $a^2=-\tfrac{115}{144}$. ∎

For (C1) these are the solutions of Q14(b) at type-(i) points. For class (A) of L8 they are all the solutions (N12(b)). This answers N.4's "for the other members of (A), (B) existence was not examined" (REVIEW8-GU R10-3).

**Corollary Q16 (the family near $s(X)$ in Z1).** Assume (H1). Together with N12:

(a) **(7,7):** for every member with $(x,y)\neq0$, $\Upsilon=0$ has no solution on any connected open $U$ meeting $s(X)$.

(b) **(14,0):** for every member with $(x,y)\neq0$, a solution forces $\lambda=\tfrac15$, a Ricci-flat base on $\pi(U\cap s(X))$, and no type-(ii) points. For real-analytic $t,u$ (and for class (A) always) it is a flat-base solution with constant $(t,u)$ and $T_e\equiv0$, as in Q15; for $c_V>0$ these exist iff $c_V\kappa^2\ge9360(x'^2+y'^2)$.

(c) **$(x,y)=0$:** only $T=0$, for every metric.

*[proven; uses N6(f), the certificate of Q8 and α, α″]*

*Proof.* The classification of M10: L8(A), L8(B) and the multiples of $W_A$ up to M1 (N12); (C1)–(C3) (Q14). In (7,7) the multiples of $W'_A$ (the $N<0$ cone with $z=0$) are multiples of $W_A$ on $(U,-g_\aleph)$ by M1(b), so N12(a) applies through Q1(b). ∎

**Answers.**
1. *Transfer.* The reduced system and M1 are metric-free. The swap metric $-g_\aleph$ is the Z1 metric with $(-c_H,-c_V)$ (Q1). Z1 has no nonconstant concircular function at all (Q2); this replaces L3/L4 and makes M8's trace step unnecessary.
2. *The twisted cone $x(W_A+aW'_C\mathrm{Ad}_\Gamma)$, $a\neq0$.*
   - The γΓ-channel still kills $u$, but now Einstein pieces exist (Q5).
   - For a nonlinear μ, the degree-1 conditions at a point are those of N7 with $\mu'(r_0)$ in place of $\mu_1$ (Q6). Along Σ, N7's quantisation of λ becomes the ODE system (B), (C) for μ along $r(\Sigma)$ (Q7). For the twisted cone it has no solution for any admissible $\lambda\neq\tfrac15$ (Q8), also through the zeros of $H$ (Q9).
   - At $\lambda=\tfrac15$, the Weyl tensor is removed by N9's pairing where $H\neq0$ (Q10) and by a new pairing α″ at the zeros of $H$ (Q11). Over a flat base nothing is left (Q12).
   - Result: the same as for $W_A$ (Q14(a), (b)).
3. *The null members and $x=\pm y$ in (7,7).* No solution near $s(X)$: the first by F6 alone (Q4), the second by N7 with $\mu_1=-\tfrac5{12}$, which allows only $\lambda\in\{\tfrac15,\tfrac{23}{100}\}$ (Q13).
4. *$(x,y)=0$.* $T=0$, for every metric (Q3).
5. *Existence.* The flat-base constant solutions exist for every (14,0) member iff $c_V\kappa^2\ge9360(x'^2+y'^2)$ (Q15).
6. With N12, the Z1 item of J.3/L.3 is settled near $s(X)$ for every member, up to type-(iii) points without analyticity (Q16).

## Q.2 What is verified, precisely

| Statement | How (section of `Q_z1_remaining.out`) | Scope |
|---|---|---|
| Q1(a), inputs | §Q0: M's sections M0–M4, M6, M7 rerun silently (22 asserted checks: L's J1/J2 inputs, M1(a), M1(b), M2, M3, M4, M7, M8) | both real forms |
| Q1(b) | §Q1: 2-jets of Z1$(\pm2,\pm3)$ at a point off the section, Ric, N1(d), the constants, signatures at a point of $s(X)$ | curved (1,3) base in general coordinates, λ=½ |
| Q2 | proof; §Q2: the $G_\lambda$ pairing identity | (4,0), (1,3), (2,2), symbolic λ |
| Q3, Q4 | proof from Q2, F6 and the identities of §Q0 | — |
| Q5 | proof; §Q2: $G=\tfrac6cg$, $c_1=6x/c$ | — |
| Q6 | proof; §Q3: exact ranks and solution vectors of the 1274×28 system with $\mu_1\to k$ | (4,0) λ=1/7, 1/5; (1,3) λ=½; constant curvature, flat, $\sigma_g=4b$ |
| Q7 | proof; §Q4: (B), (C) ⇔ (C′); the derivative identity giving (D); N7(c) at $k=\mu_1$; the linear values of λ | symbolic |
| Q8 | §Q5: μ′ and μμ″ from M's tables; $\tilde R$, the two resultants, their gcd, the candidates | both real forms |
| Q9 | proof; §Q6: $H$, $q$, $q'$, $\det M$ in closed form; the resultants | both |
| Q10 | proof (N9; α from N §N7) | (14,0) |
| Q11 | proof; §Q7: $\tilde T$ along the fibre to order 2 (N's truncated series), the pairing α″, homogeneity, ∇R-independence, the traceless pairing | (4,0); $(c_H,c_V)=(2,3),(1,1),(3,\tfrac12)$ |
| Q12, Q13 | proof; §Q8 (sympy) | — |
| Q14, Q16 | proof | (H1) |
| Q15 | proof; §Q9 (sympy) | (14,0) |

Run: `cd computation/geometric-unity && python3 Q_z1_remaining.py > Q_z1_remaining.out 2>&1`. Result: 31 checks, 31 True, each asserted; 258 s on one core, a second run 259 s, with identical output apart from the timing annotations.

## Q.3 Limits (layers b, d)

- **Computer-assisted inputs.**
  - N6(f), as in N.3. Q6, Q7, Q10, Q11, Q13, Q14(a)–(c) and Q16 rest on it; Q1–Q5, Q12 and Q15 do not.
  - The certificate of Q8 (exact resultants and gcd), the resultants of Q9 and α″ of Q11 (exact evaluations; Schur and the parity step do the rest) are computations. Q15's $\det M$ statement is a computed resultant.
- **Type-(iii) points.** As in N.3: in (14,0), without analyticity, only $\mathrm{Ric}_g=0$ is proven at limit points of Einstein pieces, not flatness.
- **Away from $s(X)$.** Nothing is proven for (C1)–(C3) on open sets that do not meet $s(X)$. For (C2) it would follow from F9 for indefinite bases, which is open. (C4) is settled on every open set.
- **Definite base in Q11.** α″ is used through $|W|^2$, so only in (14,0); in (7,7), $\lambda=\tfrac15$ does not occur.
- **M1** is a symmetry of the reduced system of J1–J2, not claimed for the full Υ (as M.3).
- **Hypotheses kept:** $n=4$; the real forms of J1–J2; $\kappa\neq0$; $T_e\in\mathcal V$ slaved with $\det M\neq0$; $T_o$ in the invariant sector; for existence, $c_V>0$.
- **Reconstruction:** reading R-Z of F with Z1. No layer-(b) choice is changed.

## Q.4 Consequences for gaps and claims

- **GU-G14.**
  - Under (H1), for every member of the family: in (7,7) no solution on a neighbourhood of a point of $s(X)$; in (14,0) a solution forces $\lambda=\tfrac15$ and a Ricci-flat base, and for analytic fields it is a flat-base solution with constant $(t,u)$ and $T_e\equiv0$; for $c_V>0$ these exist iff $c_V\kappa^2\ge9360(x'^2+y'^2)$. For $(x,y)=0$ only $T=0$ (Q16).
  - So the Z1 question near $s(X)$ in this sector is settled for every member, as round 9 settled the Z2 question (M10(a)). The two metrics differ only in (14,0) at $\lambda=\tfrac15$: Z2 has no solution, Z1 has the flat-base solutions with $\Lambda_4=0$.
  - Open in Z1: type-(iii) points without analyticity; open sets away from $s(X)$; other sectors.
- **GU-G16.** In Z1, near $s(X)$, the cone $\ell=0$ and $N=0$ are no longer exceptional. The nonlinear μ of the twisted cone (M4) is excluded by its persistence along Σ (Q8); in Z2 it is still the one exception away from $s(X)$ (M10(b)).
- **GU-054.** Status unchanged (C).
- **N.3, N.4, J.3, L.3.** The "members" item is settled near $s(X)$. N.4's "for the other members of (A), (B) existence was not examined" is answered by Q15. Q2 shows that N3(b) holds without its Ricci hypothesis. No statement made before this note is contradicted.
