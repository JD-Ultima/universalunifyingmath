# O. The ℷ-variation of (9.4): fibre weights and the 4D equation

Claims: GU-054 (context), GU-052 (the ℷ-part of (9.5)). Gaps: GU-G15 (ℷ-part), GU-G14. Code: `computation/geometric-unity/O_gimel.py` → `.out` (sections O0–O9; exact: sympy rationals and rational-function fields, Fractions, Gaussian rationals; fixed seeds; 158 s on one core, a second run 157 s, identical apart from the timing annotations; 59 checks). It imports `F_jets.py` and `I_elliptic_reduction.py` (and through it `G_symbols.py`, `H_ellipticity.py`, `shiab_checks.py`), all unmodified. Record: `metrics/geometric-unity/O-gimel-variation.yaml`. Builds on F (R-Z0–R-Z3, F1–F5, F7, L-F1, L-F6), G (G1, G.2.5), D (D2, D3.4, D5), I (I1, I2, I5/I6, I7(b)), K (Theorem K), L (L.0 inputs) and N (N12(c)). Layers: (b) reconstruction, (c) proven, (d) open.

## O.0 Question and setting (layer b)

The draft treats ℷ as a field (G, g1) but writes only $dI\vert_Y$ (g2). G.2.5 left the ℷ-equation "formal; not well-posed without a fibre weight". This note makes it precise.

**Notation.** $g=ℷ$ (R-Z0), $n=\dim X$ ($=4$ unless stated), $N=n(n+1)/2$. $F_x=\mathrm{Met}_{(p,q)}(T_xX)$ with $c_VG_\lambda$, $\cong\mathbb R\times S$, $S=SL(n)/SO(p,q)$ (F1). $u:=g_x^{-1}y$, $s:=\tfrac1n\log|\det u|$. $K_x:=O(g_x)$ acts on $F_x$ by $y\mapsto P^TyP$. $dvol_F$ is the Riemannian measure of $c_VG_\lambda$; it is $GL(n)$-invariant, hence the fibre-invariant ($SL$- and scaling-invariant) measure, unique up to a constant. $\rho:=dvol_{g_\aleph}/(dvol_g\otimes dvol_F)$.

**The action.** Unitary gauge of G1, normalisation of I.0 ($b=\mathrm{Tr}/128$): $I=J(S)$, $S=\pi'(\varepsilon,\varpi)$, $J(T)=\int_YL$ with
$$L=\Psi(T\wedge F_0)+\tfrac12\Psi(T\wedge d_0T)+\tfrac13\Psi(T\wedge T\wedge T)+\tfrac\kappa2\langle T,T\rangle\,\mathrm{vol},\qquad \circledast=xW_A+yW'_A+zW_C\mathrm{Ad}_\Gamma+wW'_C\mathrm{Ad}_\Gamma .$$
$\chi:=-\Gamma^2$ ($+1$ in (14,0), $-1$ in (7,7)); real forms as in I.0. Invariant sector: $T_o=t\Phi^1+u\Phi^1\Gamma$; put
$$a:=xt+\chi yu,\qquad \tilde P(t,u):=(t^2+\chi u^2)\big(\tfrac\kappa2+104\,a\big).$$
$\sigma_1=-\tfrac{23}2-\tfrac{1+\lambda}{1-4\lambda}$ (F4(b)), $s_F=-9/c_V$; $\Lambda_A^{Z2}=\tfrac{9c_H}{2c_V}$, $\Lambda_A^{Z1}=-\tfrac{c_H\sigma_1}{2c_V}=\tfrac{5c_H}{4c_V}\tfrac{5-18\lambda}{1-4\lambda}$, and $\Lambda_P$ as in Cor. F5.

**Definition O0.**
- A *fibre weight* is a family $\nu=(\nu_x)$ of Borel measures on $F_x$. It is *natural* if $\nu_x$ is obtained from $g_x$ alone by a rule that commutes with linear isomorphisms (the weights "built from $g_x$, $y$ and invariants such as $s$, $\mathrm{tr}\,u^k$"). A natural weight is $K_x$-invariant: a local diffeomorphism fixing $x$ with differential in $K_x$ preserves $g_x$.
- Examples: $w(u)\,dvol_F$ with $w$ a class function of $u$; the *ray weights* $\nu^\mu_x:=(\tau\mapsto e^\tau g_x)_*\mu$ for a measure $\mu$ on $\mathbb R$ (for $p=q$ also $(\tau\mapsto-e^\tau g_x)_*\mu_-$).
- Weighted action: $I_\nu:=\int_X\big(\int_{F_x}\ell\rho\,d\nu_x\big)dvol_g$, where $L=\ell\,dvol_{g_\aleph}$. For $\nu=w\,dvol_F$ this is $\int_YLw$.
- $\mathcal E_ℷ:=\delta I_\nu/\delta ℷ$ for compactly supported $\deltaℷ$ on $X$. Such a variation is not compactly supported on $Y$, which is why a weight is needed for ℷ and not for $\varpi$.
- Two systems: **(S1)** every field equation from $I_\nu$; **(S2)** $\Upsilon=0$ from local unweighted variations (the draft's $dI\vert_Y$, g2), together with $\mathcal E_ℷ=0$.

## O.1 Results (layer c)

**Proposition O1 (how ℷ enters; proven).**

| ingredient | Z1 | Z2 |
|---|---|---|
| $g_\aleph$ | only through $\aleph$ (Γ, the 1-jet; F2.4); invariant under $ℷ\to cℷ$ | $ℷ_x$ itself (block $c_Hℷ$) and $\aleph$ (in $D$) |
| $dvol_{g_\aleph}$ | **independent of ℷ** (O2) | through $\lvert\det ℷ\rvert^{1/2}$ |
| $*$, $\Phi^1$ (soldering $\Theta_\aleph$), hence ⊛ | through $D$ | through $D$ and $ℷ_x$ |
| $F_0$, $d_0$, $\pi'$ | through $\nabla^0$; $F_0$ involves $\nabla\mathrm{Rm}_ℷ$ (O'Neill's $\nabla A$ terms) | same |
| $s_Y$ | 2-jet only: $\sigma_1/c_V+\mathrm{tr}(y^{-1}\mathrm{Ric})/c_H+q_y$ (O8) | 2-jet only: $R/c_H+s_F-\lvert A\rvert^2$ (F3, L-F1) |
| bundle of $\varpi$ | $(C,g_C)$ ℷ-independent: "fixed $\varpi$" is canonical | $g_C$ (or $g_\aleph$) depends on $ℷ_x$: needs an identification (O12) |

Every term of $L$ depends on ℷ ($*$ in all, and vol in Z2 only; $F_0$ in the first; $d_0$ in the second; ⊛ through $\Phi^1$). In Z2 the pulled-back metric is ℷ: R-Z0 has no other metric on $X$. With a fixed background $g_0$ in the HH block (not in the draft) Z2 would depend on ℷ only through ℵ.

*Proof.* R-Z2, F2.4 and O2. The jet orders: Γ is the 1-jet, and the curvature of $g_\aleph$ needs $\partial^2g_\aleph\ni\partial^2\Gamma$. For $s_Y$ use O8 and L-F1. ∎

**Proposition O2 (volume forms; proven, computed §O2).**
- (a) The Jacobian of $y\mapsto P^TyP$ on $\mathrm{Sym}^2$ is $(\det P)^{n+1}$. Hence $\mu_{\rm can}:=|\det y|^{-n/2}d^nx\,d^Ny$ is a Diff(X)-natural measure on Met(X).
- (b) $\det G_\lambda|_y=(1-n\lambda)2^{N-n}(\det y)^{-(n+1)}$ in the basis $E_{ii}$, $E_{ij}+E_{ji}$.
- (c) $dvol_{g_\aleph}=|c_H|^{n/2}|\det h_H|^{1/2}dvol_F\,d^nx$, independent of Γ. So Z1: $dvol_{g_\aleph}=C\mu_{\rm can}$ with $C=|c_H|^{n/2}|c_V|^{N/2}|1-n\lambda|^{1/2}2^{(N-n)/2}$, independent of ℷ; Z2: $=Ce^{-ns/2}\mu_{\rm can}$. Thus $\rho_{Z1}=|c_H|^{n/2}e^{ns/2}$ and $\rho_{Z2}=|c_H|^{n/2}$.
- (d) $y\mapsto e^\tau y$ is an isometry of Z2, and maps Z1 with $c_H$ to Z1 with $e^\tau c_H$.

*Proof.* (a) Exponent $-(n+1)+1+n=0$. (b) GL-equivariance from $y=I$. (c) The coframe $(dx,Dy)$ is unitriangular with respect to $(dx,dy)$. (d) $G_{cy}(ch,ck)=G_y(h,k)$ and $D(cy)=c\,Dy$. ∎

**Proposition O3 (T = 0; proven, §O9).**
- (a) $L\vert_{T=0}\equiv0$ for every ℷ, reading and ⊛. In the unitary gauge $\mathcal E_ℷ=0$ at $T=0$ for every weight: **the ℷ-equation at $T=0$ is empty.**
- (b) At fixed $(\varepsilon,\varpi)$, $\partial_ℷ\pi'=-(1-\mathrm{Ad}_\varepsilon)\delta A_0$. So $\delta_ℷI(\varepsilon,\varpi)=\partial_ℷJ(S)-p_\nu\big((1-\mathrm{Ad}_\varepsilon)\delta A_0,E_J(S)\big)$. Off-shell this changes along a $\tau'$-orbit, e.g. at $T=0$, $\varepsilon\neq1$ it is $-p_\nu((1-\mathrm{Ad}_\varepsilon)\delta A_0,\circledast F_0)$. It is orbit-independent on $\Upsilon=0$. Below, $\mathcal E_ℷ$ is taken at fixed $S$.
- (c) $\frac{d}{d\tau}\big|_0\mathcal E_ℷ(\tau\hat T)=\mathcal E_ℷ\big[p_\nu(\hat T,\circledast F_0)\big]$. For $\hat T=t\Phi^1+u\Phi^1\Gamma$ this is $a$ times the ℷ-variation of the weighted Einstein–Hilbert integral $\int6s_Y\rho\,d\nu$ (O4).

*Proof.* (a) Every term of $L$ has a factor $T$ (G, R-G0). (b) $d_0\varepsilon=d\varepsilon+A_0\varepsilon-\varepsilon A_0$, then the chain rule with $E_J$ the gradient (G2). (c) $L=\Psi(T\wedge F_0)+O(T^2)$. ∎

**Proposition O4 (the density on the invariant sector; proven, computed §O0–§O1).** Let $T_o=t\Phi^1+u\Phi^1\Gamma$ and $T_e=0$, with $t,u$ functions. Pointwise
$$\ell=6a\,s_Y+14\tilde P(t,u),\qquad \partial_t\ell=6xs_Y+14c_1,\quad \partial_u\ell=\chi(6ys_Y+14c_{13}),$$
with $c_1,c_{13}$ of I6. The kinetic term $\tfrac12\Psi(T\wedge d_0T)$ vanishes pointwise.

*Proof.*
- $\Psi(T\wedge F_0)=-\langle T,xG\,e\otimes\gamma+yG\,e\otimes\gamma\Gamma\rangle$ (D3.4) $=6a\,s_Y$, since $\mathrm{tr}\,G=-6s_Y$ and $\langle\Phi^1\Gamma,\Phi^1\Gamma\rangle=14\chi$.
- $\langle T,T\rangle=14(t^2+\chi u^2)$, and $\tfrac13\Psi(T\wedge T\wedge T)=14\cdot104\,a(t^2+\chi u^2)$; the twisted operators contribute 0.
- $d_0T=dt\wedge\Phi^1+du\wedge\Phi^1\Gamma$, whose image under ⊛ lies in grades 2, 4, 12 (I1, I5); it pairs to zero with grades 1, 13.
- Checks: $c_1,c_{13}$ are recomputed by I's `i5_invariant`; $(c_1,\chi c_{13})=\nabla\tilde P$; direct Clifford evaluation of each term at 6 points per operator in both real forms, the curvature term for a random algebraic curvature tensor. ∎

**Proposition O5 (divergence for the fibre-invariant measure; proven).** Fix $x$ and constant $(t,u)$. Take $\nu_x=w\,dvol_F$ with $w=1$, i.e. the unweighted $\int_YL$, or, in Z1, $w=e^{-ns/2}$, which puts the fibre-invariant measure in place of Z1's own volume.
- (a) Along each scaling line $\tau\mapsto e^\tau y_1$, $\ell\rho w$ is constant (Z2) or an exponential polynomial $\sum_kc_ke^{k\tau}$ (Z1, $n=4$), with $k\in\{0,1,2\}$ for $w=1$ and $k\in\{0,-1,-2\}$ for $w=e^{-2s}$. Hence $\int_{F_x}|\ell\rho|\,w\,dvol_F<\infty$ iff $\ell\equiv0$ on $F_x$.
- (b) If $a\neq0$, $\ell\equiv0$ on $F_x$ requires:
  - Z2: $|A|^2(x,\cdot)\equiv0$ (in definite signature: $\mathrm{Rm}(x)=0$) and $6a(R/c_H+s_F)+14\tilde P=0$;
  - Z1: $\mathrm{Ric}(x)=0$, $q_y(\mathrm{Rm}(x))\equiv0$ and $6a\sigma_1/c_V+14\tilde P=0$.

  If $a=0$, it requires $\tilde P=0$. So for generic $(t,u)$ the fibre integral diverges at every $x$.
- (c) At $T=0$ the action is $0$ for every weight. Its $T$-derivative, the weighted EH integral $\int_Fs_Y\rho\,dvol_F$, diverges at every $x$ in Z2: $s_Y\equiv0$ on $F_x$ would need $|A|^2\equiv0$, which for $n=4$ forces $\sum\epsilon_i\epsilon_jR_{ij}\otimes R_{ij}=0$ (the kernel of $B\mapsto f_B$ is $\mathbb R\star$, and Bianchi excludes $c\star$, $c\neq0$), hence $R(x)=0$, and $s_F\neq0$ (REVIEW9-GU R11-1). In Z1 it diverges unless $\lambda=5/18$ ($\sigma_1=0$), $\mathrm{Ric}(x)=0$ and $q\equiv0$.
- (d) The same holds for $t,u$ depending on $(x,s)$ only (the radial ansatz of J5, with $T_e=0$): $\int_Sdvol_S=\infty$. Proven for flat bases, where the integrand is constant on $S$, and in Z2 for definite bases, where O9 and the triangle inequality on $K$-orbits apply. In indefinite signature $|A|^2$ can vanish identically for $\mathrm{Rm}\neq0$ (e.g. $R=\omega\otimes\omega$ with ω a null simple bivector), and for general curved indefinite bases the divergence is not shown. (Corrected after review 2026-10-30, R11-1.)

*Proof.*
- (a) Z2 by O2(d). Z1 by O2(d) and O8, using $q_{cy}=c^{-2}q_y$. Then Fubini on $\mathbb R\times S$, and $\int_{\mathbb R}|\sum_kc_ke^{k\tau}|\,d\tau=\infty$ unless every $c_k=0$.
- (b) Z2: $|A|^2(x,g_x)=0$ and $|A|^2$ is constant on $F_x$. In definite signature $|A|^2\equiv0$ forces $R\cdot y=0$ for all $y$, i.e. $\mathrm{Rm}=0$ (as in L2). Z1: by O8, $\mathrm{tr}(y^{-1}\mathrm{Ric})\equiv0$ forces $\mathrm{Ric}=0$.
- (c), (d) follow. ∎

**Theorem O6 (indefinite signature: no natural convergent weight; proven).** Let $pq>0$, and let $\nu$ be a signed Borel measure of finite total variation on $F_x$, invariant under $SO_0(g_x)$. Then $|\nu|$ is concentrated on $\mathrm{Fix}:=\{cg_x:c\in\mathbb R^\times\}\cap F_x$: the ray $\mathbb R_{>0}g_x$, plus $\mathbb R_{<0}g_x$ when $p=q$. Consequently:
- (a) Every $SO_0(g_x)$-invariant $w\in L^1(F_x,dvol_F)$ vanishes a.e. In particular **no natural weight function** (built from $g_x$, $y$, e.g. from $s$ and $\mathrm{tr}\,u^k$) has finite fibre integral, unless it is 0 a.e.
- (b) Suppose a natural $\nu$ gives an absolutely convergent fibre integral $\int_{F_x}|\ell\rho|\,d|\nu_x|<\infty$ at one point of one flat base and one constant $(t,u)$ with $6as_F+14\tilde P\neq0$ (Z2), resp. $6a\sigma_1/c_V+14\tilde P\neq0$ (Z1). Then $\nu$ is a ray weight.

*Proof.*
- Let $X\in\mathfrak{so}(g_x)$ be a boost. Its action $y\mapsto X^Ty+yX$ on $\mathrm{Sym}^2$ is diagonalisable with real eigenvalues in $\{0,\pm\beta,\pm2\beta\}$ (§O6). Write $\mathrm{Sym}^2=V_+\oplus V_0\oplus V_-$ accordingly.
- $b:=\exp X$ preserves $F_x$ and $|\nu|$ (uniqueness of the Jordan decomposition). Put $h(y):=|\pi_+y|$; then $h(b^ky)\ge e^{k\beta}h(y)$.
- A point of $A_c:=\{c\le h<2c\}$ leaves $A_c$ for good after $\lceil\log2/\beta\rceil$ steps. By Poincaré recurrence $|\nu|(A_c)=0$, so $|\nu|\{\pi_+y\neq0\}=0$. With $b^{-1}$, also $|\nu|\{\pi_-y\neq0\}=0$.
- So $|\nu|$ lives on $\ker X$ for every boost. The boosts span $\mathfrak p$; a $y$ killed by $\mathfrak p$ is killed by $[\mathfrak p,\mathfrak p]=\mathfrak k$, hence by $\mathfrak{so}(g)$, and $\mathrm{Sym}^2_0$ has no invariants. So the joint kernel is $\mathbb Rg_x$ (also computed, §O6, for (1,1), (2,1), (1,2), (1,3), (3,1), (2,2), (2,3)).
- (a) $|w|\,dvol_F$ is finite and invariant, and Fix is $dvol_F$-null.
- (b) For a flat base, $\ell\rho$ is a nonzero constant (Z2) or a nonzero multiple of $e^{ns/2}$ (Z1), hence $K_x$-invariant. So $|\ell\rho|\,|\nu|$ is a finite invariant measure. ∎

The mechanism: $K_x$ is non-compact, and regular $K_x$-orbits in $F_x$ have discrete (in fact finite) stabiliser (§O6), so invariance spreads any mass over infinite Haar volume. Gaussian or Wishart weights in $\mathrm{tr}\,u$ fail for a second reason: on $F_x$ for (1,3), $\mathrm{tr}\,u$ is unbounded below and $u$ can have complex eigenvalues (§O6).

**Corollary O7 (classification of convergent natural weights; proven).**
- (a) **Indefinite base**, i.e. every (7,7) case: (1,3), (3,1), (2,2). The natural weights for which the invariant-sector action converges for all bases and all constant $(t,u)$ are exactly the ray weights with $\mu(\mathbb R)<\infty$ (Z2), resp. $\int(e^\tau+e^{2\tau})\,d|\mu|<\infty$ (Z1, $n=4$).
- (b) **Definite base**, the (14,0) case. Natural means $O(4)$-invariant; an absolutely continuous natural weight is $\nu=W(a)\prod_{k<l}|\sinh\tfrac{a_k-a_l}2|\,da\otimes dk$, where $u=k\,\mathrm{diag}(e^{a})k^{-1}$ (§O6 for $n=2$).
  - Z2: convergent iff $\int(1+\mathrm{tr}\,u\,\mathrm{tr}\,u^{-1})\,d|\nu|<\infty$.
  - Z1: necessary $\int e^{2s}(1+\mathrm{tr}\,u^{-1})\,d|\nu|<\infty$; sufficient $\int e^{2s}(\mathrm{tr}\,u+\mathrm{tr}\,u^{-1})^kd|\nu|<\infty$ for all $k$.
  - Examples: compactly supported $W$; $e^{-\beta(\mathrm{tr}\,u+\mathrm{tr}\,u^{-1})}$; $e^{-|a|^2/\sigma^2}$.
  - Never: weights in $s$ alone ($\int_Sdvol_S=\infty$).

*Proof.*
- (a) By O6(b). On the ray, $\ell\rho$ is $c_H^2(6a(R/c_H+s_F)+14\tilde P)$ in Z2 and $c_H^2\big(6ae^{\tau}R/c_H+e^{2\tau}(6a\sigma_1/c_V+14\tilde P)\big)$ in Z1 (O4, O2(c), F3, F4(b), O8).
- (b) Z2: $\ell\rho=c_H^2(6a(R/c_H+s_F-|A|^2)+14\tilde P)$ with $|A|^2\ge0$; then O9 and Tonelli.
- Z1: $s_Y$ is rational in $y$ with denominators powers of $\det y$, so $|s_Y|\le C(\mathrm{tr}\,u+\mathrm{tr}\,u^{-1})^k$. Necessity follows from flat bases and from the $\varepsilon$-linear coefficient under $\mathrm{Rm}\to\varepsilon\mathrm{Rm}$ (O8, with $\mathrm{Ric}\propto g$ at $x$). ∎

**Lemma O8 (Z1 scalar curvature off the section; proven, computed §O3).** For every $n$, signature and $\lambda\neq1/n$, at $(x,y)\in Y$:
$$s_Y^{Z1}(x,y)=\frac{\sigma_1}{c_V}+\frac{\mathrm{tr}(y^{-1}\mathrm{Ric}_ℷ(x))}{c_H}+q_y\big(\mathrm{Rm}_ℷ(x)/c_H\big).$$
Here $q_y$ is a quadratic form depending on $(y,c_V,\lambda)$ only, homogeneous of degree $-2$ in $y$, and $q_{cg_x}=0$.

*Proof.*
- In ℵ-normal coordinates at $x$, $\Gamma(x)=0$ and $\partial\Gamma(x)$ is linear in $\mathrm{Rm}(x)$.
- The HH block $c_Hy$ and $\det g_\aleph$ (O2(c)) do not depend on $x$, and at $\Gamma=0$ the metric is block-diagonal. So no term of the scalar curvature contains two $x$-derivatives of a mixed component, i.e. $\partial^2\Gamma$ does not occur.
- Hence $s_Y$ is a polynomial of degree $\le2$ in $\partial\Gamma(x)$, GL-invariant in $(y,\mathrm{Rm})$.
- Degree 0: for flat ℵ every constant $y$ is a parallel section, so F4(b) gives $\sigma_1/c_V$.
- Degree 1: the only GL-invariant linear function of $R^a{}_{bcd}$ (with $R^a{}_{acd}=0$) built from $y$ is a multiple of $\mathrm{tr}(y^{-1}\mathrm{Ric})$, and F4(b) at $y=g$ fixes the coefficient $1/c_H$.
- Rescaling $x$ by $\sqrt{c_H}$ gives the dependence on $\mathrm{Rm}/c_H$, and O2(d) gives $q_{cy}=c^{-2}q_y$.
- On the ray, $cg$ is a parallel section, and F4(b) for the base metric $cg$ gives $\sigma_1/c_V+R/(cc_H)$. ∎

Checks: symbolic ε in $g=g_0+\varepsilon(\text{quad}+\text{cub})$, for $n=3$ in (3,0), (2,1), (1,2) and $n=4$ in (1,3). The result has degree 2 exactly, constant $\sigma_1/c_V$, linear coefficient $\mathrm{tr}(y_0^{-1}\mathrm{Ric})/c_H$, and quadratic coefficient nonzero off the ray and zero on it. The cubic terms of $g$ do not matter, and doubling $c_H$ halves the linear coefficient and quarters the quadratic one.

**Lemma O9 (Z2, definite base: the A-tensor on the fibre; proven, computed §O5).** At $u=\mathrm{diag}(e^{a_k})$ in a $g$-orthonormal frame,
$$|A|^2=\frac{c_V}{4c_H^2}\sum_{ijkl}R_{klij}^2\cdot4\sinh^2\tfrac{a_k-a_l}2 .$$
For $n=4$ and every $O(4)$-invariant measure ν,
$$\int|A|^2d\nu=\frac{c_V}{24c_H^2}\,|\mathrm{Rm}|^2\int(\mathrm{tr}\,u\,\mathrm{tr}\,u^{-1}-16)\,d\nu .$$

*Proof.*
- $A_XY=-\tfrac12(R(X,Y)\cdot y)^{\rm vert}$ (L-F5). With $R\in\mathfrak{so}(g)$, $y^{-1}R^Ty=-u^{-1}Ru$, which gives $G_y(R\cdot y,R\cdot y)=\sum_{kl}R_{kl}^2(u_k/u_l+u_l/u_k-2)$.
- Averaging: Weyl, $\mathrm{Ric}_0$ and scalar parts are pairwise inequivalent irreducible $O(4)$-modules. So the $O(4)$-average of $Q_u$ is $\sum_V\frac{\mathrm{tr}\,Q_u|_V}{\dim V}|P_VR|^2$, and the three ratios equal $\sum_{k\neq l}f_{kl}/12$, with $\sum_{k\ne l}f_{kl}=2(\mathrm{tr}\,u\,\mathrm{tr}\,u^{-1}-16)$.
- Equivalently, the average of $\xi\otimes\xi$ over unit simple bivectors is $\mathrm{id}/6$; the Hodge part drops because $\xi\wedge\xi=0$. ∎

Checked: brute force at $n=3$; exact Schur traces at three random $u$ for $n=4$.

**Theorem O10 (the ℷ-equation on the invariant sector; system (S2); proven, computed §O7).** Take $n=4$, $T=t\Phi^1+u\Phi^1\Gamma$ with constant $t,u$ and $T_e=0$, varied at fixed $(t,u)$, with $a\neq0$. Put
$$\Lambda_O:=\Lambda_A-\tfrac{7c_H}6\,\frac{\tilde P}a=\Lambda_A-c_H(t^2+\chi u^2)\Big(\frac{7\kappa}{12a}+\frac{364}3\Big).$$
- **(a) Ray weights** (forced in every (7,7) case by O7(a)). $I_\nu=\int_X(\alpha R+\beta)\,dvol_g$ with $\alpha=6ac_Hm_1$ and $\beta=c_H^2m_2(6a\sigma/c_V+14\tilde P)$:
  - Z2: $m_1=m_2=\mu(\mathbb R)$, $\sigma=-9$;
  - Z1: $m_k=\int e^{k\tau}d\mu$, $\sigma=\sigma_1$.

  So $\mathcal E_ℷ=0$ is **exactly the 4D vacuum Einstein equation** $G(ℷ)+\Lambda ℷ=0$, with
  $$\text{Z2: }\Lambda=\Lambda_O^{Z2}\ \text{(independent of }\mu),\qquad \text{Z1: }\Lambda=r_\mu\Lambda_O^{Z1},\quad r_\mu=m_2/m_1 .$$
  This needs $\alpha\neq0$, i.e. $m_1\neq0$; if $m_1=0$ there is no Einstein term and $\mathcal E_ℷ=\tfrac\beta2g$. For $\mu\ge0$ on the positive ray, $r_\mu$ takes every value in $(0,\infty)$ ($\mu=\delta_{\tau_0}$ gives $e^{\tau_0}$). For $p=q$, in Z1, $m_1=\int e^\tau d\mu_+-\int e^\tau d\mu_-$ and $m_2=\int e^{2\tau}d(\mu_++\mu_-)$ (the negative ray has $s_Y=-e^{-\tau}R/c_H+\sigma_1/c_V$). So $r_\mu$ takes every nonzero real value ($\mu_-=\delta_{\tau_0}$, $\mu_+=0$ gives $-e^{\tau_0}$). In Z2, $m_1=m_2=(\mu_++\mu_-)(\mathbb R)$ (R11-2).
- **(b) Smooth natural weights, definite base.**
  - Z2: $\mathcal E_ℷ=0$ reads
    $$G+\Lambda_O^{Z2}g+\kappa_wE_{|\mathrm{Ric}_0|^2}=0,\qquad \kappa_w=\frac{c_V}{6c_H}\big\langle\mathrm{tr}\,u\,\mathrm{tr}\,u^{-1}-16\big\rangle_w\ge0,$$
    where $E_{|\mathrm{Ric}_0|^2}$ is the Euler–Lagrange tensor of $\int|\mathrm{Ric}_0|^2$ (fourth order) and $\langle\cdot\rangle_w$ the $w$-average. $\kappa_w>0$ for every positive absolutely continuous weight ($\mathrm{tr}\,u\,\mathrm{tr}\,u^{-1}\ge16$, with equality only on the ray), and $\kappa_w\to0$ under localisation onto the ray. $\kappa_w$ does not depend on $(t,u)$.
  - Z1: Einstein metrics solve it iff $\Lambda=r_w\Lambda_O^{Z1}$, with $r_w=\int e^{2s}d\nu\big/\int\tfrac14\mathrm{tr}(u^{-1})e^{2s}d\nu$. The quadratic-curvature term is not computed.
- **(c)** If $a=0$, $\mathcal E_ℷ=0$ reads $\tilde P\,g=0$ (ray weights, Z2).

In (a) and (b) every Einstein metric with that Λ is a solution; in (a) these are all the solutions.

*Proof.*
- (a) On the ray, $s_Y=R/c_H+s_F$ in Z2 ($A=0$ there, by F3 and the isometry O2(d)) and $s_Y=e^{-\tau}R/c_H+\sigma_1/c_V$ in Z1 (F4(b) for the parallel section $e^\tau g$; O8). Take ρ from O2(c) and the density from O4. The composite is a functional of ℷ alone, including the moving ray. L-F6: $\delta\int(\alpha R+\beta)=-\alpha\int\langle G+\Lambda g,\delta g\rangle$ with $\Lambda=-\beta/(2\alpha)$.
- (b) Use O9 and, in dimension 4, $|\mathrm{Rm}|^2=E_{GB}+4|\mathrm{Ric}_0|^2$; the Gauss–Bonnet integrand has identically vanishing EL tensor. $E_{|\mathrm{Ric}_0|^2}$ has a factor $\mathrm{Ric}_0$ or a derivative of it in every term, so it vanishes on Einstein metrics.
- Z1: by O8 and $K_x$-invariance, $\int u^{-1}\rho\,d\nu=\tfrac14\int\mathrm{tr}(u^{-1})\rho\,d\nu\cdot I$. The $O(4)$-invariant quadratic forms on curvature tensors are spanned by $|W|^2$, $|\mathrm{Ric}_0|^2$, $R^2$, and in dimension 4 their EL tensors vanish on Einstein metrics (Bach tensor; $\delta\int R^2=2R(\mathrm{Ric}-\tfrac14Rg)+$ derivatives of $R$; Besse, *Einstein Manifolds*, §4.H).
- (c) is immediate. ∎

**Corollary O11 (leading order, localisation, comparison with F; proven, computed §O7).**
- (a) For every admissible weight, $\mathcal E_ℷ(\tau T)=\tau a\,\mathcal E_ℷ^{EH}+O(\tau^2)$. The **linearised ℷ-equation is the weighted Einstein–Hilbert equation of $g_\aleph$.** As $T\to0$ (with $a\neq0$), $\Lambda_O\to\Lambda_A$, with first-order correction $-\tfrac{7c_H\kappa}{12}\tfrac{t^2+\chi u^2}a$.
- (b) *Localisation at $s(X)$* ($\mu=\delta_0$; or smooth $w\to\delta_{g_x}$ in definite signature, where $\kappa_w\to0$). The leading-order ℷ-equation is **F's reading (A)**, $G+\Lambda_Ag=0$ with $\Lambda_A$ of Cor. F5, in Z1 and in Z2. In Z1, $\Lambda_A\neq\Lambda_P$ (λ = ½: $5$ against $21/4$, in units of $c_H/c_V$). So the ℷ-variation selects (A), not (P).
- (c) *Weight-independent statements.*
  - Z2: $\Lambda=\Lambda_O^{Z2}$ for every ray weight, and for the Einstein solutions under every smooth definite weight.
  - Z1: the weight enters only through $r_\mu$. This is equivalent to the fibre scale $c_H/c_V\to r_\mu c_H/c_V$, which the draft leaves free (F.0). For $\mu\ge0$ on the positive ray (for $p\neq q$, every $\mu\ge0$), the sign of Λ does not depend on the weight. For $p=q$ a weight on the negative ray can reverse it (R11-2). The zero $\lambda=5/18$ (at leading order) does not depend on the weight whenever $m_1\neq0$.
- (d) *Against (P) and (E).*
  - (P) of I7(a): $\Lambda_4=\Lambda_P-c_Hc_1/x$, first order $-c_H\kappa t/x$. The two differ at first order in $T$, also in Z2 where $\Lambda_P=\Lambda_A$.
  - (E) requires $g_\aleph$ Einstein on $s(X)$ (F6). The ℷ-equation neither implies (E) nor follows from it.
- (e) $\tilde P=0$, hence $\Lambda_O=\Lambda_A$ exactly, on $a=-\kappa/208$, which is I6's $W_A$ branch $t=-\kappa/(208x)$ for $y=0$, and in (7,7) also on $t=\pm u$.

**Proposition O12 (identifications; proven, computed §O4, §O7).** The ℷ-variation at fixed ϖ needs an identification ψ of the configuration spaces for different ℷ. Then $\mathcal E_ℷ^\psi=\mathcal E_ℷ^{(t,u)}+p_\nu(\delta_\psi T,\Upsilon)$, and the difference vanishes where $\Upsilon=0$.
- (a) Z1, canonical identification through $\Theta_\aleph$ into the ℷ-independent $(C,g_C)$ (or Bourguignon–Gauduchon for $g_\aleph$): $\delta\Theta=(\delta D,0)$ is HV-mixed. It pairs with $E^{(i)}_{HV}=x_i\mathrm{Ric}^Y_{HV}$, which is $0$ on the ray (§O4; Z1 at $cg$ is the section $cg$, F4(c)). So for ray weights, fixed ϖ = fixed $(t,u)$.
- (b) Z2, Bourguignon–Gauduchon: $\delta_\psi T=-\tfrac12\delta g_\aleph\cdot(t\,e\otimes\gamma+u\,e\otimes\gamma\Gamma)$. On the ray $\delta g_\aleph$ has HH block $c_H\deltaℷ$, and $\mathrm{Ric}^Y_{HH}=\mathrm{Ric}_ℷ$ (§O4). Then
  $$G+\Lambda'g=0,\qquad \Lambda'=\Lambda_A^{Z2}-\frac{c_H(14\tilde P-D)}{11a}=\Lambda_A-c_H(t^2+\chi u^2)\Big(\frac{6\kappa}{11a}+104\Big),\quad D=tc_1+\chi uc_{13}.$$
  $\Lambda'=\Lambda_O$ at leading order ($\Lambda_A$), but not at first order.

**Proposition O13 (system (S1): the weight in the ϖ-equation; proven, computed §O8).**
- (a) For a smooth weight $w$, $\Upsilon_w:=\delta I_w/\delta\varpi=w\Upsilon+\tfrac12\circledast(dw\wedge T)$. Equivalently $\tfrac12w\Psi(T\wedge d_0T)=\tfrac12\Psi(\tilde T\wedge d_0\tilde T)$ with $\tilde T=w^{1/2}T$. The ε-equation stays redundant, since $I_w=J_w\circ\pi'$ pointwise.
- (b) On $T_o=t\Phi^1+u\Phi^1\Gamma$, $T_e=0$, the grade-2 part of $\Upsilon_w=0$ is $d(w^{1/2}a)=0$ (unweighted: $da=0$). For constant $(t,u)$ and $dw\neq0$, this forces $a=0$.
- (c) Ray weights. For every conormal $n$ of the cone $C=\{e^\tau g_x\}$ (dimension 5, codimension 9), the ϖ-equation of $I_\nu$ contains $\circledast(n^\flat\wedge T)=0$ on $C$. For $T_o=t\Phi^1+u\Phi^1\Gamma$ and any even $T_e$, its grade-2 part is $a\,W_A(n^\flat\wedge\Phi^1)_2\neq0$. So $a\vert_C=0$, the coefficient of $s_Y$ vanishes on $C$, and **within (S1) the invariant sector gives no 4D gravitational equation.**

*Proof.*
- (a) The kinetic term is the only derivative of $T$ in $L$. $\Psi(T\wedge\xi\wedge T)=0$ because the symbol of a self-adjoint ⊛ is antisymmetric (D2; checked for all 48 types, §O8). Then Leibniz; a toy Euler–Lagrange check is in §O8.
- (b) With $T_e=0$, the grade-2 part of $\Upsilon_{\rm ev}$ (I2) is $\circledast(d_0T_o)_2$. By §O8 the grade-2 part of $\circledast(\xi_1\wedge\Phi^1+\xi_2\wedge\Phi^1\Gamma)$ is $W_A((x\xi_1+\chi y\xi_2)\wedge\Phi^1)_2$, and $\xi\mapsto W_A(\xi\wedge\Phi^1)_2$ is injective (I5). With the weight, $(dt,du)$ becomes $w^{-1/2}d(w^{1/2}(t,u))$.
- (c) The variation of $\tfrac12\int\Psi(T\wedge d_0T)\,d\nu$ contains $\tfrac12\int\Psi(T\wedge n^\flat\wedge\nabla_n\alpha)\,d\nu$ with $\nabla_n\alpha\vert_C$ free, and $\Psi(T\wedge n\wedge\beta)=-p(\beta,\circledast(n\wedge T))$.
- $T_e$ has grades 0, 2, 12, 14, which the four operators send to grades other than 2 (I1(a)).
- The conormals include non-null covectors, since $g_\aleph\vert_{V_0}$ is nondegenerate. ∎

**Remark O14 (constant-mode system; computed §O7).** In (S2) the constant modes are not determined: $\Upsilon=0$ fails in this sector at every point of $s(X)$ except in Z1 at $\lambda=\tfrac15$ over a Ricci-flat base (F6), on every open set in Z2 (L2(a)), and near $s(X)$ in Z1 except for the flat-base solutions of O15 (N10, N12). Suppose one adds the minisuperspace equations $\partial_tI_\nu=\partial_uI_\nu=0$ (ray weight, Z2, $y=0$). Their solutions are:
- $u=0$, $t^2=-27/(728c_V)$, $\Lambda=9c_H/c_V-7c_H\kappa t/(12x)$;
- $t=-\kappa/(208x)$ with $u^2=\kappa^2/(43264x^2)+27/(728c_V)$ in (7,7), resp. $u^2=-\kappa^2/(43264x^2)-27/(728c_V)$ in (14,0), and $\Lambda=\Lambda_A^{Z2}$.

For $c_H,c_V>0$ (the (7,7) sign rule A4 with $c_H>0$), only the second branch is real in (7,7), so $\Lambda=9c_H/(2c_V)$, the anti-de Sitter sign (F5). This truncation is not the restriction of either system: symmetric criticality is not available here.

**Corollary O15 (an on-shell test of (S2); proven, computed §O7).** The flat-base solutions of round 10 are genuine solutions of $\Upsilon=0$ on all of $Y$: Z1, base (4,0), $\lambda=\tfrac15$, $c_V\kappa^2\ge9360x'^2$, $u=0$, $T_e\equiv0$, $312x's^2+\kappa s+\tfrac{15x'}{2c_V}=0$ (N12(c), I7(b); for $\circledast_\sharp$ and $x W_A$, i.e. $y=0$; $t=is$, $x=-ix'$). There O12 makes $\mathcal E_ℷ$ identification- and orbit-independent. At a flat base it reduces, for every admissible weight (ray or smooth), to $\beta=0$, i.e. $\Lambda_O^{Z1}=0$, with $\Lambda_A^{Z1}(\tfrac15)=\tfrac{35c_H}{4c_V}$. Together with $\Upsilon=0$ this holds **iff $c_V\kappa^2=12480\,x'^2$ and $s=-\kappa/(416x')$**, where $\det M=-\tfrac{24}{169}\kappa^4\neq0$. So among the N12(c) solutions exactly those with this one value of $c_V\kappa^2$ solve (S2) on all of $Y$. They are flat ($\Lambda_4=0$) and Euclidean. In (7,7) none is known: Z1 has none near $s(X)$ (N12), and in Z2 this sector has none on any open set for $(x,y)\neq0$ (L2(a): Z2 is Einstein on no open set).

*Proof.* On a flat base $s_Y\equiv\sigma_1/c_V$ on $Y$, so the constant term of $I_\nu$ is $c_H^2(6a\sigma_1/c_V+14\tilde P)\int\rho\,d\nu$. The Einstein–Hilbert and quadratic terms have vanishing EL tensors at a flat metric. Eliminating between $\kappa t+312xt^2=\tfrac{15x}{2c_V}$ and $\kappa t+208xt^2=\tfrac{15x}{c_V}$ gives the result. ∎

**Answers.**
1. *Dependence.* Z1 depends on ℷ through ℵ only, and its volume form does not depend on ℷ at all. Z2 depends on $ℷ_x$ and ℵ. Every term of (9.4) carries ℷ, but $L$ vanishes at $T=0$ (O1–O3).
2. *Integration.*
   - With the fibre-invariant measure, every nonzero fibre integrand of this sector diverges (O5).
   - In every (7,7) case no natural weight function converges; the convergent natural weights are measures on the scaling ray through $s(X)$ (O6, O7(a)).
   - In (14,0) a large class converges (O7(b)).
3. *The equation.*
   - At $T=0$ it is empty.
   - At first order in $T$ it is the weighted EH equation of $g_\aleph$.
   - On the invariant sector, for ray weights, the ℷ-equation of (S2) (off-shell, at fixed $(t,u)$) is exactly 4D vacuum GR with $\Lambda=\Lambda_O^{Z2}$ (weight-independent) or $r_\mu\Lambda_O^{Z1}$, tending to $\Lambda_A$ (Z2) and $r_\mu\Lambda_A$ (Z1), i.e. reading (A). In this sector $\Upsilon=0$ itself holds near $s(X)$ only on round 10's flat-base solutions (O15, where $\Lambda_4=0$), and $T\to0$ is off-shell ($\Upsilon\vert_{T=0}\neq0$) (R11-3).
   - Under (S1) the weight also changes the ϖ-equation and forces $a=0$, which removes the Einstein term (O10, O11, O13).
   - On-shell, the only known solutions of (S2) are in Z1, (14,0), over a flat base, with $c_V\kappa^2=12480x'^2$ (O15).
4. *Weight-independent.* Z2: Λ. Z1: Λ up to the free fibre scale (for $p=q$ also its sign); localisation onto $s(X)$ gives exactly (A) in both readings.

## O.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| O4 (§O0, §O1) | $c_1,c_{13}$ from I's `i5_invariant`; gradient and potential (sympy); direct Clifford evaluation of the cubic, κ and curvature terms at 6 points; the kinetic term on random $(dt,du)$ | both real forms, all four operators / all 48 types |
| O2 (§O2) | Sym² Jacobian ($n=2$ symbolic, $n=3,4$ random); $\det G_\lambda$ (λ symbolic, $n=2,3,4$); $\det g_\aleph$ at points with $\Gamma\neq0$; scaling identities by brute-force scalar curvature | $n=2,3$, three signatures; scaling $n=3$ |
| O8 (§O3) | brute-force scalar curvature of Z1 (F_jets) with symbolic ε, off and on the ray; cubic terms; $c_H$-scaling | $n=3$: (3,0), (2,1), (1,2); $n=4$: (1,3) |
| ray values (§O4) | brute-force $s_Y$, Ricci HH and HV blocks at $y=\pm2g(x_0)$, general coordinates | $n=4$: Z1, Z2 in (1,3); Z1 at $-2g$ in (2,2) |
| O9 (§O5) | brute force against the eigen-formula ($n=3$); 20-dim space of curvature tensors, exact Schur traces | $n=3$ (3,0); $n=4$ Euclidean, three random $u$ |
| O6, O7 inputs (§O6) | boost spectra on Sym², joint kernel, orbit rank; definite case; a (1,3) point with $\mathrm{tr}\,u=-18$ and complex eigenvalues; $n=2$ Weyl Jacobian | (1,1), (2,1), (1,2), (1,3), (3,1), (2,2), (2,3); (4,0) |
| O10–O12, O14, O15 (§O7) | sympy: Λ for ray weights (Z1, Z2), smooth Z2 weights ($\kappa_w$), BG identification, closed forms, limits and first-order terms, (P) comparison, $\tilde P$ on the branch, constant-mode solutions, (14,0) reality; the on-shell elimination of O15 and $\det M$ there | both real forms; O15 in (14,0) |
| O13 (§O8) | toy EL operator (sympy); $\Psi_t(T\wedge\xi\wedge T)=0$ for all 48 types on random dense $T$; grade-2 part of $\circledast(\xi\wedge T)$ equals $a\,W_A(\xi\wedge\Phi^1)_2$ | both real forms, three members, three ξ |
| O3(b) (§O9) | 2×2 symbolic | — |
| O1, O3(a),(c), O5, O6, O7, O10–O13, O15 | proofs above | as stated |

Run: `cd computation/geometric-unity && python3 O_gimel.py > O_gimel.out 2>&1`. Result: 59 checks, 59 True, each asserted; 158 s on one core, a second run 157 s, with identical output apart from the timing annotations.

## O.3 Limits (layers b, d)

- **Sector.** The explicit ℷ-equation covers only $T_o=t\Phi^1+u\Phi^1\Gamma$ with constant $t,u$ and $T_e=0$ (O5(d) also covers $t,u$ depending on $(x,s)$, over flat bases, and in Z2 over definite bases). Not treated:
  - nonconstant $t,u$;
  - the slaved $T_e$ of rounds 6–9;
  - the symmetric sector θ (I4), whose density couples $\mathrm{Ric}^Y$ to θ.
- **Z1 off the ray.** $q_y$ is not computed, and neither is the quadratic term for smooth definite weights in Z1.
- **Class of weights.** Natural means built from $g_x$ and $y$. Weights using further structure (a timelike unit field, a background Riemannian metric, the curvature of ℷ) are outside the class: they break Diff(X)-naturality or add curvature terms. Conditionally convergent signed weights and Wick rotation between (4,0) and (1,3) are not considered.
- **The draft fixes neither** the weight nor the system ((S1) or (S2)), nor the identification of configuration spaces for different ℷ. In Z2 the ℷ-equation off-shell depends on that identification beyond leading order (O12).
- **Joint solutions.** In Z2, with $T_e=0$, $\Upsilon=0$ has no solution on open sets for $(x,y)\neq0$ (L2(a)), so (S2) has none there. In the slaved vector sector the same holds except for M10(b)'s members away from $s(X)$ (L, M). In Z1 the only known ones are those of O15: (14,0), a flat base, $c_V\kappa^2=12480x'^2$. O10 itself is an off-shell statement about the ℷ-equation at fixed $(t,u)$.
- O9 and the $\kappa_w$ formula are stated for $n=4$. For general $n$ the factor $1/24$ becomes $1/(2n(n-1))$ (same bivector argument), which is not run.

## O.4 Consequences for gaps and claims

- **GU-G15 (ℷ-part).** Now precise, and partly closed.
  - The ℷ-variation needs a weight, because compactly supported variations of ℷ are not compactly supported on $Y$.
  - In every (7,7) case naturality and convergence force ray weights (O6, O7).
  - The equation is empty at $T=0$, and its linearisation is the weighted EH equation.
  - On the invariant sector the ℷ-equation of (S2) is, off-shell, exactly 4D GR with $\Lambda=\Lambda_O$; under (S1) it has no Einstein term (O10, O13). Υ = 0 itself holds in this sector near $s(X)$ only on the solutions of O15.
  - On-shell: among the round-10 solutions (Z1, (14,0), flat base) exactly those with $c_V\kappa^2=12480x'^2$ solve (S2) (O15). In (7,7) no joint solution is known.
  - Open: other sectors, and the draft's choice of weight and system.
- **GU-G14.** The ℷ-route is F's reading (A) at leading order, not (P). In Z1 the two differ ($\Lambda_A\neq\Lambda_P$), so if the 4D equation is to come from varying ℷ, (A) is the candidate. Beyond leading order the ℷ-route gives $\Lambda_O(t,u)$, different from I7(a)'s (P) value.
- **GU-054.** Status unchanged (C). The "$\Lambda g$" of (9.10) appears in the ℷ-equation as $\Lambda_A$ at leading order: weight-independent in Z2 ($9c_H/2c_V$, anti-de Sitter sign); in Z1 up to the free fibre scale.
- **G.2.5.** Its "vol" in the list of ℷ-dependent ingredients is wrong for Z1: $dvol_{g_\aleph}$ does not depend on ℷ (O2(c)). Its "formal only" is replaced by O3–O13.
- **F.1** ("$\int_Ys_Y\,dvol$ is not defined") is made precise by O5(c) and O6. No other earlier statement is contradicted.
