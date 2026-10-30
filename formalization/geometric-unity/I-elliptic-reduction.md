# I. The elliptic twist $W_C\circ\mathrm{Ad}_\Gamma$ in $\Upsilon=0$ and in the 14→4 reduction

Claims: GU-054 (and context for GU-057, GU-072). Gaps: GU-G14, GU-G15, GU-G16. Code: `computation/geometric-unity/I_elliptic_reduction.py` → `.out` (sections I0–I6; exact: integers, rationals, Gaussian rationals; ~100 s on one core). It imports `G_symbols.py`, `H_ellipticity.py` and `shiab_checks.py` without modifying them. Independent check: `computation/geometric-unity/xcheck_round5.py` → `.out` (own Clifford algebra, no repo imports; sections X0–X11). Record: `metrics/geometric-unity/I-elliptic-reduction.yaml`. Builds on D (D2, D3, D5), G (G1, G2, G.2.5), H (H2, H3, H7) and F (F2–F6, F8). Layers: (b) reconstruction, (c) proven, (d) open.

## I.0 Question and setting (layer b)

Round 4 (H7) found that $\circledast_\sharp=W_A+W_C\circ\mathrm{Ad}_\Gamma$ is self-adjoint, elliptic in (14,0) and equal to $W_A$ on Levi-Civita curvature. Rounds 2–3 reduced $\Upsilon=0$ along a section only at $T=0$ (F5, F8.1) or at $T=t\Phi^1$ with the printed Υ (F8.2, constants open). Since $T\neq0$ on solutions near $s(X)$ (G.2.5, Cor. I2.1), the question is whether the twisted part changes the equations.

**Setting.** Everything below is in the unitary gauge of G1: $T$ stands for the $\tau'$-invariant field $S=\varepsilon T\varepsilon^{-1}$, with $T=\varpi-\varepsilon^{-1}d_0\varepsilon$ as in G R-G0. Then $A=\nabla^0+T$ and $F_A=F_0+d_0T+T\wedge T$, where $F_0$ is the Levi-Civita spin curvature of $g_\aleph$, $d_0=d_{\nabla^0}$ and $(T\wedge T)_{ab}=T_aT_b-T_bT_a$.

- For $\circledast\in\mathcal S_{\rm sa}$ (parallel coefficients) the true field equation (D5/G2) is
  $$\Upsilon=\circledast F_0+\circledast(d_0T)+\kappa{*}T+Q(T),\qquad Q(T):=\circledast(T\wedge T)-\tfrac13h(T),$$
  equivalently $p(\alpha,Q(T))=\tfrac13\Psi(\alpha\wedge T\wedge T)+\tfrac13\Psi(T\wedge[T\wedge\alpha])$, the gradient of $\tfrac13\Psi(T\wedge T\wedge T)$. $Q$ is **quadratic** in $T$; $B$ denotes its polarisation, $Q(T)=B(T,T)$.
- Family (H7), in the normalisation of `G_symbols.py`:
  $$\circledast=x\,W_A+y\,W'_A+z\,W_C\mathrm{Ad}_\Gamma+w\,W'_C\mathrm{Ad}_\Gamma,\qquad \circledast_\sharp:\ (x,y,z,w)=(1,0,1,0).$$
  That normalisation has two parts: $\Phi^3=\sum e^{abc}\otimes i\gamma_{abc}$, and the factor $i$ in $[a\wedge b]_A=i\sum e^I\wedge e^J\otimes(a_Ib_J+b_Ja_I)$ (C, R-S; `G_symbols.br`). The real forms below depend on both. With these conventions $W_A(F_0)=-{*}(G_{ab}e^a\otimes\gamma^b)$ and $W'_A(F_0)=-{*}(G_{ab}e^a\otimes\gamma^b\Gamma)$, with $G=G(g_\aleph)$ (D3.4; re-checked in output §I0).
- *Real forms.*
  - (7,7), $\mathfrak u(64,64)$: real $x,y,z,w$, real $T$, real κ. Here $\gamma_a$ and $\gamma_a\Gamma$ lie in $\mathfrak u(64,64)$, and the four operators preserve $\mathfrak u(64,64)$ (output §I0, explicit matrices).
  - (14,0), $\mathfrak u(128)$: $(x,y,z,w)=-i(x',y',z',w')$ (H7) and $T=i(\ldots)$ with real coefficients, κ real.
  - The §I0 checks are spot checks (3 samples for $\mathfrak u(64,64)$). An exact criterion, $p_t\in i^{\epsilon_j+\epsilon_l}\mathbb R$ blockwise, confirms both real forms, the (7,7) one given $\gamma_a\in\mathfrak u(64,64)$ (`xcheck_round5.out` X3b).
- *Parity.* $\mathfrak g=\mathfrak g_{\rm ev}\oplus\mathfrak g_{\rm odd}$ by Clifford parity, so $\mathrm{Ad}_\Gamma=\pm1$ on the two parts; write $T=T_e+T_o$. $d_0$ preserves grades: the spin connection acts by commutators with grade-2 elements.
- *Readings along a section* $s(X)$ are those of F (R-Z3): (P) is the Hodge-dual pull-back $s^*({*}\Upsilon)=0$; (E) is $\Upsilon=0$ at the points of $s(X)$, or on a neighbourhood. Z1 and Z2 are the two horizontal metrics of F (R-Z2).

The computational method: every quantity is computed for each of the 48 basis types of $\mathcal S_{\rm sa}$ (H3) through $\Psi$ (D2: $b(\alpha\wedge\circledast\zeta)=\Psi(\alpha\wedge\zeta)$), hence for all of $\mathcal S_{\rm sa}$ at once. "Output §I$k$" refers to a section of `I_elliptic_reduction.out`; I1–I7 without that prefix are the results below.

- Output §I0 checks that this reproduces the operator code of `G_symbols.py` and `H_ellipticity.py` exactly on random inputs, for all four operators and both signatures.
- It also checks the functional $\alpha\mapsto\Psi(T\wedge[T\wedge\alpha])$ against direct evaluation, and the gradient identity $C(T+\alpha)-C(T-\alpha)=2\,dC_T(\alpha)+2C(\alpha)$.

## I.1 Results (layer c)

**Proposition I1 (what the twisted operators see; exact, output §I1).**

(a) *Grade rules.* A type $(j,l)$ maps input grade $l$ to output grade $j$, and $j$ to $l$. From the H3 coordinates:

| input grade | 0 | 1 | 2 | 12 | 13 | 14 |
|---|---|---|---|---|---|---|
| $W_A$ | 3 | 2 | 1, 5 | 11 | 10 | – |
| $W'_A$ | 11 | 10 | 9, 13 | 3 | 2 | – |
| $W_C\mathrm{Ad}_\Gamma$ | – | 4 | 3 | 9, 13 | 12 | 11 |
| $W'_C\mathrm{Ad}_\Gamma$ | – | 12 | 11 | 1, 5 | 4 | 3 |

The output prints all 15 grades.

(b) *Constituents* (exact, both signatures). Each listed $\mathrm{Spin}(14)$-constituent (in (7,7): $\mathrm{Spin}(7,7)$) of $\Lambda^2\otimes\Lambda^k$ occurs once, so an operator is either zero or injective on it (Schur); one test vector per constituent decides which.
- *Grade 2* ($\mathfrak{so}(14)$-valued 2-forms; $F_0$, the curvature of a metric connection, $T_o\wedge T_o$ for $T_o$ of grade 1):
  - $W_C\mathrm{Ad}_\Gamma$ and $W'_C\mathrm{Ad}_\Gamma$ vanish on Weyl $\oplus\,S^2_0\oplus\mathbb R$, i.e. on every algebraic curvature tensor.
  - They are injective on $\Lambda^4$ (the first-Bianchi defect) and on $\Lambda^2\oplus[2,1,1]_0$ (the pair-antisymmetric part), with outputs in grades 3 and 11.
  - $W_A$ sees $S^2_0,\mathbb R,\Lambda^2$ (output grade 1) and $\Lambda^4$ (output grade 5). It kills Weyl and $[2,1,1]_0$.
- *Grade 1* ($\Lambda^2\otimes\Lambda^1=\Lambda^3\oplus\Lambda^1\oplus[2,1]_0$, where $d^\nabla\theta$ lives):
  - $W_A$ is injective (output grade 2).
  - $W_C\mathrm{Ad}_\Gamma$ sees only the totally antisymmetric part $\Lambda^3$ (output grade 4).
  - $W'_C\mathrm{Ad}_\Gamma$ is injective (output grade 12).
  - $W'_A$ sees only $\Lambda^3$ (output grade 10).
- *Grade 13:* the roles are exchanged ($W'_A$ and $W_C\mathrm{Ad}_\Gamma$ are injective). On grade 12 the twisted pair acts on the curvature-type constituents as $W_A,W'_A$ do on grade 2, with Einstein-type outputs in grades 13 and 1.
- *Grades 0 and 14:* only the untwisted pair sees grade 0, and only the twisted pair sees grade 14.

So in $\circledast_\sharp F_A=W_A F_A+W_C\mathrm{Ad}_\Gamma F_A$, the twisted part never sees $F_0$. On the $\mathfrak{so}(14)$ part of $d_0T+T\wedge T$ it sees exactly the components that violate the algebraic Bianchi identity or pair symmetry.

**Proposition I2 (parity split; every $\circledast\in\mathcal S_{\rm sa}$).**
$$\Upsilon_{\rm ev}=\circledast(d_0T_o)+\kappa{*}T_e+2B(T_e,T_o),\qquad \Upsilon_{\rm odd}=\circledast F_0+\circledast(d_0T_e)+\kappa{*}T_o+Q(T_o)+Q(T_e).$$
That is, $Q(T_e+T_o)=[Q(T_e)+Q(T_o)]_{\rm odd}+[2B(T_e,T_o)]_{\rm ev}$: the cross term is purely even.

*Proof.*
- Every type pairs an even with an odd grade (H2). So $\circledast$ maps $\mathfrak g_{\rm ev}$-valued forms to $\mathfrak g_{\rm odd}$-valued ones and conversely, and $\circledast F_0$ is odd. Equivalently, $\Psi(\alpha\wedge\zeta)$ vanishes unless α and ζ have opposite parity.
- $Q(T_o)$ is odd. In $\Psi(\alpha_e\wedge T_o\wedge T_o)$ an even α meets the even $T_o\wedge T_o$. In $\Psi(T_o\wedge[T_o\wedge\alpha_e])$ the odd $T_o$ meets the odd bracket. Both vanish by parity.
- The same argument shows $Q(T_e)$ is odd.
- $B(T_e,T_o)$ is even. $6\,p(\alpha,B(T_e,T_o))=\Psi(\alpha\wedge(T_e\wedge T_o+T_o\wedge T_e))+\Psi(T_e\wedge[T_o\wedge\alpha])+\Psi(T_o\wedge[T_e\wedge\alpha])$. For odd α the three terms pair α with an odd product, $T_e$ with an even bracket and $T_o$ with an odd bracket, so all vanish.

Exact check (output §I2): $Q_t(T_e)_{\rm ev}=Q_t(T_o)_{\rm ev}=0$ and $2B_t(T_e,T_o)_{\rm odd}=0$ for all 48 types on random exact $T_e,T_o$, and $B_t(T_e,T_o)\neq0$ for all 48 (also `xcheck_round5.out` X5). ∎

*Corollary I2.1.* Let $\kappa\neq0$ and $(x_1,x_{13})\neq0$ (true for every elliptic member, H6.2). Let $\mathcal C:=\{(x,e^\sigma g(x)):x\in X,\ \sigma\in\mathbb R\}$ be the fibre-scaling cone through $s(X)$. It contains $s(X)$ and is unbounded in every fibre. Then, in Z1 or Z2, no solution has $T_o\equiv0$ on an open set $U$ meeting $\mathcal C$.

*Proof.* With $T_o=0$ the even equation gives $\kappa{*}T_e=0$, so $T=0$. Then $\Upsilon=\circledast F_0={*}(x_1G\otimes\gamma+x_{13}G\otimes\gamma\Gamma)=0$ (D3.2, where $W_A$ has $x_1=-1$) forces $G(g_\aleph)=0$ on $U$. But $G(g_\aleph)\neq0$ on $\mathcal C$:
- On $s(X)$ this is F6, which holds for every $g,\lambda,c$. In Z1, F6 rests on F4(c), now proven for every $n$ (K, Cor. K4(b); before 2026-10-30 graded "proven for $n=1,2$; exact verification for $n=3,4,5$; general-$n$ sketch"). The Z2 half of F6 is proven.
- In Z2, $\varphi_\sigma(y)=e^\sigma y$ is an isometry, since $G_\lambda|_y(Dy,Dy)$ is scale invariant and the horizontal block is $c_Hg(x)$.
- In Z1, $\varphi_\sigma(s(X))$ is the section of the base metric $e^\sigma g$. That metric has the same $\aleph$, hence the same $g_\aleph$ (F2.4), so F6 applies to it. ∎

In particular no solution is $\mathfrak g_{\rm ev}$-valued near $\mathcal C$. Since $\mathcal C$ is unbounded in the fibres, $T_o$ cannot have compact support along the fibres either. This excludes a pure contorsion (an $\mathfrak{so}(14)$-valued $T$, i.e. a metric connection with torsion), a pure $\mathfrak u(1)$ field, and any combination of these. The displacement must have an odd part (grades $1,3,\dots,13$).

*Corollary I2.2 (slaving).* Let $\kappa\neq0$. Then $L_{T_o}:=\kappa{*}+2B(\cdot,T_o)$ is invertible at $T_o=0$, where it is $\kappa{*}$, hence for $|T_o|$ small compared with κ. Where it is invertible, $\Upsilon_{\rm ev}=0$ is equivalent to
$$T_e=-L_{T_o}^{-1}\circledast(d_0T_o).$$
So $T_e$ is a first-order expression in $T_o$, and the ansatz $T_e=0$ is consistent iff $\circledast(d_0T_o)=0$.

**Proposition I3 (the twist drops out of the algebraic part on odd fields).** For $\circledast\in\{W_C\mathrm{Ad}_\Gamma,\ W'_C\mathrm{Ad}_\Gamma\}$ and every $T$:
$$\Psi_\circledast(T\wedge T\wedge T)=2\,\Psi_\circledast(T_o\wedge T_e\wedge T_e).$$
Hence:
- $Q_\circledast(T_o)=0$ for every odd $T_o$;
- $Q_\circledast(T_e)=\tfrac23\circledast(T_e\wedge T_e)$;
- yet $\circledast(T_o\wedge T_o)\neq0$ in general.

So in $\Upsilon$ the term $\circledast(T_o\wedge T_o)$ coming from $F_A$ is cancelled exactly by $-\tfrac13h$. Moreover
$$\{\circledast\in\mathcal S_{\rm sa}:\ Q_\circledast(T)=0\ \text{for all odd }T\}=\mathrm{span}\big(O_{0,3},\,E_{0,11},\,W_C\mathrm{Ad}_\Gamma,\,W'_C\mathrm{Ad}_\Gamma\big).$$

*Proof.*
- *Untwisted identity.* For $W_C=[\Phi_{11}\wedge\cdot]_C$ (any $\Phi$), trace cyclicity gives
  $$p(T,W_C(T\wedge T))=\pm\!\int\mathrm{Tr}\big(\Phi_{11}\wedge((T\wedge T)\wedge T-T\wedge(T\wedge T))\big)=0$$
  for every $T$.
- *Twist.* Put $\tau(A,B,C)=p(A,W_C(B\wedge C))$. Since $\mathrm{Ad}_\Gamma$ is an automorphism, $\Psi_{W_C\mathrm{Ad}_\Gamma}(T\wedge T\wedge T)=\tau(T,\hat T,\hat T)$ with $\hat T=T_e-T_o$.
- *Parity.* $W_C$ shifts parity by one, and $p$ pairs equal parities. This leaves $-\tau(T_e,T_e,T_o)-\tau(T_e,T_o,T_e)+\tau(T_o,T_e,T_e)+\tau(T_o,T_o,T_o)$.
- *Conclusion.* The untwisted identity at $T_o$ gives $\tau(T_o,T_o,T_o)=0$. At $T=T_e+T_o$, where the same parity selection applies, it says the first two terms sum to $-\tau(T_o,T_e,T_e)$. Hence the total is $2\tau(T_o,T_e,T_e)$. The statements on $Q$ follow by differentiation, since $Q=\tfrac13\nabla C$ (output §I3 checks them exactly).
- *The span.* "⊇": $O_{0,3}$ and $E_{0,11}$ see only the scalar parts of $T\wedge T$ and $[T\wedge\alpha]$, which vanish. "⊆": the exact null space of $t\mapsto Q_t(T_o)$ on three random odd samples has dimension 4 in both signatures (output §I3). ∎

For the printed equation (9.7), which has no $h$-term, the cancellation fails, and $W_C\mathrm{Ad}_\Gamma(T_o\wedge T_o)$ survives.

**Theorem I4 (Υ = 0 for the elliptic family on odd fields).** Let $T_e\equiv0$. Then $\Upsilon_{x,y,z,w}=0$ is equivalent to the two equations
$$\text{(odd)}\quad -{*}\big(xG\otimes\gamma+yG\otimes\gamma\Gamma\big)+\kappa{*}T+x\,Q_{W_A}(T)+y\,Q_{W'_A}(T)=0,\qquad\text{(even)}\quad \circledast_{x,y,z,w}(d_0T)=0.$$
The algebraic ("Einstein") equation does not contain $(z,w)$. The twist enters only the first-order linear constraint. For $\circledast_\sharp$:
$$-{*}G(g_\aleph)\otimes\gamma+\kappa{*}T+Q_{W_A}(T)=0,\qquad W_A(d_0T)+W_C(\mathrm{Ad}_\Gamma d_0T)=W_A(d_0T)-W_C(d_0T)=0.$$

*Proof.* I2 with $T_e=0$, together with I3 and $W_C\mathrm{Ad}_\Gamma F_0=W'_C\mathrm{Ad}_\Gamma F_0=0$ (output §I0; H6.1). ∎

**Proposition I5 (soldering sector).** Let $T=\theta_{ab}e^a\otimes\gamma^b+\theta'_{ab}e^a\otimes\gamma^b\Gamma$, with $\theta,\theta'$ arbitrary 2-tensors and $T_e=0$.

(a) For symmetric θ and $\theta'=0$ (exact fit, 392 equations, both signatures; complete by the argument below):
- $Q_{W_A}(T)=W_A(T\wedge T)=-{*}(G(R_\theta)\otimes\gamma)$, so $h_{W_A}$ vanishes on this sector;
- $Q_{W'_A}(T)=-\tfrac13{*}(G(R_\theta)\otimes\gamma\Gamma)$;

where
$$(R_\theta)_{abcd}=4(\theta_{ac}\theta_{bd}-\theta_{ad}\theta_{bc}),\qquad G(R_\theta)=4\big[(\mathrm{tr}\,\theta)\theta-\theta^2-\tfrac12((\mathrm{tr}\,\theta)^2-\mathrm{tr}\,\theta^2)g\big].$$
$R_\theta$ is an algebraic curvature tensor, and $T\wedge T$ is its spin curvature. So
$$\Upsilon_{\rm odd}=-{*}\Big(\big[x(G(g_\aleph)+G(R_\theta))-\kappa\theta\big]_{ab}e^a\otimes\gamma^b+y\big[G(g_\aleph)+\tfrac13G(R_\theta)\big]_{ab}e^a\otimes\gamma^b\Gamma\Big).$$
For general θ (with an antisymmetric part) $Q_{W_A}$ also has grade 5 and $Q_{W'_A}$ grade 9. The twisted pair contributes nothing (I3).

*Completeness of the fit.* By I1(a), for grade-1 $T$ the output grades of $Q_{W_A}(T)$ lie in $\{1,3,5\}$, and those of $Q_{W'_A}(T)$ in $\{9,11,13\}$. As an $SO(14)$-module (in (7,7) after complexification), $\mathrm{Sym}^2(\mathrm{Sym}^2V)=[4]\oplus[2,2]\oplus2[2]\oplus2[0]$. It has no constituent in common with $V\otimes\Lambda^3$ or $V\otimes\Lambda^5$ ($\cong V\otimes\Lambda^{11}$, $V\otimes\Lambda^9$), so grades 3, 5, 9 and 11 vanish. Since $V\otimes V=[2]\oplus[0]\oplus\Lambda^2$, the grade-1 (resp. 13) part has exactly $2+2=4$ equivariant parameters, $\theta^2$, $\mathrm{tr}\,\theta\,\theta$, $\mathrm{tr}\,\theta^2g$ and $(\mathrm{tr}\,\theta)^2g$. So the 4-term ansatz is complete, and the unique exact fit proves (a) for every symmetric θ. `xcheck_round5.out` X10 reproduces it as an exact identity.

(b) *Constraint.* Decompose $(d^\nabla\theta,d^\nabla\theta')$ by constituents $\rho\in\{\Lambda^3,\Lambda^1,[2,1]_0\}$ of $\Lambda^2\otimes\Lambda^1$, with coefficients $(a,b)$ on $(v,v\Gamma)$. The even equation becomes a linear system in $(a,b)$ per ρ (output §I4). Its determinants and minors are given up to a constant, which depends on the normalisation of the reference patterns (`xcheck_round5.out` X9 gives them without the factor 4):
- on $\Lambda^1$ and $[2,1]_0$ its determinant is $4i(xz+yw)$ in (14,0) and $4i(xz-yw)$ in (7,7);
- on $\Lambda^3$ (four output grades) it has rank 2 unless all its 2×2 minors vanish. In (14,0) the minors include $4(x^2+y^2)$ and $4(z^2+w^2)$, so the rank is 2 for every nonzero real member, $(x,y,z,w)=-i(x',y',z',w')$. In (7,7) the minors are $4(x^2-y^2)$, $4(z^2-w^2)$, $4i(xz-yw)$ and $4i(wx-yz)$.

The $\Lambda^1$ and $[2,1]_0$ factor is H7's ellipticity factor, as expected. On these constituents $\circledast(d_0T)$ is the principal symbol $\circledast(\xi\wedge\cdot)$ at $\xi=v$, applied to the $\mathrm{Spin}(13)$-invariant odd copies $a\Phi^1+b\Phi^1\Gamma$ (for $\Lambda^1$) and to a hook vector (for $[2,1]_0$), with outputs in the even copies. These are H4's blocks $C_{\Lambda^0}$ and $C_{H_1}$, so the agreement holds by construction; it is not an independent link between the constraint and ellipticity.

Hence:
- if $xz+yw\neq0$ in (14,0), which is H7's ellipticity condition, then $d^\nabla\theta=d^\nabla\theta'=0$;
- for $W_A$ alone the constraint gives only $d^\nabla\theta=0$ and $(d^\nabla\theta')_{\Lambda^3}=0$; the $\Lambda^1\oplus[2,1]_0$ part of $d^\nabla\theta'$ is unconstrained.

(c) For $\circledast_\sharp$ (and every $x\neq0$, $y=0$) in this sector:
$$\Upsilon=0\iff G(g_\aleph)+G(R_\theta)=\tfrac{\kappa}{x}\,\theta\quad\text{and}\quad d^\nabla\theta=0\ \ (\theta\ \text{Codazzi}),$$
with $\theta'=0$. The γΓ-channel vanishes identically. $W_A$ gives the same system.

*Proof.* Propositions I1(b), I2, I3 and Theorem I4, with the exact fits and ranks of output §I4. ∎

**Proposition I6 (invariant sector; closes F8's open constants).** Let $T=t\,\Phi^1+u\,\Phi^1\Gamma$ with $t,u$ functions on $Y$. This is the most general pointwise $\mathrm{Spin}(14)$-invariant (in (7,7): $\mathrm{Spin}(7,7)$-invariant) $T$ (D1). Then exactly:
$$\Upsilon=-{*}\Big([xG(g_\aleph)-c_1\,g]_{ab}e^a\otimes\gamma^b+[yG(g_\aleph)-c_{13}\,g]_{ab}e^a\otimes\gamma^b\Gamma\Big)+\circledast(dt\wedge\Phi^1+du\wedge\Phi^1\Gamma),$$
$$(14,0):\ \ c_1=\kappa t+x(312t^2+104u^2)+208\,y\,tu,\qquad c_{13}=\kappa u+208\,x\,tu+y(104t^2+312u^2);$$
$$(7,7):\ \ c_1=\kappa t+x(312t^2-104u^2)-208\,y\,tu,\qquad c_{13}=\kappa u+208\,x\,tu+y(104t^2-312u^2).$$
- $(z,w)$ do not occur.
- The derivative terms occupy grades 2 and 12. In a fixed pattern they are $[-2x,\mp2y]\cdot(dt,du)$ and $[2iw,-2iz]\cdot(dt,du)$, with determinant $4i(xz\pm yw)$ up to a constant (upper signs (14,0), lower (7,7)). Since $dt\wedge\Phi^1$ is the $\Lambda^1$ constituent with $v=dt$, this is the $\Lambda^1$ block of I5(b) again: I5(b) and I6 are one 2×2 block, not two confirmations.
- So for $\circledast_\sharp$ the constraint forces $dt=du=0$. For $W_A$ alone it forces only $dt=0$; there $u$ is fixed by $c_{13}=u(\kappa+208t)=0$, or, on the branch $t=-\kappa/208$, by a contracted Bianchi identity. For solutions of the full equation on an open set in $Y$, the 14D identity makes $c_1$, hence $u^2$ and $u$, constant there. **The two solution sets on such open sets coincide.** Under (P), the 4D identity for $G(g)+\Lambda_4(x)g=0$ makes $\Lambda_4=\Lambda_P-c_Hc_1/x$ constant only along $s(X)$. It fixes $u|_{s(X)}$, but leaves normal derivatives of $u$ free. The twist still requires those derivatives to vanish, so **the projected solution sets differ**, although their induced 4D algebraic equations and allowed cosmological constants agree (counterexample I6.1).
- For $T=t\Phi^1$ one has $h_{W_A}(\Phi^1)=0$, so true and printed Υ agree. F8.2's $\Lambda_Y(t)=\alpha t^2+\beta\kappa t$ (with $G+\Lambda_Yg=0$) has $\alpha=-312$, $\beta=-1$ in this normalisation.
  - This also holds for F8's R-S operator with the printed Υ. Here $F_A=F_0+t^2\Phi^1\wedge\Phi^1$ is the spin curvature of the algebraic curvature tensor $R_0+R_{tg}$, and R-S $=W_A$ on such tensors (D3.4).
  - It also holds under reading (P). The pulled-back derivative terms still force $dt=0$, and $du=0$ for the twist: the component with form index $b$ involves every $\partial_ct$ ($c\ne b$) through distinct monomials, so two tangential indices suffice (output §I5).

*Proof.* Exact computation, output §I5: every T-term was checked componentwise to be of the invariant form. The derivative pattern is the same for the directions $e^0,e^5,e^9$, and for all directions by equivariance. ∎

**Counterexample I6.1 (normal derivatives under (P); added 2026-10-30).** Take a flat Lorentzian 4D base, Z2, $\lambda=\tfrac12$, $c_H=1$, $c_V=\tfrac3{52}$, $\kappa=208$, $t=-1$ and $u|_{s(X)}=\tfrac12$. The observerse signature is $(7,7)$. Extend $u$ smoothly off the section with a nonzero normal derivative and no tangential derivative; such an extension exists locally. Then
$$c_1=208(-1)+312-104/4=78=\Lambda_P,\qquad c_{13}=u(208-208)=0.$$
The projected algebraic equation for $W_A$ is $G(g)+(\Lambda_P-c_1)g=0$, which holds for the flat base. Its derivative term also vanishes: $dt=0$ and $W_A(du\wedge\Phi^1\Gamma)=0$. For $\circledast_\sharp$, the same algebraic equation holds but the projected grade-12 constraint is nonzero. Hodge dualising and pulling back keeps tangential **form** indices while retaining the Clifford indices, so a normal derivative is still detected. Thus the same field can satisfy (P) for $W_A$ and fail it for $\circledast_\sharp$.

`computation/geometric-unity/projection_regression.py` → `.out` checks the rational coefficients and this normal-derivative term directly from the Clifford definitions, with integer arithmetic. It does not assert a solution of $\Upsilon=0$ on a neighbourhood in $Y$. The former equivalence claim used the 4D Bianchi identity to control normal derivatives, which it cannot do.

**Theorem I7 (reduction along $s(X)$, $n=4$).** Use the F closed forms (output §I6 re-derives $\Lambda_P$):
- $\Lambda_P(\mathrm{Z2})=\tfrac{9c_H}{2c_V}$ and $\Lambda_P(\mathrm{Z1})=\tfrac{3c_H}{2c_V}\tfrac{4-15\lambda}{1-4\lambda}$;
- in the F6 case (Z1, $\lambda=\tfrac15$, $\mathrm{Ric}(g)=0$), $G^Y=-\Lambda_Yg_\aleph$ on $s(X)$ with $\Lambda_Y=-\tfrac{15}{2c_V}$.

(a) *Invariant sector, reading (P).* By I6, $t$ and $u$ are constant along $s(X)$ (for $\circledast_\sharp$ by the constraint; for $W_A$ with the 4D Bianchi identity, I6). The pull-back is 4D vacuum Einstein, $G(g)+\Lambda_4g=0$, with
$$(7,7):\ \Lambda_4=\Lambda_P-\tfrac{c_H}{x}c_1(t,u),\qquad (14,0):\ \Lambda_4=\Lambda_P+\tfrac{c_H}{x'}\big(\kappa s+x'(312s^2+104r^2)+208y'sr\big),$$
on the curve $yc_1=xc_{13}$ (resp. $y'e_1=x'e_{13}$, with $e_1,e_{13}$ the brackets of output §I6). For $\circledast_\sharp$ and $W_A$ the common algebraic curve is $u(\kappa+208t)=0$, resp. $r(\kappa+208s)=0$. Their projected constraints differ on normal derivatives (I6.1).
- For $\circledast_\sharp$ and $c_H>0$, $\Lambda_4$ takes **every real value** in (7,7). In (14,0) it takes exactly the values in $[\Lambda_P-c_H\kappa^2/1248,\infty)$. For $x'W_A+z'W_C\mathrm{Ad}_\Gamma$ ($x'>0$) the bound is $\Lambda_P-c_H\kappa^2/(1248x'^2)$ (`xcheck_round5.out` X11).
- This holds in Z1 and Z2 alike, since only $\Lambda_P$ differs.
- So **reading (P) is underdetermined**. Besides $T=0$, with $\Lambda_4=\Lambda_P$ (F5, F8.1), it admits invariant-sector $T\neq0$ with any of these $\Lambda_4$, and the fibre geometry no longer fixes Λ under (P). Solutions of $\Upsilon=0$ near $s(X)$ are a different class, see (b).
- The cosmological term of (9.10) is reproduced, with $\Lambda_4$ undetermined. Only $\kappa t$ comes from $T_\omega$ ($\kappa{*}T$). $\Lambda_P$ and $312t^2$ come from $S_\omega=\circledast F_A$, and the $u$-dependence of the (7,7) range from $-\tfrac13h$, which the draft does not have. With the printed (9.7) the same sector gives $c_1^{\rm pr}=\kappa t+312x(t^2\pm u^2)$ and $c_{13}^{\rm pr}=\kappa u+312y(t^2\pm u^2)$ (upper signs (14,0); `xcheck_round5.out` X8p). For $\circledast_\sharp$ and $\kappa\neq0$ the curve then forces $u=0$, and (7,7) gives only $\Lambda_4\le\Lambda_P+c_H\kappa^2/1248$ (X11). "Every real value" is a property of the true Υ (D5), not of the draft's equation.

(b) *Invariant sector, reading (E) on $s(X)$.* It requires $g_\aleph$ Einstein on $s(X)$. By F6 this means Z1, $\lambda=\tfrac15$ and $\mathrm{Ric}(g)=0$, in every member with $(x,y)\neq0$ and for every $(t,u)$. For $x=y=0$, $\circledast F_0=0$, $c_1=\kappa t$, $c_{13}=\kappa u$, and $T=0$ solves $\Upsilon=0$ identically (`xcheck_round5.out` X11); every elliptic member has $(x,y)\neq0$ (H6.2).
- In (7,7) there is therefore **no solution**: $\lambda=\tfrac15$ gives signature (8,6) or (10,4).
- In (14,0), with base (4,0) and $c_V>0$, $\circledast_\sharp$ (and $W_A$) needs $312s^2+\kappa s+104r^2+\tfrac{15}{2c_V}=0$ together with $r(\kappa+208s)=0$. A real solution exists **iff $c_V\kappa^2\ge9360$**, with $r=0$; the branch $s=-\kappa/208$ has none.
- These are the configurations Cor. I2.1 concerns: a solution of $\Upsilon=0$ on a neighbourhood of a point of $s(X)$ satisfies (E) there, and then I2.1 gives $T_o\neq0$. Their (P) value is fixed: (E) at $s(X)$ implies (P), and $\mathrm{Ric}(g)=0$ gives $\Lambda_4=0$. For $\circledast_\sharp$, $e_1=-\tfrac{15}{2c_V}$ and $\Lambda_4=\Lambda_P(\mathrm{Z1},\tfrac15)-\tfrac{15c_H}{2c_V}=\tfrac{15c_H}{2c_V}-\tfrac{15c_H}{2c_V}=0$ (X11). So the free $\Lambda_4$ of (a) and the $T_o\neq0$ of I2.1 apply to disjoint classes.
- Off $s(X)$ the Einstein condition is not available. F's E2(v) shows it failing for one Weyl base, in base signature (1,3), and Conjecture F9 says it fails for every curved ℵ. So (E) on a neighbourhood is not established in this sector.

(c) *Symmetric displacement (I5(c)), reading (P).* With θ block-diagonal on $s(X)$, write $h:=s^*\theta$, $\tau_1=\mathrm{tr}_{g_V}\theta_V$, $\tau_2=\mathrm{tr}_{g_V}\theta_V^2$. Then
$$G(g)+\Lambda_Pg+4\Big[\big(\tfrac{\mathrm{tr}_gh}{c_H}+\tau_1\big)h-\tfrac{h^2}{c_H}-\tfrac{c_H}2\Big(\big(\tfrac{\mathrm{tr}_gh}{c_H}+\tau_1\big)^2-\tfrac{\mathrm{tr}_gh^2}{c_H^2}-\tau_2\Big)g\Big]-\tfrac{\kappa}{x}h=0 .$$
This is the algebraic part of (P): 4D Einstein + $\Lambda_P$ plus a **displacement stress**, quadratic and linear ($-\tfrac\kappa xh$) in a new symmetric field $h$ and two fibre scalars. (P) does not determine $h,\tau_1,\tau_2$.
- The (P)-projection of the constraint $\circledast(d_0T)=0$ was not computed. It keeps $4\cdot91=364$ of the $14\cdot91$ grade-2 output components. For symmetric θ, $d^\nabla\theta$ lies in $\Lambda^1\oplus[2,1]_0$, of dimension 910, so the projected constraint cannot force $d^\nabla\theta=0$. Which 4D conditions on $(g,h,\tau_1,\tau_2)$ it imposes is open (I.3).
- For solutions of $\Upsilon=0$ near $s(X)$, (E) at $s(X)$ implies (P). The equation above then holds automatically once $h$ is the solution of (d), and carries nothing beyond (E).

(d) *Symmetric displacement, reading (E).* Let $F(\theta):=\kappa\theta-xG(R_\theta)$, so $dF_0=\kappa\,\mathrm{Id}$.
- At each point of $s(X)$, the algebraic equation $F(\theta)=xG^Y$ has a unique small solution whenever $|G^Y|$ is small compared with $\kappa^2$ (inverse function theorem). On $s(X)$, $G^Y$ contains the fibre part $G^Y_{VV}=O(1/c_V)$ (F3/F4), so with $x=O(1)$ this needs $c_V\kappa^2\gg1$, together with $|G(g)+\Lambda_Pg|\ll\kappa^2$; it is the analogue of the invariant-sector condition $c_V\kappa^2\ge9360$. That solution is block-diagonal: $F$ is $O(T_yY,g_\aleph)$-equivariant, $G^Y_{HV}=0$ (F3/F4), and uniqueness gives invariance under $\sigma=\mathrm{diag}(1_H,-1_V)\in O(T_yY,g_\aleph)$, since $H\perp V$ (F2). **Pointwise, for $|G^Y|\ll\kappa^2$, (E) imposes no condition on $g$.**
- The constraint $d^\nabla\theta=0$ on a neighbourhood has HHH components equal to $d^{\nabla_g}(s^*\theta)$. Write $\nabla_YY'=\nabla^s_YY'+\mathrm{II}(Y,Y')$ for $Y,Y'$ tangent to $s(X)$; $\theta(\mathrm{II},\cdot)$ vanishes because $\theta_{VH}=0$. The induced connection is Levi-Civita of $c_Hg$.
- So the HHH constraint requires $s^*\theta$, the implicit function above, to be Codazzi. With $s^*\theta=\tfrac x\kappa(G(g)+\Lambda_Pg)+O(\kappa^{-3})$, this is **$G(g)$ Codazzi to leading order in $1/\kappa$**; that truncation is not itself a necessary condition. Contracting gives $dR=0$, so $G(g)$ Codazzi is equivalent to $\mathrm{Ric}(g)$ being Codazzi, i.e. to harmonic curvature $\delta\mathrm{Rm}=0$. Every Einstein metric satisfies it, but it is strictly weaker: products of Einstein metrics with different constants have parallel, non-Einstein Ricci tensor.
- The HHV, HVV and VVV components involve normal derivatives of θ and were not computed.

(e) In the algebraic equations of (a)–(d) the twisted coefficients $z,w$ appear nowhere. Their effect in these sectors is the constraint on the Γ-partner: $du=0$, $d^\nabla\theta'=0$. Under (P), $du=0$ includes normal derivatives and changes the set of allowed fields off the section (I6.1), while leaving the induced 4D algebraic equation and the range of $\Lambda_4$ unchanged. Equality of full-neighbourhood solution sets in the invariant sector follows instead from the 14D Bianchi identity (I6).

*Proof.* (a), (b): Proposition I6 together with F3–F6; the extrema are completed-square identities, and the discriminant is computed exactly in output §I6. (c): substitute I5(a), using $s^*g_\aleph=c_Hg$ and $s^*G^Y=G(g)+\Lambda_Pg$ (F5). (d): as stated. (e): I4–I6. ∎

**Answers to the three questions.**
1. The twisted part never acts on $F_0$. On the $\mathfrak{so}(14)$-valued part of $F_A$ it acts exactly on the non-curvature constituents ($\Lambda^4$, $\Lambda^2$, $[2,1,1]_0$); on the other grades it acts as listed in I1. In $\Upsilon$ its $T\wedge T$ contribution is cancelled by $-h/3$ on odd fields. What remains of it is $W_C\mathrm{Ad}_\Gamma(d_0T)$ and $\tfrac13\nabla[2\Psi_{W_C\mathrm{Ad}_\Gamma}(T_o\wedge T_e\wedge T_e)]$ (I3).
2. Relative to F3/F4 the 4D equations do change once $T\neq0$. Under (P), which is underdetermined, Λ becomes a free constant (invariant sector), or a displacement stress appears (symmetric sector). Solutions of $\Upsilon=0$ near $s(X)$ satisfy (E) there. Under (E), with $(x,y)\neq0$, (7,7) has no solutions in the invariant sector and (14,0) has $\Lambda_4=0$; in the symmetric sector, for $|G^Y|\ll\kappa^2$, no pointwise 4D algebraic equation arises at all. **Relative to $W_A$, $\circledast_\sharp$ preserves these 4D algebraic equations and cosmological ranges; under (P) its constraint excludes additional normal derivatives** (I7(e), I6.1).
3. See I.3.

## I.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| output §I0 (conventions) | exact comparison with the existing operator code; functional and gradient identities; curvature normalisation. Spot checks (25 random inputs; 3 $\mathfrak u(64,64)$ samples): the proofs rest on the H3 coordinates, which `xcheck_round5.out` X3 rebuilds from the definitions | random exact inputs; both signatures; all 4 operators and all 48 types |
| I1 | grade rules from the H3 coordinates; one exact test vector per constituent plus Schur | input grades 0, 1, 2, 12, 13, 14 (constituents); all grades (rules) |
| I2, I3 | proof; exact checks on random exact fields; exact null space (29516 × 48 and 26864 × 48) | all of $\mathcal S_{\rm sa}$; both signatures |
| I4 | proof from I2, I3 | all odd $T$ |
| I5 | exact fits ($a_1,\dots,a_4$, 392 equations), complete by the $SO(14)$ multiplicities of $\mathrm{Sym}^2(\mathrm{Sym}^2V)$; exact per-constituent systems and minors | soldering sector; both signatures |
| I6 | exact; invariant form checked componentwise; derivative pattern checked in 3 directions; I6.1 checks the projected normal-derivative counterexample directly from the Clifford definitions | invariant sector; full-neighbourhood equivalence and projected inequivalence distinguished |
| I7 | sympy identities from F3–F6 closed forms; proof for (c), (d) | $n=4$; Z1, Z2; (P), (E) |

## I.3 Limits and ill-posedness (layers b, d)

- **Not covered.**
  - Odd $T$ outside the soldering sector (grades 3–11): Theorem I4 still applies to the algebraic part, but the constraint was not analysed.
  - The soldering sector with $\theta'\neq0$: $Q$ was fitted in closed form only for symmetric θ with $\theta'=0$; for $\theta'\neq0$ only its grades were computed. Whether $W_A$ and $\circledast_\sharp$ have the same full-neighbourhood solution set there is not determined: $W_A$ leaves the $\Lambda^1\oplus[2,1]_0$ part of $d^\nabla\theta'$ free, $\circledast_\sharp$ does not. In the invariant sector the full-neighbourhood solution sets coincide (I6); the projected solution sets already differ because of normal derivatives (I6.1).
  - Solutions with slaved $T_e\neq0$ (I2.2). There the odd equation becomes second order in $T_o$, and $z,w$ do enter it, through $\circledast_\sharp(d_0T_e)$ and $\tfrac23W_C\mathrm{Ad}_\Gamma(T_e\wedge T_e)$ (I3), and through $L_{T_o}$ in the slaving map. $B_{W_C\mathrm{Ad}_\Gamma}(\cdot,T_o)$ is even (I2), so it enters only the even equation. Not reduced here.
  - In the symmetric sector, the (P)-projection of the constraint $\circledast(d_0T)=0$ (I7(c)). It cannot force $d^\nabla\theta=0$ (364 against 910 components); which 4D conditions on $(g,h,\tau_1,\tau_2)$ it imposes, for example whether it contains $d^{\nabla_g}h=0$, is open.
  - The reduced-action reading (A) with $T\neq0$.
  - A chimeric-invariant ansatz using $\Gamma_X,\Gamma_V$.
  - The signatures (8,6) and (10,4) of the F6 case with base (1,3): the identities were computed in (14,0) and (7,7) only, and the real-root analysis of (E) was done only for (14,0)/$\mathfrak u(128)$.
  - The ℷ-variation (G.2.5) stays open.
- **Precise negative statements.**
  1. Reading (P) does not determine $T$ along $s(X)$. In the invariant sector the 4D cosmological constant is then a free parameter (I7(a)).
  2. Reading (E) at points of $s(X)$ is either impossible ((7,7), invariant sector, $(x,y)\neq0$), fixes $\Lambda_4=0$ ((14,0), invariant sector), or is vacuous pointwise (symmetric sector, for $|G^Y|\ll\kappa^2$). Its 4D content comes only from the constraint equations, whose remaining components involve normal derivatives of $T$ at $s(X)$; data on $s(X)$ do not fix these.
- **Extra structure needed** to make the 14→4 reduction well-posed for $T\neq0$:
  - a rule fixing $T$ off $s(X)$, either a declared ansatz or $\Upsilon=0$ on a neighbourhood together with boundary or decay conditions along the non-compact fibre (the fibre weight of G.2.5);
  - and a choice among (P), (A), (E).
  The draft supplies neither.
- **Scope of G.2.5 (logged in `reports/CORRECTIONS.md`).** G.2.5 says $T$ "cannot vanish on any open set". The argument it cites (F6) gives $G(g_\aleph)\ne0$ only on $s(X)$, and nothing in the repo shows $G(g_\aleph)\neq0$ on all of $Y$.
  - Corollary I2.1 proves the statement for open sets meeting the scaling cone $\mathcal C\supset s(X)$, with $T_o$ in place of $T$. In Z1 it uses F4(c), proven for every $n$ in K.
  - That scope still excludes compact support along the fibres, so G.2.5's conclusion is unaffected.
  - *Remark (2026-10-30, hypothesis: $X$ connected and $g$ real-analytic).* Then G.2.5's "any open set" holds as printed, with $T_o$ in place of $T$. In Z1 and in Z2, $g_\aleph$ is real-analytic on $Y$, since $y^{-1}$, $D$ and the Christoffel symbols of $g$ are. $Y=\mathrm{Met}_{(p,q)}(X)$ is connected, because $O(p,q)$ meets both components of $GL(4)$. So $G(g_\aleph)$ cannot vanish on an open set without vanishing on all of $Y$, including $s(X)$, which contradicts F6. The proof of I2.1 then gives $T_o\not\equiv0$ on every open set (under I2.1's hypotheses; in Z1 via F4(c), proven in K). For smooth non-analytic $g$ the scope remains that of I2.1.

## I.4 Consequences for gaps and claims

- **GU-G14:** refined.
  - (P) is underdetermined. Besides $T=0$ ($\Lambda_4=\Lambda_P$) it admits invariant-sector $T\neq0$ with any $\Lambda_4$ (every real value in (7,7)), and in the symmetric sector a displacement stress (algebraic part; the projected constraint is open).
  - Solutions of $\Upsilon=0$ near $s(X)$, for which I2.1 gives $T_o\neq0$, satisfy (E) on $s(X)$. In the invariant sector none exists in (7,7), and in (14,0) they have $\Lambda_4=0$ (both for $(x,y)\neq0$). In the symmetric sector (E) yields no pointwise 4D equation for $|G^Y|\ll\kappa^2$, and its HHH constraint requires $s^*\theta$ Codazzi.
  - The elliptic twist preserves these algebraic equations and cosmological ranges. Its first-order constraint on the Γ-partner changes the projected solution set through normal derivatives (I6.1).
- **GU-G16:** H7's ellipticity factor $xz+yw$ is, up to a constant, the determinant of the reduced constraint system on the $\Lambda^1$ and $[2,1]_0$ parts of $(d^\nabla\theta,d^\nabla\theta')$ (I5(b)); the $(dt,du)$ system of I6 is the same $\Lambda^1$ block. This holds by construction: there the constraint is the principal symbol of $\circledast\circ d$, and the blocks are H4's $C_{\Lambda^0}$ and $C_{H_1}$.
- **GU-G15 / G.2.5:** for $\kappa\neq0$ and $(x_1,x_{13})\neq0$, $T_o\not\equiv0$ on every open set meeting $\mathcal C$ (on every open set if $g$ is real-analytic, I.3), and $T_e$ is slaved to $T_o$ where $L_{T_o}$ is invertible (I2).
- **GU-054:** status C unchanged; the constants are now explicit for $T\neq0$ (amended in `claims/geometric-unity.md`).
- **F8:** its open constants are now determined (I6), and its real-root condition is made explicit in (14,0) (I7(b)).
