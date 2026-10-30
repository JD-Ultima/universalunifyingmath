# J. The slaved even sector: does the elliptic twist change the reduced equations once $T_e\neq0$?

Claims: GU-054 (context for GU-057, GU-072). Gaps: GU-G14, GU-G15, GU-G16. Code: `computation/geometric-unity/J_slaved_sector.py` → `.out` (sections J0–J5 and J4b; exact: integers, rationals, Gaussian rationals and integers, sympy; 567 s (~10 min) on one core after the review additions; a second run gave identical output apart from the timing annotations, 560 s). It imports `I_elliptic_reduction.py` (and through it `G_symbols.py`, `H_ellipticity.py`, `shiab_checks.py`) without modifying it. Record: `metrics/geometric-unity/J-slaved-sector.yaml`. Builds on I (I2, I2.2, I3, I6, I7), H (H3, H4, H7), D (D3.4, D5), G (G1) and F (F2–F6). Layers: (b) reconstruction, (c) proven, (d) open.

## J.0 Question and setting (layer b)

Round 5 settled the reduction only for $T_e=0$. There the twisted coefficients $z,w$ of $\circledast_{x,y,z,w}$ drop out of every reduced equation (I7(e)). I.3 left open the solutions with $T_e\neq0$, which the even equation forces as soon as $\circledast(d_0T_o)\neq0$ (Cor. I2.2). This note settles that case on the largest sector in which the even equation closes exactly.

**Setting.** As in I.0: the unitary gauge of G1, $A=\nabla^0+T$, $F_A=F_0+d_0T+T\wedge T$, the true field equation (D5) $\Upsilon=\circledast F_0+\circledast(d_0T)+\kappa{*}T+Q(T)$, and the family (H7, code normalisation of `G_symbols.py`)
$$\circledast=x\,W_A+y\,W'_A+z\,W_C\mathrm{Ad}_\Gamma+w\,W'_C\mathrm{Ad}_\Gamma,\qquad \circledast_\sharp=(1,0,1,0).$$
- Real forms are those of I.0: (7,7), $\mathfrak u(64,64)$, real $x,y,z,w,T,\kappa$; (14,0), $\mathfrak u(128)$, $(x,y,z,w)=-i(x',y',z',w')$ and $T=i(\ldots)$.
- We write $\Upsilon={*}Y$, i.e. $Y_a=\eta_{aa}D_a$ with $\Upsilon(\alpha)=\sum_a\mathrm{Tr}(D_a\alpha_a)$. Then the κ-term is $Y=\kappa T$, and the I6 part is $Y=-[xG-c_1g]_{ab}e^a\otimes\gamma^b-[yG-c_{13}g]_{ab}e^a\otimes\gamma^b\Gamma$.
- By I2, $Y_{\rm ev}=\circledast(d_0T_o)+\kappa T_e+2B(T_e,T_o)$ and $Y_{\rm odd}=\circledast F_0+\circledast(d_0T_e)+\kappa T_o+Q(T_o)+Q(T_e)$ (dualised).
- Readings (P), (E) and metrics Z1, Z2 are those of F (R-Z2, R-Z3). Base dimension $n=4$; $V=\mathrm{Sym}^2T^*_xX$ along $s(X)$, with $\tau(h)=\mathrm{tr}(g^{-1}h)$ and $G_0(h,k)=\mathrm{tr}(g^{-1}hg^{-1}k)$; $\mathfrak h:=g_\aleph|_H=c_Hg$ on $s(X)$.

**The sector $\mathcal V$ (declared ansatz).**
- $T_o=t\,\Phi^1+u\,\Phi^1\Gamma$ with $t,u$ functions on $Y$: the invariant sector of I6, now with $dt,du\neq0$ allowed.
- $T_e=P_0(A)+P_2(C)+P_{12}(D)+P_{14}(B)$ with four covector fields $A,C,D,B$, where
$$P_0(v)_a=v_a\,1,\quad P_{14}(v)_a=v_a\,\Gamma,\quad P_2(v)_a=\sum_{b\neq a}v_b\,\gamma_a\gamma^b,\quad P_{12}(v)=P_2(v)\Gamma .$$
  These are the four Spin(14)-equivariant copies of the vector representation in $\Lambda^1\otimes\mathfrak g_{\rm ev}$ (grades 0, 2, 12, 14; $V\otimes\Lambda^k\supset V$ only for $k\in\{0,2,12,14\}$). They have disjoint monomial supports, so they are linearly independent for $v\neq0$. The representation theory is motivation only: closure is proven by computation (J1).
- Why this sector: J1(a) shows that $\circledast(d_0T_o)$ lies in it for every invariant $T_o$, and J1(b) that $L_{T_o}$ preserves it. So it contains the slaved $T_e$ of every invariant $T_o$ wherever $L_{T_o}$ is invertible; for $\circledast_\sharp$ with $u\neq0$ all four components occur (J1(d)). Nonconstant $t,u$ force $T_e\neq0$ whenever the source matrix $S$ of J1 has full rank, i.e. $xz\pm yw\neq0$ (for $z=w=0$, i.e. $\mathrm{span}(W_A,W'_A)$: whenever $x\,dt\pm y\,du\neq0$; for $W_A$: whenever $dt\neq0$).

**Rules fixing $(t,u)$ near $s(X)$ (declared where used).**
- (B) *basic*: $t=\pi^*\tau$, $u=\pi^*\upsilon$.
- (J) *isotropy-invariant 2-jet along $s(X)$*: at each $p\in s(X)$ the 2-jet of $(t,u)$ is invariant under the isotropy $O(g_x)$ acting on $T_pY=H\oplus V$. $O(g_x)$ contains $-1$, which acts as $-1$ on $H$ and $+1$ on $V$; this kills $H^*$ and $H\otimes V$. The remaining statements are $\mathfrak{so}(g_x)$-representation theory: $V^*$ contains only τ, $S^2H$ only $g$, and $S^2V$ only $G_0$ and $\tau\otimes\tau$. (Corrected after review 2026-10-30, R6-6; see `reports/CORRECTIONS.md`.) The tangential derivatives along $s(X)$ are part of the 2-jet. So (for $X$ connected): $t,u$ constant on $s(X)$; $d\phi=a_\phi\tau$ with constants $a_\phi$; $\mathrm{Hess}\,\phi=\beta_\phi\,\mathfrak h+p_\phi G_0+q_\phi\,\tau\otimes\tau$ ($\phi\in\{t,u\}$), with $p_\phi,q_\phi$ free (second normal derivatives) and $\beta_\phi$ forced (J4).
- (R) *radial*: $t=f(s)$, $u=h(s)$ near $s(X)$, with $s:=\tfrac14\log|\det(g(x)^{-1}y)|$, the fibre-scaling coordinate relative to the section. Its 2-jet on $s(X)$ is of type (J) with $p_\phi=0$.

**Hypotheses used throughout:** $\kappa\neq0$, and $\det M(t,u)\neq0$ for the slaving block $M$ of J1 (true for $|T_o|$ small compared with κ, where $M\approx\kappa\,1$).

## J.1 Results (layer c)

**Proposition J1 (closure of the even equation; explicit slaving).** On $\mathcal V$, for both signatures and every $(x,y,z,w)$:

(a) $Y_{\rm ev}[\circledast(d_0T_o)]=\sum_kP_k\big(S_{k,t}\,dt+S_{k,u}\,du\big)$, with rows $k=(0,2,12,14)$
$$(14,0):\ S=\begin{pmatrix}0&0\\-24x&-24y\\24iw&-24iz\\0&0\end{pmatrix},\qquad (7,7):\ S=\begin{pmatrix}0&0\\-24x&24y\\24iw&-24iz\\0&0\end{pmatrix}.$$

(b) $Y_{\rm ev}[\kappa T_e+2B(T_e,T_o)]=\sum_jP_j\big((M\vec A)_j\big)$ with $\vec A=(A,C,D,B)$ and $M=\kappa1+tM_1+uM_{13}$:
$$M=\begin{pmatrix}\kappa&0&208\epsilon\delta&0\\0&\kappa+224\sigma&384i\epsilon\theta&16\epsilon\delta\\16\delta&-384i\theta&\kappa+544\sigma&0\\0&208\delta&0&\kappa\end{pmatrix},\qquad \sigma=tx+\epsilon uy,\ \theta=tw-uz,\ \delta=ty-ux,$$
with $\epsilon=+1$ in (14,0) and $\epsilon=-1$ in (7,7). Its determinant is
$$\det M=(\rho+224\kappa\sigma)(\rho+544\kappa\sigma)-\epsilon\,(384\kappa\theta)^2,\qquad \rho=\kappa^2-3328\epsilon\,\delta^2 ,$$
so $\det M=\kappa^4+O(\kappa^3|T_o|)$. The grade-0 and grade-14 components $A,B$ couple to $C,D$ only through $\delta=ty-ux$, and $z,w$ enter $M$ only through θ.

(c) Hence on $\mathcal V$ the even equation $Y_{\rm ev}=0$ is exactly the linear system $M\vec A=-S\,(dt,du)^\top$. Where $\det M\neq0$ its solution is
$$\vec A=\Sigma\,(dt,du)^\top,\qquad \Sigma=-M^{-1}S .$$
Where $L_{T_o}$ is invertible on all of $\Omega^1(\mathfrak g_{\rm ev})$ (for instance $|T_o|$ small, I2.2), this is *the* slaved $T_e$ of I2.2.

(d) To leading order in $1/\kappa$ (14,0): $C=\tfrac{24}{\kappa}(x\,dt+y\,du)$, $D=\tfrac{24i}{\kappa}(-w\,dt+z\,du)$, and $A,B=O(\kappa^{-2})$. Exactly, for $\circledast_\sharp$ all four components are nonzero for $u\neq0$; for $W_A$, $\Sigma$ has only the entries $\Sigma_{C,t}=24\kappa/(\kappa^2+224\kappa t-3328\epsilon u^2)$ and $\Sigma_{B,t}=4992u/(\kappa^2+224\kappa t-3328\epsilon u^2)$ (output §J1; $\epsilon=+1$ in (14,0), $-1$ in (7,7)). **$z$ and $w$ enter $T_e$**, through the $D$-row of $S$ and through θ in $M$.

*Proof.* (a), (b): the maps $v\mapsto Y_{\rm ev}[\circledast(v\wedge\Phi^1)]$, $v\mapsto Y_{\rm ev}[\circledast(v\wedge\Phi^1\Gamma)]$ and $v\mapsto Y_{\rm ev}[2B(P_k(v),\Phi^1)]$, $Y_{\rm ev}[2B(P_k(v),\Phi^1\Gamma)]$ are linear in $v$. Output §J1 evaluates each on all 14 basis covectors $e^c$, for each of the four operators, and finds it exactly equal to $\sum_js_jP_j(e^c)$ with the same $s_j$ for every $c$. By linearity this holds for every $v$. Since $\Phi^1$ and Γ are $\nabla^0$-parallel, $d_0T_o=dt\wedge\Phi^1+du\wedge\Phi^1\Gamma$, and $B$ is bilinear, so (a) and (b) follow pointwise. (c): the $P_j(v)$ are linearly independent for $v\ne0$, so $Y_{\rm ev}=0$ on $\mathcal V$ is equivalent to $S(dt,du)^\top+M\vec A=0$. (d): Neumann series of $M^{-1}$. ∎

**Theorem J2 (the odd equation on $\mathcal V$).** For arbitrary covector fields $A,C,D,B$ (slaved or not) and $T_o=t\Phi^1+u\Phi^1\Gamma$,
$$Y_{\rm odd}=\sum_bE^{(1)}_{ab}\gamma^b+\sum_bE^{(13)}_{ab}\gamma^b\Gamma+\sum_{b,c\neq a}\omega^{(3)}_{bc}\gamma_a\gamma^b\gamma^c+\sum_{b,c\neq a}\omega^{(11)}_{bc}\gamma_a\gamma^b\gamma^c\Gamma\quad(\text{component }a),$$
the ω-sums running over $b,c\neq a$, i.e. the grade-3 and grade-11 parts. (Corrected after review 2026-10-30, R6-7; see `reports/CORRECTIONS.md`.)
with ($x_1=x$, $x_{13}=y$; $N_k=\nabla A_k$, $(N_k)_{ab}=\nabla_aA_{k,b}$; $w_{kk}=\tfrac12$, $w_{kl}=1$ for $k<l$)
$$E^{(i)}=-x_iG+c_ig+\sum_k\ell^{(i)}_k\big(N_k^\top-(\nabla\!\cdot\!A_k)\,g\big)+\sum_{k\le l}w_{kl}\big[b^{(i)}_{kl,1}A_k\otimes A_l+b^{(i)}_{kl,2}A_l\otimes A_k+b^{(i)}_{kl,3}\langle A_k,A_l\rangle g\big],$$
$$\omega^{(j)}=\sum_k\ell^{(j)}_k\,dA_k+\sum_{k<l}b^{(j)}_{kl}\,A_k\wedge A_l\qquad(j=3,11),$$
where $G=G(g_\aleph)$, $c_1,c_{13}$ are those of I6 (re-derived in §J0), and the coefficients are:

*Derivative coefficients* $\ell^{(\cdot)}_k$ (coefficient of $N_k^\top-\mathrm{tr}N_k\,g$ in $E^{(1)},E^{(13)}$; of $dA_k$ in $\omega^{(3)},\omega^{(11)}$):

| $k$ | (14,0): $E^{(1)}$ | $E^{(13)}$ | $\omega^{(3)}$ | $\omega^{(11)}$ | (7,7): $E^{(1)}$ | $E^{(13)}$ | $\omega^{(3)}$ | $\omega^{(11)}$ |
|---|---|---|---|---|---|---|---|---|
| $A$ ($P_0$) | 0 | 0 | $-x$ | $-y$ | 0 | 0 | $-x$ | $-y$ |
| $C$ ($P_2$) | $24x$ | $24y$ | $11iz$ | $11iw$ | $24x$ | $24y$ | $11iz$ | $11iw$ |
| $D$ ($P_{12}$) | $24iw$ | $-24iz$ | $11y$ | $-11x$ | $-24iw$ | $-24iz$ | $-11y$ | $-11x$ |
| $B$ ($P_{14}$) | 0 | 0 | $-iw$ | $iz$ | 0 | 0 | $iw$ | $iz$ |

*Quadratic coefficients* $(b_1,b_2,b_3)$ in $E^{(1)}$, $E^{(13)}$ and $b_4$ in $\omega^{(3)}$, $\omega^{(11)}$ (pairs not listed, including $(A,A)$ and $(B,B)$, are 0):

| $(k,l)$ | (14,0) | (7,7) |
|---|---|---|
| $(A,C)$ | $\omega^{(3)}:-\tfrac{22}3x$; $\omega^{(11)}:-\tfrac{22}3y$ | same |
| $(A,D)$ | $E^{(1)}:(0,-16y,16y)$; $E^{(13)}:(0,16x,-16x)$ | $E^{(1)}:(0,16y,-16y)$; $E^{(13)}:(0,16x,-16x)$ |
| $(A,B)$ | $\omega^{(3)}:\tfrac23y$; $\omega^{(11)}:-\tfrac23x$ | $\omega^{(3)}:-\tfrac23y$; $\omega^{(11)}:-\tfrac23x$ |
| $(C,C)$ | $E^{(1)}:(-224x,-224x,-176x)$; $E^{(13)}:(-224y,-224y,-176y)$ | same |
| $(C,D)$ | $E^{(1)}:(-32iw,-32iw,-352iw)$; $E^{(13)}:(32iz,32iz,352iz)$; $\omega^{(3)}:-\tfrac{266}3y$; $\omega^{(11)}:\tfrac{266}3x$ | $E^{(1)}:(32iw,32iw,352iw)$; $E^{(13)}:(32iz,32iz,352iz)$; $\omega^{(3)}:\tfrac{266}3y$; $\omega^{(11)}:\tfrac{266}3x$ |
| $(C,B)$ | $E^{(1)}:(16y,0,-16y)$; $E^{(13)}:(-16x,0,16x)$ | $E^{(1)}:(-16y,0,16y)$; $E^{(13)}:(-16x,0,16x)$ |
| $(D,D)$ | $E^{(1)}:(-160x,-160x,528x)$; $E^{(13)}:(-160y,-160y,528y)$ | $E^{(1)}:(160x,160x,-528x)$; $E^{(13)}:(160y,160y,-528y)$ |
| $(D,B)$ | $\omega^{(3)}:\tfrac{22}3x$; $\omega^{(11)}:\tfrac{22}3y$ | $\omega^{(3)}:-\tfrac{22}3x$; $\omega^{(11)}:-\tfrac{22}3y$ |

The twisted pair ($z,w$) occurs only in the rows $D$, $B$ of the first table, in $\ell^{(3)}_C,\ell^{(11)}_C$, and in the $(C,D)$ entries of $E^{(1)},E^{(13)}$ of the second.

In particular no other grade occurs, and the symmetric part of the derivative term is always $\ell(N^\top-\mathrm{tr}N\,g)$: the coefficient of $N$ is 0 and that of $\mathrm{tr}N\,g$ is $-\ell$.

*Proof.*
- The $T_o$-terms: $\circledast F_0+\kappa T_o+Q(T_o)$ is I6 (re-derived exactly in §J0), and $Q(T_e+T_o)$ has no odd cross term (I2).
- *Derivative term.* $P_k$ is built from $\nabla^0$-parallel objects (the soldering form, the Clifford product, Γ), so $\nabla(P_k(A))=P_k(\nabla A)$ and $(d_0P_k(A))_{ca}=P_k(N_{c\cdot})_a-P_k(N_{a\cdot})_c$ ($\nabla^0$ is torsion-free). The map $N\mapsto Y_{\rm odd}[\circledast(d_0P_k[N])]$ is linear in the 2-tensor $N$. Output §J2 evaluates it on all 196 basis tensors $e^p\otimes e^q$, for $k=0,2,12,14$ and each operator, and finds exactly the stated form. By linearity it holds for every $N$.
- *Quadratic term.* $(A,A')\mapsto Y_{\rm odd}[2B(P_k(A),P_l(A'))]$ is bilinear. Output §J2 evaluates it on every basis pair $(e^p,e^q)$ ($p\le q$ when $k=l$, using the symmetry of $B$) and finds exactly the stated form. By bilinearity it holds for all $(A,A')$, and $Q(T_e)=\sum_{k\le l}w_{kl}\,2B(P_k(A_k),P_l(A_l))$.
- §J3 confirms the assembled formula end to end on a random 2-jet. ∎

**Substituting the slaving.** With $A_k=\Sigma_{k,t}dt+\Sigma_{k,u}du$ one has $N_k^\top=\sum_\phi\Sigma_{k\phi}\nabla^2\phi+\sum_{\phi,\psi}\partial_\psi\Sigma_{k\phi}\,d\phi\otimes d\psi$ and $dA_k=(\partial_t\Sigma_{ku}-\partial_u\Sigma_{kt})\,dt\wedge du$. So
$$E^{(i)}=-x_iG+c_ig+\sum_\phi H_{i\phi}\big(\nabla^2\phi-\Box\phi\,g\big)+\sum_{\phi,\psi}\big[\mathcal Q^{(i)}_{\phi\psi}\,d\phi\otimes d\psi+\mathcal Q'^{(i)}_{\phi\psi}\langle d\phi,d\psi\rangle g\big],\qquad \omega^{(j)}=\varpi^{(j)}\,dt\wedge du,$$
with the **kinetic matrix** $H=\ell_E\Sigma=-\ell_EM^{-1}S$ ($\ell_E$ the $2\times4$ matrix of $\ell^{(1)}_k,\ell^{(13)}_k$). "Kinetic" means only that $H$ is the principal (second-derivative) coefficient of the reduced equations; it multiplies $\nabla^2\phi-\Box\phi\,g$, the nonminimal-coupling structure of the Remark below, and no reduced action is derived. Further, $\mathcal Q^{(i)}_{\phi\psi}=\sum_k\ell^{(i)}_k\partial_\psi\Sigma_{k\phi}+\sum_{k\le l}w_{kl}(b_1\Sigma_{k\phi}\Sigma_{l\psi}+b_2\Sigma_{l\phi}\Sigma_{k\psi})$, $\mathcal Q'^{(i)}_{\phi\psi}=-\sum_k\ell^{(i)}_k\partial_\psi\Sigma_{k\phi}+\sum_{k\le l}w_{kl}b_3\Sigma_{k\phi}\Sigma_{l\psi}$, and $\varpi$ collecting the antisymmetric parts. The equation is second order in $(t,u)$.

**Corollary J2.1 (where $z,w$ enter, and where they do not).**

(a) *Zeroth order.* $z,w$ do not occur in $-x_iG+c_ig$ (I6): the zeroth-order term $c_i(t,u)\,g$ is unchanged. The effective 4D cosmological constant of the reductions can still depend on $(z,w)$: under (J) with a normal gradient $a\neq0$ through $\mathfrak B_i(\alpha)$ and, in Z1, through $13\beta(Ha)_i$ (J4(c′)), and under (P) through the $\mathcal Q'\langle d\tau,d\tau\rangle g$ terms.

(b) *Kinetic matrix.* To leading order,
$$(14,0):\ H=\frac{576}{\kappa}\begin{pmatrix}x^2+w^2&xy-wz\\xy-wz&y^2+z^2\end{pmatrix},\qquad (7,7):\ H=\frac{576}{\kappa}\begin{pmatrix}x^2-w^2&wz-xy\\xy-wz&z^2-y^2\end{pmatrix},$$
with $\det H=576^2(xz\pm yw)^2/\kappa^2$ (upper sign (14,0)). In the $\mathfrak u(128)$ real form, after multiplying the equations by $-i$, the Hessian terms read $-\sum_\phi K_{i\phi}(\nabla^2\phi'-\Box\phi'\,g)$ with
$$K=\frac{576}{\kappa}\begin{pmatrix}x'^2+w'^2&x'y'-w'z'\\x'y'-w'z'&y'^2+z'^2\end{pmatrix}=\frac{576}{\kappa}\,\mathrm{Gram}\big((x',w'),(y',-z')\big),$$
positive semidefinite for $\kappa>0$ (negative semidefinite for $\kappa<0$) (Corrected after review 2026-10-30, R6-8; see `reports/CORRECTIONS.md`.), and definite iff $x'z'+y'w'\neq0$, which is H7's ellipticity condition. In (7,7), weighting the γΓ-channel by its Gram sign gives a symmetric form with determinant $-576^2(xz-yw)^2/\kappa^2$, indefinite whenever it is nondegenerate. The factor is H7's by construction: $\ell_E$ and $S$ are the $\Lambda^1$ blocks of the principal symbol of $\circledast\circ d$ (H4's $C_{\Lambda^0}$), so $H\approx-\ell_ES/\kappa$ is the symbol of $\circledast d\,\kappa^{-1}{*}^{-1}\circledast d$ on the invariant odd copies.

(b′) *Exactly, for every member* (after review 2026-10-30, R6-2). $\ell_E$ vanishes on the $A,B$ columns and $S$ on the $A,B$ rows, so $H=-\ell_{CD}(M^{-1})_{CD}S_{CD}$. Jacobi's complementary-minor identity gives $\det(M^{-1})_{CD}=\det M_{AB}/\det M=\kappa^2/\det M$, and $\det\ell_{CD}\det S_{CD}=576^2(xz\pm yw)^2$. Hence, wherever $\det M\neq0$,
$$\det H=\frac{576^2\,(xz\pm yw)^2\,\kappa^2}{\det M}\qquad(\text{upper sign (14,0)}).$$
Also, column $C$ of $\ell_E$ is $24(x,y)$, and $\det[\text{column }D,(x,y)]=24i(xz\pm yw)$. So whenever $xz\pm yw=0$, $\mathrm{im}\,H(t,u)\subset\mathrm{span}(x,y)$ exactly, for every $(t,u)$ (output §J2).

(c) *Exactly*, for $\circledast_\sharp$ (14,0):
$$H_\sharp=\frac{576\kappa}{\det M}\begin{pmatrix}\kappa^2+544\kappa t-3328u^2&-384\kappa u\\-384\kappa u&\kappa^2+224\kappa t-3328u^2\end{pmatrix},\qquad \det H_\sharp=\frac{576^2\kappa^2}{\det M},$$
nonzero wherever the slaving is defined. In (7,7) the same holds, with $-3328u^2\to+3328u^2$ and off-diagonal entries $\pm384\kappa u$ (output §J2), again with $\det H_\sharp=576^2\kappa^2/\det M$. For $W_A$: $H_{W_A}=\dfrac{576\kappa}{\kappa^2+224\kappa t-3328\epsilon u^2}\begin{pmatrix}1&0\\0&0\end{pmatrix}$, of rank 1: only $t$ has a second-order term, and $u$ enters only algebraically.

(d) *Gradient terms* are $O(\kappa^{-2})$ and depend on $z,w$ (leading coefficients in output §J2). The twisted pair enters $Q(T_e)$ only through the $(C,D)$ cross term, i.e. through $\tfrac23W_C\mathrm{Ad}_\Gamma(T_e\wedge T_e)$ (I3).

(e) *$\Lambda^2$-type equations.* To leading order, in (14,0),
$$\varpi^{(1)}=-46080\,w\,\frac{xz+yw}{\kappa^2},\quad \varpi^{(13)}=46080\,z\,\frac{xz+yw}{\kappa^2},\quad \varpi^{(3)}=-88320\,i\,y\,\frac{xz+yw}{\kappa^2},\quad \varpi^{(11)}=88320\,i\,x\,\frac{xz+yw}{\kappa^2}.$$
In (7,7) the same four read $46080\,w$, $46080\,z$, $88320\,i\,y$, $88320\,i\,x$ times $(xz-yw)/\kappa^2$. All four vanish to this order for $W_A$. Exactly, however, $\varpi^{(1)}_{W_A}=\tfrac12\ell^{(1)}_C\partial_u\Sigma_{C,t}=\pm1916928\,\kappa u/(\kappa^2+224\kappa t\mp3328u^2)^2$ (upper signs (14,0)), which is $O(\kappa^{-3})$ and nonzero for $u\neq0$; the other three vanish (output §J2). For $\circledast_\sharp$, exactly, $\varpi^{(13)}=9216\kappa^2\,(5\kappa^4-1600\kappa^3t-\ldots)/(\det M)^2$ (output §J2), which is nonzero near $T_o=0$. The antisymmetric part of $E^{(i)}$ at a point is $\varpi^{(i)}\,dt\wedge du$ (the Hessian, $G$ and $g$ are symmetric), and $\omega^{(j)}=\varpi^{(j)}dt\wedge du$. So at every point where $\Upsilon=0$ and some $\varpi\neq0$, $dt\wedge du=0$: for $\circledast_\sharp$ near $T_o=0$, for $W_A$ wherever $u\neq0$. (Corrected after review 2026-10-30, R6-5; see `reports/CORRECTIONS.md`.)

(f) *Hence $\circledast_\sharp$ and $W_A$ give different reduced equations once $T_e\neq0$*: for $\circledast_\sharp$ both scalars are second order (rank $H=2$); for $W_A$ only $t$ is (rank $H=1$). The $\Lambda^2$-constraint $dt\wedge du=0$ distinguishes them only at leading order and at $u=0$ (e). (Corrected after review 2026-10-30, R6-5; see `reports/CORRECTIONS.md`.) The answer of I7(e), "$z,w$ appear nowhere", is specific to $T_e=0$.

*Proof.* Substitution of J1 into J2. (b) and (e) to leading order: the Neumann series $\Sigma=-S/\kappa+(tM_1+uM_{13})S/\kappa^2+O(\kappa^{-3})$ (exact sympy). (b′): the vanishing pattern of $\ell_E$ and $S$, Jacobi's identity $\det\mathrm{adj}(M)_{CD}=\kappa^2\det M$ and $\det\ell_{CD}\det S_{CD}$, all as exact polynomial identities (output §J2). (c), (e) and the exact $\varpi^{(13)}_\sharp$, $\varpi_{W_A}$: exact sympy on $\Sigma=-M^{-1}S$. Real form: the substitution of I6. ∎

*Remark.* The combination $\nabla^2\phi-\Box\phi\,g$ is the one that the variation of a nonminimal coupling $f(\phi)R$ produces in a scalar–tensor theory. Here it arises from $\circledast(d_0T_e)$ with $T_e$ slaved; there is no $f(\phi)G$ term at this order.

**Proposition J3 (basic scalars (B)).** Let $t=\pi^*\tau$, $u=\pi^*\upsilon$, with the slaved $T_e$ of J1. On $s(X)$:
$\mathrm{Hess}\,\phi|_{HH}=\mathrm{Hess}_g\tau_\phi$, $\mathrm{Hess}\,\phi|_{VV}=0$, and $\mathrm{Hess}\,\phi(X,U)=0$ in Z2, $=g_\aleph(\mathrm{II}(X,\mathrm{grad}\,\phi),U)$ in Z1, with II of F4(a); $d\phi$ is horizontal. Then:

(a) *(P), Z2.* $s^*({*}\Upsilon)=0$ is the 4D system ($i=1,13$; all operations for $g$)
$$-x_i\big(G(g)+\Lambda_Pg\big)+c_Hc_i(\tau,\upsilon)\,g+\sum_\phi H_{i\phi}\big(\nabla^2\tau_\phi-\Box\tau_\phi\,g\big)+\sum_{\phi,\psi}\big[\mathcal Q^{(i)}_{\phi\psi}d\tau_\phi\otimes d\tau_\psi+\mathcal Q'^{(i)}_{\phi\psi}\langle d\tau_\phi,d\tau_\psi\rangle g\big]=0,$$
together with $\varpi^{(j)}\,d\tau\wedge d\upsilon=0$. The coefficients are those of J2, unchanged, and independent of $c_H$.

(b) *(P), Z1.* (P) keeps the $\gamma^A$-components $E^{(i)}_{\mu A}$. These are $\sum_\phi H_{i\phi}\,g_\aleph(\mathrm{II}(\partial_\mu,\mathrm{grad}\,\phi),\partial_A)$, so (P) contains $\mathrm{II}(X,\sum_\phi H_{i\phi}\mathrm{grad}\,\tau_\phi)=0$ for all $X$, i.e. $H\,(d\tau,d\upsilon)^\top=0$. If $\det H\neq0$, then $d\tau=d\upsilon=0$, so $T_e=0$ and (P) reduces to I7(a). By Cor. J2.1(b′) this holds for **every** member with $xz\pm yw\neq0$ wherever $\det M\neq0$, in particular $\circledast_\sharp$. For $W_A$ it gives $d\tau=0$, hence again $T_e=0$, with $\upsilon$ as in I6.

(c) *(E) on $s(X)$.* The VV block is $-x_iG^Y_{VV}+(\text{scalar})\,g_{VV}$. For $(x,y)\neq0$ it requires $\mathrm{Ric}^Y_{VV}\propto g_{VV}$, i.e. Z1 with $\lambda=\tfrac15$ (F3, F4(c)); the HV block is that of (b). So for every member with $xz\pm yw\neq0$ (where $\det M\neq0$), and for $W_A$, the solutions are exactly those of I7(b): none in (7,7), and in (14,0) $\mathrm{Ric}(g)=0$, $\Lambda_4=0$. For the other members with $xz\pm yw=0$ this is not shown.

*Proof.*
- Hessian of a basic function on $s(X)$. HH: $\nabla_XY=\nabla^s_XY+\mathrm{II}(X,Y)$ with II vertical, and vertical derivatives of $\pi^*\tau$ vanish. VV: the fibres are totally geodesic (F3, F4(a)). HV: $\mathrm{Hess}(X,U)=-(\nabla_XU)\phi=g(\nabla_X\mathrm{grad}\,\phi,U)$, which is $g(A_X\mathrm{grad}\,\phi,U)=0$ in Z2 (A vanishes on $s(X)$, F3) and $g(\mathrm{II}(X,\mathrm{grad}\,\phi),U)$ in Z1. Also $\Box_Y\phi=c_H^{-1}\Box_g\tau_\phi$ and $\langle d\phi,d\psi\rangle_{g_\aleph}=c_H^{-1}\langle d\tau_\phi,d\tau_\psi\rangle_g$, and $s^*G^Y=G(g)+\Lambda_Pg$ (F5); this gives (a).
- (b): $\mathrm{II}(Z,Z)=-\tfrac{c_H}{2c_V}[Z^\flat\otimes Z^\flat+\mu|Z|^2g]\neq0$ for $Z\neq0$ (rank count, also for null $Z$). For $W_A$, $H$ has only the entry $H_{1t}\neq0$, and $\Sigma$ has no $du$ column.
- (c): with $\mathrm{Hess}_{VV}=0$ and $d\phi$ horizontal, every $(t,u)$-term in the VV block is a multiple of $g_{VV}$. F6's argument then applies. ∎

**Proposition J4 ((E) on $s(X)$ for the invariant 2-jet (J)).** Assume (J) at the points of $s(X)$ and the slaved $T_e$ of J1 near $s(X)$.

(a) At $p\in s(X)$ all $(t,u)$-terms of $Y_{\rm odd}$ lie in $\mathrm{span}(\mathfrak h,G_0,\tau\otimes\tau)$, and the antisymmetric and $\Lambda^2$ parts vanish ($dt\wedge du=0$). In Z1, $\beta_\phi=a_\phi/(2c_V(1-4\lambda))$, from $\mathrm{Hess}\,\phi|_{HH}=-a_\phi\,\tau(\mathrm{II})$; in Z2, $\beta_\phi=0$. The HV block vanishes. For $(x,y)\neq0$ the traceless HH part forces $g$ Einstein, $\mathrm{Ric}(g)=\Lambda_gg$, i.e. $G(g)+\Lambda_gg=0$. So $\Upsilon=0$ at the points of $s(X)$ is equivalent to six scalar equations: the $\mathfrak h$, $G_0$ and $\tau\otimes\tau$ coefficients of $E^{(1)}$ and $E^{(13)}$.

(b) *Balance identity.* For $i=1,13$ and any slaving matrix Σ,
$$x_i\Big(\frac{\Lambda_g}{c_H}+\rho\Big)=\sum_\phi H_{i\phi}\Big(\beta_\phi-\frac{p_\phi}{c_V}\Big),\qquad \rho_{\rm Z2}=\frac1{c_V},\quad \rho_{\rm Z1}=\frac{5\lambda-1}{c_V(4\lambda-1)} .$$
It is the difference of the $\mathfrak h$ coefficient and $c_V^{-1}$ times the $G_0$ coefficient of $E^{(i)}$. With $H=0$ ($T_e\equiv0$) it gives $\Lambda_g=-c_H\rho$, which vanishes exactly at Z1, $\lambda=\tfrac15$, consistent with F6 (there the other equations then force that case).

(c) *Branch $a_t=a_u=0$* ($T_e=0$ on $s(X)$, $\nabla T_e\neq0$). The six equations are equivalent to
$$c_1(t,u)=x\,\Theta,\qquad c_{13}(t,u)=y\,\Theta,\qquad \Theta=\mathrm{tr}_VG^Y_{VV}-9\,G^Y_{\mathfrak h}=\Theta_0-\frac{11}{c_H}\Lambda_g,$$
together with $\sum_\phi H_{i\phi}p_\phi=x_i\big(G^Y_{G_0}-c_VG^Y_{\mathfrak h}\big)$ and $\sum_\phi H_{i\phi}q_\phi=x_i\big(G^Y_{\tau\tau}+\lambda c_VG^Y_{\mathfrak h}\big)$, where $\Theta_0=-\tfrac{9}{2c_V}$ in Z2 and $\Theta_0=\tfrac{3(2-15\lambda)}{2c_V(4\lambda-1)}$ in Z1. The last two are solvable for $(p,q)$ whenever $(x,y)^\top\in\mathrm{im}\,H(t,u)$. By Cor. J2.1(b′) this holds exactly: if $xz\pm yw\neq0$, $H$ is invertible wherever $\det M\neq0$; if $xz\pm yw=0$, $\mathrm{im}\,H(t,u)\subset\mathrm{span}(x,y)$ for every $(t,u)$, so the condition is $H(t,u)\neq0$. In particular it holds for $\circledast_\sharp$ wherever $\det M\neq0$, for $W_A$ always, and near $t=u=0$ for every member with $H_0\neq0$ (in the $\mathfrak u(128)$ real form: every nonzero member). (This replaces a leading-order argument; after review 2026-10-30, R6-2.)
So, for $x\neq0$: $(t,u)$ lies on the curve $yc_1=xc_{13}$ of I7(a), and
$$\Lambda_g=\frac{c_H}{11}\Big(\Theta_0-\frac{c_1(t,u)}{x}\Big).$$
This branch contains the solutions of I7(b): at Z1, $\lambda=\tfrac15$, $\Lambda_g=0$ one has $\Theta_0=\tfrac{15}{2c_V}$, which is I7(b)'s condition $e_1=-\tfrac{15}{2c_V}$ ($c_1/x=-e_1$); there $G^Y\propto g_\aleph$, so the right-hand sides vanish and $p=q=0$ is allowed (constant $t,u$, $T_e\equiv0$).

(c′) *Every normal gradient $a$* (after review 2026-10-30, R6-1). Put $\alpha_k=\sum_\phi\Sigma_{k\phi}a_\phi$, so that $A_k=\alpha_k\tau$ on $s(X)$, and $\mathfrak B_i(\alpha)=\sum_{k\le l}w_{kl}\,(b_1+b_2+b_3)^{(i)}_{kl}\,\alpha_k\alpha_l$. Then for every Σ, $\partial\Sigma$, $a$, $p$, $q$:
$$-9\,[\mathfrak h]+\mathrm{tr}_{g_\aleph}G_0\,[G_0]+|\tau|^2\,[\tau\tau]=c_i(t,u)-x_i\Theta-13\beta\,(Ha)_i+|\tau|^2\,\mathfrak B_i(\alpha),\qquad |\tau|^2=\frac{4}{c_V(1-4\lambda)},$$
where $[\cdot]$ are the three coefficients of $E^{(i)}$. The other two coefficients fix $\big(\sum_\phi H_{i\phi}p_\phi,\ \sum_\phi H_{i\phi}q_\phi\big)$ uniquely, through a scalar $2\times2$ system of determinant $-9$. So, where $\det H\neq0$, the six equations are equivalent to
$$c_i(t,u)=x_i\Theta+13\beta\,(Ha)_i-|\tau|^2\,\mathfrak B_i(\alpha)\qquad(i=1,13),$$
with $(p,q)$ then determined. $\partial\Sigma$ drops out; $13=9+4$; at $a=0$ this is (c). This closed form agrees with the review's (`xcheck_round6.out` X4b).

(d) *Consequences.*
- (E) at the points of $s(X)$ has solutions with $T_e\not\equiv0$ near $s(X)$ in **both signatures, in Z1 and Z2, for every $\lambda\neq\tfrac14$ on the side the signature requires** ($\lambda<\tfrac14$ for definite (14,0), $\lambda>\tfrac14$ for (7,7), A5), for every member with $(x,y)\neq0$ satisfying the image condition (in particular $\circledast_\sharp$ and $W_A$). (Corrected after review 2026-10-30, R6-8; see `reports/CORRECTIONS.md`.)
- The 4D cosmological constant $\Lambda_4:=\Lambda_g$ is **not fixed**. For $\circledast_\sharp$ it takes every real value in (7,7). In (14,0), with $c_H>0$, **on the branch $a=0$** it takes exactly the values in $[\tfrac{c_H}{11}(\Theta_0-\tfrac{\kappa^2}{1248}),\infty)$, since $c_1/x=-e_1$ and $e_1\ge-\kappa^2/1248$ (I7(a)). With $a\neq0$ this is the full range in Z2, while in Z1 the full range is larger (d′). (Corrected after review 2026-10-30, R6-1; the range had been stated for all of (J), without the restriction $a=0$; see `reports/CORRECTIONS.md`.)
- Exact sample solutions (Z1 and Z2; (14,0) with $\lambda=0$ in the $\mathfrak u(128)$ real form, (7,7) with $\lambda=\tfrac12$) are verified in all six components in output §J4.
- So I7(b)'s "no solution in (7,7)" and "$\Lambda_4=0$ in (14,0)" are properties of $T_e\equiv0$, not of (E) at $s(X)$.

(d′) *(14,0) range for $\circledast_\sharp$, every $a$* (after review 2026-10-30, R6-1). In the $\mathfrak u(128)$ real form the slaved $T_e$ is $\mathfrak u(128)$-valued (the slaving preserves the real form), so $\alpha_C\in\mathbb R$ and $\alpha_D=i\delta$ with $\delta\in\mathbb R$ (also $\alpha_A\in i\mathbb R$, $\alpha_B\in\mathbb R$; output §J4b at random exact points). For $y=w=0$ (any $z$) the tables give $\mathfrak B_1=x(-312\alpha_C^2+104\alpha_D^2)$ and $(Ha)_1=24x\,\alpha_C$, and $c_1/x=-e_1$ with $e_1=\kappa s+312s^2+104r^2\ge-\kappa^2/1248$ for all real $s,r$. Hence, completing the square,
$$\Theta=-e_1+\frac{78\beta^2}{|\tau|^2}-312|\tau|^2\Big(\alpha_C+\frac{\beta}{2|\tau|^2}\Big)^2-104|\tau|^2\delta^2\ \le\ \frac{\kappa^2}{1248}+\frac{78\beta^2}{|\tau|^2},$$
for every $a$ and every $(t,u)$. With $\Lambda_g=\tfrac{c_H}{11}(\Theta_0-\Theta)$ and $c_H>0$:
- **Z2** ($\beta=0$): $\Lambda_g\ge\tfrac{c_H}{11}(\Theta_0-\tfrac{\kappa^2}{1248})$ for every $a$, so with (d) the range is exactly this half-line.
- **Z1** ($78\beta^2/|\tau|^2=\tfrac{39}{8c_V(1-4\lambda)}$): the range is exactly
$$\Big[\tfrac{c_H}{11}\Big(\Theta_0-\tfrac{\kappa^2}{1248}-\tfrac{39}{8c_V(1-4\lambda)}\Big),\ \infty\Big).$$
  The bound holds for every $a$. It is attained on the slice $u=a_u=0$, where the channel-13 equation holds identically and $\alpha_C=24a_s/(\kappa+224s)$; the slice is connected, has $\det M\neq0$ for $s\ge-\kappa/624$ ($\kappa>0$), attains the minimum at $s=-\kappa/624$, $\alpha_C=-\beta/(2|\tau|^2)$, and is unbounded above.
- Exact configurations (Z1, $\lambda=0$, $c_H=c_V=\kappa=1$, real form; output §J4b): the review's, with $t=i/1000$, $a_t=-51i/16000$ and $\Lambda_g=-984211/1375000\approx-0.7158$, below the $a=0$ bound $-3745/13728\approx-0.2728$ (all six equations vanish, $\det M=29529/15625$); and one at the bound, $t=-i/624$, $a_t=-25i/14976$, $\Lambda_g=-9829/13728$.

The review left the full Z1 range open; it is settled here for $\circledast_\sharp$.

(e) *Radial sub-case* ($p_\phi=0$). In Z2 the balance identity gives **$\Lambda_g=-c_H/c_V$** for every member with $(x,y)\neq0$, every $(z,w)$ and every $\lambda\neq\tfrac14$. In Z1 it gives $\Lambda_g=-c_H\rho_{\rm Z1}+\tfrac{c_H}{x}\beta\sum_\phi H_{1\phi}a_\phi$ (for $x\neq0$), i.e. $\Lambda_g=\tfrac{c_H(1-5\lambda)}{c_V(4\lambda-1)}$ at $a=0$ (zero exactly at $\lambda=\tfrac15$, the F6 value), shifted by the normal gradient $a$. In (14,0) Z2, $\Lambda_g=-c_H/c_V$ lies in the range of (d′) iff $c_V\kappa^2\ge8112$: a necessary condition for a radial 2-jet with $\circledast_\sharp$, already pointwise at $s(X)$ (output §J4b).

*Proof.*
- (a): by the definition of (J). The HV part of $\mathrm{Hess}\,\phi$ is $a_\phi\nabla^\perp\tau$. The normal connection on $V|_{s(X)}$ is the one induced by ℵ: Koszul gives $2g(\nabla_XU,W)=Xg(U,W)+g([X,U],W)-g([X,W],U)$ for $U,W$ vertical, which involves only the vertical metric and the brackets, the same in Z1 and Z2. τ is ℵ-parallel, so the HV part vanishes for constant $a_\phi$. $G^Y$ on $s(X)$ is block-diagonal with the F3/F4 blocks (output §J4 prints them).
- (b), (c), (c′): polynomial identities in generic symbols for the entries of $\Sigma$ and $\partial\Sigma$ (output §J4, §J4b), hence valid for every slaving matrix.
- (d): substitution of the actual Σ, exact sympy. The ranges need $\det M\neq0$ on the branches used.
  - On $u=0$ (resp. $r=0$), $\det M=\kappa^2(\kappa+224\sigma)(\kappa+544\sigma)$ with $\sigma=t$ (resp. $s$). For $\kappa>0$ this is nonzero for $\sigma\ge-\kappa/624$, which already gives all values $c_1\ge-\kappa^2/1248$ (resp. $e_1\ge-\kappa^2/1248$).
  - On the (7,7) branch $t=-\kappa/208$, $\det M=v^2+\tfrac{554}{13}\kappa^2v+\tfrac{21}{169}\kappa^4>0$ with $v=3328u^2$. This branch gives $c_1=\kappa^2/416-104u^2$, i.e. all values $\le\kappa^2/416$.
  - $\kappa<0$ is symmetric.
- (d′): the reality of $\alpha$ follows from the real structure (I.0) and is confirmed at random exact points. The table identities for $\mathfrak B_1$, $(Ha)_1$ and $e_1$ and the completed square are exact sympy identities. The two configurations are solved and checked exactly (output §J4b).
- (e): (b) with $p=0$; the $c_V\kappa^2\ge8112$ condition is (d′) with $\Lambda_g=-c_H/c_V$ and $\Theta_0=-\tfrac9{2c_V}$. ∎

**Theorem J5 (radial ansatz on a neighbourhood, Z2).** Let Z2 hold with $\lambda\neq\tfrac14$, let $(x,y)\neq0$ and $\kappa\neq0$. Let $U\subset Y$ be a connected open set meeting $s(X)$, on which $t=f(s)$ and $u=h(s)$, $\det M(t,u)\neq0$, $T_e$ is the slaved field of J1, and $\Upsilon=0$.

(a) $\nabla ds=0$ on $Y$ (for every base metric $g$).

(b) There are functions $\alpha_i,\beta_i$ on $U$ with $x_i\,\mathrm{Ric}(g_\aleph)=\alpha_ig_\aleph+\beta_i\,ds\otimes ds$ on $U$. The leaves $\{s=\mathrm{const}\}\cap U$ are Einstein. At $p\in U\cap s(X)$ the Einstein constant is $-1/c_V$, and $\mathrm{Ric}(g)=-\tfrac{c_H}{c_V}\,g$ at $\pi(p)$.

(c) **There is no such solution**, in any signature of $g_\aleph$ (in particular (14,0) and (7,7)), for any $(z,w)$. (The proof below covers definite $g_\aleph$; the indefinite case is settled after round 6 in `J5d-indefinite-leaf.md`, Cor. L2.)

(d) *Indefinite signature* (settled 2026-10-30 in `J5d-indefinite-leaf.md`). A solution needs, besides (b), $\sum_i\epsilon_i\,g_\aleph(A_XE_i,A_YE_i)=0$ on the leaf through $p$ (O'Neill's A-tensor, $E_i$ an orthonormal horizontal frame). At second order in the fibre direction this reads $F(X,Y;h):=\sum_i\epsilon_iG_0(R(X,E_i)\cdot h,R(Y,E_i)\cdot h)=0$ for all $X,Y$ and all $g$-traceless $h$. Two traces of $F$ satisfy $8C_2-C_1=28|\mathrm{Ric}|^2$ for every algebraic curvature tensor in dimension 4 (J5d, Theorem L1), so $\Lambda_g=0$, against $\Lambda_g=-c_H/c_V$ from (b). For a Lorentzian base a single evaluation already forces $R_g(x_0)=0$ (J5d, Prop. L3). This covers (7,7) from a Lorentzian base and, locally, from a (2,2) base with $\lambda>\tfrac14$ (A5; A7.6 excludes a global Spin(7,7) structure over $\mathrm{Met}_{(2,2)}$). The constant-curvature failures of output §J5 are special cases. (Corrected after review 2026-10-30, R6-9, for the (2,2) base; until J5d this item said that existence in (7,7) was open.)

*Proof.*
- (a): the Euler field $E=y\,\partial_y$ generates $y\mapsto e^\sigma y$, an isometry of Z2 (I2.1's proof), so $E$ is Killing. Its dual is $E^\flat=c_VG_\lambda(y,Dy)=c_V(1-4\lambda)\,\mathrm{tr}(y^{-1}Dy)$. Also $\mathrm{tr}(y^{-1}Dy)=d\log|\det y|-2\,\mathrm{tr}\,\Gamma_\mu dx^\mu=d\log|\det(g^{-1}y)|=4\,ds$. So $E^\flat=4c_V(1-4\lambda)\,ds$ is closed, and a Killing field with closed dual is parallel. $|ds|^2=1/(4c_V(1-4\lambda))\neq0$.
- (b): with $t=f(s)$, $u=h(s)$ the slaved $A_k=F_k(s)\,ds$, so $N_k=F_k'(s)\,ds\otimes ds$ by (a). Every $(t,u)$-term of J2 lies in $\mathrm{span}(g_\aleph,ds\otimes ds)$, and $dt\wedge du=dA_k=A_k\wedge A_l=0$. So $E^{(i)}=0$ reads $x_i\mathrm{Ric}^Y=\alpha_ig_\aleph+\beta_ids\otimes ds$, and since $(x,y)\neq0$ some $x_i\neq0$. By (a), $\nabla s$ is a parallel non-null field, so $U$ is locally a product $I\times Y'$ and $\mathrm{Ric}^Y|_{TY'}=\mathrm{Ric}^{Y'}$ is a multiple of $g_{Y'}$. By Schur ($\dim Y'=13\ge3$) the factor is constant on each connected leaf. At $p\in s(X)$, $T_pY'=H\oplus V_0$ ($V_0=\ker\tau$). By F3, $\mathrm{Ric}^Y=\mathrm{Ric}_g$ on $H$ and $-G_0=-c_V^{-1}g_\aleph$ on $V_0$. Hence the constant is $-1/c_V$ and $\mathrm{Ric}_g=-\tfrac{c_H}{c_V}g$ at $\pi(p)$.
- (c): $Y'\to X$ is a pseudo-Riemannian submersion onto $(X,c_Hg)$ whose fibres $\{\det(g^{-1}y)=1\}$ are totally geodesic (level sets of the affine function $s$ inside the totally geodesic fibres of Z2). Its horizontal space is $H$, since horizontal transport preserves $\det(g^{-1}y)$, and its A-tensor is that of $Y$. O'Neill (Besse 9.36, $T=0$) gives $\mathrm{Ric}^{Y'}(X,X)=\mathrm{Ric}_g(X,X)-2(A_X,A_X)$. On the leaf through $p$, at the points over $x_0=\pi(p)$, both $\mathrm{Ric}^{Y'}$ and $\mathrm{Ric}_g(x_0)$ equal $-\tfrac{c_H}{c_V}g$, so $(A_X,A_X)=0$ at all points $(x_0,y)$ of an open subset $W$ of the unimodular slice. If $g_\aleph$ is definite, $A=0$ there. By L-F5, $A_XY=-\tfrac12(R_g(X,Y)\cdot y)^{\rm vert}$, so $R_g(X,Y)\cdot y=0$ for $y\in W$. $W$ spans $\mathrm{Sym}^2$: it contains $g$, and its tangent space at $g$ is the traceless part. So $B^\top y+yB=0$ for all symmetric $y$, with $B=R_g(X,Y)$; taking $y=e_i\otimes e_i$ gives $B_{ik}=0$ for all $k$, so $B=0$. Hence $R_g(x_0)=0$, contradicting $\mathrm{Ric}_g(x_0)=-\tfrac{c_H}{c_V}g\neq0$.
- (d): expand $(A_X,A_Y)$ at $y=g+\varepsilon h$. $R\cdot g=0$ and $R\cdot h$ is $g$-traceless, so the $\lambda$-term of $G_\lambda$ drops. The constant-curvature computation is exact in output §J5; the quadratic form scales with $k^2$. The general case is `J5d-indefinite-leaf.md` (L0 re-derives the expansion; L1, L2), which also gives (c) in indefinite signature. ∎

**Answers to the three questions.**
1. **Yes, $z$ and $w$ enter once $T_e\neq0$; they do not cancel.** They enter through the source $S$ ($D$-row), the slaving block $M$ (via $\theta=tw-uz$), the derivative map ($\ell_D$, and $\ell_C,\ell_B$ in grades 3, 11) and the $(C,D)$ cross term of $Q(T_e)$. In the reduced odd equation they appear in the principal (Hessian) coefficient matrix $H=-\ell_EM^{-1}S$ (to leading order $\tfrac{576}{\kappa}$ times a Gram matrix; exactly $\det H=576^2(xz\pm yw)^2\kappa^2/\det M$), in the $O(\kappa^{-2})$ gradient terms, and in the $\Lambda^2$-type equations ($\propto xz\pm yw$ to leading order). They do not enter the zeroth-order term $c_ig$; the effective $\Lambda_4$ of J4 depends on them once $a\neq0$ (J4(c′)). For $\circledast_\sharp$ both scalars carry second derivatives; for $W_A$ only one does.
2. **Existence near $s(X)$.**
   - (E) at the points of $s(X)$, invariant 2-jet: solutions exist in (14,0) and in (7,7), in Z1 and Z2, for every member with $(x,y)\neq0$ and the image condition of J4(c) (exactly: always if $xz\pm yw\neq0$ and $\det M\neq0$, otherwise iff $H(t,u)\neq0$; Cor. J2.1(b′)). This includes $\circledast_\sharp$ and $W_A$, so here $z,w$ change only how the second normal derivatives realise the solution.
   - Basic scalars: (E) gives back exactly I7(b) for every member with $xz\pm yw\neq0$ (where $\det M\neq0$) and for $W_A$.
   - $\Upsilon=0$ on a neighbourhood, radial ansatz, Z2 ($\kappa\neq0$, $T_e$ the slaved field in $\mathcal V$, $\det M\neq0$): no solution in any signature, (14,0) and (7,7) included, for any $(x,y,z,w)$ with $(x,y)\neq0$ (J5(c); indefinite case `J5d-indefinite-leaf.md`).
3. **$\Lambda_4$ and 4D field content.**
   - (P), basic scalars, Z2: a 4D system in $(g,\tau,\upsilon)$: an Einstein-type equation with Hessian terms $H_{i\phi}(\nabla^2\tau_\phi-\Box\tau_\phi g)$ and a field-dependent cosmological term $c_Hc_i(\tau,\upsilon)-x_i\Lambda_P$; for $y=0$ (including $\circledast_\sharp$) a curvature-free tensor constraint on the scalars (the γΓ-channel); and $\varpi\,d\tau\wedge d\upsilon=0$. Its solutions are not classified; nonconstant scalars restrict $g$ strongly (J.3). $H$ has rank 2 for $\circledast_\sharp$ and rank 1 for $W_A$. In Z1, (P) forces the scalars constant ($\det H\neq0$), which is I7(a) again. (Corrected after review 2026-10-30, R6-4; the text had "the 4D theory is Einstein gravity coupled to two scalars"; see `reports/CORRECTIONS.md`.)
   - (E) at the points of $s(X)$, 2-jet (J): $\Lambda_4$ is free. For $\circledast_\sharp$ it takes every real value in (7,7); in (14,0) ($c_H>0$) it takes the half-line $[\tfrac{c_H}{11}(\Theta_0-\tfrac{\kappa^2}{1248}),\infty)$ in Z2 and the larger half-line $[\tfrac{c_H}{11}(\Theta_0-\tfrac{\kappa^2}{1248}-\tfrac{39}{8c_V(1-4\lambda)}),\infty)$ in Z1 (J4(d), (d′)). (Corrected after review 2026-10-30, R6-1; the text had "a half-line in (14,0)" for both metrics; see `reports/CORRECTIONS.md`.)
   - Radial: $\Lambda_4=-c_H/c_V$ in Z2 (in (14,0) for $\circledast_\sharp$ this needs $c_V\kappa^2\ge8112$).

## J.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| §J0 | $c_1,c_{13}$ re-derived with the round-5 evaluator, equal to I6 | both signatures |
| J1 | exact evaluation on all 14 basis covectors, both maps, four operators; linearity | all $v$; both signatures |
| J2 | exact evaluation on all 196 basis tensors (derivative map) and all basis pairs (quadratic map), four operators; (bi)linearity | all $A_k$; both signatures |
| Cor. J2.1 | exact sympy: Neumann series to $O(\kappa^{-2})$; exact $H_\sharp$, $H_{W_A}$, $\varpi^{(13)}_\sharp$, $\varpi_{W_A}$; exact $\det H$ for every member (block structure plus Jacobi, polynomial identities) | family; both signatures |
| §J3 | end-to-end on a random 2-jet (seed fixed), member $(2,-1,3,1)$, $\kappa=5$, $(t,u)=(\tfrac17,-\tfrac13)$: $Y_{\rm ev}=0$ exactly; odd part equals the assembled formula exactly (scaled Gaussian integers) | spot check of the assembly |
| J3 | proof from F3, F4(a), F5, F6 | Z1, Z2; (P), (E) |
| J4 | proof; balance identity, branch $a=0$ and the general-$a$ closed form (c′) as polynomial identities for generic Σ, $\partial\Sigma$; exact sample solutions | Z1, Z2; both signatures |
| J4(d′) | reality of α at 4 random exact points plus the real structure; table identities and completed square exact; two exact Z1 configurations (the review's, and one at the bound) | (14,0), $\circledast_\sharp$ |
| J5 | proof; constant-curvature computation exact (sympy), bases (4,0), (1,3), (2,2); every signature: `J5d_leaf.py` | Z2 |
| Review | independent cross-check `xcheck_round6.py` (REVIEW4-GU), 90/90 | all of J1–J5 |

## J.3 Limits (layers b, d)

- **Not covered.**
  - $T_o$ outside the invariant sector (e.g. the soldering sector of I5 with $T_e\neq0$). There $L_{T_o}$ is not equivariant and the even equation does not close on a small sector.
  - Components of $T_e$ outside $\mathcal V$. The slaved solution has none where $L_{T_o}$ is invertible on all of $\Omega^1(\mathfrak g_{\rm ev})$ (e.g. $|T_o|$ small). Only the $\mathcal V$-block $M$ of $L_{T_o}$ was computed, so for large $|T_o|$ there may be further branches.
  - The exact gradient coefficients $\mathcal Q,\mathcal Q'$ and $\varpi$ are printed only to $O(\kappa^{-2})$ (and $\varpi^{(13)}$ exactly for $\circledast_\sharp$). The exact ones are used numerically in §J3.
  - $\Upsilon=0$ on a neighbourhood beyond the radial Z2 ansatz. In particular: Z1 (the Euler field is not Killing there, so $s$ is not parallel); non-radial $t,u$ (one only knows $dt\wedge du=0$ at points where some $\varpi\ne0$).
  - The (14,0) range of J4(d′) for members other than $\circledast_\sharp$ (it uses $y=w=0$), and the J3(b), (c) conclusions for members with $xz\pm yw=0$ other than $W_A$.
  - A classification of the solutions of the 4D system of J3(a). For $\circledast_\sharp$ its γΓ-channel has no curvature term, and its traceless part makes a nonconstant scalar concircular after reparametrisation, unless its reduced Hessian coefficient vanishes. Not analysed further.
  - Reading (A); the signatures (8,6), (10,4) of the F6 case; the ℷ-variation.
  - In Z1 every statement uses F4(c), which is now proven for every $n$ (`K-z1-ricci-general-n.md`, 2026-10-30).
- **Precise statements replacing I.3's "not covered: slaved $T_e\neq0$".**
  1. The twist changes the reduced odd equation once $T_e\neq0$ (Cor. J2.1), but not its zeroth-order part.
  2. Reading (E) at the points of $s(X)$ no longer obstructs (7,7) and no longer fixes $\Lambda_4$ (J4). With free second normal derivatives it is as underdetermined as (P).
  3. A genuine neighbourhood condition can still exclude solutions: radial Z2, every signature (J5, J5d).
- **What would make the reduction well-posed.** As in I.3, a rule for $T$ off $s(X)$, now including its second normal derivatives (which feed $\nabla T_e$), or $\Upsilon=0$ on a neighbourhood with fibre boundary conditions. J4 and J5 show that the answer depends on which.

## J.4 Consequences for gaps and claims

- **GU-G14.** Refined. Round 5's (E) results ("none in (7,7)", "$\Lambda_4=0$ in (14,0)") hold for $T_e\equiv0$. With the slaved $T_e$ and an invariant 2-jet along $s(X)$, (E) at $s(X)$ has solutions in both signatures with $\Lambda_4$ free (J4; in (14,0) for $\circledast_\sharp$ a half-line whose end depends on Z1/Z2, J4(d′)). In the radial Z2 sector ($\kappa\ne0$, slaved $T_e\in\mathcal V$, $\det M\ne0$) no solution on a neighbourhood exists in any signature (J5; indefinite case J5d). Under (P) with basic scalars in Z2 one gets a 4D system with two scalars, Hessian terms and, for $y=0$, a curvature-free constraint on the scalars (J3(a)); its solutions are not classified.
- **GU-G16.** The ellipticity factor $xz\pm yw$ reappears in $\det H=576^2(xz\pm yw)^2\kappa^2/\det M$ (exact, every member) and as a factor of every leading $\Lambda^2$ coefficient. For $H$ this holds by construction (H4's $C_{\Lambda^0}$ blocks).
- **GU-G15.** The slaving of I2.2 is explicit on $\mathcal V$ (J1).
- **GU-054.** Status C unchanged. The (P) and (E) cosmological constants of round 5 are not predictions once $T_e\neq0$ (J4). No statement made before round 6 is contradicted: round 5's (E) statements are all scoped to the invariant sector $T=t\Phi^1+u\Phi^1\Gamma$, which has $T_e\equiv0$ by definition (confirmed by REVIEW4-GU, R6-11). Those summaries do not say "$T_e\equiv0$" explicitly; adding it is left to the files that carry them. This note's own J4(d) range (round 6) was too broad and is corrected in place (R6-1; `reports/CORRECTIONS.md`).
