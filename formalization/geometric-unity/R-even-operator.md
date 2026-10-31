# R. The even operator $L_{T_o}$ on all of $\Lambda^1\otimes\mathfrak g_{\rm ev}$: where the slaved field leaves $\mathcal V$

Claims: GU-054 (context). Gaps: GU-G15, GU-G14 (and GU-G16 for the role of $\ell=xz\pm yw$). Code: `computation/geometric-unity/R_even_operator.py` → `.out` (sections R0–R4 with R1b–R1d; exact: Python Fractions, sympy rationals and Gaussian rationals; no randomness; 342 s on one core, a second run 353 s, identical apart from the timing annotations). It imports `J_slaved_sector.py` (and through it `I_elliptic_reduction.py`, `G_symbols.py`, `H_ellipticity.py`, `shiab_checks.py`), unmodified, and calls J's own `j1`. Output sections: §R0 constituents, §R1 blocks, §R1b equivariance, §R1c source, §R1d traces, §R2 determinants, §R3 singular sets and kernels, §R4 consequences. Record: `metrics/geometric-unity/R-even-operator.yaml`. Builds on I (I2, Cor. I2.2), J (J.0, J1), M (M1), H (H7) and, for R9, I7(b), N12(c), Q15. Layers: (b) reconstruction (that of J, unchanged), (c) proven, (d) open.

## R.0 Question and setting (layer b)

Rounds 6–13 take $T_e$ to be the slaved field $\Sigma(dt,du)^\top$ in the vector sector $\mathcal V$, with $\det M\neq0$. J.3 lists "components of $T_e$ outside $\mathcal V$" as not covered, because only the $\mathcal V$-block $M$ of $L_{T_o}$ was computed. This note computes $L_{T_o}$ on all of $\Lambda^1\otimes\mathfrak g_{\rm ev}$ for invariant $T_o$.

**Setting.** That of J.0: unitary gauge, the true Υ (D5), the family $\circledast=xW_A+yW'_A+zW_C\mathrm{Ad}_\Gamma+wW'_C\mathrm{Ad}_\Gamma$, the real forms (14,0) and (7,7), and $Y_a=\eta_{aa}D_a$.
- $T_o=t\Phi^1+u\Phi^1\Gamma$. By I2 the even equation reads, pointwise,
$$L_{T_o}T_e=-Y_{\rm ev}[\circledast(d_0T_o)],\qquad L_{T_o}T_e:=Y_{\rm ev}[\kappa T_e+2B(T_e,T_o)].$$
  It is algebraic in $T_e$. $L_{T_o}$ acts on the fibre $\Lambda^1\otimes\mathfrak g_{\rm ev}$, of complex dimension $14\cdot2^{13}=114688$, and depends only on the values of $(t,u)$ at the point.
- *Equivariance.* $\Phi^1$, Γ and the four operators are invariant under Spin(14) (in (7,7): Spin(7,7)); $\mathrm{Ad}_g\Gamma=\Gamma$ for even $g$. So $L_{T_o}$ commutes with $\mathfrak{so}(14,\mathbb C)$, the common complexification.
- *Notation* (J1): $\sigma=tx+\epsilon uy$, $\theta=tw-uz$, $\delta=ty-ux$, with $\epsilon=+1$ in (14,0) and $-1$ in (7,7); $\ell=xz+\epsilon yw$ (H7). Real forms as in L.0: in (14,0), $t,u,x,\dots$ stand for the real $t',u',x',\dots$ of I.0, and σ, θ, δ have the same form in them (output §R2).
- *Copies.* $\Lambda^m$ is identified with the grade-$m$ part of Cl (Chevalley map, equivariant). For $\omega\in\Lambda^m$ put $c(\omega)_a=[\gamma_a\omega]_{m-1}$ and $w(\omega)_a=[\gamma_a\omega]_{m+1}$ (the two grade parts of $\Phi^1\omega$), and let $R_\Gamma$ be right multiplication by Γ. All three are equivariant ($\mathrm{Ad}_g(X\Gamma)=\mathrm{Ad}_g(X)\Gamma$). $H_k:=[2,1^{k-1}]$ is the common kernel of $X\mapsto[\sum_a\gamma^aX_a]_{k\pm1}$ on $V\otimes\Lambda^k$. $\Lambda^7_\pm$ are the $R_\Gamma$-eigenspaces in $\Lambda^7$, eigenvalue $\pm i$ in (14,0) and $\pm1$ in (7,7) ($\Gamma^2=\mp1$).
- *Singular sets.* $Z:=\{(t,u)\in\mathbb R^2:\det L_{T_o}=0\}$; $Z_{\mathcal V}:=\{\det M=0\}$; $Z_\perp:=\bigcup_{\rho\neq\Lambda^1}\{\det B_\rho=0\}$, with the blocks $B_\rho$ of R2. So $Z=Z_{\mathcal V}\cup Z_\perp$.

## R.1 Results (layer c)

**Proposition R1 (constituents; proven).** As $\mathfrak{so}(14,\mathbb C)$-modules
$$\Lambda^1\otimes\mathfrak g_{\rm ev}=4\Lambda^1\oplus4\Lambda^3\oplus4\Lambda^5\oplus2\Lambda^7_+\oplus2\Lambda^7_-\oplus2H_2\oplus2H_4\oplus2H_6,$$
with dimensions 14, 364, 2002, 1716, 1716, 896, 11648, 36608. The copies:

| ρ | copies (grades) |
|---|---|
| $\Lambda^1$ | $c,\ w,\ wR_\Gamma,\ cR_\Gamma$ (0, 2, 12, 14): J's $P_0,P_2,P_{12},P_{14}$, with $\omega=\gamma^c$ for $v=e^c$ |
| $\Lambda^3$ | $c,\ w,\ wR_\Gamma,\ cR_\Gamma$ (2, 4, 10, 12) |
| $\Lambda^5$ | $c,\ w,\ wR_\Gamma,\ cR_\Gamma$ (4, 6, 8, 10) |
| $\Lambda^7_\pm$ | $c,\ w$ on $\Lambda^7_\pm$ (6, 8) |
| $H_k$, $k=2,4,6$ | $X,\ XR_\Gamma$ ($k$, $14-k$) |

So $\mathcal V$ is the $\Lambda^1$-isotypic component. Its complement $\mathcal V^\perp$, the sum of the other seven isotypic components, has dimension 114632.

*Proof.*
- For $\mathfrak{so}(14)$, $V\otimes\Lambda^k=\Lambda^{k+1}\oplus\Lambda^{k-1}\oplus H_k$ for $1\le k\le5$, and $V\otimes\Lambda^6=\Lambda^7_+\oplus\Lambda^7_-\oplus\Lambda^5\oplus H_6$. This is the rule $V\otimes V_\lambda=\bigoplus V_{\lambda\pm\varepsilon_i}$ over dominant weights, at $\lambda=\varepsilon_1+\dots+\varepsilon_k$. $H_k$ is irreducible with highest weight $2\varepsilon_1+\varepsilon_2+\dots+\varepsilon_k$. $R_\Gamma$ maps grade $k$ isomorphically onto grade $14-k$.
- The Weyl dimensions (output §R0) add up to $14\binom{14}k$ in each grade and to 114688 in all. Also $\dim H_k=14\binom{14}k-\binom{14}{k+1}-\binom{14}{k-1}$, so $H_k$ is the whole common kernel.
- $c$ and $w$ are nonzero equivariant maps on irreducibles, hence injective. Copies of one ρ lie in distinct grades, except on $\Lambda^7$. There $c(\omega)\Gamma=w(\omega\Gamma)$ and $w(\omega)\Gamma=c(\omega\Gamma)$, since $(\iota_a\omega)\Gamma=e_a\wedge(\omega\Gamma)$ (output §R1), so each of $\Lambda^7_\pm$ has exactly the two copies $c$, $w$. ∎

Test vectors: $\omega=E_K$ for two index sets per $m$, one inside a signature block and one across it; hook vectors $X_p=\gamma_qE_K$, $X_q=c\,\gamma_pE_K$ ($p,q\notin K$, $|K|=k-1$, $c$ fixed by $\sum_a\gamma^aX_a=0$).

**Proposition R2 (the blocks; computed, exact).** By Schur, $L_{T_o}=\bigoplus_\rho B_\rho\otimes1_\rho$, i.e. $L_{T_o}\circ\phi_i=\sum_j(B_\rho)_{ji}\phi_j$ for the copies $\phi_i$ of R1, with $B_\rho=\kappa1+tB^{(1)}_\rho+uB^{(13)}_\rho$. Every entry depends on $(t,u)$ and the member only through σ, θ, δ. In the order of R1:
$$B_{\Lambda^1}=M=\begin{pmatrix}\kappa&0&208\epsilon\delta&0\\0&\kappa+224\sigma&384i\epsilon\theta&16\epsilon\delta\\16\delta&-384i\theta&\kappa+544\sigma&0\\0&208\delta&0&\kappa\end{pmatrix}\ \ (\text{J1(b)}),$$
$$B_{\Lambda^3}=\begin{pmatrix}\kappa+\tfrac{112}3\sigma&0&\tfrac{352}3\epsilon\delta&64i\epsilon\theta\\0&\kappa+240\sigma&320i\epsilon\theta&32\epsilon\delta\\32\delta&-320i\theta&\kappa+400\sigma&0\\-64i\theta&\tfrac{352}3\delta&0&\kappa+\tfrac{272}3\sigma\end{pmatrix},\qquad B_{\Lambda^5}=\begin{pmatrix}\kappa+96\sigma&0&48\epsilon\delta&128i\epsilon\theta\\0&\kappa+\tfrac{704}3\sigma&256i\epsilon\theta&\tfrac{80}3\epsilon\delta\\\tfrac{80}3\delta&-256i\theta&\kappa+\tfrac{832}3\sigma&0\\-128i\theta&48\delta&0&\kappa+160\sigma\end{pmatrix},$$
$$B_{\Lambda^7_\pm}=\begin{pmatrix}\kappa+176\sigma&\pm192\theta\\\pm192\theta&\kappa+208\sigma\end{pmatrix}\ (14,0),\qquad B_{\Lambda^7_\pm}=\begin{pmatrix}\kappa+176\sigma&\mp192i\theta\\\mp192i\theta&\kappa+208\sigma\end{pmatrix}\ (7,7),$$
$$B_{H_k}=\begin{pmatrix}\kappa+a_k\sigma&-32i\epsilon\theta\\32i\theta&\kappa+b_k\sigma\end{pmatrix},\qquad(a_k,b_k)=\big(-\tfrac{56}3,-\tfrac{136}3\big),\ (-24,-40),\ \big(-\tfrac{88}3,-\tfrac{104}3\big)\ \ (k=2,4,6).$$
- (a) The twisted coefficients enter only through θ, and $\Phi^1\Gamma$ mixes the γ- and γΓ-channels only through δ. $\Lambda^7$ and the hooks have no δ-term. In the 4×4 blocks δ couples $c\leftrightarrow wR_\Gamma$ and $w\leftrightarrow cR_\Gamma$, and θ couples $c\leftrightarrow cR_\Gamma$ (absent for $\Lambda^1$) and $w\leftrightarrow wR_\Gamma$.
- (b) σ, θ, δ, hence every block, are invariant under M1(a) and M1(b) (output §R3). M1 was proven a symmetry of the reduced system of J1–J2; for $L_{T_o}$ it holds on all of $\Lambda^1\otimes\mathfrak g_{\rm ev}$.
- (c) $\mathrm{tr}\,L_{T_o}=114688\,\kappa$: the σ-traces cancel, $768(14+364+2002)+384\cdot3432-64(896+11648+36608)=0$.

*Verification* (output §R1, §R1b, §R1d).
- For every constituent, both test vectors, both maps ($\Phi^1$ and $\Phi^1\Gamma$) and all four operators, $L(\phi_i(v))=\sum_jB_{ji}\phi_j(v)$ with zero residual, and the two test vectors give the same matrices. The $\Lambda^1$ block equals J's $M$ from J's own `j1`.
- Equivariance spot check: $L(\delta X)=\delta L(X)$ exactly for the generators $(0,1)$ and $(2,9)$ (the latter across the signature split in (7,7)), for every constituent and all operators.
- Traces over the whole fibre, by orbit representatives of $e^a\otimes E_K$ under the signed permutations in $SO(\eta)$ (14 types in (14,0), 112 in (7,7)):
  - $\mathrm{tr}\,L_{{\rm op},\phi}=0$ for all eight pairs (operator, map), as predicted;
  - in (14,0), all 64 $\mathrm{tr}(L_aL_b)$, computed by applying $L$ twice, equal $\sum_\rho\dim\rho\,\mathrm{tr}(B_{\rho,a}B_{\rho,b})$;
  - in (7,7), $\mathrm{tr}((L-\kappa)^2)$ at one numeric member equals its prediction.

  These test the decomposition, the multiplicities and all blocks at once. ∎

**Theorem R3 (the determinant; computed).**
$$\det L_{T_o}=(\det M)^{14}\,(\det B_{\Lambda^3})^{364}\,(\det B_{\Lambda^5})^{2002}\,(\det B_{\Lambda^7})^{3432}\,(\det B_{H_2})^{896}\,(\det B_{H_4})^{11648}\,(\det B_{H_6})^{36608},$$
with $\det B_{\Lambda^7_+}=\det B_{\Lambda^7_-}=:\det B_{\Lambda^7}$ and:
- $\det M=(\rho+224\kappa\sigma)(\rho+544\kappa\sigma)-\epsilon(384\kappa\theta)^2$, $\rho=\kappa^2-3328\epsilon\delta^2$ (J1(b));
- $\det B_{\Lambda^7}=(\kappa+176\sigma)(\kappa+208\sigma)-36864\,\epsilon\theta^2$ and $\det B_{H_k}=(\kappa+a_k\sigma)(\kappa+b_k\sigma)-1024\,\epsilon\theta^2$;
- $9\det B_{\Lambda^3}=9\kappa^4+6912\kappa^3\sigma+1631744\kappa^2\sigma^2+130088960\kappa\sigma^3+2924544000\sigma^4-\epsilon(67584\kappa^2+25952256\kappa\sigma+1239941120\sigma^2)\delta^2+126877696\delta^4-\epsilon(958464\kappa^2+141557760\kappa\sigma+6658457600\sigma^2)\theta^2+1384120320\delta^2\theta^2+3774873600\theta^4$;
- $9\det B_{\Lambda^5}=9\kappa^4+6912\kappa^3\sigma+1903616\kappa^2\sigma^2+220725248\kappa\sigma^3+8996782080\sigma^4-\epsilon(23040\kappa^2+8847360\kappa\sigma+739246080\sigma^2)\delta^2+14745600\delta^4-\epsilon(737280\kappa^2+226492416\kappa\sigma+18656264192\sigma^2)\theta^2+754974720\delta^2\theta^2+9663676416\theta^4$.

Structure:
- At θ = 0 each 4×4 determinant factors as $[d_cd_{wR}-\epsilon P\delta^2][d_wd_{cR}-\epsilon P\delta^2]$, with $P=3328,\tfrac{11264}3,1280$ and $d$ the diagonal entries. At δ = 0 it factors as $[d_cd_{cR}-\epsilon q^2\theta^2][d_wd_{wR}-\epsilon r^2\theta^2]$, $(q,r)=(0,384),(64,320),(128,256)$.
- As polynomials in $(\kappa,\sigma,\theta,\delta)$ the seven determinants have real (rational) coefficients, each is irreducible over ℚ, and $\det M$ is coprime to the other six.
- $\det L_{T_o}=\kappa^{114688}$ at $T_o=0$, and every row of every block is strictly diagonally dominant when
$$544|\sigma|+384|\theta|+208|\delta|<|\kappa| ,$$
  so $L_{T_o}$ is invertible there. For $\circledast_\sharp$, $(\sigma,\theta,\delta)=(t,-u,-u)$ and the condition reads $544|t|+592|u|<|\kappa|$; for $W_A$, $(t,0,-u)$ and $544|t|+208|u|<|\kappa|$. This makes I2.2's "$|T_o|$ small compared with κ" explicit.

*Proof.* Schur and R2. The determinants and factorisations are exact sympy identities, and the dominance bound is a coefficient count (output §R2). ∎

**Proposition R4 (the source and the even equation; proven).**
- (a) For every member, $Y_{\rm ev}[\circledast(d_0T_o)]\in\mathcal V$ at every point. Indeed $d_0T_o=dt\wedge\Phi^1+du\wedge\Phi^1\Gamma$ (J1), and $v\mapsto Y_{\rm ev}[\circledast(v\wedge\Phi^1)]$, $v\mapsto Y_{\rm ev}[\circledast(v\wedge\Phi^1\Gamma)]$ are equivariant maps on the irreducible $V$. So their images lie in the $\Lambda^1$-isotypic component, which is $\mathcal V$ (R1). Exact check: J1(a) recomputed with zero residual on all 14 basis covectors, with J's $S$ (output §R1c).
- (b) At a point with values $(t,u)$, the solutions of $\Upsilon_{\rm ev}=0$ are exactly
$$T_e=T_e^{\mathcal V}+k,\qquad M\,T_e^{\mathcal V}=-S\,(dt,du)^\top,\qquad k\in K(t,u):=\bigoplus_{\rho\neq\Lambda^1}\ker B_\rho(t,u)\otimes\rho\ \subset\mathcal V^\perp .$$
  If $(t,u)\notin Z$, the unique solution is $T_e=\Sigma(dt,du)^\top$ (J1(c)). The same holds if $\det M\neq0$ and $(t,u)\notin Z_\perp$.

*Proof.* $L_{T_o}$ preserves the isotypic decomposition, and by (a) the source has no component in $\mathcal V^\perp$. ∎

**Corollary R5 (what the hypothesis "$T_e\in\mathcal V$" amounts to; proven).** Let $U$ be open, $t,u$ smooth on $U$ with $\det M(t,u)\neq0$, and $T_e$ any continuous even 1-form with $\Upsilon_{\rm ev}=0$ on $U$. Put $U_\perp:=\{q\in U:(t,u)(q)\in Z_\perp\}$.
- (a) If $U_\perp$ has empty interior, then $T_e=\Sigma(dt,du)^\top$ on $U$. This holds in particular if $(t,u)(U)\cap Z=\emptyset$, i.e. if $L_{T_o}$ is invertible on $U$.
- (b) If $U_\perp$ contains a nonempty open $W$, then $dt\wedge du=0$ on $W$. For every smooth section $k$ of $K(t,u)$ over $W$, $T_e=\Sigma(dt,du)^\top+k$ solves $\Upsilon_{\rm ev}=0$ on $W$, and nonzero such $k$ exist.

*Proof.*
- (a) $k:=T_e-\Sigma(dt,du)^\top$ is continuous, and $k(q)\in K((t,u)(q))=0$ for $q\notin U_\perp$, a dense set.
- (b) Each $\det B_\rho$ equals $\kappa^{n}\neq0$ at the origin, so $Z_\perp$ is a proper algebraic subset of $\mathbb R^2$ and has empty interior. Hence $(t,u)|_W$ has rank at most 1. $\dim K$ is upper semicontinuous, so it is locally constant on an open dense subset of $W$; there $K$ is a smooth bundle of rank at least 1. A local section times a bump function is a section over $W$. ∎

So in rounds 6–13 the hypothesis "$T_e=\Sigma(dt,du)^\top$, the slaved field in $\mathcal V$, with $\det M\neq0$ on $U$" can be replaced by

> **(R)** $\det M(t,u)\neq0$ on $U$, and $\{q\in U:(t,u)(q)\in Z_\perp\}$ has empty interior,

with $T_e$ an arbitrary smooth even 1-form. (R) is a condition on $(t,u)$ alone. A sufficient form is $(t,u)(U)\cap Z=\emptyset$.

**Proposition R6 ($\circledast_\sharp$ and $W_A$; computed).**

(a) *The t-axis.* On $u=0$ the members $\circledast_\sharp$ and $W_A$ have the same $L_{T_o}$ in both real forms ($\sigma=t$, $\theta=\delta=0$). Every block is then diagonal with entries $\kappa+at$, so $Z\cap\{u=0\}=\{-\kappa/a\}$, 18 points. In order of distance from 0 (κ = 1; the copy that spans the kernel, with its grade):
- $t<0$: $-\tfrac1{544}$ ($\Lambda^1$, $P_{12}$; $\det M=0$), $-\tfrac1{400}$ ($\Lambda^3$, $wR_\Gamma$, 10), $-\tfrac3{832}$ ($\Lambda^5$, $wR_\Gamma$, 8), $-\tfrac1{240}$ ($\Lambda^3$, $w$, 4), $-\tfrac3{704}$ ($\Lambda^5$, $w$, 6), $-\tfrac1{224}$ ($\Lambda^1$, $P_2$; $\det M=0$), $-\tfrac1{208}$ ($\Lambda^7$, $w$, 8, both ±), $-\tfrac1{176}$ ($\Lambda^7$, $c$, 6, both ±), $-\tfrac1{160}$ ($\Lambda^5$, $cR_\Gamma$, 10), $-\tfrac1{96}$ ($\Lambda^5$, $c$, 4), $-\tfrac3{272}$ ($\Lambda^3$, $cR_\Gamma$, 12), $-\tfrac3{112}$ ($\Lambda^3$, $c$, 2);
- $t>0$: $\tfrac3{136}$ ($H_2$, $XR_\Gamma$, 12), $\tfrac1{40}$ ($H_4$, $XR_\Gamma$, 10), $\tfrac3{104}$ ($H_6$, $XR_\Gamma$, 8), $\tfrac3{88}$ ($H_6$, $X$, 6), $\tfrac1{24}$ ($H_4$, $X$, 4), $\tfrac3{56}$ ($H_2$, $X$, 2).

The component of $\{u=0\}\setminus Z$ containing 0 is $-\kappa/544<t<3\kappa/136$ (κ > 0), while $\det M\neq0$ on all of $t>-\kappa/544$. On the positive t-axis the first singularity of $L_{T_o}$ is a hook block, not $M$.

(b) *Full determinants* (both real forms; output §R3).
- $W_A$:
  - $\det M=(\kappa^2+224\kappa t-3328\epsilon u^2)(\kappa^2+544\kappa t-3328\epsilon u^2)$;
  - $9\det B_{\Lambda^3}=(3\kappa^2+992\kappa t+65280t^2-11264\epsilon u^2)(3\kappa^2+1312\kappa t+44800t^2-11264\epsilon u^2)$;
  - $9\det B_{\Lambda^5}=(3\kappa^2+1120\kappa t+79872t^2-3840\epsilon u^2)(3\kappa^2+1184\kappa t+112640t^2-3840\epsilon u^2)$;
  - $\det B_{\Lambda^7}=(\kappa+176t)(\kappa+208t)$ and $\det B_{H_k}=(\kappa+a_kt)(\kappa+b_kt)$.

  So $Z_{\mathcal V}(W_A)$ is two parabolas, and $Z_\perp(W_A)$ is eight lines $t={\rm const}$ and four conics (hyperbolas in (14,0), ellipses in (7,7)).
- $\circledast_\sharp$:
  - $\det M=\kappa^4+768\kappa^3t+121856\kappa^2t^2-\epsilon(154112\kappa^2+2555904\kappa t)u^2+11075584u^4$;
  - $\det B_{\Lambda^7}=(\kappa+176t)(\kappa+208t)-36864\epsilon u^2$ and $\det B_{H_k}=(\kappa+a_kt)(\kappa+b_kt)-1024\epsilon u^2$;
  - $\det B_{\Lambda^3}$ and $\det B_{\Lambda^5}$ are those of R3 at $(\sigma,\theta,\delta)=(t,-u,-u)$.

  The hook and $\Lambda^7$ curves are hyperbolas in (14,0) and ellipses in (7,7). For example, the $H_2$ ellipse of (7,7) meets the axis at $t=3\kappa/136$ and $t=3\kappa/56$.

(c) Hence $Z_{\mathcal V}$ and $Z_\perp$ are nonempty for both members in both real forms. Exact kernel elements, checked on the Clifford algebra with a non-kernel copy as control (κ = 1; output §R3):
- (14,0), $\circledast_\sharp$:
  - at $(\tfrac3{136},0)$, $XR_\Gamma\in H_2$;
  - at $(-\tfrac1{400},0)$, $w(\omega)R_\Gamma\in\Lambda^3$;
  - at $(0,\tfrac1{32})$, $X+iXR_\Gamma$ in each of $H_2,H_4,H_6$. There the three hook blocks have rank 1 and the other four blocks are invertible, so $\ker L_{T_o}\cong H_2\oplus H_4\oplus H_6$, of dimension 49152, while $\det M\neq0$.
- (7,7): $\circledast_\sharp$ at $(\tfrac3{136},0)$, $XR_\Gamma\in H_2$; $W_A$ at $(-\tfrac1{176},0)$, $c(\omega)\in\Lambda^7$, and at $(-\tfrac3{832},0)$, $w(\omega)R_\Gamma\in\Lambda^5$.
- (14,0), $W_A$:
  - at $(-\tfrac{1831}{52400},-\tfrac{17}{524})$ on the conic $(1+\tfrac{112}3t)(1+400t)=\tfrac{11264}3u^2$ of $B_{\Lambda^3}$, the δ-coupled vector $\tfrac{25}2c(\omega)+w(\omega)R_\Gamma$;
  - at $(-\tfrac1{544},0)$, $P_{12}$, which lies in $\mathcal V$ (here $\det M=0$).

**Proposition R7 (when $Z$ is empty; proven, the $\Lambda^3$, $\Lambda^5$ step by an exact certificate).**
- (a) (14,0): for every nonzero member, each hook block and the $\Lambda^7$ block vanish somewhere. So $Z_\perp\neq\emptyset$.
- (b) (7,7): $Z=\emptyset$ iff $(x,y)=0$, or $(z,w)=c\,(y,x)$ with $c^2>\tfrac{25}{144}$. For every other member $Z_\perp\neq\emptyset$, already through the $H_2$ block. Up to M1 the second family is $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, $y(W'_A+cW_C\mathrm{Ad}_\Gamma)$ and the null members $(x,\pm x,\pm cx,cx)$, all with $|c|>\tfrac5{12}$. It lies in the cone $\ell=0$.
- (c) $\partial(\sigma,\theta)/\partial(t,u)=-\ell$. So $\ell=0$ is exactly the case where σ and θ are dependent on the $(t,u)$-plane.
- (d) For the (7,7) members of (b), $L_{T_o}$ is invertible for all real $(t,u)$. Every solution of $\Upsilon_{\rm ev}=0$ then has $\det M\neq0$ and $T_e=\Sigma(dt,du)^\top\in\mathcal V$, with no hypothesis on $T_e$. These members fall under M8 and P7 (Z2) and Q3, Q14–Q16 (Z1), whose (7,7) conclusions therefore hold for every smooth even $T_e$.

*Proof.*
- *2×2 blocks.* Write $\det=(\kappa+a\sigma)(\kappa+b\sigma)-\epsilon g^2\theta^2$, with $g=32$ (hooks) or 192 ($\Lambda^7$).
  - If σ and θ are independent ($\ell\neq0$), take $\theta=0$ and $\sigma=-\kappa/a$.
  - If $\theta\equiv0$ ($z=w=0$), then σ ≢ 0 unless the member is 0, and $\sigma=-\kappa/a$ is attained.
  - If $\sigma\equiv0$ ($(x,y)=0$), the determinant is $\kappa^2-\epsilon g^2\theta^2$: zero at $\theta=\pm\kappa/g$ in (14,0), positive in (7,7).
  - If $\sigma=\lambda\theta$ with $\lambda\neq0$: in (14,0) the discriminant in θ is $\kappa^2((a-b)^2\lambda^2+4g^2)>0$. This proves (a).
  - In (7,7), $\ell=0$ with $(x,y)\neq0\neq(z,w)$ means $(z,w)=c(y,x)$, so $\theta=c\sigma$. The block is then a binary quadratic form in $(\kappa,\sigma)$, positive definite iff $c^2>(a-b)^2/4g^2$: $\tfrac{25}{144}$ ($H_2$), $\tfrac1{16}$ ($H_4$), $\tfrac1{144}$ ($H_6$ and $\Lambda^7$). At $c^2=\tfrac{25}{144}$, $\det B_{H_2}=(\kappa-32\sigma)^2$.
- *$\Lambda^1$*, (7,7), $\theta=c\sigma$: as a quadratic in $\rho=\kappa^2+3328\delta^2$, $\det M$ has minimum $\kappa^2\sigma^2(384^2c^2-160^2)$. This is positive for $c^2>\tfrac{25}{144}$ and σ ≠ 0, and at σ = 0 the determinant is $\rho^2>0$.
- *$\Lambda^3$, $\Lambda^5$*, (7,7), $\theta=c\sigma$ (exact certificate, output §R3). Both determinants are even in θ and in δ. With $m=c^2$, $r=\delta^2\ge0$ and κ = 1, $9\det=Ar^2+B(\sigma,m)r+C(\sigma,m)$.
  - $\partial_m(9\det)=\sigma^2(\alpha m\sigma^2+\beta r+Q(\sigma))$, with $\alpha,\beta\ge0$ and $Q$ positive definite.
  - At $m=\tfrac{25}{144}$: $C\ge0$ (for $\Lambda^3$, $C=(64\sigma+1)^2(1024000\sigma^2+5760\sigma+9)$), and $4AC-B^2>0$ wherever $B<0$. So $9\det\ge0$ there.
  - Hence $9\det>0$ for $m>\tfrac{25}{144}$ and σ ≠ 0. At σ = 0 all coefficients are positive.
- *$(x,y)=0$*, (7,7): σ = δ = 0, and every determinant is a polynomial in $\kappa^2,\theta^2$ with positive coefficients.
- M1 preserves σ, θ, δ (R2(b)), which gives the orbit description in (b). (c) is a determinant (output §R2). (d) follows from (b) and R5(a). ∎

**Proposition R8 (comparison with $\det M$; proven).** $\det M$ is a factor of $\det L_{T_o}$, with exponent 14. There are six further factors, $\det B_{\Lambda^3}$, $\det B_{\Lambda^5}$, $\det B_{\Lambda^7}$ and $\det B_{H_k}$, and each is coprime to $\det M$ (R3). So $Z_{\mathcal V}\subset Z$, strictly for $\circledast_\sharp$ and $W_A$: the six positive points of R6(a), and for $\circledast_\sharp$ in (14,0) the point $(0,\kappa/32)$, lie in $Z_\perp\setminus Z_{\mathcal V}$. "$\det M\neq0$" does not imply that $L_{T_o}$ is invertible.

**Remark R9 (a flat-base configuration in $Z_\perp$; computed).** Take (14,0) and $\circledast_\sharp$ ($x'=z'=1$).
- The flat-base solutions of I7(b), N12(c) and Q15 (Z1, $\lambda=\tfrac15$) have $u=0$ and constant $s=t$ with $312s^2+\kappa s+\tfrac{15}{2c_V}=0$, so $s\in(-\kappa/312,0)$.
- That interval meets $Z$ only at two points. One is $s=-\kappa/400$: $\Lambda^3$, copy $wR_\Gamma$ (grade 10), $c_V\kappa^2=\tfrac{150000}{11}\ge9360$, $\det M\neq0$. The other is $s=-\kappa/544$: $\det M=0$, $c_V\kappa^2=\tfrac{277440}{29}$. At $c_V\kappa^2=\tfrac{150000}{11}$ the other root, $s=-\tfrac{11}{15600}\kappa$, is regular.
- So for this one value of $c_V\kappa^2$, at the root $s=-\kappa/400$ the even equation is solved by $T_e=k$ for every section $k$ of $w(\Lambda^3)\Gamma$ (grade 10, rank 364). $(t,u)$ is constant there, so (R) fails on all of $U$.
- Whether some $k\neq0$ also solves the odd equation was not examined. N12's and Q15's "$T_e\equiv0$" at that root rests on their hypothesis $T_e\in\mathcal V$.

**Answers.**
1. *$L_{T_o}$ on all of $\Lambda^1\otimes\mathfrak g_{\rm ev}$.* It has seven blocks besides the 𝒱-block $M$: 4×4 on $\Lambda^3$ and $\Lambda^5$, and 2×2 on $\Lambda^7_\pm$, $H_2$, $H_4$, $H_6$ (R1, R2).
2. *Where it is invertible.* Exactly off $Z$, with $\det L_{T_o}$ as in R3. $Z\neq\emptyset$ for $\circledast_\sharp$ and $W_A$ in both real forms (R6), for every nonzero member in (14,0), and in (7,7) for every member except $(x,y)=0$ and $(z,w)=c(y,x)$ with $c^2>\tfrac{25}{144}$ (R7).
3. *Consequences.*
   - (a) Under (R), in particular off $Z$, the even equation has the unique solution $T_e=\Sigma(dt,du)^\top\in\mathcal V$, because the source lies in $\mathcal V$ (R4, R5). The hypothesis "$T_e\in\mathcal V$" of rounds 6–13 reduces to (R), a condition on $(t,u)$.
   - (b) Over $Z_\perp$ the kernel is $\bigoplus_\rho\ker B_\rho\otimes\rho$. Generically this is one copy of one constituent ($\Lambda^3$, $\Lambda^5$, $\Lambda^7$ or a hook); at special points it has several (R6(c)). Additional solution branches of the even equation exist exactly on open sets where $(t,u)$ takes values in $Z_\perp$, and there $dt\wedge du=0$ (R5(b)). Whether they extend to solutions of $\Upsilon=0$ is open (R.3). Affected results: R.4.
4. *Comparison.* $\det M$ is a factor of $\det L_{T_o}$ (exponent 14), and there are six further factors, coprime to it (R8).

## R.2 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| R1 (§R0, §R1) | Weyl dimensions of all constituents ($D_7$); dimension count per grade and in total; copies nonzero and in the stated grades; hook vectors in the common kernel; $c(\omega)\Gamma=w(\omega\Gamma)$ on $\Lambda^7$; J's $P_k$ are the $\Lambda^1$ copies | both signatures |
| R2 (§R1, §R1b, §R1d) | exact fits with zero residual, two test vectors per constituent, both maps, four operators; templates in σ, θ, δ; $\Lambda^1$ block = J's `j1`; $\Lambda^7$ split by $R_\Gamma$ (commutation, invariance, equal determinants); equivariance for two generators; traces over the whole fibre: all eight $\mathrm{tr}\,L_{{\rm op},\phi}=0$, all 64 $\mathrm{tr}(L_aL_b)$ in (14,0), $\mathrm{tr}((L-\kappa)^2)$ at a numeric member in (7,7) | both signatures |
| R3 (§R2) | exact determinants; degree count 114688; reality; J1(b)'s form of $\det M$; the 2×2 forms; factorisations at θ = 0 and δ = 0; coprimality with $\det M$; irreducibility; dominance bound row by row; σ, θ, δ under the (14,0) real-form substitution | symbolic member; both signatures |
| R4 (§R1c) | proof (Schur); J1(a) recomputed with zero residual on all 14 $e^c$, $S$ = J's | both signatures |
| R5 | proof | both |
| R6 (§R3) | factored determinants for $\circledast_\sharp$, $W_A$; t-axis diagonal entries and the 18 points; ten exact kernel elements evaluated on the Clifford algebra with controls; ranks at $(0,\tfrac1{32})$ | both signatures |
| R7 (§R2, §R3) | proof; discriminant identity; 2×2 thresholds; $\Lambda^1$ minimum; $\Lambda^3$, $\Lambda^5$ certificate (polynomial identities, real-root isolation); $(x,y)=0$ positivity; M1 invariance; $H_2$ at the threshold; the Jacobian $-\ell$ | (14,0) (a); (7,7) (b) |
| R8 | from R3, R6 | both |
| R9 (§R4) | the I7(b) roots on $u=0$ against $Z$; both values of $c_V\kappa^2$; regularity of the other root | (14,0), $\circledast_\sharp$ |

Run: `cd computation/geometric-unity && python3 R_even_operator.py > R_even_operator.out 2>&1`. Result: 48 checks, 48 True, each asserted; 342 s on one core, a second run 353 s, with identical output apart from the timing annotations.

## R.3 Limits (layers b, d)

- **The extra branches are not analysed.** On an open set where $(t,u)\in Z_\perp$, either $(t,u)$ is constant at a point of $Z_\perp$ or it moves along an arc of $Z_\perp$. There $T_e=\Sigma(dt,du)^\top+k$ with $k$ a section of $K(t,u)$. The odd equation then gains $\circledast(d_0k)+2B(\Sigma(dt,du)^\top,k)+Q(k)$, all odd by I2, which was not computed. These branches are outside every hypothesis list of rounds 6–13, and R does not decide whether they carry solutions of $\Upsilon=0$.
- **Representation theory.** The multiplicities use the standard $\mathfrak{so}(14)$ rule for $V\otimes\Lambda^k$. They are confirmed by the dimension count and by the trace identities over the whole fibre, not re-derived from characters. Schur needs equivariance, which holds by construction and is spot-checked for two generators.
- **R7(b)** for $\Lambda^3$ and $\Lambda^5$ rests on an exact computer certificate (sympy polynomial identities and real-root isolation). The 2×2 blocks and $\Lambda^1$ are done by hand.
- **Real forms.** In (14,0), statements about real $(t,u)$ refer to the primed variables of I.0. The Z2 signatures (8,6), (10,4), … are not set up.
- **Sector.** $T_o$ is in the invariant sector only; outside it $L_{T_o}$ is not equivariant (J.3). Solutions crossing $Z_{\mathcal V}$ are not covered, as before.
- **Reconstruction.** No layer-(b) choice of J is changed.

## R.4 Consequences for gaps and claims

- **GU-G15.**
  - The slaving of I2.2 holds exactly off $Z$, where $T_e=\Sigma(dt,du)^\top\in\mathcal V$ is the unique solution of the even equation (R4).
  - "$|T_o|$ small compared with κ" is quantified: $544|\sigma|+384|\theta|+208|\delta|<|\kappa|$ suffices (R3).
  - In (14,0), $Z_\perp\neq\emptyset$ for every nonzero member, so the slaving fails somewhere for every member. In (7,7), $L_{T_o}$ is invertible everywhere for $(x,y)=0$ and for $(z,w)=c(y,x)$ with $c^2>\tfrac{25}{144}$ (R7).
- **GU-G14; the hypothesis lists of J, J5d, L, M, N, P, Q.** Replace "$T_e=\Sigma(dt,du)^\top$ the slaved field in $\mathcal V$ with $\det M\neq0$ on $U$" by (R), with $T_e$ arbitrary. Every result of rounds 6–13 then holds as stated, because under (R) every solution of $\Upsilon_{\rm ev}=0$ has $T_e=\Sigma(dt,du)^\top$ (R5(a)). In detail:
  - J1(c)'s "this is *the* slaved $T_e$" holds off $Z$, not only for small $|T_o|$.
  - J3 and J4 are pointwise at $s(X)$; their existence statements are unaffected, and the "exactly" ranges of J4(d), (d′) are exact under (R) near $s(X)$. J4(d′)'s slice $u=a_u=0$, $s\ge-\kappa/624$, contains the six positive points of R6(a); at them the slaved field is one solution of the even equation among many.
  - J5, J5d (radial Z2), L8, L9, M6–M10, P7 (non-radial Z2), N12 and Q14–Q16 (Z1) are unchanged under (R).
  - None of them covers open sets on which $(t,u)\in Z_\perp$ (R5(b)). For $W_A$ these include the radial configurations with $t\equiv3\kappa/136$ (or one of the other seven lines of R6(b)) and $u=h(s)$ arbitrary; there $T_e^{\mathcal V}=0$, since Σ has no $du$ column for $W_A$ (J1(d)), and $T_e=k$, e.g. $k\in H_2R_\Gamma$. For constant $(t,u)$ they include the flat-base root of R9.
  - For the (7,7) members of R7(b) no hypothesis on $T_e$ is needed at all (R7(d)).
- **GU-G16.** $\ell=xz\pm yw$ is, up to sign, the Jacobian of $(\sigma,\theta)$ in $(t,u)$. In (7,7), $\ell\neq0$ forces $Z_\perp\neq\emptyset$; the members with globally invertible $L_{T_o}$ lie in $\ell=0$ or $(x,y)=0$ (R7).
- **GU-054.** Status unchanged (C).
  - No committed statement is contradicted. Every statement of rounds 6–13 that involves $T_e$ carries the hypothesis "slaved in $\mathcal V$" (or $T_e\equiv0$), and J.3 lists the complement as not covered.
  - J.3's "for large $|T_o|$ there may be further branches" is answered: there are, exactly over $Z_\perp$. For $\circledast_\sharp$ and $W_A$, $Z_\perp$ already contains $t=3\kappa/136$ on the t-axis, where $\det M\neq0$ (R6).
