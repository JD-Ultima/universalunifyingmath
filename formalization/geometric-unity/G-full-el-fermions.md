# G. The complete Euler–Lagrange system, the fermion operator and the deformation complex

Claims: GU-051, 052, 053, 055, 056, 057, 072. Gaps: GU-G07, G12, G13, G15, G16, G22. Code (outputs committed):

- `computation/geometric-unity/G_el_noether.py` → `.out` (N1–N4; exact Euler operator in a non-abelian model);
- `G_fermion.py` → `.out` (F0–F6; exact, explicit real $\mathrm{Cl}(7,7)$);
- `G_symbols.py` → `.out` (Y0–Y3; exact symbol ranks in $d=14$, $\mathfrak g=\mathfrak u(64,64)$).

Lean: `proofs/lean/geometric-unity/GStueckelberg.lean` (24 theorems, no `sorry`; checked with the installed `lean` 4.24.0, core only; output `.lean.out`; `#print axioms` lists at most `propext`, `Quot.sound`). Records: `metrics/geometric-unity/G-*.yaml`. It builds on C (R-C1′, C4), D (H1–H6, D2–D5), F (Z1/Z2, F6, F8) and B (B1–B4). Layers: (a) source, (b) reconstruction, (c) proven, (d) conjectured or open, (e) testable.

## G.0 Source assertions (layer a)

Quotations are from [W21], with page = PDF page. Formulas were read from 110–260 dpi renders of pp. 44 and 46–49.

| # | Location | Text (abridged) |
|---|---|---|
| g1 | (9.1), (9.4), p. 43–44 | $I_1^B:\mathcal G\times\mathrm{MET}(X^{1,3})\to\mathbb R$, "for only the Bosonic fields of integral spin ($\varepsilon,\varpi$ on $Y$ and ℷ on $X$)" |
| g2 | (9.5)–(9.6), p. 44 | "$dI_1^B(ℷ,\omega)\vert_{Y^{14}}=(\Upsilon_\omega\oplus\Xi_\omega)\in\Omega^{d-1}(\mathrm{ad})\oplus\Omega^d(\mathrm{ad})$ where generally $\Xi=\mathfrak D_\omega\Upsilon_\omega$ for some operator differential $\mathfrak D_\omega$ so that the vanishing of $\Xi_\omega=0$ need not be considered if $\Upsilon_\omega=0$" |
| g3 | fn. 10, p. 43 | the lost Shiab operator "was chosen for its properties relative to the Bianchi identity" |
| g4 | (9.16), p. 46 | $(\bar\zeta_-\ \bar\zeta_+\ \bar\nu_-\ \bar\nu_+)_{\rho(\varepsilon)}\cdot\mathcal D_\omega\cdot(\zeta_+\ \zeta_-\ \nu_+\ \nu_-)^T_{\rho(\varepsilon^{-1})}$, with rows of $\mathcal D_\omega$: $[{*}{\circledast}\varpi_{++},\ {*}{\circledast}(d_0+\varpi_{+-}),\ \varpi_{++},\ d_0+\varpi_{+-}]$, $[{*}{\circledast}(d_0+\varpi_{-+}),\ {*}{\circledast}\varpi_{--},\ d_0+\varpi_{-+},\ \varpi_{--}]$, $[-\bar\varpi^*_{++},\ -d_0^*-\bar\varpi^*_{+-},\ 0,\ 0]$, $[-d_0^*-\bar\varpi^*_{-+},\ -\bar\varpi^*_{--},\ 0,\ 0]$; "other versions of the theory exist including one with a non-trivial map in the lower right quadrant" |
| g5 | §9.3, p. 46 | "we can at the classical level take fields $\bar\nu,\bar\zeta,\nu,\zeta$ on $Y$ to be four distinct fields when ultimately we will wish to integrate out the Dirac like operator to take a Berezinian 'integral'" |
| g6 | (9.21)–(9.22), p. 47 | $\sqrt{\Upsilon_\omega}=\delta^\omega$, "$\Upsilon_\omega=(\delta^\omega)^2=\delta_2^\omega\circ\delta_1^\omega=0$" |
| g7 | (10.1), (10.3), p. 47 | $0\to\Omega^0(\mathrm{ad})\xrightarrow{\delta_1}\Omega^1(\mathrm{ad})\oplus\Omega^0(\mathrm{ad})\xrightarrow{\delta_2}\Omega^{d-1}(\mathrm{ad})\to0$ |
| g8 | (10.6), p. 48 | $\frac d{ds}(\varepsilon,\varpi)\cdot(\exp(s\gamma),\exp(s\gamma)^{-1}d_{A_0}\exp(s\gamma))=(DL_\varepsilon\gamma,\ d_{A_0}\gamma-[\gamma,\varpi])=(DL_\varepsilon\gamma,\ d_{A_\varpi}\gamma)$ |
| g9 | (10.7)–(10.9), p. 48 | $\delta_1=(d_{A_\omega}, DL_{\varepsilon_\omega})$; $\delta_{2,a}={\circledast}_\omega\circ d_{A_\omega}(\cdot)+\kappa_1{*}(\cdot)$, $\delta_{2,b}=F_{A_\omega}\wedge{\circledast}(\cdot)-\kappa_1{*}d_{B_\omega}(\cdot)$ |
| g10 | (10.10), p. 49 | diagram "carried over from an older version and may contain some inconsistancies until it is stabilized. Caveat Emptor." |
| g11 | §12 item (x), p. 65 | "$\Upsilon=0$ carries an elliptic deformation complex in Euclidean signature once the redundant Euler-Lagrange equations are discarded" (fn. 14: Seiberg–Witten) |

**Observation (supports R-C1′).** In (10.6) the tilt is $\tau_{A_0}(e^{s\gamma})=(e^{s\gamma},\,+e^{-s\gamma}d_{A_0}e^{s\gamma})$. That is the repaired $\tau'$ of Prop. C3′, not the printed (6.4) $\tau(h)=(h,-h^{-1}d_{A_0}h)$. Its result $(DL_\varepsilon\gamma,d_{A_\varpi}\gamma)$ is the infinitesimal $\tau'$-action. Sections 6 and 10 of the draft therefore use opposite signs, and §10 uses the one that makes Lemmas 6.5 and 7.2 true.

**Not specified by the draft:**

- the $\varepsilon$-component $\Xi$ and the operator $\mathfrak D$;
- the variation in ℷ;
- ⊛ on spinor-valued forms;
- the spinor pairing and the reality and chirality of $\nu,\zeta$;
- the meaning of $\bar\varpi^*$ and of the $\rho(\varepsilon)$ placements;
- the principal symbols in §10.

## G.1 Reconstruction (layer b)

- **R-G0 (bosons).** We use R-C1′ and D.1 H1–H5. The action is
  $$I(\varepsilon,\varpi)=p\big(T,\circledast_\varepsilon(F_B+k_2d_BT+k_3T\wedge T)\big)+\tfrac\kappa2p(T,{*}T),$$
  with:
  - $p(\alpha,X)=\int b(\alpha\wedge X)$ and $b$ Ad-invariant;
  - $T=\varpi-\varepsilon^{-1}d_0\varepsilon$, $B=A_0\cdot\varepsilon$, $A=A_0+\varpi=B+T$, and $A_0=\nabla^0$;
  - $\circledast_\varepsilon=\mathrm{Ad}_\varepsilon^{-1}\circ\circledast\circ\mathrm{Ad}_\varepsilon$.

  For G1–G3, $\circledast$ is **any** pointwise map. The printed coefficients are $k_2=\tfrac12$ and $k_3=\tfrac13$, with $T\wedge T$ (G12).
- **Variations.** $\delta\varpi=\alpha$ and $\delta\varepsilon=\varepsilon\gamma$, where $\gamma\in\Omega^0(\mathrm{ad})$ is the left trivialisation, the draft's $DL_\varepsilon$. $\Upsilon$ and $\Xi$ are defined by $\delta I=p(\alpha,\Upsilon)+\int b(\gamma\,\Xi)$. So $\Upsilon\in\Omega^{d-1}(\mathrm{ad})$ and $\Xi\in\Omega^d(\mathrm{ad})$, as in g2.
- **ℷ.** It is a background for G1–G4; see G.2.5.
- **R-F (fermions).** $\psi=(\zeta,\nu)\in(\Omega^1\oplus\Omega^0)(Y,S)$, $S=\mathbb C^{128}$, with chirality given by $\Gamma$. The operator (9.16) is read as
  $$\mathcal D_\omega\begin{pmatrix}\zeta\\\nu\end{pmatrix}=\begin{pmatrix}{*}\circledast^S_\varepsilon d_A\zeta+d_A\nu\\-d_A^{*}\zeta\end{pmatrix},\qquad I_F=\int\langle\bar\zeta,{*}\circledast^S_\varepsilon d_A\zeta+d_A\nu\rangle-\langle\bar\nu,d_A^*\zeta\rangle .$$
  - Here $d_A=d_0+\varpi$ acts on spinor-valued forms by matrix multiplication ($\mathfrak u(64,64)\subset\mathrm{End}\,S$).
  - $\bar\varpi^*$ is the adjoint, so that $-d_0^*-\bar\varpi^*=-(d_A)^*$.
  - $\circledast^S_\varepsilon=\rho(\varepsilon)^{-1}\circledast^S\rho(\varepsilon)$; this is the $\rho(\varepsilon)$ sandwich.
  - $\langle\cdot,\cdot\rangle$ is the metric contraction tensored with a spinor pairing $C\in\{B,B'\}$ (real) or $H(\psi,\chi)=\psi^\dagger B\chi$ (complex, $U(64,64)$-invariant).

  The $\pm$ blocks of (9.16) are the chirality components.
- **R-F2 (⊛ on spinor-valued forms).** The natural Clifford contraction is $\circledast^S_\Phi\xi=\Phi_{11}\wedge\xi$, with $\Phi_{11}={*}\sum e^{abc}\otimes\gamma_{abc}$ acting by Clifford multiplication. This is the spinor version of $W_A$, with or without the factor $i$. The alternatives are:
  - $\circledast^S_{\Phi'}$, built with $\gamma_{abc}\Gamma$;
  - $\circledast^S_\gamma\xi={*}(e^a\otimes\gamma^b\xi_{ab})$.

  By G5 these exhaust the equivariant maps.

## G.2 The bosonic Euler–Lagrange system (layer c)

**Theorem G1 (Stueckelberg factorisation).** Assume R-G0, for arbitrary $\circledast$, Ad-invariant $b$, and $k_2,k_3,\kappa$.

1. $\tau'(\mathcal H)$ acts freely on $\mathcal G$. The fibres of $\pi'(\varepsilon,\varpi)=\varepsilon\varpi\varepsilon^{-1}-(d_0\varepsilon)\varepsilon^{-1}$ are exactly the $\tau'$-orbits. Every $g$ decomposes as $g=(1,\pi'g)\cdot\tau'(\varepsilon)$ (the "unitary gauge").
2. Pointwise, as densities, $L_I(\varepsilon,\varpi)=L_J(\pi'(\varepsilon,\varpi))$. Here
   $$L_J(S)=b\big(S\wedge\circledast(F_0+k_2d_0S+k_3S\wedge S)\big)+\tfrac\kappa2\,b(S\wedge{*}S),$$
   with $F_0=F_{\nabla^0}$, $d_0=d_{\nabla^0}$ and the *untilted* $\circledast$. Hence $I=J\circ\pi'$.

*Proof.*
- (1) This is pure group theory. Lean proves `decomp`, `fiber`, `free` and `pi'_unit` in the abstract model of E (group, module, 1-cocycle).
- (2) Conjugation by $\varepsilon$ gives $T=\varepsilon^{-1}S\varepsilon$, $F_B=\varepsilon^{-1}F_0\varepsilon$, $d_B(\varepsilon^{-1}X\varepsilon)=\varepsilon^{-1}(d_0X)\varepsilon$, $T\wedge T=\varepsilon^{-1}(S\wedge S)\varepsilon$ and $\circledast_\varepsilon(\varepsilon^{-1}Y\varepsilon)=\varepsilon^{-1}\circledast(Y)\varepsilon$. Here $*$ commutes with Ad and $b$ is Ad-invariant (L-C1).
- Equivalently: the density is $\tau'$-invariant (C4.4), and Lean `factor` gives $I(g)=I(1,\pi'g)$. ∎

*Exact check (model M3).* $d=3$, $\eta=\mathrm{diag}(1,1,-1)$, $SL(2)$, and a constant non-flat $A_0$. Three operators are used: $\circledast=\mathrm{ad}_Q^2$ (Tr-symmetric), $\mathrm{ad}_Q\mathrm{ad}_R$ (not Tr-symmetric) and $\mathrm{Id}$ (Chern–Simons).
- N1: on random integer polynomial fields, the 20 Taylor coefficients of degree ≤ 3 of $L_I$ and $L_J$ agree exactly, for 2 configurations × 3 operators.
- N4: the density is $\tau'$-invariant for a random polynomial $h$.

*Consequences.*
- The "gauge invariance bought back" by the tilt (GU-051) is of **Stueckelberg type**: $\varepsilon$ is pure gauge.
- In the unitary gauge $\varepsilon=1$ ($S=\varpi$) the theory is $J$. This is a theory of one field $S\in\Omega^1(Y,\mathrm{ad})$ on the fixed background $(Y,g_\aleph,\nabla^0)$, with the untilted ⊛.
- $J$ is not invariant under ordinary gauge transformations of $S$: $\Xi\not\equiv0$ off-shell, by G2 and N2.

**Theorem G2 (complete EL system; Noether identity).** Under the hypotheses of G1:
$$\delta I=p(\alpha,\Upsilon_\omega)+\int b(\gamma\,\Xi_\omega),\qquad \Upsilon_\omega=\mathrm{Ad}_\varepsilon^{-1}E_J(\pi'\omega),\qquad \Xi_\omega\equiv d_{A_\omega}\Upsilon_\omega\ \ \text{(off-shell identity)},$$
where $E_J$ is the $p$-gradient of $J$. The complete EL system of (9.4) in $(\varepsilon,\varpi)$ is therefore $\Upsilon_\omega=0$. The $\varepsilon$-equation is redundant.

Explicitly, $E_J(S)=\circledast F_0+k_2(L+L^\dagger)S+\kappa{*}S+k_3c(S)$, where $L=\circledast\circ d_0$ and $c$ is the gradient of $p(S,\circledast(S\wedge S))$.
- For $\circledast\in\mathcal S_{\rm sa}$ with parallel coefficients (the GU case: the $\Phi$'s are $\nabla^0$-parallel) and $k_2=\tfrac12$, $k_3=\tfrac13$, this is D5: $\Upsilon=\circledast_\varepsilon F_A+\kappa{*}T-\tfrac13h_\varepsilon(T)$.
- If $\circledast_\varepsilon$ is not $d_B$-parallel, a term $+\tfrac12\big(d_B(\circledast_\varepsilon T)-\circledast_\varepsilon(d_BT)\big)$ is added.

*Proof.*
- $\delta(\varepsilon\varpi\varepsilon^{-1})=\varepsilon[\gamma,\varpi]\varepsilon^{-1}$ and $\delta((d_0\varepsilon)\varepsilon^{-1})=\varepsilon(d_0\gamma)\varepsilon^{-1}$. So $\delta S=\varepsilon(\alpha-d_A\gamma)\varepsilon^{-1}$, and $\delta I=p(\alpha,\Upsilon)-p(d_A\gamma,\Upsilon)$.
- Then use $d\,b(\gamma\Upsilon)=b(d_A\gamma\wedge\Upsilon)+b(\gamma\,d_A\Upsilon)$ and Stokes.
- The infinitesimal $\tau'$-action is $(\varepsilon\gamma,d_A\gamma)$, which is exactly g8; along it $\delta S=0$. The formula for $E_J$ is D2/D5; the non-parallel term comes from one integration by parts. ∎

*Exact check (M3, N2/N3).* The Euler operator $E_u=\sum_N(-D)^N\partial L/\partial u_N$ is evaluated exactly at $x=0$. The tool is truncated Taylor arithmetic with 66 dual jet directions.
- $\Xi=d_A\Upsilon$ holds exactly for all three operators and both configurations, including the non-self-adjoint one.
- $\Xi\neq0$, so the identity is not vacuous.
- The general formula for $\Upsilon$ holds, including the non-parallelism term, which is non-zero for $\mathrm{ad}_Q^2$ because $[A_0,Q]\neq0$.
- For the Chern–Simons control, $h=D=0$ and $\Upsilon=\Upsilon^{\rm pr}$.
- This independently re-derives D5 by an Euler-operator method.

**Corollary G3 ((9.5)–(9.6)).** For the true $\Upsilon$ of (9.4), (9.6) holds with $\mathfrak D_\omega=d_{A_\omega}:\Omega^{d-1}(\mathrm{ad})\to\Omega^d(\mathrm{ad})$. This matches the draft's degree bookkeeping, and the redundancy claim of g2 is correct.

For the printed $\Upsilon^{\rm pr}=\circledast_\omega F_A+\kappa_1{*}T$ of (9.7), $\Xi=d_A(\Upsilon^{\rm pr}-\tfrac13h(T))\ne d_A\Upsilon^{\rm pr}$ (N3: unequal, except for the Chern–Simons control). If (9.7) is imposed as the field equation, the $\varepsilon$-equation adds the independent condition $d_{A}h_\varepsilon(T)=0$.

Since D4 shows $\Upsilon^{\rm pr}\ne\Upsilon$ for every nonzero admissible ⊛, g2 applies to D5's equation and not to (9.7). In the unitary gauge, $\Upsilon=0$ is a determined system: $14\dim\mathfrak g$ equations for $14\dim\mathfrak g$ unknowns.

**Proposition G4 (the Bianchi-type identity; fn. 10).** Let $\circledast\in\mathcal S$ (all 200 dimensions), with $(x_1,y_1,x_{13},y_{13})$ as in D3.2. For the Levi-Civita connection $\nabla$ of any metric on $Y^{14}$,
$$d_\nabla(\circledast F_\nabla)=\big(\tfrac{x_1}2+y_1\big)\,\partial_b\mathrm{Scal}\ \mathrm{vol}\otimes\gamma^b+\big(\tfrac{x_{13}}2+y_{13}\big)\,\partial_b\mathrm{Scal}\ \mathrm{vol}\otimes\gamma^b\Gamma$$
(up to the overall sign convention of $*$). So $d_\nabla(\circledast F_\nabla)\equiv0$ for all metrics iff $y_i=-x_i/2$ for $i=1,13$, the Einstein ratio in both channels.
- These are two independent linear conditions. R-S$(c)$ has $(x_1,y_1)=(-1,c/2)$, with $x_{13}=y_{13}=0$. Its $\Phi^1\Gamma$-variant at $c=0$ has $(x_{13},y_{13})=(-1,0)$.
- So the subspace $\mathcal S_B$ has dimension 198. It contains $\mathcal S_{\rm sa}$ (48, by D3.3) and R-S$(c{=}1)$.
- For $\circledast\in\mathcal S_B$ and $T=0$, the Noether identity reads $\Xi=\mathrm{Ad}_\varepsilon^{-1}d_0(\circledast F_0)\equiv0$. There (9.6) is the contracted Bianchi identity $\nabla^aG_{ab}=0$ of $g_\aleph$, and $\Upsilon=0$ is $G(g_\aleph)=0$ (F8.1).

A "property relative to the Bianchi identity" thus selects the Einstein ratio and nothing more.

*Proof.* By D3.2, $\circledast F_\nabla={*}(K_{ab}e^a\otimes\gamma^b+K'_{ab}e^a\otimes\gamma^b\Gamma)$ with $K=x_1\mathrm{Ric}+y_1\mathrm{Scal}\,g$, and similarly $K'$. The soldering forms are parallel, so $d_\nabla$ acts as the divergence. The contracted Bianchi identity $\nabla^a\mathrm{Ric}_{ab}=\tfrac12\partial_b\mathrm{Scal}$ finishes the computation, and metrics with $d\,\mathrm{Scal}\neq0$ exist. The $(x,y)$ values are those of C5 and D3. ∎

**G.2.5 The ℷ-variation (layers b, d).** The draft records only $dI\vert_{Y}$ (g2).
- Under Z1, $g_\aleph$ depends on ℷ only through $\aleph$ (F2.4). It enters $\pi'$ (via $d_0$), $F_0$, $*$, $\Phi$ and vol. So $\delta_ℷI$ is a pairing of $\partial_\Gamma g_\aleph$ with a fibre integral over $GL(4)/O(1,3)\cong\mathbb R\times SL(4)/SO(1,3)$, which is non-compact (F1).
- On solutions of $\Upsilon=0$, $T$ cannot vanish on any open set, because $\Upsilon\vert_{T=0}=-{*}G(g_\aleph)\otimes\gamma\neq0$ (F6/F8). So no compact support or decay along the fibres is available, and the X-equation is formal.
- **Open**, and not well-posed without a fibre weight.

## G.3 The fermionic operator (layers b, c)

**F0 (identification).** ${*}^{-1}\circledast^S_\Phi(\xi\wedge\zeta)_a=\gamma_a\gamma^b\gamma^c\xi_b\zeta_c$ for distinct $a,b,c$, exactly, with coefficient $+1$. So the ζζ block of (9.16) is the **Rarita–Schwinger operator** $\gamma_a{}^{bc}\nabla_b\zeta_c$, and $\mathcal D$ couples it to a spin-½ field ν through $d_A\nu$ and $d_A^*\zeta$.

**Proposition G5 (chirality).** $\mathrm{Hom}_{\mathrm{Spin}(7,7)}(\Lambda^2\otimes S^\epsilon,\Lambda^{13}\otimes S^{\epsilon'})$ is $0$ for $\epsilon'=\epsilon$. For $\epsilon'=-\epsilon$ it has dimension 2, spanned by $\circledast^S_\Phi$ and $\circledast^S_\gamma$; on $S^\epsilon$, $\circledast^S_{\Phi'}=\epsilon\circledast^S_\Phi$. So **every equivariant spinor Shiab flips chirality**, while $d_0$ (grade 2) preserves it.

Hence the labels of (9.16) are inconsistent for every equivariant $\circledast^S$:
- In row $\bar\zeta_-$, the ⊛-entry containing $d_0$ acts on $\zeta_-$ and returns $S^+$.
- The $d_0$-entry in the same row acts on $\nu_-$ and returns $S^-$.
- No invariant pairing with $\bar\zeta_-$ accepts both outputs, because $B$, $B'$, $H$ and the dual pairing each pair a fixed chirality with a fixed chirality.

**Repair R-F3:** exchange the columns $\zeta_+\leftrightarrow\zeta_-$ of the ⊛-block.

*Proof.*
- $\mathrm{Hom}(S^\epsilon,S^{\epsilon'})\cong\bigoplus\Lambda^k$, over $k\le7$ with $k$ even for $\epsilon'=\epsilon$ and odd otherwise. This follows from B1 (even monomials preserve chirality, odd ones flip it), Schur, and $64^2=4096$.
- $\Lambda^2\otimes\Lambda^1=\Lambda^3\oplus\Lambda^1\oplus\Gamma_{21}$ (as in D1). So invariants of $\Lambda^2\otimes\Lambda^1\otimes\Lambda^k$ exist only for $k\in\{1,3\}$, once each.
- F1 checks the flips and the independence of the two maps exactly. ∎

**Lemma G6 (symbol criterion for fermions).** Take a real action $\int C(\psi,\mathcal D\psi)$ with Grassmann ψ (Majorana) and parallel coefficients. Only the antisymmetric part of the kernel $K(\psi_1,\psi_2)=\int C(\psi_1,\mathcal D\psi_2)$ contributes. The EL equation is $\mathcal D\psi=0$ iff $K$ is antisymmetric, iff every symbol form $k^b$ is **symmetric**.

For the ⊛-block this means $\theta(a,b,c;\psi,\chi)=C(e^a\psi\wedge\circledast(e^b\wedge e^c\chi))_{\rm top}$ is totally antisymmetric in $(a,b,c)$ and **antisymmetric** in $(\psi,\chi)$. That is, $\circledast^S=\Psi\wedge$ with $\Psi\in\Lambda^3V\otimes\Lambda^2S^*$: this is D2 with the opposite spinor symmetry.

For complex fields with $H$, reality holds iff $K$ is Hermitian, iff $k^b$ is $H$-anti-Hermitian.

*Proof.* The $S_3\times S_2$ argument of D2 applies, with the character $\sigma\mapsto-1$, $\tau\mapsto+1$. For Grassmann variables, $\psi^TK\psi=\tfrac12\psi^T(K-K^T)\psi$, and $(\psi^\dagger K\psi)^*=\psi^\dagger K^\dagger\psi$. ∎

A check against known physics: in $d=4$, $C\gamma^{(3)}$ is antisymmetric, so the Majorana gravitino term exists, as in supergravity.

**Theorem G7 (no Majorana(–Weyl) Shiab term in (7,7)).** Take real spinor fields and $C\in\{B,B'\}$, the only invariant forms (B1.3).

1. $\Lambda^2(S^\epsilon)^*\cong\Lambda^1\oplus\Lambda^5$ and $\mathrm{Sym}^2(S^\epsilon)^*\cong\Lambda^3\oplus\Lambda^7_\epsilon$; $C$ pairs $S^+$ with $S^-$ only through even forms. Hence $\mathrm{Hom}_{\mathrm{Spin}(7,7)}(\Lambda^3,\Lambda^2S^*)=0$, and no nonzero equivariant ⊛-block is Grassmann-consistent. On the 2-dimensional space of G5 the Grassmann defect has exact rank 2 (F2).
2. The natural contraction $\circledast^S_\Phi$, and $\circledast^S_{\Phi'}$, have a $C$-kernel of the bosonic, D2 type for both $C$: exactly $\beta_b+\beta_b^t=0$ (F2). In a Majorana action the Rarita–Schwinger kinetic term of ζ therefore vanishes identically.
3. The ζν block with the printed signs $(+d_0,-d_0^*)$ is Grassmann-consistent iff $C$ is symmetric, i.e. $C=B$ (F4).
4. In the MW reading of B2–B4 ($\nu\in S^+$, $\zeta\in T^*\otimes S^-$), the action reduces to the mixing $\int B(\zeta,d_A\nu)-B(\nu,d_A^*\zeta)$, up to zeroth-order terms. Its symbol has rank 128 of 960 at **every** covector (F6: 256 of 1920 over both Weyl sectors). The $13\cdot64$ components of ζ transverse to ξ have no kinetic term, so the operator is non-elliptic in every signature.
5. The minimal real repair is two MW copies with $\epsilon_{ij}$ (symplectic Majorana–Weyl) and $C=B'$. Both blocks are then consistent: (sym ⊗ antisym) and (antisym ⊗ antisym). This is equivalent to one complex Weyl field with $H$, and the content doubles to $6\times16_+\oplus2\times144_-$ (the complex-Weyl row of B4).

*Proof.*
- (1) F3 gives the symmetry of $C\gamma_I\vert_{S^\epsilon\times S^\epsilon}$ for every monomial of grades 1, 3, 5, 7: antisymmetric, symmetric, antisymmetric, symmetric, for both $C$. Even grades vanish there.
- The equivariant maps $\Lambda^k\to\mathrm{Bil}(S^\epsilon)$, $e^I\mapsto C\gamma_I\vert$, are injective on the pairwise non-isomorphic irreducibles $\Lambda^1,\Lambda^3,\Lambda^5,\Lambda^7_\epsilon$ (Schur).
- Their dimensions are $14+364+2002+1716=4096=\dim\mathrm{Bil}$. The split is $364+1716=2080=\dim\mathrm{Sym}^2$ and $14+2002=2016=\dim\Lambda^2$.
- (2)–(4) are exact (F2, F4, F6) together with G6. (5) follows from the symmetry of tensor products. ∎

*Remark (alternative, layer c; elementary count, not machine-checked).* For the chimeric group $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$, Weyl's first fundamental theorem gives $\dim\mathrm{Hom}(\Lambda^3V,\Lambda^1V\oplus\Lambda^5V)=1+1$. Both invariants use the 4-volume form $\varepsilon_X$. Majorana-consistent ζζ terms therefore exist only if Spin(7,7)-covariance is given up and $X$'s orientation is used. Whether they give a non-degenerate kinetic term is open.

**Proposition G8 (Hermitian readings).** Take complex fields and $H(\psi,\chi)=\psi^\dagger B\chi$.
1. Build $\circledast^S_\Phi$ from the $\mathfrak u(64,64)$-elements $i\gamma_{abc}$, as in $W_A$, and take the ζν block with the printed signs. Both blocks are then anti-Hermitian (F5, exact), so $\int H(\psi,i\mathcal D\psi)$ is real and the printed relative phases are consistent. With real $\gamma_{abc}$ the ⊛-block is Hermitian, and a relative factor $i$ is needed.
2. Complex Weyl ($\nu\in S^+_{\mathbb C}$, $\zeta\in T^*\otimes S^-_{\mathbb C}$) is consistent with G5: $H$ pairs $S^+$ with $S^-$, ⊛ flips chirality and $d_0$ preserves it.
   - The chirality-changing components of ϖ drop out. Only $\mathfrak u(64,64)_{\rm ev}\cong\mathfrak{gl}(64,\mathbb C)$ couples, i.e. 8192 of 16384 generators.
   - The content is the complex-Weyl row of B4.
3. The Dirac reading is consistent and vector-like (B4).

*Proof.*
- (1) F5.
- (2) Use $H(S^\epsilon,S^\epsilon)=0$ ($\Gamma$ is real) and chirality bookkeeping.
- For $\mathfrak u(64,64)_{\rm ev}\cong\mathfrak{gl}(64,\mathbb C)$: an element $X=\mathrm{diag}(X_+,X_-)$ satisfies $X^\dagger B+BX=0$ with $B$ off-diagonal in chirality iff $X_-$ is determined by $X_+$, and $X_+$ is free. ∎

**Theorem G9 (fermion symbol).** Consider $\mathcal D$ with $\circledast^S_\Phi$, in the Dirac or complex-Weyl reading.
- $\sigma(\mathcal D)(\xi)$ is invertible iff $\eta(\xi,\xi)\ne0$.
- In (7,7) the rank is $1920/1920$ (Weyl sector $960/960$) at $e^0$, $e^7$ and $e^0+e^1$.
- At the null covectors $e^0+e^7$ and $e^0+e^1+e^7+e^8$ the rank is $1024$, a deficit of $896=7\cdot128$; per Weyl sector it is $512/960$.
- In (14,0) the rank is $1920/1920$ at every tested ξ, and at all ξ ≠ 0 by the proof.

So the fermion operator is **elliptic in Euclidean signature**. In (7,7) it is singular exactly on the null cone.

*Proof.*
- For non-null ξ, choose an orthonormal coframe with $e^0\parallel\xi$. Then $\sigma(\zeta,\nu)=0$ forces $\nu=0$, $\zeta_0=0$ and $\zeta_a=\gamma_a(\gamma\cdot\zeta)$, hence $\gamma\cdot\zeta=13\,\gamma\cdot\zeta$, so $\zeta=0$. The system is square.
- For null ξ, $(\xi\otimes\psi,0)$ lies in the kernel. The exact deficit is F6. ∎

## G.4 The deformation complex (layer c)

**Proposition G10 (structure).** Under R-C1′:
1. $\delta_1\gamma=(d_{A_\omega}\gamma,\gamma)$ (g9, left-trivialised). $\delta_2$, the linearisation of $\Upsilon$ in $(\varpi,\varepsilon)$, is
   $$\delta_2(\alpha,\gamma')=\mathcal L_\omega(\alpha-d_A\gamma')+[\Upsilon_\omega,\gamma'],\qquad \mathcal L_\omega=\circledast_\varepsilon d_A+\kappa{*}-\tfrac13Dh_{\varepsilon,T}.$$
   For $\Upsilon^{\rm pr}$, $\delta_2$ equals the draft's (10.9) term by term:
   - $\delta_{2,a}=\circledast d_A+\kappa_1{*}$;
   - $\delta_{2,b}\gamma'=\circledast_\varepsilon[\gamma',F_A]-[\gamma',\circledast_\varepsilon F_A]-\kappa_1{*}d_B\gamma'$. The first two terms are the variation of the tilt; reading them as the draft's "$F_A\wedge\circledast(\cdot)$" is our interpretation of that symbol. The $-\kappa_1{*}d_B$ term matches the draft exactly.

   In both cases $\delta_2\circ\delta_1=\mathrm{ad}(\Upsilon_\omega)$, which is zeroth order. So (9.22)/(10.5) "$\Upsilon=\delta_2\circ\delta_1$" holds in this sense, and (10.1) is a complex exactly on-shell.
2. *Stueckelberg elimination.* Substitute $\beta=\alpha-d_A\gamma'$. At $\Upsilon=0$ the complex then becomes $(\Omega^0\xrightarrow{\rm id}\Omega^0)\oplus(\Omega^1\xrightarrow{\mathcal L}\Omega^{13})$, so formally $H^0=0$, $H^1\cong\ker\mathcal L$ and $H^2\cong\mathrm{coker}\,\mathcal L$.
3. *Symbols (Atiyah–Bott, uniform order 1).* Since $d_Ad_A=[F_A,\cdot]$, $\mathcal L(d_A\gamma')$ is first order. Then
   - $\sigma\delta_1(\xi)\gamma=(\xi\otimes\gamma,0)$;
   - $\sigma\delta_2(\xi)(\alpha,\gamma')=\circledast_\varepsilon(\xi\wedge\alpha)-M_T(\xi\otimes\gamma')$, with $M_T=\kappa{*}-\tfrac13Dh_T$ for the true $\Upsilon$ and $\kappa{*}$ for the printed one;
   - $\sigma\delta_2\circ\sigma\delta_1=0$, and $\sigma\delta_1$ is injective;
   - the ranks give $1-15+14=0$.

   So the symbol sequence is exact iff $\mathrm{rank}\,\sigma\delta_2(\xi)=14\dim\mathfrak g$.

*Proof.*
- Both Υ's are $\tau'$-covariant: G1 for the true one, C4 for the printed one. Differentiate $\Upsilon_{\omega\tau'(h)}=h^{-1}\Upsilon_\omega h$ along $(d_A\gamma,\gamma)$.
- The (10.9) identification is a direct computation with $\delta T=-d_B\gamma'$ and $\delta(\circledast_\varepsilon X)=\circledast_\varepsilon[\gamma',X]-[\gamma',\circledast_\varepsilon X]$. ∎

**Theorem G11 (the deformation complex is not elliptic).**

1. *Null-cone obstruction.* Take any $\circledast\in\mathcal S_{\rm sa}$ (48-dimensional), any κ and any background $(\varepsilon,T)$, with the true or the printed Υ. For every null $\xi\ne0$, $\mathrm{im}\,\sigma\delta_2(\xi)\subset\mathrm{Ann}(\xi\otimes\mathfrak g)$, which has codimension $\dim\mathfrak g$. So **in signature (7,7) the deformation complex of (9.4) is never elliptic.**
2. *Exact ranks* ($d=14$, $\mathfrak u(64,64)$, $\dim\Lambda^1\otimes\mathfrak g=229376$; Y1–Y3). The entry "rank ⊛(ξ∧·)" is the rank of the symbol of $\circledast\circ d$; ellipticity modulo the gauge directions would need $229376-16384=212992$.

| ⊛ (signature, algebra) | rank ⊛(ξ∧·): non-null / null ξ | rank σδ₂ (κ=1), non-null ξ | rank σδ₂ (κ=1), null ξ |
|---|---|---|---|
| $W_A$, $W'_A$, $W_A+2W'_A$ — (7,7), $\mathfrak u$ | 122850 / 122850 | 139234 (deficit 90142) | 126960 (deficit 102416) |
| $W_A\pm W'_A$ — (7,7), $\mathfrak u$ | 106496 / 86016 | 122880 (deficit 106496) | 95232 (deficit 134144) |
| R-S $(C;A,A)$, $c=1$ — (7,7), $\mathfrak u$ | 122864 / 122864 | 131057 (deficit 98319) | 126974 (deficit 102402) |
| $W_A$, $W'_A$, $W_A+W'_A$, $W_A+2W'_A$ — (14,0), $\mathfrak u$ | 122850 / — | 139234 (deficit 90142) | — (no real null ξ) |
| R-S — (14,0), $\mathfrak u$ | 122877 / — | 131070 (deficit 98306) | — |
| $W_A$, $W'_A$ — (7,7), $\mathfrak{so}(64,64)$ (of 113792; need 105664) | 60580 / 60580 | 68708 (deficit 45084) | 62609 (deficit 51183) |
| $W_A+W'_A$ — (7,7), $\mathfrak{so}(64,64)$ | 52416 / 42560 | 60544 (deficit 53248) | 47152 (deficit 66640) |
| $W_A$, $W'_A$, $W_A+W'_A$ — (14,0), $\mathfrak{so}(64,64)$ | 60580 / — | 68708 (deficit 45084) | — |

   Tested covectors: in (7,7), non-null $e^0$, $e^7$, $e^0+e^1$ and null $e^0+e^7$, $e^0+e^1+e^7+e^8$; in (14,0), $e^0$ and $e^0+e^1$ for $\mathfrak u$, and $e^0$ and $e^0+e^7$ for $\mathfrak{so}$. All listed operators and $\kappa{*}$ are Spin(7,7)-equivariant, and the tilt acts at a point by an invertible conjugation. $SO_0(7,7)$ is transitive on each of $\{\eta(\xi,\xi)=c\}$, $c\ne0$, and on the nonzero null cone, and $SO(14)$ is transitive on spheres. The ranks therefore hold at **every** nonzero covector of each class; the equal ranks at $e^0$ and $e^0+e^1$ confirm this.
   Moreover the γ′-columns add at most $\dim\mathfrak g$ to the rank, so $\mathrm{rank}\,\sigma\delta_2(\xi)\le122877+16384<229376$. This holds for **every background** $(\varepsilon,T)$, every κ, and the true or printed Υ, since $M_T$ only enters those columns. The symbol sequence is **exact at no covector, in either signature or algebra**.
3. *The whole pencil.* On a chirality block $\mathrm{Hom}(S^\tau,S^\sigma)$, write $L=\gamma_J\cdot$ and $R=\cdot\gamma_J$ for the two halves of $W_A$'s anticommutator. Then $W'_A=\sigma L-\tau R$, so $xW_A+yW'_A=(x+y\sigma)L+(x-y\tau)R$.
   - The block graph (source $(\tau,\sigma)\to$ targets $(\tau,-\sigma)$ via $L$ and $(-\tau,\sigma)$ via $R$) consists of two 4-cycles with zero holonomy.
   - So for $x+y\sigma\ne0\ne x-y\tau$, invertible block rescalings of source and target conjugate $xW_A+yW'_A$ to $W_A$. Rescaling the γ′ columns by the target factor makes this compatible with the κ-term.
   - Hence every member of $\mathrm{span}(W_A,W'_A)$ has the ranks of $W_A$, except, in (7,7) ($\sigma,\tau=\pm1$), the two lines $W_A\pm W'_A$, which are computed separately.
   - In (14,0), $\Gamma^2=-1$ and $\sigma,\tau=\pm i$, so every real member has $W_A$'s ranks. The table confirms this for $W_A+W'_A$, which in (7,7) is degenerate.

   So the conclusion covers **every operator of the draft's self-adjoint family (8.1)** ($\mathrm{span}(W_A,W'_A)$, D3.5), and R-S, at every covector and every background.

   The source of the degeneracy: $W_A$'s coefficient is an anticommutator with grade-3 Clifford elements, which vanishes on every monomial that anticommutes with them. For example, $W_A(\xi\wedge e^a\otimes\Gamma)=0$ for all $a$.

*Proof.*
- (1) $\circledast_\varepsilon=\Psi_\varepsilon\wedge$ (D2), so $b(\xi\gamma\wedge\circledast_\varepsilon(\xi\wedge\alpha))=\Psi_\varepsilon(\xi\gamma,\xi,\alpha)=0$.
- $b(\xi\gamma\wedge{*}(\xi\gamma'))=\eta(\xi,\xi)\,b(\gamma,\gamma')=0$.
- Every term of $p(\xi\gamma,Dh_T(\xi\gamma'))$ evaluates $\Psi$ with ξ in two of its three antisymmetric form slots.
- (2) The computation is exact: monomial Clifford algebra, connected blocks, integer elimination. (3) is the rescaling argument above, checked by the equal ranks of $W_A$, $W'_A$ and $W_A+2W'_A$ in (7,7) and of all real members in (14,0).
- It is cross-checked against the explicit $128\times128$ matrices: products; $W_A$, $W'_A$ and R-S on random inputs; assembled columns; and an independent matrix-route rank on a 210-column subspace, which agrees at 195. ∎

So the claim of g11, "an elliptic deformation complex in Euclidean signature", **fails for the draft's self-adjoint family and for R-S**. In (7,7) it fails for all of $\mathcal S_{\rm sa}$. The fermionic block is elliptic in Euclidean signature (G9).

## G.5 What is verified, precisely

| Statement | Verified how | Scope |
|---|---|---|
| G1 (1) | Lean, abstract model | all groups, modules, cocycles |
| G1 (2), G2, G3 | proof; exact jet / Euler-operator evaluation | proof: any $d$, ⊛, $b$; check: model M3 ($d=3$, $SL(2)$, non-flat $A_0$, 3 operators, 2 random configurations, at one point) |
| G4 | proof from D3.2 and contracted Bianchi | all of $\mathcal S$ |
| G5–G9 | proofs; exact integer arithmetic in the $\mathrm{Cl}(7,7)$ model | ranks at the listed covectors; symmetry types on all monomials; Grassmann defects on random integer data (any nonzero defect is a disproof; zero defects are backed by the symmetry-type proofs) |
| G10 | proof | — |
| G11 (1) | proof | all of $\mathcal S_{\rm sa}$, all κ, $T$ |
| G11 (2)–(3) | exact ranks at orbit representatives; equivariance; rescaling lemma; rank bound for the γ′ columns | span$(W_A,W'_A)$ and R-S; all covectors and backgrounds; both signatures; $\mathfrak u(64,64)$ and $\mathfrak{so}(64,64)$; not the remaining 46 dimensions of $\mathcal S_{\rm sa}$ |

## G.6 Grading and claims

| Item | Status |
|---|---|
| G1, G2, G3, G4, G5, G6, G7, G8, G9 (non-null part by proof, null deficit computed), G10, G11 (1) | **Proven** |
| G11 (2)–(3) | **Proven** (exact computation, equivariance, rescaling lemma) for all of $\mathrm{span}(W_A,W'_A)$ and for R-S, at every covector and background |
| Chimeric Majorana remark | Proven count, not machine-checked; non-degeneracy open |
| ℷ-equation | Open; formal only (G.2.5) |
| Euclidean ellipticity for the other 46 dimensions of $\mathcal S_{\rm sa}$ | Open |

Claims:
- **GU-051** C: true for $\tau'$, and the symmetry is of Stueckelberg type (G1).
- **GU-052** C, from O: (9.5)–(9.6) hold with $\mathfrak D=d_{A_\omega}$ for the true Υ, but not for the printed (9.7).
- **GU-053** F, unchanged. D5's equation is now known to be the complete EL system.
- **GU-055** unchanged.
- **GU-056** C/F, from O: well-defined under R-F/R-F3. The printed chirality labels are inconsistent (G5). In Majorana(–Weyl) form the Shiab block vanishes (G7).
- **GU-057** split. "$\Upsilon=\delta_2\circ\delta_1$" holds as $\mathrm{ad}(\Upsilon)$ (V). Ellipticity is F (G11).
- **GU-072** bosonic part F as for GU-057. The fermionic part is V in Euclidean signature, conditional on the Hermitian reading (G9).

Gaps:
- **G15** closed for $(\varepsilon,\varpi)$; the ℷ-part stays open.
- **G16** closed: negative for the specified operators, and in (7,7) for all of $\mathcal S_{\rm sa}$.
- **G22** closed under R-F, with the inconsistencies above.
- **G07** refined: the MW reading is incompatible with the ζ kinetic term (G7). See `reports/CORRECTIONS.md`.

**Empirical (layer e).** Nothing new is testable. G7 bears on the generation count of B4: with a Rarita–Schwinger kinetic term for ζ, the content is the doubled, complex-Weyl one.


## Independent cross-check (coordinator)

`computation/geometric-unity/xcheck_mw_kinetic.py` uses its own real construction of $\mathrm{Cl}(7,7)$ on $\mathbb R^{128}$, built as a graded tensor power of $\mathrm{Cl}(1,1)$; it does not use `clifford77.py`. It finds:

- the two Spin-invariant bilinears are $C_+=\gamma_1\cdots\gamma_7$, which is antisymmetric, and $C_-=\gamma_8\cdots\gamma_{14}$, which is symmetric;
- both pair $S^+$ with $S^-$ and vanish on $S^+\times S^+$;
- every Dirac kinetic block $C\gamma^\mu|_{S^+\times S^+}$ is exactly antisymmetric.

For a Grassmann field, $\int\psi^TM\partial\psi$ only sees the symmetric part of $M$, so a single Majorana–Weyl spinor has no Dirac kinetic term in $(7,7)$. It also finds that the full Rarita–Schwinger kinetic matrix $(C\gamma^{\mu\nu\rho})$ on $(T^*\otimes S^-)^{\otimes2}$, including its vector-index antisymmetry, has exactly zero symmetric part for both $C_\pm$. So a single Majorana–Weyl $\zeta$ also has no Rarita–Schwinger kinetic term. This supports Theorem G7 and the doubling to complex Weyl.
