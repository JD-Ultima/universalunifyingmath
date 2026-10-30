# D. The space of Shiab operators: self-adjointness, the Einstein property and whether (9.7) is variational

Claims: GU-045, GU-047, GU-050, GU-053. Gaps: GU-G10, GU-G13. Pattern: GU-P09. Code: `computation/geometric-unity/D_rep_counts.py` → `.out` (exact integers; brute-force ranks by SVD with the gap printed) and `D_selfadjoint_shiab.py` → `.out` (exact: integer and dyadic data, `==` tests, magnitudes asserted below $2^{53}$). Metric records: `metrics/geometric-unity/D-*.yaml`. Builds on C.1–C.4 (`C-gauge-shiab.md`).

## D.0 Source assertions (layer a)

Def. 8.1: $\{\Phi_i\}_{i=0}^{14}$ is "a basis for the invariant subspaces of $[\Lambda^i(\mathbb R^{7,7})\otimes\mathfrak u(64,64)]^{\mathrm{Spin}(7,7)}$". (8.1): a Shiab operator "might look like" $*^n i^{\frac{1\pm1}2}[(\varepsilon^{-1}\Phi_r\varepsilon)\wedge*^m\eta]_\pm$ "or perhaps a sum of such terms". (9.3): the chosen operator kills Weyl and is "Ricci like − ½ $g$-like · scalar like". (9.4) is the action and (9.7) its claimed Euler–Lagrange equation $\Upsilon=\circledast_\omega F_{A_\omega}+\kappa_1*T_\omega$. Fn. 10: the original operator was "chosen for its properties relative to the Bianchi identity" and is lost.

## D.1 Hypotheses and reconstruction choices (layer b)

- **H1 (structure group).** $K=\mathrm{Spin}(7,7)$ acts on $V=T_yY\cong\mathbb R^{7,7}$ through $SO(7,7)$, and on $\mathfrak g$ by conjugation through the spin representation in $\mathrm{Cl}(7,7)\otimes\mathbb C=M_{128}(\mathbb C)$.
  - *Justification.* Def. 8.1 takes $\mathrm{Spin}(7,7)$-invariants. A pointwise map built from $g_\aleph$ and the spin structure alone must be $K$-equivariant, and such a map is parallel for the Levi-Civita connection, as Lemma C6a requires.
  - *Alternatives.* Theorem D4 also covers the chimeric reduction $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$, and D4(a) holds with no equivariance at all. The dimension counts are for $K$ only.
- **H2 (gauge algebra).** The draft allows $\mathfrak u(64,64)$ (Def. 8.1) and $\mathfrak{so}(64,64)$ ((8.4)–(8.5)). As test sectors we also use two Clifford-graded subalgebras: $\mathfrak{so}(8,8)=$ grades $\{1,2,13,14\}$ and $\mathfrak{so}(7,8)=$ grades $\{1,2\}$. They are $\Lambda^2(V\oplus\mathbb R^{1,1})$ and $\Lambda^2(V\oplus\mathbb R^{0,1})$, the smallest graded algebras containing $\mathfrak{spin}(7,7)$ and a copy of $V$. We also treat $\mathfrak{so}(7,7)$ itself.
  - As a real $K$-module, $\mathfrak g\cong\bigoplus_{k\in\mathrm{grades}}\Lambda^kV$. The quantisation map $\Lambda^\bullet V\to\mathrm{Cl}(V)$ is $O(V)$-equivariant, and the factors $i$ in (8.6) do not change the real module.
  - In signature $(7,7)$, $*^2=+1$ on $\Lambda^7$. So $\Lambda^7=\Lambda^7_+\oplus\Lambda^7_-$ over $\mathbb R$, and all constituents are absolutely irreducible. Real dimensions below therefore equal complex dimensions for $D_7=\mathfrak{so}(14,\mathbb C)$, and every result also holds for the complexification $\mathfrak{gl}(128,\mathbb C)$.
- **H3 (operator class).**
  $$\mathcal S:=\mathrm{Hom}_K(\Lambda^2V^*\otimes\mathfrak g,\Lambda^{13}V^*\otimes\mathfrak g),$$
  real-linear maps, tilted as $\circledast_\varepsilon=\mathrm{Ad}_\varepsilon^{-1}\circ\circledast\circ\mathrm{Ad}_\varepsilon$ (C4.2). This class contains (8.1), all sums of such terms, and R-S. Following (9.3) and the integral line of (9.4), $\circledast_\omega:=\circledast_\varepsilon$ depends on $\varepsilon$ only.
- **H4 (pairing).** $p(T,X)=\int b(T\wedge X)$ with $b$ a nondegenerate symmetric bilinear form on $\mathfrak g$ that is $\mathrm{Ad}$-invariant. For $\mathfrak u(64,64)$ this means $b=a\,\mathrm{Tr}(XY)+c\,\mathrm{Tr}X\,\mathrm{Tr}Y$; the printed choice is $\mathrm{Tr}$.
  - Writing $b=\mathrm{Tr}(\cdot\,M\cdot)$, the map $\circledast\mapsto M\circledast$ is a bijection of $\mathcal S$ that preserves the Einstein output (D3). So D2, D3 and D5 do not depend on $b$.
  - D4 holds for any $K$-invariant $b$.
  - Reading $\langle T,*\circledast(\cdot)\rangle_g$ instead of $\int\mathrm{Tr}(T\wedge\circledast(\cdot))$ changes $p$ by a constant sign.
- **H5 (variation).** Vary $\varpi$ with $\varepsilon$ fixed, as in (9.7). Set $T=\varpi-\varepsilon^{-1}d_0\varepsilon$ and $F_{A}=F_B+d_BT+T\wedge T$ (where $T\wedge T=\tfrac12[T\wedge T]$). The coefficients $k_2$ of $d_BT$ and $k_3$ of the cubic term in (9.4) are left free.
- **H6 (Einstein-type).** $\circledast(F_R)=*\bigl((\mathrm{Ric}-\lambda\,\mathrm{Scal}\,g)_{ab}\,e^a\otimes\iota(e^b)\bigr)$ for all algebraic curvature tensors $R$, with $\iota:V\to\mathfrak g$ nonzero and $K$-equivariant. *Einstein* means $\lambda=\tfrac12$. Here $F_R=\sum_{a<b}e^{ab}\otimes\tfrac14R_{abcd}\gamma^c\gamma^d$, as in Thm C5.

Definitions:
- $\Phi^3=\sum_{a<b<c}e^{abc}\otimes i\gamma_a\gamma_b\gamma_c$ and $\Phi^{3\prime}=\sum e^{abc}\otimes i\gamma_{abc}\Gamma$, with $\Gamma=\gamma_0\cdots\gamma_{13}$.
- $\Phi_{11}=*\Phi^3$ and $\Phi_{11}'=*\Phi^{3\prime}$.
- $W_s\xi=[\Phi_{11}\wedge\xi]_s$ and $W'_s\xi=[\Phi'_{11}\wedge\xi]_s$ for $s\in\{A,C\}$.

$\Phi^3$ is the antisymmetrised cube of the invariant $\Phi^1=e^m\otimes\gamma_m$, so it is $K$-invariant. $W_s$ is the member $r=11$, $m=n=0$ of (8.1).

## D.2 Results (layer c)

**Proposition D1 (size of the space).** Exact multiplicities, from the Weyl character formula for $D_7$; every decomposition is checked against the Weyl dimension formula (D_rep_counts D1–D3):

| $\mathfrak g$ | $\dim\mathfrak g$ | $\dim\mathcal S$ | $\dim\mathcal S_{\rm sa}=\dim\mathrm{Hom}_K(\Lambda^3,S^2\mathfrak g)$ | $\dim\mathcal S_{\rm skew}=\dim\mathrm{Hom}_K(\Lambda^3,\Lambda^2\mathfrak g)$ |
|---|---|---|---|---|
| $\mathfrak u(64,64)$ | 16384 | **200** | **48** | 48 |
| $\mathfrak{so}(64,64)$ | 8128 | 52 | 12 | 12 |
| $\mathfrak{so}(8,8)$ | 120 | 16 | 2 | 2 |
| $\mathfrak{so}(7,8)$ | 105 | 6 | 1 | 1 |
| $\mathfrak{so}(7,7)$ | 91 | 0 | 0 | 0 |

Two further facts:
- $\dim[\Lambda^i\otimes\mathfrak u(64,64)]^K=2$ for every $i=0,\dots,14$, spanned by $\Phi_i$ and $\Phi_i\Gamma$. Def. 8.1's "basis $\{\Phi_i\}_{i=0}^{14}$" therefore names 15 of the 30 invariants.
- For $\mathfrak g=\mathfrak{so}(7,7)$ no Shiab operator exists, because $\Lambda^2\otimes\Lambda^2$ and $\Lambda^1\otimes\Lambda^2$ share no constituent.

*Checks.*
- The identity $\dim\mathcal S=\dim\mathcal S_{\rm sa}+\dim\mathcal S_{\rm skew}+\dim\mathrm{Hom}_K(V\oplus\Gamma_{21},\mathfrak g\otimes\mathfrak g)$, which follows from $V\otimes\Lambda^2=\Lambda^3\oplus V\oplus\Gamma_{21}$, holds in every row.
- For $D_2$ and $D_3$ ($n=4,6$, $\mathfrak g=\Lambda^\bullet$), brute force gives $\dim\mathcal S=40,72$ and $\dim\mathcal S_{\rm sa}=8,16$. Brute force here means zero-weight vectors killed by the simple raising operators, with the symmetry conditions of D2 imposed as linear constraints. All 86 blocks agree with the characters. The smallest kept singular value is $0.73$; the largest dropped is $1.2\times10^{-15}$. This tests D2 itself, not only the character code. ∎

**Theorem D2 (symbol criterion; any $K$, any $b$).** For a pointwise map $\circledast:\Lambda^2\otimes\mathfrak g\to\Lambda^{13}\otimes\mathfrak g$, set
$$\theta(a,b,c;X,Y):=b\bigl(e^a X\wedge\circledast(e^b\wedge e^c\,Y)\bigr)_{\rm top}.$$
The principal symbol of $\circledast\circ d_B$ is formally $p$-self-adjoint, i.e. $\beta_b(T,\alpha)=-\beta_b(\alpha,T)$ for all $b$ (Lemma C6a), **iff** $\theta$ is totally antisymmetric in $(a,b,c)$ and symmetric in $(X,Y)$. Equivalently, there is $\Psi\in\Lambda^3V\otimes S^2\mathfrak g^*$ with $b(\alpha\wedge\circledast\xi)=\Psi(\alpha\wedge\xi)$; in short, "$\circledast=\Psi\wedge$". For parallel $\circledast$ (H1), this is equivalent to formal self-adjointness of $\circledast\circ d_B$. Hence $\mathcal S_{\rm sa}\cong\mathrm{Hom}_K(\Lambda^3V,S^2\mathfrak g)$.

*Proof.* Always $\theta(a,c,b;X,Y)=-\theta(a,b,c;X,Y)$; call this $\sigma$, since $\theta$ depends on $e^b\wedge e^c$. Self-adjointness is $\theta(c,b,a;Y,X)=-\theta(a,b,c;X,Y)$; call this $\tau$.
- $\tau\sigma\tau$ is the transposition $(ab)$ of the form slots alone, so $\langle\sigma,\tau\rangle\supset S_3$ acting on $(a,b,c)$ with $(X,Y)$ fixed.
- $\tau\circ(ac)$ is the swap $X\leftrightarrow Y$ alone, so $\langle\sigma,\tau\rangle=S_3\times S_2$.
- The character with $\sigma,\tau\mapsto-1$ is $\mathrm{sgn}\otimes1$.

The converse direction and the parallel case follow from Lemma C6a and the nondegeneracy of $p$. ∎

Sanity check: Chern–Simons ($d=3$, $\circledast=\mathrm{Id}$) has $\theta=\epsilon_{abc}\mathrm{Tr}(XY)$.

**Theorem D3 (Weyl and Einstein; $K=\mathrm{Spin}(7,7)$).**
1. **Every** $\circledast\in\mathcal S$ annihilates Weyl curvature. *Proof:* By Schur's lemma, since $\Lambda^1\otimes\Lambda^l$ contains no Weyl-type constituent $(2,2,0,\dots)$ for any $l$ (D_rep_counts D4). This is confirmed exactly on $312\,W$ for $W_A$, $W'_A$, $W_C$, $W'_C$ and all six R-S terms (E2). The draft's "the Weyl tensor is annihilated" thus holds on all 200 dimensions and does not select an operator.
2. For every $\circledast\in\mathcal S$ there are $x_1,y_1,x_{13},y_{13}$, linear in $\circledast$, with
   $$\circledast(F_R)=*\Bigl(e^a\otimes\bigl[(x_1\mathrm{Ric}+y_1\mathrm{Scal}\,\eta)_{ab}\gamma^b+(x_{13}\mathrm{Ric}+y_{13}\mathrm{Scal}\,\eta)_{ab}\gamma^b\Gamma\bigr]\Bigr).$$
   *Proof:* The trivial and $S^2_0$ constituents occur in $\Lambda^1\otimes\Lambda^l$ only for $l=1,13$, once each. So every admissible operator is "Ricci like + scalar like" on curvature. The Einstein condition H6 amounts to two linear conditions, $y_i=-\tfrac12x_i$, plus $(x_1,x_{13})\neq0$.
3. **Self-adjointness forces the Einstein ratio.** If $\circledast\in\mathcal S_{\rm sa}$, then $y_1=-\tfrac12x_1$ and $y_{13}=-\tfrac12x_{13}$, giving the Einstein tensor $G=\mathrm{Ric}-\tfrac12\mathrm{Scal}\,g$ of $Y^{14}$. *Proof:* $\mathrm{Hom}_K(\Lambda^3,\Lambda^1\otimes\Lambda^2)=\mathrm{Hom}_K(\Lambda^3,\Lambda^2\otimes\Lambda^{13})=1$, and the generators $W_A$, $W'_A$ (item 4) give $(x,y)=(-1,\tfrac12)$ exactly, in rational arithmetic (E2). Equivalently, $\Psi$ contracts $R$ with $\delta^{apq}_{mcd}$, and $\delta^{apq}_{mcd}R^{cd}{}_{pq}=-4G^a{}_m$ is the Cartan–Palatini Einstein $(n-1)$-form.
4. **Existence.** $W_A$ and $W'_A$ lie in $\mathcal S_{\rm sa}$: the form factor of $\mathrm{Tr}(\alpha\wedge[\Phi_{11}\wedge\xi]_A)$ is $\epsilon(a,I,b,c)$ and its algebra factor $i\mathrm{Tr}(X\{\Phi_I,Y\})$ is symmetric. They satisfy $W_A(F_R)=-*(G_{ab}e^a\otimes\gamma^b)$ and $W'_A(F_R)=-*(G_{ab}e^a\otimes\gamma^b\Gamma)$.
   - This is exact (E1, E2): the self-adjointness defect is $0$ on random integer $\mathfrak u(64,64)$ data with $|\beta|\sim10^4$, and $W_C$, $W'_C$ are exactly skew.
   - On curvature $W_A$ coincides with R-S$(C;A,A;c=1)$ (E2), but it is a different operator: the difference on a generic grade-2 2-form is nonzero. At the round-1 counterexample, $W_A$ gives $\beta_0=0$ both ways, where R-S gives $256$ and $0$.
5. **Within the draft's own family (8.1).** The single-term members mapping $\Omega^2\to\Omega^{13}$ are $r=11$ with $m=n=0$, and $r=1$ with $m=1$. Allowing $\Phi_r$ or $\Phi_r\Gamma$ and either bracket gives 8 operators. Their self-adjoint linear combinations are exactly $\mathrm{span}(W_A,W'_A)$: the defect matrix has exact rank 6 (E1). Every nonzero member of this span reproduces the Einstein tensor.

So **GU-P09 as stated is false**: Einstein-reproducing, Weyl-annihilating $\circledast$ with $\circledast\circ d_B$ formally self-adjoint exist, and in the draft's own form. ∎

**Theorem D4 (no-go for (9.7)).** Assume H4–H5, and let $\circledast_\varepsilon$ be *any* bundle map $\Lambda^2\otimes\mathrm{ad}\to\Lambda^{13}\otimes\mathrm{ad}$. No equivariance or parallelism is assumed. Suppose $\Upsilon(\varpi)=\circledast_\varepsilon F_{A_\omega}+\kappa_1*T_\omega$ is the $p$-gradient of some functional of $\varpi$. Then, at each point,
$$b(\alpha\wedge\circledast\xi)=\textstyle\sum_k\Omega_k\wedge I_k(\alpha\wedge\xi),$$
where $\{I_k\}$ is a basis of the $\mathrm{ad}$-invariant symmetric forms on $\mathfrak g$ ($\mathrm{Tr}(XY)$ and $\mathrm{Tr}X\,\mathrm{Tr}Y$ for $\mathfrak u(64,64)$; $\mathrm{Tr}$ for the simple algebras of H2) and the $\Omega_k$ are 11-forms. Consequently:
- **(a)** $\circledast\xi=\Omega_0\wedge\xi+\Omega_1\wedge\mathrm{Tr}(\xi)\mathbf 1$ preserves each Clifford grade. So $\circledast(F_R)$ is $\mathfrak{spin}(7,7)$-valued: $x_1=x_{13}=0$, and **the Einstein property and variationality of (9.7) are incompatible for every $\circledast$**. This holds for every $\mathrm{Ad}$- or $K$-invariant $b$.
- **(b)** If $\circledast$ is $K$-equivariant for a group with $(\Lambda^3V)^K=0$, then $\circledast=0$. The groups $\mathrm{Spin}(7,7)$, its identity component, and the chimeric $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$ all qualify; for the last, every summand of $\Lambda^3(V_4\oplus V_{10})$ is a nontrivial irreducible.

*Proof.*
1. A gradient has a symmetric second variation, so $\int b(\alpha\wedge D\Upsilon_T\beta)$ is symmetric for compactly supported $\alpha,\beta$ and every $T$.
2. The linearisation is $D\Upsilon_T\beta=\circledast d_B\beta+\circledast[T\wedge\beta]+\kappa_1*\beta$, where $[T\wedge\beta]=T\wedge\beta+\beta\wedge T$ has coefficients $[T_i,\beta_j]$.
3. At $T=0$ the principal symbol must be symmetric, so by D2, $\circledast=\Psi\wedge$. This step uses no parallelism, because the principal symbol of an adjoint is $-$(transpose).
4. The difference between $T$ and $0$ is zeroth order, so $P(\alpha,\beta;T)=b(\alpha\wedge\circledast[T\wedge\beta])$ is pointwise symmetric.
5. Take $\alpha=e^1A$, $T=e^2X$, $\beta=e^3C$. Using the total antisymmetry of $\Psi$, symmetry gives $\psi^{123}(A,[X,C])+\psi^{123}([X,A],C)=0$. So each $\psi^{abc}$ is an $\mathrm{ad}$-invariant symmetric form, which proves the displayed formula.
6. (a) For $b=\mathrm{Tr}$ (or $M$-equivalent), read off $\circledast$; $F_R$ is grade 2. A $K$-invariant $b$ cannot map grade 2 to grades 1 or 13, since $\Lambda^2\not\cong\Lambda^1$.
7. (b) Conjugation fixes the $I_k$, so $\Psi\in(\Lambda^3V)^K\otimes\mathrm{span}\{I_k\}=0$. ∎

*Sharpness.* A closed invariant 11-form would give $\int\Omega\wedge\mathrm{CS}$, which is variational. In $d=3$, $\Omega=1$ is Chern–Simons. Exact controls (E3): for $d=3$ Chern–Simons, $H=0$. For $d=14$ with the non-invariant $\Omega=e^{3\cdots13}\otimes\mathbf1$, both the self-adjointness defect and $H$ are $0$. For $W_A$ and $W'_A$, $H\neq0$ on $\mathfrak{so}(8,8)$ data. The $2\times2$ determinant is $-1.49\times10^{12}$, so no combination $xW_A+yW'_A\neq0$ passes.

**Corollary D5 (the true first variation of (9.4)).**
- For $\circledast\notin\mathcal S_{\rm sa}$, no $k_2$ gives $\circledast d_BT$: the term is $k_2(L+L^\dagger)T$, as in C6, now for all of $\mathcal S$.
- For $\circledast\in\mathcal S_{\rm sa}$, $k_2=\tfrac12$ is correct. Taking $\alpha=T$ forces $k_3=\tfrac13$ whenever $C(T)=p(T,\circledast(T\wedge T))\not\equiv0$, and the $\varpi$-EL equation of (9.4) is
$$\circledast F_{A}+\kappa_1*T-\tfrac13h(T)=0,\qquad p(\alpha,h(T))=H(\alpha,T;T):=p(\alpha,\circledast[T\wedge T])-p(T,\circledast[T\wedge\alpha]).$$
  By D4(b), $h\not\equiv0$ for every nonzero $\circledast\in\mathcal S_{\rm sa}$. For $\circledast=W_A$, tilted ($\Phi\to\varepsilon^{-1}\Phi\varepsilon$),
$$h(T)=\bigl[T\wedge[T\wedge\Phi_{11}]_C\bigr]_A .$$
  E4 checks this exactly: the derivative identity $C(T+\alpha)-C(T-\alpha)=2\,dC_T(\alpha)+2C(\alpha)$; the residual identity; ratio $H/p(\alpha,[T\wedge[T\wedge\Phi_{11}]_C]_A)=1$ on two samples.
- Under the other reading $[T,T]=[T\wedge T]_C$ (G12), the $T\wedge T$ coefficient is also off, by a factor of 2.

**Theorem D (answer to GU-P09; the two layers).** Under H1–H6, for each $\mathfrak g$ in H2 containing a copy of $V$ and each pairing in H4:
1. The *quadratic* part of (9.4) can be made variational together with the Einstein and Weyl properties. The possible operators form a 48-dimensional space ($\mathfrak u(64,64)$), and on curvature they reduce to $xW_A+yW'_A$ with the exact Einstein ratio.
2. The *equation* (9.7) is variational for **no** nonzero admissible $\circledast$, for any coefficients in (9.4) and any $K$-invariant pairing. Without equivariance it is variational only for generalized Chern–Simons operators, which have no Einstein output.

The answer does not depend on the interpretive choices listed in D.1. Outside the scope of D4: field-dependent multipliers, $\varpi$-dependent $\circledast$, and the $\varepsilon$-variation.

## D.3 Changes relative to round 1 (dated in `reports/CORRECTIONS.md`)

- GU-P09 (conjecture) is refuted (D3.4) and replaced by the stronger D4.
- GU-G13 is now proved for all of $\mathcal S$, not only R-S. Remedy (a) of the gap register is closed in the negative for (9.7).
- Prop. C6′ said the cubic term "requires $\tau(T,T,\alpha)=\kappa\tau(\alpha,T,T)$". That condition is not necessary, since the variation also involves $\tau(T,\alpha,T)$. The exact criterion is the Helmholtz condition of D4. C6′'s conclusion (the cubic term fails) stands and is now proved.
- Round-1 conclusions that are unchanged: C6 (for R-S), C5 and C4.

## D.4 Open (layer d)

- The $\varepsilon$-variation, and the Noether identity $\Xi=\mathfrak D\Upsilon$ (G15).
- Whether some $\varpi$-dependent multiplier or $\circledast_\omega$ rescues (9.7).
- Which Bianchi-identity property (fn. 10) singles out an operator. Weyl annihilation cannot be it: it is automatic (D3.1).
- The 14→4 reduction of $G(g_\aleph)$ (G14).
