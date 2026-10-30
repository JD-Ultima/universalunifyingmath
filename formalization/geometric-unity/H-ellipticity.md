# H. Euclidean ellipticity of the deformation complex for every self-adjoint Shiab operator

Claims: GU-057, GU-072. Gap: GU-G16. Code: `computation/geometric-unity/H_ellipticity.py` → `.out` (exact: integers, rationals, Gaussian rationals; exact sparse ranks by integer elimination). Records: `metrics/geometric-unity/H-euclidean-ellipticity.yaml`, `H-twisted-shiab.yaml`. Builds on D (D1–D3, H1–H6 of D.1) and G (Prop. G10, Thm G11). Layers: (a) source, (b) reconstruction, (c) proven, (d) open.

## H.0 Source and question (layer a)

- g11 (§12 (x), p. 65): "$\Upsilon=0$ carries an elliptic deformation complex in Euclidean signature once the redundant Euler–Lagrange equations are discarded" (fn. 14: Seiberg–Witten).
- Round 3 (G11): the symbol sequence is exact at no covector for $\mathrm{span}(W_A,W'_A)$ (the draft's self-adjoint family (8.1), D3.5) and for R-S, in (7,7) and (14,0). In (7,7) no $\circledast\in\mathcal S_{\rm sa}$ is elliptic (null cone). **Open:** the other 46 dimensions of $\mathcal S_{\rm sa}$ for $\mathfrak u$, and 10 of 12 for $\mathfrak{so}$, in (14,0).

## H.1 Hypotheses (layer b)

- **E1 (Euclidean real form).** $Y$ has signature (14,0), $K=\mathrm{Spin}(14)$, $\gamma_a$ Hermitian with $\gamma_a^2=1$. The spinor module $S=\mathbb C^{128}$ carries a positive-definite $K$-invariant Hermitian form, so the Wick rotation of the draft's construction (gauge algebra = Lie algebra of the unitary group of the spinor Hermitian structure, §8) has $\mathfrak g=\mathfrak u(128)$. We take this as the reading of "Euclidean signature".
  - *Alternative real forms.* $\mathfrak u(64,64)=\mathfrak u(S,h(i\Gamma\cdot,\cdot))$ is also $K$-stable, since $i\Gamma$ is Hermitian and central in $\mathrm{Cl}_{\rm ev}$. For $\mathfrak{so}(64,64)$ the compact analogue is $\mathfrak{so}(128)$ (Majorana, $\mathrm{Cl}(0,14)\cong M_{128}(\mathbb R)$), with the same Clifford grades $\{1,2,5,6,9,10,13,14\}$.
  - *Why the choice does not matter.* All of these are real forms of $\mathfrak{gl}(128,\mathbb C)$ or $\mathfrak{so}(128,\mathbb C)$, and as real $K$-modules each is $\bigoplus_k\Lambda^kV$ over its grades. The parameter space $\mathcal S_{\rm sa}(\mathfrak g_{\mathbb R})$ is a totally real form of $\mathcal S_{\rm sa}(\mathfrak g_{\mathbb C})$ of full dimension, hence Zariski dense. Every condition below is the non-vanishing of polynomials, so it holds on an open dense subset of the real form iff it holds somewhere on the complexification.
- **E2 (complex).** As in G10: Atiyah–Bott symbols of uniform order 1,
  $$0\to\mathfrak g\xrightarrow{\ \xi\otimes\cdot\ }\Lambda^1\mathfrak g\oplus\mathfrak g\xrightarrow{\ \sigma\delta_2(\xi)\ }\Lambda^{13}\mathfrak g\to0,\qquad\sigma\delta_2(\xi)(\alpha,\gamma')=\circledast_\varepsilon(\xi\wedge\alpha)-M_T(\xi\otimes\gamma'),$$
  with $M_T=\kappa{*}-\tfrac13Dh_T$ for the true Υ (D5) and $M_T=\kappa{*}$ for the printed one. *Elliptic* means exact for every $\xi\neq0$.
- **E3 (operators).** $\circledast\in\mathcal S_{\rm sa}$, i.e. $b(\alpha\wedge\circledast\zeta)=\Psi(\alpha\wedge\zeta)$ with $\Psi\in\mathrm{Hom}_K(\Lambda^3V,S^2\mathfrak g)$ (D2). Here $b$ is any nondegenerate Ad-invariant pairing (D.1 H4), and the background $(\varepsilon,T)$ and $\kappa$ are arbitrary.

Notation: Clifford monomials $E_K$ ($K\subset\{0,\dots,13\}$); $X_j$ is the grade-$j$ part; $q(e^A)=\gamma^A$ is the equivariant quantization; $\Gamma=E_{0\cdots13}$; $\mathrm{Tr}$ is normalised by $\mathrm{Tr}\,1=1$. $\mathrm{Ad}_\Gamma X=\Gamma X\Gamma^{-1}=(-1)^{\rm grade}X$ is the grade involution. $\xi=e^0$ and $W=\xi^\perp=\langle e^1,\dots,e^{13}\rangle$.

## H.2 Results (layer c)

**Lemma H1 (reduction).** Under E1–E3, in (14,0), the symbol sequence at $\xi\ne0$ is exact iff
1. $\kappa\neq0$, and
2. the antisymmetric form $\beta_\xi(\alpha,\alpha')=\Psi(\xi\wedge\alpha\wedge\alpha')$ on $W\otimes\mathfrak g$ is nondegenerate.

The criterion does not depend on the background, on $b$, or on whether the true or the printed Υ is used. It is the same for every $\xi\ne0$.

*Proof.*
- By D2, $b(\xi\gamma\wedge\circledast_\varepsilon(\xi\wedge\alpha))=\Psi_\varepsilon(\xi,\xi,\ldots)=0$. So $\mathrm{im}\,\circledast_\varepsilon(\xi\wedge\cdot)\subset A:=\mathrm{Ann}(\xi\otimes\mathfrak g)$, which has codimension $\dim\mathfrak g$.
- The γ′-columns pair with $\xi\otimes\mathfrak g$ through $\kappa\,\eta(\xi,\xi)\,b$; the $Dh_T$ part contains ξ twice (G11 proof). Hence $\mathrm{rank}\,\sigma\delta_2=\mathrm{rank}\,\circledast_\varepsilon(\xi\wedge\cdot)+\dim\mathfrak g$ when $\kappa\ne0$, and exactness needs $\mathrm{rank}\,\circledast_\varepsilon(\xi\wedge\cdot)=13\dim\mathfrak g$.
- That rank equals the rank of $(\alpha,\alpha')\mapsto\Psi_\varepsilon(\alpha\wedge\xi\wedge\alpha')$ on $\Lambda^1\mathfrak g/\xi\otimes\mathfrak g\cong W\otimes\mathfrak g$.
- $\Psi_\varepsilon=\Psi\circ(\mathrm{Ad}_\varepsilon\times\mathrm{Ad}_\varepsilon)$ has the same rank.
- $SO(14)$ is transitive on spheres and the symbol is homogeneous in ξ. ∎

*Check.* G11's (14,0) ranks satisfy $139234=122850+16384$.

**Lemma H2 (parity).** $\Gamma\in\mathrm{Spin}(14)$ covers $-1\in SO(14)$, so $\Psi(\omega;\hat X,\hat Y)=-\Psi(\omega;X,Y)$. Hence $\beta_\xi$ pairs $W\otimes\mathfrak g_{\rm ev}$ only with $W\otimes\mathfrak g_{\rm odd}$, and every block of $\mathcal S_{\rm sa}$ has $j+l$ odd. ∎

**Proposition H3 (exact basis of $\mathcal S_{\rm sa}$).** For $j<l$ define
$$O_{(j,l)}(\omega;X,Y)=\mathrm{Tr}\bigl(q(\omega)(X_jY_l+Y_jX_l)\bigr),\qquad E_{(j,l)}=\mathrm{Tr}\bigl(\Gamma q(\omega)(X_jY_l+Y_jX_l)\bigr).$$
The nonzero ones are:
- 24 $O$-types, with $|j-l|\in\{1,3\}$ and $3\le j+l\le25$;
- 24 $E$-types, with $E_{(j,l)}\leftrightarrow O_{(j,14-l)}$.

Together they form a basis of $\mathcal S_{\rm sa}$ (complex dimension 48). The coordinates are $p_0,\dots,p_{47}$, listed in H0 of the output.
- The block multiplicities agree block by block with the $D_7$ characters of D1 (no $S^2\Lambda^j$ blocks).
- *Real form $\mathfrak u(128)$:* $\Psi=\sum p_t\Psi_t$ is real iff $p_t\in i^{\epsilon_j+\epsilon_l}\mathbb R$, where $\epsilon_k=1$ for $k\equiv0,1\pmod 4$ and $0$ otherwise.
- *The draft's operators, exact coordinates:*

| Operator | Coordinates |
|---|---|
| $W_A$ | $-2$ on the 12 $O$-blocks $(0,3),(1,2),(2,5),(3,4),(4,7),(5,6),(6,9),(7,8),(8,11),(9,10),(10,13),(11,12)$, i.e. the blocks where $\gamma_J\gamma_L=+\gamma_L\gamma_J$ |
| $W'_A$ | $+2$ on 12 $E$-blocks |
| $W_C\circ\mathrm{Ad}_\Gamma$ | $\pm2i$ on the complementary 12 $O$-blocks |
| $W'_C\circ\mathrm{Ad}_\Gamma$ | $\pm2i$ on the complementary 12 $E$-blocks |
| $W_C$, $W'_C$ | not in $\mathcal S_{\rm sa}$ (the exact system is inconsistent) |

*Proof.*
- Invariance follows from trace cyclicity, equivariance of $q$, and $K$-invariance of the grade projections. This is checked exactly on 4000 samples per signature.
- Symmetry in $(X,Y)$ holds by construction.
- Distinct types have disjoint monomial supports ($A\triangle J\triangle L=\emptyset$ or everything, with fixed grades), so they are independent. Each type is nonzero on an explicit witness.
- 48 equals $\dim\mathcal S_{\rm sa}$ (D1).
- The coordinates of $W_A$ and $W'_A$ are the unique exact solutions of 4000 equations generated by the G_symbols operator code. They reproduce G_symbols' symbol matrix at $e^0$ entrywise (1 277 952 nonzeros), in both signatures. ∎

**Theorem H4 (Spin(13) block reduction).** As a $\mathrm{Spin}(13)$-module,
$$W\otimes\mathfrak g\cong4\Lambda^0\oplus8(\Lambda^1\oplus\cdots\oplus\Lambda^6)\oplus4(H_1\oplus\cdots\oplus H_6),$$
where $H_j$ is the hook $[2,1^{j-1}]$ (dims 90, 715, 2925, 7722, 14300, 19305).
- *Copies.* The copies are labelled by Clifford data:
  - the maps $\omega\mapsto\sum e^a\otimes e^a\wedge\omega$ and $\omega\mapsto\sum e^a\otimes\iota_a\omega$;
  - the Hodge star on $W$;
  - left multiplication by $\gamma^0$.

  This uses $\mathfrak g|_{\mathrm{Spin}(13)}=\bigoplus_m\Lambda^mW\otimes\{1,\gamma^0\}$, $W\otimes\Lambda^m=\Lambda^{m+1}\oplus\Lambda^{m-1}\oplus H_m$, $\Lambda^m\cong\Lambda^{13-m}$ and $H_m\cong H_{13-m}$.
- *Block form.* By Schur, $\beta_\xi=\bigoplus_\rho\begin{psmallmatrix}0&C_\rho\\-C_\rho^{T}&0\end{psmallmatrix}\otimes\langle,\rangle_\rho$, where $C_\rho$ is the odd-copies × even-copies matrix. It is $2\times2$ for $\Lambda^0$ and $H_j$, and $4\times4$ for $\Lambda^1,\dots,\Lambda^6$. Its entries are listed in H2 of the output.
- *Rank.* $\mathrm{rank}\,\beta_\xi=\sum_\rho2\,\mathrm{rank}(C_\rho)\dim\rho$.

*Checks (exact).*
- Antisymmetry and parity of every block.
- The same $C_\rho$ from three independent test vectors.
- Orthogonality between different test vectors and between different ρ.
- $\sum m_\rho\dim\rho=212992$.
- The predicted rank equals the exact rank of the full $212992\times212992$ matrix of $\beta_\xi$ at six points: $W_A$, $W'_A$ and $W_A+W'_A$ give 122850, equal to G_symbols; a random point gives 212992; one point on each of two components gives the predicted deficits 1430 and 572. In (7,7) the Schur ranks of the non-null pencil also reproduce G_symbols (106496 for $W_A\pm W'_A$). ∎

**Theorem H5 (main; Euclidean, $\mathfrak g=\mathfrak u(128)$ or any real form of $\mathfrak{gl}(128,\mathbb C)$).** Under E1–E3 the deformation complex of (9.4) is elliptic in (14,0) **iff $\kappa\ne0$ and $\Delta(\Psi)\ne0$**, where
$$\Delta=Q_0\prod_{j=1}^{6}Q_j\prod_{i=2}^{6}D_i .$$
The factors are:
- $Q_j=p_{O(j,j+1)}p_{O(13-j,14-j)}+(-1)^jp_{E(j,13-j)}p_{E(j+1,14-j)}$ for $j=1,\dots,6$, from $H_j$. $\Lambda^0$ also gives $Q_1$.
- $Q_0=p_{O(0,3)}p_{O(11,14)}+p_{E(0,11)}p_{E(3,14)}$. The $\Lambda^1$ determinant is $Q_0Q_1$.
- $D_i=\det C_{\Lambda^i}$, $i=2,\dots,6$: quartics in 16, 16, 16, 16 and 12 variables, printed in H3 of the output.

Consequently:
- **(a)** The elliptic set $V=\{\Delta\ne0\}$ is nonempty and Zariski-open, so it is open and dense in the 48-dimensional $\mathcal S_{\rm sa}$. Its complement is a hypersurface of degree 34, the union of 12 hypersurfaces.
  - The 12 factors are irreducible over $\mathbb Q(i)$. The quadrics are absolutely irreducible, being of rank 4.
  - On the real form $\mathfrak u(128)$, each factor is a real polynomial up to a unit and takes both signs. So $V\cap\mathcal S_{\rm sa}(\mathfrak u(128))$ is the complement of 12 real hypersurfaces.
  - All 48 coordinates enter Δ.
- **(b)** Not elliptic: every member of $\mathrm{span}(W_A,W'_A)$, since all six $\det C_{H_j}$ and $\det C_{\Lambda^0},\det C_{\Lambda^1}$ vanish identically on it. The $\Lambda^{2..6}$ determinants are $c_i(x^2+y^2)^2$ in (14,0) and $c_i(x-y)^2(x+y)^2$ in (7,7). This explains G11's rescaling lemma and its degenerate lines.

*Proof.* Lemma H1, Lemma H2 and Theorem H4. Exactness at one ξ is exactness at all. The determinants are exact symbolic computations. Nonemptiness is witnessed by an exact rank of 212992 at an integer point (H3), which is a $\mathfrak u(128)$-real point by H3's real structure. ∎

**Proposition H6 (Einstein ratio; the draft's form (8.1)).**
1. On Levi-Civita curvature only $p_{O(1,2)}$ (grade-1 output) and $p_{E(2,13)}$ (grade-13 output) contribute, each with the exact Einstein tensor. Outputs of grade 3, 5, 9 and 11 vanish. This is checked exactly for every type and agrees with D3.2.
2. Every elliptic Ψ has a nonzero Einstein output: $Q_1\ne0$ forces $(p_{O(1,2)},p_{E(2,13)})\neq0$. The output has the exact ratio $y=-x/2$, forced by self-adjointness (D3.3). **So ellipticity implies the Einstein property.**
3. **No elliptic operator is of the draft's form (8.1).** Its self-adjoint members are exactly $\mathrm{span}(W_A,W'_A)$ (D3.5), and that span lies in $\{\Delta=0\}$. ∎

**Proposition H7 (a minimal elliptic extension of (8.1)).** Let $\Phi^E=\sum e^{abc}\otimes\gamma_{abc}$ be the $\mathfrak u(128)$-valued Φ³. Consider the real four-parameter family
$$\circledast_{x,y,z,w}=x\,W_A+y\,W'_A+z\,W_C\circ\mathrm{Ad}_\Gamma+w\,W'_C\circ\mathrm{Ad}_\Gamma .$$
Here all operators are built from $*\Phi^E$ (and its Γ-partner), and the coefficients are real.
- Every member lies in $\mathcal S_{\rm sa}$, is $\mathfrak u(128)$-real, annihilates Weyl, has the Bianchi property (G4), and on curvature equals $xW_A+yW'_A$.
- On the family:
  - $\det C_{H_j}$, $\det C_{\Lambda^0}$ and $\det C_{\Lambda^1}$ are powers of $(xz+yw)$ up to a unit;
  - $\det C_{\Lambda^i}=a_i(x^2+y^2-z^2-w^2)^2+b_i(xz+yw)^2$ with $a_ib_i>0$ for $i=2,\dots,6$ (exact identities; $(a_i,b_i)\propto(55,256),(15,64),(6,25),(70,289),(35,144)$).
- Hence **$\circledast_{x,y,z,w}$ is elliptic iff $\kappa\ne0$ and $xz+yw\ne0$.**
- *Example.* $\circledast_\sharp=W_A+W_C\circ\mathrm{Ad}_\Gamma$, i.e. $\circledast_\sharp(\eta)=[\Phi_{11}\wedge\eta]_A+[\Phi_{11}\wedge\Gamma\eta\Gamma^{-1}]_C$.
  - Exact rank of $\beta_\xi$: 212992/212992.
  - Exact rank of $\sigma\delta_2(e^0)$ at $\kappa=1$: $229376/229376$, so the sequence is exact.
  - $\circledast_\sharp(F_R)=W_A(F_R)$, the Einstein tensor of $Y$.
  - It is not of form (8.1), because of the grade involution $\mathrm{Ad}_\Gamma$, the action of the lift of $-1\in O(14)$.

  *Proof.* Proposition H3 (coordinates), Theorem H4 and exact sympy identities. ∎

**Theorem H8 ($\mathfrak{so}(64,64)$ / $\mathfrak{so}(128)$: never elliptic).** For every Ψ in the 12-dimensional $\mathcal S_{\rm sa}(\mathfrak{so})$, every κ and every background, $\mathrm{rank}\,\beta_\xi\le60580<105664=13\cdot8128$. The bound is attained generically and by $W_A$.

The corank $45084=45058+26$ comes from two sources:
- 45058: the odd and even multiplicities differ for $\Lambda^0$ and all six $H_j$. The odd part of $\mathfrak{so}$ ($\Lambda^2S^+\oplus\Lambda^2S^-$, dim 4032) is not isomorphic, as a Spin(13)-module, to the even part ($\mathfrak{gl}(64)$, dim 4096).
- 26: the grade-14 copy $W\otimes\Gamma$ has no partner, because its only couplings, $O_{(11,14)}$ and $E_{(3,14)}$, need grades 3 and 11, which are not in $\mathfrak{so}$.

So G11's $\mathfrak{so}$ figure for $W_A$ (60580 of 105664) is the maximum over all of $\mathcal S_{\rm sa}(\mathfrak{so})$. ∎

**Remark H9 ((7,7)).** The same blocks over $\mathbb C$ give exactness at non-null ξ for generic Ψ: a random point has exact rank 212992 at $e^0$ in (7,7). In split signature the only obstruction is therefore the null cone (G11.1). The degenerate lines $W_A\pm W'_A$ of G11(3) are the factors $(x\mp y)^2$ of the pencil determinants.

## H.3 What is verified, precisely

| Statement | How | Scope |
|---|---|---|
| H1, H2 | proof; numeric consistency with G11 ($139234=122850+16384$) | all Ψ, κ, $(\varepsilon,T)$, $b$, true/printed Υ |
| H3 | proof plus exact checks (character multiplicities, invariance and symmetry on 4000 samples, witnesses, overdetermined exact fits, entrywise identity with G_symbols) | all 48 types; (14,0) and (7,7) |
| H4 | Schur's lemma plus exact block matrices; full exact sparse ranks at 6 points (14,0) and 1 point (7,7) agree with the Schur prediction | all Ψ |
| H5 | exact symbolic determinants; factorisation over $\mathbb Q$ and $\mathbb Q(i)$; sign changes at real integer points | all Ψ ∈ $\mathcal S_{\rm sa}$ |
| H6–H7 | exact rational curvature contraction; exact sympy identities; exact ranks of $\beta$ and $\sigma\delta_2$ for $\circledast_\sharp$ | all Ψ; the 4-parameter family |
| H8 | exact Schur blocks restricted to the $\mathfrak{so}$ grades; full exact ranks at a random point and at $W_A$ | all 12 dimensions |

## H.4 Grading and changes to earlier conclusions

| Item | Status |
|---|---|
| H1–H8 | **Proven** (exact) |
| Euclidean ellipticity of the bosonic complex | **True for an open dense set of self-adjoint Shiab operators** ($\kappa\ne0$, $\Delta\ne0$). **False** for the draft's family (8.1), for R-S, and for every $\mathfrak{so}(64,64)$ operator. In (7,7) it fails always (null cone) |
| Einstein + self-adjoint + elliptic | coexist (H6.2, H7). Not in the draft's form (8.1) |

Claims:
- **GU-072** (bosonic): F → **C**. It holds for generic $\circledast\in\mathcal S_{\rm sa}(\mathfrak u)$ and fails for the operators the draft writes down.
- **GU-057**, ellipticity part: F → **C** on the same terms. The part $\delta_2\delta_1=\mathrm{ad}\Upsilon$ stays V.
- **GU-G16**: the negative verdict of round 3 is restricted to (8.1), R-S, $\mathfrak{so}$ and (7,7).

These are logged in `reports/CORRECTIONS.md`. G11 itself is unchanged: every number in it is reproduced here.

## H.5 Open (layer d)

- Whether fn. 10's lost operator, or a principle stated in the draft, selects a point of $V$. H6.2 shows ellipticity forces the Einstein output, but Weyl annihilation and the Bianchi property select nothing (D3.1, G4).
- Ellipticity in any sense for non-self-adjoint ⊛ beyond R-S. Such operators are outside the variational class anyway.
- Analytic consequences (Fredholm index on compact $Y^{14}$, moduli dimension), which need global data the draft does not give.
- The equation (9.7) is still not variational for any ⊛ (D4). Ellipticity here refers to the complex of the true EL equation (D5/G10).

**Coordinator reproduction (2026-10-30).** `H_ellipticity.py` was re-run and its output is identical to the committed `.out` apart from timing annotations.
