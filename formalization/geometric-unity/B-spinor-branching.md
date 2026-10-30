# B. Spinors of $\mathrm{Spin}(7,7)$, their branching to $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$, and the generation count

Claims: GU-016…020, GU-025…027, GU-029, GU-058…066. Code: `clifford77.py` → `clifford77.out` (explicit real model) and `spinor_branching.py` → `spinor_branching.out` (weights). Metric record: `metrics/geometric-unity/B-spinor-branching.yaml`. Imported lemmas: L-B1…L-B4.

## B.0 Source assertions (layer a)

- (3.18) $\mathrm{Cl}_{7,7}\cong\mathbb R(128)$. (3.19), (3.25), (3.30): $\mathrm{Spin}(7,7)\to GL(64,\mathbb R)^{L}\times GL(64,\mathbb R)^{R}$, $\to SO(64,64)\to U(64,64)$, with the left and right real Weyl spinors dual to each other. (3.20) "$\rho_{\mathrm{Spin}}:\mathrm{Spin}(14,\mathbb C)\to U(128,\mathbb C)$". (3.21), (3.22), (3.26), (3.27), (8.6): grade lists for $\mathfrak{so}(7,7)$, $\mathfrak{gl}(64,\mathbb R)$, $\mathfrak{so}(64,64)$, $\mathfrak u(64,64)$ inside $\mathrm{Cl}$. Fn. 3: "should be checked".
- (4.1)–(4.3) Under observation a chimeric spinor appears as $\$_x(TX)\otimes\$_x(N)$. The internal dimension is $2^{\dim N/2}$ (Dirac) or $2^{\dim N/2-1}=16$ (Weyl).
- (4.6), (4.9), §11.4: $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)\subset\mathrm{Spin}(7,7)$. The maximal compact subgroup $\mathrm{Spin}(6)\times\mathrm{Spin}(4)\cong SU(4)\times SU(2)\times SU(2)$ is Pati–Salam.
- (11.3) $R(U\oplus V)=R(U)\otimes S(V)\oplus S(U)\otimes R(V)\oplus S(U)\otimes S(V)$. (11.5)–(11.6) is the decomposition of the "rolled up Fermionic complex" into $\mathcal F_{1/2}^\pm$, $\mathcal Q_{3/2}^\pm$, $\mathcal Z_{1/2}^\pm$, built from $2_\mp\otimes16_\pm$, $6_\mp\otimes16_\pm$ and $2_\mp\otimes144_\pm$ with $64^\pm$ and $832^\pm$. The tables on pp. 52–53 give the SM quantum numbers.
- (12.10), (12.22), (12.28): $\nu\in\Omega^0(\$_L)$ gives the first generation and the $\$$-part of $\zeta\in\Omega^1(\$_R)$ gives the second. The third generation is an "imposter" in the Rarita–Schwinger part. (12.9), (12.20): the fundamental theory is non-chiral, and chirality is "merely effective".

## B.1 Reconstruction (layer b)

**R-B1.** $\eta=(+^7,-^7)$. $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)\subset\mathrm{Spin}(7,7)$ is the block subgroup for $\mathbb R^{7,7}=\mathbb R^{1,3}\oplus\mathbb R^{6,4}$. This is Theorem A4's $(7,7)$ case with $V\cong\mathbb R^{6,4}$ and $H\cong\mathbb R^{1,3}$. **R-B2.** The field reality condition is left open in the source. The draft uses Dirac/$U(64,64)$ in §3.6 and real Weyl/$GL(64,\mathbb R)$ in (3.25). We compute three readings: Majorana–Weyl (MW, real), complex Weyl, and Dirac. **R-B3.** SM embedding (L-B4): $Y=\tfrac12(B-L)+T_{3R}$, and $n:=6Y$ as in the draft's tables.

## B.2 Results (layer c)

**Theorem B1 (Clifford structure; all verified exactly in `clifford77.out`).**
1. Explicit real matrices $\gamma_0,\dots,\gamma_{13}\in M_{128}(\{0,\pm1\})$ satisfy $\{\gamma_a,\gamma_b\}=2\eta_{ab}$. All $2^{14}$ monomials are trace-orthogonal, so $\mathrm{Cl}(7,7)\xrightarrow{\ \sim\ }\mathbb R(128)$. This is (3.18).
2. $\Gamma=\gamma_0\cdots\gamma_{13}$ is real with $\Gamma^2=+1$, so real Majorana–Weyl modules $S^\pm_{\mathbb R}$ of dimension 64 exist, as required by L-B1, since $p-q=0$.
3. There is a unique, up to scale, form $B$ with $B$-adjoint $\gamma_a^*=-\gamma_a$. It is symmetric, has signature $(64,64)$, and both $S^\pm$ are isotropic. Its adjoint acts on grade $k$ by $(-1)^{k(k+1)/2}$, which is Clifford conjugation. Hence $\mathfrak{so}(B)$ consists of the grades $k\equiv1,2\pmod 4$, $\dim\mathfrak{so}(B)=8128$, and this is exactly (3.26). $\mathfrak u(B_{\mathbb C})=\mathfrak{so}(B)\oplus i\cdot\{k\equiv0,3\}$ is exactly (3.27)/(8.6). The even skew grades $\{2,6,10,14\}$ have dimension $4096=64^2$, which is $\mathfrak{gl}(S^+)$ as in (3.22). The odd skew grades have dimension $4032=2\binom{64}{2}$, which is $\Lambda^2S^+\oplus\Lambda^2S^-$. The other form $B'$ ($\gamma_a^*=+\gamma_a$) is antisymmetric, with $\mathfrak{sp}$ of dimension 8256.
4. Consequently $\mathrm{Spin}(7,7)\to GL(S^+_{\mathbb R})\subset SO(B)\subset U(B_{\mathbb C})=U(64,64)$, which confirms (3.19), (3.25) and (3.30).
5. **(3.20) is ill-typed as written.** A holomorphic homomorphism from $\mathrm{Spin}(14,\mathbb C)$, a complex semisimple group, into a compact group is trivial. Read charitably it means $\mathrm{Spin}(14,\mathbb C)\to GL(128,\mathbb C)$. The same holds for (12.23).

**Theorem B2 (a Majorana–Weyl spinor of $\mathrm{Spin}(7,7)$ is exactly one generation).** As a real $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$-module,
$$S^+_{\mathbb R}(7,7)\;\cong\;\bigl(S^+(1,3)\otimes_{\mathbb C}S^+(6,4)\bigr)_{\mathbb R}=\bigl((2,1)\otimes 16\bigr)_{\mathbb R},$$
irreducible and of complex type. Restricted to $SL(2,\mathbb C)\times SU(4)\times SU(2)_L\times SU(2)_R$ it is $(2,1)\otimes[(4,2,1)\oplus(\bar4,1,2)]$. Restricted to the SM it is one generation of left-handed Weyl fermions including $\nu^c$, with hypercharges $n=6Y\in\{+1,+2,-4,-3,+6,0\}$. This equals the draft's p. 52 table exactly (B5: `True`).
*Proof.* Over $\mathbb C$, $S^+(14)|=(2_+\otimes16_+)\oplus(2_-\otimes16_-)$ (B1). By L-B1 the Weyl modules of $\mathrm{Spin}(1,3)$ and $\mathrm{Spin}(6,4)$ ($p-q\equiv\mp2$) are swapped by complex conjugation, so the real structure of $S^+_{\mathbb R}$ swaps the two summands. Explicitly, $\omega_4=\gamma_0\gamma_7\gamma_8\gamma_9$ has $\omega_4^2=-1$, commutes with all 51 generators of $\mathfrak{spin}(1,3)\oplus\mathfrak{spin}(6,4)$, and on $S^+$ satisfies $\omega_{10}=-\omega_4$: the 4d and 10d chiralities are locked. The real commutant is 2-dimensional (SVD gap: smallest retained singular value 1.58, tolerance $2.3\times10^{-8}$), so it is $\mathbb C=\mathrm{span}\{1,\omega_4\}$, and L-B3 applies. The SM content is B3/B4 of `spinor_branching.out`. ∎

**Theorem B3 (chirality-refined Rarita–Schwinger branching).** With $R^\epsilon:=V\otimes S^\epsilon\ominus S^{-\epsilon}$ (dimensions $6$, $144$, $832$ for $D_2,D_5,D_7$), and with $S^\epsilon(14)|=\bigoplus_{\sigma=\pm}2_\sigma\otimes16_{\sigma\epsilon}$:
$$R^\epsilon(14)|=\bigoplus_{\sigma=\pm}\Bigl[(6_\sigma\otimes16_{\sigma\epsilon})\oplus(2_\sigma\otimes144_{\sigma\epsilon})\oplus(2_\sigma\otimes16_{-\sigma\epsilon})\Bigr],\qquad \bigoplus_\sigma 2_\sigma\otimes16_{-\sigma\epsilon}=S^{-\epsilon}(14)|$$
(B1). This is (11.3) with chirality labels. The draft's $\mathcal F,\mathcal Q,\mathcal Z$ blocks (11.6) are these summands. The extra spin-½ summand of $R^\epsilon(14)$ restricts like $S^{-\epsilon}(14)$, the spinor of the *opposite* chirality.

**Theorem B4 (field content of $\nu\oplus\zeta$).** Take $\nu\in S^{+}$ and $\zeta\in T^*Y\otimes S^{-}=S^+\oplus R^-$, the chiralities of (12.28). The left-handed spin-½ content (B8, B9) is:

| reading | LH spin-½ $\mathrm{Spin}(10)$ content | net $\#16_+-\#16_-$ |
|---|---|---|
| Majorana–Weyl (real) | $3\times16_+\ \oplus\ 144_-$ | $3$ |
| complex Weyl | $6\times16_+\ \oplus\ 2\times144_-$ | $6$ |
| Dirac ($S^+\oplus S^-$) | $6\times16_+\oplus6\times16_-\oplus2\times144_+\oplus2\times144_-$ | $0$ |
| MW with $\zeta\in T^*\otimes S^+$ (same chirality) | $16_+\oplus2\times16_-\oplus144_+$ | $-1$ |

*Note (2026-10-30, round 3, `G-full-el-fermions.md` Theorem G7).* The table is representation theory and is unchanged. The Majorana–Weyl rows are not dynamically consistent with the fermion operator (9.16) once its ζζ (Shiab) block is kept: in (7,7) a single Majorana–Weyl ζ has no Rarita–Schwinger kinetic term, and a single Majorana–Weyl spinor no Dirac kinetic term. The minimal consistent real reading is two MW copies, equivalent to complex Weyl (second row). Statements that use the MW row (e.g. "one generation from ν", GU-027, 063, 064) are conditional on dropping the ζζ block or adding non-Spin(7,7)-covariant structure.

In the MW reading the spin-3/2 content is $(3,2)\otimes16_-$, which is the draft's $\mathcal Q^+_{3/2}$. Its SM content equals the draft's second p. 52 table (B5: `True`).
1. **The p. 53 table is exactly the SM content of $144_-$** (B5: `True`, with SU(3) conjugations transcribed from a 220-dpi render), including all electric charges, e.g. $Q\in\{\pm\tfrac43,\pm\tfrac53,\pm2\}$ in exotic multiplets.
2. Pati–Salam: $144_-=(20,1,2)\oplus(\overline{20},2,1)\oplus(4,2,1)\oplus(4,2,3)\oplus(\bar4,1,2)\oplus(\bar4,3,2)$ (B3). It contains the Pati–Salam content of one generation, $(4,2,1)\oplus(\bar4,1,2)$.
3. **Generation count.** In the MW reading, the spin-½ sector contains exactly three $\mathrm{Spin}(6,4)$-$16_+$ multiplets. Two sit in the 64-dimensional $\mathrm{Spin}(7,7)$ modules ($\nu$ and $\gamma$-trace$(\zeta)$) and one in the 832 ($R^-(14)$). This matches "2+1" in (12.22). At the level of $SU(3)\times SU(2)\times U(1)$ quantum numbers, however, the content contains **four** copies of a generation, because $144_-$ contains one more, as maximal multiset containment shows (B6, B8). Removing all four leaves 128 exotic spin-½ states.
4. The draft labels the "imposter" quarks and leptons *inside* the 144 table (p. 53). In (12.22) the imposter is the summand $\$(TX)\otimes\$(N)\subset R(TY)$. These are different subspaces: the first is in $2\otimes144$, the second in $2\otimes16$. See gap GU-G17.

*Proofs.* All are weight-multiset identities with exact integer arithmetic. Decompositions are by highest-weight peeling against explicitly constructed characters, and every subtraction is checked to be non-negative and to exhaust the module. ∎

**Proposition B5 ((4.2)–(4.3)).** For $\dim N=n(n+1)/2$ even, the Dirac and Weyl dimensions are $2^{n(n+1)/4}$ and $2^{n(n+1)/4-1}$. For $n=4$ this gives $32$ and $16$ (B7). The arithmetic is correct. That $N$ is even requires $n\equiv0,3\pmod4$.

## B.3 Assessment

- **Holds (layer c):** GU-016, 017, 019, 020, 025, 026, 058, 059, 060, 065, 066, and the representation content of 061 and 062. GU-027, 063 and 064 hold under the choices below. The rep-theoretic core of the fermion sector is internally consistent and reproduces the draft's tables exactly, provided one fixes the MW reading, the opposite Spin(7,7) chiralities for $\nu$ and $\zeta$, and $Y=\tfrac12(B-L)+T_{3R}$.
- **Fails or ill-typed:** (3.20) and (12.23) as written.
- **Depends on an unstated choice (gap GU-G07):** the reality condition. The draft's main bundle is complex Dirac ($U(64,64)$). In that reading the content is vector-like with net chirality 0, and $\nu$ alone gives two generations plus two mirror generations, not one (viii of (12.11)). The "one generation from $\nu$" statement holds only for MW.
- **Not established (layer d):** the dynamical claims that the "dark/looking-glass" sector decouples (12.12–12.20), that the 144 and the spin-3/2 matter are heavy, and that the three $16_+$ differ in mass. No mechanism, scale or coupling is given.
- **Empirically testable (layer e), in principle:** the predicted extra spin-½ states (SM content of $144_-$ minus one generation: 128 states, including colour sextets and octets and charges $\pm\tfrac43,\pm\tfrac53,\pm2$) and spin-3/2 states with conjugate generation quantum numbers. With no mass scale these are not yet falsifiable predictions.
- **Physical caveat (not mathematics):** the internal symmetry acting on the 16 is the non-compact $\mathrm{Spin}(6,4)$. Only its maximal compact subgroup is Pati–Salam, and the draft acknowledges the unitarity issue (§4.2).

## Independent cross-check (coordinator)

`computation/geometric-unity/xcheck_144.py` checks the representation content by a different method from `spinor_branching.py`. It does not use weight multisets. Instead it builds an explicit $\mathrm{Cl}(10)$ on $\mathbb C^{32}$ and realises $144_\pm=\ker(V\otimes S_\pm\to S_\mp)$ concretely. It then reads off SM labels by joint diagonalisation of the SU(3) and SU(2)$_L$ Casimirs, $Y$ and $B-L$, and separates $3$ from $\bar3$ by triality $3(B-L)\bmod 3$.

- $16_+$ is exactly one left-handed generation with $\nu^c$, which fixes the sign of $B-L$.
- **$144_-$ contains the SM content of one $16_+$ generation; $144_+$ contains that of its conjugate.**
- Both have 23 SM irreducibles and electric charges $\{0,\pm\tfrac13,\pm\tfrac23,\pm1,\pm\tfrac43,\pm\tfrac53,\pm2\}$.

This confirms the representation-theoretic part of B6/B8. The count "four generation-equivalents" additionally uses the chirality assignment $3\times16_+\oplus144_-$ of Theorem B3, which this cross-check does not re-derive.
