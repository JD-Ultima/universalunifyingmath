# Geometric Unity: Phase 1 report, with Phase 2 results on three claims

**Scope.** E. Weinstein, *Geometric Unity: Author's Working Draft v1.0* (1 Apr 2021). Hypothesis-neutral: critiques are treated as leads, and every verdict below comes from our own derivation or computation. Files: `sources/`, `claims/`, `gaps/`, `patterns/` (`geometric-unity.md` in each), `formalization/geometric-unity/{A,…,J}-*.md`, `computation/geometric-unity/` (29 scripts, outputs committed; see its README), and `metrics/geometric-unity/*.yaml`.

## 1. Sources

- **Obtained: the 2021 draft**, via a third-party GitHub mirror. sha256 `3f28d742…86d4`, 69 pp., PDF timestamps 2021-04-01. The original host was unreachable here, so the copy is unverified against the original.
- **Not obtained:** the 2013 Oxford lectures, the 2020 "First Look" video, all transcripts, and Nguyen–Polya (2021). We attribute no content to any of them.

## 2. Phase 1

- **Claim Register.** 72 claims (GU-001…072). Primary types: 40 theorem-claims, 17 definitions, 6 physical predictions, 6 analogies, and 1 each numerical observation, conjecture and source statement.
- **Status after Phase 2.** 34 verified; 19 conditional, i.e. holding under a stated reading or repair; 6 open or not well-posed; 6 false as stated (GU-018, 042, 044, 045, 046, 053); 2 trivially true; 5 not mathematical. A re-graded claim is counted by its current status; GU-057 (V/F/C) and GU-072 (F/C) are counted as C. (Corrected 2026-10-30: the counts had taken each claim's pre-Phase-2 status; see `reports/CORRECTIONS.md`.) Counts come from `computation/geometric-unity/register_stats.py`.
- **Specification gaps.** 22 gaps (GU-G01…G22). Three are proved inconsistencies:
  - G11: tilt sign;
  - G13: the field equation is not variational;
  - G17: "imposter" generation has two definitions.
  The underspecified core is the Zorro metric (G06), the Shiab operator (G10), the 14→4 reduction (G14), the reality conditions (G07) and symmetry breaking (G18).

The draft's finite-dimensional algebra is largely correct and checkable. Its analytic core, from the action to the field equations, is not specified precisely enough to be well-posed. Under the most natural reconstruction it is inconsistent.

## 3. Phase 2 results

| # | Result | Status |
|---|---|---|
| A | **Chimeric signature.** All $GL(4)$-natural fiber metrics form the DeWitt family $c\,G_\lambda$. The fiber signature is $(7,3)$ for $\lambda<\tfrac14$ and $(6,4)$ for $\lambda>\tfrac14$ with base $(1,3)$ or $(3,1)$, and $(10,0)$ or $(9,1)$ with base $(4,0)$. Among bases $(1,3)$, $(3,1)$, $(4,0)$, the total on $V\oplus H^*$ is $(7,7)$ **iff** the base is Lorentzian, $\lambda>\tfrac14$, and a fixed relative block sign holds. A base $(2,2)$ with $\lambda>\tfrac14$ also gives $(7,7)$, for either sign. A Riemannian base never gives $(7,7)$. For $n\ge2$, neutral signatures occur only for bases $(1,3)$, $(3,1)$, $(2,2)$ [$(7,7)$] and $(2,3)$, $(3,2)$ [$(10,10)$]. $w_2(C)=\pi^*w_2(X)$, but for an indefinite metric the $\mathrm{Spin}(7,7)$ lifting obstruction is $w_2(C_+)+w_2(C_-)=\pi^*w_2(X)+w_1(\ell)^2$, with $\ell$ the time line (Karoubi). Hence the draft's $\mathrm{Spin}(7,7)$ bundle **does not exist over $\mathrm{Met}(X)$** for any $X$. Over the time-oriented double cover of (3.33) it exists iff $w_2(X)=0$, i.e. iff $X$ is spin when $X$ is orientable. The draft's fiber in (3.33) is the double cover of the fiber of $\mathrm{Met}(X)$. | **Proven** (closed form, exact computation). Spin statement corrected 2026-10-30 (`A_spin_obstruction.out`) |
| B | **Spinors.** $\mathrm{Cl}(7,7)\cong\mathbb R(128)$. Majorana–Weyl spinors exist. All of the draft's Clifford grade lists for $\mathfrak{so}(64,64)$ and $\mathfrak u(64,64)$ are exact. One MW spinor restricted to $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$ is $((2,1)\otimes16)_\mathbb R$: exactly one SM generation with $\nu^c$, Pati–Salam $(4,2,1)\oplus(\bar4,1,2)$. With $\nu\in S^+$, $\zeta\in T^*\otimes S^-$, in the Majorana–Weyl reading, the left-handed spin-½ content is $3\times16\oplus144$, plus spin-3/2 $(3,2)\otimes\overline{16}$. Round 3 (G6–G7) shows that this reading has no ζ kinetic term. All three of the draft's SM tables (pp. 52–53), including the 144 table with its charges, are reproduced exactly. The SU(3) bars of the p. 53 table were transcribed from a 220-dpi render. At the SM level, however, the 144 contains a **fourth** generation-equivalent plus 128 exotics. In the complex Weyl reading the content doubles. In the Dirac/$U(64,64)$ reading it is vector-like. "Spin(14,ℂ)→U(128)" is ill-typed. | **Proven** (exact weights and an explicit real model) |
| C1 | **Inhomogeneous gauge group.** $\mathcal G=\mathcal H\ltimes\mathcal N$ is a group with the stated actions. $\tau$ is a homomorphism and equals the right stabiliser of $A_0$. **Lemmas 6.5 and 7.2 are false as printed**, with residual $-2h^{-1}d_{A_0}h$. They hold for the opposite-sign tilt $\tau'$, the *left* stabiliser of $A_0$, but then Lemma 6.3 fails. One global sign choice makes Secs. 6–7 consistent, but it changes at least two printed formulas: $\tau'$ with $\pi'$ (changing (6.4), (6.13), (6.11), and (6.15) for the induced left action), or $\tau$ with $\varpi$ flipped in $\mu$ (6.18)–(6.19), $A_\omega$ (6.22) and $T$ (7.3), which also conflicts with the draft's (10.6) (`C_tilt_repair.out` R1–R7; draft fn. 7 warns of conflicting conventions). | **Proven** (exact symbolic) |
| C2 | **Shiab operator** (reconstruction R-S: soldering invariants $\Phi^1,\Phi^2$). The degrees are consistent. It is gauge covariant only when tilted by $\varepsilon$; untilted, only for the centre $U(1)$. On Levi-Civita curvature, the operator with a commutator in the first term and $i$-anticommutators in the second gives $-\varepsilon^{-1}{*}G\,\varepsilon$ for $c=1$, the Einstein tensor of $Y^{14}$ (not $X^4$). The Weyl part is annihilated. | **Proven** for R-S |
| C3 | **Field equation.** $\circledast\circ d_B$ is not formally self-adjoint: there is an exact counterexample for every $c$ for each of the 8 bracket choices (`C_rs_all_members.out`), and $O(1)$ defects for all variants. So the printed EL equation (9.7) ("Swervature = Displasion") is not the variation of the printed action (9.4), and no rescaling fixes this. The cubic term also lacks the needed symmetry. | **Proven** for R-S; open beyond R-S |

**Empirical.** No claim is presently testable. The representation content predicts exotic spin-½ states (including charges $\pm\tfrac43,\pm\tfrac53,\pm2$, and colour sextets and octets) and spin-3/2 states, but no masses, couplings or breaking mechanism are given.

## 4. Proposed formalization order

This is the plan as written after Phase 1. All five items have since been done: item 1 in round 2 (D), item 2 in round 2 (F), items 3 and 4 in rounds 3–4 (G, H), and item 5 in `proofs/lean/geometric-unity/GUTilt.lean` (E). What remains open is listed in `reports/SYNTHESIS.md` §4.

1. **GU-P09: does a self-adjoint Shiab exist?** Decide by finite linear algebra over $\mathrm{Hom}_{\mathrm{Spin}(7,7)}(\Lambda^2\otimes\mathfrak u,\Lambda^{13}\otimes\mathfrak u)$, organised by Clifford grade. This settles whether any reading of (9.4)/(9.7) is consistent. It is the highest-value and lowest-cost item.
2. **Zorro metric (G06), then the 14→4 reduction (G14, GU-P10).** Write $g_\aleph$ explicitly. Compute $\gimel^*G(g_\aleph)$ by O'Neill formulas with the DeWitt fiber, and determine whether 4D Einstein gravity appears, and with what $\lambda$-dependent corrections.
3. **Full EL system** (G15) under repair R-C1′: the $\varepsilon$- and $\gimel$-variations, and the Noether identity for $\tau'(\mathcal H)$ as the candidate for $\Xi=\mathfrak D\Upsilon$.
4. **Fermion operator (G22)** and its self-adjointness. Then the deformation complex and its symbol and ellipticity (G16).
5. **Optional Lean 4 formalization** of Props. C2/C3/C3′. They are pure group theory on a semidirect product with a 1-cocycle and need no analysis.

Items 1–2 are well-posed now. Items 3–4 depend on 1.

## Round 2 — reduction (GU-G06, GU-G14)

Write-up: `formalization/geometric-unity/F-zorro-reduction.md`. Code: `computation/geometric-unity/F_{jets,fiber_curvature,zorro_reduction,zorro_einstein}.py`, outputs committed. Records: `metrics/geometric-unity/F-*.yaml`. All results are exact over $\mathbb Q(\lambda,c_H,c_V)$.

**Reconstruction.** No natural horizontal distribution exists on $\mathrm{Met}(X)$ (Lemma F0), so a connection on $TX$ must be supplied, here the Levi-Civita connection $\aleph$ of $g=\gimel$. The metric on $Y$ is
$$g_\aleph=c_H\,h_H(dx,dx)+c_V\,G_\lambda(Dy,Dy),\qquad Dy=dy-(\Gamma^Ty+y\Gamma)dx .$$
The draft supports two horizontal metrics:
- **Z1:** $h_H=y$, the point itself. This equals $\Theta_\aleph^*g_C$, the chimeric metric transported by the LC splitting. It is supported by (3.7), (3.15), p. 36 and p. 56.
- **Z2:** $h_H=g$, the base metric. This is a Riemannian submersion, supported by p. 22 "pulled up from the base".

**Results.**

| # | Statement | Status |
|---|---|---|
| F1 | The fiber $(GL(n)/O(p,q),G_\lambda)$ is the symmetric space $\mathbb R\times SL(n)/SO(p,q)$. Its curvature is independent of $\lambda$: $\mathrm{Ric}_F=-\tfrac n4G_{1/n}$ and $s_F=-n(n-1)(n+2)/(8c_V)$ ($=-9/c_V$ for $n=4$). | Proven |
| F2 | $s^*g_\aleph=c_Hg$; $s(X)$ is horizontal; $N_s=V$ isometrically. In Z1, $g_\aleph$ depends on $g$ only through $\aleph$, as the draft claims (p. 22). | Proven |
| F3 | **Z2** along $s(X)$: $s(X)$ is totally geodesic; $\mathrm{Ric}_{HH}=\mathrm{Ric}_g$, $\mathrm{Ric}_{HV}=0$, $\mathrm{Ric}_{VV}=\mathrm{Ric}_F$; $s_Y=R/c_H+s_F$. Hence $s^*G^Y=G(g)+\Lambda g$ with $\Lambda=9c_H/(2c_V)$ ($n=4$). The reduced action gives the same $\Lambda$, independent of $\lambda$. Under $(7,7)$ this has the anti-de Sitter sign. | Proven (O'Neill) |
| F4 | **Z1** along $s(X)$: $\mathrm{II}\neq0$ and $\propto c_H/c_V$; $\mathrm{Ric}_{HH}=\mathrm{Ric}_g+\frac{c_H}{4c_V(n\lambda-1)}g$; $\mathrm{Ric}_{VV}=-\frac{n+1}4G_{1/(n+1)}$; $s_Y=R/c_H+\sigma_1(\lambda)/c_V$. The cosmological term depends on $\lambda$ and also on the route. For $n=4$: projection gives $\Lambda_P=\frac{3c_H}{2c_V}\frac{4-15\lambda}{1-4\lambda}$ and the reduced action gives $\Lambda_A=\frac{5c_H}{4c_V}\frac{5-18\lambda}{1-4\lambda}$. | Proven ($\mathrm{II}$, $s_Y$, and the Ricci blocks for every $n$: K, 2026-10-30) |
| F6 | $g_\aleph$ is Einstein along $s(X)$ **iff** Z1, $\lambda=1/(n+1)$ and $\mathrm{Ric}(g)=0$. Then the 14D Einstein equation with $\Lambda_Y=-15/(2c_V)$ is exactly 4D vacuum GR. But $\lambda=\tfrac15<\tfrac14$ gives signature $(8,6)$ or $(10,4)$, not $(7,7)$. With the draft's $(7,7)$, and in Z2 for every $\lambda$, the condition never holds. It also fails off $s(X)$. For flat $\aleph$, Z1 at $\lambda=1/(n+1)$ is locally the symmetric space $SL(n+1)/SO$; for $n=4$ this is $SL(5)/SO(3,2)$ (Prop. F7). | Proven (Z1 half via K, 2026-10-30; also checked on a Weyl base at $n=4$); F7 proven |
| F8 | Under R-S/C5, at $T_\omega=0$ the pulled-back $*\Upsilon$ is exactly $-(G+\Lambda_Pg)\otimes\gamma$, a precise version of (9.10). Still $\Upsilon\neq0$, because its vertical–vertical block is non-zero. $T_\omega=t\Phi^1$ adds a 14D cosmological term, so F6 decides solvability. | Conditional; constants open |

**Verdict on GU-G14.** Under R-Z, 4D Einstein gravity with a cosmological term does arise from $Y$, in the projection and reduced-action senses. Its $\Lambda$ is fixed by the fiber geometry: $|\Lambda|\sim1/\ell_F^2$, with no small parameter. Matching the observed value needs $\ell_F\approx2\times10^{26}$ m.

The answer depends on three choices the draft leaves open:
- the horizontal metric (Z1 or Z2);
- $\lambda$ (Z1 only);
- the readout (projection, reduced action, or the field equation on $Y$).

An exact 4D vacuum reduction exists only at $\lambda=1/(n+1)$, which is incompatible with the draft's $(7,7)$ choice. Claims: GU-006 C, GU-021 C, GU-024 V, GU-054 C (previously O).

## Round 2 — Shiab operators and the field equation (GU-P09, GU-G10, GU-G13)

Write-ups: `formalization/geometric-unity/D-selfadjoint-shiab.md` and `E-lean-tilt.md`. Code: `computation/geometric-unity/D_{rep_counts,selfadjoint_shiab}.py`, outputs committed. All results are exact (integer characters; integer or dyadic Clifford data). Lean: `proofs/lean/geometric-unity/GUTilt.lean`. Records: `metrics/geometric-unity/D-*.yaml`.

**Hypotheses.** $\mathcal S=\mathrm{Hom}_{\mathrm{Spin}(7,7)}(\Lambda^2\otimes\mathfrak g,\Lambda^{13}\otimes\mathfrak g)$ (Def. 8.1), with:
- $\mathfrak g$ one of $\mathfrak u(64,64)$, $\mathfrak{so}(64,64)$, or the subalgebras $\mathfrak{so}(8,8)$ and $\mathfrak{so}(7,8)$;
- an $\mathrm{Ad}$-invariant pairing;
- $\varpi$ varied with $\varepsilon$ fixed, as in (9.7).

| # | Statement | Status |
|---|---|---|
| D1 | $\dim\mathcal S=200$ and $\dim\mathcal S_{\rm sa}=48$ for $\mathfrak u(64,64)$. For $\mathfrak{so}(64,64)$ the values are 52 and 12. For $\mathfrak{so}(7,7)$, $\mathcal S=0$, so no Shiab operator exists. $[\Lambda^i\otimes\mathfrak u]^{\mathrm{Spin}(7,7)}$ is 2-dimensional for every $i$, so Def. 8.1's 15-element basis is half of the 30. | Proven (characters; brute force in $D_2$, $D_3$) |
| D2 | $\circledast\circ d_B$ is formally self-adjoint iff $\circledast=\Psi\wedge$ with $\Psi\in\Lambda^3\otimes S^2\mathfrak g^*$. This holds for any structure group. | Proven |
| D3 | Every $\circledast\in\mathcal S$ kills Weyl (Schur), so the draft's Weyl claim selects nothing. Every $\circledast\in\mathcal S$ sends Levi-Civita curvature to $x\,\mathrm{Ric}+y\,\mathrm{Scal}\,g$. Self-adjointness forces $y=-x/2$, the exact Einstein tensor of $Y^{14}$. $W_A=[*\Phi^3\wedge\cdot]_A$ (draft form (8.1), $r=11$) and $W'_A$ are self-adjoint with $W_A(F_R)=-*G$. Among the 8 single-term members of (8.1), the self-adjoint combinations are exactly $\mathrm{span}(W_A,W'_A)$. **GU-P09 as stated is false.** | Proven (exact) |
| D4 | Vary $\varpi$ at fixed $\varepsilon$, with $\circledast$ pointwise and independent of $\varpi$, and a constant Ad-invariant pairing $b$. Then $\Upsilon=\circledast F_A+\kappa_1*T$ is a $\varpi$-gradient only if $b(\alpha\wedge\circledast\xi)=\sum_k\Omega_k\wedge I_k(\alpha\wedge\xi)$, with $\Omega_k$ 11-forms and $I_k$ the invariant symmetric forms. For $\mathfrak u(64,64)$ and $b=\mathrm{Tr}$ this is $\circledast\xi=\Omega_0\wedge\xi+\Omega_1\wedge\mathrm{Tr}(\xi)\mathbf 1$ (generalized Chern–Simons). Such an operator has no Einstein output, whatever its equivariance. If moreover $\circledast$ is $\mathrm{Spin}(7,7)$-equivariant, or $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$-equivariant, then $\circledast=0$. | Proven (Helmholtz); exact illustrations and controls |
| D5 | For $\circledast\in\mathcal S_{\rm sa}$, the $\varpi$-EL equation of (9.4) is $\circledast F_A+\kappa_1*T-\tfrac13h(T)=0$. For $W_A$, $h(T)=[T\wedge[T\wedge\Phi_{11}]_C]_A\neq0$. | Proven; exact check |
| E | Props. C1–C3′ are machine-checked in Lean 4 core: 42 theorems, no `sorry`. These include the $-2c(h)$ residual, "printed Lemma 7.2 holds at $h$ iff $2c(h)=0$", and the $\tau'$ repair. The model is abstract (a group, a module and a 1-cocycle); L-C1 is an axiom of the model. | Proven (Lean) |

**Verdict on GU-P09 and GU-G13.** The quadratic part of (9.4) can be made variational in a way compatible with the Einstein and Weyl properties. The field equation (9.7) cannot be made variational with any nonzero admissible Shiab operator, any coefficients, or any invariant pairing. The obstruction is the $T\wedge T$ term of $F_A$, and it is equivalent to the absence of a $\mathrm{Spin}(7,7)$-invariant 11-form. GU-G13 is now proved for the whole admissible space. Remedy (a) is closed in the negative; remedies (b) (take $\Upsilon=0$ as the definition) and (c) (use D5's equation) remain.

**Claims.**
- GU-045: Def. 8.1 undercounts the invariants.
- GU-047 C: the Weyl part is automatic, and the Einstein ratio is forced by self-adjointness.
- GU-053 F: for all of $\mathcal S$.

Round-1 Prop. C6′'s criterion is superseded; see `reports/CORRECTIONS.md`.

## Round 3 — full EL system, fermion operator, deformation complex (GU-G15, GU-G22, GU-G16)

Write-up: `formalization/geometric-unity/G-full-el-fermions.md`. Code: `computation/geometric-unity/G_{el_noether,fermion,symbols}.py`, outputs committed. Lean: `proofs/lean/geometric-unity/GStueckelberg.lean` (24 theorems, no `sorry`). Records: `metrics/geometric-unity/G-*.yaml`. Everything is exact (integers, rationals, Gaussian integers).

**Reconstruction.**
- Bosons: R-C1′ and the D.1 hypotheses, with ⊛ arbitrary for G1–G3 and $W_A$ for the Einstein case.
- Fermions: (9.16) is read as $\mathcal D=\begin{psmallmatrix}{*}\circledast^S d_A& d_A\\-d_A^*&0\end{psmallmatrix}$ on $(\zeta,\nu)$, tilted by $\rho(\varepsilon)$, with $\circledast^S=\Phi_{11}\wedge$ acting by Clifford multiplication. The alternatives are the other Spin(7,7)-equivariant maps.
- Source observation: the draft's own (10.6) uses the repaired sign $\tau'$ (C3′), not the printed (6.4).

| # | Statement | Status |
|---|---|---|
| G1 | **Stueckelberg factorisation.** $\tau'(\mathcal H)$ acts freely, with orbits equal to the fibres of $\pi'$. $I_1^B=J\circ\pi'$, where $J(S)$ is (9.4) at $\varepsilon=1$ with the untilted ⊛. This holds for every ⊛, every pairing and every set of coefficients. The "bought-back" gauge invariance makes ε pure gauge: the theory is one field $S\in\Omega^1(\mathrm{ad})$ on the fixed $(Y,g_\aleph,\nabla^0)$. | Proven; Lean core; exact density check |
| G2–G3 | **Complete EL system.** The system is $\Upsilon=0$, with the off-shell Noether identity $\Xi\equiv d_{A_\omega}\Upsilon_\omega$. So (9.5)–(9.6) hold with $\mathfrak D=d_{A_\omega}$, and D5's equation is the full EL system. They fail for the printed (9.7), where $\Xi\ne d_A\Upsilon^{\rm pr}$. | Proven; exact Euler-operator check in a non-abelian $d=3$ model (incl. non-self-adjoint and Chern–Simons control) |
| G4 | **Bianchi property** (fn. 10). $d_\nabla(\circledast F_\nabla)\equiv0$ for all metrics iff the Einstein ratio holds. This cuts out 198 of the 200 dimensions and contains $\mathcal S_{\rm sa}$. At $T=0$ the Noether identity is the contracted Bianchi identity of $g_\aleph$. | Proven |
| G5 | Every equivariant spinor Shiab **flips chirality**. The natural one is the Rarita–Schwinger operator. The chirality labels of (9.16) are inconsistent, and a one-swap relabelling repairs them. | Proven; exact |
| G6–G7 | **Fermionic symbol criterion** (D2 with the opposite spinor symmetry). In (7,7), $\Lambda^3\subset\mathrm{Sym}^2S^\pm$ and $\Lambda^3\not\subset\Lambda^2S$. So **no Majorana(–Weyl) action has a Shiab/RS kinetic term for ζ**. The MW reading of B4 ($3\times16\oplus144$) leaves ζ's transverse components non-propagating (symbol rank 128/960). The minimal real repair (two MW copies ≅ complex Weyl) doubles the content to $6\times16\oplus2\times144$. | Proven; exact |
| G8 | The Hermitian $U(64,64)$ readings (complex Weyl, Dirac) are consistent with the printed phases. Only $\mathfrak u(64,64)_{\rm ev}\cong\mathfrak{gl}(64,\mathbb C)$ couples to Weyl fermions. | Proven; exact |
| G9 | The fermion symbol is invertible iff $\eta(\xi,\xi)\ne0$: rank 1920/1920 off the null cone, 1024/1920 on it, and 1920/1920 in (14,0). So the fermion operator is **elliptic in Euclidean signature**. | Proven + exact ranks |
| G10 | $\delta_2\circ\delta_1=\mathrm{ad}(\Upsilon)$, which confirms (9.22) in this sense. The draft's (10.9) is the linearisation of the printed Υ. The complex reduces to $\Omega^1\xrightarrow{\mathcal L}\Omega^{13}$. | Proven |
| G11 | **The bosonic deformation complex is not elliptic.** In (7,7) this holds for every self-adjoint Shiab, because of a null-cone obstruction. For the whole self-adjoint family of (8.1), $\mathrm{span}(W_A,W'_A)$ (via a block-rescaling lemma), and for R-S, the symbol sequence is exact at no covector, for any background, in (7,7) **and (14,0)**. The exact rank of $\sigma(\circledast\circ d)$ is at most 122877 of the 212992 needed on $\mathfrak u(64,64)$, and 60580 of the 105664 needed on $\mathfrak{so}(64,64)$. | Proven + exact ranks |

**Verdicts.**
- GU-G15 is closed for $(\varepsilon,\varpi)$: the ε-equation is redundant, as the draft says, but only for the true Υ. The ℷ-variation is formal, with no fibre decay, and stays open.
- GU-G22 is closed under R-F. The printed chirality labels fail, and the B4 Majorana–Weyl reading is incompatible with the ζ kinetic term.
- GU-G16 is closed in the negative. The "elliptic deformation complex in Euclidean signature" (§12 (x)) fails for every Shiab of the draft's self-adjoint family and for R-S. The fermion block alone is elliptic in (14,0).
- Claims: GU-052 O→C; GU-056 O→C/F; GU-057 split (V for $\delta_2\delta_1=\mathrm{ad}\Upsilon$, F for ellipticity); GU-072 F (bosonic) / C (fermionic).
- Earlier conclusions qualified: see `reports/CORRECTIONS.md` (B4/GU-G07; GU-052; GU-057/072).

## Round 4 — Euclidean ellipticity over all self-adjoint Shiab operators (GU-072, GU-057, GU-G16)

Write-up: `formalization/geometric-unity/H-ellipticity.md`. Code: `computation/geometric-unity/H_ellipticity.py`, output committed. It imports `G_symbols.py`, `D_rep_counts.py` and `shiab_checks.py` without modifying them. Records: `metrics/geometric-unity/H-{euclidean-ellipticity,twisted-shiab}.yaml`. Everything is exact (integers, rationals, Gaussian rationals; exact sparse ranks).

**Question.** Round 3 (G11) left open whether *some* $\circledast\in\mathcal S_{\rm sa}$ (48-dimensional for $\mathfrak u$, 12 for $\mathfrak{so}$) has an elliptic deformation complex in (14,0). The round-3 negative covered only $\mathrm{span}(W_A,W'_A)$ and R-S.

**Hypotheses.**
- *Euclidean real form:* $\mathrm{Spin}(14)$ and $\mathfrak g=\mathfrak u(128)$ (the unitary algebra of the positive spinor Hermitian form). The answer is the same for $\mathfrak u(64,64)$ and for any real form of $\mathfrak{gl}(128,\mathbb C)$, since the conditions are polynomial and real forms are Zariski dense.
- *Complex:* as in G10.
- Any pairing, κ, background $(\varepsilon,T)$, true or printed Υ.

| # | Statement | Status |
|---|---|---|
| H1 | **Reduction.** Exact at ξ iff $\kappa\ne0$ and $\beta_\xi(\alpha,\alpha')=\Psi(\xi\wedge\alpha\wedge\alpha')$ is nondegenerate on $\xi^\perp\otimes\mathfrak g$. This is independent of the background, of $b$ and of the Υ version. One ξ suffices. | Proven |
| H2–H3 | **Exact basis.** The basis of $\mathcal S_{\rm sa}$ is the 48 Clifford trace invariants $O_{(j,l)}=\mathrm{Tr}(q(\omega)(X_jY_l+Y_jX_l))$ and $E_{(j,l)}$ (with Γ). They agree block by block with the $D_7$ characters, and every block pairs even with odd grades. Coordinates: $W_A$ is $-2$ on 12 $O$-blocks and $W'_A$ is $2$ on 12 $E$-blocks. $W_C\circ\mathrm{Ad}_\Gamma$ and $W'_C\circ\mathrm{Ad}_\Gamma$ fill the other 24 blocks. These reproduce G_symbols' symbol matrix entrywise. | Proven; exact |
| H4 | **Spin(13) reduction.** $W\otimes\mathfrak g=4\Lambda^0+8(\Lambda^1..\Lambda^6)+4(H_1..H_6)$. β splits into 13 matrices $C_\rho$ ($2\times2$ or $4\times4$) with entries linear in the coordinates. The Schur-predicted rank equals the exact $212992^2$ sparse rank at 6 points. This includes $W_A$, $W'_A$ and $W_A+W'_A$ at 122850, the round-3 values. | Proven; exact |
| H5 | **Main.** For $\mathfrak u(128)$ (any real form of $\mathfrak{gl}(128)$) the complex is elliptic in (14,0) **iff $\kappa\ne0$ and $\Delta(\Psi)\ne0$**. Here Δ is the product of 7 quadrics and 5 quartics (degree 34), all irreducible over $\mathbb Q(i)$ and sign-changing on the real form. So the elliptic set is **open and dense** in $\mathcal S_{\rm sa}$. $\mathrm{span}(W_A,W'_A)$ lies entirely in $\{\Delta=0\}$. | Proven; exact |
| H6 | **Ellipticity forces Einstein.** $Q_1\ne0$ forces a nonzero $O_{(1,2)}$ or $E_{(2,13)}$ coefficient, which are the only coefficients curvature sees. So every elliptic member reproduces the Einstein tensor of $Y$ with the exact ratio (D3.3). **No elliptic member is of the draft's form (8.1).** | Proven |
| H7 | **Minimal extension of (8.1).** $xW_A+yW'_A+zW_C\mathrm{Ad}_\Gamma+wW'_C\mathrm{Ad}_\Gamma$, with $\mathrm{Ad}_\Gamma$ the grade involution and real coefficients, is elliptic iff $\kappa\ne0$ and $xz+yw\ne0$. The quartic determinants are $a(x^2+y^2-z^2-w^2)^2+b(xz+yw)^2$ with $ab>0$. Example: $\circledast_\sharp=W_A+W_C\circ\mathrm{Ad}_\Gamma$ gives rank $\sigma\delta_2=229376/229376$ and $\circledast_\sharp(F_R)=W_A(F_R)$. | Proven; exact |
| H8 | **$\mathfrak{so}(64,64)$ / $\mathfrak{so}(128)$: never elliptic.** $\mathrm{rank}\le60580<105664$ for all 12 dimensions, attained by $W_A$. The corank $45084=45058+26$ comes from unequal odd/even Spin(13)-multiplicities ($\mathfrak{gl}(64)$ vs $\Lambda^2S^+\oplus\Lambda^2S^-$) plus an uncoupled grade-14 copy. | Proven; exact |
| H9 | **(7,7).** Generic Ψ is exact at non-null ξ, so the null cone (G11.1) is the only obstruction. The pencil determinants $\propto(x-y)^2(x+y)^2$ explain G11's degenerate lines $W_A\pm W'_A$. | Proven; exact |

**Verdicts.**
- §12 (x) "elliptic deformation complex in Euclidean signature" holds for a **generic** self-adjoint Shiab operator on $\mathfrak u$, with $\kappa\ne0$. It fails for every operator the draft actually writes: the self-adjoint (8.1) family, R-S, and anything on $\mathfrak{so}(64,64)$. In (7,7) it fails always.
- Elliptic, self-adjoint and Einstein are compatible, but only outside (8.1), and nothing in the draft selects such a member.
- Unchanged: (9.7) is still not variational (D4); ellipticity refers to the complex of the true EL equation (D5).
- Claims:
  - GU-072 bosonic: F → **C**;
  - GU-057 ellipticity: F → **C**; $\delta_2\delta_1=\mathrm{ad}\Upsilon$ stays V;
  - GU-G16: the negative verdict is restricted accordingly.
- Changes are logged in `reports/CORRECTIONS.md`. G11's numbers are all reproduced and unchanged.

## Round 5 — the elliptic twist and the 14→4 reduction (GU-054, GU-G14, GU-G16)

Write-up: `formalization/geometric-unity/I-elliptic-reduction.md`. Code: `computation/geometric-unity/I_elliptic_reduction.py`, output committed, exact, ~100 s. Record: `metrics/geometric-unity/I-elliptic-reduction.yaml`. Revised the same day after `reports/reviews/REVIEW3-GU-2026-10-30.md` (independent check `computation/geometric-unity/xcheck_round5.py`).

Round 4's elliptic operator $\circledast_\sharp=W_A+W_C\circ\mathrm{Ad}_\Gamma$ differs from $W_A$ off curvature. The question was whether this changes the reduced equations, since $T\neq0$ on solutions near $s(X)$.

- **Which parts it sees (I1).** The twisted pair never sees $F_0$, or any algebraic curvature tensor. On $\mathfrak{so}(14)$-valued 2-forms it sees exactly the first-Bianchi defect and the pair-antisymmetric parts.
- **Parity (I2).** For every self-adjoint Shiab operator, $\Upsilon_{\rm ev}=\circledast(d_0T_o)+\kappa{*}T_e+2B(T_e,T_o)_{\rm ev}$, and $B(T_e,T_o)$ is even, so $\Upsilon_{\rm odd}$ has no cross term. For $\kappa\neq0$, $T_e$ is slaved to $T_o$ where $L_{T_o}$ is invertible ($|T_o|$ small compared with κ). If also $(x_1,x_{13})\neq0$, the odd part $T_o$ cannot vanish on any open set meeting the fibre-scaling cone through $s(X)$; so a pure contorsion is not a solution on any such open set.
- **Twist identity (I3).** The twisted cubic term is $2\Psi(T_o\wedge T_e\wedge T_e)$. On odd fields its quadratic contribution to Υ therefore vanishes; the $T_o\wedge T_o$ part of $F_A$ is cancelled exactly by $-h/3$. So for $T_e=0$ the algebraic ("Einstein") part of $\Upsilon=0$ is the same for the whole elliptic family as for $\mathrm{span}(W_A,W'_A)$ (I4). The twist enters only the first-order constraint $\circledast(d_0T)=0$.
- **Constraint system (I5, I6).** In the soldering and invariant sectors (one 2×2 block) this constraint has determinant $4i(xz+yw)$ in (14,0), up to a constant. That is H7's ellipticity factor, as expected: there the constraint is the principal symbol of $\circledast\circ d$ (H4's blocks $C_{\Lambda^0}$, $C_{H_1}$). For $\circledast_\sharp$ it forces the Γ-partner to be closed ($du=0$, $d^\nabla\theta'=0$); $W_A$ alone forces only $(d^\nabla\theta')_{\Lambda^3}=0$, and not $du=0$. In the invariant sector full-neighbourhood solution sets coincide by the 14D Bianchi identity. Projected solution sets differ: the 4D identity controls only tangential derivatives, whereas the twist also excludes normal derivatives of $u$ (I6.1, `projection_regression.py`). With general $\theta'\neq0$ full-neighbourhood equivalence was not checked.
- **Exact Υ in the two sectors.**
  - Invariant sector: $\Upsilon=-{*}([xG-c_1g]\,e\otimes\gamma+[yG-c_{13}g]\,e\otimes\gamma\Gamma)$ with $c_1=\kappa t+x(312t^2\pm104u^2)\pm208\,ytu$; the upper signs hold in (14,0), the lower in (7,7). This also fixes F8's open constants.
  - Symmetric displacement sector: $\Upsilon=0\iff G(g_\aleph)+G(R_\theta)=\kappa\theta/x$ and θ is Codazzi, where $R_\theta=4(\theta_{ac}\theta_{bd}-\theta_{ad}\theta_{bc})$.
- **Reduction (I7), Z1 and Z2 alike.**
  - (P) is underdetermined: it admits $T=0$ ($\Lambda_4=\Lambda_P$) and, in the invariant sector, 4D Einstein with $\Lambda_4=\Lambda_P-c_Hc_1/x$, a **free constant**. For $\circledast_\sharp$, $\Lambda_4$ takes every real value in (7,7), and exactly the values $\ge\Lambda_P-c_H\kappa^2/1248$ in (14,0) (the bound is $\Lambda_P-c_H\kappa^2/(1248x'^2)$ for $x'W_A+z'W_C\mathrm{Ad}_\Gamma$, $x'>0$). Only $\kappa t$ comes from $T_\omega$; with the printed (9.7), (7,7) gives only $\Lambda_4\le\Lambda_P+c_H\kappa^2/1248$.
  - (P), symmetric sector, algebraic part: Einstein plus a displacement stress, quadratic and linear in $h=s^*\theta$. The projected constraint was not computed.
  - (E), which solutions of $\Upsilon=0$ near $s(X)$ satisfy (for these I2.1 gives $T_o\neq0$): for $(x,y)\neq0$, no solution in (7,7) in the invariant sector. In (14,0) a solution needs Z1, $\lambda=\tfrac15$, $\mathrm{Ric}(g)=0$ and $c_V\kappa^2\ge9360$, and has $\Lambda_4=0$. In the symmetric sector (E) is solvable pointwise for $|G^Y|\ll\kappa^2$ (on $s(X)$ this needs $c_V\kappa^2\gg1$), so it imposes no *pointwise* 4D condition. The HHH constraint requires $s^*\theta$ to be Codazzi, i.e. $G(g)$ Codazzi to leading order in $1/\kappa$.
  - **The 4D algebraic equations and cosmological ranges agree for $\circledast_\sharp$ and $W_A$. Their projected normal-derivative constraints differ.**
- **Verdict.** The elliptic twist preserves the 4D algebraic content in these sectors while changing which normal extensions satisfy the projected equation. The source gives no rule fixing $T$ off $s(X)$ and no choice of reading; those choices remain necessary to assess the reduction. G.2.5's "any open set" is narrowed to open sets meeting the fibre-scaling cone through $s(X)$. That still excludes compact support along the fibres. For a real-analytic base metric the original statement holds (`reports/CORRECTIONS.md`).

## Round 6 — the slaved even sector (GU-054, GU-G14, GU-G16)

Write-up: `formalization/geometric-unity/J-slaved-sector.md`. Code: `computation/geometric-unity/J_slaved_sector.py`, output committed, exact, ~10 min (567 s; a second run 560 s, identical apart from the timing lines) on one core. Record: `metrics/geometric-unity/J-slaved-sector.yaml`. Revised the same day after `reports/reviews/REVIEW4-GU-2026-10-30.md` (independent check `computation/geometric-unity/xcheck_round6.py`, 90/90).

Round 5 settled the reduction only for $T_e=0$, where the twisted coefficients $z,w$ drop out (I7(e)). The question was whether they change the reduced equations once the even equation forces $T_e\neq0$ (Cor. I2.2).

- **Sector.** $T_o=t\Phi^1+u\Phi^1\Gamma$ with nonconstant $t,u$, and $T_e$ in the four copies of the vector representation in $\Lambda^1\otimes\mathfrak g_{\rm ev}$ (grades 0, 2, 12, 14). The even equation closes exactly on this sector. It is a 4×4 linear system $M\vec A=-S(dt,du)^\top$, so the slaved field is explicit, $T_e=-M^{-1}S\,(dt,du)$ (J1). $z,w$ enter both the source $S$ and $M$.
- **Odd equation (J2).** On this sector it is second order in $(t,u)$:
  - two Einstein-type channels, $-x_iG+c_ig+\sum_\phi H_{i\phi}(\nabla^2\phi-\Box\phi\,g)+(\text{gradient terms})=0$ for $i=1,13$;
  - $\Lambda^2$-type equations $\varpi\,dt\wedge du=0$.

  Every coefficient is exact. The linear maps were checked on full bases, the bilinear ones on all basis pairs, and the assembled formula end to end on a random 2-jet.
- **Where the twist enters (Cor. J2.1).**
  - It enters the principal (Hessian) coefficient matrix $H=-\ell M^{-1}S$ ("kinetic matrix" in the note), to leading order $\tfrac{576}{\kappa}$ times a Gram matrix. In the $\mathfrak u(128)$ real form this is $\mathrm{Gram}((x',w'),(y',-z'))$, semidefinite of the sign of κ and definite iff $x'z'+y'w'\neq0$: H7's ellipticity factor, by construction. Exactly, for every member, $\det H=576^2(xz\pm yw)^2\kappa^2/\det M$ (block structure plus Jacobi's identity; added after review, R6-2).
  - It enters the $O(\kappa^{-2})$ gradient terms, and the $\Lambda^2$ coefficients, which are all $\propto xz\pm yw$ at leading order.
  - It does not enter the zeroth-order term $c_i(t,u)\,g$.
  - For $\circledast_\sharp$, $\det H=576^2\kappa^2/\det M\neq0$; for $W_A$, $H$ has rank 1 and only $t$ carries second derivatives. Both have a Λ²-constraint: $dt\wedge du=0$ at every point where $\Upsilon=0$ and some $\varpi\neq0$ — for $\circledast_\sharp$ near $T_o=0$, for $W_A$ (exactly, at $O(\kappa^{-3})$) wherever $u\neq0$. (Corrected after review 2026-10-30, R6-5; the text had "$dt\wedge du=0$ is forced on open solution sets" for $\circledast_\sharp$ only; see `reports/CORRECTIONS.md`.)
  - **So $\circledast_\sharp$ and $W_A$ give different reduced equations once $T_e\neq0$**: they differ in the rank of $H$; the Λ²-constraint distinguishes them only at leading order and at $u=0$.
- **Reduction along $s(X)$.**
  - *Basic scalars* ($t,u$ pulled back from $X$):
    - (P) in Z2 gives a 4D system in $(g,\tau,\upsilon)$: an Einstein-type equation with Hessian terms $H_{i\phi}(\nabla^2\tau_\phi-\Box\tau_\phi g)$ and a field-dependent cosmological term, for $y=0$ (including $\circledast_\sharp$) a curvature-free tensor constraint on the scalars, and $\varpi\,d\tau\wedge d\upsilon=0$. Its solutions are not classified. (Corrected after review 2026-10-30, R6-4; the text had "4D Einstein gravity coupled to two scalars"; see `reports/CORRECTIONS.md`.)
    - (P) in Z1 contains $H(d\tau,d\upsilon)=0$, so the scalars are constant wherever $\det H\neq0$, i.e. for every member with $xz\pm yw\neq0$ ($\det M\neq0$), and also for $W_A$: round 5 again.
    - (E) gives exactly round 5's solutions (J3) for every member with $xz\pm yw\neq0$ ($\det M\neq0$) and for $W_A$; for the other members with $xz\pm yw=0$ it is not shown. (Corrected after review 2026-10-30, R6-10; see `reports/CORRECTIONS.md`.)
  - *Isotropy-invariant 2-jet along $s(X)$:* (E) at the points of $s(X)$ reduces to six scalar equations and a balance identity.
    - It has solutions in **(7,7) and (14,0), in Z1 and Z2, for every $\lambda\neq\tfrac14$ on the side the signature requires** ($\lambda<\tfrac14$ for definite (14,0), $\lambda>\tfrac14$ for (7,7)), with $T_e=0$ on $s(X)$ but $\nabla T_e\neq0$. Exact samples exist for $\circledast_\sharp$ and $W_A$.
    - On that branch ($a=0$) the 4D cosmological constant $\Lambda_g=\tfrac{c_H}{11}(\Theta_0-c_1/x)$ is free: every real value in (7,7) for $\circledast_\sharp$, a half-line in (14,0) ($c_H>0$) (J4).
    - With a normal gradient $a\neq0$ (still the invariant 2-jet) the six equations reduce to $c_i=x_i\Theta+13\beta(Ha)_i-|\tau|^2\mathfrak B_i(\alpha)$ (J4(c′)), so $\Lambda_g$ then depends on $(z,w)$. For $\circledast_\sharp$ in (14,0) the full range is the same half-line in Z2, but in Z1 it extends down to $\tfrac{c_H}{11}(\Theta_0-\tfrac{\kappa^2}{1248}-\tfrac{39}{8c_V(1-4\lambda)})$, attained exactly (J4(d′); after review, R6-1).
    - So round 5's "(E): none in (7,7), $\Lambda_4=0$ in (14,0)" is a property of $T_e\equiv0$.
  - *Radial ansatz* $t=f(s)$, $u=h(s)$, with $s$ the fibre-scaling coordinate, in Z2:
    - (E) at $s(X)$ forces $\Lambda_g=-c_H/c_V$ (in (14,0) for $\circledast_\sharp$ this needs $c_V\kappa^2\ge8112$);
    - with $\kappa\neq0$ and $T_e$ the slaved field in the vector sector ($\det M\neq0$), $\Upsilon=0$ on a neighbourhood of a point of $s(X)$ has **no solution in either signature of the setting, (14,0) or (7,7)**, for every member with $(x,y)\neq0$ and every $(z,w)$. The proof uses a parallel $ds$, Einstein leaves and O'Neill's A-tensor; in indefinite signature two traces of the order-$\varepsilon^2$ leaf condition (the second symmetrised in $h$) give $4C_2^{\rm sym}-C_1=32\Lambda_g^2$, so $\Lambda_g=0$, a contradiction (`J5d-indefinite-leaf.md`). (Extended 2026-10-30 by J5d from definite signature to (7,7); its trace step corrected after review, R7-1, R7-2; see `reports/CORRECTIONS.md`.) (Corrected after review 2026-10-30, R6-10; the hypotheses had been omitted; see `reports/CORRECTIONS.md`.)
    - In particular there is none in (7,7), from a base (1,3), (3,1) or (2,2) (J5d; this settles J5(d), which had left it open). The signatures (8,6), (10,4) of Z2 are outside the setting of J1–J2 and not covered.
- **Verdict.** The twist does change the reduced equations once $T_e\neq0$. It acts on the derivative structure and never on the zeroth-order term $c_i(t,u)\,g$; the effective $\Lambda_4$ of the 2-jet reduction depends on $(z,w)$ once $a\neq0$. (Corrected after review 2026-10-30, R6-3; the text had "never on the cosmological term"; see `reports/CORRECTIONS.md`.)
  - The physically relevant output still depends on the undeclared rule for $T$ off $s(X)$. With free second normal derivatives, (E) at $s(X)$ becomes as underdetermined as (P). A genuine neighbourhood condition can exclude solutions altogether (radial Z2, (14,0) and (7,7); in Z2 also without the radial ansatz for $\circledast_\sharp$, $W_A$, every elliptic member in (14,0) and every member with $xz-yw\neq0$, $x\neq\pm y$ in (7,7), round 8).
  - No statement made before round 6 is contradicted: round 5's (E) statements are all scoped to the invariant sector, which has $T_e\equiv0$ by definition (confirmed by REVIEW4-GU, R6-11). This round's own J4(d) range was too broad and is corrected (R6-1, `reports/CORRECTIONS.md`).
  - Not covered:
    - $T_o$ outside the invariant sector;
    - $T_e$ outside the vector sector for large $|T_o|$;
    - neighbourhood solutions in Z1 beyond Round 10 (limit points of Einstein pieces without analyticity; the members left open by L.3); in Z2 beyond the radial ansatz, only the members with $xz\pm yw=0$, $x^2\pm y^2\neq0$, $(z,w)\neq0$ (up to symmetry $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, $c\neq0$, $c^2\neq\tfrac{11}6$ in (14,0)), away from $s(X)$ (Rounds 8, 9);
    - reading (A).

## Round 8 — $\Upsilon=0$ beyond the radial Z2 ansatz (GU-054, GU-G14, GU-G16)

Write-up: `formalization/geometric-unity/L-nonradial-z2.md`. Code: `computation/geometric-unity/L_nonradial.py`, output committed, exact, 73 asserted checks, ~6.5 min (380 s) on one core; a second run was identical apart from the timing lines. Record: `metrics/geometric-unity/L-nonradial-z2.yaml`. Refereed in `reports/reviews/REVIEW6-GU-2026-10-30.md` (no errors or gaps; two summary overstatements and four presentation items corrected; independent check `computation/geometric-unity/xcheck_round8.py`, 70/70).

Hypotheses of J5 without the radial ansatz: Z2, $n=4$, $\lambda\neq\tfrac14$, the real forms (14,0) and (7,7), $\kappa\neq0$, $t,u$ arbitrary, $T_e$ the slaved field in $\mathcal V$ with $\det M\neq0$, $U$ connected and open.

- $dt\wedge du=0$ on every open solution set, for every member with some $\varpi\not\equiv0$ (L1; this includes $\circledast_\sharp$, $xW_A$ and every member with $xz\pm yw\neq0$).
- Z2 has no open set with $\mathrm{Ric}\in\mathrm{span}(g,ds\otimes ds)$; in particular it is Einstein on no open set (L2). A concircular function on Z2 has parallel gradient (L3), so the concircular branch reduces to excluded cases (L4).
- Where the curvature-free channel $yE^{(1)}-xE^{(13)}$ has no Hessian term, its residual $R'=-1916928\kappa^3(x^2\pm y^2)(xz\pm yw)^2(ty-ux)$ forces $(t,u)\parallel(x,y)$; this line is invariant iff $xw=yz$ (L5).
- On the line one equation remains, $\nabla^2r=-\tfrac7{18}r\,\mathrm{Ric}$, with $Er=0$ (L6). It fails at order $\varepsilon^2$ in the fibre at a point of $s(X)$ (L7).
- **Result.** No solution on any open subset of $Y$ for $\circledast_\sharp$, $W_A$ and its multiples, every elliptic member in (14,0), and every member with $xz-yw\neq0$, $x\neq\pm y$ in (7,7) (L8, L9).
- **Open.** Z1; the other non-elliptic members; $x=\pm y$ in (7,7); the other sectors.
- No statement made before round 8 is contradicted.

## Round 9 — the members round 8 left open (GU-054, GU-G14, GU-G16)

Write-up: `formalization/geometric-unity/M-z2-remaining-members.md`. Code: `computation/geometric-unity/M_remaining.py`, output committed, exact, 35 asserted checks, ~4.5 min (270 s) on one core; a second run was identical apart from the timing lines. Record: `metrics/geometric-unity/M-z2-remaining-members.yaml`. Refereed in `reports/reviews/REVIEW7-GU-2026-10-30.md` (independent check `computation/geometric-unity/xcheck_round9.py`, 47/47; one side claim corrected, R9-1, which settles the member $c^2=\tfrac{11}6$ in (14,0) on every open set).

Hypotheses of round 8; for $(x,y)=0$ the same without $(x,y)\neq0$.

- Symmetries of the reduced system: a rotation (14,0) / boost (7,7) of $(x,y)$, $(z,w)$, $(t,u)$, and in (7,7) a swap with $g_\aleph\to-g_\aleph$ (M1). Up to them the cone $xz\pm yw=0$ is $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, the (7,7) null members and $(x,y)=0$ (M2).
- For $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$ the γΓ-channel forces $u\equiv0$; one scalar equation remains, nonlinear after reparametrisation for $c\neq0$ except $c^2=\tfrac{11}6$ in (14,0), which is then settled on every open set (M3, M4, M6(c)). At a point of $s(X)$ the fibre equation forces, for every $\mu$, $\mu'=0$ and the order-1 condition (C1)$_0$, which no Ricci-flat curvature satisfies (Gröbner certificate; M5). In Lorentzian bases non-Ricci-flat curvature does satisfy it (REVIEW7-GU R9-5), and M6 excludes that case by a separate step.
- **Results.** No solution near $s(X)$ for any member with $xz\pm yw=0$, $(x,y)\neq0$ (M6(a), M7(a)); none on any open set for $\mathrm{span}(W_A,W'_A)$, including $W'_A$ (M6(b)), and for $x=\pm y$ in (7,7) (M7; invariant line $\sigma=-\kappa/208$ with $\nabla^2r=-\tfrac5{12}r\,\mathrm{Ric}$). For $(x,y)=0$, $\Upsilon=0$ iff $T=0$, for every metric (M8). With round 8: no solution near $s(X)$ for any member with $(x,y)\neq0$, and none on any open set except possibly the members with $xz\pm yw=0$, $x^2\pm y^2\neq0$, $(z,w)\neq0$ (up to M1 $x(W_A+cW'_C\mathrm{Ad}_\Gamma)$, $c\neq0$, $c^2\neq\tfrac{11}6$ in (14,0)), away from $s(X)$ (M10).
- **Open.** That exception (off the section the first integrability condition fails at every sampled point, so a solution would need special base curvature, M9, computed); Z1; other sectors.
- No statement made before round 9 is contradicted.

## Round 10 — $\Upsilon=0$ on a neighbourhood in Z1 (GU-054, GU-G14)

Write-up: `formalization/geometric-unity/N-z1-neighbourhood.md`. Code: `computation/geometric-unity/N_z1_neighbourhood.py`, output committed, exact, 44 asserted checks, ~7.5 min (450 s) on one core; a second run was identical apart from the timing lines. Record: `metrics/geometric-unity/N-z1-neighbourhood.yaml`. Not yet internally reviewed.

- Z1 geometry: the Euler field is a gradient with $\nabla E=\tfrac12P_H$, not Killing, and a Ricci eigenvector ($-1/c$). The flat model is $\mathrm{Ric}=\beta_Sg+4bc\,e^{-s/2}\nabla^2e^{s/2}$, Einstein iff $\lambda=\tfrac15$.
- Neighbourhood solutions reduce to Einstein pieces or the line equation $\nabla^2r=-\tfrac7{18}r\,\mathrm{Ric}+\nu g$. Its 3-jet conditions along the radial set force a Ricci-flat base and $\lambda\in\{\tfrac15,\tfrac{13}{56}\}$; the trace equation removes $\tfrac{13}{56}$; at $\tfrac15$ an order-$\varepsilon^2$ pairing gives $\tfrac{45}{8c_H^2c_V}|W|^2=0$.
- **Result.** For $\circledast_\sharp$, $W_A$ and L8's classes (A), (B): none in (7,7). In (14,0), only $\lambda=\tfrac15$ over a Ricci-flat base; for analytic fields exactly the flat-base, $T_e\equiv0$ solutions (F7 + I7(b), $c_V\kappa^2\ge9360$), which hold on all of $Y$. F9 for Z1 is proven over flat bases and near $s(X)$ over definite bases.
- One input ($\nabla_U\mathrm{Ric}$ at $s(X)$, N6(f)) is computed, with a Schur reduction. Open: limit points without analyticity, indefinite bases in F9, the members left open by L.3.
- F's sentence after F6 ("(E) holds only along $s(X)$") was too broad for flat ℵ, and I7(b)'s "(E) on a neighbourhood is not established in this sector" did not cover the flat base; both are qualified (`reports/CORRECTIONS.md`). Otherwise no earlier statement is contradicted.

## Round 11 — the ℷ-variation (GU-G15, GU-G14, GU-054)

Write-up: `formalization/geometric-unity/O-gimel-variation.md`. Code: `computation/geometric-unity/O_gimel.py`, output committed, exact, 59 asserted checks, ~2.5 min (158 s) on one core; a second run was identical apart from the timing lines. Record: `metrics/geometric-unity/O-gimel-variation.yaml`. Not yet internally reviewed.

- How ℷ enters: Z1 only through ℵ, with a ℷ-independent volume form (the canonical natural measure $|\det y|^{-n/2}d^nx\,d^Ny$ up to a constant); Z2 through $ℷ_x$ and ℵ. A compactly supported variation of ℷ is not compactly supported on $Y$, so the ℷ-equation needs a fibre weight; the ϖ-equation does not.
- Weights. With the fibre-invariant measure every nonzero fibre integrand diverges. Over every indefinite base (all (7,7) cases) a finite $SO_0(g_x)$-invariant measure lives on $\mathbb R^\times g_x$ (Poincaré recurrence for boosts on $\mathrm{Sym}^2$). So no natural weight function converges, and the convergent natural weights are measures on the scaling ray through $s(X)$. Over (4,0) an explicit class converges.
- The equation. Empty at $T=0$; at first order, the weighted Einstein–Hilbert equation of $g_\aleph$. On the invariant sector with constant $t,u$ and $\Upsilon=0$ kept unweighted, it is exactly 4D vacuum GR with $\Lambda=\Lambda_A-\tfrac{7c_H}6\tilde P/a$. In Z2 this is weight-independent; in Z1 it is multiplied by a moment ratio of the weight, equivalent to a change of the free fibre scale. As $T\to0$, $\Lambda\to\Lambda_A$: reading (A) of F, not (P). A smooth (4,0) weight adds a positive $|\mathrm{Ric}_0|^2$ term in Z2. This is the ℷ-equation alone; in this sector $\Upsilon=0$ itself holds only for the flat-base solutions of round 10.
- If the weight is used for all fields, the ϖ-equation becomes $w\Upsilon+\tfrac12\circledast(dw\wedge T)$, and for ray weights it forces $a=0$: no Einstein term.
- On-shell: among round 10's flat-base Z1 solutions in (14,0), exactly those with $c_V\kappa^2=12480$ (for $x'=1$) also solve the ℷ-equation. None is known in (7,7).
- G.2.5's list of ℷ-dependent ingredients included "vol", which is ℷ-independent in Z1; corrected there. No other earlier statement is contradicted.

## F4(c) for every $n$ (round 2 F, completed)

Write-up: `formalization/geometric-unity/K-z1-ricci-general-n.md`. Code: `computation/geometric-unity/K_z1_ricci.py`, output committed, exact, ~8 min on one core, no repository imports. Record: `metrics/geometric-unity/K-z1-ricci.yaml`.

- **Result** (refereed in `reports/reviews/REVIEW5-GU-2026-10-30.md`, correct as written). Theorem F4(c), the Ricci blocks of Z1 along $s(X)$, is proven for every $n\ge1$, every signature of $g$, every $\lambda\neq1/n$ and every $c_H,c_V\neq0$. The closed form is the one F states.
- **Proof.**
  - K1: an exact formula for $\mathrm{Ric}(h+\delta)-\mathrm{Ric}(h)$ at a zero of $\delta$, linear in $\nabla^2\delta$ and quadratic in $\nabla\delta$.
  - K2: the jets of $\delta=\mathrm{Z1}-\mathrm{Z2}$ on $s(X)$ for a curved base. They depend on $g_x$ alone, because the curvature of the base could enter only through the A-tensor, which vanishes on $s(X)$.
  - K3: an evaluation of K1 on those jets for symbolic $n$ (Sherman–Morrison, and an $O(g)$-invariant diagram calculus with loop weight $n$).
  - The HV block vanishes by the isotropy element $-1\in O(g_x)$ (equivalently the $h_p$-isometry $P=\mathrm{id}_H\oplus(-\mathrm{id}_V)$, which reverses $\nabla\delta$ and preserves $h$ and $\nabla^2\delta$). This replaces the sketch's parity argument; Φ itself is not an isometry of Z1.
- **Checks.** Independent brute-force Ricci tensors of Z1 and Z2 at a point of $s(X)$ for $n=2,\dots,6$, in two or three signatures each, on random curved base metrics in non-normal coordinates; 94/94 exact checks.
- **Consequences.** Every Z1 statement that carried "the status of F4(c)" is now proven under its own other hypotheses: $\Lambda_P$ (Z1) of Cor. F5, the Z1 half of F6, Cor. I2.1 and I7 in Z1, and J3–J4 in Z1. No earlier statement changes.
