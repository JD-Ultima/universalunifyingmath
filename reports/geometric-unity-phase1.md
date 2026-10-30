# Geometric Unity: Phase 1 report, with Phase 2 results on three claims

**Scope.** E. Weinstein, *Geometric Unity: Author's Working Draft v1.0* (1 Apr 2021). Hypothesis-neutral: critiques are treated as leads, and every verdict below comes from our own derivation or computation. Files: `sources/`, `claims/`, `gaps/`, `patterns/` (`geometric-unity.md` in each), `formalization/geometric-unity/{A,B,C}-*.md`, `computation/geometric-unity/` (5 scripts, outputs committed), and `metrics/geometric-unity/*.yaml`.

## 1. Sources

- **Obtained: the 2021 draft**, via a third-party GitHub mirror. sha256 `3f28d742…86d4`, 69 pp., PDF timestamps 2021-04-01. The original host was unreachable here, so the copy is unverified against the original.
- **Not obtained:** the 2013 Oxford lectures, the 2020 "First Look" video, all transcripts, and Nguyen–Polya (2021). We attribute no content to any of them.

## 2. Phase 1

- **Claim Register.** 72 claims (GU-001…072). Primary types: 40 theorem-claims, 17 definitions, 6 physical predictions, 6 analogies, and 1 each numerological, conjecture and source statement.
- **Status after Phase 2.** 34 verified; 13 conditional, i.e. holding under a stated reading or repair; 13 open or not well-posed; 5 false as stated (GU-018, 042, 044, 046, 053); 2 trivially true; 5 not mathematical. Counts come from `computation/geometric-unity/register_stats.py`.
- **Specification gaps.** 22 gaps (GU-G01…G22). Three are proved inconsistencies:
  - G11: tilt sign;
  - G13: the field equation is not variational;
  - G17: "imposter" generation has two definitions.
  The underspecified core is the Zorro metric (G06), the Shiab operator (G10), the 14→4 reduction (G14), the reality conditions (G07) and symmetry breaking (G18).

The draft's finite-dimensional algebra is largely correct and checkable. Its analytic core, from the action to the field equations, is not specified precisely enough to be well-posed. Under the most natural reconstruction it is inconsistent.

## 3. Phase 2 results

| # | Result | Status |
|---|---|---|
| A | **Chimeric signature.** All $GL(4)$-natural fiber metrics form the DeWitt family $c\,G_\lambda$. The fiber signature is $(7,3)$ for $\lambda<\tfrac14$ and $(6,4)$ for $\lambda>\tfrac14$ with base $(1,3)$ or $(3,1)$, and $(10,0)$ or $(9,1)$ with base $(4,0)$. The total on $V\oplus H^*$ is $(7,7)$ **iff** the base is Lorentzian, $\lambda>\tfrac14$, and a fixed relative block sign holds. A Riemannian base never gives $(7,7)$. Over all $n$, neutral signatures occur only for bases $(1,3)$, $(3,1)$, $(2,2)$ [$(7,7)$] and $(2,3)$, $(3,2)$ [$(10,10)$]. Also $w_2(C)=\pi^*w_2(X)$, so the draft's $\mathrm{Spin}(7,7)$ bundle exists iff $X$ is spin. The draft's fiber in (3.33) is the double cover of the fiber of $\mathrm{Met}(X)$. | **Proven** (closed form, exact computation) |
| B | **Spinors.** $\mathrm{Cl}(7,7)\cong\mathbb R(128)$. Majorana–Weyl spinors exist. All of the draft's Clifford grade lists for $\mathfrak{so}(64,64)$ and $\mathfrak u(64,64)$ are exact. One MW spinor restricted to $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$ is $((2,1)\otimes16)_\mathbb R$: exactly one SM generation with $\nu^c$, Pati–Salam $(4,2,1)\oplus(\bar4,1,2)$. With $\nu\in S^+$, $\zeta\in T^*\otimes S^-$ the left-handed spin-½ content is $3\times16\oplus144$, plus spin-3/2 $(3,2)\otimes\overline{16}$. All three of the draft's SM tables (pp. 52–53), including the 144 table with its charges, are reproduced exactly. At the SM level, however, the 144 contains a **fourth** generation-equivalent plus 128 exotics. In the complex Weyl reading the content doubles. In the Dirac/$U(64,64)$ reading it is vector-like. "Spin(14,ℂ)→U(128)" is ill-typed. | **Proven** (exact weights and an explicit real model) |
| C1 | **Inhomogeneous gauge group.** $\mathcal G=\mathcal H\ltimes\mathcal N$ is a group with the stated actions. $\tau$ is a homomorphism and equals the right stabiliser of $A_0$. **Lemmas 6.5 and 7.2 are false as printed**, with residual $-2h^{-1}d_{A_0}h. They hold for the opposite-sign tilt $\tau'$, the *left* stabiliser of $A_0$. A one-sign repair makes Secs. 6–7 consistent. | **Proven** (exact symbolic) |
| C2 | **Shiab operator** (reconstruction R-S: soldering invariants $\Phi^1,\Phi^2$). The degrees are consistent. It is gauge covariant only when tilted by $\varepsilon$; untilted, only for the centre $U(1)$. On Levi-Civita curvature, the operator with a commutator in the first term and $i$-anticommutators in the second gives $-\varepsilon^{-1}{*}G\,\varepsilon$ for $c=1$, the Einstein tensor of $Y^{14}$ (not $X^4$). The Weyl part is annihilated. | **Proven** for R-S |
| C3 | **Field equation.** $\circledast\circ d_B$ is not formally self-adjoint: there is an exact counterexample for the Einstein-reproducing member for every $c$, and $O(1)$ defects for all variants. So the printed EL equation (9.7) ("Swervature = Displasion") is not the variation of the printed action (9.4), and no rescaling fixes this. The cubic term also lacks the needed symmetry. | **Proven** for R-S; open beyond R-S |

**Empirical.** No claim is presently testable. The representation content predicts exotic spin-½ states (charges $\pm\tfrac43,\pm\tfrac53,\pm2$, colour sextets and octets) and spin-3/2 states, but no masses, couplings or breaking mechanism are given.

## 4. Proposed formalization order

1. **GU-P09: does a self-adjoint Shiab exist?** Decide by finite linear algebra over $\mathrm{Hom}_{\mathrm{Spin}(7,7)}(\Lambda^2\otimes\mathfrak u,\Lambda^{13}\otimes\mathfrak u)$, organised by Clifford grade. This settles whether any reading of (9.4)/(9.7) is consistent. It is the highest-value and lowest-cost item.
2. **Zorro metric (G06), then the 14→4 reduction (G14, GU-P10).** Write $g_\aleph$ explicitly. Compute $\gimel^*G(g_\aleph)$ by O'Neill formulas with the DeWitt fiber, and determine whether 4D Einstein gravity appears, and with what $\lambda$-dependent corrections.
3. **Full EL system** (G15) under repair R-C1′: the $\varepsilon$- and $\gimel$-variations, and the Noether identity for $\tau'(\mathcal H)$ as the candidate for $\Xi=\mathfrak D\Upsilon$.
4. **Fermion operator (G22)** and its self-adjointness. Then the deformation complex and its symbol and ellipticity (G16).
5. **Optional Lean 4 formalization** of Props. C2/C3/C3′. They are pure group theory on a semidirect product with a 1-cocycle and need no analysis.

Items 1–2 are well-posed now. Items 3–4 depend on 1.
