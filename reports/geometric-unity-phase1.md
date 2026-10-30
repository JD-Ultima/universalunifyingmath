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
| F4 | **Z1** along $s(X)$: $\mathrm{II}\neq0$ and $\propto c_H/c_V$; $\mathrm{Ric}_{HH}=\mathrm{Ric}_g+\frac{c_H}{4c_V(n\lambda-1)}g$; $\mathrm{Ric}_{VV}=-\frac{n+1}4G_{1/(n+1)}$; $s_Y=R/c_H+\sigma_1(\lambda)/c_V$. The cosmological term depends on $\lambda$ and also on the route. For $n=4$: projection gives $\Lambda_P=\frac{3c_H}{2c_V}\frac{4-15\lambda}{1-4\lambda}$ and the reduced action gives $\Lambda_A=\frac{5c_H}{4c_V}\frac{5-18\lambda}{1-4\lambda}$. | Proven ($\mathrm{II}$, $s_Y$); Ricci blocks proven for $n\le2$, exact for $n=3,4,5$ |
| F6 | $g_\aleph$ is Einstein along $s(X)$ **iff** Z1, $\lambda=1/(n+1)$ and $\mathrm{Ric}(g)=0$. Then the 14D Einstein equation with $\Lambda_Y=-15/(2c_V)$ is exactly 4D vacuum GR. But $\lambda=\tfrac15<\tfrac14$ gives signature $(8,6)$ or $(10,4)$, not $(7,7)$. With the draft's $(7,7)$, and in Z2 for every $\lambda$, the condition never holds. It also fails off $s(X)$. For flat $\aleph$, Z1 at $\lambda=1/(n+1)$ is locally the symmetric space $SL(n+1)/SO$; for $n=4$ this is $SL(5)/SO(3,2)$ (Prop. F7). | Z2 half and F7 proven; Z1 half has the status of F4's Ricci blocks (exact for $n=4$, also on a Weyl base) |
| F8 | Under R-S/C5, at $T_\omega=0$ the pulled-back $*\Upsilon$ is exactly $-(G+\Lambda_Pg)\otimes\gamma$, a precise version of (9.10). Still $\Upsilon\neq0$, because its vertical–vertical block is non-zero. $T_\omega=t\Phi^1$ adds a 14D cosmological term, so F6 decides solvability. | Conditional; constants open |

**Verdict on GU-G14.** Under R-Z, 4D Einstein gravity with a cosmological term does arise from $Y$, in the projection and reduced-action senses. Its $\Lambda$ is fixed by the fiber geometry: $|\Lambda|\sim1/\ell_F^2$, with no small parameter. Matching the observed value needs $\ell_F\approx2\times10^{26}$ m.

The answer depends on three choices the draft leaves open:
- the horizontal metric (Z1 or Z2);
- $\lambda$ (Z1 only);
- the readout (projection, reduced action, or the field equation on $Y$).

An exact 4D vacuum reduction exists only at $\lambda=1/(n+1)$, which is incompatible with the draft's $(7,7)$ choice. Claims: GU-006 C, GU-021 C, GU-024 V, GU-054 C (previously O).
