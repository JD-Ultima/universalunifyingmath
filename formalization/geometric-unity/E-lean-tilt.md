# E. Lean 4 formalization of Props. C1–C3′ (inhomogeneous gauge group and the sign of the tilt)

File: `proofs/lean/geometric-unity/GUTilt.lean` (Lean 4.35.0, core only, no Mathlib, no imports). Recorded output: `GUTilt.lean.out`. It contains 42 theorems and no `sorry`. `#print axioms` lists only `propext` and `Quot.sound`, or nothing. Check with `lean proofs/lean/geometric-unity/GUTilt.lean`.

## What the abstract model is

- $H$ is a group and $N$ is an abelian group.
- $x\mapsto x^h$ is a right action of $H$ on $N$ by additive automorphisms. In the draft this is $\varpi\mapsto h^{-1}\varpi h$ on $\Omega^1(Y,\mathrm{ad})$.
- $c:H\to N$ is a 1-cocycle, $c(h_1h_2)=c(h_1)^{h_2}+c(h_2)$. In the draft, $c(h)=h^{-1}d_{A_0}h$; this identity is lemma L-C1.
- A connection $A_0+a$ is represented by its offset $a\in N$, with gauge action $a\cdot h=a^h+c(h)$.

The draft's (5.9)–(5.11), (6.4), (6.11), (6.15), (6.18) and (7.3) are transcribed under reading R-C1, as the definitions `gmul`, `actR`, `actL`, `tau`, `pi`, `actN`, `mu` and `T`.

## What is proved

| Lean theorem(s) | Statement | Round-1 item |
|---|---|---|
| `gmul_assoc`, `gmul_gone`, `gone_gmul`, `gmul_ginv`, `ginv_gmul` | $\mathcal G=H\ltimes N$ is a group | C1 (G1) |
| `actR_mul`, `actL_mul` | (5.10) is a right action and (5.11) a left action | C1 (G2, G3) |
| `tau_hom`, `stab_right` | $\tau$ is a homomorphism; $\mathrm{Stab}_{\rm right}(A_0)=\tau(H)$, as an iff | C2 (G4′, G5) |
| `pi_tau`, `pi_mul` | Lemma 6.3: $\pi(g\tau(h))=\pi(g)$; $\pi(g_1g_2)=g_1\cdot\pi(g_2)$ | C2 (G6, G7′) |
| `mu_tau_residual`, `T_tau_residual` | $\mu(g\tau(h))=(\mu(g)_1\cdot h-2c(h),\ \mu(g)_2\cdot h)$ and $T(g\tau(h))=T(g)^h-2c(h)$ | **C3** (G8b, G9a) |
| `T_tau_iff` | printed Lemma 7.2 holds at $h$ **iff** $2c(h)=0$ | C3 (new: sharp form) |
| `tau'_hom`, `stab_left` | $\tau'$ is a homomorphism; $\mathrm{Stab}_{\rm left}(A_0)=\tau'(H)$, as an iff | C3′ (P1, P2) |
| `mu_tau'`, `T_tau'`, `B_equivariant` | Lemmas 6.5 and 7.2 hold for $\tau'$; $B=A_0\cdot\varepsilon$ is equivariant for any tilt | C3′ (P4–P6) |
| `pi'_tau'` | $\pi'(g)=(\varpi-c(\varepsilon))^{\varepsilon^{-1}}$ is $\tau'$-invariant | C3′ (P8) |
| `phi_hom`, `phi_tau`, `phi_invol` | $(\varepsilon,\varpi)\mapsto(\varepsilon,-\varpi)$ is an involutive automorphism with $\varphi\circ\tau=\tau'$ | C3′ (P9, P10) |
| `lemma72_as_printed_false`, `lemma65_as_printed_false`, `tau'_not_right_stabiliser`, `pi_not_tau'_invariant` | in the counter-model $H=N=(\mathbb Z,+)$, trivial action, $c=\mathrm{id}$, the printed Lemmas 6.5 and 7.2, P3 and P7 fail | C3, C3′ (P3, P7) |

## What it does not capture

- **Inputs taken as axioms, not proved here.** The cocycle identity for $h^{-1}d_{A_0}h$ (L-C1), and the fact that conjugation is a right action by automorphisms. Both are one-line calculations, checked symbolically in round 1 (`gauge_group_identities.py`, G1–G9).
- **The meaning of "false".** The counter-model shows the printed lemmas do not follow from the group structure. Failure in the actual gauge group is the separate fact `T_tau_iff` + $c(h)=h^{-1}d_{A_0}h\neq0$: the printed lemmas fail at every $h$ that is not $A_0$-parallel.
- **What is out of scope.** No analysis, smoothness or infinite-dimensionality, and nothing about the action (C4.4), the Shiab operator, or D.
