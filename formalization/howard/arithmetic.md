# Howard arithmetic: "1 × 1 = 2" (TH-001 – TH-004, TH-007)

**Layer key.** (a) source assertion · (b) our reconstruction · (c) proven · (d) conjectured · (e) testable.
Lean: `proofs/lean/howard/HowardArith.lean` (Lean 4.35.0; the theorems that predate the 2026-10-30 revision were also checked with 4.24.0; core library only, no `sorry`). Full axiom audit of all 80 theorems in `HowardArith.lean.out`: the generic theorems use no axioms, the ℤ, ℤ/4 and Euler theorems at most `propext` and `Quot.sound`, and the ℚ instances also `Classical.choice`, inherited from core's `Rat` lemmas.
Code: `computation/howard/arith_candidates.py`, `axiom_lattice.py` (outputs `*.out` beside them).

## 0. What the sources assert (layer a)

| Rule | Verbatim | Location | Our formal reading (b) |
|---|---|---|---|
| H1 | "One times one equals two" | RS 2015 (TH-S1); JRE #2152 39:03 (TH-S2) | $1\otimes1=2$ |
| SQ | "because the square root of four is two, so what's the square root of two? Should be one" | RS 2015 | $\sqrt[\otimes]{4}=2$ and $\sqrt[\otimes]{2}=1$, where $\sqrt[\otimes]{y}$ is the positive solution of $x\otimes x=y$ |
| H2 | (implied by SQ) | RS 2015 | $2\otimes2=4$ |
| EFF | "If one times one equals one that means that two is of no value because one times itself has no effect" | RS 2015 | weak: $1\otimes1\neq1$. Strong: no element is a $\otimes$-identity |
| INC | "how is it that multiplication, if it means to make more and increase in number, how is one times one equaling one part of the multiplication table?" | JRE #2152 40:00 | $a\otimes b>\max(a,b)$ for all positive $a,b$ in the carrier |
| BAL | "that equal sign is supposed to show that there's a balance … What happened to the other one in this equation?" | JRE #2152 40:30 | $a\otimes b\ge a+b$ for positive $a,b$ |

**Source quality.** The JRE #2152 quotes come from a third-party machine transcript (TH-S2) that was not checked against the audio. INC and BAL rest on it alone; H1 is also in RS 2015.

**Interpretive choices.** (I1) The carrier is $\mathbb Z_{>0}$, $\mathbb Z$, $\mathbb Q$ or $\mathbb R$; results are stated per carrier. (I2) "+" is ordinary addition, except in §5, where we relax it. (I3) Numerals denote the usual integers. (I4) LD/RD mean left/right distributivity over +, AS associativity, CM commutativity.

**Caveat on status.** A secondary report says that on JRE #2171 Howard called 1×1=2 "just a metaphor for challenging the status quo" (TH-S3; transcript not obtained). If so, TH-001 is not a mathematical assertion. We evaluate the mathematical readings regardless.

## 1. The literal reading (c)

**Theorem A0.** In Peano arithmetic, $a\cdot 1=a$ (from $a\cdot S0=a\cdot0+a$). So $1\cdot1=1$. *Lean:* `std_peano`, `std_one_mul_one`.
TH-001 is therefore false as a statement about standard multiplication. It can only be read as a proposal for a different operation $\otimes$. H1 already rules out 1 as a left or right $\otimes$-identity (`h1_not_left_identity`, `h1_not_right_identity`).

## 2. Candidate operations (b, c)

These are derived from the rules in §0, not from critics. Every universal property below was decided by z3 (negation UNSAT, which is a proof), and every "N" comes with a counterexample. Source: `arith_candidates.out`.

| $a\otimes b$ | motivated by | H1 | H2 | LD | RD | AS | CM | identity ($\mathbb R$) | INC ($\mathbb R_{>0}$) | INC ($\mathbb Z_{\ge1}$) | BAL ($\mathbb Z_{\ge1}$) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| $ab$ | baseline | N | Y | Y | Y | Y | Y | 1 | N | N | N |
| $a+b$ | BAL, SQ | Y | Y | N | N | Y | Y | 0 | Y | Y | Y |
| $2ab$ | forced by Thm A3 | Y | N | Y | Y | Y | Y | ½ (none in $\mathbb Z$) | N | Y | Y |
| $ab+1$ | EFF ("must have an effect") | Y | N | N | N | N | Y | none | N | Y | Y |
| $a(b+1)$ | "b times more than a" | Y | N | N | Y | N | N | none | N | Y | Y |
| $ab+\min(a,b)$ | min-type patch | Y | N | N | N | N | Y | none | N | Y | Y |
| $2b$ | LD witness | Y | Y | Y | N | N | N | none | N | N | N |
| $ab$ except $1\otimes1:=2$ | edit the table | Y | Y | N | N | N | Y | none | N | N | N |

$\otimes$-square roots (exact, sympy). Here $\sqrt[\otimes]2$ and $\sqrt[\otimes]4$ are, respectively:
- $a+b$: $1$, $2$
- $2ab$: $1$, $\sqrt2$
- $ab+1$: $1$, $\sqrt3$
- $ab+a$: $1$, $(\sqrt{17}-1)/2$
- $2b$: $1$, $2$

Only $a+b$ and $2b$ reproduce both values in SQ. The patched table also gives $\sqrt[\otimes]4=2$, but there $\sqrt[\otimes]2$ is not unique: both $1$ and $\sqrt2$ solve $x\otimes x=2$.

## 3. No-go and classification theorems (c)

**Theorem A1 (two-sided distributivity).** Let $(S,+)$ be any magma with $1\in S$. Put $\bar2=1+1$, $\bar4=\bar2+\bar2$ and $\bar8=\bar4+\bar4$, and let $\otimes$ satisfy LD and RD. Then H1 ($1\otimes1=\bar2$) implies $\bar2\otimes\bar2=\bar8$. If also H2 ($\bar2\otimes\bar2=\bar4$) holds, then $\bar4+\bar4=\bar4$. In a cancellative monoid this gives $\bar4=0$. So H1+H2+LD+RD is impossible in $\mathbb N,\mathbb Z,\mathbb Q,\mathbb R$, or in any ring whose characteristic does not divide 4.
*Proof.* $(1+1)\otimes(1+1)=1\otimes(1+1)+1\otimes(1+1)=4\cdot(1\otimes1)$. ∎
*Lean:* `two_two_eq_eight`, `nogo_generic` (no axioms used), `nogo_int`.
*Sharpness:* $\mathbb Z/4$ with $a\otimes b=2ab$ satisfies H1, H2, LD, RD (`char4_*`). Among $\mathbb Z/n$, $n\le12$, the only nondegenerate model is $n=4$ (`axiom_lattice.py` Part 3).
*Prior remarks:* a search-result summary reports that a 2018 blog post by Lipton & Regan, and a commenter there, discuss the tension $2\otimes2=8$ vs. 4 (TH-S9). We could not read the post (host unreachable), so nothing is attributed to it; the derivation here does not depend on it.

**Theorem A2 (one distributive law plus one instance of AS or CM).** In any magma $(S,+)$ with $1$ (numerals as in A1), each of the sets {LD, AS, H1, H2}, {LD, CM, H1, H2}, {RD, AS, H1, H2}, {RD, CM, H1, H2} forces $\bar8=\bar4$. Only the single instances $(1\otimes1)\otimes1=1\otimes(1\otimes1)$ and $\bar2\otimes1=1\otimes\bar2$ are needed. So each set is inconsistent in $\mathbb N,\mathbb Z,\mathbb Q,\mathbb R$ and in any ring whose characteristic does not divide 4; $\mathbb Z/4$ with $2ab$ (associative and commutative) shows the condition is needed.
*Proof (LD, AS).* LD gives $1\otimes\bar2=\bar2+\bar2=\bar4$. AS with H1 gives $\bar2\otimes1=1\otimes\bar2=\bar4$. LD then gives $\bar2\otimes\bar2=\bar2\otimes1+\bar2\otimes1=\bar8$. The RD cases are the mirror image. ∎ *Lean:* generic `nogo_generic_ld_assoc`, `nogo_generic_ld_comm`, `nogo_generic_rd_assoc`, `nogo_generic_rd_comm` (no axioms); ℤ: `nogo_int_ld_assoc`, `nogo_int_ld_comm`, `nogo_int_rd_assoc`, `nogo_int_rd_comm`; ℚ: `nogo_rat`, `nogo_rat_ld_assoc`, `nogo_rat_ld_comm`, `nogo_rat_rd_assoc`, `nogo_rat_rd_comm`.
**Corollary (near-rings).** No left or right near-ring structure $(\mathbb Z,+,\otimes)$, $(\mathbb Q,+,\otimes)$ or $(\mathbb R,+,\otimes)$ satisfies H1 and H2 (ℤ, ℚ in Lean; ℝ by A2, since $4\ne0$ in ℝ). Near-ring as in Pilz (1983).

**Theorem A3 (classification).** Every bi-additive $\otimes$ on $\mathbb Z$ has the form $a\otimes b=ab\,(1\otimes1)$. With H1, $\otimes$ is forced to be $a\otimes b=2ab=:a\otimes_2 b$.
*Lean:* `right_additive_linear`, `biadditive_classification`, `forced_two_ab`.
- The same proof holds on $\mathbb Q$: additive maps are $\mathbb Q$-linear.
- On $\mathbb R$, the conclusion holds under any of these regularity conditions on $b\mapsto a\otimes b$ and $a\mapsto a\otimes b$: continuity, measurability, or boundedness on an interval. This is the imported Cauchy theorem (Aczél 1966, Ch. 2).
- Without regularity there are non-measurable $\mathbb Q$-bilinear solutions (Hamel 1905), also with H1: e.g. $a\otimes b=2ab+g(a)g(b)$ with $g$ additive, $g(1)=0$, $g\ne0$. They still satisfy $2\otimes2=4(1\otimes1)$, so A1 is unaffected. Uniqueness of $\otimes_2$ on ℝ therefore needs the regularity condition.
- Brute force over all $3^9$ tables on $\mathbb Z/3$ finds only $k\,ab$ (Part 3).

**Theorem A4 (the forced operation).** The map $\varphi(x)=2x$ is an injective ring homomorphism $(\mathbb Z,+,\otimes_2)\to(\mathbb Z,+,\cdot)$ with image $2\mathbb Z$ (`phi_add`, `phi_mul`, `phi_inj`, `phi_image`). Hence $(\mathbb Z,+,\otimes_2)\cong(2\mathbb Z,+,\cdot)$, a commutative ring without identity. Its properties:
- no $\otimes_2$-identity (`op2_no_identity`);
- no nonzero idempotent (`op2_idempotent`);
- INC and BAL hold on $\mathbb Z_{\ge1}$ (`op2_increasing`, `op2_ge_add`).

On $\mathbb R$, $\varphi$ is a ring isomorphism $(\mathbb R,+,\otimes_2)\cong(\mathbb R,+,\cdot)$. The identity is $\tfrac12$, and $\varphi(\mathbb Q)=\mathbb Q$. So every algebraic and arithmetic property of $\mathbb R$ transfers unchanged, including irrationality and the existence of $\sqrt2$.

**Theorem A5.** If $e$ is a two-sided identity of $\otimes$, then $e\otimes e=e$ (`identity_idempotent`). So EFF-strong ("nothing times itself is without effect") is incompatible with any unital reconstruction. The objection Howard raises against 1 is moved to $e$, not removed.

**Theorem A6 (Archimedean no-go for "multiplication makes more").** Let $S\subseteq\mathbb R$ contain a positive element, and let its positive elements be closed under $+$ and under division by positive integers, for example $S=\mathbb Q,\mathbb R,\mathbb Q_{>0},\mathbb R_{>0}$. Then no $\otimes:S\times S\to\mathbb R$ satisfies LD (used only at positive arguments) together with even the weak form $a\otimes b\ge a$ for all positive $a,b\in S$. The same holds with RD and $a\otimes b\ge b$ for all positive $a,b$. No regularity assumption is needed.
*Proof.* Fix positive $a,b$ and set $f(x)=a\otimes x$. LD gives $f(b)=n\,f(b/n)\ge na$ for every $n\ge1$, which contradicts the Archimedean property of ℝ. ∎
- *Hypotheses corrected 2026-10-30.* The earlier statement had neither the positive element nor "for positive $a,b$". Then $S=\{0\}$, or $S=\mathbb R_{\le0}$ with $\otimes\equiv0$, was a counterexample (`reports/CORRECTIONS.md`).
- *Lean:* `archimedean_ld_generic`, `archimedean_rd_generic` (no axioms) prove A6 for any carrier with `add`, `le`, `pos` such that positives are closed under `add`, `add` is monotone, positives have positive $n$-th parts, and `le` is Archimedean on positives. Instances over ℚ: `archimedean_ld_rat`, `archimedean_rd_rat`, and the lattice cores `core_ld_inc_rat`, `core_rd_inc_rat`. ℝ is not in Lean core; for ℝ the four hypotheses are the standard ordered-field and Archimedean properties (paper).
- *Finite certificates exist only with H1 or H2:* H1+LD force $1\otimes\tfrac14=\tfrac12<1$, and H2+LD force $2\otimes\tfrac14=\tfrac12<2$; RD is the mirror image. z3 returns UNSAT on these ground instances (`axiom_lattice.py`, cores H1+LD+INC, H1+RD+INC, H2+LD+INC, H2+RD+INC).
- *Without H1 and H2 there is no finite certificate.* Let $K=\mathbb R((t^{\mathbb Q}))$ (Hahn series; a real closed field in which $t$ is a positive infinitesimal). Let $\phi:\mathbb Q\to\mathbb Q_{<0}$ be an order isomorphism, e.g. $\phi(q)=q-1$ for $q\le0$ and $\phi(q)=-1/(1+q)$ for $q>0$. Let $\Phi(\sum r_qt^q)=\sum r_qt^{\phi(q)}$; this is an order-preserving ℝ-linear bijection from $K$ onto the series supported in $\mathbb Q_{<0}$, a set closed under $+$ and $\cdot$. Then $a\otimes b:=\Phi^{-1}(\Phi(a)\Phi(b))$ is bi-additive, associative and commutative. For positive $a,b$, $\Phi(a)$ and $\Phi(b)$ are infinitely large, so $a\otimes b$ is infinitely larger than $a$ and $b$, and $a\otimes b>a+b>\max(a,b)$. So LD, RD, AS, CM, INC and BAL all hold in $K$. A finite set of ground instances that holds in $K$ is satisfiable over ℝ (Tarski transfer). Hence the 12 lattice subsets of {LD,RD,AS,CM,INC} that contain INC and LD or RD have no z3 certificate, and rest on A6 alone.
- On $\mathbb Z_{\ge1}$ the theorem fails: $\otimes_2$ satisfies LD, RD and INC. "Multiplication always increases" is consistent with distributivity exactly when the carrier has no arbitrarily small elements.

**Theorem A7 (regular "increasing multiplications" are additions).** Let $I\subseteq\mathbb R$ be an interval and $\otimes:I^2\to I$ continuous, associative and strictly increasing in each variable. Then (imported: Aczél 1949/1966 §6.2; Craigen & Páles 1989) there is a continuous strictly monotone $\varphi$ with $a\otimes b=\varphi^{-1}(\varphi(a)+\varphi(b))$.
- If $\varphi$ is increasing, INC is equivalent to $\varphi>0$ on $I$.
- So every such $\otimes$ satisfying INC is isomorphic, as an ordered semigroup, to a subsemigroup of $(\mathbb R_{>0},+)$.
- H1 and H2 add only $\varphi(2)=2\varphi(1)$ and $\varphi(4)=2\varphi(2)$. The case $\varphi=\mathrm{id}$, i.e. $\otimes=+$, satisfies all of H1, H2, SQ, EFF-weak, INC, BAL, AS, CM.
- On $I=(0,\infty)$, A6 shows that none of these operations distributes over $+$. On $I=[1,\infty)$ the operation $\otimes_2$ is regular, increasing and distributive, which is consistent with A6, because $[1,\infty)$ has no arbitrarily small elements.

**Complete axiom lattice over $\mathbb R$** (`axiom_lattice.out`). The seven rules are H1, H2, LD, RD, AS, CM and INC. EFF and BAL are not among them (see below). All $2^7$ subsets are decided: 58 inconsistent, 70 consistent, 0 open.
- *Certification.* z3 certifies nine cores by UNSAT of their own ground instances: the five cores containing H1 and H2 listed below, and the non-minimal H1+LD+INC, H1+RD+INC, H2+LD+INC, H2+RD+INC. They cover 46 inconsistent subsets. The other 12 are the subsets of {LD,RD,AS,CM,INC} containing INC and LD or RD; they rest on A6 alone (Lean generic and ℚ; ℝ on paper). Before 2026-10-30 only six z3 cores were used, which covered 29 subsets, and the remaining 29 were decided by the paper proof.
- Minimal inconsistent cores:
  - {H1,H2,LD,RD}
  - {H1,H2,LD,AS}
  - {H1,H2,LD,CM}
  - {H1,H2,RD,AS}
  - {H1,H2,RD,CM}
  - {LD,INC} (A6; no finite ground certificate)
  - {RD,INC} (A6; no finite ground certificate)
- Maximal consistent sets:
  - {H1,H2,LD} (witness $2b$)
  - {H1,H2,RD} (witness $2a$)
  - {H1,H2,AS,CM,INC} (witness $+$)
  - {H1,LD,RD,AS,CM} (witness $2ab$)
  - {H2,LD,RD,AS,CM} (witness $ab$)
- *BAL* ($a\otimes b\ge a+b$ for positive $a,b$) added as an eighth rule (Part 2b): 256 subsets, 154 inconsistent (118 by z3, 36 by A6 alone), 102 consistent, 0 open. Each verdict equals that of the subset with BAL replaced by INC, because BAL implies INC on positives and $a+b$, the only witness with INC, satisfies BAL.
- *EFF.* The weak form follows from H1. The strong form (no identity) is not in the lattice, and A5 is the only result about it here. Its status depends on regularity. For {H2,LD,RD,AS,CM}, continuity and A3 force $ab$, which has the identity 1. Without regularity, choose a $\mathbb Q$-complement $V$ of $\mathbb Q\cdot1$ in ℝ, let $\pi$ be the projection onto $\mathbb Q$ along $V$, and set $a\otimes b=\pi(a)\pi(b)$. This is bi-additive, associative and commutative, satisfies $2\otimes2=4$, and has no identity, because $e\otimes v=0\ne v$ for $v\in V\setminus\{0\}$.

**Verdict (c).** No operation on $\mathbb N$, $\mathbb Z$, $\mathbb Q$ or $\mathbb R$ satisfies Howard's two explicit data points (H1, H2) and is a ring, semiring or near-ring over ordinary + (where such a structure exists on the carrier). The same holds in any ring whose characteristic does not divide 4; $\mathbb Z/4$ with $2ab$ is a ring model. If H2 is dropped, the unique bi-additive option on $\mathbb Z$ and $\mathbb Q$ (on ℝ, under regularity) is $\otimes_2$. Via $x\mapsto2x$ it is isomorphic to ordinary multiplication on ℚ and ℝ, and to $(2\mathbb Z,+,\cdot)$ on ℤ.

## 4. Consequences for his other statements under $\otimes_2$ (c)

| Statement | Under $\otimes_2$ |
|---|---|
| $\sqrt[\otimes]2=1$ (SQ) | holds |
| $\sqrt[\otimes]4=2$ (SQ) | fails: $\sqrt[\otimes]4=\sqrt2$ |
| EFF (1 is not idle) | holds. On $\mathbb R$, $\tfrac12$ is idle (A5). On $\mathbb Z$ there is no idle element |
| INC, BAL | hold on $\mathbb Z_{\ge1}$; INC fails on $\mathbb R_{>0}$ (A6) |
| "$x^3=2x=x+x$ is unnatural" (TH-005) | The $\otimes_2$ analogue $x\otimes x\otimes x=x+x$ gives $4x^3=2x$, so $x=1/\sqrt2$. This is irrational, and $\varphi$ maps it to $\sqrt2$. The "loop" is invariant under the relabelling, so it is not caused by $1\cdot1=1$ |

## 5. The idea inside the claim (b, c), and what it costs

1. **Pairs.** $(\mathbb Z,+,\otimes_2)\cong(2\mathbb Z,+,\cdot)$. "1⊗1=2" is ordinary multiplication of even numbers counted in pairs: 1 pair × 1 pair = 4 units = 2 pairs. This fits Howard's "equal and opposite" ontology (every unit is a dyad).
   *Cost:* $2\otimes2=8$, so SQ fails. It is not a new arithmetic.
2. **Addition.** $\otimes=+$ satisfies every verbatim rule in §0.
   *Cost:* no distributivity (A6). The "multiplication" is addition (A7).
3. **Exponents and log scale.** Take the semiring $(\mathbb R,\oplus,\otimes)$ with $a\oplus b=\log_2(2^a+2^b)$ and $\otimes=+$ (Maslov dequantization; Litvinov 2007).
   - It is distributive. H1 and H2 hold, and $1\oplus1=2$.
   - $x\mapsto2^x$ is an isomorphism onto $(\mathbb R_{>0},+,\cdot)$. Under it, "1⊗1=2" reads $2\cdot2=4$.
   - *Cost:* $2\oplus2=3$, so numerals no longer count: the $n$-fold $\oplus$-sum of 1 is $\bar n=1+\log_2 n$. Theorem A1 still applies: $\bar2\otimes\bar2=\bar8$, and $\bar8$ is the real number 4, while $\bar4=3$ (`semirings.out`).
4. **Tropical.** $(\mathbb R\cup\{-\infty\},\max,+)$ also gives H1, H2 and distributivity.
   *Cost:* $1\oplus1=1$, i.e. idempotent addition (Maclagan & Sturmfels 2015).

Each of these four reconstructions is isomorphic to a standard structure, and none yields a new number system. This is not a classification of all consistent readings. The only classification results are A3 (bi-additive operations; on ℝ under regularity) and A7 (continuous, strictly monotone, associative operations). Other consistent rule sets have witnesses that are not relabellings of a standard multiplication: $2b$ and $2a$, the witnesses for {H1,H2,LD} and {H1,H2,RD}, are neither associative nor commutative.

## 6. Failure conditions (what would change these verdicts)

- If the primary book (TH-S5, not obtained) defines $2\times2\ne4$, then H2 is removed. The forced distributive reading is then $\otimes_2$, and the lattice should be re-run with the new data.
- If the book defines "+" differently, §5 applies.
- If a stated rule contradicts the $\otimes=+$ reading (e.g. an explicit $2\times3=7$), A7's reconstruction must be revisited.

## References (imported lemmas)

- Aczél, J. (1966). *Lectures on Functional Equations and Their Applications*. Academic Press (Cauchy equation Ch. 2; associativity equation §6.2).
- Craigen, R. & Páles, Zs. (1989). The associativity equation revisited. *Aequationes Math.* 37, 306–312.
- Hamel, G. (1905). *Math. Ann.* 60, 459–462.
- Pilz, G. (1983). *Near-Rings*, 2nd ed. North-Holland.
- Litvinov, G. L. (2007). Maslov dequantization, idempotent and tropical mathematics. *J. Math. Sci.* 140, 426–444.
- Maclagan, D. & Sturmfels, B. (2015). *Introduction to Tropical Geometry*. AMS GSM 161.

Not read, and not relied on: Lipton, R. J. & Regan, K. W. (2018-05-22). "Does 1 Times 1 Equal 2?", *Gödel's Lost Letter and P=NP* (blog). Known only from a search-result summary (TH-S9).
