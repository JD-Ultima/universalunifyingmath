# Howard arithmetic: "1 × 1 = 2" (TH-001 – TH-004, TH-007)

**Layer key.** (a) source assertion · (b) our reconstruction · (c) proven · (d) conjectured · (e) testable.
Lean: `proofs/lean/howard/HowardArith.lean` (Lean 4.24.0, core library only, no `sorry`; axiom audit in `HowardArith.lean.out`: at most `propext`, `Quot.sound`).
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
*Prior art:* the tension $2\otimes2=8$ vs. 4 was noted by Lipton & Regan (2018) and a commenter (TH-S9). We re-derive it here, and the derivation does not depend on theirs.

**Theorem A2 (one distributive law plus one instance of AS or CM).** Over $\mathbb Z$, each of the following sets is inconsistent: {LD, AS, H1, H2}, {LD, CM, H1, H2}, {RD, AS, H1, H2}, {RD, CM, H1, H2}. Only the single instances $(1\otimes1)\otimes1=1\otimes(1\otimes1)$ and $2\otimes1=1\otimes2$ are needed.
*Proof (LD, AS).* LD gives $1\otimes2=2(1\otimes1)=4$. AS with H1 gives $2\otimes1=1\otimes2=4$. LD then gives $2\otimes2=2(2\otimes1)=8\ne4$. ∎ *Lean:* `nogo_int_ld_assoc`, `nogo_int_ld_comm`, `nogo_int_rd_assoc`, `nogo_int_rd_comm`.
**Corollary (near-rings).** No left or right near-ring structure $(\mathbb Z,+,\otimes)$ satisfies H1 and H2. The same holds over $\mathbb Q$ and $\mathbb R$, by the identical computation. Near-ring as in Pilz (1983).

**Theorem A3 (classification).** Every bi-additive $\otimes$ on $\mathbb Z$ has the form $a\otimes b=ab\,(1\otimes1)$. With H1, $\otimes$ is forced to be $a\otimes b=2ab=:a\otimes_2 b$.
*Lean:* `right_additive_linear`, `biadditive_classification`, `forced_two_ab`.
- The same proof holds on $\mathbb Q$: additive maps are $\mathbb Q$-linear.
- On $\mathbb R$, the conclusion holds under any of these regularity conditions on $b\mapsto a\otimes b$ and $a\mapsto a\otimes b$: continuity, measurability, or boundedness on an interval. This is the imported Cauchy theorem (Aczél 1966, Ch. 2).
- Without regularity there are non-measurable $\mathbb Q$-bilinear solutions (Hamel 1905). They still satisfy $2\otimes2=4(1\otimes1)$, so A1 is unaffected.
- Brute force over all $3^9$ tables on $\mathbb Z/3$ finds only $k\,ab$ (Part 3).

**Theorem A4 (the forced operation).** The map $\varphi(x)=2x$ is an injective ring homomorphism $(\mathbb Z,+,\otimes_2)\to(\mathbb Z,+,\cdot)$ with image $2\mathbb Z$ (`phi_add`, `phi_mul`, `phi_inj`, `phi_image`). Hence $(\mathbb Z,+,\otimes_2)\cong(2\mathbb Z,+,\cdot)$, a commutative ring without identity. Its properties:
- no $\otimes_2$-identity (`op2_no_identity`);
- no nonzero idempotent (`op2_idempotent`);
- INC and BAL hold on $\mathbb Z_{\ge1}$ (`op2_increasing`, `op2_ge_add`).

On $\mathbb R$, $\varphi$ is a ring isomorphism $(\mathbb R,+,\otimes_2)\cong(\mathbb R,+,\cdot)$. The identity is $\tfrac12$, and $\varphi(\mathbb Q)=\mathbb Q$. So every algebraic and arithmetic property of $\mathbb R$ transfers unchanged, including irrationality and the existence of $\sqrt2$.

**Theorem A5.** If $e$ is a two-sided identity of $\otimes$, then $e\otimes e=e$ (`identity_idempotent`). So EFF-strong ("nothing times itself is without effect") is incompatible with any unital reconstruction. The objection Howard raises against 1 is moved to $e$, not removed.

**Theorem A6 (Archimedean no-go for "multiplication makes more").** Let the carrier $S$ be closed under $+$ and under division by positive integers, for example $\mathbb Q_{>0}$ or $\mathbb R_{>0}$. Then no $\otimes:S\times S\to S$ satisfies LD together with even the weak form $a\otimes b\ge a$. The same holds with RD and $a\otimes b\ge b$. No regularity assumption is needed.
*Proof.* Fix $a$ and set $f(b)=a\otimes b$. LD gives $f(b)=n\,f(b/n)\ge na$ for every $n\in\mathbb N$, which contradicts the Archimedean property. ∎
- With H1, a finite certificate exists: $1\otimes\tfrac14=\tfrac12<1$. z3 returns UNSAT on the ground instances (`axiom_lattice.py`, core H1+LD+INC).
- On $\mathbb Z_{\ge1}$ the theorem fails: $\otimes_2$ satisfies LD, RD and INC. "Multiplication always increases" is consistent with distributivity exactly when the carrier has no arbitrarily small elements.

**Theorem A7 (regular "increasing multiplications" are additions).** Let $I\subseteq\mathbb R$ be an interval and $\otimes:I^2\to I$ continuous, associative and strictly increasing in each variable. Then (imported: Aczél 1949/1966 §6.2; Craigen & Páles 1989) there is a continuous strictly monotone $\varphi$ with $a\otimes b=\varphi^{-1}(\varphi(a)+\varphi(b))$.
- If $\varphi$ is increasing, INC is equivalent to $\varphi>0$ on $I$.
- So every such $\otimes$ satisfying INC is isomorphic, as an ordered semigroup, to a subsemigroup of $(\mathbb R_{>0},+)$.
- H1 and H2 add only $\varphi(2)=2\varphi(1)$ and $\varphi(4)=2\varphi(2)$. The case $\varphi=\mathrm{id}$, i.e. $\otimes=+$, satisfies all of H1, H2, SQ, EFF-weak, INC, BAL, AS, CM.
- On $I=(0,\infty)$, A6 shows that none of these operations distributes over $+$. On $I=[1,\infty)$ the operation $\otimes_2$ is regular, increasing and distributive, which is consistent with A6, because $[1,\infty)$ has no arbitrarily small elements.

**Complete axiom lattice over $\mathbb R$** (`axiom_lattice.out`). All $2^7$ subsets of {H1, H2, LD, RD, AS, CM, INC} are decided: 58 inconsistent, 70 consistent, 0 open.
- Certified minimal inconsistent cores:
  - {H1,H2,LD,RD}
  - {H1,H2,LD,AS}
  - {H1,H2,LD,CM}
  - {H1,H2,RD,AS}
  - {H1,H2,RD,CM}
  - {LD,INC} (A6)
  - {RD,INC} (A6)
- Maximal consistent sets:
  - {H1,H2,LD} (witness $2b$)
  - {H1,H2,RD} (witness $2a$)
  - {H1,H2,AS,CM,INC} (witness $+$)
  - {H1,LD,RD,AS,CM} (witness $2ab$)
  - {H2,LD,RD,AS,CM} (witness $ab$)

**Verdict (c).** No operation on $\mathbb Z$, $\mathbb Q$ or $\mathbb R$ satisfies Howard's two explicit data points (H1, H2) and is a ring, semiring or near-ring over ordinary +. If H2 is dropped, the unique distributive option is $\otimes_2$, which is ordinary multiplication up to the relabelling $x\mapsto2x$.

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

Each consistent reconstruction is isomorphic to a standard structure. None yields a new number system.

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
