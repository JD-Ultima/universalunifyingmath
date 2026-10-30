# Howard arithmetic and geometry — Phase 1 report and first Phase 2 results

**Scope.** We tested the mathematical content of Terrence Howard's public claims: arithmetic, numerical statements, and the Lynchpin geometry. Physical claims are registered but not evaluated.

**Files.**
- Registers: `claims/howard.md` (31 claims: TH-001–014 and TH-020–036), `gaps/howard.md` (16 gaps), `patterns/howard.md`, `sources/howard.md`.
- Proofs and write-ups: `formalization/howard/{arithmetic,numerics,geometry}.md`, Lean `proofs/lean/howard/HowardArith.lean`.
- Code: `computation/howard/`. Metrics: `metrics/howard/`.

## 1. Sources and acquisition limits

**Obtained**
- *Rolling Stone* 2015 quotes, taken verbatim from the CBS News reprint, which is preserved in a public news corpus (TH-S1).
- A complete timestamped third-party machine transcript of JRE #2152, 2024-05-18, 3 h 08 min (TH-S2; GitHub `zreese/jre-excerpts`). It was not checked against the audio, so speaker attribution and digits may be wrong. INC and BAL rest on it alone.
- Patent metadata and abstract-level summaries for US 11,117,065 B2 (Lynchpin), D1,002,750, and the building-block family (TH-S6–S8).

**Not obtained, flagged rather than guessed**
- The original *Rolling Stone* text.
- The book *One Times One Equals Two* (tcotlc.com). This is the decisive gap, TH-G01.
- The 36-page treatise sent to Tyson.
- The JRE #2171 transcript. A secondary report says Howard there called 1×1=2 "a metaphor".
- Any Oxford Union 2017 recording.
- The full patent texts.

Most hosts were blocked by the sandbox proxy. Note: the coordinator's "JRE #2126" is #2152 (first appearance) and #2171 (second).

## 2. Main results (Phase 2)

**(a) Arithmetic.** These results are proven. The arithmetic theorems on ℤ, the generic (magma) forms of A1 and A2, A6 in generic form, the ℚ instances of A1, A2 and A6, and the Euler bound G1 are machine-checked in Lean 4 core (4.35.0). A `#print axioms` audit covers all 80 theorems. The generic theorems use no axioms, the ℤ ones at most `propext` and `Quot.sound`, and the ℚ instances also `Classical.choice`, inherited from core's `Rat` library.

Howard's stated data are:
- H1: $1\otimes1=2$;
- H2: $2\otimes2=4$, inferred from "the square root of four is two";
- the rules "multiplication makes more" (INC), "one times itself must have an effect" (EFF) and "balance" (BAL).

| Result | Status |
|---|---|
| A0. $1\cdot1=1$ in Peano arithmetic, so TH-001 can only define a new $\otimes$. | Lean |
| A1/A2. **No-go.** In any additive magma, H1+H2 together with {LD,RD}, {LD,AS}, {LD,CM}, {RD,AS} or {RD,CM} forces $\bar8=\bar4$. So these sets are inconsistent on ℕ, ℤ, ℚ, ℝ and in every ring whose characteristic does not divide 4. Hence **no ring or semiring on these carriers with ordinary +, and no near-ring on ℤ, ℚ or ℝ, realises Howard's data.** The condition is needed: $\mathbb Z/4$ with $2ab$ is a ring model. | Lean (magma, ℤ, ℚ) + z3 |
| A3/A4. **Classification.** If H2 is dropped, the only bi-additive option on ℤ and ℚ is $a\otimes b=2ab$; on ℝ this needs a regularity condition (continuity, say). $x\mapsto2x$ maps $(\mathbb Z,+,\otimes)$ onto $(2\mathbb Z,+,\cdot)$ (not onto $(\mathbb Z,+,\cdot)$). This is ordinary arithmetic counted in pairs. It satisfies H1, EFF, INC and BAL on positive integers, but gives $2\otimes2=8$. | Lean (ℤ); paper (ℚ, ℝ) |
| A5. Any identity element is idempotent, so the strong form of EFF excludes every unital reconstruction. | Lean |
| A6. **Archimedean no-go.** "Multiplication always makes more" contradicts either distributive law on $\mathbb Q$, $\mathbb R$, $\mathbb Q_{>0}$ and $\mathbb R_{>0}$. The carrier must contain a positive element, and the inequality is quantified over positive $a,b$ (hypotheses corrected 2026-10-30). It is consistent on $\mathbb Z_{\ge1}$. | Lean (generic, ℚ); ℝ instance on paper; z3 with H1 or H2 |
| A7. Regular increasing $\otimes$ are additions in disguise (imported Aczél theorem). $\otimes=+$ satisfies **every** verbatim rule, but is not distributive. | paper |
| Lattice. All 128 subsets of the seven rules {H1,H2,LD,RD,AS,CM,INC} are decided: 58 inconsistent, 70 consistent, 0 open. EFF and BAL are not among the seven. Adding BAL changes no verdict; EFF-strong is not decided. | 46 inconsistent by z3 ground certificates, 12 by A6 alone; 70 by witnesses |
| §5. With a non-standard addition (log₂ or tropical semiring), H1, H2 and distributivity all hold. These structures are isomorphic to standard ones, and in them the numerals stop counting. | sympy |

**Mathematical reconstructions examined.** The four readings we examine are relabellings of standard structures: pairs, addition, exponents (where $2^1\cdot2^1=2^2$), and the tropical semiring. None of them yields a new number system. This is not a classification of all consistent readings. The only classifications are A3 (bi-additive) and A7 (continuous, strictly monotone, associative), and each has hypotheses. Other consistent sets have witnesses that are not relabellings, e.g. $2b$ and $2a$ for {H1,H2,LD} and {H1,H2,RD}. For each reading that keeps ordinary +, the analysis identifies the compatible and incompatible sets of the explicitly stated rules.

**(b) Numerical claims.** These are exact unless a precision is stated.
- The √2 "loop" is exactly true arithmetic. It is the fixed point $\sqrt c$ of $x^3/c$, and it is repelling: in float64 it escapes after 28 steps.
- "most numbers satisfy $2\sqrt n>n$" is false: the condition holds only for $n\in(0,4)$.
- His speed-of-light figures are off by $1.3\times10^{-4}$ and $1.8\times10^{-3}$ relative. He quotes 186 624 mi/s, which is exactly $432^2$.
- The element "tones" 40.5/81/162/324 Hz are exactly the Pythagorean E over C = 256 Hz. No tested element property doubles.
- 64 is not a Flower-of-Life lattice count.

**(c) Lynchpin.** The sources alone do not define it well. Six pentagons cannot close a polyhedron (Lean). Under reconstruction L-T (six walls, four pyramidal chambers, tetrahedral support), the resulting mathematical object is the four-film junction of Plateau's laws, i.e. Taylor's T-cone: rays at 109.47°, walls meeting at 120°. That object is well posed but classical. Regular pentagons miss the required angle by 1.4712°. An explicit pentagonal instance is given. "Wave conjugation", "super symmetrical system" and "poles" are undefined, so they cannot be formalized.

## 3. Assessment

| Layer | Content |
|---|---|
| (a) asserted | 31 claims. 13 are mathematically formalizable or partially so (TH-001–010, 020, 022, 024). The rest are physical, historical or undefined. |
| (c) proven | A0–A7, G1, T1–T3, N1–N5. |
| (d) conjectured | P-TH-1 (prime-factor observation in the stated units), P-TH-3 (octave reading of 1×1=2), P-TH-5 (whether measured models satisfy Plateau's laws), and whether primary drawings support reconstruction L-T. |
| (e) testable | TH-009 (once "tone" is defined); TH-020 (patent angles); TH-035 (curricula in 2035). |

This phase makes no claim of mathematical novelty. Its reproducible contribution is a complete consistency analysis of the seven specified rules, with certificates: z3 ground certificates for 46 inconsistent subsets, and Lean proofs of A6 in generic form and over ℚ for the other 12.

## 4. Proposed formalization order

1. **Obtain TH-S5 (book) and TH-S3 (JRE #2171).** Extract every stated product and re-run `axiom_lattice.py`. This would determine whether the primary text supports the inferred rule H2.
2. **Port A3 (ℚ, ℝ) and the ℝ instances of A1, A2, A6 to Lean + Mathlib** once the cache is reachable. Already done in core Lean: the ℤ results, A1/A2 for any magma and over ℚ, A6 in generic form and over ℚ, and G1.
3. **Obtain the patent full text (TH-S6).** Check L-T against the figure angles; prediction: 109.47° and 120°. If it matches, formalize T1 in Mathlib's Euclidean geometry.
4. **Define "tone" (TH-G06)** and test TH-009 against NIST spectral lines.
5. Only once definitions exist: TH-022 (poles), TH-023 (wave conjugations) and TH-025 (star prediction).
