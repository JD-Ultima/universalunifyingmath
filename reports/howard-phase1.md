# Howard ("Terryology") — Phase 1 report and first Phase 2 results

**Scope.** We tested the mathematical content of Terrence Howard's public claims: arithmetic, numerical statements, and the Lynchpin geometry. Physical claims are registered but not evaluated.

**Files.**
- Registers: `claims/howard.md` (31 claims: TH-001–014 and TH-020–036), `gaps/howard.md` (16 gaps), `patterns/howard.md`, `sources/howard.md`.
- Proofs and write-ups: `formalization/howard/{arithmetic,numerics,geometry}.md`, Lean `proofs/lean/howard/HowardArith.lean`.
- Code: `computation/howard/`. Metrics: `metrics/howard/`.

## 1. Sources (all informal)

**Obtained**
- *Rolling Stone* 2015 quotes, taken verbatim from the CBS News reprint, which is preserved in a public news corpus (TH-S1).
- A complete timestamped machine transcript of JRE #2152, 2024-05-18, 3 h 08 min (TH-S2).
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

**(a) Arithmetic.** These results are proven. The arithmetic theorems on ℤ, plus the Euler bound G1, are machine-checked in Lean 4 core, which used only `propext` and `Quot.sound`.

Howard's stated data are:
- H1: $1\otimes1=2$;
- H2: $2\otimes2=4$, inferred from "the square root of four is two";
- the rules "multiplication makes more" (INC), "one times itself must have an effect" (EFF) and "balance" (BAL).

| Result | Status |
|---|---|
| A0. $1\cdot1=1$ in Peano arithmetic, so TH-001 can only define a new $\otimes$. | Lean |
| A1/A2. **No-go.** H1+H2 is inconsistent with {LD,RD}, {LD,AS}, {LD,CM}, {RD,AS} and {RD,CM} over any char-0 structure. So **no ring, semiring or near-ring over ordinary + realises Howard's data.** The bound is sharp: $\mathbb Z/4$ is a model. | Lean + z3 |
| A3/A4. **Classification.** If H2 is dropped, the only distributive option is $a\otimes b=2ab$, and $(\mathbb Z,+,\otimes)\cong(2\mathbb Z,+,\cdot)$. This is ordinary arithmetic counted in pairs. It satisfies H1, EFF, INC and BAL on positive integers, but gives $2\otimes2=8$. | Lean (ℤ); paper (ℚ, ℝ) |
| A5. Any identity element is idempotent, so the strong form of EFF excludes every unital reconstruction. | Lean |
| A6. **Archimedean no-go.** "Multiplication always makes more" contradicts either distributive law on $\mathbb Q_{>0}$ and $\mathbb R_{>0}$. It is consistent on $\mathbb Z_{\ge1}$. | paper + z3 |
| A7. Regular increasing $\otimes$ are additions in disguise (imported Aczél theorem). $\otimes=+$ satisfies **every** verbatim rule, but is not distributive. | paper |
| Lattice. All 128 subsets of {H1,H2,LD,RD,AS,CM,INC} are decided: 58 inconsistent, 70 consistent, 0 open. | z3 + paper |
| §5. With a non-standard addition (log₂ or tropical semiring), H1, H2 and distributivity all hold. These structures are isomorphic to standard ones, and in them the numerals stop counting. | sympy |

**The idea inside the claim.** Every consistent reading is a relabelling of standard arithmetic: pairs, addition, or exponents (where $2^1\cdot2^1=2^2$). None yields a new number system. Once we fix any reading that keeps ordinary +, the failing parts of Howard's text are determined.

**(b) Numerical claims.** These are exact unless a precision is stated.
- The √2 "loop" is exactly true arithmetic. It is the fixed point $\sqrt c$ of $x^3/c$, and it is repelling: in float64 it escapes after 28 steps.
- "most numbers satisfy $2\sqrt n>n$" is false: the condition holds only for $n\in(0,4)$.
- His speed-of-light figures are off by $1.3\times10^{-4}$ and $1.8\times10^{-3}$ relative. He quotes 186 624 mi/s, which is exactly $432^2$.
- The element "tones" 40.5/81/162/324 Hz are exactly the Pythagorean E over C = 256 Hz. No tested element property doubles.
- 64 is not a Flower-of-Life lattice count.

**(c) Lynchpin.** The sources alone do not define it well. Six pentagons cannot close a polyhedron (Lean). Under the best-supported reconstruction (six walls, four pyramidal chambers, tetrahedral support), it is the four-film junction of Plateau's laws, i.e. Taylor's T-cone: rays at 109.47°, walls meeting at 120°. That object is well posed but classical. Regular pentagons miss the required angle by 1.4712°. An explicit pentagonal instance is given. "Wave conjugation", "super symmetrical system" and "poles" are undefined, so they cannot be formalized.

## 3. Assessment

| Layer | Content |
|---|---|
| (a) asserted | 31 claims. 13 are mathematically formalizable or partially so (TH-001–010, 020, 022, 024). The rest are physical, historical or undefined. |
| (c) proven | A0–A7, G1, T1–T3, N1–N5. |
| (d) conjectured | P-TH-1 (3-limit numerology), P-TH-3 (octave reading of 1×1=2), P-TH-5 (models are Plateau-compliant), L-T as the patent's intent. |
| (e) testable | TH-009 (once "tone" is defined); TH-020 (patent angles); TH-035 (curricula in 2035). |

No result in this phase is new mathematics. The new content is the complete consistency analysis of Howard's rule set, and its machine-checked certificates.

## 4. Proposed formalization order

1. **Obtain TH-S5 (book) and TH-S3 (JRE #2171).** Extract every stated product and re-run `axiom_lattice.py`. This is the cheapest step and decides whether H2 stands.
2. **Port A3 (ℚ, ℝ), A6 and G1 to Lean + Mathlib** once the cache is reachable. The ℤ results are already done.
3. **Obtain the patent full text (TH-S6).** Check L-T against the figure angles; prediction: 109.47° and 120°. If it matches, formalize T1 in Mathlib's Euclidean geometry.
4. **Define "tone" (TH-G06)** and test TH-009 against NIST spectral lines.
5. Only once definitions exist: TH-022 (poles), TH-023 (wave conjugations) and TH-025 (star prediction).
