# Audit report: Chain Theory compendium, release 0.1

Source: A. D. Wissner-Gross, *Chain Theory ... A Compendium of Preliminary Findings*, release 0.1,
`chain-theory` commit `a28db899` (2026-10-24). Scope: the compendium only; the alpha Letter is audited
separately and enters here only as a dependency. Register: `claims/chain-theory-compendium.md`
(CT-C001-CT-C037). Notes: `formalization/chain-theory/compendium/00-04`. Code and recorded outputs:
`computation/chain-theory/compendium/`. Metrics: `metrics/chain-theory/compendium/*.yaml`.

## Findings

**F1. Census reproduced** (c). 100 THEOREM / 109 CONDITIONAL / 198 OPEN boxes; 64 theorems tied to
Lean; Parts VI-XI 17/68/78. README's "Tables 56" is not reproduced (54).

**F2. "Tied to Lean" is a section-level label; for the audited claims the Lean has no physics or
arithmetic content** (c). The `lean` column is inherited by every statement of a section. Of the 64
tied theorems, 11 cite only `opaque : Prop` declarations or open-result records, 37 at best `decide`
proofs over hand-entered tables (137 of 150 data modules have only `String` fields), 16 a restated
AlphaDerivation theorem or kernel module. For 141.1, 143.1, 64.1, 64.2, 65.1 the Lean proves counts
of labels such as `"EQUAL"`, `"ZERO"`, `"YES"`; none of the printed numbers (309, 67/(2^28 3^15), Z6) occurs
in the declarations these chapters name, and the E1 denominators occur nowhere in the Lean tree. (Build not re-run.)

**F3. Grades violate the compendium's own chain rule** (c). 66 of 100 THEOREMs list an OPEN premise
directly (robust to parsing choices); 84 list an OPEN or CONDITIONAL one. Premise lists cite whole
chapters and label CONDITIONAL results "Theorem". Under the index's designator mapping, 124.1 (alpha)
lies on premise cycles (with 16.1 and 58.1; 22-node component). Only 43.4% of premise edges are echoed
in `used_by`.

**F4. Load-bearing unproven premises** (c). Five readings of the axiom clauses in section 2 (chapter 1-1),
2.2 (C) and 2.4, 2.5, 2.7, 2.8 (O; "One act clause discharge", "Record clause ...", "The act alphabet", "The
unknowability clause"), underlie 183 of 209 graded claims; the readability typing 9.1 (C) with 9.7 (C),
9.8 (O) underlies 175; the readout law of chapter 18 (18.1/18.3/18.4 C, 18.5/18.6 O) 152.
124.1 (alpha) has 22 direct users and 70 transitive dependents, including 141.1, 143.1, 146.1, 147.1.

**F5. Electron anomalous moment (141.1, C): no QED coefficient is derived or inserted** (a, b, c, e).
The chapter prints T0 = 1, T1 = 309, T2 = 309^2 + closure (closure not printed), T3 = 0 and states "no
digit of any constant displayed here". T3 = 0 is the general truncation of 146.1, not an anomaly
computation. Under the literal identification T_n ~ C_n (alpha/pi)^n with any normalisation, T3 = 0
contradicts C3 = 1.1812414566 (recomputed from the Laporta-Remiddi closed form); dropping C3 moves a_e
by 9.3e4 sigma. Under the source's own "correspondence by name", the box is an internal equality at one
adopted reading (W1; three other readings part) and makes no contact with a_e. Digits matched: 0.

**F6. Hydrogen (143.1, C): 1/n^2 for n <= 8; the scale is off by ~9 orders** (a, b, c, e). Printed
E1 = 67/(2^28 3^15) or 67/(2^29 3^15) "in the electron unit" (two faces, factor 2, none selected); with
that unit = m_e c^2, E1/(alpha^2 mu/2) = 6.5e-10 / 3.3e-10, and E1 contains no alpha. Under our
reconstruction (twirled symmetric power, two-vertex read) 1/n^2 follows iff the register is
two-dimensional; with dimension 3 or 4 the ladder is not Bohr. Fine structure and Lamb shift: no values
(Lamb sign "positive" agrees); hyperfine OPEN with computed reading 0.

**F7. Gauge group (64.2, 65.1, T): the mathematics is correct and standard; the physics rests on choices**
(a, c, e). Verified exactly: image of U(2) x U(3) on C^2 (x) C^3 = S(U(2) x U(3)) =
(SU(2) x SU(3) x U(1))/Z6; all six genre-matrix (dimension, central scalar) entries; the u(1)
representative = 6Y on the 5. The (3,2) genre is singled out only under read cap d1 d2 <= 6
(9.1, CONDITIONAL) plus d1 != d2; cap 8 admits (4,2) as well. The "dim 12, rank 4" uniqueness cell uses
the Standard Model's numbers as targets. The Z6 global form is consistent with every observed field but
not experimentally established. Representations and hypercharges are not derived (66.x OPEN).

**F8. Other register checks** (c, e). The weak-angle ratio 1/2 (164.5), which the source identifies
with cos^2 theta_W, is 1067 sigma from (M_W/M_Z)^2 = 0.77679(26). The printed running (173.1) is
geometric in an unmapped count metric. The "conservation laws" theorem 74.1 proves the conserved-functional
kernel is trivial. The "Chern number" table is titled "awaited". The source itself states that m_p/m_e
"stands outside the comparison by orders" and that decay-width ratios "part at every pair by many orders".
Signature (1,3) (54.1) follows from pairing the Hermitian unit with anti-Hermitian directions under
tr(XY); with Hermitian directions the same form is Euclidean.

**Summary.** No audited chapter outputs a number that agrees with a measured value beyond integer
structure (the exponent of n, a sign, the su(3) + su(2) + u(1) dimensions once (3,2) is chosen). Where the
source's own stated correspondences give numbers (cos^2 theta_W, hydrogen scale, T3), they disagree with
experiment. The abstract's physics list should be read as correspondences by name, which is how the
chapters themselves describe them.

## Proposed order of further work

1. **Statement-level premises** (cheap, decisive for F3/F4): ask for, or reconstruct, per-statement
   logical premises for 2.2-2.8, 9.1, 18.x, 124.1, 200.x; re-grade under the chain rule; break the
   124.1 cycles.
2. **Formalize G1-G3 in Lean/Mathlib** (`03-gauge-group.md`): the only audited content that is precise
   enough; it fixes exactly which inputs (cap 6, d1 != d2, alphabet) select (3,2).
3. **Pin down the read cap** (9.1 and the Letter's read cap): it carries both the gauge-group selection
   and, per the README, the alpha result; coordinate with the alpha audit.
4. **Request the identification maps** (G1 in `gaps/`): record order -> alpha/pi, electron unit,
   count metric -> energy, tick -> s. Without them CT-C015/16/20/22 are untestable; with them several
   are already refuted.
5. **Re-derive the hydrogen body registers** 2^30 3^18, 2^31 3^17 from the stated construction, and test
   reconstruction R against the source's ladder beyond n = 8.
6. **Lean audit at scale**: extend `lean_classify.py` to resolve restated AlphaDerivation theorems and
   the three kernel modules; only those can certify any computation.
7. Lower priority: Parts XIII-XV theorems (second law, distillation, bands), which appear to hold by
   construction (e.g. fidelity 1 "by the Axiom") and need only a tautology check.

## Round 2 — numeric scorecard

Note: `formalization/chain-theory/compendium/05-numeric-scorecard.md`. Code and outputs:
`computation/chain-theory/compendium/scorecard_extract.py`, `scorecard_data.py`, `scorecard.py`,
`scorecard_freedom.py`, `outputs/scorecard*.{txt,tsv}`. Metrics: `metrics/chain-theory/compendium/scorecard.yaml`.
References: CODATA 2022 (`scipy.constants` 1.18.1) and PDG 2024 (MC table `mass_width_2024.mcd`, with some values
transcribed and flagged in the note).

**F9. Coverage** (c). A mechanical screen of the 468 boxes found 242 that name a physical quantity, 174 of them
with numerals. The curated scorecard covers 53 named quantities and 124 printed numbers ("faces"). All 210
anchors (numbers, identification sentences and verdicts) are verbatim. Each of the 114 screened boxes without a row has a reason code. The source states a
normalisation or unit for 9 quantities: α(0) at zero momentum transfer, and masses, binding and the vacuum value
"in the electron unit". The other 43 are correspondences "by name" or "nearest", and for one more the pairing is
ours. The release defers its comparisons 25 times to "the comparison table", and the release contains no such
table.

**F10. Tally** (c, e). 34 quantities can be compared numerically. Classified at their best face: 2 lie within 3σ,
5 within 1%, 9 within a factor 10, and 18 disagree beyond a factor 10 (or are zero where the measured value is
nonzero). Only 16 of the 34 have a face the source adopts; at that face the count is 1 / 1 / 3 / 11. Of the
qualitative checks, 10 agree and 7 disagree. The 10 agreements are the Tsirelson and CGLMP values, the signature,
the gauge dimensions, the 3 generations, the photon, the n⁻², the Lamb sign, the asymptotic-freedom sign and the
Chern integrality. Most of them rest on round-1 choices or on standard quantum mechanics. The 7 disagreements are
a_e at order 3, the hyperfine interval, the Casimir form, the direction of QED running, the electron decay cell,
the rate ordering and the periodic table. Large failures include:
- the W, Z and H masses and v in the electron unit: factors 10³ to 10⁶;
- G_F: factor 10¹⁰;
- α_s at its adopted face: 6.5e-8 against 0.118;
- CKM |J| = √3/18, the maximum possible value: factor 3×10³;
- trimaximal PMNS: factor 15 at |U_e3|²;
- the heavy-quark masses: in the wrong order;
- the decay-board ratios: off by 7 to 15 orders, as the source itself states;
- G: factor 10³⁵.

**F11. Agreements besides α(0)** (b, c, e). Five of the seven agreements within 1% rest on a choice among faces the
source prints without selecting one.
- **sin²θ_W = (17/36)² = 289/1296.** This face is printed at 164.2 (OPEN). We re-derived it: 17 of the 36 Weyl
  pair letters are nonzero at the (3,2) pair ket. It lies −0.83σ from the on-shell value
  1 − (M_W/M_Z)² = 0.223209(259), but 3.6% from the MS-bar value. The source refuses this identification step,
  and its adopted value, cos²θ_W = 1/2 (164.5), is 1070σ off. Illustrative reference counts: 0.019 under a
  uniform-value model for the 11 distinct printed values; 0.14 from local densities in a declared family
  of 192 variants. These are not calibrated p-values. Five discrete choices lie behind the agreement.
- **Lepton and proton mass ratios.** The "derived cells" of 186.1 (THEOREM) are e = 67, μ = 963, τ = 3969 and
  p = 2871. Combined with mass = (first power)² (185.9), they give:
  - m_μ/m_e = 206.587 (−0.088%);
  - m_τ/m_e = 3509.24 (+0.92%);
  - Koide = 0.667660 (+0.15%);
  - m_p/m_e = (2871/67)² = 1836.18646 (+1.84×10⁻⁵, about 10⁶σ).

  The source displays none of these. It prints the hub face of the lepton ratios instead, (963/201)² = 22.95
  and (3969/201)² = 389.9, both a factor 9 low. It calls 2871 "the candidate cell" (20.1, 175.4) and reports it
  "refuted at the bounded closures" (175.3).

  The agreement needs the electron undressed and the muon and tau dressed. That violates the source's own face
  dissolution theorem 155.1 by the factor 3² = 9.

  Declared-space sensitivity: density-window reference counts 0.007 (μ), 0.07 (τ)
  and 0.0025 (p); one lepton combination in 504 puts both ratios within 1%; one of the 9 printed proton cells lies
  within 1%; at least 4 or 5 discrete choices each. As with α, the order in which cells and targets were fixed is
  not documented; this establishes no construction search history or author intent. The text itself
  marks the proton cell as a candidate.
- **m_n/m_p = 1 (−0.14%).** This is the isospin value. The printed splitting, 5.4×10⁻⁸, is too small by a
  factor 2.6×10⁴.

The numerical agreements and declared-family hit counts can be compared directly. Their density references
use different illustrative models and do not support a ranking of calibrated improbabilities. The α(0)
audit's 4.4×10⁻⁶ is an expected count only in its stated uniform-random-target model.

**F12. Internal inconsistencies exposed by scoring** (c).
- The source uses two mass conventions:
  - mass = S² (185.9 and the Koide arithmetic);
  - office = S² = mass² (the 181.5 and 200.6 bridges; the Higgs "mass displays" 16/67 and 12/67).

  The lepton and proton agreements exist only under the first. "Stands outside the comparison by orders" (176.1)
  is true only under the second.
- 176.1 says the lepton ratios are displayed "nowhere", while 176.2 displays them.
- 211.1 derives the rate ordering electromagnetic > weak > strong, opposite to the ordering its own bridge cites.

**Changes to round 1.** These are logged in `reports/CORRECTIONS.md`. CT-C027 ("No model number" for Koide) is
wrong: 157.1 and 154.1 print four Koide values. CT-C026 and the round-1 summary need the qualification of F11.
The round-1 summary sentence still holds for chapters 141, 143, 64 and 65.

## Round 3 — mass cells

Note: `formalization/chain-theory/compendium/06-mass-cells.md`. Code: `computation/chain-theory/compendium/masscells.py`,
`masscells_data.py`. Outputs: `outputs/masscells.txt`, `outputs/masscells_readings.tsv`. Metrics:
`metrics/chain-theory/compendium/mass-cells.yaml`. The question is whether m_p/m_e ≈ (2871/67)² and
m_μ,τ/m_e ≈ (963/67)², (3969/67)² are derivations or selections. All 102 quoted strings are checked verbatim.

**F13. Trace of the integers** (a, c).
- Re-derived from the construction:
  - 67 = 79 − 12: the electron web is parsed from the 152.1 figure.
  - 963 = 107·9 and 3969 = 147·27: the webs are built from 152.1 and 195.1, and site, link, measure, strike, rank and
    dressing all equal the printed values.
  - 319 = 79 + 2·48 + 3·42 + 6·3: from a combinatorial (2,3) record census. This is fixed once the α weight law is
    fixed; 19 of 4608 squared-factor rules give 319.
- Entered:
  - The tau web is the "fixed library" tau, which no sink of the saturation reaches (189.1).
  - 2871 = 319·9 is introduced as "the candidate cell" and "the proton figure of the comparison table" (175.3). It is
    the integer nearest 67√(m_p/m_e) = 2870.97.
  - 186.1's rule, which "returns" 2871, is tested against it.
- Choices among printed readings:
  - The electron first power: 67 at the path reading of 186.1/195.1, or 201 under the dressing law 198.3.
  - The muon tie: "both passed" (20.1) or "both unpassed" (186.1).
  - The proton typing: capped at 319, or uncapped (747).
  - The proton strike: 0 (186.1) or 18 (2709, 20.1).
  - The proton dressing: 9 or 27.
  - The mass law: S, S², S² + closure, or S²(1 + G₂²).

**F14. Distribution and declared-space sensitivity** (c, e). (2871/67)² is the closest value in every space considered.
The following are density-window reference counts, or uniform-phase lattice coverage, at 1.84×10⁻⁵:

| space | illustrative reference count / coverage |
|---|---|
| printed proton cells × printed electron cells × p (120 variants) | 2–3×10⁻⁴ |
| 72-reading grid, cells chosen per body | 2–3×10⁻⁴ |
| round-2 declared space | 2.5×10⁻³ |
| narrow family (319·3^u/67·3^v)^p | 2–5×10⁻⁵ |
| any integer over 67 | 0.053 |

The density references assume local uniformity in log-value; lattice coverage assumes a uniform target phase
modulo the spacing. The corresponding Poisson probability in the output requires an additional illustrative
homogeneous-process assumption. Fixing a cell before comparison alone does not establish these models, and
the figures are not calibrated p-values. The source's own mass laws
spread m_p/m_e over a relative width of 1.2×10⁻⁴: +1.84×10⁻⁵ at S², +6.2×10⁻⁶ at the 195.1 offices, and −1.0×10⁻⁴
under 176.2's second-order law without G₂ at the proton. CODATA lies inside this band, so the quoted precision is not a
prediction with a stated error budget. At that spread the uniform-phase integer reference coverage is 0.34.
For the muon the density reference counts are 3–7×10⁻³; 963 is also the integer nearest 67√(m_μ/m_e).
For the tau the reference counts are 0.04–0.23.

**F15. Consistency** (c). No face of the mass law gives all three agreements.
- The proton ratio holds at both faces that 155.1 equates (2871/67 = 8613/201).
- The lepton ratios hold only at the "mixed face" (the hub numerator over the dock unit). At the definitional dressing
  of 198.3, 185.6 and 152.1 they are 9× low (22.95 and 389.9, as 176.2 prints).
- In a 72-reading grid that reproduces 6 printed tuples, the four cells meet together at exactly one reading, that of
  186.1. Without 186.1's clause "struck at an uncapped entry alone" they meet at none, as 20.1 states.
- The THEOREM 186.1 and the checkpoint of THEOREM 195.1 (electron dressing 1) contradict the CONDITIONALs 198.3, 155.1
  and 185.12 (dressing 3); 195.1's own mass-register row uses dressing 3.
- The same segment walk gives the top sink the proton's cell 2871 (181.5), so m_t = m_p, against a measured
  m_t/m_p = 184.

**F16. Grade** (b, e).
- m_p/m_e: selection of the reading by the source's stated equality test. The relation 2871 = 9 × 319 has a
  density-window reference count about 2.5×10⁻³ in the declared space; this is not a calibrated probability
  and establishes no history for the original candidate.
- m_μ/m_e and m_τ/m_e: conditional on a face choice that the compendium's own definitions exclude.
- The only content independent of the reading is m_τ/m_μ = (441/107)² = 16.987, +1.0% (198σ).
- As exact statements all three are excluded: p at ≥3.6×10⁵σ, μ at 3.9×10⁴σ, τ at 181σ.
- As approximate statements they carry no error budget.

**Changes to round 2** are logged in `reports/CORRECTIONS.md`:
- SC-38's precision is convention-dependent.
- There are 15 printed proton cells, not 9.
- The face remark of F11 applies to the leptons only.
- The τ reference is now PDG 2024: +0.917%.

## Round 4 — statement-level premises

Note: `formalization/chain-theory/compendium/07-statement-premises.md`. Code: `computation/chain-theory/compendium/premises.py`.
Outputs: `outputs/premises.txt`, `outputs/premises_items.tsv` (1069 premise items, each with its verbatim quote and location).
Metrics: `metrics/chain-theory/compendium/statement-premises.yaml`. The task is proposed further work item 1: re-derive
F3 and F4 from what each derivation says it rests on.

**F17. Where the premises are** (a, c).
- The index premise line is chapter-level by the source's own description: it is "generated mechanically from the
  ancestry of the chapters" (`01-framework.tex:109`).
- Statement-level premises appear only in the opening sentence of each derivation. That sentence cites a chapter, a
  description and a grade word, e.g. "the stabiliser of the pair ket (chapter 5-1, an open result)".
- The 297 derivations contain 0 citations by statement number.
- 2.1–2.9, 9.1–9.8 and 18.1–18.6 each share one derivation, so the source gives no per-box premises for 2.2–2.8, 9.7,
  9.8 or 18.1–18.6. The three procedures conclude 2.7, 9.1 and 18.1.
- The mechanical rule (layer M) gives 1069 items, with candidates filtered by grade word: 237 items have one
  candidate, 750 have several. A separate judgment layer (J) covers 21 derivation groups: 101 items, 30 narrowed by a
  verbatim title or text.
- In derivation text, chapter designators follow the module numbering. Where the two readings differ (171 items),
  description–title overlap favours that reading 129 to 3.

**F18. Chain rule at statement level** (c). The rule is at `front.tex:34`: "A theorem is a statement derived from the
axiom, the definitions and the results proved ahead of it; ... a composed claim carrying the lowest status in its
chain."
- **The 14 requested THEOREMs.** Index: 10 list an OPEN premise and 12 an OPEN or CONDITIONAL one. Derivations (J): 8
  declare an OPEN premise in their own text (16.1, 18.2, 58.1, 8.1, 27.1, 39.1, 69.1, 70.1), 10 violate, 4 are
  undetermined (2.3, 2.6, 74.1, 74.4) and 0 are consistent.
- **Verdicts that move.** 74.1 and 74.4 (index: 14 CONDITIONAL premises) become undetermined. 9.2–9.6 have no OPEN
  premise in their derivation (index: 4); under J they violate through 66.4 (CONDITIONAL).
- **All 100 THEOREMs.** 41 declare an OPEN premise, the same under every reading. With J: 35 more violate definitely,
  19 are undetermined, 3 are consistent (2.10, I.7.1, I.52.1) and 2 have no premise sentence. (Corrected after review
  REVIEW3-CT, R3-CT-8: the first run read 217.1, 217.2 and 223.2 as undetermined, giving 38 / 35 / 22; see
  `reports/CORRECTIONS.md`.)
- **Against the index's 66.** Of the 66 THEOREMs the index gives an OPEN premise, 34 declare one in their derivation;
  7 more declare one that the index does not list.
- **124.1 (CONDITIONAL).** Its derivation names two OPEN premises (an open result of section 2, and 57.2). The
  "lowest status" clause would grade it OPEN; its box treats these premises as conditions "each at its own standing".

**F19. The 124.1 cycles** (c).
- **The index cycles do not survive.** None of the four index edges 124.1↔16.1 and 124.1↔58.1 is cited in a
  derivation. 124.1's derivation cites sections 1, 2, 9, 18, 28, 35, 40, 50, 53 and 57. 16.1 and 58.1 cite "an open
  result" of chapter 5-1: the grade word excludes 124.1 (CONDITIONAL), and the description names 131.1.
- **A different cycle does survive.** 124.1 rests on "the theorem of chapter 18" (one of 18.1–18.4). The shared
  derivation of 18.1–18.6 rests on "the explicit two-point function with its records (chapter 5-1, a theorem)" (one of
  124.1, 134.1, 135.1), and 134.1 and 135.1 each cite 124.1 alone. So 124.1 rests on itself whichever box each
  citation means.
- The cycle passes through 35.1 (state space) and 40.1 (Born rule): the component forced under every resolution is
  {35.1, 40.1, 124.1}. With J, the strict graph has one cycle component, {124.1, 18.2, 18.3}.
- The cycle needs the module-numbering reading of "chapter 5-1" in the derivation sentence at `01-framework.tex:691`,
  and "theorem" read as admitting CONDITIONAL. Under the `designators.tsv` reading, 124.1 is on no forced cycle. The
  source says that its running text cites by the module numbering (`back.tex:1024`). But the premise line printed
  beneath 18.1–18.6 (`01-framework.tex:657`–`684`) cites neither chapter 5-1 nor 124.1. So the cycle sets that
  derivation sentence, under its better-supported but not forced reading, against the boxes' own premise line.
- 16.1 and 58.1 lie on no cycle in the strict graph.

**F20. Load-bearing premises at statement level** (c). Transitive THEOREM+CONDITIONAL dependents, of 209, index against
statement level. The lower bound holds whatever box each chapter citation means; the upper bound is reachability in the
any-of graph.

| premises | index | statement level |
|---|---|---|
| 2.4, 2.5, 2.7, 2.8 (as a set) | 183 | ≥ 153 |
| 9.1, 9.7 | 175 | ≥ 151 (with J) |
| 2.2 | 183 | ≥ 1 (only 13.1 names it) |
| 9.8 | 175 | ≥ 0 |
| 18.1/18.3/18.4 | 152 | ≥ 0 |
| 18.5/18.6 | 152 | ≥ 0 |
| 124.1 | 42 (70 of all grades) | ≥ 152 (313 of all grades) |

For 124.1 the direct users are 25 with J and 14 without, against 22 in the index. The upper bounds are 182–185 for
every row.

**Changes to rounds 1–3.**
- Unchanged:
  - All index numbers, recomputed: 8367 edges, an SCC of 22, 66 / 84 THEOREMs, and 183 / 175 / 152.
  - The conclusion that most THEOREMs violate the chain rule: ≥ 76 at statement level (≥ 73 before R3-CT-8).
- Changed:
  - The 66 is an index count; at statement level the count is 41.
  - The F4 loads of 2.2, 9.8 and 18.x are not supported box by box.
  - 124.1's load is larger than F4 states.
  - The F3/F4 cycle statement holds for the index only. At statement level the cycle runs through chapter 18, not
    through 16.1 or 58.1.
- Corrected in `reports/CORRECTIONS.md`: `00-dependency-structure.md` §2 says that `designators.tsv` resolves 4-5
  uniquely. It does not, and the index two-cycle 124.1↔58.1 depends on reading 4-5 as III-5.

## Round 5 — hydrogen body registers and the ladder beyond n = 8

Note: `formalization/chain-theory/compendium/08-hydrogen-registers.md`. Code: `computation/chain-theory/compendium/hydrogen_registers.py`.
Output: `outputs/hydrogen_registers.txt`. Metrics: `metrics/chain-theory/compendium/hydrogen-registers.yaml`. The task is
proposed further work item 5. All 53 quoted strings are checked verbatim.

**F21. The body registers are re-derived** (a, b, c).
- The rules are stated in other chapters:
  - the per-link whole face (fc)² and cut face (f²−1)(c²−1) (16.1 T, 103.1 C);
  - "the register of a body at a face the product over its links of the per-link face" (103.1);
  - the W word's product register (147.1).

  They reproduce all 16 printed genre faces.
- 2^30·3^18 is the product of (d_a d_b)² over the 10 links of the 152.1 electron web. Given the web (entered, round 3),
  it follows from these rules, provided the struck dock links are counted. Without them the product is 2^26·3^14.
- 2^31·3^17 is the same product with only the port B–w1 at its cut face 24. This is a reading:
  - The 103.1 sentence applied at every link gives 2^27·3^11·5^4.
  - Cutting both docks, or reading the port at its support face, gives 2^32·3^16.
  - The port-only reading also produces the cut values the source prints for the Z body and the gluon word (3 of 3
    numbers). The other readings produce none of them.
  - The source states the port-only rule nowhere.
- The units 67 and 201 come from the web. The text pairs the whole face with 201 and the cut face with 67 but gives no
  reason; the cross pairings give 67/(2^28·3^16) and 67/(2^29·3^14).
- The factor 2 between the faces is the hub dressing 3 times the (2,3) cut/whole ratio 2/3.
- The mass face 36 has no formula:
  - The displayed held chain (6 sites, 5 links) has measure 30. The text also gives it 8 links (measure 48), and the
    glossary calls the bond "the link".
  - 36 equals the measure of the docked carrier as the docking act describes it (6 links), and the whole face of one
    (2,3) link.
- Every input of E1 except the mass face is a cell of the electron body alone.

**F22. Nothing beyond n = 8 to test** (a, b, c).
- A mechanical scan covered every ladder sentence: 15 sentences, of which 6 state the range "eight" or "1, …, 8". None
  names a larger rank. No value beyond 8 is printed, and the per-record ladder has no printed value at any rank.
- Under R the two-vertex ladder is exactly 1/n² at every n: the sl2 commutant is scalar, checked for n ≤ 200. With
  support dimension 3, 4 or 6 it equals 1/n² only at n = 1. R agrees with all 8 printed values.
- Two sentences go beyond the values, and under some readings R parts from them:
  - "a pure power of the rank at none" (the per-record ladder) fails for a single read (1/n, at every n) and for the
    second moment ((n²−1)/12, at n = 2).
  - The 45.1 prose, "the order-k term ... over n squared", holds under R only for terms that read the state twice.

  The source defines neither "record" nor "order-k term" at the held pair, so neither sentence decides R.

**F23. The comparison with experiment lacks the source's unit** (a, c, e).
- E1 is identified by name, as "the Rydberg constant in the electron unit". The source itself leaves the mass unit open:
  - 192.1: "the mass unit stays a question of which cell is meant".
  - The appendix of conditions lists "Mass unit as the rung register (an open question)" among the premises of the
    Rydberg-constant result.
- 176.1 uses "in the electron unit" to mean a ratio to m_e. Under that reading E1 is 6.54×10⁻¹⁰ (whole) and 3.27×10⁻¹⁰
  (cut) of the measured value, 9.18 and 9.49 orders low, as round 1 found.
- We declared a family of 36 conversions built from the source's unit spellings. None lies within a factor 10; the
  nearest is 62× off.
- The ratio E_2/E_1 = 1/4 needs no unit. It differs from the measured B(2S)/B(1S) by 4.5×10⁻⁶ relative:
  - α²/16 = 3.3×10⁻⁶ of this comes from the Dirac energies, and the rest from QED and recoil terms beyond them (not
    evaluated here).
  - This is the size of the fine-structure and Lamb terms that 143.1 states without values.

**Changes to rounds 1–4.**
- No earlier statement is wrong, and nothing is added to `reports/CORRECTIONS.md`.
- `02-hydrogen-spectrum.md` §1 says the body registers are "not derived on the page". That holds for the page. The
  registers follow from rules printed elsewhere, given the web and two readings. Its §4 choice (5) is now those two
  readings.
- F6's "none selected" holds for 143.1. 229.1 later resolves the unit to the dock, at a register (16) that matches
  neither printed face.
- G1 and G2 in `gaps/chain-theory-compendium.md` get dated notes.

## Round 6 — the Lean at scale: restated theorems, kernel modules, build and axioms

Note: `formalization/chain-theory/compendium/09-lean-audit-at-scale.md`. Code, in `computation/chain-theory/compendium/`:
`lean_resolve.py` and `lean_resolve_data.py`; `lean_build.py` with `lean_build_patch.tsv`; `lean_audit.sh`;
`lean_audit/{ResolveChecks,Diagnose,AxiomTrace}.lean`. Outputs: `outputs/lean_resolve.{txt,tsv}`,
`outputs/lean_build_{A,B}.{txt,tsv}`, `outputs/lean_axioms.txt`, `outputs/lean_resolve_checks.txt`,
`outputs/lean_diagnose.txt` and `outputs/lean_axiom_trace.txt`. Metrics:
`metrics/chain-theory/compendium/lean-audit-at-scale.yaml`. The task is proposed further work item 6. All 44 Lean
quotes and all 52 box quotes are checked verbatim.

**F24. What the 16 "restated or kernel" theorems cite** (c).
- They cite 8 declarations: 7 restated `AlphaDerivation` conjunctions, with 26 conjuncts in all, and one kernel check.
- All of them compute. 18 conjuncts are proved by `decide +kernel`, 4 by an anonymous constructor and 4 by short
  tactic proofs. None uses `sorry`, `native_decide` or a declared axiom.
- Part of the computation is genuine linear algebra: the 12-matrix su(2) + su(3) + u(1) algebra, its closure, its
  one-dimensional centre, TU = i·Id_6 on C²⊗C³, and a scalar commutant.
- Where the computation stops:
  - The (2,2) centre is not computed. `Expected.centreDims` is entered and no theorem uses it.
  - "Abelian", "conserved", "massless" and "spin one" are Boolean definitions. "Conserved" reduces to d1 ≠ d2,
    because one of its two clauses is true at every pair (checked here).
  - Five opaque hypotheses (`SoftPhotonTheorem` and four others) are not used. The "on-shell coupling" equals the
    alpha closed form at (2,3) by a case split in its definition.
  - The zero momentum transfer is entered: (0,0,0) in 172 of 172 rows.
  - The 172-record weight check is an identity over {±1}², so it does not depend on the records.
  - The "Ward identity" is a set of telescoping-sum lemmas.
  - Six conjunct parts are entered constants equal to literals.
- Some parts are never restated by any chapter declaration:
  - LemmaA, CollapseG1–G3, Generation, Span and SpanRank. Six of these are among the 8 modules at 4 GB or more.
  - GainLemma, which states 18.2's end-set sums and Born count, OneLinkWord and Correspondence.

**F25. For 10 of the 16 rows, the computation belongs to another box or to none** (c).
- Clause by clause:
  - MATCH: 1 (9.4).
  - PARTIAL: 5 (2.1, 9.2, 9.3, 18.2, 196.2).
  - NONE: 10 (2.3, 2.6, 2.9, 2.10, 9.5, 9.6, 56.1, 74.1, 74.4, 78.3).
- Five of the 8 declarations state content of a CONDITIONAL box in their section: ORKC 2.2, GWDC 18.3, MDFN 74.2,
  WDCL 78.1, WSWC 78.2. Where legible, their four-letter codes are those boxes' initials. They reach THEOREM rows only
  because the `lean` column is per section.
- 56.1 says "the deciding table re-derived in Lean at the forty-eight matchings". Its numbers 322714 and 319332 occur
  nowhere in the Lean tree.
- The three kernel modules compute their rows. All 147 data modules enter theirs as literal lists.
- **Coverage re-stated.** By kind, the tiers are unchanged at 16 / 37 / 11. By content, 6 rows cite a computation that
  states their box, in whole or in part. The other 10 fall back to table checks (9) or to opaque Props (1). That gives
  6 / 46 / 12, with the table tier not checked for content.

**F26. Build and axioms** (c).
- The tree is core-only: no Mathlib, requirements by path only, pinned to Lean 4.32.1.
- **Build A.** Under Lean 4.35.0 the published sources do not build. `AxiomDeclaration/Theorems.lean:42` and `:52`
  fail because `decide` finds no `Decidable` instance: 4.35.0 does not unify `Fin worldThree.qudits` with `Fin 3`.
  277 of the 493 modules depend on that module.
- **Build B.** Five such leaf proofs (the two above and `AxiomDeclarationTwo/Theorems.lean:42`, `:52`, `:56`) were
  changed to `rfl`, with no statement changed and no other library declaration using them. Then 493 of 493 modules
  compile, with 0 errors, 0 warnings and 0 `sorry`. Cost: 5398 s summed, 2494 s wall, peak 7.7 GB.
- 215 modules, including all of `AlphaDerivation` and ORKC, AGUQ and GWDC, do not depend on the patch.
- `#print axioms` on the source's full list of 6662 names:
  - `sorryAx`: 0.
  - No axiom other than `propext`, `Quot.sound` and `Classical.choice`.
- Compared with the source's 4.32.1 report, 6259 axiom sets are equal. The other 403 gain `Classical.choice`. It enters
  through core 4.35.0 order instances (`Int`/`Fin.instLinearOrderPackage` → `Classical.propDecidable`), so the profile
  depends on the toolchain.
- The pinned 4.32.1 could not be fetched (the download was refused in this session), so the source's own record is
  not re-verified on its toolchain.

**Changes to rounds 1–5.**
- No number is wrong; one gloss is corrected (`reports/CORRECTIONS.md`, 2026-10-30: the 04 "hand-entered" gloss and
  "Build not re-run").
- F2's 11/37/16 stands as a classification by kind. By content, its top tier is 6.
- "Build not re-run" (F2, G6, CT-C003) is superseded by F26.
- `04-lean-coverage.md` §1 calls these lists "hand-entered". That fails for 2 of the 137 string-only modules, the kernel
  modules LGIS and ESRL, whose rows are computed.
- G6 in `gaps/chain-theory-compendium.md` gets a dated note.
- Item 6 is done as far as core Lean 4.35.0 allows. Still open: a build on 4.32.1, and the content of the 37 table-check
  rows.

## Round 7 — the table-check and opaque-only rows

Note: `formalization/chain-theory/compendium/10-lean-table-rows.md`. Code, in `computation/chain-theory/compendium/`:
`lean_tables.py` and `lean_tables_data.py`, and `lean_audit/TableChecks.lean`. `lean_classify.py`, `lean_resolve.py`,
`lean_resolve_data.py` and `boxes.py` are imported unchanged. Outputs: `outputs/lean_tables.{txt,tsv}`,
`outputs/lean_tables_clauses.tsv` and `outputs/lean_table_checks.txt`. Metrics:
`metrics/chain-theory/compendium/lean-table-rows.yaml`. The task is the rest of item 6: the 48 rows that round 6 left.
All 248 box quotes and every check anchor are checked verbatim.

**F27. What the 37 table checks and the 94 other names prove** (c).
- The 48 rows cite 131 distinct names: 37 table checks (13 under an opaque hypothesis), 38 opaque Props, 53
  open-result records and 3 aliases.
- **Table checks.** Each restates a data module's `checks`. All 37 modules enter their rows as literal lists, and all
  618 check lemmas are `by decide`.
  - By form: 560 say that a key/value string pair occurs, 13 count rows with string fields equal to literals, and 37
    give the row count.
  - Only 8 test `Int` fields, in three modules: WAMR (28.x), SLBC (217.x) and NCRE (234.x).
- **Hypotheses.** The 13 hypotheses are opaque and never used (`fun _ =>`). No other theorem mentions them.
- **The other 94 names** (opaque Props, records whose `statement` is an opaque constant, and aliases) are not
  theorems, and no library theorem mentions them. **For the 11 opaque-only rows the Lean asserts nothing.** 170.1 and
  229.1 say "the deciding table is re-derived in Lean", and their tied declarations are records.

**F28. Of the 48 rows, 2 state a clause of their box** (c).
- Round 6's rule is used. A check that an entered string records the clause's verdict, count or letter is not a
  counterpart; it gives a label-inclusive bound.
- **PARTIAL: 2.**
  - 217.1, 4 of 8 clauses: SLBC checks `W_target - W_source = W_diff >= 0` at the 228 entered reads.
  - 234.1, 1 of 7: NCRE checks that orders and kinds sum to the records at the 246 entered rows.
- **NONE: 46.**
- **Clauses:** 248 in all; 5 stated, 44 recorded only, 199 with nothing.
- **Label-inclusive:** MATCH 1 (28.2: `"5 of 5"`, `"10 of 10"`, `BUILT`), PARTIAL 22, NONE 25.
- **Sensitivity.** Counting WAMR's `wash_moved_acts == 0` as a counterpart, as round 6 counted the entered zero
  transfer, makes 28.1 PARTIAL. Nothing else changes.
- **Notable.**
  - 223.2 says the asymmetric Onsager pairs "are exactly the pairs at a non-zero production difference at every one of
    the 43". Its table enters and checks `onsg.asym.pairs = "24"`. The key is not defined, so the two may count
    different things.
  - 145.1's table stops at 60 entries and covers one chain at six separations, against "all forty-four chains". 225.1's
    table covers 18 of the 20 networks.
  - The box labels of 181.5 and 200.6 name data modules (QARM, HCOS) that exist; the index does not tie them.
  - 3.1's sentence on the reads Z, X, XZ, XZ² is entered in a table tied to 57.1. 3.1 itself says its "Lean proof
    [is] still to be written".
  - In 217.1's table the checked W count never changes (0 at all 228 reads). Its printed 114, 78, 150 and 14 follow
    from the entered rows but are not checked (`Round7.slbc_*`). 28.1's 228 is likewise the unchecked sum of WAMR's
    entered acts.
  - Printed numbers of 5 or more: 135 (row, number) pairs, 23 in a tied `checks`. 15 occur nowhere in the tree, among
    them 101761, 17689, 270584, 322516 and 46563.

**F29. The 64-row table by content** (c).
- Round 6's 16 rows also cite 9 table checks. None has an `Int` relation, so they add no counterpart and round 6's
  verdicts hold over all tied declarations. Recomputed, they agree with round 6 at 16 of 16.
- **All 64:** MATCH 1, PARTIAL 7, NONE 56. **8 rows** have a tied declaration that states part or all of their box:
  6 by computation (round 6) and 2 by relations over entered tables.
- **NONE: 56.**
  - 10: the computation belongs to another box.
  - 35: table checks. In 21 of these the labels record some clause; in 14 they record nothing of the box.
  - 11: opaque/records only.
- **Label-inclusive:** 2 / 32 / 30.

**Changes to rounds 1–6.**
- The tiers 16 / 37 / 11 by kind stand.
- 09's open table tier is resolved. "37 (+9 = 46)" becomes 2 PARTIAL and 44 NONE, and the opaque tier becomes 12 NONE.
  The README's "coverage by content: 6 / 46 / 12" becomes 8 rows stating a clause and 56 stating none.
- Two earlier statements are wrong in detail (`reports/CORRECTIONS.md`, 2026-10-30, round 7):
  - Round 6's printed-number list for 2.3 (`[5]`, a parser artefact; no verdict depends on it).
  - CT-C003's "37 at best string-label table checks". 7 of the 37 rows cite relation checks over entered `Int` fields.
- G6 gets a dated note. Item 6 is done as far as the tied declarations go. Still open:
  - a build on 4.32.1;
  - Lean statements about the mathematical objects rather than entered tables.
