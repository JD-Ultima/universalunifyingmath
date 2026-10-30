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
  and its adopted value, cos²θ_W = 1/2 (164.5), is 1070σ off. Expected chance hits: 0.019 among the source's 11
  distinct printed values; 0.14 in a declared family of 192 variants. Five discrete choices lie behind it.
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

  Freedom in declared spaces built from the source's own alternatives: expected chance hits 0.007 (μ), 0.07 (τ)
  and 0.0025 (p); one lepton combination in 504 puts both ratios within 1%; one of the 9 printed proton cells lies
  within 1%; at least 4 or 5 discrete choices each. As with α, the order in which cells and targets were fixed is
  not documented. Here, unlike for α, the text itself marks the proton cell as a candidate.
- **m_n/m_p = 1 (−0.14%).** This is the isospin value. The printed splitting, 5.4×10⁻⁸, is too small by a
  factor 2.6×10⁴.

None of these approaches the improbability of the α(0) agreement: F-CT-4 gives 4.4×10⁻⁶ expected chance hits
within 1σ for α(0).

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

**F14. Distribution and look-elsewhere** (c, e). (2871/67)² is the closest value in every space considered. The
expected chance hits at 1.84×10⁻⁵ are:

| space | expected chance hits |
|---|---|
| printed proton cells × printed electron cells × p (120 variants) | 2–3×10⁻⁴ |
| 72-reading grid, cells chosen per body | 2–3×10⁻⁴ |
| round-2 declared space | 2.5×10⁻³ |
| narrow family (319·3^u/67·3^v)^p | 2–5×10⁻⁵ |
| any integer over 67 | 0.053 |

These assume the cell was fixed without the target, which the text does not document. The source's own mass laws
spread m_p/m_e over a relative width of 1.2×10⁻⁴: +1.84×10⁻⁵ at S², +6.2×10⁻⁶ at the 195.1 offices, and −1.0×10⁻⁴
under 176.2's second-order law without G₂ at the proton. CODATA lies inside this band, so the quoted precision is not a
prediction of the source. At that spread the any-integer chance is 0.34. For the muon the expected chance hits are 3–7×10⁻³;
963 is also the integer nearest 67√(m_μ/m_e). For the tau they are 0.04–0.23.

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
- m_p/m_e: selection of the reading. One coincidence remains, 2871 = 9 × 319, with an expected chance of about 2.5×10⁻³
  in the declared space if the cell had been fixed beforehand.
- m_μ/m_e and m_τ/m_e: conditional on a face choice that the compendium's own definitions exclude.
- The only content independent of the reading is m_τ/m_μ = (441/107)² = 16.987, +1.0% (198σ).
- As exact statements all three are excluded: p at ≥3.6×10⁵σ, μ at 3.9×10⁴σ, τ at 181σ.
- As approximate statements they carry no error budget.

**Changes to round 2** are logged in `reports/CORRECTIONS.md`:
- SC-38's precision is convention-dependent.
- There are 15 printed proton cells, not 9.
- The face remark of F11 applies to the leptons only.
- The τ reference is now PDG 2024: +0.917%.
