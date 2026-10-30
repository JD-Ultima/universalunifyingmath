# Numeric scorecard of the compendium (round 2)

Source: compendium 0.1, commit `a28db899`. Scripts and recorded outputs:
`computation/chain-theory/compendium/scorecard_extract.py`, `scorecard_data.py`, `scorecard.py`,
`scorecard_freedom.py` → `outputs/scorecard_extract.txt`, `scorecard_candidates.tsv`, `scorecard.txt`,
`scorecard.tsv`, `scorecard_freedom.txt`. Every number below is copied from those outputs. Metrics:
`metrics/chain-theory/compendium/scorecard.yaml`. The α(0) Letter is audited in
`formalization/chain-theory/alpha/audit.md` and enters here as one row.

## 1. Method

**Screen (mechanical).** The 468 boxed statements are segmented with the index numbering (0 mismatches
against `statements.tsv`). A statement is a candidate when its box, or its section's bridge with
citations, cross-reference lists and the boilerplate "the record supplies …" clauses removed, names a
physical quantity from a fixed 42-entry lexicon. Result: 242 candidates, 174 of them with numerals in the
box or display equations. Of the 174, 60 carry a scorecard row. The other 114 each get a reason code:
68 STRUCT (construction counts or registers, no value for the named quantity), 29 DUP (a face already
in a row), 6 NOVAL, 5 MATH, 2 NOREF, 2 AMBIG, 2 UNMAPPED. Four statements were added by reading: 171.1,
175.1, 175.3, 175.4.

**Rows and faces.** A row is one named physical quantity. Its faces are the numbers the text prints for
that quantity; many quantities get several faces ("both faces printed, none selected"). Each face is
anchored by a verbatim string, and so are the identification sentence and, where there is one, the
source's own verdict. All 210 strings (125 numbers, 53 identification sentences, 32 verdicts) are checked
against the source, with 0 missing. The sentences are listed in §7.

**Codes.**
- Identification:
  - E: the number is stated as the quantity, with a normalisation or unit (`under the correspondence
    the fine-structure constant at zero momentum transfer`; `in the electron unit`).
  - N: the source names the quantity only ("the nearest object"; "The correspondence is by name").
  - P: the source names only a generic class, and the pairing is ours.
- Derivation:
  - R: re-derived here exactly.
  - A: we checked the arithmetic from the printed inputs; the inputs themselves cannot be derived from
    the page.
  - T: the value is target-like. Either the source calls the cell a candidate, or a rule is fitted so
    that it returns the value.
  - N: the value is stated without derivation.

**Comparison.** Each face is converted to the measured quantity under its literal reading. A
dimensionless face is taken as printed. A mass in "the electron unit" is taken as the ratio to m_e. The
exponent is a documented choice: mass = (first power)² follows 185.9 and the Koide arithmetic of
154.1/157.1, and mass = first power follows the 181.5 and 200.6 bridges ("office" = squared mass). Bins:
≤3σ; ≤1% (beyond 3σ); ≤×10 (beyond 1%); disagree (beyond ×10, or zero against nonzero). For integers,
signs and forms the bin is agree or disagree. A face is "not comparable" when its unit is never mapped,
for example the tick. Each row is classified at its best face. The adopted face, where the source has
one, is reported beside it.

**References.**
- CODATA 2022 values are read at run time from `scipy.constants` 1.18.1 and asserted against the
  transcription.
- PDG 2024 masses and widths come from PDG's Monte Carlo table `mass_width_2024.mcd` (RPP 2024; shipped
  in scikit-hep `particle` 0.25.2). The run checks the transcription against that file. M_W uses the
  summary value 80.3692(133) GeV.
- The following values are transcribed and not machine-checked: α_s(M_Z) = 0.1180(9),
  α⁻¹(M_Z) = 127.951(9), |V_ud| = 0.97367(32), J = 3.08(14)×10⁻⁵, sin²θ₁₃ = 0.0220(7), N_ν = 2.9963(74),
  and the Planck 2018 parameters as quoted by PDG 2024.
- The weak angle is compared with the on-shell value 1 − (M_W/M_Z)² = 0.223209(259). Its scheme
  dependence is treated in §4.

## 2. Summary (task 3)

| | rows | faces |
|---|---|---|
| Named physical quantities / printed numbers | **53** | **124** |
| · source states a normalisation or unit (E) | 9 | 28 |
| · name or nearest object only (N), or pairing ours (P) | 44 | 96 |
| Numerical comparison possible | 34 | 101 |
| · within 3σ | **2** | 2 |
| · within 1% (beyond 3σ) | **5** | 5 |
| · within a factor 10 (beyond 1%) | 9 | 39 |
| · beyond a factor 10, or zero against nonzero | 18 | 55 |
| Theorem of standard QM reproduced (CHSH 2√2, CGLMP) | 1 | 2 |
| Integer, sign or form: agrees | 9 | 9 |
| Integer, sign, form or zero: disagrees | 7 | 8 |
| Not comparable (tick or count metric unmapped) | 2 | 4 |

Five of the 7 numerical agreements within 1% need a choice among faces the source prints without
selecting one: SC-18, SC-35, SC-36, SC-37 and SC-38. The two that need no such choice are α(0) and
m_n/m_p = 1, which is the isospin value.

Only 16 numerical rows have a face the source adopts. At that face the count is: 1 within 3σ (α(0)),
1 within 1% (m_n/m_p = 1), 3 within ×10, and 11 that disagree.

The release refers 25 times to "the comparison table" (for example "stands at the comparison table
alone"). No boxed statement or table of the release is that table. Part XVI, "The Scorecard", contains
tallies only.

## 3. The scorecard

Columns:
- st/gr: statements and their grades (T THEOREM, C CONDITIONAL, O OPEN, Tb TABLE).
- f: number of faces.
- id: identification code.
- d: derivation code.
- "best face": the face closest to the measurement.
- "verdict": the source's own words, verbatim.
- Values are rounded.

The full per-face table, with n_σ, relative deviation and factor for all 124 faces, is in
`outputs/scorecard.txt`.

| ID | st/gr | quantity | f | id | d | best face (as measured quantity) | measured | result | source's verdict |
|---|---|---|---|---|---|---|---|---|---|
| SC-01 | 36.1 O | CHSH max at a max-entangled pair; CGLMP (qutrits) | 2 | N | R | 2√2; (12+8√3)/9 | 2√2; 4(3+2√3)/9 | theorem reproduced | — |
| SC-02 | 54.1 C | spacetime signature | 1 | N | R | (1,3) | (1,3) | agree (rests on the split choice, round 1) | — |
| SC-03 | 64.2 T, 65.1 T | gauge algebra dims; global Z6 | 1 | N | R | 8+3+1 = 12; Z6 | SM | agree (rests on read cap 6, round 1) | — |
| SC-04 | 83.1 T | number of generations | 1 | N | A | 3 | N_ν = 2.9963(74) | agree (1, 2, 3 hubs built into the roster; 197.5 O: "Face A … is unmet") | — |
| SC-05 | 124.1 C | photon mass, polarisations | 1 | N | A | massless, two transverse states | 0, 2 | agree | — |
| SC-06 | 124.1 C | 1/α(0) | 1 | **E** | R | 137.035999190 | 137.035999177(21) | **≤3σ (+0.60σ)** | — |
| SC-07 | 78.3 T | a_e: order-one coefficient | 1 | N | A | 67 | C1 = 1/2 | disagree (×134) | "compared with nothing here" |
| SC-08 | 141.1 C | a_e: order-three coefficient | 1 | N | A | 0 | C3 = 1.1812 | disagree (round 1: 9e4σ) | "no digit of any constant displayed here" |
| SC-09 | 143.1 C | H level-one binding / m_e c² | 2 | **E** | A | 67/(2²⁸3¹⁵) = 1.74e-14 | 2.661e-5 | disagree (×1.5e9) | — |
| SC-10 | 143.1 C | exponent of n in E_n | 1 | N | R | 2 (n ≤ 8) | 2 | agree (iff 2-dim register, round 1) | — |
| SC-11 | 143.1 C | sign of the Lamb shift | 1 | N | N | + | +1057.8 MHz | agree (sign only) | — |
| SC-12 | 144.1 O | H hyperfine interval | 1 | N | A | 0 | 1420.4 MHz | disagree | "an output of nothing as yet" |
| SC-13 | 145.1 T | Casimir force vs separation | 1 | N | N | "a power at no chain" | ∝ a⁻⁴ | disagree (form) | — |
| SC-14 | 147.1 C, 148.1 C | α at the Z scale | 3 | N | A | α_joined(Z) = 1156.0 | 1/127.95 | disagree (×1.5e5) | "no value of the coupling at any scale is displayed at this chapter" |
| SC-15 | 148.1 C | direction of α running | 1 | N | A | decreasing (×1/319 per rung) | increasing | disagree (under rung↑ = Q²↑) | — |
| SC-16 | 149.1 C, 172.2 O, 173.1 C, 193.1 O | α_s | 7 | N | A | rung-2 ladder 0.0198; 172 dock 0.0174; adopted 149.1: 6.5e-8 | 0.1180(9) | ≤×10 (×6.0) at a free rung; adopted ×1.8e6 | "the value of the coupling stands in the comparison table alone" |
| SC-17 | 173.1 C | direction of α_s running | 1 | N | A | decreasing (×1/36 per rung) | decreasing | agree (same reading as SC-15, so SC-15 and SC-17 cannot both agree) | — |
| SC-18 | 164.2–164.5, 160.1 O, I.38.1 O, I.40.1 T | sin²θ_W | 13 | N | R | (17/36)² = 289/1296 = 0.222994 (164.2 O) | 0.223209(259) on-shell | **≤3σ (−0.83σ)**; adopted 164.5 C: 1/2, +1070σ; MS-bar −3.6% | "the step … to the square of the angle is refused at its first step" |
| SC-19 | I.43.1 O | G_F m_e² | 2 | N | A | 784/4489 = 0.175 | 3.05e-12 | disagree (×6e10) | "equal across the bodies at no reading of the grid" |
| SC-20 | 164.4 (164.2) O | M_W/m_e | 2 | **E** | A | 28/67 = 0.418 | 157279 | disagree (×3.8e5) | "no value of either mass is displayed in this chapter" |
| SC-21 | 200.5 O | M_H/m_e | 2 | **E** | A | 16/67 = 0.239 ("the mass displays") | 245010 | disagree (×1.0e6) | "compared at the comparison table" |
| SC-22 | 200.6 T | M_H/M_W | 1 | N | A | √(Higgs office/carrier office) = 0.5714 | 1.5578 | ≤×10 (×2.7) | — |
| SC-23 | 202.1 T, 86.1 T, 200.5 O | Higgs quartic λ | 5 | N | A | 9/16 (200.5); adopted 202.1: 5.9e-6 | 0.1293 (tree) | ≤×10 (×4.4); adopted ×2.2e4 | — |
| SC-24 | 200.3 C, 200.5 O | v/m_e | 3 | **E** | A | 6561/67 = 97.9 | 481840 | disagree (×4.9e3) | "every mass figure stands at the comparison table" |
| SC-25 | 199.1 O | top Yukawa y_t | 2 | N | A | 319/9 = 35.4 | 0.991 | disagree (×36) | — |
| SC-26 | 170.2 O | Γ_Z/Γ_W | 1 | N | A | 1 (readings beside) | 1.200(23) | ≤×10 (−17%) | "no width and no rate is displayed in this chapter" |
| SC-27 | 170.2 O | Γ_H/Γ_W | 1 | N | A | 5929/4320 = 1.37 | 1.8e-3 | disagree (×770) | same |
| SC-28 | 205.1 Tb | Γ_τ/Γ_μ (decay board) | 1 | N | N | 0.816 | 7.57e6 | disagree (×9.3e6) | "part at every pair by many orders" (reproduced) |
| SC-29 | 205.1 Tb | Γ_H/Γ_μ (decay board) | 1 | N | N | 26.9 | 1.23e16 | disagree (×4.6e14) | same |
| SC-30 | 205.1 Tb | electron decay rate | 1 | N | N | 1108/4489 per tick | 0 (stable) | disagree | — |
| SC-31 | 211.1 C | rate ordering by interaction | 1 | N | A | EM > weak > strong (forced by the registers) | strong > EM > weak | disagree | — |
| SC-32 | 166.4/166.6/168.1 O | CKM \|J\| | 6 | N | R | √3/18 (= max over all 3×3 unitaries) … 0 | 3.08(14)e-5 | disagree (≥×1.4e3) | "made at the comparison table alone" |
| SC-33 | 166.4 O | CKM \|V_ud\|² | 2 | N | R | 1 (identity block); flat 1/3 | 0.9480(6) | ≤×10 (+5.5%) | "no element of the block is displayed here" |
| SC-34 | 166.1 C | PMNS \|U_e3\|² | 1 | N | R | 1/3 (order-3 Fourier) | 0.0220(7) | disagree (×15) | "no angle, phase or figure is displayed at this page" |
| SC-35 | 176.2 O, 186.1 T | m_μ/m_e | 4 | **E** | T | (963/67)² = 206.587 (derived cells); printed (963/201)² = 22.95 | 206.7682827(46) | **≤1% (−0.088%)**; printed face ×9.0 low | "this compendium displays them nowhere" (176.1) |
| SC-36 | 176.2 O, 186.1 T | m_τ/m_e | 3 | **E** | T | (3969/67)² = 3509.24; printed (3969/201)² = 389.9 | 3477.23(23) | **≤1% (+0.92%)**; printed ×8.9 low | same |
| SC-37 | 157.1 O, 154.1 T | Koide Q | 4 | N | A | dock 16684819/24990001 = 0.667660; hub 0.634619 | 0.6666645(51) | **≤1% (+0.15%, +196σ)** | "equals two thirds in the field at no reading" |
| SC-38 | 186.1 T; 175.1/175.3/175.4/180.1/214.6 O | m_p/m_e | 9 | **E** | T | (2871/67)² = 1836.18646 (2871 = 319·9) | 1836.152673426(32) | **≤1% (+1.84e-5, 1.1e6σ)**; 1 of 9 faces | "the candidate is refuted at the bounded closures" (175.3); "stands outside the comparison by orders" (176.1) |
| SC-39 | 186.1 T, O.1 O, 171.1 C | m_n/m_p | 3 | N | A | 1 (adopted: neutron through the deuteron 2871 = proton) | 1.00137842 | ≤1% (−0.14%; isospin value) | — |
| SC-40 | O.1 O | (m_n − m_p)/m_p | 3 | N | A | 5.4e-8 (mass-row closure); adopted 0 | 1.378e-3 | disagree (×2.6e4); sign + agrees | — |
| SC-41 | 181.5 T | m_u/m_d | 1 | N | A | 1 | 0.46(2) | ≤×10 (×2.2) | "no quark mass and no mass ratio is displayed here" |
| SC-42 | 178.1 O, 181.4 C | m_t/m_c | 3 | N | A | 2.24 (181.4, mass = S²); 1/3 (178.1) | 135.6 | disagree (×60); t lightest of the heavy at 178.1 | "compared with experiment at the comparison table alone" |
| SC-43 | 178.1 O, 181.4 C | m_s/m_c | 2 | N | A | 3; 3.51 | 0.073 | disagree (×41; s > c) | same |
| SC-44 | 181.4 C | m_t/M_W | 1 | N | N | 891 | 2.147 | disagree (×415) | same |
| SC-45 | 214.6 O | B_d/(m_e c²) | 2 | **E** | A | 2X = 56/27 (port); 21/2 (dock) | 4.353 | ≤×10 (×2.1); the office arithmetic adds the binding (deuteron heavier than p + n) | "the deuteron office parting at both register faces" |
| SC-46 | 91.1 O, 95.1 O | G m_e²/(ħc) | 2 | N | A | 1/1088391168 = 9.2e-10 | 1.75e-45 | disagree (×5e35) | "no comparison with observation is made here" (113.1) |
| SC-47 | 209.1 O, 115.1 C, 96.1 O | Λ | 3 | N | A | 2/36⁶, 36⁻⁶³, 1/23328 "in tick units" | Λl_P² = 2.85e-122 | not comparable (tick unmapped) | "no value of the term is displayed here" (115.1) |
| SC-48 | 122.1 O | Hubble rate | 1 | N | A | growth 36 per wave | H0 t_P = 1.2e-61 | not comparable | "no value of the scale or of the rate is displayed here" |
| SC-49 | 209.1 O | Bekenstein–Hawking S/A | 3 | N | R | S/hidden = (2ln2+10ln3)/3 = 4.12 (rungs 0, 1); 2.06 (rung 2) | 1/4 | disagree (×16.5); the printed inverse 0.2425 lies 3% from 1/4 | "Open is the area coefficient" |
| SC-50 | I.51.2 T, I.51.3 T | Chern / Hall integrality | 1 | N | N | integer winding (12/24 of a turn per period) | integers | agree (structural) | — |
| SC-51 | 104.1 T, U.1 T | dark fraction of matter | 4 | **P** | A | dock-free share 0.961 (also 0.651, 0.698) | 0.843(2) | ≤×10 (+14%) | — |
| SC-52 | U.1 T, 211.2 C | baryon asymmetry ≠ 0 | 2 | N | N | "a definite excess" (211.2); 0, census symmetric (U.1) | ≠ 0 | agree at one of two faces | — |
| SC-53 | 228.1 O | periods of Z = 1, Z = 10 | 1 | N | N | 2, 11 (hub count) | 1, 2 | disagree | — |

## 4. Agreements besides α(0), and their freedom (`scorecard_freedom.py`)

The method follows F-CT-4:
- (A) the alternatives the source itself prints;
- (B) a declared rule space built only from choices the source prints or names. For this space we
  report the number of variants, the hits at the achieved precision, and the expected number of chance
  hits from the local density of values;
- (C) the list of discrete choices behind the agreement.

For comparison, F-CT-4 for α(0) gives 4.4×10⁻⁶ expected chance hits within 1σ in a declared space of
46 080 rules.

**W. sin²θ_W = (17/36)² = 289/1296 (SC-18; 164.2, OPEN).**
- *Re-derived (R).* Take the (3,2) pair ket on its two-dimensional common support. Of the 36 Weyl
  letter pairs XᵃZᵇ ⊗ XᶜZᵈ, 17 have nonzero amplitude and 14 have no identity letter at either end. This
  reproduces the printed 17/36 and 14/36. With Gell-Mann × Pauli letters the counts are 5 and 3.
- *Deviation by scheme.*

  | scheme | sin²θ_W | deviation of 289/1296 |
  |---|---|---|
  | on-shell, PDG 2024 masses | 0.223209(259) | −0.83σ |
  | on-shell, CODATA 2022 | 0.22305(23) | −0.24σ |
  | MS-bar at M_Z | ≈ 0.23122 | −3.6% |
  | effective leptonic | ≈ 0.23153 | −3.7% |

- *(A) printed faces.* The text prints 13 faces for the angle, 11 of them distinct (as sin²): 1/2 (the
  CONDITIONAL 164.5), 17/36, 7/18, 289/1296, 49/324, 3/8, 2/5, 2/3, 1/3, 0, and 0.612 (I.38.1). Only
  289/1296 lands within 2.2×10⁻⁴ of any scheme. The expected chance count, with 11 values uniform on
  [0, 1] and 4 schemes, is 0.019.
- *(B) declared space.* The space covers the count-face family: six genres × {Weyl, Gell-Mann/Pauli} ×
  {whole, cut} numerator × {all, cut} denominator × {face, Born square} × {sin², cos²}. That gives 192
  variants and 120 distinct values. The expected number of chance hits summed over the four schemes is
  0.14.
- *(C) choices.*
  1. The source takes the Born square of the count face, which is the office-rule route, not the count
     face that the bridge names. The source itself refuses this step: "refused at its first step".
  2. It takes the whole count (17), not the cut count (14).
  3. It uses Weyl letters.
  4. It reads the value as sin², not cos².
  5. The agreement needs the on-shell scheme.

  Against this agreement stand two facts. The source's adopted value for the angle is cos²θ_W = 1/2,
  which is ≈1.07×10³σ away. And its bridge for 164.5 names the office ratio, not the count face, as the
  angle.
- *Verdict.* The agreement is real at the on-shell scheme and re-derivable. It sits at an OPEN box
  whose identification step the source rejects, and chance odds of a few percent cannot be excluded.

**L. Lepton mass ratios from the "derived cells" (SC-35, SC-36, SC-37).**
- *The cells.* 186.1 (THEOREM) "returns the derived cells" 67, 963, 3969 and 2871 for e, μ, τ and p.
  These follow from first power = rank × 3^(unpassed in-legs), with ranks 27 + 40k = 67, 107, 147
  (re-derived). With mass = S² (185.9) they give:

  | ratio | value | deviation |
  |---|---|---|
  | m_μ/m_e | 206.587 | −0.088% |
  | m_τ/m_e | 3509.24 | +0.92% |
  | Koide Q | 0.667660 | +0.15% |

- *What the source prints.* None of these ratios is printed. 176.2 prints the hub face
  (963/201)² = 22.954 and (3969/201)² = 389.9 as "the muon and the tau over the electron". Both are low
  by a factor of 9.0. 176.1 states the ratios are displayed "nowhere".
- *The source's own inconsistency.* The near-agreement needs the electron at its dock-struck cell 67
  ("its in-leg passed", dressing 1), while the muon and the tau keep the dressings 9 and 27. 155.1 (the
  face dissolution theorem) states that the dock and hub readings are one law, with unit-independent
  ratios. The hybrid dock face used in 154.1 and 157.1 violates this by exactly 3² = 9.
- *(B) declared space.* Each measure is struck or not, each hub in-leg is passed or not, the
  stop-clause cells are included, e ∈ {67, 79, 201, 237} and p ∈ {1, 2}. Results:

  | ratio | variants | expected chance hits at the achieved precision |
  |---|---|---|
  | m_μ/m_e | 56 | 0.007 |
  | m_τ/m_e | 72 | 0.07 |
  | both within 1% | 1 of 504 combinations | — |

  The one combination is (67, 963, 3969, p = 2).
- *(C) choices.*
  1. The electron is undressed while μ and τ are dressed.
  2. Mass = S², not office = mass².
  3. The strikes are 12, 10 and 8.
  4. The roster webs themselves.

**P. m_p/m_e = (2871/67)² = 1836.18646 (SC-38).**
- *Deviation.* +1.84×10⁻⁵, which is the most precise agreement after α(0). In σ it is 1.1×10⁶.
- *The cell.* 2871 = 319·9: the (2,3) Born count times dressing 9. It is also the integer nearest
  67√(m_p/m_e) = 2870.974, off by 0.026.
- *What the source says about 2871.* The source calls 2871 "the candidate cell" (20.1, 175.4). It
  reports that no sink reads it: "the candidate is refuted at the bounded closures" (175.3). It says the
  proton in the electron unit "stands outside the comparison by orders at every reading the construction
  has produced" (176.1). That last statement holds only under mass = S, which gives 42.85, 1.6 orders
  low. Under the source's own mass law 185.9 it fails.
- *(A) printed proton cells.* Over the electron unit, 1 of 9 printed proton cells lies within 1%. The
  others include the seed (2727/201)² = 184.1, the first catalogued sink (2187/201)² = 118.4, the capped
  strike (2709/67)² = 1634.8 at −11%, and (25839/201)² = 16526.
- *(B) declared space.* The space has 330 cells and 1320 variants. Records come from the six genre Born
  counts or the printed proton measures, the strike is 0 or 18, dressings are 3ᵃ4ᵇ, e ∈ {67, 201}, and
  p ∈ {1, 2}. The expected number of chance hits at 1.84×10⁻⁵ is 0.0025. If any integer can be built,
  the chance of landing within 0.026 of the target integer is 0.053; if only multiples of 9 can, it is
  0.006.
- *Interpretation.* The text calls 2871 a candidate and calls its use a refutation. 186.1 then fits a
  segment rule that "meets every cell at 8 of its 36 readings". As with α (CT-G10), the order in which
  the cell and the target were fixed is not documented. Unlike α, the source's own words mark the cell
  as target-like.

**Others within 1%.**
- m_n/m_p = 1 (−0.14%) is the isospin value. The mass-row closure gives (m_n − m_p)/m_p = 5.4×10⁻⁸
  against 1.378×10⁻³. The sign is right; the size is too small by a factor of 2.6×10⁴.
- The Koide value at the dock face (+0.15%) follows from the cells of L. It is not independent of them.

**Structural agreements.** These are not independent evidence:
- SC-01 consists of theorems of standard QM.
- SC-02, SC-03 and SC-10 rest on the choices found in round 1 (the split, read cap 6, a 2-dim
  register).
- In SC-04, the generation count is built into the roster.
- In SC-50, the "integer" winding is half a turn per period.
- SC-15 and SC-17 use the same rung towers, so no single orientation of rung against scale makes both
  the QED and the QCD running directions right.

## 5. Internal inconsistencies found while scoring

1. **Two mass conventions.**
   - mass = S² (first power squared): 185.9, 176.2, and the Koide arithmetic of 154.1, 157.1 and 176.2.
   - office = S² = mass²: the 181.5 and 200.6 bridges, and the Higgs "mass displays" 16/67 and 12/67
     (200.5).

   The lepton and proton agreements exist only under the first. Every boson mass in the electron unit
   is 10⁵ to 10⁶ off under either.
2. **Face dissolution against the dock face.** 155.1 proves the ratio carries no reading. 157.1 and
   154.1 nonetheless replace only the electron's cell, which changes every lepton ratio by a factor
   of 9.
3. **"Displays nowhere" against what is displayed.**
   - 176.1 says the lepton ratios are displayed nowhere, yet 176.2 displays (963/201)² and
     (3969/201)².
   - The bridge of 170 says "no width and no rate is displayed in this chapter", yet the box of 170.2
     prints 270/49 and 121/16.
4. **Comparisons deferred to a table that does not exist.** The text defers comparisons to "the
   comparison table" 25 times, and the release has no such table. The only comparison the release
   makes itself is on the decay board (205.1): "part at every pair by six to fifteen orders". It is
   reproduced here at 7.0 orders (τ/μ) and 14.7 orders (H/μ).
5. **Rate ordering against the cited ordering.** The ordering 211.1 derives from the registers
   (electromagnetic > weak > strong) contradicts the ordering its own bridge cites ("the strong decays
   the fastest and the weak the slowest").

## 6. Effect on round-1 conclusions

These are logged in `reports/CORRECTIONS.md`, 2026-10-30.

- **CT-C027.** The round-1 statement "No model number" for Koide is wrong. 157.1 prints K₁ = 16684819/24990001
  and 1857859/2927521, and 154.1 prints two more values.
- **CT-C026.** The round-1 statement "No mass compared successfully" is true of the source's own
  comparisons, but needs a qualification. Under 185.9, the cells of 186.1 imply m_μ/m_e, m_τ/m_e and
  m_p/m_e within 0.09%, 0.92% and 1.8×10⁻⁵.
- **The round-1 summary sentence** ("No audited chapter outputs a number that agrees … beyond integer
  structure") stands for chapters 141, 143, 64 and 65. Compendium-wide, one OPEN face agrees within 1σ
  (SC-18), and the cell agreements of §4 exist.
- **F8's quotation for m_p/m_e** is accurate, but it holds only under the mass = S convention.

## 7. Identification sentences and verdicts (verbatim)

This list is copied from `outputs/scorecard.txt` block [IDENTIFICATION]. LaTeX is kept as printed.
[E/N/P] is the identification code. The statement in parentheses is the unit where the sentence stands.

- **SC-01** [N] (I.20.1): "The pair ket is a maximally entangled state and the bound its Bell value, both tags nearest."
- **SC-02** [N] (54.1): "it equals the Lorentz algebra of the literature by name"
- **SC-03** [N] (64.2): "isomorphic to $(SU(2)$ x $SU(3)$ x $U(1))/\mathbb{Z}_6$ by an explicit covering whose kernel is the cyclic group of order six"
- **SC-04** [N] (83.1): "the three generations are named by correspondence"
- **SC-05** [N] (124.1): "abelian, conserved and massless, with two transverse states at its open read"
- **SC-06** [E] (124.1): "under the correspondence the fine-structure constant at zero momentum transfer"
- **SC-07** [N] (78.3): "the Ward identity, the anomaly, and the dipole factor are named by correspondence" — verdict (78.3): "every number is an exact object of the derivation, compared with nothing here"
- **SC-08** [N] (141.1): "and the order-three term is zero" — verdict (141.1): "no digit of any constant displayed here"
- **SC-09** [E] (143.1): "its level-one binding through the electron body is an exact element of the field in the electron unit"
- **SC-10** [N] (143.1): "E_{n} = \frac{E_{1}}{n^{2}}, \qquad n = 1, \dots, 8"
- **SC-11** [N] (143.1): "the Lamb interval at the hydrogen sink, the s self-term less the p self-term at rung two, is one exact element, positive and split"
- **SC-12** [N] (144.1): "the hyperfine interval is no output of it as yet" — verdict (144.1): "the hyperfine interval an output of nothing as yet"
- **SC-13** [N] (145.1): "so the Casimir effect is a composed cell of the second-law and chain results"
- **SC-14** [N] (148.1): "the coupling through the joined web is an exact nonzero element of the field at all three" — verdict (148.1): "no value of the coupling at any scale is displayed at this chapter"
- **SC-15** [N] (148.1): "and the running is its register growth, the readout at the cap register negative at every rung"
- **SC-16** [N] (172.2): "What this chapter calls the gluon closed read the literature calls the strong coupling, the closure sector its value at the hadron scale" — verdict (172.2): "the value of the coupling stands in the comparison table alone"
- **SC-17** [N] (173.1): "the running, the ratio of successive rungs, is one over thirty-six"
- **SC-18** [N] (164.4): "What this chapter calls the dock share, the count face and the office ratio the literature calls the squared sine of the weak mixing angle, each the nearest object." — verdict (164.5): "The comparison of the cell with the measured angle stands at the comparison table alone."
- **SC-19** [N] (I.43.1): "What this chapter calls the office of the W word the literature calls the Fermi constant" — verdict (I.43.1): "they are equal across the bodies at no reading of the grid"
- **SC-20** [E] (I.43.1): "the first power the mass of the body" — verdict (159.2): "no value of either mass is displayed in this chapter"
- **SC-21** [E] (200.5): "The mass displays at the two faces 16 over 67 and 12 over 67." — verdict (200.5): "the three displays are compared at the comparison table"
- **SC-22** [N] (200.6): "What this chapter calls the Higgs office the literature calls the squared Higgs mass in the unit of the construction, the carrier offices the squared boson masses"
- **SC-23** [N] (202.1): "The Higgs quartic is the ratio of two offices at every Higgs-form sink" — verdict (200.5): "the three displays are compared at the comparison table"
- **SC-24** [E] (200.5): "the vacuum value the register over the dock-struck electron unit 67" — verdict (200.4): "every mass figure stands at the comparison table"
- **SC-25** [N] (199.1): "the top Yukawa cell of the Higgs-and-carrier-offices entry is re-met as the proton office over the register $101761/81$ with the implied per-leg cell $319/9$"
- **SC-26** [N] (170.2): "the W and the Z at 270 over 49 each and the Higgs at 121 over 16" — verdict (170.2): "no width and no rate is displayed in this chapter"
- **SC-27** [N] (170.2): "the W and the Z at 270 over 49 each and the Higgs at 121 over 16" — verdict (170.2): "no width and no rate is displayed in this chapter"
- **SC-28** [N] (205.1): "every cell an exact count per tick with the tick a name" — verdict (205.1): "the ratios of the cells against the ratios of the observed widths part at every pair by many orders"
- **SC-29** [N] (205.1): "every cell an exact count per tick with the tick a name" — verdict (205.1): "the ratios of the cells against the ratios of the observed widths part at every pair by many orders"
- **SC-30** [N] (205.1): "every cell an exact count per tick with the tick a name"
- **SC-31** [N] (211.1): "the (2,3) kind above the hub-leg kind at (2,4) above the (3,3) kind"
- **SC-32** [N] (166.4): "its Jarlskog cell the Jarlskog determinant, the nearest objects" — verdict (166.4): "the comparison with the measured entries is made at the comparison table alone"
- **SC-33** [N] (166.6): "What this chapter calls the block at the read sites the literature calls the quark mixing matrix" — verdict (166.6): "no element of the block is displayed here"
- **SC-34** [N] (166.1): "What this chapter calls the overlap of the corner frame with the depth eigenbasis the literature calls the lepton mixing matrix" — verdict (166.1): "no angle, phase or figure is displayed at this page"
- **SC-35** [E] (176.2): "What this chapter calls the ratio at the second-order read the literature calls the ratio of two lepton masses" — verdict (176.1): "this compendium displays them nowhere"
- **SC-36** [E] (176.2): "What this chapter calls the ratio at the second-order read the literature calls the ratio of two lepton masses" — verdict (176.1): "this compendium displays them nowhere"
- **SC-37** [N] (157.1): "The Koide relation states that the sum of the three charged-lepton masses divided by the square of the sum of their square roots is two thirds" — verdict (157.1): "so it equals two thirds in the field at no reading"
- **SC-38** [E] (176.1): "The proton in the electron unit stands outside the comparison by orders at every reading the construction has produced" — verdict (175.3): "so the candidate is refuted at the bounded closures"
- **SC-39** [N] (O.1): "the mass difference is named by correspondence"
- **SC-40** [N] (O.1): "the mass difference is named by correspondence" — verdict (O.1): "the neutron office parts from the proton office by that closure"
- **SC-41** [N] (181.5): "the ratios the quark mass ratios, each the nearest object" — verdict (181.5): "no quark mass and no mass ratio is displayed here"
- **SC-42** [N] (181.5): "the heavy sink offices the squared masses of the strange, charm, bottom and top" — verdict (181.4): "the quark rows are compared with experiment at the comparison table alone"
- **SC-43** [N] (181.5): "the heavy sink offices the squared masses of the strange, charm, bottom and top" — verdict (181.4): "the quark rows are compared with experiment at the comparison table alone"
- **SC-44** [N] (181.4): "the top-to-W ratio the ratio of the top and W masses" — verdict (181.4): "the quark rows are compared with experiment at the comparison table alone"
- **SC-45** [E] (214.5): "The deuteron and the alpha are exact objects at the register of the common support in the electron unit at both readings of the binding" — verdict (214.5): "the deuteron office parting at both register faces"
- **SC-46** [N] (91.1): "the coefficient that couples the curvature to its source is the gravitational constant" — verdict (113.1): "no comparison with observation is made here"
- **SC-47** [N] (209.1): "the cosmological constant is the vacuum energy density of the field equations" — verdict (115.1): "no value of the term is displayed here"
- **SC-48** [N] (122.1): "the Hubble ratio as the growth thirty-six per wave" — verdict (122.1): "no value of the scale or of the rate is displayed here"
- **SC-49** [N] (209.1): "The entropy of a black hole is proportional to the area of its horizon with a coefficient of one quarter in Planck units" — verdict (209.1): "Open is the area coefficient"
- **SC-50** [N] (I.51.3): "the quantum Hall integer a name beside"
- **SC-51** [P] (104.1): "In cosmology the abundances of the species are fixed at an epoch by the expansion rate and the reaction rates"
- **SC-52** [N] (U.1): "the baryon asymmetry is named by correspondence"
- **SC-53** [N] (228.1): "the objects named the periodic table of Mendeleev"


## Round 3 amendments (from `06-mass-cells.md`; see `reports/CORRECTIONS.md`)

- **SC-38.** The $+1.84\times10^{-5}$ figure depends on the mass-law convention. The source's own printed laws give $+1.84\times10^{-5}$, $+6.2\times10^{-6}$ (195.1 offices) or $-1.0\times10^{-4}$ (176.2). CODATA lies inside this $1.2\times10^{-4}$ band, so the precision is not a prediction of the source.
- **Proton cells.** There are 15 printed proton cells, not 9.
- **Face violation.** The remark that a reading violates face dissolution (155.1) applies to the leptons only. The proton ratio holds at both faces that 155.1 equates ($2871/67=8613/201$).
- **τ value.** Now evaluated with PDG 2024: $+0.917\%$.
- **Grading.** 2871 is *entered* as a candidate, the integer nearest $67\sqrt{m_p/m_e}$. The proton agreement is a selection of reading, and the lepton agreements are conditional on a face choice the compendium's own definitions exclude. The only face-free content is $m_\tau/m_\mu=16.987$, which is $+1.0\%$ ($198\sigma$).
