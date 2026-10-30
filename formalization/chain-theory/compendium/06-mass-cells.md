# The mass cells: derivation or selection? (round 3)

Source: compendium 0.1, commit `a28db899`. Script: `computation/chain-theory/compendium/masscells.py` with the
transcription module `masscells_data.py`. Recorded outputs: `outputs/masscells.txt` (all numbers below) and
`outputs/masscells_readings.tsv` (the 72-reading grid). Metrics: `metrics/chain-theory/compendium/mass-cells.yaml`.
Object: scorecard rows SC-35, SC-36 and SC-38 of `05-numeric-scorecard.md`, i.e. m_p/m_e ≈ (2871/67)² (+1.84×10⁻⁵),
m_μ/m_e ≈ (963/67)² (−0.088%) and m_τ/m_e ≈ (3969/67)² (+0.92%).

References: CODATA 2022 m_p/m_e = 1836.152673426(32) and m_μ/m_e = 206.7682827(46), read from `scipy.constants` 1.18.1
and asserted against the transcription. PDG 2024 lepton masses come from the MC table `mass_width_2024.mcd`, checked at
run time. They give m_τ/m_e = 3477.365(176) (CODATA 2022 has 3477.23(23)), m_τ/m_μ = 16.81769(85) and the Koide ratio
Q = 0.6666645(51).

Every quoted string is checked verbatim against its statement unit by the script: 102 strings, 0 missing. The grades
of the statements cited are T THEOREM, C CONDITIONAL and O OPEN, as printed.

## 0. Verdict

| item | grade | one-line reason |
|---|---|---|
| integers 67, 963, 3969 | derived, given the webs | rank × 3^(in-legs) of the electron, muon and tau webs, re-derived from the 152.1 figure and the printed descriptions; the tau web is a "fixed library" web that the saturation never reaches (189.1 O) |
| integer 319 | derived, given the α weight law | it is the (2,3) Born count, re-derived combinatorially: 79·1 + 48·2 + 42·3 + 3·6 |
| cell 2871 = 319·9 | entered as a candidate | the text calls it "the candidate cell" and "the proton figure of the comparison table"; the rule that yields it (186.1) is tested against it; it is the integer nearest 67√(m_p/m_e) |
| m_μ/m_e, m_τ/m_e agreements | conditional on one choice the compendium's own definitions exclude | they need the electron read at the dock face and μ, τ at the hub face; under 198.3, 185.6 and 155.1 both ratios are 9× low, as 176.2 prints |
| m_p/m_e agreement | selection of the reading | the readings that meet 2871 were adopted because they meet it (20.1, 186.1); the claimed precision lies inside the spread of the source's own mass laws |
| all three at once | no single definitional reading | 1 of 72 uniform readings gives all three; it is the reading of 186.1, whose test is equality with these cells |

## 1. Trace of each integer (task 1)

Status codes: **F**, forced by a stated definition; **R**, one of several readings the text prints; **E**, entered
(an identification, or a value introduced without derivation).

| # | step | value | where | re-derived here | status |
|---|---|---|---|---|---|
| 1 | electron web | 9 sites, 10 links | 152.1 (C) figure; 195.1 | parsed from the picture: every genre label equals (dim, dim) of its ends | **E**, the roster's electron |
| 2 | electron measure | 79 = 27 + 40 + 12 | 152.1; measure = Σ d₁d₂ (185.12) | 79 | **F** given 1 |
| 3 | dock strike | 12 (two (2,3) docks) | 152.1, 185.12 | 12 | printed per body. 155.1 (C) states that the strike rule "reproduces the struck set at no tabulated lepton result"; 32.1 (O) derives it "at the pass line". **R/E** |
| 4 | rank | 67 = 79 − 12 | 152.1 | 67 | **F** given 2, 3 |
| 5 | electron first power | 67 (dressing 1) or 201 (dressing 3) | 186.1 (T) and 195.1 (T) path reading: 67. 185.12 (C), 152.1 (C) and 198.3 (C) dressing law: 201 | breadth-first path: the one in-leg lies on the geodesic | **R**; the definitional law 198.3 gives 201 |
| 6 | muon | 963 = (117 − 10)·9 | 152.1 (rank 107), 185.12, 186.1 | web built from 152.1 and 195.1: 13 sites, 15 links, measure 117, strike 10, dressing 9 | **F** under 198.3. Under the path reading, the tie needs "both in-legs unpassed" (186.1); 20.1 resolves the same tie "both passed" and gets 69 |
| 7 | tau | 3969 = (155 − 8)·27 | 152.1, 185.12, 186.1 | 17 sites, 20 links, 155, 8, 147, 27 | **F** given the web. The web is **E**: "the tau of the library is a web the census reaches at no hub-count-two and no hub-count-three sink" (189.1 O). The census tau sink stands at 4185 (154.1); 20.1's adopted reading gives 432 |
| 8 | 319 | Born count of the (2,3) link | 176.2: "319 as the sum over the 172 record letters" | 172 records; closed 79 = 48 + 28 + 3 (as printed); Σk = 319 | **F** given the α Letter's weight law. Over the α audit's 4608 squared-factor rules the Born count takes 517 values (169 to 1068); 319 occurs at 19 rules |
| 9 | proton records | 319 (capped entry, counted at its cap) | 32.1 (O), 186.1, 175.4 | proton seed: 10 sites, 12 links, 101 = 27+18+36+16+4 | **R**: the uncapped typing gives (101 − 18)·9 = 747 (179.1). 32.1: "each still owing a derivation" |
| 10 | proton strike | 0 | 186.1: "struck at an uncapped entry alone" | — | **R**: 20.1's adopted reading strikes 18 and gives 2709 (−11%). 186.1: "every meeting reading at the strike and the beyond records adopted here" |
| 11 | proton dressing | 9: two of three in-legs unpassed | 20.1, 186.1 | breadth-first path: 1 of 3 in-legs on the geodesic | **R**: 27 under 198.3 and the 21.1 stop clause. The ratio is unchanged if the electron is read at the same face |
| 12 | proton cell | 2871 = 319·9 | 186.1 (T) "returns the derived cells". 20.1, 21.1, 175.4: "the candidate cell". 175.3: "no sink … reads the proton figure of the comparison table at any reading tried, so the candidate is refuted". 175.4: "the two proton cells of the construction". Chapter VII-5: "a candidate stated as a candidate" | 67√(m_p/m_e) = 2870.9736; nearest integer 2871 | **E**; its provenance is not documented |
| 13 | mass law | m = S² | 185.9 (C) leading term; with its body term, = 195.1 office S² + closure | — | **R**: S (181.5 and 200.6 bridges; 200.3 Higgs exponent one), S², S² + closure, S²(1 + G₂²) (176.2) |
| 14 | electron unit in the ratio | 67 | — | — | **R**, from the printed set {67, 201, 79, 237} |

The source states the split itself. 20.1 (O): "Over the grid of sixty readings … no one reading meets the four derived
cells at once … The proton is met at six readings, every one at the capped register of its entered link with nothing
removed, … so the passed records at the leptons and the passed leg at the proton are objects of two rules." 186.1
(THEOREM) then adds two clauses, "the dock set is struck at an uncapped entry alone" and "a capped link beyond the
first cap reads nothing". Its "test applied to each read is the equality of this sum with the derived cell of the
body"; it "meets every cell at 8 of its 36 readings".

## 2. The space of readings and the distribution (task 2)

Notation:
- "Variant": one choice of (numerator cell, electron cell, exponent p ∈ {1, 2}).
- δ: the achieved |ln(pred/meas)|. It is 1.84×10⁻⁵ for p, 8.77×10⁻⁴ for μ and 9.12×10⁻³ for τ.
- ρ: distinct predicted values per unit ln within |ln(v/meas)| < w, for windows w = 0.5 and w = 1.0.
- λ = 2δρ, the expected number of chance values that close. P = 1 − e^(−λ). This is the F-CT-4 / round-2 method.

The spaces:
1. **Printed.** The 15 proton first powers the text prints for the proton's own webs: 2871, 8613, 2709, 2727, 101,
   2187, 3108, 25839, 729, 747, 1701, 255, 85, 201, and S_p of 214.6. The 4 printed electron cells: 67, 201, 79, 237.
   Each is anchored verbatim. A mechanical sweep of the integers in proton first-power contexts leaves 0 unclassified;
   the exclusions are listed with reasons in the output.
2. **Reading grid.** We reconstruct the axes the text names, and S = (records − strike)·3^u:
   - 4 records typings: the reach (32.1/186.1), the register (20.1), and the measure with every link uncapped (seed or
     sink).
   - 3 strike rules: every entry (20.1), the uncapped entry alone (186.1), none.
   - 6 dressing rules: the path with the tie passed (20.1), the path with the tie unpassed (186.1), all in-legs
     (198.3, the hub face), all but one (155.1, the dock face), the stop clause (21.1), none.

   The grid has 72 readings. It reproduces all 6 tuples the text prints: 20.1 (67, 69, 432, 2709), 21.1 (67, 621, 432,
   8613), 186.1 (67, 963, 3969, 2871), the 185.12 hub leptons (201, 963, 3969), 179.1 (747) and the 195.1 mass-register
   row (201, 963, 3969, 1701). It also reproduces both of 20.1's counting statements: 6 readings meet 2871 under the
   20.1 strike rules, and none meets all four cells.
3. **Round-2 declared space.** 330 cells × {67, 201} × p, recomputed with the same code.
4. **Narrow family.** (319·3^u/(67·3^v))^p, with 319 and 67 fixed upstream.
5. **Saturation censuses.** 154.1: 63 distinct tau-candidate first powers between 369 and 4617. 189.1: 92 two-hub
   values between 95 and 3321. Only the counts and ranges are printed, so log-uniform spacing is assumed.
6. **Any integer.** The numerator over e = 67 may be any integer, or any multiple of 9.

**m_p/m_e.** (2871/67)² ranks first in every space.

| space | variants (distinct) | next closest, \|ln dev\| | λ at δ = 1.84e-5 (w = 0.5 / 1.0) | λ at the convention spread 1.2e-4 |
|---|---|---|---|---|
| printed | 120 (103) | 0.103 ((2727/67)², −9.8%) | 2.9e-4 / 1.8e-4 | 1.9e-3 / 1.2e-3 |
| reading grid, cell chosen per body | 144 (120) | 0.103 | 2.6e-4 / 1.8e-4 | — |
| reading grid, uniform readings | 144 (34) | 0.116 | 1.1e-4 / 5.5e-5 | — |
| round-2 declared | 1320 (874) | 0.0125 | 2.5e-3 / 2.3e-3 | — |
| narrow family | — | — | 1.7e-5 (p = 2 fixed); 5.0e-5 (p free) | — |
| any integer over 67 | — | — | 0.053 | 0.34 |
| any multiple of 9 over 67 | — | — | 0.0059 | 0.038 |

In the printed space only 2 of 120 variants lie within 10% of the measured ratio, and both are the same ratio
(2871/67 = 8613/201). The median variant lies a factor e^4.91 ≈ 136 away. With the electron cell free in the printed
set, the any-integer chance at δ = 1.84×10⁻⁵ rises to 0.39. The scorecard compares 34 quantities numerically, so a
trials factor of up to 34 applies across quantities.

**Leptons.**

| ratio | printed (λ, w = 0.5 / 1.0) | per body over the grid | census | integer level |
|---|---|---|---|---|
| m_μ/m_e | 3.5e-3 / 2.6e-3 (24 variants) | 7.0e-3 / 6.1e-3 | 0.023 | 963 is the integer nearest 67√(m_μ/m_e) = 963.42; for any integer P = 0.84, for a multiple of 9 P = 0.094 |
| m_τ/m_e | 0.15 / 0.073 (40 variants) | 0.073 / 0.037 | 0.23 | 3969 lies 18.1 above the target cell 3950.9; the nearest multiple of 27 is 3942 |

In the printed space both lepton ratios fall within 1% for 1 of the 60 (e, μ, τ) combinations: (67, 963, 3969).

**The precision depends on the mass law.** The table gives the agreeing cells under each printed law
(`[MASSLAW]`).

| mass law | m_p/m_e | rel. | m_μ/m_e | m_τ/m_e |
|---|---|---|---|---|
| S² | 1836.186456 | +1.84e-5 | 206.5870 | 3509.24 |
| S² + closure (195.1 office, T; 185.9 with its body term) | 1836.164052 | +6.20e-6 | 206.5864 | 3509.22 |
| as above, with the proton closure at the mass-row weight (O.1) | 1836.164054 | +6.20e-6 | 206.5864 | 3509.22 |
| S²(1 + G₂²) of 176.2, with G₂ at e, μ and the (2,3) proton entry | 1836.186456 | +1.84e-5 | 206.5870 | 3508.82 |
| S²(1 + G₂²), with G₂ at e and μ only | 1835.966478 | −1.01e-4 | 206.5870 | 3508.82 |
| S (office = m²) | 42.850746 | −97.7% | 14.37 | 59.24 |

The printed second-order laws spread m_p/m_e over a relative width of 1.20×10⁻⁴, and the CODATA value lies inside that
band. The source's own next terms therefore exceed the claimed agreement by up to a factor 6.5, and nothing in the text
selects among them. As exact statements, the S² and office readings miss by +1.06×10⁶σ and +3.56×10⁵σ (p),
−3.9×10⁴σ (μ) and +181σ (τ).

## 3. Internal consistency (task 3)

**Faces** (`[FACES]`; 155.1 names "the dock reading", "the hub reading" and "the mixed face the hub numerator over the
dock unit").

| face | numerators (e, μ, τ, p) | m_μ/m_e | m_τ/m_e | m_p/m_e |
|---|---|---|---|---|
| dock / dock | 67, 321, 1323, 2871 | 22.954 | 389.92 | 1836.19 |
| hub / hub (198.3) | 201, 963, 3969, 8613 | 22.954 | 389.92 | 1836.19 |
| mixed: hub / dock unit | 67 unit; 963, 3969, 8613 | 206.587 | 3509.24 | 16525.7 (×9) |
| 186.1 path reading | 67, 963, 3969, 2871 | 206.587 | 3509.24 | 1836.19 |

1. **No face gives all three agreements.** The proton ratio is face-invariant (2871/67 = 8613/201), so it holds at
   both readings that 155.1 proves equivalent. The lepton agreements hold only at the mixed face, where the electron
   absorbs a leg and μ and τ do not. The 186.1 path reading absorbs one in-leg at e (1 of 1) and at p (1 of 3), and none
   at μ (0 of 2) or τ (0 of 3). That is, it is the dock face for e and p and the hub face for μ and τ. The dock-face
   numerators 321 and 1323 are printed nowhere in the release.
2. **Grid** (`[READINGS]`, p = 2).
   - The four cells together occur at 1 of 72 readings: the reach, the 186.1 strike clause and the tie broken as
     unpassed. Without the 186.1 strike clause they occur at none, as 20.1 itself says.
   - Ratios within 1%: μ and τ together at 12 readings, every one with e = 67 and μ, τ at hub numerators. p/e at 8
     readings, every one with the 186.1 strike clause, including the 198.3 hub face. All three at 1 reading. At p = 1,
     none.
   - The 195.1 mass-register row is the uniform reading that applies the 198.3 dressing law, one strike rule and one
     typing to every body. It prints (201, 963, 3969, 1701) and gives 22.95, 389.9 and 71.6: no agreement at all.
3. **Conflicts with stated results.**
   - (a) The lepton agreements need an electron dressing of 1. 198.3 (C) defines the dressing as "a function of its
     hub in-legs alone"; that gives 3 at the electron, as do 185.6 (C), 185.12 (C, S₁ = 201) and 152.1 (C). 152.1 says
     "the two offices parted by the factor nine at every body". 155.1 (C) says "the in-leg count is one function of the
     web alone, invariant under every reading" and that "the ratio form carries no reading". Under these the ratios are
     22.95 and 389.9, which are the values 176.2 prints. The THEOREM 186.1 and the "natural checkpoint" of THEOREM 195.1 use the path dressing (67); 195.1's own
     mass-register row uses the 198.3 dressing (201). The compendium's theorems and its conditional definitions therefore
     disagree on the electron's first power, and the agreements need the path version.
   - (b) The muon tie is resolved oppositely in 20.1 ("both passed") and in 186.1 ("both in-legs unpassed").
   - (c) The proton agreement needs the capped-entry typing (32.1 O, "owing a derivation") and the strike exemption at
     a capped entry. 20.1's adopted reading strikes the proton and gives 2709.
   - (d) 176.1 (C), "stands outside the comparison by orders at every reading the construction has produced", holds
     only under m ∝ S. Under 185.9 it contradicts 186.1 and 195.1.
   - (e) **Cross-check: the same rule at other bodies** (`[CROSS]`). 181.5 (THEOREM) applies the same segment walk
     with a capped entry to the quark sinks. It gives the top sink 2871 = 319·9, so m_t = m_p at either exponent,
     against a measured m_t/m_p = 183.9. It also orders the heavy quarks s > c = b > t.

## 4. Grade and falsifiable content (task 4)

- **m_p/m_e: selection** of the reading.
  - The readings that produce 2871 were adopted because they meet it: 20.1 says "the proton is met at six readings",
    and 186.1 says "every meeting reading at the strike and the beyond records adopted here".
  - The cell was introduced as "the candidate" and "the proton figure of the comparison table". The order in which the
    cell and the target were fixed is not documented, as with CT-G10 for α.
  - The 1.84×10⁻⁵ is smaller than the source's own spread of mass laws (1.2×10⁻⁴).
  - One coincidence remains. The integer nearest 67√(m_p/m_e) equals 9 × 319, where 319 is fixed upstream by the α
    weight law. Its chance is about 2.5×10⁻³ in the round-2 ingredient space and about 2×10⁻⁴ in the printed space,
    before the trials factor across quantities. It would be notable if the cell had been fixed before the comparison;
    the text gives no evidence either way.
- **m_μ/m_e and m_τ/m_e: conditional on choices.** The integers are derived, but the agreement needs the mixed face,
  which 198.3 and 155.1 exclude. Under the definitional law both ratios are 9× low (8.9× for τ). The only
  reading-independent content is m_τ/m_μ = (441/107)² = 16.987. That is +1.0% from 16.81769(85), or +198σ. The muon cell
  is the integer nearest the target at e = 67, so its 0.088% agreement is what an integer numerator gives typically
  (P = 0.84). The chance that some power-of-3 rescaling brings a fixed prediction this close is 1.6×10⁻³ (8.0×10⁻⁴ for powers of 9).
- **Joint:** no reading consistent with the compendium's definitions gives the lepton and the proton agreements
  together. Each agreement uses a different face.
- **Falsifiable content.**
  - As exact equalities, every one is excluded: p at ≥ 3.6×10⁵σ, μ at 3.9×10⁴σ, τ at 181σ, τ/μ at 198σ and Koide at
    196σ.
  - As approximate relations they carry no error budget. The source's own printed next terms move m_p/m_e by up to
    1.2×10⁻⁴, and none is stated to be the last.
  - Two sharp predictions are independent of the reading: m_τ/m_μ = (441/107)², which fails at 198σ, and, through
    181.5, m_t = m_p, which fails by a factor 184.
  - The proposal would become testable with a rule that fixes the face, the typing and the strike before comparison,
    together with a bound on the uncomputed terms.

## 5. Choices made in this audit

1. Segmentation and normalisation are those of round 2 (`scorecard_extract.units`, `norm`).
2. The muon, tau and proton-seed webs are built from the printed descriptions (152.1, 175.3) and the 195.1 rows
   (`masscells_data.WEB_NOTES`). The placement of the (2,2) pair inside each hub assembly is ours. Measure, site and
   link counts and the paths do not depend on it, and all match the printed values.
3. The 72-reading grid is our reconstruction of the axes the text names. 20.1's own 60-reading grid is not specified
   in enough detail to be reproduced. Our grid is validated against 6 printed tuples and 2 printed counts.
4. Membership of the printed space is curated, with verbatim anchors. The neutron, the deuteron, the 116.1
   insertion readings and the quark sinks are excluded, each with a stated reason.
5. λ uses the local log-density, as in F-CT-4 and round 2. These values are chance-hit expectations under the
   assumption that the cells were fixed without the target. They are not p-values for a derivation.
6. The "convention spread" is the width of the printed second-order mass laws at the agreeing cells.

## 6. Effect on round 2

These are logged in `reports/CORRECTIONS.md` (2026-10-30, round 3).

- SC-38's "+1.84×10⁻⁵" is the value under m = S² only. The source's own THEOREM 195.1 offices give +6.2×10⁻⁶, and
  176.2's second-order law can give −1.0×10⁻⁴.
- The printed proton cells number 15, not 9. One distinct cell is within 1%, as before.
- Round 2's face-violation remark (L, F11) applies to the leptons only. The proton ratio does not depend on the face.
- The τ reference is now PDG 2024: +0.917%, 181σ.
- The 3969 web is the library tau, which the saturation does not reach (189.1).
- Round 2's declared-space figure 0.0025 is reproduced (2.5×10⁻³ at w = 0.5).
