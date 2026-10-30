# The hydrogen body registers and the ladder beyond n = 8 (round 5)

Source: compendium 0.1, commit `a28db899`. Script: `computation/chain-theory/compendium/hydrogen_registers.py`.
Recorded output: `outputs/hydrogen_registers.txt`; block names in brackets below refer to it. Metrics:
`metrics/chain-theory/compendium/hydrogen-registers.yaml`. Object: item 5 of the proposed further work in
`reports/chain-theory-compendium.md`. Round 1 is `02-hydrogen-spectrum.md` (H1–H4 and the reconstruction R).

The script checks all 53 quoted strings verbatim against their statement units or files [ANCHORS]. Grades are as printed:
T THEOREM, C CONDITIONAL, O OPEN, Tab TABLE. 143.1 is C.

Status codes used below:
- **F**: follows from a definition or rule the source states.
- **R**: one of several readings. Each reading is listed with the number it gives.
- **E**: entered, with no formula in the source.

## 0. Verdict

| item | status | one-line reason |
|---|---|---|
| whole register 2^30·3^18 | F given the electron web, plus one reading | It is the product of (d_a d_b)² over the 10 links of the 152.1 electron web: the 16.1/103.1 whole face combined with the 103.1/147.1 product rule. The struck docks must be counted. |
| cut register 2^31·3^17 | R | It is the same product with only the port B–w1 at its cut face (2²−1)(3²−1) = 24. Read literally at the cut face, the 103.1 rule gives 2^27·3^11·5^4. The port-only reading is the one behind every printed cut register we could test (2 bodies, 3 numbers). |
| units 67 and 201 | F given the web | 67 is the electron rank 79 − 12, and 201 is the rank times the hub product 3 (152.1; round 3). The source prints which unit goes with which face but gives no reason. |
| mass face M_L = 36 | E | No formula is printed. Two stated quantities equal 36; the displayed held chain has measure 30. |
| ladder beyond n = 8 | nothing printed | No value, formula or rank beyond 8 appears anywhere in the release. R gives exactly 1/n² at every n. |
| comparison with experiment | identification missing | The source itself says the mass unit is open (192.1). Under X/m_e c², E1 is 9.18 or 9.49 orders low. E_2/E_1 = 1/4 agrees with the measured ratio to 4.5×10⁻⁶. |

## 1. The statements (task a)

The binding (143.1 derivation, `07-cell-3-2-photon.tex:504-505`): "by the read through the electron body at its dock
port $B$-$w_{1}$ (the mass face multiplied by the unit over the register of the body at that face)". It continues: "the
whole face with the hub unit two hundred and one over the body register two to the thirtieth multiplied by three to the
eighteenth, and the cut face with the dock unit sixty-seven over the body register two to the thirty-first multiplied by
three to the seventeenth (the port at its cut twenty-four)". The page adds (line 512): "The two body registers of the
level-one binding are cells of the read tables of this result, printed in their prime factorizations and re-derived on
the page not".

The source states the rules that produce the registers in other chapters:
- **16.1 (T), derivation** (`01-framework.tex:582`): "From the whole register $f^{2} c^{2}$ and the cut register
  $(f^{2} - 1)(c^{2} - 1)$". The derivation prints the count register cut·whole² at (2,2) and (4,4): 2304 and 14745600.
- **103.1 (C)** (`06-cell-2-2-metric-and-horizon.tex:814-845`): "the register of a body at a face the product over its
  links of the per-link face". It also defines four per-link faces: "the whole, the square of the product of the ends,
  the cut, the product of the squared ends less one each, the rank, the product of the ends, and the support". It prints
  them for four genres, 16 numbers in all.
- **147.1 (C)** (W word): "The product register of the body at the whole port face is the product over its four links of
  the squared measures", with $R = (8 \cdot 6 \cdot 6 \cdot 8)^{2} = 5308416$. The cut chain sum is 22, from a measure of
  28.
- **149.1 (C)**:
  - Z body: "the registers are 5308416 at the whole face and 3538944 at the cut face, the chain sums 28 and 22, and the
    factors $7/1327104$ and $11/1769472$".
  - Gluon word, as a counter-test: "moves the factors to $1/19683$ and $1/23328$".
- **122.1 (O)**: "the whole register $36=4\cdot 9$ and the cut register $24=3\cdot 8$". **164.4 (O)**: "At the
  qutrit-qubit genre the whole 36, the cut 24, the common-support register 16". The 12-6 derivation (unit 214.6) gives
  the support as "the smaller end of the port genre to the fourth power".

**Result A0 (per-link faces; c).** whole (fc)², cut (f²−1)(c²−1), rank fc and support min(f,c)⁴ reproduce all 16 faces
printed in 103.1. They also reproduce 16.1's count registers: 2304 at (2,2), 14745600 at (4,4), and the 147.1 star 31104
= 24·36² at (2,3) [FACES].

## 2. The counting steps

| # | step | formula | value | status |
|---|---|---|---|---|
| 1 | electron web | parsed from the 152.1 figure | 9 sites, 10 links: core (3,3)×3, in-leg (3,4), hub legs (2,4)×3, qubit pair (2,2), docks (2,3)×2 (B–w1, C–w2) | **E**: the roster's electron (round 3, `06-mass-cells.md` §1 row 1) |
| 2 | per-link whole face | (d_a d_b)² | 81 ×3, 144, 64 ×3, 16, 36 ×2 | **F** (16.1, 103.1) |
| 3 | whole body register | R_w = ∏_ℓ (d_a d_b)² = ∏_s d_s^{2 deg s} | **2^30·3^18** | **F** given 1 (103.1, 147.1). One reading: the struck dock links are counted (W0, printed). Counting the 8 unstruck links only (W1) gives 2^26·3^14. The port B–w1 is itself a struck dock, so W1 has no port to read. |
| 4 | per-link cut face | (f²−1)(c²−1) | 24 at (2,3) | **F** (16.1, 103.1, 122.1, 164.4) |
| 5 | cut body register | R_c = R_w · cut(port)/whole(port) = R_w·24/36 | **2^31·3^17** | **R**: C1 (port only, printed); see Result A2 for the others |
| 6 | units | dock 67 = 79 − 12; hub 201 = 67·3 | 67, 201 | **F** given 1 (152.1; round 3 rows 2–5) |
| 6′ | face–unit pairing | whole ↔ 201, cut ↔ 67 | — | **R**: printed without a reason; the cross pairings give 67/(2^28·3^16) and 67/(2^29·3^14) |
| 7 | mass face | M_L | 36 | **E** (Result A4) |
| 8 | binding | E1 = M_L · u / R | 67/(2^28·3^15), 67/(2^29·3^15); port readings 67 and 201/2 | **F** (143.1's read rule); 4 of 4 printed equalities reproduced exactly [E1] |

**Result A1 (whole register; c).** With the 152.1 web, 81³·144·64³·16·36² = 2^30·3^18 [ELECTRON]. Given the web and the
one reading of step 3, the number is forced by rules printed in 16.1, 103.1 and 147.1. It equals ∏_s d_s^{2 deg s}, the
dimension of the tensor product over the links of End(C^{d_a} ⊗ C^{d_b}).

**Result A2 (cut register; b, c).** The printed 2^31·3^17 is the whole register with the port alone at its cut face.
The electron cut-face readings give:

| reading | rule | register |
|---|---|---|
| C1 | the port at its cut, every other link whole | 2^31·3^17 (printed) |
| C2 | every link at its cut (the 103.1 sentence read literally) | 2^27·3^11·5^4 |
| C3 | both (2,3) docks at their cut | 2^32·3^16 |
| C4 | the port at its support face 16 | 2^32·3^16 (equal to C3, since (24/36)² = 16/36) |
| C5 | the port at its rank face 6 | 2^29·3^17 |

C1 is also the rule behind the cut values the source prints for the two other bodies we could test [BODIES]:
- **Z body (4,2,3,2,4).** 5308416·24/36 = 3538944, and 22/3538944 = 11/1769472. These are the printed register and
  factor.
- **Gluon word (four qutrits).** 18/(64·81·81) = 1/23328, the printed factor. This needs the cut chain sum
  S_cut = 27 − 9 = 18, the 147.1 pattern (28 → 22). Among integer pairs (S_cut ≤ 27, port factor ≤ 81) only
  (9, 32) and (18, 64) fit, and 32 is not a face of 103.1.

C2 and C3 reproduce none of these three numbers. The source states the port-only rule nowhere; we infer it from 2 bodies.

**Result A3 (the factor 2 between the faces; c).** E1(whole)/E1(cut) = (201/67)·(24/36) = 3·(2/3) = 2. It is the
electron's hub dressing times the (2,3) cut/whole ratio. Neither factor involves the photon word, the proton or n. More
generally, every input of E1 other than M_L is a cell of the electron body alone.

**Result A4 (mass face; a, c).** The source prints "the mass face of the held chain thirty-six and the exchanged mass face
six root six" with no formula for either [CARRIER]:
- The displayed held chain ("6 sites, 5 links of type (2,3)", the 143.1 figure) has measure 30.
- "The eight links of the held chain" (lines 488 and 511) would give 48.
- The docking act as described ("one link made at each end and one contraction at one link") leaves 6 links of type
  (2,3), measure 36.
- The whole face of one (2,3) link is 36. 192.1 lists "the whole rung register as thirty-six" as one of its mass-unit
  spellings.

So 36 equals two stated quantities, and the text names neither. The text counts the links of the held chain or the bond
three ways: 5 (display), 8 (fine-structure bullet, census) and 1 (glossary, "the link that joins the electron body to
the proton body").

**Result A5 (face selection elsewhere; a).** 229.1 (T) says "the hub-or-dock question of the previous releases resolves
to the dock wherever the joined read opens at all", with "the register 16". That selects the dock unit, but at a
register that matches neither printed face:
- 36·67/16 = 603/4 as a port reading;
- 67/(2^30·3^14) with the port at its support face in the body register [E1].

## 3. The ladder beyond n = 8 (task b)

**Result B1 (what the source prints; c).** We scanned every sentence of the tex that states the hydrogen ladder
[LADDER (1)].
- There are 15 such sentences, counting section titles and the appendix of conditions. 6 state the range "eight" or "1, …, 8", and none names a larger rank.
- No ladder value or formula is printed for n > 8.
- The per-record ladder has no printed value at any rank.
- In Lean the ladder appears only as a record name, "RYDBERG-LADDER-FROM-THE-PER-RECORD-READ"
  (`ConditionUnions.lean:119`).

So R cannot be compared with the source beyond n = 8. The source enumerates "the eight ranks of the held chain" and does
not say whether a ninth rank exists.

**Result B2 (R at every n; c, proven).** Take m = 2.
- The commutant of the sl2 action on Sym^(n−1)(C²) is C·I. The diagonal of H is distinct and the support of E is
  connected. This is checked exactly for n ≤ 200 and by brute force over Q for n ≤ 5, and the argument holds at every n.
- So the twirled state is I/n, one read gives 1/n, and the two-vertex read gives exactly 1/n² at every n.

For m ≥ 3, C(n+m−2, m−1) > n when n ≥ 2. The ladder then equals 1/n² only at n = 1 (checked for m = 3, 4, 6 and
n ≤ 200) [LADDER (2), (3)]. R agrees with the 8 printed values and continues them with no cap.

**Result B3 (where R and the source can part; b, c).** Two sentences say more than the two-vertex values.
1. **The per-record sentence.** 143.1 says "the per-record ladder is exact at every rank and no pure power of the rank",
   and the census says "a pure power of the rank at none". Under R [LADDER (5)]:

   | reading | value under R | pure power of n | verdict |
   |---|---|---|---|
   | P1, one read | 1/n | at every rank | parts |
   | P2, second moment ⟨J_z²⟩ of one traceless letter | (n²−1)/12 | at n = 2 only (1/4 = 2⁻²) | parts at that rank |
   | P0, first moment | 0 | never | agrees, trivially |

2. **The all-orders sentence.** The 3-7 prose after 45.1 says "the order-$k$ term at rank $n$ the term at rank one over
   $n$ squared at all eight ranks". Under R, a term that reads the state r times scales as n^(−r). So R agrees for every
   k only if every order-k term at the held pair reads the state twice [LADDER (4)].

The source defines neither "record" nor "order-k term" at the held pair, so neither sentence decides R.

## 4. Comparison with experiment (task c)

**Result C1 (the identification; a).** The source identifies E1 by name only:
- the box: E1 is "an exact element of the field in the electron unit";
- the bullet: "the level-one binding is the two-vertex term at the bond, the Rydberg constant in the electron unit the
  ladder cell of the held carrier";
- the bridge: "The correspondence is by name".

The unit itself is left open by the source:
- 192.1 (Tab) ends "the mass unit stays a question of which cell is meant". Its five spellings are 1, 36, 16, 67 and 67.
- The appendix of conditions lists "The Rydberg constant as the two vertex term at the bond" as resting on, among others,
  "Mass unit as the rung register (an open question)" (`back.tex:63`).

So the map needed to compare E1 numerically is missing. The one usage that fixes a meaning is 176.1's "The proton in the
electron unit", a ratio to m_e.

**Result C2 (under X/m_e c²; c, e).** We use CODATA 2022 and the NIST ionization energy 13.598434599702 eV [COMPARE].

| face | E1 / (IE/m_e c²) | log10 gap | E1 / (α²/2) |
|---|---|---|---|
| whole | 6.5365e-10 | 9.185 | 6.5330e-10 |
| cut | 3.2683e-10 | 9.486 | 3.2665e-10 |

This reproduces round-1 H2 with new code.

**Result C3 (declared family; c).** The family takes the four printed binding numbers {E1 whole, E1 cut, 67, 201/2}. Each
is divided by u^p, with u ∈ {1, 16, 36, 67, 201} and p ∈ {1, 2}, or left undivided: 36 variants. The u are the 192.1
spellings plus the hub unit.
- None is within a factor 10 of IE/m_e c².
- The nearest is 67/201², 62 times the measured value.
- The family spans 10^−18.7 to 10^2.

No conversion built from the source's own unit spellings brings E1 near the measured scale.

**Result C4 (the unit-free ratio; c, e).** E_2/E_1 = 1/4 does not involve the unit. We compare it with
B(2S)/B(1S) = 1 − ν(1S–2S)/IE, using the 1S–2S hyperfine centroid 2466061413187035(10) Hz:
- The result is 4B(2S)/B(1S) − 1 = +4.537×10⁻⁶, with a measurement uncertainty of 2.7×10⁻¹².
- The Dirac energies alone give α²/16 = 3.33×10⁻⁶. The remaining 1.2×10⁻⁶ comes from the QED and recoil terms beyond
  them (the Lamb shifts among them; not evaluated here).

So the exact 1/4 differs from the measured ratio by an amount of the size of the fine-structure and Lamb terms, which
143.1 states without values. At the level the source evaluates, its ladder agrees with experiment as the Bohr formula
does, and it cannot be tested further. The Lamb sign and the hyperfine reading 0 are as in round 1 (H4).

## 5. Limits

- The electron web is the roster's (round 3). A different web moves both registers.
- The port-only cut (C1) is inferred from 2 bodies and 3 printed numbers. The gluon check uses S_cut = 27 − 9, itself
  inferred from the 147.1 pattern.
- R is our reconstruction (02 H3). B2 and B3 are void if the source's ladder is a different object. The per-record
  readings P0–P2 are ours.
- C2 and C3 rest on the reading "in the electron unit" = X/m_e c². C3's family is ours and declared above. C4 compares
  only n = 2, the one level with a sub-ppm measurement used here.

## 6. Changes to earlier rounds

- `02-hydrogen-spectrum.md` §1 lists the body registers among the "inputs not derived on the page". That stays true of
  the page. The registers do follow from rules the release states elsewhere, given the web and the readings of steps 3
  and 5.
- Choice (5) of §4, "the two body registers", should now be read as those two readings.
- Choice (3), "none selected": 229.1 later resolves the unit to the dock, at a register that matches neither printed
  face. So no printed face is selected.
- No earlier statement is wrong, so nothing is added to `reports/CORRECTIONS.md`. G1 and G2 in `gaps/` get dated notes.
