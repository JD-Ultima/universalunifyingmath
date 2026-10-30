# CT-C016: "the hydrogen spectrum" (statement 143.1, with 144.1)

Source: compendium 0.1, commit `a28db899`. Abstract (`tex/front.tex:25`): "the hydrogen spectrum" as a
conditional theorem. Box 143.1 "Hydrogen from the held carrier", `tex/parts/07-cell-3-2-photon.tex:481-484`;
derivation and bullets lines 488-541; open result 144.1 (hyperfine) lines 551-568. Computations:
`computation/chain-theory/compendium/hydrogen_audit.py` (output `outputs/hydrogen_audit.txt`).

## 1. What the source asserts (a)

Box: "Hydrogen is the held carrier: the photon word held at both ends, the qutrit end at the electron
dock and the qubit end at the up core of the proton, is a closed read with a mass reading, and its
level-one binding through the electron body is an exact element of the field in the electron unit."
Grade CONDITIONAL. Direct premises include 124.1 (alpha) and 141.1 (CONDITIONAL), open results 18.5,
18.6, 28.3, 28.4, 57.2 and H.1, H.2, H.4-H.8, and CONDITIONAL H.9 (statements.tsv).

What is actually derived on the page, item by item:

| Item | Printed result | Status on the page |
|---|---|---|
| level-one binding | E1(whole) = 36 * 201/(2^30 3^18) = 67/(2^28 3^15); E1(cut) = 36 * 67/(2^31 3^17) = 67/(2^29 3^15) | "both faces printed"; the two body registers "re-derived on the page not" |
| Bohr ladder | E_n = E1/n^2, n = 1..8 ("two-vertex ladder"); the per-record ladder is "a pure power of the rank at none" | enumerated to n = 8 |
| fine structure | "the two-vertex term multiplied by the square of the cut gain, conditional on the admission of the second-order read at the (2,3) link, a posed step" | no value; cut gain not given |
| Lamb shift | "one exact element, positive and split" | no value |
| hyperfine | open result 144.1: interval = 0 at the standing transport, undefined at the ray-cycle reading | OPEN |
| m_p/m_e | "an exact rational ... compared at the comparison table"; 176.1: "the proton in the electron unit stands outside the comparison by orders" | not displayed |

Inputs not derived on the page: the photon word (6-site (2,3) chain), the docking act at both ends, the
mass face 36, the hub unit 201 and dock unit 67, and the body registers 2^30 3^18 and 2^31 3^17.

## 2. Results

**Proposition H1 (arithmetic; proven, c).** 36 * 201/(2^30 3^18) = 67/(2^28 3^15) = 1.739467e-14 and
36 * 67/(2^31 3^17) = 67/(2^29 3^15) = 8.697333e-15 exactly; the two faces differ by the factor 2;
36 * 67/36 = 67 and 36 * 67/24 = 201/2. All printed identities hold.

**Proposition H2 (magnitude; c, e).** Read "electron unit" as m_e c^2 (the reading 176.1 uses when it
places m_p in the electron unit). Physical: E1/m_e c^2 = (alpha^2/2) mu/m_e = 2.661118e-5 (CODATA 2022
alpha and m_p/m_e), measured ionization energy / m_e c^2 = 2.661147e-5. Then E1(whole)/physical =
6.54e-10 and E1(cut)/physical = 3.27e-10: the printed binding is low by 9.18 and 9.49 orders of
magnitude. Writing 2 E1 = alpha^k gives k = 6.30 (whole) and 6.44 (cut), against k = 2. E1 is a rational
with only the primes 2, 3, 67; it contains no factor of the model's alpha (124.1). Matching the Bohr
relation E1 = alpha^2 mu/2 would require alpha = 1/5.36e6 (whole) or 1/7.58e6 (cut).
*Failure condition of this proposition:* if "electron unit" denotes a quantity other than m_e c^2, H2
does not apply; the source defines no other unit for this chapter.

**Proposition H3 (the 1/n^2 law and the choice that fixes it; reconstruction b, proven under R).**
Reconstruction R: "the washed symmetric power of n less one passages along the bond", with "the register
of the ladder equal to the dimension of the register of the pair", is read as the SU(m)-twirled state on
Sym^(n-1)(C^m), which by Schur's lemma is I/dim (dim = C(n+m-2, m-1)); the two-vertex term reads it
twice and gives (1/dim)^2. Then:
- m = 2 (the rank-2 common support of a (2,3) pair): dim = n and the ladder is 1/n^2 for all n.
- m = 3: 4/(n^2 (n+1)^2); m = 4: 36/(n^2 (n+1)^2 (n+2)^2); neither is 1/n^2 for n = 1..8 (checked exactly).
The twirl identity was checked by Monte Carlo for n = 1..4 (eigenvalues 1/n to 3 decimals).
Hence, under R, the Bohr exponent is fixed by the choice of a two-dimensional register, not by any
Coulomb dynamics; and the level-n register has dimension n, whereas the hydrogen level n has orbital
degeneracy n^2. *Failure condition:* R is our reconstruction; if the source's ladder is a different
object the proposition is void (the source prints no formula from which to decide).

**Proposition H4 (fine structure, Lamb shift, hyperfine; a).** The page prints no number for any of the
three. The only testable bit is the sign of the Lamb interval ("positive"), which agrees with
2S1/2 - 2P1/2 = +1057.845(9) MHz (Lundeen & Pipkin 1981). The hyperfine reading actually computed is
zero (144.1), against 1420.405751768 MHz; the source grades it OPEN.

## 3. Derived or inserted?

The n-dependence 1/n^2 is derived from a construction whose register dimension is chosen (H3). The
scale E1 is computed from register cells that the page does not re-derive, and disagrees with hydrogen by
about nine orders under the only stated unit (H2). No fine-structure, Lamb or hyperfine number is derived.

## 4. Choices vs digits matched

Choices visible on the page: (1) docking at both ends (the counter-test docks at one end and the mass
face becomes undefined); (2) qutrit end at the electron dock, qubit end at the up core; (3) whole vs cut
closure face (factor 2; both printed, none selected); (4) hub unit 201 vs dock unit 67; (5) the two body
registers; (6) the two-vertex term vs the per-record ladder; (7) the range n <= 8. At least 7 choices.
Matched: one integer (the exponent 2 of n) and one sign (Lamb). Digits of R_infinity matched: 0.

## 5. Lean (task 4)

The section names eight declarations. Seven are `opaque ... : Prop` with no content
(`CompositeMatter.Ch14.row_524296_HPMR`, `CompositeMatter.Ch15.row_519522_HYHC`, `_519828_RYHR`,
`_519937_FSHL`, `_520015_RLRR`, `_520248_RLPQ`, `_520663_RLTV`). The eighth,
`CompositeMatter.Ch15.row_531968_AMBL : cond_531968_AMBL -> D531968_AMBL.checks`, takes an opaque
hypothesis and proves that a one-row `String` table has one row with `key = "C-EQLA"`, `E1 = "YES"`,
`E2 = "YES"`, `letter = "COMPOSED"`. No rational, no ladder and no 1/n^2 appears in Lean. The Lean does
not match the prose claim, nor any arithmetic sub-step of it.

## 6. Failure conditions

- "Hydrogen spectrum" as a physical claim fails at the scale (H2) unless a unit other than m_e c^2 is
  defined and justified.
- The 1/n^2 law would be falsified, as a derivation, if the choice of a rank-2 register is shown not to
  follow from the axiom (H3), or if extended beyond n = 8 the two-vertex ladder departs from 1/n^2.
- Fine structure: testable once the cut gain is given; it must be O(alpha^2) relative to E_n and depend on j.
