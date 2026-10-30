# CT-C015: "the electron anomalous moment through third order" (statement 141.1)

Source: compendium 0.1, commit `a28db899`. Abstract (`tex/front.tex:25`) lists "the electron anomalous
moment through third order" among the conditional theorems. The boxed statement is 141.1,
`tex/parts/07-cell-3-2-photon.tex:426-428`; derivation lines 438-444; data table
`lean/ChainTheory/ChainTheory/Data/D526870_PAET.lean`. Computations: `computation/chain-theory/compendium/ae_audit.py`
(output `outputs/ae_audit.txt`).

## 1. What the source asserts (a)

Box 141.1, verbatim core: "At the anomaly web with the measuring body in place the readout at order two
equals the sum of the tabulated self-terms at all eight own-link records, and the order-three term is
zero, the direct record sum graded by the surviving exchange count terminating at order three at all
eighteen joined webs." Grade: CONDITIONAL. Direct premises (statements.tsv): Definitions 1.1-1.6;
18.1, 18.3, 18.4, 124.1 (CONDITIONAL); 18.2, 28.1, 28.2, 28.5 (THEOREM); 18.5, 18.6, 28.3, 28.4 (OPEN).

Numbers printed in the derivation (line 441): T0 = 1, T1 = 3(67+30+6) = 309, T2 = 309^2 + closure
(closure "an exact element of the field of the two roots", value not printed), T3 = 0,
G2(2,2) = (5/216)(1+sqrt2). The derivation states: "no digit of any constant displayed here".
The QED coefficient 1/2 appears only in the literature paragraph (`chbridge`, line 435, citing
Schwinger 1948); no expression in alpha/pi appears in the chapter.

## 2. Reconstructed definitions (b)

- Link measure mu(d1,d2) = d1 d2. (Checked: the printed measures 9, 12, 8, 4, 6 equal d1 d2 for
  (3,3), (3,4), (2,4), (2,2), (2,3).)
- Anomaly web W: a two-dock body (8 unstruck links, measure sum 67, dressing 3 at its one hub) joined
  to the photon word (6 sites, 5 links of type (2,3)) by one made (2,3) link; 15 sites, 16 links,
  2 dock links struck.
- Record sum R(W) = sum_n T_n, graded by "surviving exchange count" n; T1 = dressing x sum of measures of
  unstruck links.
- The chapter does not give the body's link list; 13 multisets of eight measures from {4,6,8,9,12}
  sum to 67 (computed), so T1 = 309 cannot be re-derived from the page alone, only re-added.

## 3. Results

**Proposition A (arithmetic; proven, c).** With the printed inputs, T1 = 3(67+30+6) = 309, T1^2 = 95481,
14 unstruck links, and 105 = C(15,2) = number of degree-2 monomials in 14 variables (also the number of
pairs of 15 sites; the page does not say which). All verified exactly.

**Proposition B (T3 = 0 is structural, not a computation of the anomaly; reconstruction b).**
Statement 146.1 (CONDITIONAL; `07-cell-3-2-photon.tex:619`) asserts that at every one of the sixteen
webs "the record sum graded by order carries exactly the three rungs zero, one and two". Granting 146.1,
T3 = 0 holds at every web; the order-three clause of 141.1 is an instance of that truncation and carries
no information specific to the electron.

**Proposition C (no QED coefficient is produced; c).** The chapter, its data table and its Lean
declarations contain no coefficient of (alpha/pi)^n. C1 = 1/2, C2 = 197/144 + pi^2/12 - (pi^2/2) ln 2 +
(3/4) zeta(3) = -0.328478965579194 and C3 = 1.1812414565872 (Laporta-Remiddi closed form, evaluated to
40 digits and asserted against the reference values) are neither derived nor inserted; they are cited
in the literature paragraph only.

**Proposition D (literal identification is refuted; c, e).** Let hypothesis L be: there is lambda > 0
with T_n lambda^n = C_n for n = 1, 2, 3 (the record-sum order identified with the power of alpha/pi up
to one free normalisation). Then L is false: T3 = 0 while C3 != 0, for every lambda. Fixing lambda by
n = 1 gives lambda = 1/618 and would require closure = -220935.00045 at n = 2 (untestable: not printed).
Empirical weight: C3 (alpha/pi)^3 = 1.48042e-8, i.e. 1.28e-5 of a_e; against
a_e = 1.15965218059(13)e-3 (Fan et al. 2023) combined with alpha(Rb 2020), dropping it moves a_e by
9.26e4 sigma (5.63e4 sigma with alpha(Cs 2018)).

**Proposition E (charitable reading; b).** Under the reading the source itself states ("the
correspondence is by name"), 141.1 asserts an internal equality: at the adopted reading W1, the
order-two readout of W equals the sum of eight tabulated self-terms. The data table shows this equality
at reading W1 only; readings W2, W3 and PWONCE "part" (0 of 8). The exact W1 self-term at link L1 is
11989/6879707136 + (1675/1719926784) sqrt3 = 3.4295e-6 (denominators 2^20 3^8 and 2^18 3^8); the source
makes no identification of it with any physical quantity. Under this reading the statement makes no
contact with the measured a_e.

## 4. Derived or inserted?

Neither. No anomalous-moment coefficient is derived; none is inserted. The word "third order" refers to
the record-sum grading, which by 146.1 terminates after order two.

## 5. Choices vs digits matched

Discrete modelling choices visible on the page: (1) the thirteen bodies and the depth rule; (2) the
per-web placement with the measuring body in place; (3) which two dock links are struck; (4) the type of
the made link, (2,3); (5) the dressing 3; (6) the reading W1 among W1/W2/W3/PWONCE; (7) the eighteen
joined webs. At least 7 choices; digits of C1, C2, C3 matched: 0.

## 6. Lean (task 4)

The section names three declarations: `Fields.Ch09.row_523203_AMSP` is `opaque ... : Prop` (no content);
`Fields.Ch09.wall_526122_EQAJ_open` is an `OpenResult` record; `Fields.Ch09.row_526870_PAET` restates
`D526870_PAET.checks`, which proves by `decide` that an 8-row table of `String` fields has 8 rows, one
row with `test = "ORDER-TWO-READOUT"`, `reading = "W1"`, `of_record = "YES"`, `letter = "EQUAL"`, one row
with `letter = "ZERO"`, three with `letter = "PARTING"`, and one `FACE` row with `count = "8 of 8"`.
The labels "EQUAL" and "ZERO" are entered as data. The Lean statement does not match the prose claim
and does not contain even the arithmetic sub-step T1 = 309.

## 7. Failure conditions and what would upgrade the claim

- The claim as an a_e result fails if T3 = 0 is read as the (alpha/pi)^3 coefficient (it is, by 9e4 sigma).
- Upgrade path: (i) define the map from record-sum order to powers of alpha; (ii) print the order-1 term
  in that normalisation and compare with 1/2; (iii) print the closure term and compare with C2;
  (iv) explain how a truncation at order two is compatible with C3, C4, C5 != 0.
