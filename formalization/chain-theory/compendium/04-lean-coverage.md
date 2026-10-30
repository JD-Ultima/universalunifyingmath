# What the compendium's Lean checks, for the audited claims (task 4)

Source: `compendium/lean/` (six core-Lean 4.32.1 packages, 915 `.lean` files). Read as text; not built
(memory). The build record and `axioms_report.txt` were therefore not re-verified. Script:
`computation/chain-theory/compendium/lean_classify.py` (output `outputs/lean_classify.txt`).

## 1. What "64 theorems tied to Lean" means

- The `lean` column of `statements.tsv` is filled per section from `chapters.tsv`; every statement of a
  section inherits the section's list (0 of 284 sections differ internally). "Tied" means the section
  names at least one declaration that exists.
- Classification of every declaration named by the 64 Lean-tied THEOREM rows (195 distinct names):
  69 `opaque X : Prop` (uninterpreted, no content), 69 `OpenResult` records, 41 restatements of a data
  module's `checks` (25 + 16 under an opaque hypothesis), 5 other table checks, 7 restatements of
  `AlphaDerivation` statements, 3 `abbrev` aliases of Props, 1 kernel-module check.
- By strongest evidence per row: 11 of 64 cite only opaque Props, aliases or open-result records;
  37 have at best data-table checks; 16 cite a restated `AlphaDerivation` theorem or a kernel module.
- 137 of the 150 data/kernel modules that define a `Row` structure have only `String` fields: their
  `checks` are `decide` proofs that given string labels ("YES", "EQUAL", "ZERO", "8 of 8") occur in a
  hand-entered list a given number of times.

This is consistent with, and sharper than, the source README's own caveat ("A chapter's named theorem
therefore usually certifies the arithmetic of that chapter's tables. The physical reading of the result
lives in the text, not in the Lean statement."): for the three audited claims, the Lean does not
certify the arithmetic either.

## 2. The three audited claims

| Claim | Declarations named | What they prove | Matches the prose? |
|---|---|---|---|
| 141.1 (a_e) | `Fields.Ch09.row_523203_AMSP` | `opaque : Prop` | no content |
| | `Fields.Ch09.row_526870_PAET` | 8-row `String` table: one row (`ORDER-TWO-READOUT`, `W1`, `YES`, `EQUAL`), one `ZERO`, three `PARTING`, one `FACE` with `8 of 8` | no; T1 = 309 is not in Lean |
| | `Fields.Ch09.wall_526122_EQAJ_open` | `OpenResult` record | n/a |
| 143.1 (hydrogen) | seven `CompositeMatter.Ch14/Ch15.row_*` | `opaque : Prop` | no content |
| | `CompositeMatter.Ch15.row_531968_AMBL` | under an opaque hypothesis, a 1-row `String` table has key `C-EQLA`, `E1 = YES`, `E2 = YES`, `letter = COMPOSED` | no; neither E1 nor 1/n^2 is in Lean (the denominators 2^28 3^15, 2^29 3^15 occur nowhere in the tree) |
| 64.1, 64.2 (gauge) | `Foundations.Ch02.row_521259_GIKC`, `StandardModel.Ch10.row_510060_GGRL` | `opaque : Prop` | no content |
| 65.1 (gauge) | `Foundations.Ch05.row_530016_GSTB` | 8-row key/value `String` table (`cell.C-EQGS = YES`, `eqgs.letter = EQUAL`, ...) | no; no group in Lean |
| | `Foundations.Ch05.row_530426_UB4H` | same shape, under an opaque hypothesis | no |

Verbatim, the strongest of these (from `ChainTheory/Data/D526870_PAET.lean`):

```lean
theorem check_3 : (rows.filter (fun r => r.test == "ORDER-THREE-TERM" && r.letter == "ZERO")).length = 1 := by decide
```

i.e. "the order-three term is zero" is a string entered in the table, and the theorem counts it.

## 3. Formalizability of the audited content

- G1-G3 of `03-gauge-group.md` (S(U(2) x U(3)) = (SU(2) x SU(3) x U(1))/Z6; dimensions and centre
  action for all six genres; uniqueness of (3,2) under cap 6 and d1 != d2) are finite or standard Lie
  theory and are formalizable with Mathlib in a few hundred lines.
- A1 (T1 = 309) and H1 (the E1 rationals) are one-line `norm_num` facts; they would certify arithmetic
  only, not the physical readings.
- Propositions D (a_e) and H2 (hydrogen scale) are comparisons with measurement and are not Lean targets.
- The constructions themselves (anomaly web, held carrier, "cut gain", "dark cut", "closure") are not
  defined precisely enough in the text to formalize.
