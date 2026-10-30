# The Lean at the table-check and opaque-only rows (round 7)

Source: compendium 0.1, commit `a28db899`, directory `compendium/lean`. Object: the 48 of the 64 Lean-tied THEOREMs that
round 6 (`09-lean-audit-at-scale.md`) left unresolved. By `04-lean-coverage.md`'s tiers, 37 cite at best a data-table
check and 11 only opaque Props, aliases or open-result records. Round 6 resolved the other 16 clause by clause. This
round does the same for the 48, then re-states the 64-row table.

Scripts (all in `computation/chain-theory/compendium/`):
- `lean_tables.py` does the resolution. It imports `lean_classify.py`, `lean_resolve.py`, `lean_resolve_data.py` and
  `boxes.py` unchanged.
- `lean_tables_data.py` is the curated layer: 248 clauses, with their counterpart checks and notes. Every clause and
  every check anchor is checked verbatim against the tex and the Lean.
- `lean_audit/TableChecks.lean` holds 7 theorems of ours about the source's entered rows, compiled against round 6's
  Build B.

Recorded outputs (`outputs/`):
- `lean_tables.txt`, sections [1]–[10], holds every number below.
- `lean_tables.tsv` has one row per tied THEOREM (64).
- `lean_tables_clauses.tsv` has one row per clause of the 48 rows.
- `lean_table_checks.txt` is the Lean output of `TableChecks.lean`.

Metrics: `metrics/chain-theory/compendium/lean-table-rows.yaml`.

## 1. The 48 rows and their 131 declarations ([1])

The 48 rows lie in 33 sections and cite 208 names, 131 of them distinct:
- 53 open-result records;
- 38 `opaque X : Prop`;
- 37 table checks: 24 unconditional and 13 under an opaque hypothesis;
- 3 aliases.

The 37 table checks restate the `checks` of 37 distinct data modules, one module each. As in round 6, every statement
of a section inherits the section's names.

## 2. What the 37 table checks prove, and how ([2])

Each tied declaration is the module's `checks`, restated. Two shapes occur, quoted verbatim:

```lean
theorem row_535360_NDFE : type_of% ChainTheory.Data.D535360_NDFE.checks := ChainTheory.Data.D535360_NDFE.checks

opaque cond_530426_UB4H : Prop

theorem row_530426_UB4H : cond_530426_UB4H -> type_of% ChainTheory.Data.D530426_UB4H.checks :=
  fun _ => ChainTheory.Data.D530426_UB4H.checks
```

- **Proof.** `checks` is a conjunction proved by `And.intro` of lemmas `check_i`. All 618 lemmas are `by decide`.
- **Hypotheses.** All 13 are opaque, and every binder is unused (`fun _ =>`). No other theorem of the library mentions
  any of them. A conditional theorem is `H -> checks` with `checks` proved outright, so the hypothesis adds nothing.
- **Data.** All 37 modules enter their `rows` as literal lists.
  - 33 have only `String` fields. 32 of them are key/value tables.
  - 4 have `Int` fields: WAMR, HTOR, SLBC and NCRE.
- **The 618 conjuncts, by form:**

  | form | meaning | conjuncts |
  |---|---|---|
  | LEN | `rows.length = N` | 37 |
  | KV | a key/value pair occurs once, e.g. `(rows.filter (fun r => r.key == "eqtf.mono" && r.value == "MONOTONE")).length = 1` | 560 |
  | CAT | a count or `all` over `String` fields compared with literals, e.g. `(rows.filter (fun r => r.equal == "YES")).length = 16` | 13 |
  | REL | a relation or value test over `Int` fields | 8 |

- **The 8 REL conjuncts** lie in three modules. These are WAMR (tied to 28.1, 28.2, 28.5), SLBC (217.1, 217.2) and NCRE
  (234.1, 234.3). The whole `checks` statements of SLBC and NCRE, verbatim:

  ```lean
  theorem checks :
      (rows.length = 228) /\
      (rows.all (fun r => decide (0 <= r.W_diff_site)) = true) /\
      (rows.all (fun r => decide (0 <= r.R_diff_site)) = true) /\
      (rows.all (fun r => r.W_target_site - r.W_source_site == r.W_diff_site) = true) /\
      ((rows.filter (fun r => r.both_engines == "SAME")).length = 228) :=

  theorem checks :
      (rows.length = 246) /\
      (rows.all (fun r => r.order0 + r.order1 + r.order2 + r.order3 == r.records) = true) /\
      (rows.all (fun r => r.same_site + r.one_link + r.cross_link == r.records) = true) /\
      ((rows.filter (fun r => r.group == "KIT" && r.polarity == "SLOT" && r.defined == "YES")).length = 64) /\
      ((rows.filter (fun r => r.both_engines == "SAME")).length = 246) :=
  ```

- **Three modules have exactly 60 entered rows** (CASM, DST2, ANDB). Two of them stop inside a record:
  - CASM's last detail key is `cas.SINK.ARC.COUNT.QUTRIT.s0-s1.N4.ratios`.
  - ANDB's last key is `web.QUTRIT.s0-s1.RING.4.tiled.spectrum.eq`.
- **Axioms** ([8], from round 6's file). WAMR, SLBC and NCRE give `[propext, Quot.sound, Classical.choice]` at 4.35.0.
  The other 34 table checks, all 94 opaque, alias and record names, and every name of the 48 rows give no axiom.
  `sorryAx` occurs 0 times.

## 3. The opaque Props, aliases and records: what the Lean asserts ([3])

Verbatim samples:
- `opaque wall_517119_RULG : Prop` (I.16.2).
- `def wall_531133_EQDU_open : ChainTheory.OpenResult := { address := 531133, stem := "DIMC", decidingClass := "EQDU", firstParting := "-", chapter := "ch03", statement := wall_531133_EQDU }`
  (28.x).
- `abbrev row_528690_VOFW : Prop := cond_528690_VOFW -> claim_528690_VOFW` (200.6).

What the 94 names are:
- 38 are opaque constants of type `Prop`. They have no definition and no proof.
- 53 are records whose `statement` field is such a constant.
- 3 are aliases. Two rename an opaque constant. One is the implication between two opaque constants.
- None of the 94 is a theorem. No library theorem (outside the source's `AxiomsCheck`/`Certificate`/`Roots` scripts)
  mentions any of them.

**For the 11 opaque-only rows, the Lean asserts nothing.** There is no proved proposition, opaque hypothesis or not.
Two of these boxes say "the deciding table is re-derived in Lean": 170.1 (two records tied) and 229.1 (one record
tied). Some of their printed numbers occur nowhere in the tree: 17689 (170.1) and 46412, 46563 (232.2).

## 4. Do the declarations state the boxes? ([4], [5])

**Rule.** The rule is round 6's (09 §3 and §8), refined below for REL-form checks over entered tables; §5 regrades
round 6's 16 rows under the refinement (corrected after review 2026-10-30, R4-CT-8; see `reports/CORRECTIONS.md`).
- Each box is split into clauses, quoted verbatim. Here that gives 248 quotes, all found.
- A clause has a counterpart when a tied declaration states its finite content under the library's own definitions of
  its words. As applied here and in round 6, a clause counts when *some part* of its finite content is stated; the
  strict reading, in which all of it must be stated, is given under Sensitivity (corrected after review 2026-10-30, R4-CT-9; see `reports/CORRECTIONS.md`).
- MATCH means every clause has a counterpart, PARTIAL some, NONE none.

Applied to tables, two kinds of check are separated:
- **ENTERED (a counterpart).** The check evaluates the relation the clause asserts over entered per-instance values of
  the quantities the clause is about. This is round 6's ENTERED kind (`LemmaB_zero_transfer` evaluates the entered
  step triples), restricted to relations: an entered value that is itself the clause's conclusion is RECORDED, and §5
  applies this to 9.4 (corrected after review 2026-10-30, R4-CT-8; see `reports/CORRECTIONS.md`).
- **RECORDED (not a counterpart).** The check confirms that the entered list holds a string that records the clause's
  conclusion. That string is a verdict (`YES`, `MONOTONE`, `NONE`), a count written as a string (`"5 of 5"`, `"13"`),
  or a per-instance verdict flag. The library defines none of the words of the key. The string is the claim.
  - A RECORDED match requires a lexical link: the key or value carries a number or content word of the clause.
  - Whole-box verdict cells (`cell.C-EQ.. = YES`) are not matched to clauses. 33 of the 48 rows have one.
  - RECORDED matches give a label-inclusive upper bound.

**Results.**
- Verdicts: PARTIAL 2 (217.1, 234.1) and NONE 46 (35 table rows, 11 opaque rows). No row is MATCH.
- Clauses: 5 stated, 44 recorded only, 199 with nothing.
- Label-inclusive: MATCH 1 (28.2), PARTIAL 22, NONE 25.

| row | verdict | label-incl. | what the tied checks state or record |
|---|---|---|---|
| 217.1 T | PARTIAL 4/8 | PARTIAL 4/8 | SLBC: at all 228 entered reads `W_target_site - W_source_site = W_diff_site >= 0` (non-decreasing, no decreasing act). Missing: the opening definition, the running count along the order coordinate, the corollary, the arrow |
| 234.1 T | PARTIAL 1/7 | PARTIAL 2/7 | NCRE: orders and kinds sum to the records at every entered network row (246). Recorded: 64 KIT SLOT rows flagged `defined = YES`. Missing: the count law, the lowering depths, termination at order three |
| 28.2 T | NONE | MATCH 1/1 | MWSP labels `eqsp.built = "5 of 5"`, `"10 of 10"`, `BUILT` |
| 223.2 T | NONE | PARTIAL 5/7 | ENT2 labels `STRICTLY-INCREASING`, `"182940"`, `TWO-OBJECTS`, `fluc.function.name = NONE`, `onsg.setequal = EQUAL`. **The table enters and checks `onsg.asym.pairs = "24"`; the box prints "every one of the 43"** (the key is not defined) |
| 220.2 T | NONE | PARTIAL 4/8 | TOP2 labels `MONOTONE`, `STRICTLY-INCREASING`, `"5 of 12"`, `FIXED`. 182940, 185496 and 2556 are not entered |
| 11.1, 65.2 | NONE | PARTIAL 3/· | cells named `UNIT`, `POLE`, `NORM` / `ALGB`, `HUBE`, `DATM`, `HIGG` = `YES`. No dimension, rank or order is entered (65.2's 15, 3, 4) |
| 221.1, 225.1 | NONE | PARTIAL 3/9 | `"9377"`, `CHARGED>DARK`, `STRICTLY-INCREASING` / per-ring `YES` labels |
| 145.1, 219.1, 220.3, 224.1, 230.1, 235.3, D.3 | NONE | PARTIAL 2/· | count and verdict labels (`"16"`, `"13"`, `"1692 of 1692"`, `"36"`, `"16 of 16"`, ...) |
| 3.1, 13.1, 28.1, 224.3, 226.1, 234.3, 235.4 | NONE | PARTIAL 1/· | one label each (`record.moves.covariance`, `FINITE`, `wash_moved_acts == 0`, ...) |
| 28.5, 57.1, 65.1, 181.5, 186.1, 200.6, 203.5, 217.2, 224.2, 235.1, 235.2, D.4, H.3, I.7.1 | NONE | NONE | nothing of the box; see below |
| 62.1, 64.1, 64.2, 73.1, 170.1, 229.1, 232.2, 233.1, 233.2, 233.4, I.16.2 | NONE | NONE | opaque/records only (§3) |

**Sensitivity.** Two choices change verdicts. [4] and [10] print both.
- **Entered conclusions.** Two REL-form checks are graded RECORDED, because the entered value is the clause's
  conclusion:
  - WAMR's `wash_moved_acts == 0` (28.1: "unmoved under every change confined to a washed end");
  - SLBC's `0 <= R_diff_site`, an entered difference.

  Counting them ENTERED, as round 6's rule would (round 6 did not separate an entered conclusion from a relation over
  entered values; §5, 9.4), adds one PARTIAL (28.1, 1/5). The other verdicts are unchanged (corrected after review 2026-10-30, R4-CT-8; see `reports/CORRECTIONS.md`).
- **Partly stated clauses** (corrected after review 2026-10-30, R4-CT-9; see `reports/CORRECTIONS.md`). Three credited clauses of the 48 are stated only in part:
  - 217.1, clause 3 ("at all 228 one-act reads of the fourteen bodies (two per link at 114 links, 78 with an open far
    end, 150 washed)"): only `rows.length = 228` is checked, a LEN conjunct; the 14 bodies, 114 links, 78 open and 150
    washed are not. This is the only LEN conjunct credited in this round (230.1's `rows.length = 14` is counted a
    coincidence).
  - 217.1, clause 4 ("rises or stays at the site bound, at the free-direction bound and at the directions read"): only
    the site bound (W) is checked. The free-direction fields are not checked, and the directions read
    (`0 <= R_diff_site`) are RECORDED.
  - 234.1, clause 2 ("summing to the records at every link and network"): NCRE checks the sums at its 246
    network-and-polarity rows (82 networks × 3 polarities). No per-link row is entered.

  Under the strict reading 217.1 is PARTIAL 2/8 (clauses 2 and 5, both the W non-decrease) and 234.1 is NONE. Of the 48
  rows, 1 then states a clause, and 2 clauses are stated rather than 5. The 64-row table under each reading is in §6.

**Printed numbers** (numbers of 5 or more, literal lookup as in round 6, so coincidences count). There are 135
(row, number) pairs:
- 23 occur in a tied `checks`, 11 only in entered rows, 1 only in another tied declaration, and 100 not in the tied
  declarations.
- 15 occur nowhere in the raw text of the 915 `.lean` files:
  - 170.1: 17689;
  - 181.5: 101761 (= 319²);
  - 220.2: 185496;
  - 220.3: 322516, 138129, 570349;
  - 221.1: 101973, 44341, 218340, 93481, 93312, 9270;
  - 232.2: 46412, 46563;
  - 234.3: 270584.
- Some literal hits are coincidences: 230.1's 14 is `rows.length = 14`, and 217.1's 13 and 14 are row values.

**Notable rows.**
1. **3.1.** The box ends: "drafted as Lean text with its Lean proof still to be written". Its own module, NDFE, is named
   in the box label `box:I-1-D535360-ndfe-t1:1` and holds seven labels. The sentence on the reads Z, X, XZ, XZ² is
   entered in HTOR, which is tied to 57.1 (section 57), not to 3.1:
   - the open read Z fixes Z, with 1 fixed and 2 unfixed;
   - the open reads X, XZ, XZ² fix none, with 3 unfixed;
   - `Round7.htor_open_reads` checks this.

   57.1 itself gets nothing from HTOR or PQCD.
2. **181.5 and 200.6: the boxes' own modules are not tied.**
   - The labels name D538212-qarm and D538733-hcos. Both modules exist, and `Frame.P07.Ch7.row_538212_QARM` and
     `Frame.P06.Ch10.row_538733_HCOS` restate their checks.
   - The index ties neither. It ties QLHL (quark-ladder hub counts 9 and 124) and HGHL (`28 of 28` keyed legs).
   - The untied modules are one-row labels, `ARMS-AND-SINKS-BY-THE-RULE` and `OFFICES-AND-RATIOS-BY-THE-RULE`. Graded,
     they too would give NONE.
3. **145.1 and 225.1: the tables stop at 60 entries** ([5]).
   - CASM's detail covers one chain (`QUTRIT.s0-s1`) at six separations, N = 2, 3, 4, 10, 11, 12, in lexical key
     order. The box says "all forty-four chains".
   - ANDB's per-ring keys cover 18 rings, the last without its translation entry. The box says "20 networks built, 20
     periodic, 20 tiled".
4. **217.1: the checked count does not move** ([5], `Round7.slbc_*`).
   - At all 228 entered reads, `W_diff_site = W_diff_free = 0`. The non-decrease that the checks verify holds with
     equality at every read.
   - The R count rises at 54 reads. `R_target_site - W_source_site = R_diff_site` holds at all 228 but is not checked;
     only `0 <= R_diff_site` is.
   - The printed 114 links, 78 open, 150 washed and 14 bodies follow from the entered rows but are not checked. WAMR,
     tied to section 28, sums to the same 228, 114, 78 and 150 (`Round7.wamr_sums`). 28.1's "two hundred
     twenty-eight one-act reads" is that sum, which WAMR's checks do not state.
5. **234.1.** The 18 of "18 of 18 defined" follows from the entered rows: 13 BODY, 3 NEUTRAL and 2 STAR outside KIT.
   `order3 = 0` holds in all 246 rows. Neither is checked (`Round7.ncre_rows`).
6. **H.3.** The box is a finiteness argument ("the moves map the webs of at most $S$ sites into themselves and such webs
   are finitely many"). ERMT records 14 demonstration runs: 7 `STAYS-BOUNDED` and 7
   `GROWS-TO-THE-DEMONSTRATION-CAP`, with max sites 40 to 42.
7. **235.4.** The box: "absent at the gluon". The table's checked `NO-CYCLE` entries are at `r-HIG`. Its readings
   COMPOSED and PHASE-CONTENT are checked at `0 of 2` and `0 of 3`, which the box does not print.
8. **203.5** is a pointer: "stated once, as Theorem~170.1". 170.1 is an opaque-only row. **62.1** says that "the
   theorem rests on the weakest of its parents, an open result among them".

## 5. The 16 rows of round 6, over all their tied declarations ([6])

Round 6 graded the 16 rows against their restated/kernel declarations only. They also cite 9 table checks: TTCF,
DIM1, LTGO, RTSE, GPOE, NRGE, NVCI, AJBP and TLWN.
- None of the 9 has a REL conjunct (kinds: LEN 9, KV 22, CAT 7). Under the rule they add no counterpart, so round 6's
  verdicts (MATCH 1, PARTIAL 5, NONE 10) hold over all tied declarations. They are recomputed here from
  `lean_resolve_data.py` and equal `lean_resolve.tsv` at 16 of 16.
- **9.4 under the refinement of §4** (corrected after review 2026-10-30, R4-CT-8; see `reports/CORRECTIONS.md`). Round 6's only MATCH has an entered counterpart for its
  third clause, "zero register displacement equals zero momentum transfer at the one-link star"
  (`01-framework.tex:360`). It is credited to `LemmaB_objectOfRecord_electromagneticLike` through its field
  `zeroTransfer` = `steps.all fun s => s.2.2.2 == (0,0,0) && (first + second) == (0,0,0)`, over the 172 entered step
  triples of `AlphaDerivation/Data/Records.lean` (09 §2 item 4). That check has two parts:
  - (A) the third triple is $(0,0,0)$ in all 172. This is the conclusion "zero transfer", entered, so it is RECORDED
    under §4's refinement;
  - (B) the second triple is minus the first in all 172 (first triples $82\times(0,0,0)$, $48\times(1,-1,0)$,
    $42\times(-1,1,0)$). This is a relation over entered values, so ENTERED. It reads as "zero register displacement"
    only by the position of the triples, which carry no field names; the library's own
    `displacement (chain objectOfRecord)`, which is $(-1,1,0)$ (`LemmaB.massless_faces`), is not used.

  9.4 therefore stays MATCH 3/3 if (B) is taken as the clause's displacement half (the practice of §4), and is
  PARTIAL 2/3 under the strict reading. Round 6's treatment counts (A) as well and also gives MATCH. The other 15
  verdicts do not depend on the refinement. Under the strict reading 196.2's clause 2 is also stated only in part (12
  open sites are computed; the 17 sites lie outside `checks`), which leaves 196.2 PARTIAL.
- Label-inclusive, five NONE rows become PARTIAL:
  - 2.6: TTCF has 32 rows, each `eq = YES`.
  - 2.9: DIM1 has `C-EQD1 ... EXCLUDED`.
  - 74.1: NVCI has `C-KERN` and `C-CHRG`.
  - 74.4: NRGE has `"577 of 577"` junctions.
  - 78.3: AJBP has `eqab.order3.zero = YES` and `eqab.order2.identity = YES`, the same pattern as 141.1's `ZERO` label
    quoted in 04 §2.

  The label-inclusive bound for the 16 is MATCH 1, PARTIAL 10, NONE 5.
- Scope note: 09 §3 says that for 2.6 "32 does not occur in the tied statements". The statements there are ORKC's,
  the restated declaration. TTCF, also tied to 2.6, has `rows.length = 32`. The verdict is unchanged.

## 6. The 64-row table re-stated ([10])

| tier by kind (04) | by content | rows |
|---|---|---|
| restated/kernel (16) | states every clause | 1 (9.4) |
| | states some clause | 5 (2.1, 9.2, 9.3, 18.2, 196.2) |
| | states no clause; the computation is another box's | 10 |
| table checks (37) | states some clause (ENTERED relations over an entered table) | 2 (217.1, 234.1) |
| | states no clause; labels record some clause | 21 |
| | states no clause; nothing of the box | 14 |
| opaque/records only (11) | states no clause; the Lean asserts nothing | 11 |

- **All 64:** MATCH 1, PARTIAL 7, NONE 56. **8 of 64 rows** have a tied declaration that states at least one clause of
  the box.
- **Sensitivity** ([10]) (corrected after review 2026-10-30, R4-CT-8; see `reports/CORRECTIONS.md`). The headline depends on two choices of §4 and §5:

  | reading | 9.4 | 28.1 | 217.1 | 234.1 | MATCH / PARTIAL / NONE | rows stating a clause |
  |---|---|---|---|---|---|---|
  | round 6's: an entered conclusion counts | MATCH | PARTIAL 1/5 | PARTIAL 4/8 | PARTIAL 1/7 | 1 / 8 / 55 | 9 |
  | headline: round 7's ENTERED/RECORDED split; a clause counts when any part of its finite content is stated | MATCH, through (B) | NONE | PARTIAL 4/8 | PARTIAL 1/7 | 1 / 7 / 56 | 8 |
  | strict: round 7's split; a clause counts only when all its finite content is stated | PARTIAL 2/3 | NONE | PARTIAL 2/8 | NONE | 0 / 7 / 57 | 7 |

  The other 60 rows have the same verdict under all three. The headline uses the middle reading. It applies this
  round's ENTERED/RECORDED split to all 64 rows, and it is the clause-counting practice recorded clause by clause in
  `lean_tables_data.py` and `lean_resolve_data.py`. The strict reading needs a further judgment, of where each clause's
  finite content ends; [10] lists the five credited clauses it removes. All three readings agree that at most one row
  (9.4) states its whole box, that none does under the strict reading, and that 7 to 9 of the 64 state some clause.
- **Label-inclusive upper bound:** MATCH 2, PARTIAL 32, NONE 30.
- Of the 8, 6 rest on round 6's restated/kernel declarations, partly over entered inputs: 9.4's zero transfer is an
  entered step table (§5), and 18.2's weights are a definitional identity applied to the entered parity and kind
  fields of the 172 records (09 §2 item 5). The other 2 rest on relations checked over entered `Int` tables, whose
  values are not derived in Lean (corrected after review 2026-10-30, R4-CT-10; see `reports/CORRECTIONS.md`).

What changes:
- **04.**
  - The tiers 16 / 37 / 11 by declaration kind are unchanged. So is "at best data-table checks": all 37 modules enter
    their rows, and all their proofs are `decide`.
  - New: by content, 2 of the 37 table rows state a clause, 35 none, and the 11 opaque rows assert nothing.
  - 7 of the 37 table rows (28.1, 28.2, 28.5, 217.1, 217.2, 234.1, 234.3) cite a table with relation checks over `Int`
    fields (WAMR, SLBC, NCRE). The register's wording "37 at best string-label table checks" (CT-C003) holds for the
    other 30. See `reports/CORRECTIONS.md`.
- **09.**
  - §7's "best evidence: data-table checks (not resolved here) 37 (+9 = 46)" is now resolved: 2 PARTIAL, 44 NONE.
  - "only opaque … 11 (+1 = 12)" becomes 12 NONE.
  - The README's "coverage by content: 6 / 46 / 12" becomes 8 rows stating a clause and 56 stating none (10 / 35 /
    11).
  - §8's limit "the 37 table-check rows and the 11 opaque rows were not content-checked" is lifted.
- **Round 6's number lookup.** `outputs/lean_resolve.txt` [4] lists a printed 5 for 2.3. The box prints none:
  `numbers_in` reads "four one-pair acts" as 4 + 1. The same parser reads "two hundred twenty-eight one-act" as 229,
  "twenty-fourth" as 20, and "four, thirty-two, three hundred eighty-four and six thousand one hundred forty-four" as
  one number.
  - `lean_tables.py` ends a number phrase at punctuation and at a unit word that follows a unit word, and reads ordinal
    words. [4] lists the 7 boxes where the two parsers differ.
  - No round-6 verdict or count depends on this (2.3 is NONE, with no tied number). Recorded in
    `reports/CORRECTIONS.md`.

## 7. Own checks ([9], `lean_audit/TableChecks.lean`)

These are 7 theorems of ours, each `by decide` over the source's entered rows. All compile against Build B in under
2 s.
- `wamr_sums`: 228 / 114 / 78 / 150.
- `slbc_counts`: OPEN 78, WASHED 150, AB 114, BA 114.
- `slbc_W_constant`.
- `slbc_R`: R consistency at all 228; 54 rises.
- `ncre_rows`: order3 = 0; 18 = 13 + 3 + 2.
- `htor_open_reads`.
- `rlro_webs`.

Three use no axiom. Four use the three core axioms. All four test `Int` fields inside `Bool` expressions, and round 6
found `Classical.choice` entering at 4.35.0 through core `Int` order instances. The route for these four was not traced.

## 8. Limits

- **Judgments.** Clause boundaries and the RECORDED lexical matches are judgments. They are recorded clause by clause
  in `lean_tables_data.py`, and their quotes and anchors are checked. The ENTERED/RECORDED line for REL-form checks and
  the counting of partly stated clauses are stated, with their sensitivities (§4, §6).
- **Keys.** The keys of the key/value tables are not interpreted. A RECORDED match says only that a key's words or
  numbers are the clause's. The ENTERED grade does interpret field names: 217.1's counterparts read the entered columns
  `W_source_site`, `W_target_site` and `W_diff_site` as the box's boundary count (the entanglement entropy across a
  cut), a reading that rests on the column names and that no library definition supports. In the ENTERED case a
  relation is computed over the entered values; in the RECORDED case a string is compared (corrected after review 2026-10-30, R4-CT-11; see `reports/CORRECTIONS.md`).
- **Scope.** Only declarations tied by the index are graded. The untied modules of 181.5 and 200.6 are reported, not
  graded. The 9 table checks of the 16 rows are graded only as to whether they could add a counterpart.
- **Toolchain.** Axiom sets are those of round 6 (Lean 4.35.0, Build B). The pinned 4.32.1 was not used.
