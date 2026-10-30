# Statement-level premises (round 4)

Source: compendium 0.1, commit `a28db899`. Script: `computation/chain-theory/compendium/premises.py`. Recorded outputs:
`outputs/premises.txt` (every number below, by section [0]-[8]) and `outputs/premises_items.tsv` (one row per premise
item, 1069 rows, with its quote, location, candidates and judgment). Metrics:
`metrics/chain-theory/compendium/statement-premises.yaml`. Object: proposed further work item 1 of
`reports/chain-theory-compendium.md`, i.e. whether F3 (chain-rule violations, the 124.1 cycles) and F4 (load-bearing
unproven premises) survive when premises are taken statement by statement instead of from the index.

Grades: T THEOREM, C CONDITIONAL, O OPEN, D DEFINITION, Tab TABLE, as printed. "Index" means the `premises` column of
`index/statements.tsv` read as in `depgraph.py` (choice C3; 8367 edges, one strongly connected component of 22,
reproduced in [2]).

## 1. What the source says

- **The chain rule** (`tex/front.tex:34`, "The status words"): "A theorem is a statement derived from the axiom, the
  definitions and the results proved ahead of it; ... a theorem conditional on a premise prints the premise beside it,
  and under the chain rule of this compendium every theorem prints at its statement the union of the conditions of the
  results it rests on, a composed claim carrying the lowest status in its chain." The same paragraph: "A
  presentation, such as a table, a tool, or this compendium itself, makes no claim of its own."
- **The index premise line is chapter-level by construction** (`tex/parts/01-framework.tex:109`): "one line names the
  chapters it rests on and the chapters that rest on it, generated mechanically from the ancestry of the chapters."
  So the F3 caveat ("premise lists cite whole chapters") is the source's own description of that line.
- **Where statement-level premises are.** Each derivation (`chproof`) that names premises does so in one opening
  sentence: "The premises are ...", "The premises of this derivation are ..." or "\textit{The premises.} ...". No
  derivation cites a statement by number: 0 citations of the form "Theorem 9.1" / "Open result 2.4" and 0
  `\ref{box:...}` inside the 297 derivation environments ([1]). Premises are cited as a chapter plus a description and
  a grade word, e.g. "the stabiliser of the pair ket (chapter 5-1, an open result)" or "the theorem of
  chapter~\ref{mod:I-3-...}".

## 2. Extraction rule

Layer **M** (mechanical, all derivations; implemented in `premises.py`, rules M1-M5 in its docstring):

1. *Segmentation.* Sections and boxes numbered as in the index (0 mismatches against `statements.tsv` in number or
   grade). A box's derivation is the next `chproof` in its section; boxes between two derivations share it
   ("derivation group"). 294 groups; 61 hold more than one box. The 9 boxes of section 1 have no derivation.
2. *Premise sentence.* The paren-balanced sentence after the marker. 261 groups have one. 33 groups (37 boxes: 31 Tab,
   3 C, 1 O, 2 T, namely I.1.2 and I.16.2) have none.
3. *Items.* Every chapter reference in the sentence is an item: `chapter~\ref{mod:L}` (one section, printed as its
   number), a designator "chapter P-c", "this/the present chapter/the ... chapter above", "the framework", or a
   premise "outside the chapters"/"placed at no chapter" of the compendium. The grade word comes from the item's own
   parenthesis, "the theorem/open result of chapter", a following ", a theorem", or a plural grade closing a run of
   references. 1069 items; each quote occurs verbatim in the source (1069 of 1069).
4. *Candidates.* The boxes of the cited sections whose grade fits the grade word. "Theorem" admits T and C, because
   the source uses it for conditional theorems: 9.1 (C) is "a theorem" at `01-framework.tex:582` and "a conditional
   theorem" at `:691`. "Open result" admits O, "table"/"instrument" Tab, and "the text"/"definitions" D (Definitions
   1.1-1.6 when the description names the axiom or definitions). An ungraded item admits every box. Items with no
   candidate: 83. Items with one: 236. Items with several: 750; for these the text says only that the box rests on
   one of them ("any-of").
5. *Designators.* Two readings. **MOD**: the module numbering, i.e. the sections whose label begins
   `<roman P>-c-`. **FIF**: the `first_section..last_section` range of `designators.tsv`, which is the index's own
   mapping. The two differ at 171 of 397 designator items. At those 171, the description shares more title words with
   the MOD candidates at 129 and with the FIF candidates at 3 (39 ties). MOD is used throughout, and FIF is reported
   as a sensitivity run. In one checked case the index's own explicit list also follows MOD: 16.1's index list has
   "Theorem 124.1", "Theorems 185.1 to 185.12" and "Theorems 200.1 to 200.3". These lie in sections 124, 185 and 200,
   the first sections of the module chapters 5-1, 8-1 and 6-10 that 16.1's derivation cites. Under FIF those
   designators are the chapters IV-1, VI-1 and X-10.

Layer **J** (judgment, kept separate). It covers the 21 derivation groups that hold the requested statements and the
22-node component, plus 64.1/64.2, 66.1-66.5 and 131.1 on the paths between 124.1 and its premises: 55 boxes, 101
items. Each item is resolved to the box its description names, by a box title or box text that the script checks
verbatim. J only narrows M; it never adds a premise. Of the 101 items, 30 are narrowed, 37 already had J = M, 33 stay
unresolved (M kept), and 1 is dropped: "the gauge group ... as a comparison target and no input", which the
derivation itself marks as not an input.

**Test.** A THEOREM *violates* the rule if, for one of its items, every candidate (itself excluded) is C or O. It
violates with a *declared OPEN* premise if the item's own grade word is "open result". It is *undetermined* if an item
mixes T/D/Tab candidates with C/O candidates. Tab and instruments are neutral, per the front-matter sentence above.

## 3. Results

**R1. Shared derivations.** 2.1-2.9, 9.1-9.8 and 18.1-18.6 each have one derivation for the whole group, so the
text gives no premises per box for 2.2-2.8, 9.7, 9.8 or 18.1-18.6. Each procedure ends by naming one box ([5],
verbatim): "which is the open result" names 2.7 (`01-framework.tex:173`); "so the readability order is a typing of the
type, which is the theorem" names 9.1 (`:396`); "So the Born count is a measure and the readout its root, which is
the theorem" names 18.1 (`:695`). The other boxes have only the plain-words paragraph, which cites nothing.

**R2. Statement-level premises (J; M candidates in brackets where J narrows).**

| statement(s) | derivation | premises |
|---|---|---|
| 2.1-2.9 (2.2 C; 2.4, 2.5, 2.7, 2.8 O; others T) | `01-framework.tex:169` | Defs; 2.3 T ("the one-act clause as a theorem of the derived act") [8 T/C of chapter 1-1]; one T/C box of chapter 1-5, unresolved: H.3, P.1, R.1, U.1 T or H.9, Q.1 C. Index: Defs only |
| 9.1-9.8 | `:392` | Defs; 66.4 C ("the carrier as the closed read of its type"; title "The carriers as the closed reads of their type faces") [66.4 C, 69.1 T, 70.1 T; FIF: 54.1 C, 55.1 T]; an instrument outside the chapters. Index: Defs, 2.1-2.10 |
| 18.1-18.6 | `:691` | Defs; 9.1 C [9.1, 9.7, 12.1, all C]; **124.1 C** ("the explicit two-point function with its records (chapter 5-1, a theorem)"; 124.1's box: "the records of its two reads are 172") [124.1 C, 134.1 T, 135.1 T]; I.1 O, declared ("the sector-to-ends map ... an open result read at its standing") [15 O of chapter 1-5]. Index: Defs, 2.x, 9.x |
| 124.1 C | `07-cell-3-2-photon.tex:44` | Defs; one of 2.4, 2.5, 2.7, 2.8 (declared O); 9.1 C, 9.2 T, 9.4 T [9.1-9.7]; 18.2 T, 18.3 C ("by the gain law") [18.1-18.4]; one of 28.1, 28.2, 28.5 T; 35.1 C; 40.1 C; 50.1 C; 53.1 C; 57.2 O (declared). The cited sections are 1, 2, 9, 18, 28, 35, 40, 50, 53, 57; neither 16 nor 58 |
| 16.1 T | `01-framework.tex:582` | 9.1 C; "the count law at the six genres (chapter 8-1)", unresolved [185.1-185.12 C, 186.1 T, 191.1 T]; 200.6 T; 131.1 O, declared ("the stabiliser of the pair ket (chapter 5-1, an open result)"; title "The stabilizer of the pair ket and the unbroken charge") |
| 58.1 T | `04-spacetime-and-symmetry.tex:464` | 131.1 O (declared); 54.1 C [54.1 C, 55.1 T]; "the boost open result (chapter 3-3)", declared O, with no O box in chapter 3-3; "the dock open result and the dock-register open result" (chapter 3-5, declared O) [57.2]; "the photon row at the centre" (chapter 6-7, declared O) [7 O boxes] |
| 200.1-200.4 | `11-cell-4-4-higgs.tex:37` | chapters 1-2, 1-3, 1-4 cited as "definitions" (no box); chapters 1-6, 4-2 and 4-3, 6-3, cited with no grade word (chapter-level); 185.1 C |
| 200.5 O / 200.6 T | `:69` / `:93` | "the Higgs chapter above" (200.1-200.4 / 200.1-200.5); the framework / 195.1 T |

The other members of the 22-node index component (8.1, 27.1, 39.1, 69.1, 70.1, 74.x, 78.x, 146.1, 147.1) are in
`premises.txt` [4], item by item.

**R3. Re-grade of the requested THEOREMs.** There are 14: 2.3, 2.6, 18.2, 16.1, 58.1, 200.6, 8.1, 27.1, 39.1, 69.1,
70.1, 74.1, 74.4 and 78.3. The full table is [5].

| | index | derivation, M | derivation, J |
|---|---|---|---|
| list or declare an OPEN premise | 10 | 8 declared | 8 declared |
| violate (OPEN or CONDITIONAL) | 12 | 9 | 10 |
| undetermined | – | 5 | 4 (2.3, 2.6, 74.1, 74.4) |
| consistent | 2 | 0 | 0 |

The eight THEOREMs that declare an OPEN premise in their own derivation are 18.2 (I.1), 16.1 (131.1), 58.1 (four
items), 8.1, 27.1 (166.6), 39.1 (E.1; 200.4/200.5), 69.1 (200.5) and 70.1 (four items). 200.6 rests on 200.1-200.5
(C, C, C, O, O). 78.3 rests on 74.3 C and 124.1 C. Four verdicts differ from the index:

- 74.1 and 74.4: the index gives 14 C premises and 0 O; the derivation gives D.3/D.4 (T) and one T/C box of chapter
  1-5, so they are undetermined, not violations.
- 2.3 and 2.6: the index gives Defs only; the derivation also cites a chapter-1-5 theorem, so they are undetermined.

Over all 23 THEOREMs of the reconstructed groups, J gives 17 violations (8 declared OPEN), 6 undetermined and 0
consistent. The index gives 17 with an OPEN premise and 19 with an OPEN or CONDITIONAL one. 9.2-9.6 have 4 O premises in
the index, but none in their derivation; under J they violate through 66.4 C. 64.1 and 64.2 violate through 124.1 C and
40.1 C.

CONDITIONALs: 6 of the 17 in the reconstructed groups rest on an OPEN premise in their derivation (18.1, 18.3, 18.4,
124.1, 146.1, 66.4); the index gives 10. For 124.1 the derivation itself names two OPEN premises. The front-matter
clause "a composed claim carrying the lowest status in its chain" would therefore grade 124.1 OPEN. Its box says it is
"conditional on ... the results of the chapters named as its premises in the derivation below, each at its own
standing". The source does not say which of the two readings it intends.

**R4. The 124.1 cycles.**

1. *The index cycles do not survive.* None of the four index edges 124.1→16.1, 16.1→124.1, 124.1→58.1, 58.1→124.1 is
   cited in a derivation, in either M or J ([6]). 124.1's derivation cites sections 9 and 57, which are the first
   sections of chapters 1-2 and 4-5; the index expands these to the whole chapters, and so takes in 16.1 and 58.1.
   16.1 and 58.1 cite "an open result" of chapter 5-1. The grade word excludes 124.1 (C), and the description names
   131.1. Their index lists have "Theorem 124.1", the first box of that chapter. In the strict graph (J choices, single-candidate items elsewhere) neither 16.1 nor 58.1 lies on any cycle.
   58.1 is not a transitive premise of 124.1 even in the any-of graph (MOD, FIF; it is one under the union of the two
   readings). 16.1 is a transitive premise of 124.1 only through any-of chapter citations; the shortest such path,
   16.1 → 3.1 → 28.1 → 124.1, uses three ([6]).
2. *A different cycle through 124.1 does survive.* 124.1's derivation rests on "the theorem of chapter 18", i.e. one
   of 18.1-18.4. All six boxes of section 18 share one derivation, and it rests on "the explicit two-point function with
   its records (chapter 5-1, a theorem)", i.e. one of 124.1, 134.1, 135.1. 134.1 and 135.1 each cite "the
   fine-structure constant at the one link of chapter 124, a theorem" (a single candidate). So 124.1 rests on itself
   **whichever box the two any-of citations mean**. The script computes this as a fixpoint ([6]). Under every
   resolution of the any-of items, the component of 124.1 is {35.1, 40.1, 124.1}. 35.1 (state space) and 40.1 (Born
   rule) are single-candidate premises of 124.1, and both cite the theorem of chapter 18. Explicitly, writing p → s
   for "p is a premise of s", the cycles are 124.1 → 18.k → 124.1, 124.1 → 18.k → 35.1 → 124.1 and
   124.1 → 18.k → 40.1 → 124.1. Here 18.k is whichever of 18.1-18.4 the citation means, and the step 124.1 → 18.k may
   pass through 134.1 or 135.1. With J, the strict graph has
   exactly one cycle component, {124.1, 18.2, 18.3}; its two-cycle is 124.1 ↔ 18.2. Its two edges are the sentences at
   `01-framework.tex:691` and `07-cell-3-2-photon.tex:44`. Under J the forced component grows to {9.1, 35.1, 40.1,
   66.4, 72.1, 73.1, 124.1}.
3. *What the cycle depends on.* It needs the shared-derivation rule and the MOD reading of "chapter 5-1". Under FIF,
   5-1 is chapter IV-1 (64.1, 64.2), and 124.1 is on no forced cycle. 35.1, 40.1 and 44.1 are on one. The evidence for
   MOD is in §2.5.
4. *The any-of graph is no evidence either way.* Expanding every chapter citation to all its grade-compatible boxes
   gives strongly connected components of 261 (MOD), 208 (FIF) and 290 (union) statements. This is the same mechanism
   that produced the index cycles. Other forced cycles outside the request: 73.1, 98.1, 121.1 (MOD).

**R5. F3 and F4 recomputed.**

*F3, restricted:* see R3. *F3, all 100 THEOREMs, derivation level ([8]; J where judged, M elsewhere):*

| verdict | count |
|---|---|
| declare an OPEN premise | 38 |
| other definite violations | 35 |
| undetermined | 22 |
| consistent (2.10, I.7.1, I.52.1) | 3 |
| no premise sentence (I.1.2, I.16.2) | 2 |

Without J, under any of MOD, FIF or the union of the two: 38 / 29 / 28 / 3 / 2.

Of the index's 66 THEOREMs with an OPEN premise, 31 declare one in their derivation and 35 do not. Another 7 declare
one that the index does not list.

*F4* ([7]): transitive T+C dependents (of 209), index against derivation level. The upper bound is reachability in
the any-of graph. The lower bound counts dependence that holds whatever box each any-of citation means (fixpoint).

| premise set | index | upper (M / J) | lower, M | lower, J | direct users: index / J-definite |
|---|---|---|---|---|---|
| 2.2 C | 183 | 184 / 184 | 1 | 1 | 62 / 1 (13.1) |
| 2.4, 2.5, 2.7, 2.8 O | 183 | 185 / 185 | 153 | 153 | 58 / 14 |
| 9.1, 9.7 C | 175 | 183 / 183 | 0 | 151 | 120 / 8 |
| 9.8 O | 175 | 185 / 185 | 0 | 0 | 97 / 0 |
| 18.1, 18.3, 18.4 C | 152 | 182 / 182 | 0 | 0 | 105 / 0 |
| 18.5, 18.6 O | 152 | 185 / 185 | 0 | 0 | 89 / 0 |
| 124.1 C | 42 (70 all grades) | 184 / 184 | 152 (313 all grades) | 152 (313) | 22 / 25 (14 without J) |

## 4. Which earlier numbers change

- **Robust.**
  - Most THEOREMs violate the chain rule. The index gives 84 (OPEN or CONDITIONAL); the derivations give at least 73
    definite violations and 22 undetermined cases.
  - "Premise lists cite whole chapters" is now the source's own statement (`01-framework.tex:109`).
  - Kind words are unreliable in the derivations too: CONDITIONAL boxes are cited as "a theorem".
  - The OPEN axiom readings 2.4/2.5/2.7/2.8, taken as a set, still sit under at least 153 of 209 T+C statements: many
    derivations cite "the open result of chapter 2".
  - 9.1/9.7 sit under at least 151, given J.
- **Changed.**
  - "66 of 100 THEOREMs list an OPEN premise" is a property of the index. At statement level, 38 THEOREMs declare an
    OPEN premise in their own derivation.
  - The F4 loads of 2.2 (183), 9.8 (175) and 18.1/18.3/18.4 and 18.5/18.6 (152 each) are not supported statement by
    statement. Their lower bounds are 1, 0, 0 and 0: these boxes are load-bearing only as members of whole-chapter
    citations.
  - 124.1 carries more than the index shows. Its definite direct users are 25 (14 without J), against 22 in the
    index; its definite transitive dependents are at least 313, against 70 (152 against 42 among T+C). The reason is
    that the gain-lemma chapter, which most derivations cite, rests on 124.1.
  - The 124.1 cycles with 16.1 and 58.1 are artifacts of whole-chapter citation. The forced cycle {124.1, 35.1, 40.1}
    through chapter 18 is not.
- **Earlier statement corrected** (logged in `reports/CORRECTIONS.md`). `00-dependency-structure.md` §2 says that
  `designators.tsv` "resolves 1-2 and 4-5 uniquely, so the 124.1 cycles do not depend on that ambiguity". But
  `designators.tsv` gives 4-5 two readings (first_section III-5, module IV-5; mapped_by "fifteen"). With 4-5 read by
  the module numbering, the index two-cycle 124.1 ↔ 58.1 disappears (8366 edges, component still 22). With every
  chapter-form designator read that way, the component has 13 members and only 124.1 ↔ 16.1 remains ([2]). All index
  counts are reproduced unchanged.

## 5. Limits of the extraction rule

1. Shared derivations (61 groups) cannot be split by box. Their premises are attributed to every box of the group.
2. Derivations are themselves chapter-granular: 750 of 1069 items have several candidates. In the judged groups, 33
   items stay unresolved because no box title or text matches the description.
3. Only the premise sentence is read. Objects the derivation body uses without listing them (e.g. "the Ward identity
   holds at the constant functional" in 124.1) are not premises under this rule.
4. Qualifiers are not modelled in M: "read at its standing", "read at their faces", and "as a comparison target and
   no input" (2 items). J drops the one in a judged group (66.1-66.5); the other (84.1-84.3) stays in M.
5. 9 premises are located outside the chapters or at no chapter. Some cited chapters hold no box of the cited grade
   (e.g. "the boost open result" in chapter 3-3; "the wash move at the cap" in chapter 11-6).
6. If "theorem" were read strictly as THEOREM, the items citing 35.1, 40.1, 50.1, 53.1 or 124.1 as "a theorem" would
   point to a THEOREM that the cited section does not contain. Such a premise would be unlocated, not proven. The test
   above reads them as the CONDITIONAL boxes that are there.
7. The parser is written to this release's wording. The checks are: verbatim quotes (1069/1069), index numbering (0
   mismatches), and the derivation sentence of every judged item, printed in [4].
