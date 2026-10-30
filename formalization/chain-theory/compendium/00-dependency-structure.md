# Dependency structure of the compendium (tasks 1-2)

Data: `compendium/index/statements.tsv` (468 rows, CRLF), `designators.tsv`, section numbering re-derived
from `compendium/tex`. Scripts: `computation/chain-theory/compendium/census.py`, `depgraph.py`
(outputs `outputs/census.txt`, `outputs/depgraph.txt`, `outputs/depgraph_no_chapter_form.txt`).
Interpretive choices C1-C4 are listed in the computation README; C3 (resolution of 23 premise strings of
the form "the theorems of Chapters X-Y") is the only one that changes any number below, and both runs
are reported.

## 1. Census (reproduced)

| | THEOREM | of which Lean-tied | CONDITIONAL | OPEN | TABLE | DEFINITION |
|---|---|---|---|---|---|---|
| Parts I-XVI | 83 | 59 | 106 | 162 | 26 | 6 |
| Appendices | 6 | 3 | 2 | 23 | 0 | 0 |
| Instruments | 11 | 2 | 1 | 13 | 28 | 1 |
| **all** | **100** | **64** | **109** | **198** | **54** | **7** |

Per Part (T / C / O): I 20/9/23; II 0/0/0; III 1/8/4; IV 6/5/4; V 14/6/24; VI 2/17/30; VII 5/10/9;
VIII 1/8/17; IX 1/12/12; X 6/18/8; XI 2/3/2; XII 1/9/16; XIII 7/0/5; XIV 6/0/3; XV 11/1/5; XVI 0/0/0.
The six cells (VI-XI): 17/68/78, as the README states.
The headline 100/109/198/64 is reproduced exactly. The README's "Tables 56" is not: the index has 54.
"Tied to Lean" is a per-section property: the `lean` column is the section's list, inherited by every
statement in it (0 of 284 sections carry differing lists), so a statement is "tied" when its section
names any declaration (see `04-lean-coverage.md`).

## 2. Graph

8367 premise edges (7552 without chapter-form strings). Premise lists are chapter-granular: e.g.
"Theorems 18.1 to 18.4 and Open results 18.5 and 18.6" cites every box of section 18. Kind words are
unreliable: "Theorem 124.1" and "Theorems 18.1 to 18.4" include CONDITIONAL statements.

**Reciprocity.** Only 3632 of 8367 premise edges (43.4%) are echoed in the premise's own `used_by`
field; 3632 of 6087 `used_by` edges are echoed in the user's premises. The two directions of the
bookkeeping are inconsistent for the majority of edges.

**Dangling references.** Premise text cites 238.1, 239.1, 240.1 (appendices A, B, C), which have no boxed
statement in the index.

**Cycles.** Under C3 there is one strongly connected component of 22 statements, containing 124.1 (the
alpha theorem), 16.1, 58.1, 69.1, 74.1-74.5, 78.1-78.3, 146.1, 147.1, 200.1-200.4, 27.1, 39.1, 8.1, 70.1.
Two-cycles: 124.1 <-> 16.1 and 124.1 <-> 58.1 (16.1 and 58.1 list "Theorem 124.1" explicitly; 124.1 lists
"the theorems of Chapters 1-0 to 1-3, ... 4-5", whose ranges contain sections 16 and 58), and
16.1 <-> 200.1-200.4, 69.1 <-> 200.1-200.4. Without chapter-form strings the graph is acyclic; the
cycles are therefore carried by the chapter-level citations of 124.1 and 200.x. Separately, 82 of the
162 chapter designators have two readings (two numberings); `designators.tsv` resolves 1-2 and 4-5
uniquely, so the 124.1 cycles do not depend on that ambiguity.

## 3. The chain rule

Front matter (`tex/front.tex`, "The status words"): "every theorem prints at its statement the union of
the conditions of the results it rests on, a composed claim carrying the lowest status in its chain."

| Test | Count (C3 applied) | Count (no chapter form) |
|---|---|---|
| THEOREMs listing >= 1 OPEN premise directly | 66 / 100 | 66 / 100 |
| THEOREMs listing >= 1 OPEN or CONDITIONAL premise directly | 84 / 100 | 78 / 100 |
| THEOREMs with an OPEN statement among transitive premises | 81 / 100 | 73 / 100 |
| CONDITIONALs listing >= 1 OPEN premise directly | 92 / 109 | 92 / 109 |

Either the premise fields are not logical premises, or two thirds of the theorems are mis-graded by the
compendium's own rule. Physics theorems affected include 64.1, 64.2, 65.1 (gauge), 74.1 (conservation),
83.1 (three generations), 145.1 (Casimir; rests on 143.1), 217.1 (second law), 220.3, 221.1, 225.1.

## 4. Load-bearing unproven premises

Ranked by the number of THEOREM+CONDITIONAL statements (of 209) resting on them transitively:

| Premise(s) | Grade | Title | Transitive T+C (C3 / no chapter form) |
|---|---|---|---|
| 2.2 | CONDITIONAL | Orbit reading from the knowledge clause | 183 / 168 |
| 2.4, 2.5, 2.7, 2.8 | OPEN | One act clause discharge; Record clause from the one act read; The act alphabet; The unknowability clause at the axiom | 183 / 168 |
| 9.1, 9.7 | CONDITIONAL | Readability order at the six genres; Universe link from parity | 175 / 160 |
| 9.8 | OPEN | Mixed genre functional at the census | 175 / 160 |
| 18.1, 18.3, 18.4 | CONDITIONAL | The read-state boundary; Gain as the washed direction count; The boson as the chain readout | 152 / 136 |
| 18.5, 18.6 | OPEN | Role from the arrow passage; The readout law along the sector route of records | 152 / 136 |
| D.1, D.2, D.5-D.7 | OPEN | the move decomposition, orientation and move-graph results (appendix D) | 142 / 125 |
| 28.3, 28.4 | OPEN | The mass unit spelling as a theorem; The dimension census of the ledger | 127 / 110 |
| H.9 / H.1, H.2, H.4-H.8 | CONDITIONAL / OPEN | saturation: cut rank, closure, network-list seeds, ... (appendix H) | 121 / 102 |

The five section-2 statements (chapter 1-1, readings of the axiom clauses) sit under 88% (183/209) of all
graded claims. The
second tier is the readability typing of the six genres (9.1), which also fixes the read cap used to
select the (3,2) genre (`03-gauge-group.md`, G3), and the readout law of chapter 18.

## 5. The alpha theorem

124.1 is listed as a premise by 22 statements (6 THEOREM, 6 CONDITIONAL, 9 OPEN, 1 TABLE), among them
141.1 (anomaly), 143.1 (hydrogen), 146.1 (finiteness), 147.1 (running), 58.1 (spin-c), 61.1 (three
involutions), 16.1. 70 statements rest on it transitively (32 without chapter-form strings). These
statements inherit whatever the coordinator's audit of the Letter finds.
