# The Lean at scale: restated AlphaDerivation theorems, kernel modules, build and axioms (round 6)

Source: compendium 0.1, commit `a28db899`, directory `compendium/lean` (six packages, 915 `.lean` files). Object:
proposed further work item 6 of `reports/chain-theory-compendium.md`. `04-lean-coverage.md` classified the
declarations by kind and did not build them. Of the 64 Lean-tied THEOREMs it found that 11 cite only `opaque : Prop`
or open records, 37 at best `decide` over entered tables, and 16 a restated `AlphaDerivation` theorem or a kernel
module. This round asks three things about those 16: what the propositions say, how they are proved, and whether they
state the boxes they are tied to. It also asks whether the tree builds with core Lean 4.35.0, and what
`#print axioms` reports.

Scripts (all in `computation/chain-theory/compendium/`):
- `lean_resolve.py` and `lean_resolve_data.py`: the resolution. The data module holds the curated layer; every quote
  in it is checked verbatim against the source.
- `lean_build.py`: the build in a scratch copy. `lean_build_patch.tsv` lists the patch (five proofs, §5).
- `lean_audit/ResolveChecks.lean`, `lean_audit/Diagnose.lean` and `lean_audit/AxiomTrace.lean`: our own Lean files.
- `lean_audit.sh`: the whole pipeline.

Recorded outputs are in `outputs/`:
- `lean_resolve.txt` holds every number below, by section [1]-[10].
- `lean_resolve.tsv` has one row per tied THEOREM.
- `lean_build_A.*` and `lean_build_B.*` are the build records, and `lean_diagnose.txt` and `lean_axiom_trace.txt` the
  outputs of `Diagnose.lean` and `AxiomTrace.lean`.
- `lean_axioms.txt` and `lean_resolve_checks.txt` are the raw Lean output.

Metrics: `metrics/chain-theory/compendium/lean-audit-at-scale.yaml`.

## 1. The 16 rows and the 8 declarations

lean_classify.py is reused unchanged and reproduces the round-1 tiers, 16 / 37 / 11 ([1]). The 16 rows cite 8
distinct declarations: 7 restated `AlphaDerivation` conjunctions and 1 kernel check. Each declaration is cited by every
statement of its section, because the tie is per section.

| declaration | type | cited by THEOREM rows | box its content states ([5]) |
|---|---|---|---|
| `Foundations.Ch01.row_524038_ORKC` | 2 conjuncts (OrbitReading) | 2.1, 2.3, 2.6, 2.9, 2.10 | 2.2 (C) "Orbit reading from the knowledge clause" |
| `StandardModel.Ch10.row_518160_AGUQ` | 12 conjuncts (IdAlg, LemmaB, Theorem5) | 9.2-9.6 | 9.2 (T) "Abelian genre uniqueness" |
| `StandardModel.Ch10.row_520119_GWDC` | 1 (Census) | 18.2 | 18.3 (C) "Gain as the washed direction count at the sector" |
| `Foundations.Ch03.row_521104_OUTR` | 3 (IdAlg, Theorem1) | 56.1 | none in section 56 |
| `Foundations.Ch05.row_518830_MDFN` | 1 (LemmaB) | 74.1, 74.4 | 74.2 (C), its momentum-transfer clause |
| `Foundations.Ch05.row_521443_WDCL` | 5 (Theorem5) | 78.3 | 78.1 (C) "The Ward identity from conservation as one link" |
| `Foundations.Ch05.row_524064_WSWC` | 2 (Theorem5) | 78.3 | 78.2 (C) "The reading of the records washed" |
| `StandardModel.Ch12.row_527992_ESRL` | `Kernel.D527992_ESRL.checks` (8) | 196.2 | 196.2 (T) |

The four-letter suffixes are the initials of the boxes in the last column where these are legible (ORKC, AGUQ, GWDC,
WDCL). Five of the eight declarations state content of a CONDITIONAL box of the same section, and reach THEOREM rows
only through the section-level tie ([5]).

## 2. What the 26 conjuncts state and how they are proved ([2])

Each restated declaration is a conjunction of `type_of% AlphaDerivation.Statements.X`, proved by `And.intro` of the
same names. Each `Statements.X` re-exports a theorem of `AlphaDerivation/Proofs/*.lean`. The script quotes all 26
resolved statements verbatim. No conjunct uses `sorry`, `native_decide` or a declared `axiom`, and the library
declares no axiom. The proof methods are `decide +kernel` (18), an anonymous constructor of other theorems (4), and
short tactic proofs with `rw`/`simp`/`unfold` (4).

| kind (definition in `lean_resolve_data.py`) | conjuncts |
|---|---|
| COMPUTED: the kernel evaluates library functions on entered inputs | 13 |
| COMPUTED+LITERAL: as above, plus equalities between entered constants and literals | 6 |
| GENERAL: proved for all values of the variables | 5 |
| DEFINITIONAL+COMPUTED: true by a definition that is written to give it | 1 |
| ENTERED: a property of an entered table as entered | 1 |

The proofs compute. Unlike the data modules, the inputs here are matrices, letter lists, and web and record lists,
and the proofs evaluate functions on them. Some of the computation is real linear algebra. `IdAlg` builds the 12
Gaussian-integer 5x5 matrices of su(2) + su(3) + u(1). It shows four facts: the matrices close under the bracket; the
centre is one-dimensional (a vector in the kernel and an 11x11 minor nonzero); TU = diag(3i,3i,-2i,-2i,-2i) acts on
C^2 (x) C^3 as i·Id_6; and the commutant of the tensor representation is the scalars (a 35x35 minor nonzero among 432
constraints). At the Lie-algebra level this overlaps G1/G2 of `03-gauge-group.md` at the genre (2,3).

Where the computation stops:

1. **The (2,2) centre is not computed.** `IdAlg_centre_acts` uses an entered TU22 = diag(i,i,-i,-i). The (2,2) centre
   dimension is entered as `Expected.centreDims = [1, 1]`, and no theorem uses it.
2. **"Abelian", "conserved", "massless", "spin one" and "electromagnetic-like" are Boolean definitions in
   `LemmaB`.**
   - `abelian g := pairScalar g != 0`, where `pairScalar g = (d2 - d1)/gcd`.
   - `conserved g := generatorBlocksScalar g && centralFunctional g != (0, 0)`. Its first clause is true at every pair
     (`Round6.generatorBlocksScalar_always`, our check). `abelian`, `conserved` and `d1 != d2` coincide on all pairs
     in 1..12 (`Round6.abelian_conserved_distinct`).
   - The spin labels come from the entered frame contents 2 ↦ [1], 3 ↦ [0,0,0], 4 ↦ [0,2].
3. **Five opaque hypotheses are never used.** `Theorem5_correspondence_at_structure` takes `ReadCapWoronowiczStormer`,
   `SoftPhotonTheorem`, `Decoupling`, `OnShellDefinition` and `OnShellIdentification`, all `opaque : Prop`. None is
   used in the proof. Its conclusion holds without them (`Round6.correspondence_without_hypotheses`).
   `Theorem4.onShellCoupling o` is defined as `(registerSum, register)` when `o` is the object of record and as
   `(Q4.zero, 1)` otherwise (`Round6.onShellCoupling_elsewhere`). So "the on-shell coupling equals the closed form
   (79, 48, 42, 3)/31104 at (2,3)" holds by that case split. The numerator is computed from the 172 entered records.
4. **The zero momentum transfer is entered.** `LemmaB_zero_transfer` checks the step table of
   `Data/Records.lean`, 172 rows. The third triple is (0,0,0) in all 172 as entered. The second triple is minus the
   first in all 172 (first triples: 82 × (0,0,0), 48 × (1,-1,0), 42 × (-1,1,0)). See also
   `Round6.steps_entered_zero`.
5. **The 172-record check does not depend on the records.** `Census_weights_from_construction` checks
   ((2-x)(3-y))/2 = 3^[x=-1]·2^[y=-1] at every record. That identity holds for every (x,y) ∈ {±1}²
   (`Round6.weight_identity`). The content is the exchange rule (same-kind records: (-p, 1); mixed: (1, -p); closure:
   ClosureSign) applied to the entered parity and kind fields. It gives the histogram (1,79), (2,48), (3,42), (6,3).
   This is the end rule of the alpha audit (CT-G02).
6. **The "Ward identity" theorems are elementary.** They are telescoping sums of differences (`telescoping`),
   c + (a-b) = c ⟺ a = b (`one_link_iff`), the constant functional (`unit_washed`) and (x-y)+(y-z) = x-z. Over the 41
   enumerated readings, the identity step holds for the unit functional on every path by `unit_washed`
   (`Round6.identityAt_unit`).
7. **Literal parts.** Six equalities hold by `rfl` (`Round6.literal_parts`): `Expected.wardRecordIdentity = "172 of
   172"`, `Expected.wardLinkPaths = "48 of 48"`, `Expected.endMap = "X:centre:3;Y:far:2"`,
   `Expected.tensorCommutantDim23 = 1`, `Expected.killingRank23 = 11`, and [1,1,-1,-1] = `Expected.centreGenerator22`.
   `Theorem1_readout_is_trace_form` does not connect the two constants 1 and 11 to the ranks it computes.

Which parts of `AlphaDerivation` are restated by any chapter declaration (restatements / statements):
- Never restated: the costly group and commutant computations (LemmaA 0/6, CollapseG1-G3 and D2 0/1 each, Generation
  0/1, Span 0/2, SpanRank 0/2). CollapseG1-G3, Generation, Span and SpanRank are 6 of the 8 modules that the source's
  build record lists at 4 GB or more; the other two are IdentificationAlgebra (restated) and
  KernelDerivationC.Brackets (behind the kernel module LGIS). Also never restated: GainLemma 0/4, OneLinkWord 0/4 and
  Correspondence 0/3.
- 18.2 is the gain lemma. Its "four end-set sums" and "Born count" are stated in `GainLemma.end_set_sums` and
  `born_count_is_a_measure`, which no chapter declaration restates.

## 3. Do the propositions state the boxes? ([4])

Method:
- Each box is split into its clauses, quoted verbatim: 52 quotes, all found.
- Each clause is given the conjuncts that state it, if any.
- MATCH means every clause has a counterpart. PARTIAL means some clause does. NONE means no clause does.
- Numbers of 5 or more printed in the box are looked up in the tied statements and in the whole Lean tree. This lookup
  is literal, so it also finds coincidences.

| row | tied | verdict | what is stated / missing |
|---|---|---|---|
| 2.1 T | ORKC | PARTIAL 2/3 | the 3 re-entry webs are commutators +(x,y)−(y,x) at the dimension-2 slot; one web is fixed by every generator. "a corollary of the Axiom at both clauses" is not in Lean |
| 2.3 T | ORKC | NONE 0/3 | one-act clause, census of moves, four one-pair acts |
| 2.6 T | ORKC | NONE 0/2 | "thirty-two derived cells" (32 does not occur in the tied statements) |
| 2.9 T | ORKC | NONE 0/2 | one-dimensional end |
| 2.10 T | ORKC | NONE 0/3 | single read |
| 9.2 T | AGUQ | PARTIAL 4/5 | all clauses except "dimension one at (2,2)" (not computed; item 1 above) |
| 9.3 T | AGUQ | PARTIAL 1/3 | the (2,3) value by the case split; "printing the other types as exact faces" is not stated (value (0,1)) |
| 9.4 T | AGUQ | MATCH 3/3 | no mass face, spin one (definitions of §2 item 2), zero transfer (entered, item 4) |
| 9.5 T | AGUQ | NONE 0/2 | hub count monotone against hub dressing |
| 9.6 T | AGUQ | NONE 0/2 | one-link tables at six types (no OneLinkWord statement is restated) |
| 18.2 T | GWDC | PARTIAL 2/4 | product of passed-end dimensions 1/2/3/6 at the 172 records (by definition, item 5); end-set sums and Born count are absent |
| 56.1 T | OUTR | NONE 0/5 | 48 matchings, 0 of 322714 records, 3072 pairs, 64 webs, 319332, 384 vs 48. **322714 and 319332 occur nowhere in the Lean tree**, although the box says "the deciding table re-derived in Lean at the forty-eight matchings". The 9 and 24 found in OUTR are the registers 9·4 and 24·36·36, coincidences |
| 74.1 T | MDFN | NONE 0/3 | left kernel of the move matrix, rank zero, broken counts |
| 74.4 T | MDFN | NONE 0/4 | junction count, broken structure and triangle counts |
| 78.3 T | WDCL, WSWC | NONE 0/4 | body/cross/self-term split, 67 with the sign plus, order-three term zero (67 is not in the tied statements) |
| 196.2 T | ESRL | PARTIAL 2/4 | every site typed `openAt` (computed); "twelve of the seventeen" (check_6 gives 12; 17 appears only in `derived_bodies_counts`, outside `checks`); the joined face and its orders are absent |

Totals: MATCH 1, PARTIAL 5, NONE 10. For 10 of the 16 rows, the declaration that made them "restated or kernel" has
nothing in common with the box. In all 10 the declaration is that of another box of the section, or of none.

## 4. The kernel modules ([3], [6])

All three kernel modules compute their rows. Every one of the 147 data modules with a `Row` structure enters its rows
as a literal list, 135 of them with only `String` fields.

| module | rows | checks over | cited by |
|---|---|---|---|
| `Kernel.D524748_SCGI` | bodies × faces, registers support^4 = 16, 81 (KernelDerivationC.Objects) | computed rows, 5 conjuncts | 172.1 C, 172.2 O |
| `Kernel.D525490_LGIS` | 41 cells derived in KernelDerivationC.Brackets (a 5.8-GB module): the brackets of the basis J1-J3, K1-K3, the Jacobi count, the Killing signature 3,3 | the target list after `rw [rows_eq_record]`; the equation derived = target is computed | 54.1 C |
| `Kernel.D527992_ESRL` | 69 site rows of 5 bodies built by applying the act library to entered step lists; open = `openAt` | computed rows, 8 conjuncts | 196.1 C, 196.2 T |

The inputs are entered: `KernelDerivationC.Inputs` has 26 definitions (sector classes, arms, registers, printed
strings), and ESRL's step lists, kits and permutations are entered. `lean_classify.py` labels ESRL
`kernel_check[strings]`, because its `Row` fields are strings, but its `"YES"`/`"NO"` labels are computed.
`04-lean-coverage.md` §1 says that "their `checks` are `decide` proofs that given string labels ... occur in a
hand-entered list". That is so for 135 of the 137 string-only modules, but not for LGIS and ESRL.

## 5. Build with core Lean 4.35.0 ([7])

1. **What the tree needs.**
   - Six lakefiles. Their requirements are relative paths only: no git source, no Mathlib.
   - All six packages pin `leanprover/lean4:v4.32.1`.
   - The 897 library modules import only the six packages' own libraries. `import Lean` occurs only in six script
     files outside the libraries.
   - The safety gate of `lean_build.py` found no `#eval`, `IO.`, macro or elaborator definitions, `extern`, `unsafe` or
     `native_decide` in the 493 modules, and the tree declares no `axiom`.
   - So core Lean alone suffices in principle. 4.35.0 is not the pinned version.
2. **Build A: the published sources, unchanged, fail under 4.35.0.**
   - `AxiomDeclaration.Theorems` fails at lines 42 and 52 (`failed to synthesize Decidable (worldThree.sharedQudits =
     [1])`, and the same for `passedTo`).
   - Every package requires `AxiomDeclaration`. 277 of the 493 modules of the source's `lake build` import the failing
     module, so that build cannot complete.
   - Cause (`lean_audit/Diagnose.lean`, `outputs/lean_diagnose.txt`): 4.35.0 instance resolution does not unify
     `Fin worldThree.qudits` with `Fin 3`, because `worldThree` is a plain `def`. The instance exists at the stated type
     (`#synth DecidableEq (List (Fin worldThree.qudits))`). `fail_if_success decide; rfl` proves both statements.
3. **Build B: copies with five proofs re-tacticked.**
   - With the first module patched, `AxiomDeclarationTwo.Theorems` fails the same way at lines 42, 52 and 56, blocking
     3 more modules.
   - `lean_build_patch.tsv` replaces `by decide` by `rfl` or `⟨rfl, rfl, rfl⟩` at these five lines. No statement
     changes. Of the library modules, the five theorems appear only in the two that define them; elsewhere they appear
     only in the source's `AxiomsCheck`/`Roots` scripts.
   - **Result: 493 of 493 modules compile, with 0 errors, 0 warnings and 0 `sorry` warnings.**
   - 215 modules import no patched module: all 39 of `AlphaDerivation`, and the chapter modules Ch01 and Ch10, which hold
     ORKC, AGUQ and GWDC. For them Build B compiles the published source unchanged. OUTR, MDFN, WDCL, WSWC and ESRL, and
     125 of the 195 names cited by the 64 THEOREM rows, lie in modules that import a patched module.
   - Cost: 5398 s summed over modules (source record 3845 s). Wall time 2494 s in three resumed runs at 3 jobs, on a
     4-core machine shared with other work. Peak 7677 MB (`CollapseG2`; source 5926 MB).
4. **Not done.** The pinned 4.32.1 could not be fetched: this session's permission policy refused the download. So the
   source's build record is confirmed at 4.35.0, up to the five proofs, and not on its own toolchain.

## 6. Axioms ([8], [9])

1. **Coverage.** Our file has `import ChainTheory`, `import AlphaDerivation` and 6662 `#print axioms` lines. It is
   emitted by `lean_resolve.py --emit-axioms-file` and run against Build B (`outputs/lean_axioms.txt`).
   - Its names are exactly those of the source's own `AxiomsCheck.lean`.
   - That list already contains every name of the `lean` column (all grades), the 26 conjuncts, their Proofs-level
     theorems and every kernel-module declaration.
   - All 6662 names are reported, with no error.
2. **sorryAx: 0.** No axiom occurs outside `propext`, `Quot.sound` and `Classical.choice`. This holds for all 6662
   names, so for the 195 names of the 64 THEOREM rows and for every conjunct.
3. **The axiom profile depends on the toolchain.**
   - At 4.35.0: none 5965, [propext] 290, [propext, Quot.sound] 4, [propext, Quot.sound, Classical.choice] 403.
   - The source's `axioms_report.txt` (4.32.1) gives none 6047, [propext] 605, [propext, Quot.sound] 10, and
     `Classical.choice` 0. These figures are its own and are reproduced exactly by our parse.
   - 6259 axiom sets are equal. The other 403 differ only by gaining `Classical.choice`: 315 had [propext], 82 had
     none, 6 had [propext, Quot.sound].
   - `lean_audit/AxiomTrace.lean` (`outputs/lean_axiom_trace.txt`) gives the shortest chain from five source
     declarations to `Classical.choice`. It runs through core 4.35.0: `Int.instLinearOrderPackage` /
     `Fin.instLinearOrderPackage` → `Std.LinearOrderPackage.ofLE._proof_1` →
     `Std.instLawfulOrderMinOfIsLinearPreorderOfLawfulOrderLeftLeaningMin` → `Classical.propDecidable`.
   - The toolchain-independent findings are therefore: no `sorryAx`, and no axiom of the source's own.
4. **By class, the 195 THEOREM-tied names at 4.35.0** (at 4.32.1: 188 none, 7 [propext]):

   | class | no axiom | all three |
   |---|---|---|
   | opaque Props | 69 | 0 |
   | open records | 69 | 0 |
   | aliases | 3 | 0 |
   | conditional table checks | 16 | 0 |
   | table checks | 27 | 3 |
   | restated | 0 | 7 |
   | kernel | 0 | 1 |

   - Per row: at 4.35.0, 23 of the 64 rows cite a declaration with the three axioms and 41 cite only axiom-free
     declarations. At 4.32.1 the split is 14 [propext] and 50 none.
   - The 8 restated/kernel declarations use the three axioms at 4.35.0. At 4.32.1 they use [propext], except MDFN,
     which uses none.
5. **Opaque usage per theorem** (`outputs/lean_resolve.tsv`, columns 7-9). An `opaque X : Prop` is a constant, not an
   axiom, so `#print axioms` does not list it. Across the 64 rows:
   - 35 rows cite opaque Props (145 with multiplicity);
   - 31 rows cite table checks under an opaque hypothesis (37);
   - 42 rows cite open records (125).

   The only opaque hypotheses inside the 26 conjuncts are the five unused ones of §2 item 3.
6. **Own checks** (`outputs/lean_resolve_checks.txt`). All 8 theorems of `ResolveChecks.lean` compile.
   `literal_parts` uses no axiom. The other 7 use the three, inherited from the same core instances.

## 7. The coverage table re-stated ([10])

| tier | round 1 (`04-lean-coverage.md`) | round 6 |
|---|---|---|
| cites a restated AlphaDerivation theorem or kernel module | 16 | 16 (kind unchanged; all of them compute) |
| — of which the computation states the whole box | — | 1 (9.4) |
| — of which it states part of the box | — | 5 (2.1, 9.2, 9.3, 18.2, 196.2) |
| — of which it is unrelated to the box | — | 10; their next-best declaration: table checks for 9 (2.3, 2.6, 2.9, 2.10, 9.5, 9.6, 74.1, 74.4, 78.3), opaque Props only for 1 (56.1) |
| best evidence: data-table checks (not resolved here) | 37 | 37 (+9 = 46 if the 10 unrelated rows are re-tiered) |
| only opaque Props / aliases / open records | 11 | 11 (+1 = 12) |

What changes in `04-lean-coverage.md`:
- The 11/37/16 counts are unchanged as a classification by declaration kind.
- Read by content, the top tier is 6 rows, not 16: 1 full and 5 partial statements of the box.
- "Only those can certify any computation" (report, item 6) is right about kind: all 16 compute. But for 10 rows the
  computation certified is that of another box.
- §1's gloss "hand-entered list" fails for 2 of the 137 string-only modules (§4 above).
- "Build not re-run" is superseded by §5-6.

No count of rounds 1-5 is wrong; one gloss is corrected (`reports/CORRECTIONS.md`, 2026-10-30: the 04 "hand-entered" gloss and "Build not re-run").

## 8. Limits

- The verdicts are judgments at clause level. They are recorded clause by clause in `lean_resolve_data.py`, and their
  quotes are checked. A clause counts as stated when a conjunct states its finite content under the library's own
  definitions of its words, whether or not those definitions capture the words (§2 items 2-5).
- The 37 table-check rows and the 11 opaque rows were not content-checked here. The 10 re-tiered rows are put on their
  next-best kind, not content-checked.
- The tree was not built with its pinned Lean 4.32.1 (download refused in this session). Build B differs from the
  published sources in five proofs. Axiom sets are measured at 4.35.0 and compared with the source's 4.32.1 report; the
  report itself is not re-derived.
- `ActClosureCapFree` was not built. Its 402 library modules are not imported by `ChainTheory` and lie outside the
  source's 493-module build.
- Build B was recorded in three runs on one scratch directory: a restart to scale the memory estimates, and a resume
  after the second failing module. `outputs/lean_build_B.txt` keeps all three logs, and `lean_build_B.tsv` records one
  completed compile per module. A fresh `lean_audit.sh` run does this in one pass.
