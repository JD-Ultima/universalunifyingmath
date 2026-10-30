# A rule space generated from the Lean definitions (F-CT-4, next step 4; gap CT-G10)

**Object.** The readout of the α Letter as `ChainAlpha.lean` defines it (Wissner-Gross, v1.0.0, commit `a28db89`, file sha256 `522844e9…c331da`): `readout objectWorld objectTwoPoint = ([(1,79),(2,48),(3,42),(6,3)], 31104)`, i.e. $(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$, $1/\text{value}=137.035\,999\,189\,566\ldots$
**Question.** F-CT-4 B measures sensitivity over a rule space that we declared by hand: 46 080 rules over the record classes, 15 639 distinct values, density 103.9 per unit of $1/\text{value}$ in $(120,155)$. In its illustrative random-target model (target $U$ uniform on $(120,155)$, rule outputs fixed) the expected number of distinct outputs within $\pm1\sigma$ of $U$ is $4.4\times10^{-6}$; this is a reference count, not a calibrated p-value. How do these quantities change when the rule family is generated mechanically from the Lean definitions by a declared set of operators, rather than by our choice of which alternatives are "natural"?
**Code.** `computation/chain-theory/alpha/mechanical_mutations.py` → `mechanical_mutations.out` (generator, evaluation, statistics). `chainalpha_transcription.py` holds the Python transcription of the Lean definitions. `proofs/lean/chain-theory/mutation_spotcheck.sh` → `mutation_spotcheck.out` runs Lean itself on the generated mutants. "out §n" means section n of `mechanical_mutations.out`.
**Date.** 2026-10-30. This is a rule-family sensitivity analysis; the audited value is not changed.

## Result in brief

1. **The space.** Fifteen declared operators are applied at every site that a token scanner finds in the 44 definitions the readout depends on. This gives 818 single mutants and 331 490 pairs of mutants with disjoint sites, 332 308 in all. The operators, the pair rule and these counts are printed in out §1, before any value is computed in the run. The version history does not record the order in which the operator set and the values were fixed: the generator, the transcription, the output and this note entered the repository in one commit (corrected after review 2026-10-30, R4-CT-5; see `reports/CORRECTIONS.md`).
2. **Validation.**
   - The transcription reproduces the unmutated readout, and all 20 lines of `upstream_mutations_lean.out` (BASE and the author's 19 rewordings) character for character, including the `main` and `commutators` verdicts (out §2).
   - Lean 4.35.0, run on every single mutant, agrees with the transcription on all 818: 770 readouts identical, 47 predicted elaboration failures confirmed, and 1 resource abort. On a seeded sample of 300 pairs it also agrees on all 300: 266 readouts identical, 31 elaboration failures and 3 resource aborts (`mutation_spotcheck.out`).
3. **Values.** 272 988 mutants are well-defined. They take 10 340 distinct exact values, 1105 of which have $1/\text{value}\in(120,155)$. That is a density of 31.6 per unit. In the illustrative random-target model of F-CT-4 B the expected match count within $\pm1\sigma$ is $1.3\times10^{-6}$, a reference count and not a calibrated p-value. The same code applied to B returns F-CT-4's 103.9 per unit and $4.4\times10^{-6}$. The density depends on the mutation depth: 1.09 per unit for single mutants and 31.6 with pairs, and a seeded sample of triples gives a low estimate of 58 at depth 3 (§4). The comparison with B is therefore a comparison at depth ≤ 2 (corrected after review 2026-10-30, R4-CT-1; see `reports/CORRECTIONS.md`).
4. **Closest values.** The only mutant value within relative $10^{-4}$ of $\alpha^{-1}$ is the Letter's own; 29 292 mutants leave it unchanged. The next closest is $(149+45\sqrt3)/31104$, with $1/\text{value}=137.0569$, $9.9\times10^{5}\sigma$ away. In B the nearest other value is $1.3\times10^5\sigma$ away (`freedom.out`).
5. **Answer.** At mutation depth ≤ 2 the mechanical family gives the same picture as B: the Letter's value is the only one within relative $10^{-4}$, the next value is $9.9\times10^5\sigma$ away, and the density near $\alpha^{-1}$ is lower than B's (31.6 against 103.9 per unit). The reference counts of the illustrative model are $1.3\times10^{-6}$ for the mechanical family, $4.4\times10^{-6}$ for B and $5.6\times10^{-6}$ for their union. They are not calibrated p-values, and they do not determine an effective trial factor or how many further rules would make a match likely. The density grows with depth, so "lower than B" holds at depth ≤ 2 and is not a bound; the first version's "the F-CT-4 conclusion is robust to replacing B by the mechanical space" is withdrawn (corrected after review 2026-10-30, R4-CT-1; see `reports/CORRECTIONS.md`). Its threshold "mutants for an $O(1)$ chance hit" is removed (corrected after review 2026-10-30, R4-CT-6; see `reports/CORRECTIONS.md`). The mechanical family does not contain everything B or Table S5 contains (§6). It also says nothing about choices made outside the Lean text (CT-G10).

---

## 1. Inventory

The readout `readout objectWorld objectTwoPoint`, the last conjunct of `main`, is computed from the definitions of Part 1 of the file and the utilities of Part 0. Following identifier references transitively from `readout`, `objectWorld` and `objectTwoPoint` gives 44 of the 52 Part 1 definitions (out §1, with site counts per definition):

| Role | Definitions |
|---|---|
| The object (constants, list literals) | `objectWorld` (4 qudits; dimension 3 at qudit 0, else 2; links `[(0,1),(0,2),(0,3)]`), `objectTwoPoint` (centre 0, reads `(0,1)`, `(0,2)`, state `(0,3)`) |
| End rule (Boolean predicate, case split) | `pairKet`, `steered`, `fixesFar`, `readEnd`, `farEnd` |
| Normalisation | `readoutDim` |
| Letters and slots | `tracelessBasis`, `slotEnd`, `slotDim`, `slotBasis`, `exch`, `letterAt`, `letterSlots`, `exchAssign`, `letterCode`, `key` |
| Record enumeration | `products`, `coeffOf`, `parity`, `weylPhaseExp`, `weylIsSign`, `coeffOn`, `commAt`, `orderSign`, `commutatorCandidates`, `sameEnd`, `commutatorRecs`, `recs` |
| The two exchanges and survival | `slotQudit`, `linkSwap`, `occupies`, `sameSet`, `moved`, `sameForm`, `survives`, `sigX`, `sigY` |
| Weight function | `wX`, `wY`, `kOf` |
| Readout | `addRoot`, `readout` |

- **Held fixed.** The Part 0 utilities the cone calls are `rank`, `neZ`, `eqZ` and `Mat.mulZ`, `Mat.pairZ`, `Mat.subZ`, `Mat.tabulate`, with what these call (`rankAux`, `Mat.get`, `Mat.traceZ`, `Mat.transposeZ`). They implement exact arithmetic, not rules, so they are not mutated.
- **Outside the readout.** $52-44=8$ Part 1 definitions are outside the cone. Four are used by the other conjuncts of `main`: `TwoPoint.valid`, `records`, `World.typeOf` and `World.wellFormed`. The other four, `World.isPair`, `World.valence`, `World.shared` and `commutatorRecords`, are used by neither the readout nor `main`; some appear in other theorems of the file (`xcheck_mechanical.out` [1]). Mutating any of the 8 cannot change the readout (corrected after review 2026-10-30, R4-CT-7(b); see `reports/CORRECTIONS.md`).
- **Types.** `World`, `TwoPoint`, `Assign`, `Formal` and `Rec` are not mutated.

The scanner finds 425 sites in the 44 bodies: 83 literals, 4 unary negations, 4 Boolean literals, 50 tests, 7 quantifiers, 12 connectives, 33 if-then-else, 7 matches, 32 arithmetic operators, 180 label occurrences, 6 list literals, 4 `++` chains and 3 `::` chains (out §1).

## 2. Declared operators and the space

| Operator | Applied to every | Variants |
|---|---|---|
| LIT | integer literal (a negative literal is one literal) | value − 1, value + 1 |
| NEG | unary minus or Boolean `!`; sign of a negative literal | deleted |
| LSW | pair of literals of one type and different value within one definition | swapped |
| BLIT | Boolean literal | `true` ↔ `false` |
| REL | test `=`, `≠`, `==`, `!=`, `eqZ`, `neZ`, `isSome` | negated; `<` → `≤`, `≥` and `≤` → `<`, `>` |
| QNT | list quantifier | `all` ↔ `any` |
| LOG | connective `&&`, `\|\|`, `∧`, `∨` | dual connective; left operand dropped; right operand dropped |
| IF | if-then-else | branches swapped |
| MATCH | pair of match arms | right-hand sides swapped |
| AOR | binary `+ − * / %`; `min` | each of the other four; `max` |
| VAR | occurrence of a declared label | its partner, where the partner is in scope |
| ENDSWAP | declared label pair within one definition | all occurrences exchanged |
| LIST | entry of a list literal | dropped; duplicated |
| CAT | summand of a `++` chain | dropped; duplicated |
| CONS | element of a `::` chain | dropped; duplicated |

The declared label pairs are the end and role pairs of the file: `dA/dB`, `d1/d2`, `dc/df`, `i/j`, `a/b`, `s1/s2`, `i1/i2`, `c1/c3`, `f/g`, `xs/ys`, `odd/even`, `read1/read2`, `readEnd/farEnd`, and `.1/.2`.

A **pair** is an unordered pair of single mutations whose sites are disjoint. A site is the token a mutation rewrites, or, for operators that move or drop a subterm (IF, MATCH, LOG drops, LIST, CAT, CONS), the whole node. So no mutation is combined with one anywhere inside a node it moves or drops. Besides a mutation in a branch that an IF swap moves, this excludes a mutation in the swapped IF's own condition (152 pairs), one in the operand that a LOG drop keeps (85 pairs), and one in an entry that a list, `++` or `::` operation keeps. Of the 2663 excluded pairs, 1367 share a site (`xcheck_mechanical.out` [8], counted with the generator's own `compatible`) (corrected after review 2026-10-30, R4-CT-7(a); see `reports/CORRECTIONS.md`).

**Arithmetic operators keep the parse tree.** An AOR or LOG replacement wraps multi-token operands in parentheses, so replacing $+$ by $*$ in $a+b*c$ gives $a*(b*c)$, not $a*b*c$.

| | Count |
|---|---|
| Single mutants | 818 = 166 LIT + 12 NEG + 93 LSW + 4 BLIT + 58 REL + 7 QNT + 36 LOG + 33 IF + 7 MATCH + 125 AOR + 180 VAR + 45 ENDSWAP + 22 LIST + 22 CAT + 8 CONS |
| … predicted to fail to elaborate | 47: 14 negative Nat literals, 26 swaps of `.1/.2` on an `Int × Assign` pair, 7 match swaps that move a pattern variable out of its arm |
| Pairs with disjoint sites | 331 490 of 334 153 |
| … containing a single that fails to elaborate | 37 105 |
| Space | 818 + 331 490 = 332 308 |

All singles are evaluated, and all pairs; triples are not. MATCH contributes no well-defined mutant (all 7 arm swaps fail to elaborate, §3), so 14 of the 15 operators populate the evaluated space (corrected after review 2026-10-30, R4-CT-7(c); see `reports/CORRECTIONS.md`).

## 3. Evaluation and validation

**Transcription.** `chainalpha_transcription.py` transcribes each of the 44 definitions and the Part 0 utilities with Lean's semantics:
- `getD` defaults;
- truncated Nat subtraction;
- $x/0=0$ and $x\,\%\,0=x$ on Nat;
- Euclidean `/` and `%` on Int, checked with Lean (`mutation_spotcheck.out`): $(-7)/2=-4$, $7/(-2)=-3$, $(-7)/(-2)=4$, $(-7)\,\%\,2=1$;
- truncating `zip`.

Each site has one hook, named `<definition>:<operator>:<k>`. The script asserts that the hooks named in the transcription are exactly the 418 sites other than the 7 MATCH sites (out §2). Every match swap moves a pattern variable out of its arm, so those 7 fail to elaborate and need no hook.

Readouts are compared exactly. $\sqrt k$ is reduced to $s\sqrt m$ with $m$ squarefree, and distinct squarefree radicals are linearly independent over $\mathbb Q$. A mutant is **well-defined** when it elaborates, its readout dimension is nonzero and its numerator is nonzero.

**Checks** (out §2; `mutation_spotcheck.out`):
- The unmutated readout equals the statement of `dressing_of_the_qutrit_qubit_link`.
- With the author's 19 rewordings written as hooks, the transcription prints the 20 lines of `upstream_mutations_lean.out` character for character, including whether `main` and `commutators` are still proved. Five rewordings are elements of the space: `reads_at_qubit` is the pair `readEnd:IF:0=swap & farEnd:IF:0=swap`; `X_at_far_end`, `exch_shift1`, `drop_fixed` and `drop_even` are singles. The other 14 change several tokens or add code.
- Lean 4.35.0 was run on every single mutant: Parts 0–1 with the edits, then `#eval readout objectWorld objectTwoPoint`.
  - All 770 readouts that Python computes are identical in Lean.
  - The 47 mutants predicted to fail to elaborate all fail, and no other mutant does.
  - The one mutant beyond the transcription's resource bound (`addRoot:CONS:2=dup0`, which doubles the readout list at every record, about $2^{172}$ terms) makes Lean abort with "deep recursion".
- A seeded sample of 300 pairs (seed 20261030) agrees as well: 266 readouts identical, 31 predicted elaboration failures confirmed, and 3 pairs beyond the resource bound, on which Lean, like on the single, aborts with "deep recursion" (4 such aborts in all). No sampled pair that Python leaves unevaluated is computed by Lean (`mutation_spotcheck.out`).
- The same output prints Lean's integer conventions that the transcription follows.

**Outcome** (out §3):

| | Singles | Pairs |
|---|---|---|
| Fail to elaborate | 47 | 37 105 |
| Readout list > 4000 terms (not evaluated) | 1 | 715 |
| Readout dimension 0 (value undefined) | 24 | 17 924 |
| Numerator 0 ($1/\text{value}$ undefined) | 4 | 3 500 |
| Well-defined | 742 | 272 246 |

All 715 pairs beyond the bound contain `addRoot:CONS:2=dup0` (out §3). The longest readout list of an evaluated mutant has 2049 terms.

## 4. Results

The same metric code is applied to the mechanical space and to B, rebuilt from `mutations_classified.py` (out §4):

| | Mechanical singles | Mechanical singles + pairs | B (F-CT-4) |
|---|---|---|---|
| Mutants / rules | 818 | 332 308 | 46 080 |
| Well-defined | 742 | 272 988 | 46 080 |
| … with the Letter's value | 244 | 29 292 | 1 |
| Distinct values | 163 | 10 340 | 15 639 |
| Distinct values with $1/\text{value}\in(120,155)$ | 38 | 1105 | 3638 |
| Density per unit there | 1.09 | 31.6 | 103.9 |
| Density per unit in $(136,138)$ | 2.0 | 35.0 | 105.5 |
| Reference count within $\pm1\sigma$, $U$ uniform on $(120,155)$ (illustrative model) | $4.6\times10^{-8}$ | $1.3\times10^{-6}$ | $4.4\times10^{-6}$ |
| Distinct values within relative $10^{-2},10^{-3},10^{-4},10^{-6}$ | 6, 1, 1, 1 | 91, 9, 1, 1 | 288, 28, 4, 1 |
| Median gap in $(120,155)$ | 0.60 | 0.022 | 0.0070 |

Here $\sigma$ is the CODATA 2022 uncertainty, $2.1\times10^{-8}$ in $\alpha^{-1}$, and the window $\pm1\sigma$ is $4.2\times10^{-8}$ wide. The reference count is computed as in `freedom.py`: the total length of the distinct values' $\pm\sigma$ intervals, clipped to the band, divided by 35. It is a count in a model in which the target, not the rule, is random; it is not a calibrated p-value for the measured $\alpha$, and it is not a search-history estimate. The row "mutants for an $O(1)$ chance hit" of the first version is removed: it divided a density threshold by the in-band yield per mutant, which mixes the density with how often the family repeats a value (corrected after review 2026-10-30, R4-CT-6; see `reports/CORRECTIONS.md`).

The density is nearly flat across $(120,155)$, $(130,145)$ and $(136,138)$ (31.6, 32.9 and 35.0 per unit), so the reference count changes little with the band (out §4: $1.3$, $1.4$ and $1.5\times10^{-6}$ for $U$ uniform on each). This describes the spacing of the values at the scale of these bands; it does not establish uniformity at the $4.2\times10^{-8}$ scale of the window (F-CT-4, Limits).

**Depth dependence** (corrected after review 2026-10-30, R4-CT-1; see `reports/CORRECTIONS.md`). The family is truncated at two edits, and the density of distinct values in $(120,155)$ is set mainly by that truncation:

| depth | distinct values in $(120,155)$ | density per unit | source |
|---|---|---|---|
| 1 (singles) | 38 | 1.09 | out §4 |
| ≤ 2 (singles + pairs) | 1105 | 31.6 | out §4 |
| 3 (triples, estimated) | ≈ 2000 (low estimate) | ≈ 58 | `xcheck_mechanical.out` [6] |

The depth-3 row comes from a seeded sample of 12 000 triples with pairwise-disjoint sites: 389 distinct in-band values, a Chao1 estimate of 1139 for the whole depth-3 space, scaled by the depth-2 ratio of the true count to the Chao1 estimate from a 12 000-pair sample (1105/616). Chao1 underestimates more as the sampling fraction falls, and the depth-3 fraction is about 270 times smaller, so the estimate is low. B is a complete product over its 8 parameters and has no comparable truncation. So the statement "the density is lower than B's" is a statement at depth ≤ 2, and nothing here bounds depth 4 or beyond. In the illustrative model the depth-3 estimate corresponds to a reference count of about $2.4\times10^{-6}$.

**Closest distinct values** (singles + pairs; out §4):

| Deviation | $1/\text{value}$ | Value | Example |
|---|---|---|---|
| $0.598\sigma$ | 137.035 999 190 | $(79+48\sqrt2+42\sqrt3+3\sqrt6)/31104$ (the Letter's) | 29 292 mutants, e.g. `pairKet:IF:0=swap` |
| $9.9\times10^5\sigma$ | 137.056 872 | $(149+45\sqrt3)/31104$ | 20 pairs, e.g. `recs:REL:1=strict & linkSwap:VAR:0=read2` |
| $1.26\times10^6\sigma$ | 137.062 470 | $(41+24\sqrt2+18\sqrt3+3\sqrt6)/15552$ | 18 pairs |
| $2.9\times10^6\sigma$ | 137.097 868 | $(79+47\sqrt2+47\sqrt3)/31104$ | 1 pair |
| $3.0\times10^6\sigma$ | 137.099 922 | $(82+51\sqrt2+42\sqrt3)/31104$ | 3 pairs |

Mutants that leave the value unchanged are common: 244 of the 742 well-defined singles. By operator (out §4) these include all 7 QNT mutants, 30 of 38 ENDSWAP mutants and 68 of 161 VAR mutants. They do not add values; the metrics count distinct values.

## 5. Comparison with B, Table S5 and the author's rewordings

- **Overlap with B.** The mechanical space and B share 333 distinct values, 100 of them in $(120,155)$, against 1105 mechanical and 3638 B values there. The two spaces are largely different. Their union has density 132.7 per unit in $(120,155)$; its reference count in the illustrative model is $5.6\times10^{-6}$.
- **The author's 19 rewordings** (`upstream_mutations_lean.out`). Here "reached" means that some single or pair has the same value; membership is by value, not by mechanism (corrected after review 2026-10-30, R4-CT-4; see `reports/CORRECTIONS.md`).
  - Two rewordings, `fixesFar_always_true` and `wX_dA_equiv`, leave the value at the Letter's, which 29 292 mutants share.
  - Of the 17 that change the value, 16 values are reached, often by a different edit. For example, `Y_at_center` (127.536) is reached by `wY:AOR:0=+`, which makes `wY` equal to `at_ 0 + at_ 1` instead of `at_ 0 * at_ 1`. Not reached: `X_where_letters_sit` (139.030), which adds a new conditional to `kOf`.
  - Some reached values arise by a different mechanism. The Born value is a count of 319 records at weight 1, not a sum of squared weights (below). Table S5's "smallest readable configuration", $(30+30\sqrt2)/256$, is reached by making the centre a qubit (`objectWorld:LIT:2=2`) and dividing instead of multiplying the unread link's traceless factors (`readoutDim:AOR:5=/`), so that the unread link contributes 1.
- **Table S5.** Of its 15 closed-form alternatives, the values of 12 are reached, some by an unrelated rule: 9 by singles and 3 by pairs ("the two ends swapped", "Born weights $k_r$", "smallest readable configuration"). Not reached:
  - "exchange of reads where the letters sit" (139.030);
  - "count of read ends" (70.190);
  - "normalised by the 172 records" (0.758).
- **The Born weights** $k_r$. Their value $319/31104$ is reached by 7 mutants, all of them pairs and all containing `recs:REL:0=neg` (out §4 lists them; e.g. `recs:REL:0=neg & sameForm:REL:0=neg`). The single `recs:REL:0=neg` negates the filter `exchAssign p == p`, so the records built from the 11 exchange-fixed products are replaced by one record for each of the 158 products the exchange moves; it has $172-11+158=319$ records and the value $(237+48\sqrt2+31\sqrt3+3\sqrt6)/31104$. Each of its seven partners alone gives every record weight 1 (value $172/31104$, Table S5's "no weights"). In each pair the numerator is therefore the record count 319 at weight 1, which equals $79+2\cdot48+3\cdot42+6\cdot3$ by coincidence (corrected after review 2026-10-30, R4-CT-3; see `reports/CORRECTIONS.md`). None of the declared operators squares $k$: the square root is how the radicand list is read, not a definition in the file.

## 6. What the mechanical space does not contain

- **Multi-token rewordings.** Only single and paired token-level edits are generated. Rules that need new code, such as the end assignment "where the letters sit" or Table S5's count of read ends, are absent; so is anything that needs three or more edits. Triples ($\binom{818}{3}\approx9.1\times10^7$ before the disjointness rule) are not evaluated, and the density depends on this truncation (§4).
- **Most of B.** 15 306 of B's 15 639 distinct values are not mechanical values. B varies the per-class factor $k\in\{1,2,3,6\}$ independently on four classes, the commutator factor and count, and five normalisations. The Lean text reaches these factors only through its own structure.
- **The amplitude/Born choice and normalisation by the record count** (CT-G05, CT-G08). The $\sqrt k$ reading and $N=172$ are not definitions in the file. B contains the Born reading (it varies amplitude against Born weighting); B and the mechanical family both vary $N$, B over five values and the mechanical family through edits of `readoutDim`, a definition in the file.
- **Choices outside the readout.**
  - Selection: the read cap, the acting-charge convention of F-CT-3 and `configuration_unique` sit in Part 1.9 and Part 5 and do not enter the readout. Other configurations arise mainly through edits of `objectWorld`/`objectTwoPoint`, for example the centre dimension $3\to2$ or $3\to4$.
  - The Part 0 arithmetic is held fixed.
- **The construction history** (CT-G10). Whatever was tried before the file was written is not visible in it.

## 7. Limits

- **The operator set and label pairs are ours.** They are standard mutation-testing operators (constant replacement, relational, logical and arithmetic operator replacement, operand deletion, branch swap), with variable swap restricted to declared same-role pairs. A different set gives a different density. Here, at depth ≤ 2, the density is lower than B's. The mechanical space has fewer distinct values (10 340 against 15 639), and a smaller share of them lies in $(120,155)$: 11%, against 23%. The density grows with the depth (§4) (corrected after review 2026-10-30, R4-CT-1; see `reports/CORRECTIONS.md`).
- **Pairs are evaluated by the transcription, not by Lean.** Lean checks every single and 300 pairs. Ill-typedness of a pair is predicted as "one of its singles fails to elaborate". This held for the 31 such pairs in the sample, and no other sampled pair failed.
- **The resource bound** excludes 715 pairs (0.2%) whose readout list would exceed 4000 terms. All of them contain the list-doubling mutant, and Lean also aborts on the 3 that fall in the sample. They are counted separately and do not enter the density.
- **Sensitivity only.** Like F-CT-4 B, this describes how the values of a declared family of variants lie near $\alpha^{-1}$. The reference counts come from an artificial model in which the target, not the rule, is random. They are not probabilities for the construction, and neither family identifies the variants actually considered or the order of rule choices and evaluation (CT-G10).
