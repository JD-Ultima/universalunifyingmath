# Computation: Chain Theory compendium 0.1 audit

Source audited: `chain-theory` repository, commit `a28db8998bce5bd5ef11195c05f13c86cd9dac5c`
(2026-10-24), directory `compendium/`. The source is treated as untrusted data: it is read as
text only; none of its scripts, build files or Lean packages are executed or built. Exception (round 6):
`lean_build.py` copies the Lean *library* sources into a fresh scratch directory, refuses any module that runs code at
elaboration time (`#eval`, `run_cmd`, `IO.`, `import Lean`, macros, `extern`, `unsafe`, `native_decide`; none found),
and compiles the copies with `lean -o`; no lakefile, lake, script or `AxiomsCheck` file of the source is run.

## Environment (as run, 2026-10-30)

| Component | Version |
|---|---|
| Python | 3.13.17 |
| sympy | 1.14.0 |
| mpmath | 1.3.0 |
| numpy | 2.5.3 (only `hydrogen_audit.py`, Monte Carlo check) |
| scipy | 1.18.1 (round 2 `scorecard*.py`: CODATA 2022 via `scipy.constants`) |
| Lean | 4.35.0, core only (round 6 `lean_audit.sh`; the source pins 4.32.1) |

Install: `pip install sympy==1.14.0 mpmath==1.3.0 numpy==2.5.3 scipy==1.18.1`. Round 2 optionally checks its PDG 2024
transcription against `mass_width_2024.mcd` from `pip download particle==0.25.2` (data file read as text; the package is not imported).
All scripts run with `python3 -I` (isolated mode). No network access is needed.

## Run

```bash
SRC=/path/to/chain-theory/compendium      # checkout at a28db899
cd computation/chain-theory/compendium
python3 -I census.py        $SRC/index/statements.tsv            > outputs/census.txt
python3 -I depgraph.py      $SRC --json outputs/depgraph.json    > outputs/depgraph.txt
python3 -I depgraph.py      $SRC --no-chapter-form               > outputs/depgraph_no_chapter_form.txt
python3 -I lean_classify.py $SRC 141.1 143.1 64.1 64.2 65.1      > outputs/lean_classify.txt
python3 -I ae_audit.py                                           > outputs/ae_audit.txt
python3 -I hydrogen_audit.py                                     > outputs/hydrogen_audit.txt
python3 -I gauge_audit.py                                        > outputs/gauge_audit.txt
python3 -I register_checks.py                                    > outputs/register_checks.txt
python3 -I boxes.py         $SRC 141.1 143.1 64.2               # prints box text (helper)
# round 2: numeric scorecard
python3 -I scorecard_extract.py  $SRC --tsv outputs/scorecard_candidates.tsv > outputs/scorecard_extract.txt
python3 -I scorecard.py --src $SRC --pdg-mcd $MCD --tsv outputs/scorecard.tsv > outputs/scorecard.txt
python3 -I scorecard_freedom.py                                       > outputs/scorecard_freedom.txt
# round 3: mass cells
python3 -I masscells.py $SRC --pdg-mcd $MCD --tsv outputs/masscells_readings.tsv  > outputs/masscells.txt
# round 4: statement-level premises (about 6 s)
python3 -I premises.py $SRC --tsv outputs/premises_items.tsv > outputs/premises.txt
# round 5: hydrogen body registers, the ladder beyond n = 8 (about 3 s; uses scipy.constants and masscells.parse_figure)
python3 -I hydrogen_registers.py $SRC > outputs/hydrogen_registers.txt
# round 6: Lean at scale (build with Lean 4.35.0 in WORK, #print axioms, own checks, resolution; see below for time)
LEAN=/path/to/lean-4.35.0/bin/lean SRC=$SRC WORK=/fresh/scratch/dir ./lean_audit.sh
python3 -I lean_resolve.py $SRC          # static part alone (about 5 s), without build, axioms or checks
# MCD=.../particle/data/mass_width_2024.mcd (optional; without it the PDG transcription is not re-checked)
```

Wall time: each script < 10 s (round 6: `lean_resolve.py` about 3 s; `lean_audit.sh` about 45 min, almost all of it
Build B: 2494 s wall at 3 jobs on 4 shared cores, 5398 s summed over modules, peak 7677 MB; the axioms file 5 s). Outputs are identical across runs and across `PYTHONHASHSEED` values
(checked). `hydrogen_audit.py` uses a fixed seed (20261030); its Monte Carlo
eigenvalues are reproducible to the printed 3 decimals.

## Files

| Script | Task | What it computes |
|---|---|---|
| `census.py` | 1 | counts by grade x Part from `statements.tsv`; checks the headline 100/109/198/64 (and README's 56 tables, 7 definitions) |
| `depgraph.py` | 2 | parses `premises`/`used_by` into a directed graph; reciprocity, cycles, chain-rule audit, load-bearing unproven premises, dependents of 124.1 (alpha) |
| `lean_classify.py` | 4 | classifies every Lean name in the `lean` column (opaque Prop / open-result record / string-table check / restated theorem / kernel module) |
| `boxes.py` | 1,3 | extracts boxed text by statement number, with file:line |
| `ae_audit.py` | 3(i) | QED C1, C2, C3 (exact/40-digit); the numbers chapter 7-2 prints; literal-identification test; experimental weight of C3 |
| `hydrogen_audit.py` | 3(ii) | exact check of the printed level-one binding; comparison with alpha^2 mu/2; reconstruction of the 1/n^2 ladder and its dependence on the support dimension |
| `gauge_audit.py` | 3(iii) | for all six genres: dimension, derived algebra, centre action, kernels; Z_6 cover for (3,2); the u(1) representative; monomial group orders; Z_6 neutrality of SM fields; sensitivity of the genre selection to the read cap and alphabet |
| `register_checks.py` | 1 | weak angle 1/2 vs (M_W/M_Z)^2; form of the printed running; Koide ratio |
| `scorecard_extract.py` | R2-1 | segments the tex into 468 statement units (index numbering, grades checked), screens them with a 42-entry physical-quantity lexicon (242 candidates, 174 with numerals), verifies the 210 verbatim anchors of the curated rows (numbers, identification sentences, verdicts), assigns every non-covered candidate a reason code |
| `scorecard_data.py` | R2-1 | curated data (module): 53 named quantities, 124 printed faces, verbatim anchors, identification / derivation codes, the source's verdicts, exclusion codes |
| `scorecard.py` | R2-2,3 | exact evaluation of every face, conversion to the measured quantity, deviation (sigma, relative, factor), bins, summary table; block REDERIVE re-derives 784/1568, 3/8, sqrt3/18, the CHSH value, the lepton and proton cells, the Koide values (17/36 is re-derived in `scorecard_freedom.py`) |
| `scorecard_freedom.py` | R2-3 | F-CT-4-style freedom tests for the agreements besides alpha(0): weak angle (count-face family), lepton cells, proton cell |
| `masscells_data.py` | R3 | curated data (module): 102 verbatim anchors for the trace, the printed closures and G2, the printed proton/electron/muon/tau cells, exclusions with reasons, web descriptions |
| `masscells.py` | R3-1..4 | mass cells: (2,3) record census (319, 79 = 48+28+3); electron web parsed from the 152.1 figure, muon/tau/proton webs built and checked; measuring-read paths; 72-reading grid validated against 6 printed tuples; 155.1 faces; printed mass-law variants; look-elsewhere in printed / grid / declared / narrow / census / integer spaces; tau/mu; 181.5 cross-check; sweep audit |
| `premises.py` | R4 | statement-level premises from the derivations' premise sentences: segmentation into derivation groups, 1069 premise items with verbatim quotes, candidates by grade word under two designator readings (module numbering / `designators.tsv`), a separate judgment layer (21 groups), chain-rule re-grade, cycles (any-of, strict, forced under every resolution), F4 upper/lower bounds; writes `outputs/premises_items.tsv` |
| `hydrogen_registers.py` | R5 | 143.1 body registers: per-link faces (16.1, 103.1) against the 16 printed genre faces; the product rule at the W/Z body and gluon word (port-only cut against every-link cut); electron web parsed from the 152.1 figure, whole register 2^30 3^18 and cut register 2^31 3^17 under 9 register readings; mass face 36 against the printed link counts; E1 at every unit and register reading; mechanical scan of every ladder sentence (range n <= 8); reconstruction R exactly at n <= 200 (sl2 commutant), support dimension 2/3/4/6, order and per-record readings; E1 against R_inf hc / m_e c^2 and IE(H) (CODATA 2022, NIST), a 36-variant declared conversion family, and E_2/E_1 against the measured 1S-2S centroid (Parthey et al. 2011) |
| `lean_resolve.py` | R6 | resolves the 16 "restated/kernel" THEOREM rows of `lean_classify.py` (imported unchanged): verbatim type of each of the 8 declarations, its 26 `AlphaDerivation.Statements` conjuncts and their Proofs-level theorems (proof method, binders, opaque hypotheses and whether used, entered `Expected` constants), content kind per conjunct, the entered step table, restatement counts per Statements group; the 3 kernel modules (rows computed or entered, what `checks` is decided over, citing rows); clause-level verdicts against the box text with printed numbers looked up in the tied statements and the whole tree; the box each declaration's content states; entered vs computed `rows` in all 150 Row modules; build facts (pins, imports, patched modules and their dependents); `#print axioms` summary and comparison with the source's 4.32.1 report; own checks; the re-stated 64-row table (`outputs/lean_resolve.tsv`); `--emit-axioms-file` writes our `#print axioms` file |
| `lean_resolve_data.py` | R6 | curated data (module): per conjunct a content kind, note and verbatim Lean anchors (44); per row the box clauses (52 verbatim quotes) with their counterpart conjuncts; the box each declaration's content states, with anchors |
| `lean_build.py` | R6 | builds the six packages' library modules (import closure of the targets) in dependency order with `lean -o`, `LEAN_PATH` in a fresh scratch directory; safety gate; optional `--patch` (TSV of single-line edits to the copies, recorded); critical-path scheduling with a memory budget; per-module time, peak RSS, exit, warnings, `sorry` warnings, sha256 of the original; `--resume` |
| `lean_build_patch.tsv` | R6 | the five compatibility edits of Build B: `:= by decide` -> `:= rfl` / `⟨rfl, rfl, rfl⟩` at AxiomDeclaration/Theorems.lean:42, :52 and AxiomDeclarationTwo/Theorems.lean:42, :52, :56 (statements unchanged; the theorems are referenced by no other library declaration) |
| `lean_audit.sh` | R6 | the round-6 pipeline: Build A (published sources, AxiomDeclaration library), Build B (all 493 modules, patched), `Diagnose.lean`, our `#print axioms` file, `ResolveChecks.lean`, `lean_resolve.py` |
| `lean_audit/Diagnose.lean` | R6 | why the five proofs fail under 4.35.0: `fail_if_success decide` then `rfl` on the five statements verbatim; `#synth` of the instance at the stated type |
| `lean_audit/ResolveChecks.lean` | R6 | own theorems about the source's definitions: `generatorBlocksScalar` always true; abelian = conserved = (d1 != d2) on 1..12; the record-wise weight identity over {+-1}^2; `onShellCoupling` = (0,1) off the object of record; the correspondence conclusion without its opaque hypotheses; `identityAt` of the unit functional on every path; six literal parts by `rfl`; the entered zero transfer in 172 rows |

## Interpretive choices in the dependency graph (`depgraph.py`)

- C1: every statement number in a premise string is a premise; "a to b" expands within a section.
  Kind words are ignored because the source labels CONDITIONAL statements as "Theorem" in these lists
  (e.g. "Theorem 124.1", "Theorems 18.1 to 18.4" where 18.1, 18.3, 18.4 are CONDITIONAL).
- C2: appendix sections lettered A..W in the index are numbered 238..260 in the premise text
  (letter L <-> 237 + position). Checked by the grade pattern of D <-> 241 and H <-> 245.
- C3: the 23 premise strings of the form "the theorems of Chapters X-Y ..." are resolved through
  `designators.tsv` (first_section..last_section, else the module-numbering range), taking THEOREM
  and CONDITIONAL statements. 82 of 162 designators have two readings; the sensitivity run
  `--no-chapter-form` drops these strings.
- C4: edge p -> s means p appears in the premise list of s. Premise lists are chapter-granular
  (they cite every statement of a cited chapter), so an edge records a declared dependency, not a
  verified logical use.

## Recorded outputs (key lines; full text in `outputs/`)

```
census:   THEOREM 100, CONDITIONAL 109, OPEN 198, THEOREM tied to Lean 64  -> all MATCH
          TABLE: headline 56, computed 54 -> MISMATCH
          Lean tie is per section: 0 of 284 sections carry differing Lean lists across their statements
depgraph: edges 8367; reciprocity 3632/8367 (43.4%)
          one strongly connected component of 22 statements incl. 124.1; 2-cycles 124.1<->16.1, 124.1<->58.1
          (vanishes with --no-chapter-form)
          THEOREMs listing >=1 OPEN premise directly: 66 of 100 (66 also with --no-chapter-form)
          THEOREMs listing >=1 OPEN or CONDITIONAL premise directly: 84 (78 with --no-chapter-form)
          most load-bearing unproven premises: 2.2 (C), 2.4, 2.5, 2.7, 2.8 (O) -> 183 of 209 T+C statements
          transitively; 9.1 (C), 9.7 (C), 9.8 (O) -> 175; 18.1/18.3/18.4 (C), 18.5/18.6 (O) -> 152
          124.1 (alpha): 22 direct users, 70 transitive dependents (32 with --no-chapter-form)
lean:     of 64 Lean-tied THEOREM rows: 11 cite only opaque Props/open-result records; 37 have at best
          string-table checks; 16 cite restated AlphaDerivation theorems or a kernel module
          137 of 150 data/kernel modules have Row structures whose fields are all String
ae:       C1 = 0.5, C2 = -0.328478965579194, C3 = 1.1812414565872 (asserted against reference values)
          chapter prints T0=1, T1=309, T2=309^2+closure (closure not printed), T3=0, G2(2,2)=5/216(1+sqrt2)
          T3 = 0 vs C3 != 0: fails for every normalisation; C3 (alpha/pi)^3 = 1.48042e-8 = 9.26e4 sigma
          (Rb alpha), 5.63e4 sigma (Cs alpha)
hydrogen: E1_whole = 67/(2^28 3^15) = 1.739467e-14, E1_cut = 67/(2^29 3^15) = 8.697333e-15 (ratio 2)
          vs alpha^2 mu/2 = 2.661118e-5: model/physical = 6.54e-10 (whole), 3.27e-10 (cut)
          1/n^2 for n=1..8 reproduced iff the symmetric power is over a 2-dim support (m=2)
gauge:    dims 7,12,19,17,24,31 and |central scalar| 0,1,2,0,1,0 reproduce the genre matrix
          (3,2) unique under (dim 12, rank 4) and under (d1 d2 <= 6 and d1 != d2); cap 8 selects (3,2),(4,2)
          ker(SU2 x SU3 x U1 -> S(U2 x U3)) = Z_6; diag(3,3,-2,-2,-2) = 6 Y on the 5; |mu4 wr S_n| = 4,32,384,6144
weak:     printed cos^2(theta_W) = 1/2 vs (M_W/M_Z)^2 = 0.77679(26): 1067 sigma
scorecard (round 2):
          screen 242 candidates / 174 with numerals; 60 carry a row, 114 reason-coded; 210 anchors verbatim (125 numbers, 53 identification sentences, 32 verdicts), 0 missing
          53 named quantities, 124 printed numbers; identification E 9, N 43, P 1
          34 numerical rows by best face: <=3 sigma 2, <=1% 5, <=x10 9, disagree 18 (adopted face, 16 rows: 1/1/3/11)
          <=3 sigma: 1/alpha(0) +0.60 sigma; sin^2 theta_W = (17/36)^2 = 289/1296 (164.2 O) -0.83 sigma on-shell (-3.6% MS-bar)
          <=1%: (963/67)^2 = 206.587 (-0.088%), (3969/67)^2 = 3509.24 (+0.92%), Koide 0.667660 (+0.15%),
                (2871/67)^2 = 1836.18646 (+1.84e-5), m_n/m_p = 1 (-0.14%)
          5 of the 7 agreements rest on a choice among printed faces; the release's 'comparison table' (25 mentions) is absent
freedom:  17/36 re-derived (17 of 36 Weyl pair letters nonzero at the (3,2) pair ket; Gell-Mann/Pauli: 5)
          weak angle: 0.019 (printed faces) / 0.14 (192-variant count-face family) expected chance hits
          leptons: 0.007 (mu), 0.07 (tau); 1 of 504 combinations has both within 1%
          proton: 1 of 9 printed cells within 1%; 0.0025 in a 1320-variant space; nearest-integer 0.053 / 0.006
masscells (round 3):
          anchors 102/102; 319 = 79+96+126+18 re-derived (172 records; 19 of 4608 alpha weight rules give 319)
          webs: e (9, 10, 79, 12, 67, 3), mu (13, 15, 117, 10, 107, 9), tau (17, 20, 155, 8, 147, 27), proton seed (10, 12, 101, 18), sink (8, 9, 81, 18) = printed
          grid 72 readings, 6 printed tuples reproduced; all four cells at 1 (186.1), at 0 without its strike clause; ratios within 1%:
          mu+tau 12, p 8, all three 1; faces: dock = hub -> 22.95, 389.9, 1836.19; mixed -> 206.59, 3509.2, 16525.7
          m_p/m_e by mass law: S^2 +1.84e-5, 195.1 office +6.20e-6, 176.2 G2 at e,mu only -1.01e-4 (spread 1.2e-4)
          lambda(p): printed 2.9e-4, grid 2.6e-4, declared 2.5e-3, narrow 1.7e-5..5e-5, any integer 0.053 (0.34 at the spread)
          tau/mu = (441/107)^2 = 16.987, +1.0% (198 sigma); 181.5 top sink 2871 = proton cell
premises (round 4):
          297 derivations, 294 groups (61 shared by several boxes), 261 with a premise sentence; 1069 items, 1069 quotes verbatim;
          0 statement-number citations inside derivations; designator readings differ at 171 of 397 items (title overlap MOD 129, FIF 3)
          index two-cycles 124.1<->16.1, 124.1<->58.1: no edge cited in any derivation
          forced cycle (every resolution of any-of items, MOD): {35.1, 40.1, 124.1} via "the theorem of chapter 18" and
          "the explicit two-point function with its records (chapter 5-1)"; strict graph with judgment: one SCC {124.1, 18.2, 18.3}
          14 requested THEOREMs: index OPEN 10, OPEN/COND 12; derivations (J) declared OPEN 8, violate 10, undetermined 4, consistent 0
          all 100 THEOREMs: declared OPEN 41, other violations 35, undetermined 19, consistent 3, no premise sentence 2 (rule fixed after REVIEW3-CT R3-CT-8; was 38 / 35 / 22)
          F4 lower bounds (T+C): 2.2 -> 1; 2.4/2.5/2.7/2.8 -> 153; 9.1/9.7 -> 151 (J); 9.8 -> 0; 18.x -> 0; 124.1 -> 152 (313 all grades)
lean (round 6):
          16 rows -> 8 declarations (7 restated AlphaDerivation, 1 kernel); 26 conjuncts: COMPUTED 13, +LITERAL 6, GENERAL 5,
          DEFINITIONAL 1, ENTERED 1; decide +kernel 18, constructor 4, tactic 4; Lean anchors 44/44, box anchors 52/52
          5 opaque hypotheses (one conjunct), 0 used; Expected.centreDims used by 0 theorems; zero transfer entered 172/172
          verdicts: MATCH 1 (9.4), PARTIAL 5 (2.1 9.2 9.3 18.2 196.2), NONE 10; 56.1's 322714, 319332 nowhere in the tree
          content of a CONDITIONAL box of the section: ORKC->2.2, GWDC->18.3, MDFN->74.2, WDCL->78.1, WSWC->78.2
          Row modules: 147 data modules entered (135 String-only), 3 kernel modules computed
          Build A (published, 4.35.0): AxiomDeclaration.Theorems:42,:52 fail (Decidable not synthesized); 277 of 493 import it
          Build B (5 proofs -> rfl): 493/493, 0 errors, 0 warnings, 0 sorry; 5398 s summed, 2494 s wall, peak 7677 MB
          #print axioms (6662 = source's list): sorryAx 0; none 5965, propext 290, propext+Quot.sound 4, +Classical.choice 403
          vs source 4.32.1 report: 6259 equal, 403 gain Classical.choice (via core Int/Fin.instLinearOrderPackage)
          coverage by content: 6 / 46 / 12 (by kind unchanged 16 / 37 / 11)
```

## External values used (with sources)

- a_e = 1.15965218059(13) x 10^-3: Fan, Myers, Sukra, Gabrielse, PRL 130, 071801 (2023).
- alpha^-1 = 137.035999206(11) (Rb): Morel et al., Nature 588, 61 (2020); 137.035999046(27) (Cs):
  Parker et al., Science 360, 191 (2018); 137.035999177(21): CODATA 2022.
- C2: Petermann, Helv. Phys. Acta 30, 407 (1957); Sommerfield, Phys. Rev. 107, 328 (1957).
  C3: Laporta & Remiddi, Phys. Lett. B 379, 283 (1996).
- m_p/m_e = 1836.152673426(32), m_e c^2 = 510998.95069 eV: CODATA 2022.
  H ionization energy 13.598434599702 eV: NIST Atomic Spectra Database (H I).
- Lamb shift 1057.845(9) MHz: Lundeen & Pipkin, PRL 46, 232 (1981). H hyperfine 1420.405751768 MHz:
  Hellwig et al., IEEE Trans. Instrum. Meas. 19, 200 (1970).
- M_W = 80.3692(133) GeV, M_Z = 91.1880(20) GeV, sin^2 theta_W(MS-bar, M_Z) = 0.23122, lepton masses: PDG 2024.
- Global form of the SM gauge group: Tong, JHEP 07 (2017) 104.
- Round 2 (`scorecard.py`): CODATA 2022 from `scipy.constants` 1.18.1 (alpha, m_mu/m_e = 206.7682827(46), m_tau/m_e =
  3477.23(23), m_p/m_e, m_n/m_p, G_F, sin^2 theta_W = 0.22305(23), masses, G); PDG 2024 masses and widths from
  `mass_width_2024.mcd` (RPP 2024; M_W 80.3692(133) GeV summary value); transcribed, not machine-checked:
  alpha_s(M_Z) = 0.1180(9), alpha^-1(M_Z) = 127.951(9), |V_ud| = 0.97367(32), J = 3.08(14)e-5, sin^2 theta_13 =
  0.0220(7), N_nu = 2.9963(74), Planck 2018 Omega_c h^2 = 0.1200(12), Omega_b h^2 = 0.02237(15), Omega_Lambda =
  0.6847(73), H0 = 67.4(5); CGLMP: Collins, Gisin, Linden, Massar, Popescu, PRL 88, 040404 (2002); sin^2 = 3/8:
  Georgi, Quinn, Weinberg, PRL 33, 451 (1974).
- Round 3 (`masscells.py`): CODATA 2022 m_p/m_e, m_mu/m_e from `scipy.constants` 1.18.1; PDG 2024 m_e, m_mu, m_tau and the quark
  masses of the 181.5 cross-check from `mass_width_2024.mcd` (checked with --pdg-mcd): m_tau/m_e = 3477.365(176), m_tau/m_mu = 16.81769(85).
