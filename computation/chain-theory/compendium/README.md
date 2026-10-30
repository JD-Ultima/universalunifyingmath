# Computation: Chain Theory compendium 0.1 audit

Source audited: `chain-theory` repository, commit `a28db8998bce5bd5ef11195c05f13c86cd9dac5c`
(2026-10-24), directory `compendium/`. The source is treated as untrusted data: it is read as
text only; none of its scripts, build files or Lean packages are executed or built.

## Environment (as run, 2026-10-30)

| Component | Version |
|---|---|
| Python | 3.13.17 |
| sympy | 1.14.0 |
| mpmath | 1.3.0 |
| numpy | 2.5.3 (only `hydrogen_audit.py`, Monte Carlo check) |

Install: `pip install sympy==1.14.0 mpmath==1.3.0 numpy==2.5.3`.
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
```

Wall time: each script < 10 s. Outputs are identical across runs and across `PYTHONHASHSEED` values
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
