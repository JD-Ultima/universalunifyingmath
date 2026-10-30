# Chain Theory α audit — computations

Environment: Python 3.13.17, numpy 2.5.3, sympy 1.14.0, mpmath 1.3.0 (`pip install numpy==2.5.3 sympy==1.14.0 mpmath==1.3.0`). Randomness only with fixed seeds (`readcap_lemmas.py`: 20261030; `symmetry.py`: 7; `operator_origin.py`: 20261030; `transport.py`: 20261031; `mechanical_mutations.py`: 20261030, for the sample of pairs checked in Lean).

The complete reproduction loop below uses Bash and was run on Linux. `transport.py` requires the `fork` multiprocessing start method: its workers inherit the operator closures and large arrays built by the parent. On Windows, run it in WSL using Linux Python; it now exits immediately with that requirement on native Windows. The other Python scripts can run natively with `python -X utf8 script.py` so Unicode output uses UTF-8.

| Script | Purpose | Runtime | Output |
|---|---|---|---|
| `records_independent.py` | Operator-level rebuild of the 172 records; sectors; exact value; 60 digits (F-CT-1) | ~2 s | `records_independent.out` |
| `freedom.py` | Published alternatives; declared rule space; lattice density (F-CT-4) | ~10 s | `freedom.out` |
| `charge_conventions.py` | Acting central charge under four readings; stabilizer dims (F-CT-3) | <1 s | `charge_conventions.out` |
| `weight_functionals.py` | Exploration: vacuum functionals of records vs weight law (F-CT-8) | ~3 s | `weight_functionals.out` |
| `readcap_lemmas.py` | SDP cross-check of imported read-cap lemmas: Choi map positive and non-decomposable, also after embedding (needs `cvxpy==1.6.0`) | ~4 s | `readcap_lemmas.out` |
| `config_enumeration.py` | Independent enumeration of Thm S-config (76/42/23) and of Theorem 1's selection under readings A, A*, B | <1 s | `config_enumeration.out` |
| `symmetry.py` | Covariance of the records under the symmetry group of the joint state; isotypic decomposition; commutant (F-CT-10) | ~10 s | `symmetry.out` |
| `operator_origin.py` | Bounded search for an operator origin of the weights $\sqrt{k_r}$ (CT-G01, CT-G12): covariance, Q\|q factorization and rationality of 38 primitives; $M$ as a polynomial in natural projectors; family F1 (all 56 355 words of length ≤ 3, tested as operator, diagonal and per-isotypic traces, also twirled, $|C|$, $\sqrt C$) and family F2 (321 489 survival × attribution dressings). Imports `records_independent.py` and `symmetry.py` (`formalization/chain-theory/alpha/operator-origin.md`) | ~16 min | `operator_origin.out` |
| `xcheck_operator_form.py` | Coordinator cross-check of `operator-origin.md` Prop. 2: applies $S$, $\Pi_X$, $\Pi_Y$ (partial traces, $E'$) to each covariant record; eigenvalues $\sqrt{k_r}$, trace, and the 12 records with no qutrit letter. Imports `symmetry.py` | ~1 min | `xcheck_operator_form.out` |
| `transport.py` | Transport of the qubit letters to the qutrit copies through the read links, then credit at the join (CT-G01, CT-G02; `formalization/chain-theory/alpha/transport.md`): uniqueness and $H=SU(2)\times U(1)$ symmetry of the transport; six E′-covariant transports on the 172 records; fiber obstruction (any $D$); rank-2 support of transported letters (credits $\sqrt3$, $2/\sqrt3$, $\sqrt2$); family F3 (338 688 frame-local dressings) and family F3a (words of length ≤ 3 in 24 covariant primitives and 9 transport letters, tests of `operator_origin.py`). Loads the definitions of `operator_origin.py` without rerunning its searches; 4 worker processes | ~13 min | `transport.out` |
| `mutations_values.py` | 1/value and σ for the readouts that the author's Lean `mutations.py` printed (log `upstream_mutations_lean.out`, Lean 4.35.0) | <1 s | `mutations_values.out` |
| `mutations_classified.py` | The 19 Lean mutation readouts of `upstream_mutations_lean.out`, row by row: exact numerator, 1/value, equal to a Table S5 row?, in the declared space B? (with a witness rule whose numerator equals Lean's exactly). Also checks all 16 closed-form rows of Table S5 against their printed inverses, with the relative change and B-membership of each (F-CT-4 A, D; referee CT-2, CT-022) | <1 s | `mutations_classified.out` |
| `mechanical_mutations.py` | Rule family generated mechanically from `ChainAlpha.lean`, for rule-family sensitivity (F-CT-4 E, CT-G10; `formalization/chain-theory/alpha/mechanical-trial-factor.md`): readout cone (44 definitions), 15 declared operators at 425 sites, 818 single mutants and 331 490 pairs with disjoint sites (mutation depth ≤ 2); exact readout of each via `chainalpha_transcription.py`; validation against `upstream_mutations_lean.out` (20/20 lines); F-CT-4 metrics (distinct values, density, and the illustrative uniform-target reference count of `freedom.py`, not a calibrated p-value) and comparison with B, Table S5 and the author's 19 rewordings, with every mutant that has the Born-weight value. Reads the pinned upstream file (needs `UPSTREAM` or `--lean`); imports `mutations_classified.py`; 2 worker processes. `--emit-lean DIR` / `--compare DIR` are the halves of `proofs/lean/chain-theory/mutation_spotcheck.sh` | ~17 min | `mechanical_mutations.out` |
| `chainalpha_transcription.py` | Python transcription of Parts 0–1 of `ChainAlpha.lean` with one hook per mutation site; module used by `mechanical_mutations.py`, checked against Lean on every single mutant (`proofs/lean/chain-theory/mutation_spotcheck.out`) | — | — |
| `algebraic.py` | Minimal polynomial, Galois conjugates, $f(d_Q,d_q)$, deviations from measurements (F-CT-6, F-CT-7) | <1 s | `algebraic.out` |

Run all of them from this directory in Linux/WSL Bash:

```bash
for s in records_independent freedom charge_conventions algebraic weight_functionals config_enumeration readcap_lemmas symmetry mutations_classified operator_origin xcheck_operator_form transport; do python3 $s.py > $s.out; done
```

The mechanical mutation family needs the pinned upstream file (`UPSTREAM` is a clone of the author's repository at `a28db89`, as for `proofs/lean/chain-theory/run.sh`):

```bash
UPSTREAM=/path/to/chain-theory python3 mechanical_mutations.py > mechanical_mutations.out
```

These scripts were written from the Letter's prose. They do not import or copy the author's code. The exceptions are `mechanical_mutations.py`, which reads the pinned `ChainAlpha.lean` as text in order to generate its mutants, and `chainalpha_transcription.py`, a line-by-line Python transcription of that file's Parts 0–1; neither runs the author's code. The cross-check against the author's `records.csv` is described in `formalization/chain-theory/alpha/audit.md` (F-CT-1).

The Lean companion is in `proofs/lean/chain-theory/` (see its `run.sh`; the Lean check of the mechanical mutants is `mutation_spotcheck.sh` → `mutation_spotcheck.out`, ~11 min with 2 Lean processes).
