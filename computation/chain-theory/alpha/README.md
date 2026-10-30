# Chain Theory α audit — computations

Environment: Python 3.13.17, numpy 2.5.3, sympy 1.14.0, mpmath 1.3.0 (`pip install numpy==2.5.3 sympy==1.14.0 mpmath==1.3.0`). No randomness is used.

| Script | Purpose | Runtime | Output |
|---|---|---|---|
| `records_independent.py` | Operator-level rebuild of the 172 records; sectors; exact value; 60 digits (F-CT-1) | ~2 s | `records_independent.out` |
| `freedom.py` | Published alternatives; declared rule space; lattice density (F-CT-4) | ~10 s | `freedom.out` |
| `charge_conventions.py` | Acting central charge under four readings; stabilizer dims (F-CT-3) | <1 s | `charge_conventions.out` |
| `weight_functionals.py` | Exploration: vacuum functionals of records vs weight law (F-CT-8) | ~3 s | `weight_functionals.out` |
| `readcap_lemmas.py` | SDP cross-check of imported read-cap lemmas: Choi map positive and non-decomposable, also after embedding (needs `cvxpy==1.6.0`) | ~4 s | `readcap_lemmas.out` |
| `config_enumeration.py` | Independent enumeration of Thm S-config (76/42/23) and of Theorem 1's selection under readings A, A*, B | <1 s | `config_enumeration.out` |
| `algebraic.py` | Minimal polynomial, Galois conjugates, $f(d_Q,d_q)$, deviations from measurements (F-CT-6, F-CT-7) | <1 s | `algebraic.out` |

Run all of them:

```bash
for s in records_independent freedom charge_conventions algebraic weight_functionals config_enumeration readcap_lemmas; do python3 $s.py > $s.out; done
```

These scripts were written from the Letter's prose. They do not import or copy the author's code. The cross-check against the author's `records.csv` is described in `formalization/chain-theory/alpha/audit.md` (F-CT-1).

The Lean companion is in `proofs/lean/chain-theory/` (see its `run.sh`).
