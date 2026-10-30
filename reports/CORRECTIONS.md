# Corrections and retractions log

Entries are appended, never edited. Format: date — file — what changed — why.

- 2026-10-30 — `computation/chain-theory/alpha/freedom.py` (before first commit) — in three of the paper's alternatives ("both exchanges at the qubit/qutrit end", "X where the letters sit"), the commutator records' weight was first coded as $\sqrt d$ instead of $\sqrt{d\cdot d}=d$. The script's own check against the paper's printed values caught it. Fixed before any result was reported.
- 2026-10-30 — `computation/chain-theory/alpha/freedom.py`, F-CT-4 in `audit.md`, `reports/chain-theory-alpha.md` (commit 1 → commit 2) — Commit 1 said the 30 720-rule space "contains every alternative the author lists". It did not: three published alternatives give commutator records $k=4$ or $k=9$. The commutator factor set is now $\{1,2,3,4,6,9\}$, and a programmatic containment check has been added. Revised numbers: 46 080 rules, 15 639 distinct values, density 103.9/unit, expected chance hits within 1σ $4.4\times10^{-6}$ (was $3.2\times10^{-6}$). The qualitative conclusion is unchanged.
