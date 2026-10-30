# Corrections and retractions log

Entries are appended, never edited. Format: date — file — what changed — why.

- 2026-10-30 — `computation/chain-theory/alpha/freedom.py` (before first commit) — in three of the paper's alternatives ("both exchanges at the qubit/qutrit end", "X where the letters sit"), the commutator records' weight was first coded as $\sqrt d$ instead of $\sqrt{d\cdot d}=d$. The script's own check against the paper's printed values caught it. Fixed before any result was reported.
