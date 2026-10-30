"""Exploration for gap CT-G01: do natural functionals of the records on the
Letter's joint state reproduce the weight law sqrt(k_r) in {1, sqrt2, sqrt3, sqrt6}?

Joint state (Supplement S7, matrices S10): after the join P3 on (Q1,Q2) the
read part is |omega_3>_{Q1Q2} (x) |omega_2>_{q1q2} (we verify this).
For each record operator O_r (normalized to Hilbert-Schmidt norm 1) we compute
  f_vev  = |<Omega|O_r|Omega>|
  f_norm = ||O_r |Omega>||^2 * 36          (36 = dim H; equals 1 on average)
  f_twirl= HS norm of the U(x)Ubar-twirl of O_r over the matched pairs (Q1Q2),(q1q2)
and report, per record class, the set of values taken.  Exploratory: a match
with sqrt(k) would give the weight law an operator origin; no match is
informative too.
Run: python3 weight_functionals.py
"""
import itertools, numpy as np
from collections import defaultdict
import importlib.util, sys, io, contextlib
spec = importlib.util.spec_from_file_location("rec", __file__.replace("weight_functionals.py", "records_independent.py"))
rec = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(rec)

DIMS = rec.DIMS
def ket(idx):
    v = np.zeros(36); v[np.ravel_multi_index(idx, DIMS)] = 1; return v
# |omega2>_{Q1 q1} |omega2>_{Q2 q2}, qutrit copies use levels 0,1 only
v = sum(ket((i, i, j, j)) for i in range(2) for j in range(2)) / 2
# P3 on (Q1,Q2): project onto |omega3>
w3 = np.zeros((3, 3)); np.fill_diagonal(w3, 1 / np.sqrt(3))
P = np.zeros((36, 36))
for idx in itertools.product(*[range(d) for d in DIMS]):
    for jdx in itertools.product(*[range(d) for d in DIMS]):
        if idx[1] == jdx[1] and idx[3] == jdx[3]:
            P[np.ravel_multi_index(idx, DIMS), np.ravel_multi_index(jdx, DIMS)] = w3[idx[0], idx[2]] * w3[jdx[0], jdx[2]]
Om = P @ v
print("norm^2 of projected vector:", round(Om @ Om, 12), "(Supplement: 1/6)")
Om = Om / np.linalg.norm(Om)
target = sum(ket((i, a, i, a)) for i in range(3) for a in range(2)) / np.sqrt(6)
print("Omega == |w3>_{Q1Q2}|w2>_{q1q2}:", np.allclose(abs(Om @ target), 1))

def cls(r):
    kind, A, e, cat = r
    if cat == 'comm': return 'comm (k=6)'
    if e == -1: return 'odd (k=1)'
    return 'even same-end (k=3)' if cat == 'same' else 'even mixed (k=2)'

# twirl over U (x) Ubar on (Q1,Q2) and V (x) Vbar on (q1,q2): commutant projector
# = projection onto span{1, P_omega} on each pair  (Horodecki 1999)
def pair_proj(dA, dB, which):  # operators on Q1Q2 or q1q2 subsystem as 4-index
    pass

out = defaultdict(set)
for r in rec.records:
    A = r[1] / np.linalg.norm(r[1])
    vev = abs(Om.conj() @ A @ Om)
    nrm = np.linalg.norm(A @ Om) ** 2 * 36
    out[cls(r)].add((round(vev, 6), round(nrm, 6)))
for c in ['odd (k=1)', 'even mixed (k=2)', 'even same-end (k=3)', 'comm (k=6)']:
    vals = sorted(out[c])
    print(f"{c:22s} distinct (|vev|, 36*||O Omega||^2): {vals[:8]}{' ...' if len(vals) > 8 else ''} ({len(vals)} distinct)")

# class means (basis-independent summaries)
S = defaultdict(list)
for r in rec.records:
    A = r[1] / np.linalg.norm(r[1])
    S[cls(r)].append((np.linalg.norm(A @ Om) ** 2 * 36, abs(Om.conj() @ A @ Om)))
print("class means:")
for c in ['odd (k=1)', 'even mixed (k=2)', 'even same-end (k=3)', 'comm (k=6)']:
    v = np.array(S[c])
    print(f"  {c:22s} n={len(v):3d}  mean 36||O Omega||^2 = {v[:,0].mean():.6f}   sum |vev| = {v[:,1].sum():.6f}")
