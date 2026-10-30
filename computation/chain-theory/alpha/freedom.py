"""Freedom (look-elsewhere) analysis for the Chain Theory alpha candidate.

Part A recomputes, from the record classes alone, every alternative in the
Letter's Table S-alternatives that has a closed form.

Part B enumerates a declared rule space built from the same record classes
and asks how many rules land near alpha(0) and how dense the outputs are
near 137.036. This bounds the trial factor implied by the *stated*
alternatives; it cannot bound choices not visible in the paper.

Part C counts generic elements (a + b r2 + c r3 + d r6)/31104 with small
non-negative integer coefficients near alpha(0) (density of the number field
lattice at this denominator).

Run: python3 freedom.py
"""
import itertools
from collections import Counter
import mpmath as mp

mp.mp.dps = 50
r = {1: mp.mpf(1), 2: mp.sqrt(2), 3: mp.sqrt(3), 6: mp.sqrt(6)}

# CODATA 2022 (Mohr et al., RMP 97, 025002 (2025)): 1/alpha = 137.035 999 177(21)
INV_CODATA = mp.mpf('137.035999177')
SIG = mp.mpf('0.000000021')

# record classes (sizes from records_independent.py)
O, Sm, SQ, Sq, C = 79, 48, 36, 6, 3
N0 = 24 * 36 * 36

def inv_of(num, N=N0):
    return N / num

# ---------------- Part A ----------------
A = {
 "computed":                         O + Sm*r[2] + (SQ+Sq)*r[3] + C*r[6],
 "ends swapped":                     O + (SQ+Sq)*r[2] + Sm*r[3] + C*r[6],
 "both exchanges at qubit":          O + (Sm+SQ+Sq)*r[2] + C*2,          # k=2*2
 "both exchanges at qutrit":         O + (Sm+SQ+Sq)*r[3] + C*3,          # k=3*3
 "X where the letters sit":          O + Sm*r[2] + SQ*r[3] + Sq*r[2] + C*2,  # comm: k=2*2
 "qubit-qubit moved by Y":           O + Sm*r[2] + SQ*r[3] + Sq*r[6] + C*r[6],
 "copies as separate qudits":        O + (SQ+Sq)*r[3] + (Sm + C)*r[6],
 "survival w/o own sign":            O + C + Sm*r[2] + (SQ+Sq)*r[3],
 "commutators survive X only":       O + Sm*r[2] + (SQ+Sq+C)*r[3],
 "sums and differences kept":        O + Sm*r[2] + (SQ+Sq)*r[3] + 2*C*r[6],
 "no commutator records":            O + Sm*r[2] + (SQ+Sq)*r[3],
 "Born weights k":                   O + 2*Sm + 3*(SQ+Sq) + 6*C,
 "no weights":                       mp.mpf(O + Sm + SQ + Sq + C),
}
print("Part A: alternatives of Table S-alternatives (recomputed)")
for k, v in A.items():
    print(f"  {k:32s} 1/value = {mp.nstr(inv_of(v), 12)}")
# paper's printed inverses, for comparison
paper = {"computed": 137.036, "ends swapped": 135.894, "both exchanges at qubit": 146.524,
         "both exchanges at qutrit": 127.536, "X where the letters sit": 139.030,
         "qubit-qubit moved by Y": 134.485, "copies as separate qudits": 112.423,
         "survival w/o own sign": 139.713, "commutators survive X only": 138.348,
         "sums and differences kept": 132.739, "no commutator records": 141.621,
         "Born weights k": 97.505, "no weights": 180.837}
for k, v in A.items():
    assert abs(float(inv_of(v)) - paper[k]) < 6e-4, (k, inv_of(v), paper[k])
print("  all 13 match the paper's printed inverses to 3 decimals")
# extra normalization alternatives not listed in the paper
print("  normalization alternatives (not in the paper's table), computed numerator:")
for name, N in [("36^3 (all whole)", 36**3), ("24^3 (all traceless)", 24**3),
                ("24*35^2 (traceless part at reads)", 24*35*35), ("35^3", 35**3)]:
    print(f"    N = {name:36s} 1/value = {mp.nstr(inv_of(A['computed'], N), 10)}")

# ---------------- Part B ----------------
classes = [O, Sm, SQ, Sq]
Ns = {"24*36^2": N0, "36^3": 36**3, "24^3": 24**3, "24*35^2": 24*35*35, "35^3": 35**3}
values = {}
for ks in itertools.product([1, 2, 3, 6], repeat=4):          # factor per class
    for kc in [1, 2, 3, 6]:                                     # commutator factor
        for nc in [0, 3, 6]:                                    # commutator count
            for power in ["amp", "born"]:
                w = (lambda k: r[k]) if power == "amp" else (lambda k: mp.mpf(k))
                num = sum(n * w(k) for n, k in zip(classes, ks)) + nc * w(kc)
                for nN, N in Ns.items():
                    v = N / num
                    values[(ks, kc, nc, power, nN)] = v
distinct = {}
for key_, v in values.items():
    distinct.setdefault(mp.nstr(v, 30), key_)
print(f"\nPart B: rule space of {len(values)} rules, {len(distinct)} distinct values")
dev = sorted(((abs(v - INV_CODATA) / SIG, key_) for key_, v in values.items()), key=lambda t: t[0])
print("  five closest rules to CODATA 2022 (deviation in sigma):")
for d, key_ in dev[:5]:
    print(f"    {mp.nstr(d, 6):>12s} sigma   rule={key_}   1/value={mp.nstr(values[key_], 15)}")
win = [mp.mpf(x) for x in (1e-2, 1e-3, 1e-4, 1e-6)]
for w in win:
    n = sum(1 for v in {mp.nstr(v, 30): v for v in values.values()}.values()
            if abs(v - INV_CODATA) / INV_CODATA < w)
    print(f"  distinct values within relative {mp.nstr(w, 2)} of 1/alpha: {n}")
band = [v for v in {mp.nstr(v, 30): v for v in values.values()}.values() if 120 < v < 155]
dens = len(band) / 35
print(f"  density of distinct values in 1/alpha in (120,155): {mp.nstr(dens, 4)} per unit")
p = dens * 2 * SIG
print(f"  expected number of chance values within +-1 sigma: {mp.nstr(p, 4)}")
# relative spacing near target
band_s = sorted(band)
gaps = [b - a for a, b in zip(band_s, band_s[1:])]
print(f"  median gap between adjacent distinct values in (120,155): {mp.nstr(sorted(gaps)[len(gaps)//2], 4)}")
print(f"  1-sigma window (absolute, 1/alpha): {mp.nstr(2*SIG, 3)}")
print(f"  rules needed for O(1) chance hit at this density: ~{mp.nstr(1/(2*SIG*dens/len(values)) , 3)}")

# ---------------- Part C ----------------
# count (a,b,c,d) in [0,M]^4 with |31104/(a+b r2+c r3+d r6) - 1/alpha| < sigma
target_num = 31104 / INV_CODATA
tol = target_num * SIG / INV_CODATA         # window half-width in numerator space
print(f"\nPart C: numerator window half-width = {mp.nstr(tol, 4)} around {mp.nstr(target_num, 15)}")
import math
R2, R3, R6 = math.sqrt(2), math.sqrt(3), math.sqrt(6)
t, tl = float(target_num), float(tol)
for M in (50, 100, 200):
    hits = []
    for b in range(M + 1):
        for c in range(M + 1):
            for d in range(M + 1):
                rest = t - b * R2 - c * R3 - d * R6
                a = round(rest)
                if 0 <= a <= M and abs(rest - a) < tl:
                    hits.append((a, b, c, d))
    tot = (M + 1) ** 4
    # expected under uniform fractional parts: (#(b,c,d) with 0<=rest<=M) * 2*tol
    print(f"  coefficients in [0,{M}]: {tot:.3e} tuples, hits within 1 sigma: {len(hits)} {hits[:6]}")

# expected chance hits in Part C: count (b,c,d) admitting a in [0,M], times window width
for M in (100, 200):
    cnt = 0
    for b in range(M + 1):
        for c in range(M + 1):
            for d in range(M + 1):
                rest = t - b * R2 - c * R3 - d * R6
                if -0.5 <= rest <= M + 0.5:
                    cnt += 1
    print(f"  M={M}: (b,c,d) triples with admissible a: {cnt}; expected chance hits = {cnt*2*tl:.4f}")
