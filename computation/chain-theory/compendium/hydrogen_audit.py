#!/usr/bin/env python3
"""Audit of compendium statement 143.1, "Hydrogen from the held carrier" (abstract: "the
hydrogen spectrum"; section title "Hydrogen: the Rydberg ladder and the fine structure").

Usage: python3 -I hydrogen_audit.py
Compendium numbers are transcribed from compendium/tex/parts/07-cell-3-2-photon.tex,
lines 478-548 (chapter 7-4 / V-4) and the hyperfine open result 144.1 (lines 550-568).
Physical inputs are cited in README.md.
"""
from fractions import Fraction as F
import sympy as sp
import mpmath as mp
mp.mp.dps = 30

print("[1] Exact arithmetic of the printed level-one binding")
E1w = F(36) * F(201, 2**30 * 3**18)
E1c = F(36) * F(67, 2**31 * 3**17)
print(f"  E1_whole = 36*201/(2^30 3^18) = {E1w}  == 67/(2^28 3^15): {E1w == F(67, 2**28 * 3**15)}")
print(f"  E1_cut   = 36*67/(2^31 3^17)  = {E1c}  == 67/(2^29 3^15): {E1c == F(67, 2**29 * 3**15)}")
print(f"  E1_whole / E1_cut = {E1w / E1c}   (two 'faces' differ by an exact factor 2; both printed, none selected)")
print(f"  side readings: 36*67/36 = {F(36*67, 36)}, 36*67/24 = {F(36*67, 24)}")
print(f"  float: E1_whole = {float(E1w):.6e}, E1_cut = {float(E1c):.6e}  (in 'the electron unit')")

print("\n[2] Comparison with hydrogen, reading 'electron unit' = m_e c^2 (the reading 176.1 uses for m_p/m_e)")
ainv = mp.mpf("137.035999177")          # CODATA 2022
al = 1/ainv
mp_me = mp.mpf("1836.152673426")       # CODATA 2022 m_p/m_e
mu = 1/(1 + 1/mp_me)
bohr = al**2/2*mu
IE_H = mp.mpf("13.598434599702")       # eV, H ionization energy (NIST ASD / Jentschura et al.)
mec2 = mp.mpf("510998.95069")          # eV, CODATA 2022
meas = IE_H/mec2
for name, v in [("E1_whole", E1w), ("E1_cut", E1c)]:
    x = mp.mpf(v.numerator)/v.denominator
    print(f"  {name}: model/Bohr(reduced mass) = {mp.nstr(x/bohr, 6)}; model/measured IE = {mp.nstr(x/meas, 6)}; "
          f"log10 gap = {mp.nstr(mp.log10(meas/x), 5)}")
    print(f"     effective exponent k with 2*E1 = alpha^k: k = {mp.nstr(mp.log(2*x)/mp.log(al), 6)} (Bohr: k = 2)")
print(f"  Bohr reduced-mass value alpha^2 mu/(2 m_e) = {mp.nstr(bohr, 10)}; measured IE/m_e c^2 = {mp.nstr(meas, 10)}")
print("  E1 is a pure 2-3-67 rational: it contains no factor of the model's alpha (124.1), so the physical")
print("  relation E1 = alpha^2 m_e/2 cannot hold for any value of the model's alpha other than")
for name, v in [("E1_whole", E1w), ("E1_cut", E1c)]:
    x = mp.mpf(v.numerator)/v.denominator
    print(f"     {name}: alpha_required = sqrt(2 E1/mu) = 1/{mp.nstr(1/mp.sqrt(2*x/mu), 10)}")

print("\n[3] The 1/n^2 ladder: reconstruction and the choice that fixes the exponent")
print("  Source: 'the held pendant pair at rank n is the washed symmetric power of n less one passages")
print("  along the bond ... the register of the ladder equal to the dimension of the register of the pair';")
print("  'the two-vertex ladder equals one over the rank squared at all eight ranks'.")
print("  Reconstruction R: the washed (twirled) state on Sym^(n-1)(C^m) is I/dim; a two-vertex term reads it")
print("  twice, giving (1/dim)^2. With m = 2 (the common support of the (2,3) pair) dim = n.")
n = sp.symbols("n", positive=True, integer=True)
for m in (2, 3, 4):
    dim = sp.binomial(n + m - 2, m - 1)
    ladder = sp.simplify(1/dim**2)
    vals = [sp.Rational(1, sp.binomial(k + m - 2, m - 1)**2) for k in range(1, 9)]
    bohr_ok = all(vals[k-1] == sp.Rational(1, k*k) for k in range(1, 9))
    print(f"  m={m}: dim Sym^(n-1)(C^{m}) = {sp.factor(sp.expand_func(dim))}; ladder (1/dim)^2 = {sp.factor(sp.expand_func(ladder))}; "
          f"equals 1/n^2 for n=1..8: {bohr_ok}")
# check of R: Haar twirl of a pure state on Sym^(n-1)(C^2) is I/n (Schur's lemma; spin-(n-1)/2 irrep).
import itertools
import numpy as np
def sym_basis(k):
    """orthonormal basis of the symmetric subspace of (C^2)^{(x)k}, as columns"""
    cols = []
    for j in range(k + 1):
        v = np.zeros(2**k)
        for bits in itertools.product((0, 1), repeat=k):
            if sum(bits) == j:
                v[int("".join(map(str, bits)) or "0", 2)] = 1.0
        cols.append(v/np.linalg.norm(v))
    return np.array(cols).T
print("  check of R on Sym^(n-1)(C^2): spectrum of the Haar SU(2)-twirl of a pure state (Monte Carlo, 20000 draws):")
rng = np.random.default_rng(20261030)
for k in range(0, 4):
    B = sym_basis(k) if k else np.ones((1, 1))
    acc = np.zeros((k+1, k+1), dtype=complex)
    N = 20000
    for _ in range(N):
        q = rng.normal(size=4); q /= np.linalg.norm(q)
        U = np.array([[q[0]+1j*q[1], -q[2]+1j*q[3]], [q[2]+1j*q[3], q[0]-1j*q[1]]])
        Uk = np.ones((1, 1), dtype=complex)
        for _i in range(k): Uk = np.kron(Uk, U)
        M = B.conj().T @ Uk @ B
        v = M[:, 0]
        acc += np.outer(v, v.conj())
    ev = np.linalg.eigvalsh(acc / N)
    print(f"    n={k+1}: eigenvalues {np.round(ev, 3)}  (1/n = {1/(k+1):.3f}); two-vertex read (1/n)^2 = {1/(k+1)**2:.4f}")

print("\n[4] What a hydrogen spectrum requires vs what 143.1 prints")
rows = [
 ("E_n proportional to 1/n^2", "printed for n=1..8 (two-vertex ladder)", "matches the n-dependence (one integer exponent)"),
 ("E_1 = alpha^2 m_e/2 (13.6 eV)", "67/(2^28 3^15) or 67/(2^29 3^15) m_e", "off by ~9.2 orders; no alpha dependence"),
 ("fine structure ~ alpha^2 E_n/n, j-dependent", "'two-vertex term times the square of the cut gain'; no value; posed step", "untestable as printed"),
 ("Lamb shift 2S1/2-2P1/2 = +1057.8 MHz", "'one exact element, positive and split'; value not printed", "sign only (1 bit) matches"),
 ("hyperfine 1420.405751768 MHz", "open result 144.1: interval = 0 at the standing transport", "graded open; the computed reading is 0"),
 ("infinitely many bound levels", "ranks 1..8 enumerated", "range of the check, not a stated cap"),
]
for r in rows: print("  - " + " | ".join(r))
lamb_MHz, hf_MHz = 1057.845, 1420.405751768
print(f"  (reference: Lamb shift {lamb_MHz}(9) MHz, Lundeen & Pipkin, PRL 46, 232 (1981);"
      f" hyperfine {hf_MHz} MHz, Hellwig et al., IEEE Trans. Instrum. Meas. 19, 200 (1970))")
