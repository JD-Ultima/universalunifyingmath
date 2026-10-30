#!/usr/bin/env python3
"""Quick numeric checks for Claim-Register entries outside the three deep audits.

Usage: python3 -I register_checks.py
Compendium numbers transcribed from the cited lines; measured values cited inline.
"""
from fractions import Fraction as F
import mpmath as mp
mp.mp.dps = 25

print("[W] 164.5 'The weak angle as the ratio of offices at the Z body' (08-cell-4-2-weak.tex:436-452)")
r = F(784, 101761) / F(1568, 101761)
print(f"  printed offices 784/101761 and 1568/101761 -> ratio {r}; source: 'the literature calls [this ratio]")
print("  the squared cosine of the weak mixing angle' (correspondence by name)")
MW, sMW = mp.mpf("80.3692"), mp.mpf("0.0133")   # PDG 2024 world average W mass, GeV
MZ, sMZ = mp.mpf("91.1880"), mp.mpf("0.0020")   # PDG 2024 Z mass, GeV
c2 = (MW/MZ)**2
sc2 = c2*mp.sqrt((2*sMW/MW)**2 + (2*sMZ/MZ)**2)
print(f"  on-shell cos^2(theta_W) = (M_W/M_Z)^2 = {mp.nstr(c2, 7)} +- {mp.nstr(sc2, 2)}; printed 1/2 differs by "
      f"{mp.nstr((c2 - mp.mpf(1)/2)/sc2, 5)} sigma (M_W/M_Z predicted 1/sqrt2 = {mp.nstr(1/mp.sqrt(2), 6)} vs {mp.nstr(MW/MZ, 6)})")
print(f"  MS-bar sin^2(theta_W)(M_Z) = 0.23122 (PDG 2024) -> cos^2 = 0.76878; printed value implies sin^2 = 1/2")

print("\n[R] 173.1 'The strong coupling at the depth scale' (09-cell-3-3-strong.tex:101)")
print("  printed: coupling ratio between successive rungs = 1/36; scale at rung = 96 x depth (count metric).")
print("  => alpha(rung k) proportional to 36^(-k) with scale linear in k: alpha ~ exp(-(ln 36/96) * scale).")
print(f"     ln36/96 = {mp.nstr(mp.log(36)/96, 6)} per count-metric unit. One-loop QCD instead gives")
print("     1/alpha_s linear in ln(Q^2): alpha_s ~ 1/(b0 ln(Q^2/Lambda^2)), b0 = (33 - 2 n_f)/(12 pi).")
for nf in (3, 4, 5, 6):
    print(f"     n_f={nf}: b0 = {mp.nstr((33-2*nf)/(12*mp.pi), 6)}")
print("  The map from 'count metric' to energy is not given, so the functional form is the only testable content:")
print("  geometric (exponential) running in a scale variable vs logarithmic running in Q. Untestable as printed.")

print("\n[K] 157.1 'The Koide relation at the dock unit' (07-cell-3-2-photon.tex:1161) -- OPEN; printed only that")
print("  the Koide ratio equals 2/3 'in the field at no reading'. Measured (PDG 2024 masses, MeV):")
me, mmu, mtau = mp.mpf("0.51099895"), mp.mpf("105.6583755"), mp.mpf("1776.93")
Q = (me + mmu + mtau)/(mp.sqrt(me) + mp.sqrt(mmu) + mp.sqrt(mtau))**2
print(f"  Q = {mp.nstr(Q, 10)} (2/3 = {mp.nstr(mp.mpf(2)/3, 10)}); the source predicts no number to compare.")
