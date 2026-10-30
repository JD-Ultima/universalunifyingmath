#!/usr/bin/env python3
"""Audit of compendium statement 141.1, "The electron dipole factor at order three from the
self-terms" (abstract: "the electron anomalous moment through third order").

Usage: python3 -I ae_audit.py
Self-contained: every compendium number below is transcribed from
  compendium/tex/parts/07-cell-3-2-photon.tex, lines 423-446 (chapter 7-2 / V-2), and
  compendium/lean/ChainTheory/ChainTheory/Data/D526870_PAET.lean (the chapter's data table).
QED inputs are cited in README.md.
"""
import itertools, math
from fractions import Fraction as F
import mpmath as mp
import sympy as sp

mp.mp.dps = 40
print(f"mpmath {mp.__version__}, sympy {sp.__version__}")

# ---------------------------------------------------------------- 1. QED coefficients
pi, z3, z5, l2 = mp.pi, mp.zeta(3), mp.zeta(5), mp.log(2)
C1 = mp.mpf(1)/2                                              # Schwinger 1948
C2 = mp.mpf(197)/144 + pi**2/12 - pi**2/2*l2 + mp.mpf(3)/4*z3  # Petermann 1957, Sommerfield 1957
a4 = mp.polylog(4, mp.mpf(1)/2)
C3 = (mp.mpf(83)/72*pi**2*z3 - mp.mpf(215)/24*z5
      + mp.mpf(100)/3*((a4 + l2**4/24) - pi**2*l2**2/24)
      - mp.mpf(239)/2160*pi**4 + mp.mpf(139)/18*z3 - mp.mpf(298)/9*pi**2*l2
      + mp.mpf(17101)/810*pi**2 + mp.mpf(28259)/5184)        # Laporta & Remiddi 1996
print("\n[1] QED mass-independent coefficients of (alpha/pi)^n")
print(f"  C1 = {mp.nstr(C1, 12)}")
print(f"  C2 = {mp.nstr(C2, 15)}   (target -0.328478965...)")
print(f"  C3 = {mp.nstr(C3, 15)}   (target  1.181241456...)")
assert abs(C2 - mp.mpf("-0.328478965579193784582")) < mp.mpf("1e-18")
assert abs(C3 - mp.mpf("1.181241456587")) < mp.mpf("1e-11")

# ---------------------------------------------------------------- 2. what the chapter prints
print("\n[2] Numbers printed by chapter 7-2 (07-cell-3-2-photon.tex:440-441)")
measures = {"(3,3)": 9, "(3,4)": 12, "(2,4)": 8, "(2,2)": 4, "(2,3)": 6}
print("  link measures 9,12,8,4,6 equal d1*d2 for the five link types:",
      all(int(k[1]) * int(k[3]) == v for k, v in measures.items()))
dressing, body8, chain5, made1 = 3, 67, 5 * 6, 6
T0, T1, T3 = 1, dressing * (body8 + chain5 + made1), 0
print(f"  T0 = {T0};  T1 = 3*(67+30+6) = {T1} (printed 309: {T1 == 309});  T1^2 = {T1**2} (printed 95481: {T1**2 == 95481});  T3 = {T3}")
print(f"  unstruck links 8+5+1 = {8+5+1} (printed fourteen); 16 links - 2 struck = {16-2}")
print(f"  '105 monomials of degree two': degree-2 monomials in 14 variables = {math.comb(15, 2)}; "
      f"unordered pairs of 15 sites = {math.comb(15, 2)} (both readings give 105)")
sols = [c for c in itertools.combinations_with_replacement(sorted(measures.values()), 8) if sum(c) == 67]
print(f"  multisets of 8 link measures from {{4,6,8,9,12}} summing to 67: {len(sols)} "
      f"(the chapter does not list the body's link types; e.g. {sols[0]})")
G2 = sp.Rational(5, 216) * (1 + sp.sqrt(2))
print(f"  G2(2,2) = 5/216 (1+sqrt2) = {sp.N(G2, 12)}")
print("  The chapter prints no coefficient of (alpha/pi)^n and states 'no digit of any constant displayed here'.")

# ---------------------------------------------------------------- 3. literal identification test
print("\n[3] Hypothesis L (literal): order-n record-sum term T_n is proportional to C_n (alpha/pi)^n")
print("    with one free normalisation lambda: T_n * lambda^n = C_n.")
lam = C1 / T1
print(f"  fit lambda from n=1: lambda = 1/(2*309) = {mp.nstr(lam, 10)}")
need_T2 = C2 / lam**2
print(f"  n=2 then requires T2 = 309^2 + closure = {mp.nstr(need_T2, 12)}, i.e. closure = {mp.nstr(need_T2 - 95481, 12)};"
      f" closure is not printed -> untestable")
print(f"  n=3 requires T3 = C3/lambda^3 = {mp.nstr(C3/lam**3, 12)}; printed T3 = 0 -> FAILS for every lambda (C3 != 0)")

# ---------------------------------------------------------------- 4. experimental weight of C3
print("\n[4] Size of the third-order QED term against experiment")
a_exp, s_exp = mp.mpf("1.15965218059e-3"), mp.mpf("0.13e-12")          # Fan et al. PRL 130, 071801 (2023)
for label, ainv, sinv in [("Rb 2020", "137.035999206", "0.000000011"),   # Morel et al. Nature 588, 61 (2020)
                          ("Cs 2018", "137.035999046", "0.000000027"),   # Parker et al. Science 360, 191 (2018)
                          ("CODATA 2022", "137.035999177", "0.000000021")]:
    al = 1/mp.mpf(ainv); x = al/pi
    term3 = C3 * x**3
    s_al = a_exp * mp.mpf(sinv)/mp.mpf(ainv)
    sig = mp.sqrt(s_exp**2 + s_al**2)
    print(f"  alpha({label}): C3 (alpha/pi)^3 = {mp.nstr(term3, 6)}; sigma(exp (+) alpha) = {mp.nstr(sig, 3)}; "
          f"setting C3 = 0 shifts a_e by {mp.nstr(term3/sig, 4)} sigma")
    if label == "Rb 2020":
        print(f"     (alpha/pi) = {mp.nstr(x, 10)}; C1 x = {mp.nstr(C1*x, 10)}; C2 x^2 = {mp.nstr(C2*x**2, 6)}; "
              f"relative weight of C3 term in a_e = {mp.nstr(term3/a_exp, 4)}")

# ---------------------------------------------------------------- 5. the exact numbers in the Lean data table
print("\n[5] Exact numbers in D526870_PAET (a|b|c|d read as a + b*sqrt2 + c*sqrt3 + d*sqrt6)")
def q(a, b, c, d):
    return sp.Rational(a) + sp.Rational(b)*sp.sqrt(2) + sp.Rational(c)*sp.sqrt(3) + sp.Rational(d)*sp.sqrt(6)
W1 = q("11989/6879707136", 0, "1675/1719926784", 0)
others = {"W2": q("299725/18246703251456", 0, "41875/4561675812864", 0),
          "W3": q("53818621/72986813005824", 0, "7519075/18246703251456", 0),
          "PWONCE": q("11989/901071765504", 0, "1675/225267941376", 0)}
print(f"  W1 (the 'EQUAL' reading) self-term at link L1 = {sp.N(W1, 10)}; denominators "
      f"{sp.factorint(6879707136)}, {sp.factorint(1719926784)}")
for k, v in others.items():
    print(f"  {k}: {sp.N(v, 10)}  ratio to W1 = {sp.N(v/W1, 8)}")
print(f"  for scale only (no identification is made by the source): alpha/2pi = {mp.nstr(1/(2*pi*mp.mpf('137.035999177')), 10)}")
print("  The Lean theorem row_526870_PAET proves only: the table has 8 rows; one row has test=ORDER-TWO-READOUT,"
      " reading=W1, of_record=YES, letter=EQUAL; one row has letter=ZERO; three rows have letter=PARTING; ...")

# ---------------------------------------------------------------- 6. tally
print("\n[6] Tally: QED coefficients reproduced = 0 of 3; C3 contradicted under hypothesis L;"
      " under the charitable reading (T_n not tied to alpha/pi) the statement makes no contact with a_e.")
