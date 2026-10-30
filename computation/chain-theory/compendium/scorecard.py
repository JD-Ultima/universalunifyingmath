#!/usr/bin/env python3
"""Numeric scorecard of the Chain Theory compendium 0.1 (round 2, tasks 2-3).

Usage: python3 -I scorecard.py --src <compendium_dir> [--pdg-mcd PATH/mass_width_2024.mcd] [--tsv OUT.tsv]

For every row of scorecard_data.ROWS and every face (printed reading) it evaluates the face exactly
(sympy), converts it to the measured quantity named by the row, and compares it with the reference
value: deviation in sigma, relative deviation, factor. It then prints the summary counts of task 3.
Block REDERIVE re-derives, in exact arithmetic, the numbers that follow from the construction as the
page states it.

References
  CODATA 2022: read at run time from scipy.constants (scipy 1.18.1 ships CODATA 2022); each value used
               is asserted against the transcription below.
  PDG 2024 masses and widths: transcribed from the PDG Monte Carlo table 'mass_width_2024.mcd'
               (Review of Particle Physics 2024, file generated 31 May 2024; distributed in the
               scikit-hep 'particle' package 0.25.2). With --pdg-mcd the transcription is checked.
               M_W uses the PDG 2024 summary value 80.3692(133) GeV (the MC file rounds it to 80.369(13)).
  PDG 2024 other (transcribed, not machine-checked): alpha_s(M_Z), alpha^-1(M_Z) MS-bar, |V_ud|, J,
               sin^2 theta_13, N_nu, cosmological parameters (Planck 2018 values as quoted by PDG 2024).
"""
import math, pathlib, sys
import sympy as sp
import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))   # our own directory only
import scorecard_data as D

mp.mp.dps = 40
S = sp.sympify

# ------------------------------------------------------------------ references
import scipy, scipy.constants as sc
assert scipy.__version__.startswith("1.18"), scipy.__version__
PC = sc.physical_constants
CODATA = {  # name: (value, sigma) as transcribed; asserted against scipy.constants
    "inverse fine-structure constant": (137.035999177, 2.1e-8),
    "fine-structure constant": (7.2973525643e-3, 1.1e-12),
    "muon-electron mass ratio": (206.7682827, 4.6e-6),
    "tau-electron mass ratio": (3477.23, 0.23),
    "proton-electron mass ratio": (1836.152673426, 3.2e-8),
    "neutron-proton mass ratio": (1.00137841946, 4.0e-10),
    "Fermi coupling constant": (1.1663787e-5, 6e-12),
    "electron mass energy equivalent in MeV": (0.51099895069, 1.6e-10),
    "proton mass energy equivalent in MeV": (938.27208943, 2.9e-7),
    "neutron mass energy equivalent in MeV": (939.56542194, 4.8e-7),
    "deuteron mass energy equivalent in MeV": (1875.61294500, 5.8e-7),
    "Newtonian constant of gravitation": (6.67430e-11, 1.5e-15),
    "weak mixing angle": (0.22305, 0.00023),
    "electron mass": (9.1093837139e-31, 2.8e-40),
    "reduced Planck constant": (1.0545718176461565e-34, 0.0),
    "speed of light in vacuum": (299792458.0, 0.0),
    "Planck length": (1.616255e-35, 1.8e-40),
}
for k, (v, e) in CODATA.items():
    pv, _, pe = PC[k]
    assert abs(pv - v) <= 1e-12 * abs(v) and abs(pe - e) <= 1e-6 * max(abs(e), 1e-300) + 1e-300, (k, pv, pe)
C = {k: v for k, (v, _) in CODATA.items()}; CS = {k: e for k, (_, e) in CODATA.items()}

PDG_MC = {  # pdgid: (mass GeV, err, width GeV, err) from mass_width_2024.mcd
    24: (8.0369e1, 1.3e-2, 2.08, 4.0e-2), 23: (9.11880e1, 2.0e-3, 2.4955, 2.3e-3), 25: (1.2520e2, 1.1e-1, 3.7e-3, 1.9e-3),
    6: (1.7257e2, 2.9e-1, 1.42, 1.9e-1), 11: (5.1099895000e-4, 1.5e-13, 0.0, 0.0), 13: (1.056583755e-1, 2.3e-9, 2.9959836e-19, 3.0e-25),
    15: (1.77693, 9.0e-5, 2.267e-12, 4.0e-15), 1: (4.70e-3, 0.7e-4, None, None), 2: (2.16e-3, 0.7e-4, None, None),
    3: (9.35e-2, 0.8e-3, None, None), 4: (1.273, 5.0e-3, None, None), 5: (4.183, 7.0e-3, None, None),
}
MW, sMW = 80.3692, 0.0133                     # PDG 2024 summary table (W mass average)


def check_mcd(path):
    rows = {}
    for line in open(path, encoding="ascii", errors="replace"):
        if line.startswith("*") or len(line) < 70: continue
        ids = [int(x) for x in line[:32].split()]
        if not line[33:51].strip(): continue
        m = float(line[33:51]); me = float(line[52:60])
        w = line[70:88].strip()
        wv = float(w) if w else None
        we = float(line[89:97]) if w else None
        for i in ids: rows[i] = (m, me, wv, we)
    for k, (m, me, w, we) in PDG_MC.items():
        r = rows[k]
        assert abs(r[0] - m) <= 1e-9 * abs(m) and abs(r[1] - me) <= 1e-9 * abs(me), (k, r)
        if w is not None and w != 0.0:
            assert abs(r[2] - w) <= 1e-9 * abs(w), (k, r)
    return len(rows)


me = C["electron mass energy equivalent in MeV"] * 1e-3      # GeV
GF = C["Fermi coupling constant"]
v = (math.sqrt(2) * GF) ** -0.5                               # GeV
MZ, sMZ = PDG_MC[23][0], PDG_MC[23][1]
MH, sMH = PDG_MC[25][0], PDG_MC[25][1]
mt, smt = PDG_MC[6][0], PDG_MC[6][1]
c2 = (MW / MZ) ** 2; sc2 = c2 * math.hypot(2 * sMW / MW, 2 * sMZ / MZ)
# Koide from PDG 2024 lepton masses
mle = [PDG_MC[11][0], PDG_MC[13][0], PDG_MC[15][0]]; sle = [PDG_MC[11][1], PDG_MC[13][1], PDG_MC[15][1]]
def koide(m): return sum(m) / sum(math.sqrt(x) for x in m) ** 2
kq = koide(mle); kq_s = math.sqrt(sum(((koide([x + (h if i == j else 0) for j, x in enumerate(mle)]) - kq) / h * sle[i]) ** 2
                                     for i, h in enumerate([1e-9, 1e-7, 1e-6])))
# cosmology (Planck 2018 as quoted by PDG 2024)
H0 = 67.4e3 / 3.0856775814913673e22; OL = 0.6847
Lam = 3 * OL * H0 ** 2 / C["speed of light in vacuum"] ** 2                    # m^-2
lP = C["Planck length"]; lC = C["reduced Planck constant"] / (C["electron mass"] * C["speed of light in vacuum"])
Och, Obh, sOch, sObh = 0.1200, 0.02237, 0.0012, 0.00015
fdm = Och / (Och + Obh); sfdm = math.hypot(Obh / (Och + Obh) ** 2 * sOch, Och / (Och + Obh) ** 2 * sObh)
mn_MeV, mp_MeV, md_MeV = (C["neutron mass energy equivalent in MeV"], C["proton mass energy equivalent in MeV"],
                          C["deuteron mass energy equivalent in MeV"])
Vud, sVud = 0.97367, 0.00032
mu_md = PDG_MC[2][0] / PDG_MC[1][0]; smu_md = mu_md * math.hypot(PDG_MC[2][1] / PDG_MC[2][0], PDG_MC[1][1] / PDG_MC[1][0])
lamH = MH ** 2 / (2 * v ** 2); slamH = lamH * 2 * sMH / MH

# value, sigma (None: no meaningful sigma), citation
REF = {
    "alpha_inv": (C["inverse fine-structure constant"], CS["inverse fine-structure constant"], "CODATA 2022"),
    "C1": (0.5, 0.0, "QED one-loop coefficient (Schwinger 1948)"),
    "C3": (1.181241456587, 1e-12, "QED three-loop coefficient (Laporta & Remiddi 1996)"),
    "E1H": (C["fine-structure constant"] ** 2 / 2 / (1 + 1 / C["proton-electron mass ratio"]), None,
            "alpha^2 mu/(2 m_e), CODATA 2022 alpha and m_p/m_e (reduced-mass Bohr level)"),
    "alphaMZ": (1 / 127.951, None, "PDG 2024: alpha^-1(M_Z) = 127.951(9), MS-bar (transcribed)"),
    "alphas": (0.1180, 0.0009, "PDG 2024 world average alpha_s(M_Z) (transcribed)"),
    "s2w_os": (1 - c2, sc2, "on-shell 1 - (M_W/M_Z)^2, PDG 2024 M_W = 80.3692(133), M_Z = 91.1880(20) GeV"),
    "GFme2": (GF * me ** 2, None, "G_F m_e^2, CODATA 2022"),
    "MW_me": (MW / me, MW / me * sMW / MW, "PDG 2024 M_W / CODATA 2022 m_e"),
    "MH_me": (MH / me, MH / me * sMH / MH, "PDG 2024 M_H = 125.20(11) GeV / CODATA 2022 m_e"),
    "MH_MW": (MH / MW, MH / MW * math.hypot(sMH / MH, sMW / MW), "PDG 2024"),
    "lambdaH": (lamH, slamH, "tree level M_H^2/(2 v^2), v = (sqrt2 G_F)^(-1/2) = 246.2197 GeV"),
    "v_me": (v / me, None, "v = (sqrt2 G_F)^(-1/2), CODATA 2022"),
    "yt": (math.sqrt(2) * mt / v, math.sqrt(2) * smt / v, "sqrt2 m_t/v, PDG 2024 m_t = 172.57(29) GeV"),
    "GZ_GW": (PDG_MC[23][2] / PDG_MC[24][2], PDG_MC[23][2] / PDG_MC[24][2] * PDG_MC[24][3] / PDG_MC[24][2], "PDG 2024 widths"),
    "GH_GW": (PDG_MC[25][2] / PDG_MC[24][2], None, "PDG 2024 widths (Gamma_H = 3.7 +1.9 -1.4 MeV)"),
    "Gtau_Gmu": (PDG_MC[15][2] / PDG_MC[13][2], None, "PDG 2024 widths"),
    "GH_Gmu": (PDG_MC[25][2] / PDG_MC[13][2], None, "PDG 2024 widths"),
    "J": (3.08e-5, 0.14e-5, "PDG 2024 CKM review, J = (3.08 +0.15 -0.13)e-5 (transcribed)"),
    "Vud2": (Vud ** 2, 2 * Vud * sVud, "PDG 2024 |V_ud| = 0.97367(32) (transcribed)"),
    "Ue3": (0.0220, 0.0007, "PDG 2024 sin^2 theta_13 = 0.0220(7) (transcribed)"),
    "mmu_me": (C["muon-electron mass ratio"], CS["muon-electron mass ratio"], "CODATA 2022"),
    "mtau_me": (C["tau-electron mass ratio"], CS["tau-electron mass ratio"], "CODATA 2022"),
    "koide": (kq, kq_s, "PDG 2024 lepton masses"),
    "mp_me": (C["proton-electron mass ratio"], CS["proton-electron mass ratio"], "CODATA 2022"),
    "mn_mp": (C["neutron-proton mass ratio"], CS["neutron-proton mass ratio"], "CODATA 2022"),
    "dnp": ((mn_MeV - mp_MeV) / mp_MeV, None, "CODATA 2022"),
    "mu_md": (mu_md, smu_md, "PDG 2024 m_u = 2.16(7), m_d = 4.70(7) MeV (MS-bar, 2 GeV)"),
    "mt_mc": (mt / PDG_MC[4][0], None, "PDG 2024 m_t (pole) / m_c(m_c) MS-bar"),
    "ms_mc": (PDG_MC[3][0] / PDG_MC[4][0], None, "PDG 2024 MS-bar masses"),
    "mt_MW": (mt / MW, mt / MW * math.hypot(smt / mt, sMW / MW), "PDG 2024"),
    "Bd_me": ((mp_MeV + mn_MeV - md_MeV) / (me * 1e3), None, "CODATA 2022 masses, B_d = m_p + m_n - m_d"),
    "alphaG": (C["Newtonian constant of gravitation"] * C["electron mass"] ** 2 / (C["reduced Planck constant"] * C["speed of light in vacuum"]),
               None, "G m_e^2/(hbar c), CODATA 2022"),
    "BH": (0.25, 0.0, "Bekenstein-Hawking S = A/(4 l_P^2) (theory)"),
    "fdm": (fdm, sfdm, "Omega_c h^2 = 0.1200(12), Omega_b h^2 = 0.02237(15) (Planck 2018 via PDG 2024; transcribed)"),
    "tsirelson": (2 * math.sqrt(2), 0.0, "Tsirelson bound 2 sqrt2; CGLMP max-entangled value 4(3+2sqrt3)/9 (Collins et al. 2002)"),
}
QUAL = {  # kind int/sign/zero/form/nc: the reference statement and the comparison rule
    "signature": ("(1,3) (Lorentzian)", lambda f: f == "(1,3)"),
    "gauge": ("su(3)+su(2)+u(1), dims 8+3+1 = 12; global quotient Z6 consistent, not established", lambda f: f == "(12,6)"),
    "Nnu": ("N_nu = 2.9963(74) (PDG 2024, Z line shape; transcribed)", lambda f: f == "3"),
    "photon": ("m_gamma = 0 (< 1e-18 eV), two helicities", lambda f: f == "(0,2)"),
    "C3": ("C3 = 1.1812414566 != 0", lambda f: f != "0"),
    "bohr_exp": ("E_n proportional to 1/n^2 (Bohr; fine structure aside)", lambda f: f == "2"),
    "lamb": ("2S1/2 - 2P1/2 = +1057.845(9) MHz (Lundeen & Pipkin 1981)", lambda f: f == "1"),
    "hfs": ("1420.405751768 MHz (Hellwig et al. 1970)", lambda f: f != "0"),
    "casimir": ("pressure pi^2 hbar c/(240 a^4): a power law", lambda f: f == "power"),
    "qed_running": ("alpha increases with Q^2 (QED beta > 0)", lambda f: f == "1"),
    "qcd_running": ("alpha_s decreases with Q^2 (asymptotic freedom)", lambda f: f == "-1"),
    "Ge": ("electron stable, tau_e > 6.6e28 yr (Borexino 2015)", lambda f: f == "0"),
    "hierarchy": ("strong > electromagnetic > weak", lambda f: f == "S>EM>W"),
    "chern": ("Hall conductance quantised in integers of e^2/h", lambda f: f == "integer"),
    "baryon": ("eta = n_B/n_gamma ~ 6.1e-10 != 0 (BBN, CMB)", lambda f: f != "0"),
    "period": ("H (Z = 1) in period 1, Ne (Z = 10) in period 2", lambda f: f == "(1,2)"),
    "Lambda": ("Lambda l_P^2 = %.2e (Planck units), Lambda (hbar/m_e c)^2 = %.2e (electron units); "
               "H0 = 67.4(5) km/s/Mpc, Omega_Lambda = 0.6847(73)" % (Lam * lP ** 2, Lam * lC ** 2), None),
    "H0": ("H0 t_P = %.2e (Planck units)" % (H0 * lP / C["speed of light in vacuum"]), None),
}


def fnum(x):
    x = float(x)
    return f"{x:.6g}" if x != 0 else "0"


def classify(val, ref, sig):
    if ref == 0: return ("D", None, None, None)
    if val == 0: return ("D", None, -1.0, math.inf)
    rel = val / ref - 1
    fac = max(val / ref, ref / val) if val * ref > 0 else math.inf
    ns = (val - ref) / sig if sig else None
    if ns is not None and abs(ns) <= 3: b = "A3s"
    elif abs(rel) <= 0.01: b = "A1%"
    elif fac <= 10: b = "AOM"
    else: b = "D"
    return (b, ns, rel, fac)


BIN_ORDER = {"A3s": 0, "A1%": 1, "AOM": 2, "D": 3}


def evaluate(row):
    out = []
    for i, f in enumerate(row["faces"]):
        r = dict(face=f, i=i)
        if row["kind"] in ("num", "theory"):
            if row["ref"] == "tsirelson":
                val = sp.N(S(f["expr"]), 30)
                ref = 2 * sp.sqrt(2) if i == 0 else 4 * (3 + 2 * sp.sqrt(3)) / 9
                r.update(val=val, ref=sp.N(ref, 30), bin="T" if sp.simplify(S(f["expr"]) - ref) == 0 else "D", ns=None, rel=0.0, fac=1.0)
            else:
                val = float(sp.N(S(f["expr"]), 30))
                ref, sig, _ = REF[row["ref"]]
                b, ns, rel, fac = classify(val, ref, sig)
                r.update(val=val, ref=ref, bin=b, ns=ns, rel=rel, fac=fac)
        elif row["kind"] == "nc":
            r.update(val=float(sp.N(S(f["expr"]), 30)), ref=None, bin="NC", ns=None, rel=None, fac=None)
        else:
            ok = QUAL[row["ref"]][1](f["expr"])
            r.update(val=f["expr"], ref=QUAL[row["ref"]][0], bin="QA" if ok else "QD", ns=None, rel=None, fac=None)
        out.append(r)
    return out


def best_of(res):
    def key(r):
        if r["bin"] in BIN_ORDER:
            return (BIN_ORDER[r["bin"]], abs(r["ns"]) if r["ns"] is not None else (abs(math.log(r["fac"])) if r["fac"] not in (None, math.inf) else 1e9))
        return ({"T": -1, "QA": 0, "QD": 3, "NC": 9}[r["bin"]], 0)
    return min(res, key=key)


def rederive():
    print("[REDERIVE] exact arithmetic from the construction as printed")
    sq2, sq3, sq6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
    a = (79 + 48 * sq2 + 42 * sq3 + 3 * sq6) / 31104
    print(f"  alpha(0): records 79+48+42+3 = {79+48+42+3}; Born count 79*1+48*2+42*3+3*6 = {79+48*2+42*3+3*6}; "
          f"1/alpha = {sp.N(1/a, 15)}")
    # weak: measures and offices
    meas = {(2, 3): 6, (2, 4): 8, (3, 3): 9, (3, 4): 12, (2, 2): 4}
    SW = 2 * meas[(2, 4)] + 2 * meas[(2, 3)]
    SZ2 = SW ** 2 + SW ** 2
    print(f"  W word: two (2,4) and two (2,3) links, S_W = sum of d1 d2 = {SW}; S_Z^2 = S_W^2 + S_U^2 = {SZ2}; "
          f"(S_W/S_Z)^2 = {sp.Rational(SW**2, SZ2)} (164.5, 159.2)")
    print(f"  count face (17/36)^2 = {sp.Rational(17,36)**2}, (14/36)^2 = {sp.Rational(14,36)**2}")
    # 3/8: tr T3^2 / tr Q^2 on the 5 = (1,2)_{1/2} + (3,1)_{-1/3}
    T3 = sp.diag(sp.Rational(1, 2), -sp.Rational(1, 2), 0, 0, 0)
    Y = sp.diag(sp.Rational(1, 2), sp.Rational(1, 2), -sp.Rational(1, 3), -sp.Rational(1, 3), -sp.Rational(1, 3))
    Q = T3 + Y
    print(f"  164.3 dock figure: tr(T3^2)/tr(Q^2) on the 5 of SU(5) = {(T3*T3).trace()/(Q*Q).trace()} "
          f"(the Georgi-Quinn-Weinberg unification value); the hub figure 2/5 is not re-derived (register not specified)")
    # Fourier block: |F|^2 and J
    w = sp.exp(2 * sp.pi * sp.I / 3)
    Fm = sp.Matrix(3, 3, lambda j, k: w ** (j * k) / sp.sqrt(3))
    J = sp.im(sp.expand_complex(Fm[0, 0] * Fm[1, 1] * sp.conjugate(Fm[0, 1]) * sp.conjugate(Fm[1, 0])))
    print(f"  166.1/166.6: order-3 Fourier block, |F_jk|^2 = {sp.simplify(sp.Abs(Fm[1,2])**2)} at all nine; "
          f"J = Im(F11 F22 F12* F21*) = {sp.nsimplify(sp.simplify(J))} = 1/(6 sqrt3), the maximum of |J| over 3x3 unitaries; "
          f"J^2 = {sp.simplify(J**2)}")
    # CHSH and CGLMP
    sx = sp.Matrix([[0, 1], [1, 0]]); sz = sp.Matrix([[1, 0], [0, -1]])
    A0, A1 = sz, sx; B0, B1 = (sz + sx) / sp.sqrt(2), (sz - sx) / sp.sqrt(2)
    kp = lambda X, Y: sp.kronecker_product(X, Y)
    Bop = kp(A0, B0) + kp(A0, B1) + kp(A1, B0) - kp(A1, B1)
    phi = sp.Matrix([1, 0, 0, 1]) / sp.sqrt(2)
    print(f"  36.1: CHSH value at the two-qubit maximally entangled state, standard settings = {sp.simplify((phi.T*Bop*phi)[0])}; "
          f"(12+8 sqrt3)/9 - 4(3+2 sqrt3)/9 = {sp.simplify((12+8*sq3)/9 - 4*(3+2*sq3)/9)}")
    # leptons and proton
    ranks = [27 + 40 * k for k in (1, 2, 3)]
    print(f"  152.1 ranks 27 + 40k = {ranks}; first powers rank*3^d: hub d = (1,2,3) -> {[r*3**d for r, d in zip(ranks, (1,2,3))]}, "
          f"derived cells d = (0,2,3) -> {[r*3**d for r, d in zip(ranks, (0,2,3))]}")
    print(f"  155.1 face dissolution M 3^k/67 = M 3^(k+1)/201 at every M, k: {sp.simplify(sp.Symbol('M')*3**sp.Symbol('k')/67 - sp.Symbol('M')*3**(sp.Symbol('k')+1)/201) == 0}")
    print(f"  186.1 proton 319*9 = {319*9}; 8613 = 3*2871: {8613 == 3*2871}; 2709 = (319-18)*9: {2709 == (319-18)*9}")
    K = lambda r: sp.Rational(sum(x * x for x in r), sum(r) ** 2)
    print(f"  Koide with roots (67,963,3969) = {K((67,963,3969))}, (201,963,3969) = {K((201,963,3969))}, (67,963,4212) = {K((67,963,4212))}")
    print(f"  (2871/67)^2 = {sp.Rational(2871,67)**2} = {sp.N(sp.Rational(2871,67)**2, 12)}; (8613/201)^2 = {sp.Rational(8613,201)**2} (same: hub/hub = dock/dock)")
    print(f"  67 sqrt(m_p/m_e) = {67*math.sqrt(C['proton-electron mass ratio']):.4f}: 2871 is the integer nearest the target, off by "
          f"{2871 - 67*math.sqrt(C['proton-electron mass ratio']):.4f}")
    # hydrogen
    print(f"  143.1 E1: 36*201/(2^30 3^18) = {sp.Rational(36*201, 2**30*3**18)}, 36*67/(2^31 3^17) = {sp.Rational(36*67, 2**31*3**17)}")
    # BH register
    print(f"  209.1 cut register 81*81*36 = {81*81*36} = 2^2 3^10: {81*81*36 == 2**2*3**10}")
    # widths, hierarchy
    print(f"  170.2 (121/16)/(270/49) = {sp.Rational(121,16)/sp.Rational(270,49)}")
    h = [sp.Rational(319, 906) * sp.Rational(25, 6) ** 2, sp.Rational(319, 1036) * sp.Rational(32, 3) ** 2, sp.Rational(906, 1036) * sp.Rational(64, 25) ** 2]
    print(f"  211.1 class ratios (Born-face ratio x register ratio^2): {h} -> (2,3) > (2,4) > (3,3) forced by the registers")
    print(f"  202.1 lambda_hub = {sp.Rational(46656001681, 182250000*6561**2)}; 200.6 M_H/M_W = {sp.N(sp.sqrt(sp.Rational(46656001681,182250000)/(sp.Rational(2857690241,3645000)+sp.Rational(98,50625)*sq2)), 10)}")
    print()


def main():
    tsv = sys.argv[sys.argv.index("--tsv") + 1] if "--tsv" in sys.argv else None
    if "--pdg-mcd" in sys.argv:
        n = check_mcd(sys.argv[sys.argv.index("--pdg-mcd") + 1])
        print(f"[references] PDG 2024 transcription checked against mass_width_2024.mcd ({n} ids parsed): OK")
    else:
        print("[references] PDG 2024 MC values transcribed (not re-checked this run; pass --pdg-mcd)")
    print("[references] CODATA 2022 values asserted against scipy.constants", scipy.__version__)
    print(f"  derived: v = {v:.4f} GeV; on-shell sin^2 = {1-c2:.6f}({sc2*1e6:.0f}e-6); Koide(PDG 2024) = {kq:.8f}({kq_s*1e8:.0f}e-8); "
          f"lambda_tree = {lamH:.5f}; f_dm = {fdm:.4f}({sfdm*1e4:.0f}e-4)\n")
    rederive()
    lines, allres = [], []
    print("[SCORECARD] per face: value as the measured quantity; reference; deviation")
    hdr = f"{'id':6s} {'stmt':7s} {'gr':3s} {'idt':3s} {'der':3s} {'face':58s} {'value':>13s} {'ref':>13s} {'n_sigma':>10s} {'rel':>10s} {'factor':>9s} bin"
    print(hdr)
    for row in D.ROWS:
        res = evaluate(row); allres.append((row, res))
        best = best_of(res)
        for r in res:
            f = r["face"]; mark = ("*" if r is best else " ") + ("p" if f["primary"] else " ")
            g = GRADES.get(f["stmt"], "?")[:3]
            val = fnum(r["val"]) if not isinstance(r["val"], str) else r["val"]
            ref = fnum(r["ref"]) if isinstance(r["ref"], (int, float, sp.Basic)) else "(qual.)" if r["ref"] else "-"
            ns = f"{r['ns']:.3g}" if r.get("ns") is not None else "-"
            rel = f"{r['rel']:.3e}" if r.get("rel") not in (None,) else "-"
            fac = (f"{r['fac']:.3g}" if r["fac"] != math.inf else "inf") if r.get("fac") is not None else "-"
            print(f"{row['id']:6s} {f['stmt']:7s} {g:3s} {row['ident']:3s} {row['derived']:3s} {mark}{f['label'][:56]:56s} {val:>13s} {ref:>13s} {ns:>10s} {rel:>10s} {fac:>9s} {r['bin']}")
            lines.append([row["id"], row["qty"], f["stmt"], GRADES.get(f["stmt"], "?"), row["ident"], row["derived"], row["kind"], f["label"],
                          f["printed"], f["how"], str(val), str(ref), ns, rel, fac, r["bin"], "best" if r is best else "", "primary" if f["primary"] else "",
                          f'{row["ident_quote"][0]}: {row["ident_quote"][1]}', (f'{row["verdict"][0]}: {row["verdict"][1]}' if row.get("verdict") else "")])
    print("\n[IDENTIFICATION] per row: code, the verbatim identification sentence, the source's own verdict (both checked by scorecard_extract.py)")
    for row in D.ROWS:
        st, q = row["ident_quote"]
        vd = f' | verdict ({row["verdict"][0]}): "{row["verdict"][1]}"' if row.get("verdict") else ""
        print(f'  {row["id"]} [{row["ident"]}] ({st}) "{q}"{vd}')
    # ------------------------------------------------------------ summary (task 3)
    print("\n[SUMMARY]")
    nrows, nfaces = len(D.ROWS), sum(len(r["faces"]) for r in D.ROWS)
    print(f"  named quantities (rows): {nrows}; printed numbers attached to them (faces): {nfaces}")
    from collections import Counter
    print(f"  identification: {dict(Counter(r['ident'] for r in D.ROWS))}  (E explicit normalisation/unit, N name/nearest only, P pairing ours)")
    print(f"  derivation: {dict(Counter(r['derived'] for r in D.ROWS))}  (R re-derived here, A arithmetic of printed inputs, T target-like/candidate, N stated only)")
    print(f"  kinds: {dict(Counter(r['kind'] for r in D.ROWS))}")
    best = {row["id"]: best_of(res) for row, res in allres}
    prim = {row["id"]: next((r for r in res if r["face"]["primary"]), None) for row, res in allres}
    bb = Counter(best[r["id"]]["bin"] for r in D.ROWS)
    print(f"  rows by best face: {dict(sorted(bb.items()))}")
    numrows = [r for r in D.ROWS if r["kind"] == "num"]
    print(f"  numerical rows: {len(numrows)}; by best face {dict(Counter(best[r['id']]['bin'] for r in numrows))}")
    print(f"  numerical rows with an explicit identification (E): {sum(1 for r in numrows if r['ident']=='E')}; "
          f"by best face {dict(Counter(best[r['id']]['bin'] for r in numrows if r['ident']=='E'))}")
    pp = [r for r in numrows if prim[r["id"]] is not None]
    print(f"  numerical rows with a face the source adopts: {len(pp)}; by that face {dict(Counter(prim[r['id']]['bin'] for r in pp))}")
    fb = Counter(x["bin"] for row, res in allres if row["kind"] == "num" for x in res)
    print(f"  numerical faces: {sum(fb.values())}; {dict(sorted(fb.items()))}")
    # compact table (task 3)
    allf = [(row, x) for row, res in allres for x in res]
    def rc(pred): return sum(1 for r in D.ROWS if pred(r, best[r["id"]]))
    def fc(pred): return sum(1 for row, x in allf if pred(row, x))
    choice = {row["id"]: (len(res) > 1 and (prim[row["id"]] is None or prim[row["id"]] is not best[row["id"]])) for row, res in allres}
    print("\n  TABLE (rows classified at their best face; faces counted individually)")
    print(f"    {'':54s} {'rows':>5s} {'faces':>6s}")
    T = [("named physical quantities / printed numbers", lambda r, b: True, lambda row, x: True),
         ("  source states a normalisation or unit (E)", lambda r, b: r["ident"] == "E", lambda row, x: row["ident"] == "E"),
         ("  name / nearest object only (N), pairing ours (P)", lambda r, b: r["ident"] != "E", lambda row, x: row["ident"] != "E"),
         ("numerical comparison possible", lambda r, b: r["kind"] == "num", lambda row, x: row["kind"] == "num"),
         ("  within 3 sigma", lambda r, b: b["bin"] == "A3s", lambda row, x: x["bin"] == "A3s"),
         ("  within 1% (beyond 3 sigma)", lambda r, b: b["bin"] == "A1%", lambda row, x: x["bin"] == "A1%"),
         ("  within a factor 10 (beyond 1%)", lambda r, b: b["bin"] == "AOM", lambda row, x: x["bin"] == "AOM"),
         ("  beyond a factor 10, or zero against nonzero", lambda r, b: b["bin"] == "D", lambda row, x: x["bin"] == "D"),
         ("theorem of standard QM reproduced", lambda r, b: b["bin"] == "T", lambda row, x: x["bin"] == "T"),
         ("integer / sign / form: agrees", lambda r, b: b["bin"] == "QA", lambda row, x: x["bin"] == "QA"),
         ("integer / sign / form / zero: disagrees", lambda r, b: b["bin"] == "QD", lambda row, x: x["bin"] == "QD"),
         ("not comparable (unit unmapped)", lambda r, b: b["bin"] == "NC", lambda row, x: x["bin"] == "NC")]
    for lab, pr, pf in T:
        print(f"    {lab:54s} {rc(pr):5d} {fc(pf):6d}")
    agr = [r for r in D.ROWS if best[r["id"]]["bin"] in ("A3s", "A1%")]
    print(f"    numerical agreements within 1%: {len(agr)} rows; of these resting on a choice among printed faces: "
          f"{sum(1 for r in agr if choice[r['id']])} {[r['id'] for r in agr if choice[r['id']]]}")
    qa = [r for r in D.ROWS if best[r["id"]]["bin"] in ("QA", "T")]
    print(f"    structural agreements (QA, T): {len(qa)} rows; resting on a choice among printed faces: {sum(1 for r in qa if choice[r['id']])}")
    print("\n  agreements (best face A3s, A1%, T or QA), with the choice behind them:")
    for row, res in allres:
        b = best[row["id"]]
        if b["bin"] in ("A3s", "A1%", "T", "QA"):
            choice = len(res) > 1 and (prim[row["id"]] is None or prim[row["id"]] is not b)
            ext = f"rel {b['rel']:+.3e}" if isinstance(b.get("rel"), float) and b["bin"] in ("A3s", "A1%") else ""
            ext += f", {b['ns']:+.3g} sigma" if b.get("ns") is not None else ""
            print(f"    {row['id']} {row['qty'][:62]:62s} [{b['bin']}] face '{b['face']['label'][:60]}' {ext}; "
                  f"{len(res)} face(s), choice: {'yes' if choice else 'no'}; ident {row['ident']}; derived {row['derived']}")
    if tsv:
        import csv
        with open(tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["id", "quantity", "statement", "grade", "ident", "derived", "kind", "face", "printed", "conversion", "value",
                        "reference", "n_sigma", "rel", "factor", "bin", "best", "primary", "identification_sentence", "source_verdict"])
            w.writerows(lines)
    print("\n  references used:")
    for k in sorted({r["ref"] for r in D.ROWS}):
        if k in REF: print(f"    {k}: {REF[k][0]:.10g} +- {REF[k][1] if REF[k][1] is not None else '-'}  [{REF[k][2]}]")
        else: print(f"    {k}: {QUAL[k][0]}")


def load_grades():
    import csv
    if "--src" not in sys.argv: return {}
    src = sys.argv[sys.argv.index("--src") + 1]
    return {r["number"]: r["grade"] for r in csv.DictReader(open(pathlib.Path(src)/"index"/"statements.tsv", encoding="utf-8", newline=""), delimiter="\t")}


GRADES = load_grades()

if __name__ == "__main__":
    main()
