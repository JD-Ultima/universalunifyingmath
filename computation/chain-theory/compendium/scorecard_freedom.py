#!/usr/bin/env python3
"""Declared-family sensitivity for scorecard agreements other than alpha(0) (round 2, task 3).

Usage: python3 -I scorecard_freedom.py

Modelled on F-CT-4 of the alpha audit (formalization/chain-theory/alpha/audit.md): for each agreement,
(A) the source's own printed alternatives, (B) a declared rule space built only from choices the
source itself prints or names, with the number of variants, the number of distinct values, the hits
at the achieved precision and density-times-window reference counts,
and (C) the list of discrete choices behind the agreeing face.

  W  weak angle: sin^2(theta_W) = (17/36)^2 = 289/1296 (164.2, OPEN), on-shell -0.83 sigma
  L  lepton mass ratios (963/67)^2, (3969/67)^2 and Koide at the dock face (186.1, 157.1)
  P  proton-electron mass ratio (2871/67)^2 (186.1 'derived cell' = 21.1/175.4 'candidate cell')
Exact arithmetic (sympy) for the re-derivation of 17/36; floating point (mpmath, 30 digits) for
densities. Reference values as in scorecard.py (CODATA 2022 from scipy.constants, PDG 2024).
Reference counts use illustrative uniform-value or locally uniform log-value assumptions.
Integer reference coverage assumes a uniform phase modulo the integer spacing. These are
not calibrated p-values for the deterministic constructions or evidence of a search history.
"""
import itertools, math
import sympy as sp
import mpmath as mp
import scipy.constants as sc

mp.mp.dps = 30
PC = sc.physical_constants
MW, sMW, MZ, sMZ = 80.3692, 0.0133, 91.1880, 0.0020           # PDG 2024
S2_OS = 1 - (MW / MZ) ** 2
S2_OS_S = (1 - S2_OS) * math.hypot(2 * sMW / MW, 2 * sMZ / MZ)
SCHEMES = {  # sin^2 theta_W in four conventions; the last two transcribed (PDG 2024 electroweak review), used as targets only
    "on-shell, PDG 2024 masses": S2_OS,
    "on-shell, CODATA 2022": PC["weak mixing angle"][0],
    "MS-bar at M_Z (~0.23122)": 0.23122,
    "effective leptonic (~0.23153)": 0.23153,
}
MMU, MTAU, MP = PC["muon-electron mass ratio"][0], PC["tau-electron mass ratio"][0], PC["proton-electron mass ratio"][0]


def local_density(values, target, half_width=0.5):
    """distinct values per unit of ln within |ln(v/target)| < half_width"""
    vs = sorted({round(float(v), 12) for v in values if v > 0})
    n = sum(1 for v in vs if abs(math.log(v / target)) < half_width)
    return n / (2 * half_width), n


# ------------------------------------------------------------------ W: weak angle
def weyl(d):
    w = sp.exp(2 * sp.pi * sp.I / d)
    X = sp.Matrix(d, d, lambda i, j: 1 if i == (j + 1) % d else 0)
    Z = sp.diag(*[w ** k for k in range(d)])
    return [(a, b, X ** a * Z ** b) for a in range(d) for b in range(d)]


def gellmann_pauli(d):
    if d == 2:
        return [(0, 0, sp.eye(2)), (1, 0, sp.Matrix([[0, 1], [1, 0]])), (1, 1, sp.Matrix([[0, -sp.I], [sp.I, 0]])), (1, 2, sp.Matrix([[1, 0], [0, -1]]))]
    out = [(0, 0, sp.eye(d))]
    for i in range(d):
        for j in range(i + 1, d):
            Eij = sp.zeros(d); Eij[i, j] = 1; Eji = Eij.T
            out += [(1, 1, Eij + Eji), (1, 2, -sp.I * Eij + sp.I * Eji)]
    for k in range(1, d):
        D = sp.zeros(d)
        for i in range(k): D[i, i] = 1
        D[k, k] = -k
        out.append((1, 3, D / sp.sqrt(k * (k + 1) / 2)))
    return out


def pair_ket(d1, d2):
    m = min(d1, d2)
    v = sp.zeros(d1 * d2, 1)
    for i in range(m): v[i * d2 + i] = 1
    return v / sp.sqrt(m)


def count_faces(d1, d2, basis):
    L1, L2 = basis(d1), basis(d2); phi = pair_ket(d1, d2)
    nz_whole = nz_cut = 0
    for a1, b1, A in L1:
        for a2, b2, B in L2:
            val = sp.nsimplify(sp.simplify((phi.H * sp.kronecker_product(A, B) * phi)[0]))
            if val != 0:
                nz_whole += 1
                if (a1, b1) != (0, 0) and (a2, b2) != (0, 0): nz_cut += 1
    return nz_whole, nz_cut, len(L1) * len(L2), (len(L1) - 1) * (len(L2) - 1)


def part_W():
    print("[W] weak angle: 289/1296 = (17/36)^2 (164.2, OPEN) against sin^2(theta_W)")
    nzw, nzc, tot, cut = count_faces(3, 2, weyl)
    print(f"  re-derivation: (3,2) pair ket on its 2-dim common support, Weyl letters X^a Z^b (x) X^c Z^d: "
          f"{tot} pair letters, {nzw} with nonzero amplitude (whole), {nzc} with both letters non-identity (cut); "
          f"printed 17 and 14 -> {'reproduced' if (nzw, nzc) == (17, 14) else 'NOT reproduced'}")
    g = count_faces(3, 2, gellmann_pauli)
    print(f"  same count with Gell-Mann x Pauli letters: whole {g[0]}, cut {g[1]} (so the 17 depends on the letter basis)")
    face = sp.Rational(289, 1296)
    print("  deviation of 289/1296 = %.6f by scheme:" % float(face))
    for k, t in SCHEMES.items():
        extra = f", {(float(face) - t) / S2_OS_S:+.2f} sigma" if k.startswith("on-shell, PDG") else (
            f", {(float(face) - t) / PC['weak mixing angle'][2]:+.2f} sigma" if "CODATA" in k else "")
        print(f"    {k:32s} {t:.6f}: rel {float(face) / t - 1:+.4%}{extra}")
    # (A) printed faces
    printed = [sp.Rational(1, 2), sp.Rational(17, 36), sp.Rational(14, 36), face, sp.Rational(49, 324), sp.Rational(3, 8), sp.Rational(2, 5),
               sp.Rational(2, 3), sp.Rational(1, 3), sp.Rational(1, 2), sp.Rational(1, 2), sp.Integer(0)]
    d23 = sp.Rational(1520, 8883) + sp.Rational(16000, 26649) * sp.sqrt(2) - sp.Rational(8800, 8883) * sp.sqrt(3) + sp.Rational(8240, 8883) * sp.sqrt(6)
    printed.append(d23 / (1 + d23))
    dist = sorted({float(x) for x in printed})
    delta = abs(float(face) - S2_OS)
    print(f"  (A) printed faces for the angle (as sin^2): {len(printed)} ({len(dist)} distinct values in [0, 1])")
    for k, t in SCHEMES.items():
        hits = [x for x in dist if abs(x - t) <= delta]
        print(f"      within |delta| = {delta:.2e} of {k}: {len(hits)}")
    print(f"      illustrative uniform-value reference, values uniform on [0,1]: {len(dist)} values x {len(SCHEMES)} schemes x 2 delta = "
          f"{len(dist) * len(SCHEMES) * 2 * delta:.3f} expected matches (counts summed over schemes)")
    # (B) declared space: the count-face family at the six genres
    vals = []
    for (d1, d2) in [(2, 2), (3, 2), (4, 2), (3, 3), (4, 3), (4, 4)]:
        for basis in (weyl, gellmann_pauli):
            nzw, nzc, tot, cut = count_faces(d1, d2, basis)
            for num in (nzw, nzc):
                for den in (tot, cut):
                    for e in (1, 2):
                        x = sp.Rational(num, den) ** e
                        if x <= 1:
                            vals += [x, 1 - x]                 # read as sin^2 or as cos^2
    dvals = sorted({float(x) for x in vals})
    print(f"  (B) declared space: six genres x {{Weyl, Gell-Mann/Pauli}} letters x {{whole, cut}} numerator x {{all letters, cut letters}} "
          f"denominator x {{count face, its Born square}} x {{sin^2, cos^2}}: {len(vals)} variants, {len(dvals)} distinct values in [0,1]")
    for k, t in SCHEMES.items():
        hits = [x for x in dvals if abs(x - t) <= delta]
        dens = sum(1 for x in dvals if abs(x - t) <= 0.1) / 0.2
        print(f"      {k:32s} hits within {delta:.2e}: {len(hits)} {[round(h, 6) for h in hits]}; "
              f"local density {dens:.1f}/unit -> density-window reference count {dens * 2 * delta:.4f}")
    tot_exp = sum(sum(1 for x in dvals if abs(x - t) <= 0.1) / 0.2 * 2 * delta for t in SCHEMES.values())
    print(f"      density-window reference counts summed over the four schemes: {tot_exp:.3f}")
    print("  (C) discrete choices behind the agreeing face: (1) the Born square of the count face rather than the count face the")
    print("      bridge names (the source refuses this step: 'refused at its first step'); (2) whole (17) rather than cut (14); (3) Weyl")
    print("      letters rather than another letter basis (Gell-Mann/Pauli gives 5); (4) sin^2 rather than cos^2; (5) the on-shell")
    print("      scheme (the MS-bar value is 3.6% away). The source's own CONDITIONAL value for the angle is cos^2 = 1/2 (164.5).\n")


# ------------------------------------------------------------------ L, P: lepton and proton cells
def part_L():
    print("[L] lepton mass ratios from first powers: derived cells (67, 963, 3969) with mass = S^2")
    tgt = {"mu/e": MMU, "tau/e": MTAU}
    ach = {"mu/e": abs((963 / 67) ** 2 / MMU - 1), "tau/e": abs((3969 / 67) ** 2 / MTAU - 1)}
    print(f"  achieved: (963/67)^2 = {(963/67)**2:.4f} vs {MMU} ({(963/67)**2/MMU-1:+.3e}); (3969/67)^2 = {(3969/67)**2:.3f} vs {MTAU} ({(3969/67)**2/MTAU-1:+.3e})")
    print("  (A) printed electron cells: 67 (dock-struck, dressing 1), 201 (hub, dressing 3); muon 963, 621 (stop clause, 21.1);")
    print("      tau 3969, 432 (stop clause), 4185 (tau sink), 4212 (generation-step web). Ratios (S_l/S_e)^p:")
    for se in (67, 201):
        for p in (1, 2):
            print(f"      e = {se:3d}, p = {p}: mu {[f'{(x/se)**p:.4g}' for x in (963, 621)]}  tau {[f'{(x/se)**p:.4g}' for x in (3969, 432, 4185, 4212)]}")
    # (B) declared space: struck or unstruck measure, each hub in-leg passed or unpassed (dressing 3 per unpassed in-leg),
    #     mass exponent p, plus the printed stop-clause cells
    M = {"e": (79, 12, 1), "mu": (117, 10, 2), "tau": (155, 8, 3)}   # unstruck measure, printed strike, hub in-legs
    cells = {}
    for k, (m, s, h) in M.items():
        cells[k] = sorted({(m - st) * 3 ** dd for st in (0, s) for dd in range(h + 1)})
    cells["mu"] = sorted(set(cells["mu"]) | {621}); cells["tau"] = sorted(set(cells["tau"]) | {432})
    print(f"  (B) declared space: S = (measure - strike or 0) x 3^(unpassed in-legs), each hub in-leg passed or not; mass = S^p, p in {{1,2}};")
    print(f"      electron cells {cells['e']}; muon {cells['mu']}; tau {cells['tau']}")
    for lab, key in (("mu/e", "mu"), ("tau/e", "tau")):
        vals = [(x / e) ** p for x in cells[key] for e in cells["e"] for p in (1, 2)]
        dens, n = local_density(vals, tgt[lab])
        hits = sorted({round(v, 6) for v in vals if abs(v / tgt[lab] - 1) <= ach[lab]})
        print(f"      {lab}: {len(vals)} variants, {len({round(v, 9) for v in vals})} distinct; hits within {ach[lab]:.2e}: {hits}; "
              f"local density {dens:.2f}/unit ln -> density-window reference count {dens * 2 * ach[lab]:.4f}")
    joint = [(e, mu, ta, p) for e in cells["e"] for mu in cells["mu"] for ta in cells["tau"] for p in (1, 2)
             if abs((mu / e) ** p / MMU - 1) <= 0.01 and abs((ta / e) ** p / MTAU - 1) <= 0.01]
    ntot = len(cells["e"]) * len(cells["mu"]) * len(cells["tau"]) * 2
    print(f"      joint: mu/e and tau/e both within 1% at {len(joint)} of {ntot} (e, mu, tau, p) combinations: {joint}")
    print("  (C) discrete choices behind the agreeing faces: (1) the electron at its dock-struck cell 67 with its one in-leg 'passed'")
    print("      (dressing 1) while the muon's two and the tau's three in-legs are unpassed (dressings 9, 27); the hub face, with every")
    print("      in-leg unpassed (dressings 3, 9, 27), gives the source's printed ratios 22.95 and 389.9 (176.2); (2) mass = S^2 (185.9)")
    print("      rather than office = mass^2 (181.5, 200.6 bridges); (3) the strikes 12, 10, 8 (the docks of each roster web); (4) the")
    print("      roster webs themselves (one, two, three hub assemblies of measure 40 on a core of measure 27). The source states the")
    print("      ratios are 'compared at the direct reading of the comparison table, and this compendium displays them nowhere' (176.1);")
    print("      186.1 tests its segment rule against these cells ('the derived cells') and reports 8 of its 36 readings meet them.")
    kq = lambda r: sum(x * x for x in r) / sum(r) ** 2
    print(f"  Koide at the dock face {kq((67, 963, 3969)):.6f} and hub face {kq((201, 963, 3969)):.6f} follow from the same cells (not independent).\n")


def part_P():
    print("[P] proton-electron mass ratio: (2871/67)^2 with 2871 = 319 x 9")
    ach = abs((2871 / 67) ** 2 / MP - 1)
    s_t = 67 * math.sqrt(MP)
    print(f"  achieved: (2871/67)^2 = {(2871/67)**2:.6f} vs {MP} (rel {(2871/67)**2/MP-1:+.3e}); target cell 67 sqrt(m_p/m_e) = {s_t:.4f}")
    printed = {2871: "186.1 derived cell / 21.1, 175.4 'candidate cell'", 8613: "hub numerator; merge across the cap", 2709: "strike at the capped entry",
               2727: "175.1 seed, hub", 101: "175.1 seed, dock", 2187: "175.3 first catalogued sink, hub", 25839: "175.4 counted reading", 729: "180.1 formula"}
    print("  (A) printed proton cells over e = 67 (and over 201 for the hub-spelled ones), mass = S^2:")
    for s_, lab in printed.items():
        r67, r201 = (s_ / 67) ** 2, (s_ / 201) ** 2
        print(f"      {s_:6d} {lab:46s} (S/67)^2 = {r67:12.4f} ({r67/MP-1:+.2e})   (S/201)^2 = {r201:10.4f} ({r201/MP-1:+.2e})")
    # (B) declared space: entry-segment records = a genre Born count or a printed proton measure, strike 0 or 18, dressing 3^a 4^b
    R = [96, 319, 906, 1036, 2053, 5100, 81, 85, 101, 125, 150]
    dress = sorted({3 ** a * 4 ** b for a in range(5) for b in range(5) if a + b <= 4})
    cells = sorted({(r - s) * d for r in R for s in (0, 18) for d in dress})
    vals = [(c / e) ** p for c in cells for e in (67, 201) for p in (1, 2)]
    dens, n = local_density(vals, MP)
    hits = sorted({(c, e, p) for c in cells for e in (67, 201) for p in (1, 2) if abs((c / e) ** p / MP - 1) <= ach})
    print(f"  (B) declared space: S_p = (records - strike) x dressing, records in the six genre Born counts {R[:6]} or the printed proton")
    print(f"      measures {R[6:]}, strike in {{0, 18}}, dressing 3^a 4^b with a + b <= 4 (198.3: three per qutrit in-leg, four per ququart");
    print(f"      in-leg), e in {{67, 201}}, p in {{1, 2}}: {len(cells)} cells, {len(vals)} variants")
    print(f"      hits within {ach:.2e}: {hits}; local density {dens:.1f} distinct values per unit ln -> density-window reference count {dens * 2 * ach:.4f}")
    w = sum(1 for v in vals if abs(v / MP - 1) <= 0.01)
    print(f"      variants within 1%: {w}")
    print(f"  integer reference: nearest integer to {s_t:.4f} is 2871, off by {2871 - s_t:+.4f}; under a uniform target phase modulo the spacing,")
    print(f"  coverage this close is {2 * abs(2871 - s_t):.3f} for spacing 1 and {2 * abs(2871 - s_t) / 9:.4f} for spacing 9 (dressing 9).")
    print("  (C) discrete choices: (1) entry segment = the capped (2,3) link at its cap 319 with nothing struck (striking gives 2709,")
    print("      -11%); (2) dressing 9 from two unpassed in-legs (27 gives 8613 over 67, x9); (3) 'a capped link beyond the first cap reads")
    print("      nothing' (merging gives 8613); (4) electron at 67 (with 2871) or at 201 (with 8613): same ratio; (5) mass = S^2. The source")
    print("      names 2871 'the candidate cell' (20.1, 175.4), finds no sink that reads it ('the candidate is refuted at the bounded")
    print("      closures', 175.3) and states that the proton in the electron unit 'stands outside the comparison by orders' (176.1).\n")


if __name__ == "__main__":
    print("[STATISTICAL SCOPE] density-window counts are illustrative uniform/local-uniform references; integer coverage assumes a uniform phase.")
    print("They are not calibrated p-values, do not supply a selection model, and establish no construction search history or author intent.\n")
    part_W(); part_L(); part_P()
    print("[reference] alpha(0), F-CT-4: unique 1-sigma match in 46 080 rules; illustrative uniform-target expected count 4.4e-6, not a calibrated p-value.")
