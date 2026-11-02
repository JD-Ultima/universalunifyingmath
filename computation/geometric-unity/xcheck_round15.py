"""Independent referee checks for S, REVIEW12-GU-S (2026-11-01).

Default/--quick: derive the rational coefficient identities independently with
SymPy; guard the constant-kappa exception and the degree drop at infinity;
check the local ramification/pole calculation used in the continuation proof.
These symbolic checks support, and do not replace, the algebraic-cover proof.

--geometry additionally studies one rank-one-Ricci jet in the unresolved S15
class using only referee xcheck_round10's independent metric/curvature engine.
No imports from S, Q, N, M, L, J or F_jets. The geometry is an exact example,
not a classification of the rank-one-Ricci sector.

Run: python3 -I xcheck_round15.py --quick
     python3 -I xcheck_round15.py --geometry
"""
import argparse
import os
import pathlib
import sys
import time
from fractions import Fraction as Fr

for _var in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_var, "1")

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def quick():
    print("X0. Rational coefficients and exceptional amplitudes", flush=True)
    tau, A, z, B = sp.symbols("tau A z B")
    for ep in (1, -1):
        D = 1 + ep * A + 544 * (1 - ep * A) * tau
        Q = 1 + 768 * tau + (121856 - 147456 * ep * A) * tau**2
        N = (7 - 3 * ep * A) + (7616 - 3008 * ep * A - 768 * A**2) * tau \
            + (2071552 - 3226624 * ep * A + 1155072 * A**2) * tau**2
        K = sp.cancel(-N / (18 * D**2))
        dn = sp.together(sp.diff(K, tau)).as_numer_denom()[0]
        common = sp.gcd_list(sp.Poly(dn, tau).all_coeffs())
        check(f"epsilon={ep}: kappa constant only at A=0 or A={ep}*11/6",
              sp.degree(sp.cancel(common / (A * (6 * A - 11 * ep))), A) == 0)
        if ep == 1:
            kc = sp.cancel(K.subs(A, sp.Rational(11, 6)))
            check("A=11/6: kappa=-3/289 and its algebraic fibre equation degenerates",
                  kc == -sp.Rational(3, 289)
                  and sp.expand((18 * kc * D**2 + N).subs(A, sp.Rational(11, 6))) == 0)
            lam = sp.simplify((1 - 1 / (16 * kc)) / 5)
            check("constant kappa persistence needs lambda=337/240, outside (14,0)",
                  lam == sp.Rational(337, 240) and lam > sp.Rational(1, 4))
            eta_inf = sp.cancel(-576 * D.subs(tau, 1 / z) / (z**2 * Q.subs(tau, 1 / z)))
            special = sp.cancel(eta_inf.subs(A, sp.Rational(119, 144)))
            check("A=119/144: Q2 is linear and eta has a double pole at infinity",
                  sp.degree(Q.subs(A, sp.Rational(119, 144)), tau) == 1
                  and sp.limit(z**2 * special, z, 0) == -sp.Rational(425, 6))
            check("A=1: eta is regular at infinity (no blanket log-pole claim)",
                  sp.limit(eta_inf.subs(A, 1), z, 0) == sp.Rational(9, 200))

    print("X1. Linear-mu trace obstruction", flush=True)
    r, Nv, kap, cv = sp.symbols("r N kappa cV", nonzero=True)
    sig = (-1 / (224 * Nv * r) - kap) / (224 * Nv)
    P0 = sp.expand(sp.Rational(7, 18) * (kap * sig + 312 * Nv * sig**2))
    coefficient = sp.expand(P0 * r**2).subs(r, 0)
    check("P0(r) has nonzero r^-2 coefficient 13/(269746176 N^3)",
          sp.simplify(coefficient - sp.Rational(13, 269746176) / Nv**3) == 0)
    lam = sp.Rational(13, 56)
    c = 4 * cv * (1 - 4 * lam)
    b = (1 - 5 * lam) / (cv * (1 - 4 * lam))
    mu1 = -sp.Rational(7, 18)
    check("lambda=13/56: d(log r)/dtheta=1/2 on the continued radial orbit",
          sp.simplify(2 * mu1 * b * c) == sp.Rational(1, 2))
    q = sp.symbols("q", nonzero=True)
    check("a Laurent polynomial with this r^-2 term cannot be constant along r=r0*q",
          sp.Poly(sp.expand(r**2 * P0).subs(r, r * q), q).nth(0) == coefficient)

    print("X2. The pole calculation in the finite-cover proof (local models)", flush=True)
    u = sp.symbols("u", nonzero=True)
    for order in (1, 2, 3):
        for ramification in (1, 2, 3, 4):
            pulled = sp.cancel((u**ramification)**(-order) * sp.diff(u**ramification, u))
            exponent = ramification * (1 - order) - 1
            check(f"pole order {order}, ramification {ramification}: pullback exponent {exponent}<0",
                  pulled == ramification * u**exponent and exponent < 0)
    print("These local models illustrate the proved general order e(m-1)+1; they do not prove the global continuation lemma.", flush=True)


def geometry():
    sys.path.insert(0, str(HERE))
    import xcheck_round10 as X
    from sympy.polys.matrices import DomainMatrix

    print("X3. One rank-one-Ricci jet in the unresolved S15 sector (exact example)", flush=True)
    g4 = X.gmat((1, 1, 1, 1))
    # R = Schouten wedge g, with Ric=diag(1,0,0,0), scalar=1, Weyl=0.
    schouten = [[(Fr(int(i == j == 0)) - g4[i][j] / 6) / 2 for j in range(4)] for i in range(4)]
    R = X.kn(schouten, g4)
    base = X.weyl_base(R, (1, 1, 1, 1))
    bg = X.base_geometry(base)
    check("base jet has Ric=diag(1,0,0,0)",
          bg["Ric"] == [[Fr(int(i == j == 0)) for j in range(4)] for i in range(4)])

    def relaxed(D):
        m = 14
        gh, gih, ric, nr, ra = D["g"], D["gi"], D["Ric"], D["nRic"], D["RA"]
        rh = X.mmul(gih, ric)
        nh = {a: X.mmul(gih, nr[a]) for a in range(m)}
        rows, rhs = [], []
        for a in range(m):
            for b_ in range(a + 1, m):
                for d in range(m):
                    row = [Fr(0)] * (3 * m)
                    for cc in range(m):
                        row[cc] = ra[a][b_][d][cc]
                        row[m + cc] = -(int(a == cc) * rh[d][b_] - int(b_ == cc) * rh[d][a])
                        row[2 * m + cc] = -(gh[a][cc] * int(d == b_) - gh[b_][cc] * int(d == a))
                    rows.append(row)
                    rhs.append(nh[a][d][b_] - nh[b_][d][a])
        return rows, rhs

    def rank(rows):
        return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(sp.QQ).rank()

    outcomes = {}
    for tag, ydiag in (("section", (1, 1, 1, 1)), ("off-section", (2, 3, 4, 5))):
        yp = [[Fr(ydiag[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]
        D = X.point_data(base, yp, Fr(1, 5), Fr(2), Fr(3), riem=True)
        rows, rhs = relaxed(D)
        pair = rank(rows), rank([row + [val] for row, val in zip(rows, rhs)])
        outcomes[tag] = pair
        print(f"  {tag}: rank(A),rank([A|b])={pair}", flush=True)
    check("rank-one control: relaxed coefficient rank at the section is 41", outcomes["section"][0] == 41)
    check("this off-section rank-one-Ricci jet is inconsistent", outcomes["off-section"][1] > outcomes["off-section"][0])
    print("By S4 the full-rank inconsistent point excludes every fibre-open consistent patch for this specified metric jet; arbitrary Weyl tensors and remaining connection jets are not covered.", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--quick", action="store_true", help="coefficient and local pole checks only (the default)")
    ap.add_argument("--geometry", action="store_true", help="also evaluate the exact rank-one-Ricci example")
    args = ap.parse_args()
    start = time.time()
    quick()
    if args.geometry:
        geometry()
    print(f"\n{len(CHECKS)} checks, all True; mode={'geometry' if args.geometry else 'quick'}")
    print(f"runtime {time.time() - start:.1f} s")


if __name__ == "__main__":
    main()
