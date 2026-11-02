"""Exact obstruction for a declared family of normal quadratic base metrics.

The family is g_(rho,v)(x)=I-rho R0(v)_{iajb} x^a x^b/3, |v|=1,
R0(v)=((v tensor v-I/6)/2) wedge I. All required connection jets are fixed by
this complete metric polynomial. This does not cover arbitrary metrics with
rank-one Ricci curvature.

Only xcheck_round10, an earlier independent referee metric/curvature engine,
is imported. The 42-unknown relaxed system and its exact augmented minors are
assembled here. The two sign representatives and two nonunit scale controls
are checked with rational arithmetic; the parameter-family extension rests on
the coordinate/fibre/homothety maps proved in U-rank-one-frontier.md.
"""
import hashlib
import os
import pathlib
import sys
import time
from fractions import Fraction as Fr

for _key in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_key, "1")

import sympy as sp
from sympy.polys.matrices import DomainMatrix

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import xcheck_round10 as X

CHECKS = []
ETA = (1, 1, 1, 1)
LAMBDA = Fr(1, 5)
YSTD = (2, 3, 4, 5)


def check(label, ok):
    value = bool(ok)
    CHECKS.append(value)
    print(f"[{len(CHECKS):02d}] {label}: {value}", flush=True)
    assert value, label


def diagonal(values):
    return [[Fr(values[i]) if i == j else Fr(0) for j in range(4)] for i in range(4)]


def base_metric(rho):
    """Complete normal quadratic metric, not just its Ricci tensor."""
    identity = X.gmat(ETA)
    schouten = [[(Fr(int(i == j == 0)) - identity[i][j] / 6) / 2 for j in range(4)] for i in range(4)]
    R0 = X.kn(schouten, identity)
    return X.weyl_base({key: Fr(rho) * val for key, val in R0.items()}, ETA)


def point(rho, cH, cV, yscale=Fr(1)):
    return X.point_data(base_metric(rho), diagonal([yscale * val for val in YSTD]),
                        LAMBDA, Fr(cH), Fr(cV), riem=True)


def relaxed(data):
    """Rows indexed by (a,b,d), a<b: unknown vector is (v^i,xi_i,n^i)."""
    dim = 14
    g, gi, ric, nr, ra = data["g"], data["gi"], data["Ric"], data["nRic"], data["RA"]
    rh = X.mmul(gi, ric)
    nh = {a: X.mmul(gi, nr[a]) for a in range(dim)}
    rows, rhs = [], []
    for a in range(dim):
        for b in range(a + 1, dim):
            for d in range(dim):
                row = [Fr(0)] * (3 * dim)
                for cc in range(dim):
                    row[cc] = ra[a][b][d][cc]
                    row[dim + cc] = -(int(a == cc) * rh[d][b] - int(b == cc) * rh[d][a])
                    row[2 * dim + cc] = -(g[a][cc] * int(d == b) - g[b][cc] * int(d == a))
                rows.append(row)
                rhs.append(nh[a][d][b] - nh[b][d][a])
    return rows, rhs


def qmatrix(rows):
    return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).convert_to(sp.QQ)


def certificate(data, label):
    rows, rhs = relaxed(data)
    aug = [row + [val] for row, val in zip(rows, rhs)]
    matrix = qmatrix(aug)
    _, independent_rows = matrix.transpose().rref()
    check(f"{label}: 43 independent rows of the augmented 1274-by-43 matrix", len(independent_rows) == 43)
    selected = [aug[index] for index in independent_rows]
    determinant = sp.QQ.to_sympy(qmatrix(selected).det())
    check(f"{label}: exact augmented minor is nonzero", determinant != 0)
    print(f"  zero-based row indices: {list(independent_rows)}", flush=True)
    print(f"  determinant: {determinant}", flush=True)
    return rows, rhs


def pullback_control(canonical, rho, cH, cV, s):
    """phi(x,y)=(s*x,q*y/s^2), phi*(alpha*h_original)=h_canonical."""
    rho, cH, cV, s = map(Fr, (rho, cH, cV, s))
    alpha = Fr(3) / cV
    q = Fr(2) * cV / (Fr(3) * cH)
    yscale = q / s**2
    target = point(rho, cH, cV, yscale)
    jac = [s] * 4 + [yscale] * 10
    m = 14
    check(f"rho={rho}: base dilation gives canonical signed amplitude", abs(rho) * s**2 == 1)
    check(f"rho={rho}: exact pulled-back metric agrees with canonical metric",
          all(canonical["g"][i][j] == alpha * jac[i] * jac[j] * target["g"][i][j]
              for i in range(m) for j in range(m)))
    check(f"rho={rho}: exact pulled-back covariant Ricci agrees",
          all(canonical["Ric"][i][j] == jac[i] * jac[j] * target["Ric"][i][j]
              for i in range(m) for j in range(m)))
    keys = set(canonical["R"]) | set(target["R"])
    check(f"rho={rho}: exact pulled-back mixed Riemann agrees",
          all(canonical["R"].get(key, Fr(0)) == jac[key[1]] * jac[key[2]] * jac[key[3]] / jac[key[0]]
              * target["R"].get(key, Fr(0)) for key in keys))
    check(f"rho={rho}: exact pulled-back covariant derivative of Ricci agrees",
          all(canonical["nRic"][k][i][j] == jac[k] * jac[i] * jac[j] * target["nRic"][k][i][j]
              for k in range(m) for i in range(m) for j in range(m)))
    rows, rhs = relaxed(target)
    rka = qmatrix(rows).rank()
    rkab = qmatrix([row + [val] for row, val in zip(rows, rhs)]).rank()
    check(f"rho={rho}, cH={cH}, cV={cV}: ranks (42,43) after nonunit scaling", (rka, rkab) == (42, 43))


def main():
    start = time.time()
    print("U. Normal quadratic rank-one-Ricci family: exact sign certificates and scaling controls", flush=True)
    print(f"referee engine sha256: {hashlib.sha256((HERE / 'xcheck_round10.py').read_bytes()).hexdigest()}", flush=True)
    print(f"sympy {sp.__version__}; lambda=1/5; canonical cH=2,cV=3; y=diag{YSTD}", flush=True)
    canonical = {}
    for sign in (1, -1):
        bg = X.base_geometry(base_metric(sign))
        check(f"rho={sign}: base Ricci is diag(rho,0,0,0)", bg["Ric"] == diagonal((sign, 0, 0, 0)))
        canonical[sign] = point(sign, 2, 3)
        certificate(canonical[sign], f"rho={sign}")
    pullback_control(canonical[1], Fr(4), Fr(5), Fr(7), Fr(1, 2))
    pullback_control(canonical[-1], Fr(-9), Fr(3, 2), Fr(5, 2), Fr(1, 3))
    print("The two exact minors prove the two signed representatives; the published naturality/homothety argument covers all nonzero amplitudes and positive cH,cV in the declared metric family.", flush=True)
    print("Arbitrary rank-one-Ricci metrics, Weyl tensors and free higher connection jets are not covered.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time() - start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
