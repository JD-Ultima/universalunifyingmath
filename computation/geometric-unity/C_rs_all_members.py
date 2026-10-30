"""GU, Prop. C6 / referee issue GU-5 (reports/reviews/REVIEW-2026-10-30.md): no member of the round-1 Shiab family
R-S is formally self-adjoint, for any of the 8 bracket choices (s1; s2, s3) in {C, A}^3 and any c.

   Shiab^{(s1;s2,s3)}_c xi = [Phi1 ^ *xi]_{s1} - (c/2) *[Phi1 ^ *[Phi2 ^ *xi]_{s3}]_{s2}     (C-gauge-shiab.md, R-S)

Criterion (Lemma C6a): Shiab o d_B is formally self-adjoint iff beta_b(T, a) + beta_b(a, T) = 0 for all b and all
T, a in Lambda^1 (x) u(64,64), where beta_b(T, a) = Tr(T ^ Shiab(e^b ^ a))_top. This is P_b - (c/2) Q_b, with P from
term 1 and Q from term 2, so it is affine in c.

Part 1 (the referee's test, re-implemented with the repository seed): random u(64,64) samples; least-squares c over C
(which contains every real c); relative defect |P - (c/2) Q| / |P| at that c. Floating point (complex128).
Part 2 (exact certificates): T = X e^0, a = Y e^1, b = 0, with X, Y Clifford monomials (times 1 or i) in u(64,64).
Then beta_0(a, T) = 0 (e^0 ^ e^0 = 0), and beta_0(T, a) = p - (c/2) q with integer p, q. A certificate that no real c
works is one pair with q = 0 != p, or two pairs with p1 q2 - p2 q1 != 0. Integer-valued arithmetic in complex128
(gamma entries in {0, +-1}, all values far below 2^53); integrality is asserted.
Builds on shiab_checks.py and clifford77.py (imported unchanged).
"""
import itertools
import numpy as np
import shiab_checks as sc
from clifford77 import clifford_pq

n = 14
pos, neg = clifford_pq(7)
gam = [g.astype(complex) for g in pos + neg]
E = sc.Ext([1] * 7 + [-1] * 7)
Sh = sc.Shiab(E, gam)
B = np.eye(128, dtype=complex)
for a in range(7, 14):
    B = B @ gam[a]
Binv = np.linalg.inv(B)
Id = np.eye(128)
COMBOS = [(s1, s2, s3) for s1 in 'CA' for s2, s3 in itertools.product('CA', repeat=2)]


def in_u(X):
    return np.allclose(X.conj().T @ B + B @ X, 0)


def ext_b(b, f):
    return E.wedge({(b,): Id}, f, lambda u, v: u @ v)


def part1(rng, trials=2):
    print("=== Part 1. Random u(64,64) samples; least-squares c in C; relative defect at that c ===")

    def rand_u():
        M = rng.standard_normal((128, 128)) + 1j * rng.standard_normal((128, 128))
        return 0.5 * (M - Binv @ M.conj().T @ B)
    worst = {}
    for trial in range(trials):
        T = {(a,): rand_u() for a in range(n)}
        al = {(a,): rand_u() for a in range(n)}

        def defect(f):
            return np.array([sc.pair(E, T, f(ext_b(b, al))) + sc.pair(E, al, f(ext_b(b, T))) for b in range(n)])
        P = {s1: defect(lambda x, s1=s1: Sh.term1(x, s1)) for s1 in 'CA'}
        Q = {(s2, s3): defect(lambda x, s2=s2, s3=s3: Sh.term2(x, s2, s3)) for s2, s3 in itertools.product('CA', repeat=2)}
        for (s1, s2, s3) in COMBOS:
            p, q = P[s1], Q[(s2, s3)]
            cb = 2 * np.vdot(q, p) / np.vdot(q, q)
            res = np.linalg.norm(p - (cb / 2) * q) / np.linalg.norm(p)
            worst[(s1, s2, s3)] = min(worst.get((s1, s2, s3), np.inf), res)
            print(f"  trial {trial} ({s1};{s2},{s3}): |P| = {np.linalg.norm(p):.3e}, |Q| = {np.linalg.norm(q):.3e}, "
                  f"best c = {complex(cb):.4f}, relative defect at best c = {res:.3f}")
    print(f"  smallest relative defect at the best c, over all 8 choices and trials: "
          f"{min(worst.values()):.3f} (0 would mean some c works)")
    return worst


def monomial(K):
    M = np.eye(128, dtype=complex)
    for m in K:
        M = M @ gam[m]
    return M


def part2():
    print("\n=== Part 2. Exact certificates (T = X e^0, a = Y e^1, b = 0) ===")
    Ks = [(), (0,), (1,), (2,), (7,), (0, 1), (0, 2), (1, 2), (0, 7), (1, 7), (0, 1, 2), (0, 1, 7), (1, 2, 7)]
    cands = []
    for K in Ks:
        M = monomial(K)
        for ph, lab in ((1, ''), (1j, 'i ')):
            if in_u(ph * M):
                cands.append((lab + ('g' + ''.join(map(str, K)) if K else '1'), ph * M))
                break
    assert len(cands) == len(Ks)
    data = {}
    for (lx, X), (ly, Y) in itertools.product(cands, repeat=2):
        T = {(0,): X}
        al = {(1,): Y}
        xa, xT = ext_b(0, al), ext_b(0, T)
        rec = {}
        for s1 in 'CA':
            rec[s1] = sc.pair(E, T, Sh.term1(xa, s1))
            assert abs(sc.pair(E, al, Sh.term1(xT, s1))) == 0
        for s2, s3 in itertools.product('CA', repeat=2):
            rec[(s2, s3)] = sc.pair(E, T, Sh.term2(xa, s2, s3))
            assert abs(sc.pair(E, al, Sh.term2(xT, s2, s3))) == 0
        for k, v in rec.items():
            assert abs(v.imag) == 0 and v.real == round(v.real), (lx, ly, k, v)
        data[(lx, ly)] = {k: int(round(v.real)) for k, v in rec.items()}
    ok_all = True
    for (s1, s2, s3) in COMBOS:
        pq = [(key, d[s1], d[(s2, s3)]) for key, d in data.items()]
        cert = next(((k, p, q) for k, p, q in pq if q == 0 and p != 0), None)
        if cert:
            (lx, ly), p, q = cert
            print(f"  ({s1};{s2},{s3}): X = {lx}, Y = {ly}: beta_0(T,a) = {p} - (c/2)*{q} = {p} for every c; "
                  f"beta_0(a,T) = 0  => not self-adjoint for any c")
            continue
        nz = [(k, p, q) for k, p, q in pq if (p, q) != (0, 0)]
        cert2 = next(((u, v) for u, v in itertools.combinations(nz, 2) if u[1] * v[2] - v[1] * u[2] != 0), None)
        if cert2:
            (k1, p1, q1), (k2, p2, q2) = cert2
            print(f"  ({s1};{s2},{s3}): pairs {k1} and {k2}: beta_0 = {p1} - (c/2)*{q1} and {p2} - (c/2)*{q2}; "
                  f"det = {p1 * q2 - p2 * q1} != 0  => no common zero c")
        else:
            ok_all = False
            print(f"  ({s1};{s2},{s3}): no certificate in this candidate set")
    return ok_all


if __name__ == '__main__':
    rng = np.random.default_rng(20210401)
    print(f"sanity: gamma_m in u(64,64): {all(in_u(g) for g in gam)}")
    worst = part1(rng)
    ok = part2()
    print(f"\nall 8 bracket choices excluded exactly (Part 2): {ok}; numerically (Part 1, every defect > 0.1): "
          f"{all(r > 0.1 for r in worst.values())}")
