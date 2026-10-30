"""GU round 2 (D): explicit self-adjoint Einstein Shiab operators in the Cl(7,7) model; Helmholtz test of (9.7).

Setting as in shiab_checks.py (d = 14, eta = (+^7,-^7), g = u(64,64) inside Cl(7,7) (x) C = M_128(C), pairing
p(T, X) = Tr(top(T ^ X)) for T in Omega^1(g), X in Omega^13(g)).  Formalization: D-selfadjoint-shiab.md.

Operators (Theorem D3): Phi3 = sum_{a<b<c} e^abc (x) i gamma_a gamma_b gamma_c,  Phi3' = sum e^abc (x) i gamma_abc Gamma
   (Gamma = gamma_0...gamma_13); both are Spin(7,7)-invariant elements of Lambda^3 (x) u(64,64).
   Phi11 = *Phi3,  Phi11' = *Phi3'  (Def. 8.1, r = 11),   W_s(xi) = [Phi11 ^ xi]_s,  W'_s(xi) = [Phi11' ^ xi]_s,
   brackets s in {A, C} as in the draft's appendix: [a^b]_C = sum e^I^e^J (a_I b_J - b_J a_I), [a^b]_A = i sum (... + ...).
   These are members of the draft's family (8.1) with r = 11, m = n = 0.
Tests
  E1  formal self-adjointness of W_A o d_B, W'_A o d_B (symbol criterion of Lemma C6a: beta_b(T,a) + beta_b(a,T) = 0);
      W_C is skew (beta_b(T,a) = beta_b(a,T)).  The 8 single-term members of the draft's family (8.1) that map
      Omega^2 -> Omega^13 (r = 11, m = n = 0; r = 1, m = 1; Phi_r or Phi_r Gamma; bracket A or C): exact null space
      of the self-adjointness defect.  Round-1 counterexample re-evaluated.
  E2  on Levi-Civita spin curvature F_R: W_A(F_R) is pure grade 1, W'_A(F_R) pure grade 13, both = x * (Ric - 1/2 Scal eta)
      exactly; comparison with the round-1 member R-S (C;A,A), c = 1 (equal on curvature up to a constant, different
      operators); Weyl curvature annihilated exactly by W, W' and all six R-S terms.
  E3  Helmholtz condition for Upsilon(T) = Shiab F_{B+T} + kappa *T:  H(a, b; T) = p(a, S[T^b]) - p(b, S[T^a]) must vanish;
      [T^b] := T^b + b^T.  Values for W_A, W'_A (2x2 determinant: no nonzero combination passes); controls:
      Chern-Simons (d = 3, S = Id) and the non-invariant d = 14 operator xi -> Omega ^ xi, Omega = e^{3..13} (x) 1.
  E4  first variation of (9.4) with S = W_A: cubic C(T) = p(T, S(T^T)); exact derivative check C(T+a) - C(T-a) =
      2 dC_T(a) + 2 C(a); residual  dC_T(a)/3 - p(a, S(T^T)) = -H(a, T; T)/3 != 0; closed form
      H(a, T; T) = c * p(a, [T ^ [T ^ Phi11]_C]_A) with c constant.
Exact arithmetic: all inputs have integer (Gaussian-integer) or dyadic entries and every operation is a signed
permutation, sum or product; all intermediate magnitudes stay < 2^53 (asserted), so float64/complex128 results
are exact and equalities are tested with ==.  Fixed RNG seed.
"""
import os
for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')          # 128x128 products: threading overhead dominates
import itertools
from fractions import Fraction
import numpy as np
import sympy as sp
from clifford77 import clifford_pq
from shiab_checks import Ext, Shiab, comm, acomm, pair, kn, ricci, random_curvature, spin_curvature

rng = np.random.default_rng(20210401)
n = 14
pos, neg = clifford_pq(7)
gam = [g.astype(complex) for g in pos + neg]
eta_d = [1] * 7 + [-1] * 7
eta = np.diag(eta_d).astype(float)
gam_up = [eta_d[c] * gam[c] for c in range(n)]
E = Ext(eta_d)
Id = np.eye(128, dtype=complex)
BR = {'C': comm, 'A': acomm}


def mul(a, b):
    return a @ b


def mono(idx):
    M = Id.copy()
    for i in idx:
        M = M @ gam[i]
    return M


Gam = mono(range(n))
Bf = mono(range(7, 14))
Bfi = np.round(np.linalg.inv(Bf).real).astype(complex)
assert np.array_equal(Bfi @ Bf, Id)
BIG = 2.0 ** 52


def in_u(X):
    return np.array_equal(X.conj().T @ Bf + Bf @ X, 0 * Id)


def chk(z):
    """assert exactness bound; return z"""
    zz = np.asarray(z)
    assert np.all(np.abs(zz) < BIG)
    return z


def rand_u(r=1):
    M = rng.integers(-r, r + 1, (128, 128)) + 1j * rng.integers(-r, r + 1, (128, 128))
    X = M - Bfi @ M.conj().T @ Bf
    assert in_u(X)
    return X


KMON = {1: [mono((a,)) for a in range(n)], 2: [mono(ab) for ab in itertools.combinations(range(n), 2)],
        13: [mono((a,)) @ Gam for a in range(n)], 14: [Gam]}


def rand_k(grades=(1, 2, 13, 14), r=2):
    X = sum(int(rng.integers(-r, r + 1)) * Mm for g in grades for Mm in KMON[g])
    return X


def form1(gen):
    return {(a,): gen() for a in range(n)}


def addf(f, g, s=1):
    out = dict(f)
    for K, v in g.items():
        out[K] = out.get(K, 0) + s * v
    return out


def fnorm(f):
    return max((np.abs(v).max() for v in f.values()), default=0.0)


def gcomm1(T, b):
    """[T ^ b] = T^b + b^T for two g-valued 1-forms (coefficients [T_i, b_j])."""
    return addf(E.wedge(T, b, mul), E.wedge(b, T, mul))


def beta(S, b, T, a):
    return chk(pair(E, T, S(E.wedge({(b,): Id}, a, mul))))


def as_vec(X):
    """13-form X = *( sum_a e^a (x) M_a ) -> list M_a"""
    out = []
    for a in range(n):
        Ic, c = E.star_coeff((a,))
        out.append(X.get(Ic, 0 * Id) / c)
    return out


def main():
    Phi3 = {abc: 1j * mono(abc) for abc in itertools.combinations(range(n), 3)}
    Phi3p = {abc: 1j * mono(abc) @ Gam for abc in itertools.combinations(range(n), 3)}
    print(f"setup: Gamma^2 = 1: {np.array_equal(Gam @ Gam, Id)}; i gamma_abc in u(64,64): "
          f"{all(in_u(v) for v in Phi3.values())}; i gamma_abc Gamma in u(64,64): {all(in_u(v) for v in Phi3p.values())}; "
          f"grades 1,2,13,14 in u(64,64): {all(in_u(Mm) for g in KMON for Mm in KMON[g])}")
    Phi11, Phi11p = E.star(Phi3), E.star(Phi3p)
    W = {s: (lambda xi, s=s: E.wedge(Phi11, xi, BR[s])) for s in 'AC'}
    Wp = {s: (lambda xi, s=s: E.wedge(Phi11p, xi, BR[s])) for s in 'AC'}
    print(f"degrees: Phi11 has {len(Phi11)} components of degree {set(len(K) for K in Phi11)}; "
          f"W(e^01 (x) X) has degrees {set(len(K) for K in W['A']({(0, 1): Id}))}")

    print("\n=== E1. Formal self-adjointness (symbol criterion, Lemma C6a); exact ===")
    T, al = form1(rand_u), form1(rand_u)
    for name, S, sgn in [("W_A  = [Phi11 ^ .]_A", W['A'], +1), ("W'_A = [Phi11' ^ .]_A", Wp['A'], +1),
                         ("W_C  = [Phi11 ^ .]_C", W['C'], -1), ("W'_C = [Phi11' ^ .]_C", Wp['C'], -1)]:
        d, mag = [], []
        for b in range(n):
            x, y = beta(S, b, T, al), beta(S, b, al, T)
            d.append(abs(x + sgn * y))
            mag.append(max(abs(x), abs(y)))
        kind = "beta_b(T,a) + beta_b(a,T)" if sgn > 0 else "beta_b(T,a) - beta_b(a,T)"
        print(f"  {name}: max_b |{kind}| = {max(d):.1f}  (max_b |beta_b| = {max(mag):.3e}; T, a random integer in "
              f"Omega^1 (x) u(64,64))  => {'self-adjoint' if sgn > 0 else 'skew-adjoint'}: {max(d) == 0}")
    # the draft's single-term family (8.1) landing in Omega^13: r = 11 (m = n = 0) or r = 1 (m = 1, n = 0),
    # Phi_r in the 2-dim invariant space (Phi_r, Phi_r Gamma), bracket +-: 8 operators; which combinations are s.a.?
    Phi1 = {(m,): gam[m] for m in range(n)}
    Phi1p = {(m,): gam[m] @ Gam for m in range(n)}
    fam = [("[Phi11^.]_A", W['A']), ("[Phi11'^.]_A", Wp['A']), ("[Phi11^.]_C", W['C']), ("[Phi11'^.]_C", Wp['C'])]
    for P, lab in ((Phi1, "Phi1"), (Phi1p, "Phi1'")):
        for s_ in 'AC':
            fam.append((f"[{lab}^*.]_{s_}", (lambda xi, P=P, s_=s_: E.wedge(P, E.star(xi), BR[s_]))))
    rows = []
    for trial in range(2):
        T2, a2 = form1(rand_u), form1(rand_u)
        for b in range(n):
            rows.append([beta(S_, b, T2, a2) + beta(S_, b, a2, T2) for _, S_ in fam])
    Dm_ = sp.Matrix([[sp.Integer(int(z.real)) + sp.I * sp.Integer(int(z.imag)) for z in r] for r in rows])
    assert all(z.real == int(z.real) and z.imag == int(z.imag) for r in rows for z in r)
    ns = Dm_.nullspace()
    print(f"  (8.1) single-term members into Omega^13 ({len(fam)}: " + ", ".join(nm for nm, _ in fam) + ")")
    print(f"    exact rank of the self-adjointness defect matrix ({Dm_.shape[0]} x {Dm_.shape[1]}) = {Dm_.rank()}; "
          f"null space (self-adjoint combinations) = " + "; ".join(str(list(v)) for v in ns))
    Tx, ax = {(0,): gam[0] @ gam[1]}, {(1,): gam[0]}
    print(f"  round-1 counterexample T = gamma_0 gamma_1 e^0, a = gamma_0 e^1: W_A gives beta_0(T,a) = "
          f"{beta(W['A'], 0, Tx, ax).real:+.0f}, beta_0(a,T) = {beta(W['A'], 0, ax, Tx).real:+.0f} "
          f"(R-S (C;A,A): +256, 0)")

    print("\n=== E2. Levi-Civita curvature: Einstein tensor, Weyl annihilation; exact ===")
    R = random_curvature(n, rng)
    Ric, s = ricci(R, eta)
    F = spin_curvature(R, gam_up, n)
    res = {}
    for name, S, post in [("W_A ", W['A'], Id), ("W'_A", Wp['A'], Gam)]:
        Mv = as_vec(S(F))
        K = np.array([[np.trace(Mv[a] @ post @ gam[c]) / 128 for c in range(n)] for a in range(n)])
        resid = max(np.abs(Mv[a] - sum(K[a, b] * gam_up[b] @ post for b in range(n))).max() for a in range(n))
        assert np.array_equal(K.imag, 0 * K.real)
        Kq = [[Fraction(K[a, b].real) for b in range(n)] for a in range(n)]
        Rq = [[Fraction(Ric[a, b]) for b in range(n)] for a in range(n)]
        sq = Fraction(s)
        a0, b0 = next((a, b) for a in range(n) for b in range(n) if a != b and Ric[a, b] != 0)
        x = Kq[a0][b0] / Rq[a0][b0]
        y = (Kq[0][0] - x * Rq[0][0]) / (sq * int(eta[0, 0]))
        exact = all(Kq[a][b] == x * Rq[a][b] + y * sq * int(eta[a, b]) for a in range(n) for b in range(n))
        res[name] = x
        lab = "gamma^b" if post is Id else "gamma^b Gamma"
        print(f"  {name}(F_R) = *( e^a (x) K_ab {lab} ): residual outside this grade = {resid:.1f}; "
              f"K = x Ric + y Scal eta exactly (rational arithmetic): {exact}; x = {x}, y = {y}, y/x = {y / x}")
    RS = Shiab(E, gam)
    rsF = RS.apply(F, 'C', 'A', 'A', 1.0)
    wF = W['A'](F)
    lam = wF[next(iter(wF))][np.unravel_index(np.argmax(np.abs(wF[next(iter(wF))])), (128, 128))] / \
        rsF[next(iter(wF))][np.unravel_index(np.argmax(np.abs(wF[next(iter(wF))])), (128, 128))]
    same = all(np.array_equal(wF.get(K, 0 * Id), lam * rsF.get(K, 0 * Id)) for K in set(wF) | set(rsF))
    xi = {ab: sum(int(rng.integers(-2, 3)) * Mm for Mm in KMON[2]) for ab in itertools.combinations(range(n), 2)}
    w2, r2 = W['A'](xi), RS.apply(xi, 'C', 'A', 'A', 1.0)
    diff = fnorm(addf(w2, {K: lam * v for K, v in r2.items()}, -1))
    print(f"  W_A(F_R) = {lam.real:+.0f} x R-S(C;A,A; c=1)(F_R) exactly: {same}; on a random grade-2-valued 2-form xi "
          f"(not a curvature tensor) max|W_A(xi) - ({lam.real:+.0f}) R-S(xi)| = {diff:.0f} (different operators)")
    W312 = 312 * R - kn(26 * Ric - s * eta, eta)          # 312 * Weyl(R), integer entries
    RicW, sW = ricci(W312, eta)
    FW = spin_curvature(W312, gam_up, n)
    outs = {"W_A": W['A'](FW), "W'_A": Wp['A'](FW), "W_C": W['C'](FW), "W'_C": Wp['C'](FW)}
    for s1 in 'CA':
        outs[f"R-S term1_{s1}"] = RS.term1(FW, s1)
    for s2, s3 in itertools.product('CA', repeat=2):
        outs[f"R-S term2_({s2},{s3})"] = RS.term2(FW, s2, s3)
    print(f"  Weyl part: 312 W = 312 R - KN(26 Ric - Scal eta, eta) has Ric = 0 exactly: {np.array_equal(RicW, 0 * RicW)}, "
          f"max|W| = {np.abs(W312).max():.0f}; images (max entry): "
          + ", ".join(f"{k}: {fnorm(v):.0f}" for k, v in outs.items()))

    print("\n=== E3. Helmholtz condition for Upsilon = Shiab F_{B+T} + kappa *T (T-linear part); exact ===")

    def H(S, a, b, T):
        return chk(pair(E, a, S(gcomm1(T, b))) - pair(E, b, S(gcomm1(T, a))))
    M = np.zeros((2, 2), dtype=complex)
    for i in range(2):
        a, b, T = form1(rand_k), form1(rand_k), form1(rand_k)
        M[i] = [H(W['A'], a, b, T), H(Wp['A'], a, b, T)]
    det = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    print(f"  a, b, T random integer in Omega^1 (x) so(8,8) [grades 1,2,13,14]: (H_W, H_W') for two samples = "
          f"{[[f'{z.real:+.0f}{z.imag:+.0f}i' for z in row] for row in M]}; det = {det.real:+.3e}{det.imag:+.1e}i")
    print(f"  => x W_A + y W'_A satisfies Helmholtz only for x = y = 0: {det != 0}")
    a, b, T = form1(rand_u), form1(rand_u), form1(rand_u)
    print(f"  a, b, T random integer in Omega^1 (x) u(64,64): H_W = {H(W['A'], a, b, T).real:+.4e}")
    # controls
    E3 = Ext([1, 1, 1])

    def ru2():
        Mm = rng.integers(-3, 4, (2, 2)) + 1j * rng.integers(-3, 4, (2, 2))
        return Mm - Mm.conj().T
    a3, b3, T3 = ({(i,): ru2() for i in range(3)} for _ in range(3))
    g3 = lambda T_, b_: addf(E3.wedge(T_, b_, mul), E3.wedge(b_, T_, mul))
    h3 = pair(E3, a3, g3(T3, b3)) - pair(E3, b3, g3(T3, a3))
    print(f"  control (Chern-Simons, d = 3, Shiab = Id, u(2)): H = {abs(h3):.1f}")
    Om = {tuple(range(3, 14)): Id}
    SO = lambda xi_: E.wedge(Om, xi_, mul)
    sa = max(abs(beta(SO, bb, a, T) + beta(SO, bb, T, a)) for bb in range(n))
    print(f"  control (d = 14, non-invariant Omega = e^(3..13) (x) 1, Shiab = Omega ^ .): self-adjointness defect = "
          f"{sa:.1f}, H = {abs(H(SO, a, b, T)):.1f}  (variational; excluded only by Spin(7,7)-invariance)")

    print("\n=== E4. First variation of (9.4) with Shiab = W_A: the cubic term; exact ===")
    S = W['A']
    C = lambda X: chk(pair(E, X, S(E.wedge(X, X, mul))))
    cs = []
    for trial in range(2):
        T, a = form1(rand_u), form1(rand_u)
        TT = E.wedge(T, T, mul)
        pTT = chk(pair(E, a, S(TT)))
        dC = chk(pTT + pair(E, T, S(gcomm1(T, a))))
        lhs = C(addf(T, a)) - C(addf(T, a, -1))
        ok_d = lhs == 2 * dC + 2 * C(a)
        CT = C(T)
        Hd = H(S, a, T, T)
        ok_r = (dC - 3 * pTT) == -Hd
        Y = E.wedge(T, Phi11, comm)
        hT = E.wedge(T, Y, acomm)
        q = chk(pair(E, a, hT))
        c = Fraction(Hd.real) / Fraction(q.real) if q.real != 0 else None
        cs.append(c)
        print(f"  sample {trial}: C(T) = {CT.real:+.4e}; derivative identity exact: {ok_d}; "
              f"3*(dC/3 - p(a,S(T^T))) = -H(a,T;T) exact: {ok_r}; H(a,T;T) = {Hd.real:+.4e}{Hd.imag:+.1e}i; "
              f"H / p(a,[T^[T^Phi11]_C]_A) = {c}")
    print(f"  => with k (T,S(T^T)) in the action, k = 1/3 is forced by a = T (C(T) != 0), and the EL equation of (9.4) is\n"
          f"     W_A F_A + kappa *T - (1/3) h(T) = 0 with p(a, h(T)) = H(a,T;T) = {cs[0]} * p(a, [T^[T^Phi11]_C]_A) != 0: "
          f"constant over samples: {cs[0] == cs[1]}")


if __name__ == '__main__':
    main()
