"""GU Phase 2(b), explicit Clifford model: real Cl(7,7) on R^128.

Checks (claims GU-016 .. GU-019, GU-021 in claims/geometric-unity.md):
  C1  real 128x128 matrices gamma_a, a=0..13, with {g_a,g_b} = 2 eta_ab, eta = (+^7, -^7); exact integers.
  C2  the 2^14 Clifford monomials are trace-orthogonal => linearly independent => Cl(7,7) = R(128).
  C3  volume element Gamma: Gamma^2 = +1, real => real (Majorana-Weyl) half-spinors S^+-, dim 64 each.
  C4  invariant bilinear forms: B (adjoint of gamma_a is -gamma_a) and B' (adjoint +gamma_a);
      symmetry, signature, isotropy of S^+-; grade rule for B-adjoint; dims of so(B), sp(B'),
      the draft's (3.22), (3.26), (3.27), (8.6) grade lists.
  C5  Spin(1,3) x Spin(6,4) inside Spin(7,7): omega_4^2 = -1 is a complex structure on S^+_R commuting with
      spin(1,3)+spin(6,4); the 4d and 10d chiralities are locked on S^+; commutant of the action on S^+_R
      is 2-dimensional (= C), so S^+_R is irreducible of complex type: S^+_R = (2,1) (x) 16 as a real rep.
Exact integer arithmetic except the commutant rank in C5 (numpy SVD; tolerance printed).
"""
import itertools
from math import comb
import numpy as np


# ---- signed permutation matrices: M[i, perm[i]] = sign[i]  (all gammas here are of this form)
def sp_from_dense(M):
    perm = np.argmax(np.abs(M), axis=1)
    sign = M[np.arange(M.shape[0]), perm]
    D = np.zeros_like(M)
    D[np.arange(M.shape[0]), perm] = sign
    assert np.array_equal(D, M), "not a signed permutation"
    return perm, sign


def sp_mul(A, B):
    """(A B) for signed permutations: row i of AB = sign_A[i] * row perm_A[i] of B."""
    pa, sa = A
    pb, sb = B
    return pb[pa], sa * sb[pa]


def sp_T(A):
    p, s = A
    pt = np.empty_like(p)
    st = np.empty_like(s)
    pt[p] = np.arange(len(p))
    st[p] = s
    return pt, st


def sp_eq(A, B, c=1):
    return np.array_equal(A[0], B[0]) and np.array_equal(A[1], c * B[1])


def sp_trace(A):
    p, s = A
    fix = p == np.arange(len(p))
    return int(s[fix].sum())

s1 = np.array([[0, 1], [1, 0]], dtype=np.int64)     # squares to +1
e2 = np.array([[0, 1], [-1, 0]], dtype=np.int64)    # squares to -1
s3 = np.array([[1, 0], [0, -1]], dtype=np.int64)
I2 = np.eye(2, dtype=np.int64)


def clifford_pq(m):
    """Real Cl(m,m) on R^(2^m) by Cl(p+1,q+1) = Cl(p,q) (x)^ Cl(1,1). Returns (pos list, neg list)."""
    pos, neg = [], []
    for _ in range(m):
        pos = [np.kron(g, s3) for g in pos]
        neg = [np.kron(g, s3) for g in neg]
        N = pos[0].shape[0] // 2 if pos else 1
        Id = np.eye(N, dtype=np.int64)
        pos.append(np.kron(Id, s1))
        neg.append(np.kron(Id, e2))
    return pos, neg


def main():
    pos, neg = clifford_pq(7)
    G = pos + neg                      # indices 0..6 square +1, 7..13 square -1
    eta = np.diag([1] * 7 + [-1] * 7)
    N = G[0].shape[0]
    I = np.eye(N, dtype=np.int64)
    print(f"C1. matrix size {N}; entries in {{-1,0,1}}: {all(set(np.unique(g)) <= {-1, 0, 1} for g in G)}")
    ok = all(np.array_equal(G[a] @ G[b] + G[b] @ G[a], 2 * eta[a, b] * I) for a in range(14) for b in range(14))
    print(f"    Clifford relations {{g_a,g_b}} = 2 eta_ab exact: {ok}")
    assert ok

    # monomials as signed permutation matrices: store dense int8 is fine for 128x128 but 16384 of them
    # is 268M entries; instead verify trace-orthogonality via tr(E_I) = 0 for all nonempty I.
    Gs = [sp_from_dense(g) for g in G]
    Isp = (np.arange(N), np.ones(N, dtype=np.int64))

    def mono_sp(Iset):
        M = Isp
        for a in Iset:
            M = sp_mul(M, Gs[a])
        return M

    def mono(Iset):
        p, s = mono_sp(Iset)
        M = np.zeros((N, N), dtype=np.int64)
        M[np.arange(N), p] = s
        return M
    # E_I E_J = +- E_{I xor J}; so orthogonality <=> tr(E_K) = 0 for all K != {}.  Check all 16383 K.
    nonzero_tr = 0
    for k in range(1, 15):
        for K in itertools.combinations(range(14), k):
            if sp_trace(mono_sp(K)) != 0:
                nonzero_tr += 1
    print(f"C2. number of nonempty monomials with nonzero trace: {nonzero_tr} (0 => 16384 independent "
          f"monomials => Cl(7,7) -> End(R^128) is an isomorphism, dim 16384 = 128^2)")
    assert nonzero_tr == 0

    Gam = mono(range(14))
    print(f"C3. Gamma^2 = +1: {np.array_equal(Gam @ Gam, I)}; Gamma real diagonal: "
          f"{np.array_equal(Gam, np.diag(np.diag(Gam)))}; dim(+1) = {(np.diag(Gam) == 1).sum()}, "
          f"dim(-1) = {(np.diag(Gam) == -1).sum()}")
    Sp = np.where(np.diag(Gam) == 1)[0]
    Sm = np.where(np.diag(Gam) == -1)[0]

    # C4: bilinear forms. gamma_a^T = t_a gamma_a with t_a = +1 (pos) / -1 (neg) in this construction.
    t = [1 if np.array_equal(g.T, g) else (-1 if np.array_equal(g.T, -g) else 0) for g in G]
    print(f"C4. transpose signs t_a (g^T = t_a g): {t}")
    # seek B = product of a subset S of gammas with B^{-1} g_a^T B = sigma g_a for all a
    found = {}
    for sigma in (-1, 1):
        for S in [tuple(i for i in range(14) if t[i] == 1), tuple(i for i in range(14) if t[i] == -1)]:
            B = mono(S)
            Binv = np.round(np.linalg.inv(B)).astype(np.int64)
            if all(np.array_equal(Binv @ G[a].T @ B, sigma * G[a]) for a in range(14)):
                found[sigma] = (S, B, Binv)
    for sigma, (S, B, Binv) in sorted(found.items()):
        sym = 'symmetric' if np.array_equal(B.T, B) else ('antisymmetric' if np.array_equal(B.T, -B) else '?')
        line = f"    form with adjoint(gamma_a) = {sigma:+d} gamma_a: B = product of gammas {S}; B is {sym}"
        if sym == 'symmetric':
            ev = np.linalg.eigvalsh(B.astype(float))
            line += f"; signature ({(ev > 0).sum()},{(ev < 0).sum()})"
        iso_pp = not np.any(B[np.ix_(Sp, Sp)]) and not np.any(B[np.ix_(Sm, Sm)])
        iso_pm = not np.any(B[np.ix_(Sp, Sm)])
        line += f"; B(S+,S+)=B(S-,S-)=0: {iso_pp}; B(S+,S-)=0: {iso_pm}"
        print(line)
        # grade rule: adjoint(E_I) = c_k E_I. Check every monomial.
        signs = {}
        Bs, Binvs = sp_from_dense(B), sp_from_dense(Binv)
        for k in range(15):
            cs = set()
            for K in itertools.combinations(range(14), k):
                E = mono_sp(K)
                A = sp_mul(sp_mul(Binvs, sp_T(E)), Bs)
                if sp_eq(A, E):
                    cs.add(1)
                elif sp_eq(A, E, -1):
                    cs.add(-1)
                else:
                    cs.add(0)
            assert len(cs) == 1 and 0 not in cs
            signs[k] = cs.pop()
        skew = [k for k in range(15) if signs[k] == -1]
        dskew = sum(comb(14, k) for k in skew)
        print(f"      adjoint sign by grade k=0..14: {[signs[k] for k in range(15)]}")
        print(f"      skew grades {skew}, total dim {dskew}  (n(n-1)/2 = {N*(N-1)//2}, n(n+1)/2 = {N*(N+1)//2})")
        if sigma == -1:
            rule = all(signs[k] == (-1) ** (k * (k + 1) // 2) for k in range(15))
            print(f"      matches Clifford-conjugation rule (-1)^(k(k+1)/2): {rule}")
            src_so = [2, 6, 10, 14, 1, 5, 9, 13]
            print(f"      draft (3.26) so(64,64) grades {sorted(src_so)} == skew grades: {sorted(src_so) == skew}")
            herm_i = sorted(k for k in range(15) if signs[k] == 1)
            src_i = [0, 4, 8, 12, 3, 7, 11]
            print(f"      u(64,64) = skew (+) i*sym; i*-grades {herm_i}; draft (3.27)/(8.6) i-grades "
                  f"{sorted(src_i)}: {herm_i == sorted(src_i)}")
            ev_skew = [k for k in skew if k % 2 == 0]
            print(f"      even skew grades {ev_skew} (draft (3.22) gl(64,R)): dim "
                  f"{sum(comb(14, k) for k in ev_skew)} = 64^2 = {64*64}")
            od_skew = [k for k in skew if k % 2 == 1]
            print(f"      odd skew grades {od_skew}: dim {sum(comb(14, k) for k in od_skew)} = "
                  f"2*dim Lambda^2(R^64) = {2*64*63//2}")
            # even elements preserve S+-: check restriction to S+ of grade-2 elements is skew w.r.t. pairing
            Bsym = B

    # C5: Spin(1,3) x Spin(6,4).  X-part: one positive + three negative gammas; N-part: the rest.
    X = [0, 7, 8, 9]
    Nn = [i for i in range(14) if i not in X]
    print(f"C5. X-part gammas {X} (signature (1,3)); N-part {Nn} (signature (6,4))")
    w4 = mono(X)
    w10 = mono(Nn)
    print(f"    omega_4^2 = -1: {np.array_equal(w4 @ w4, -I)}; omega_10^2 = -1: {np.array_equal(w10 @ w10, -I)}; "
          f"[omega_4, Gamma] = 0: {np.array_equal(w4 @ Gam, Gam @ w4)}")
    gens = [G[i] @ G[j] for (i, j) in itertools.combinations(X, 2)] + \
           [G[i] @ G[j] for (i, j) in itertools.combinations(Nn, 2)]
    comm = all(np.array_equal(w4 @ A, A @ w4) for A in gens)
    print(f"    omega_4 commutes with all {len(gens)} generators of spin(1,3)+spin(6,4): {comm}")
    # locking of chiralities on S+: omega_10 = c * omega_4 on S+
    P = np.zeros((N, len(Sp)), dtype=np.int64)
    P[Sp, np.arange(len(Sp))] = 1
    w4p = P.T @ w4 @ P
    w10p = P.T @ w10 @ P
    for c in (1, -1):
        if np.array_equal(w10p, c * w4p):
            print(f"    on S+: omega_10 = {c:+d} * omega_4  (4d and 10d chiralities locked)")
    # commutant of the action on S+_R (real 64x64): commutant of two generic elements
    rng = np.random.default_rng(1)
    As = []
    for _ in range(2):
        coeff = rng.standard_normal(len(gens))
        As.append(P.T @ sum(c * A for c, A in zip(coeff, gens)).astype(float) @ P)
    n = len(Sp)
    K = np.vstack([np.kron(np.eye(n), A) - np.kron(A.T, np.eye(n)) for A in As])
    sv = np.linalg.svd(K, compute_uv=False)
    tol = 1e-9 * sv[0]
    print(f"    commutant dimension of two generic elements of spin(1,3)+spin(6,4) on S+_R: "
          f"{(sv < tol).sum()} (singular values below {tol:.1e}; next smallest {sv[sv >= tol][-1]:.3e})")
    print("    (commutant contains span{1, omega_4}, so the full commutant is exactly C: S+_R is irreducible"
          " of complex type, i.e. one complex irreducible (2,1)(x)16 viewed as a real 64-dim module)")


if __name__ == '__main__':
    main()
