"""GU round 3 (G), Task 2: the fermionic operator (9.16): Shiab on spinor-valued forms, chirality, the
symbol criterion for a real (Majorana) or Hermitian action, and the principal symbol.

Formalization: formalization/geometric-unity/G-full-el-fermions.md, Sec. G.3.
Model: the explicit real Cl(7,7) on R^128 of clifford77.py (integer gammas, Gamma real diagonal, S^+- = +-1
eigenspaces of Gamma), d = 14, eta = (+^7, -^7).  Invariant bilinear forms (Theorem B1.3):
   B  = gamma_7...gamma_13  (symmetric, gamma_a^* = -gamma_a),   B' = gamma_0...gamma_6 (antisymmetric, gamma_a^* = +gamma_a);
   both pair S^+ with S^- (S^+ and S^- isotropic).  Hermitian form of U(64,64): H(psi, chi) = psi^dagger B chi.
Spinor Shiab operators Lambda^2 (x) S -> Lambda^13 (x) S (Spin(7,7)-equivariant; Prop. G6):
   Sh_Phi  xi = Phi11 ^ xi       Phi11 = *(sum e^abc (x) gamma_abc)  (natural Clifford contraction; = W_A form)
   Sh_Phi' xi = Phi11' ^ xi      Phi11' = *(sum e^abc (x) gamma_abc Gamma)
   Sh_gam  xi = *( sum_a e^a (x) sum_b gamma^b xi_ab )                (the second equivariant structure)
Tests (all exact: integer or Gaussian-integer arithmetic, == comparisons)
  F0  Sh_Phi o d is the Rarita-Schwinger operator: *^{-1} Sh_Phi(xi ^ zeta)_a = lam * gamma_a^{bc} xi_b zeta_c.
  F1  chirality: every listed Sh maps Lambda^2 (x) S^+- into Lambda^13 (x) S^-+ (block test); d_0 (even) preserves
      chirality, so the chirality labels of (9.16) are inconsistent (row zeta-bar_-: S^+ from the Shiab block,
      S^- from the nu block); linear independence of Sh_Phi, Sh_gam on S^+.
  F2  symbol criterion (Lemma G7): for Grassmann fields a real bilinear action needs beta_b(z1,z2) = C(z1 ^ Sh(e^b ^ z2))
      SYMMETRIC in (z1, z2); for commuting fields ANTISYMMETRIC (Theorem D2).  Tested for Sh_Phi, Sh_Phi', Sh_gam with
      C = B, B'; exact null space of the Grassmann defect on span(Sh_Phi, Sh_gam).
  F3  symmetry type of C gamma_I restricted to S^eps x S^eps (all monomials of grades 1, 3, 5, 7): fixes
      Lambda^2(S^eps)* = Lambda^1 + Lambda^5, Sym^2(S^eps)* = Lambda^3 + Lambda^7_eps.
  F4  cross block (zeta-bar, nu): with the printed signs (+d_0, -d_0^*) the Grassmann condition holds iff C is symmetric.
  F5  Hermitian (complex) reading: reality condition kappa_b(z1,z2) = -conj kappa_b(z2,z1) for Sh with gamma_abc vs
      i gamma_abc, and for the cross block.
  F6  principal symbol of D = [[*Sh_Phi d, d], [-d^*, 0]] on (zeta, nu) in (Lambda^1 (+) Lambda^0) (x) S: exact ranks at
      non-null and null covectors in (7,7), and in (14,0); per Weyl sector; and with the Shiab block removed
      (what survives in a Majorana action).
"""
import itertools
import os
import sys
import time

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from clifford77 import clifford_pq          # noqa: E402
from shiab_checks import Ext, perm_sign     # noqa: E402
from G_symbols import sparse_rank           # noqa: E402

n = 14
rng = np.random.default_rng(20210401)
pos, neg = clifford_pq(7)
GAM = [g.astype(np.int64) for g in pos + neg]
ETA = [1] * 7 + [-1] * 7
I128 = np.eye(128, dtype=np.int64)


def mono(gam, idx):
    M = np.eye(gam[0].shape[0], dtype=gam[0].dtype)
    for i in idx:
        M = M @ gam[i]
    return M


Gam = mono(GAM, range(n))
B = mono(GAM, range(7, 14))
Bp = mono(GAM, range(7))
SP = np.where(np.diag(Gam) == 1)[0]
SM = np.where(np.diag(Gam) == -1)[0]
mv = lambda M, v: M @ v


def make_shiabs(gam, eta, E):
    Phi3 = {abc: mono(gam, abc) for abc in itertools.combinations(range(n), 3)}
    G_ = mono(gam, range(n))
    Phi3p = {abc: mono(gam, abc) @ G_ for abc in itertools.combinations(range(n), 3)}
    Phi11, Phi11p = E.star(Phi3), E.star(Phi3p)
    gup = [eta[b] * gam[b] for b in range(n)]

    def sh_gam(xi):
        one = {}
        for (a, b), v in xi.items():           # xi_ab e^a^e^b, a<b : contributes to a-slot with gamma^b, b-slot with -gamma^a
            one[(a,)] = one.get((a,), 0) + mv(gup[b], v)
            one[(b,)] = one.get((b,), 0) - mv(gup[a], v)
        return E.star(one)
    return {'Phi': lambda xi: E.wedge(Phi11, xi, mv), "Phi'": lambda xi: E.wedge(Phi11p, xi, mv), 'gam': sh_gam}


def rand_form1(gen):
    return {(a,): gen() for a in range(n)}


def addf(f, g, s=1):
    out = dict(f)
    for K, v in g.items():
        out[K] = out.get(K, 0) + s * v
    return out


def pair_top(E, z1, X, C, herm=False):
    """C-pairing of the top coefficient of z1 ^ X (z1 a spinor 1-form, X a spinor 13-form)."""
    f = (lambda u, w: np.conj(u) @ C @ w) if herm else (lambda u, w: u @ C @ w)
    t = E.top(E.wedge(z1, X, f))
    return t


def main():
    t0 = time.time()
    E = Ext(ETA)
    SH = make_shiabs(GAM, ETA, E)

    print("=== F0. Sh_Phi o d = Rarita-Schwinger operator ===")
    gup = [ETA[b] * GAM[b] for b in range(n)]
    xi1 = rng.integers(-2, 3, n)
    z = rand_form1(lambda: rng.integers(-2, 3, 128))
    xz = E.wedge({(b,): int(xi1[b]) for b in range(n)}, z, lambda s, v: s * v)
    out = E.star(SH['Phi'](xz))                      # 1-form (** = +1 on 13-forms in (7,7))
    rs = {}
    for a in range(n):
        acc = np.zeros(128, dtype=np.int64)
        for b in range(n):
            for c in range(n):
                if len({a, b, c}) == 3:
                    # gamma_a^{bc} = gamma_a gamma^b gamma^c for distinct indices
                    acc += xi1[b] * (GAM[a] @ gup[b] @ gup[c] @ z[(c,)])
        rs[a] = acc
    lam = None
    ok = True
    for a in range(n):
        o = out.get((a,), np.zeros(128, dtype=np.int64))
        nz = np.flatnonzero(rs[a])
        if len(nz):
            l = o[nz[0]] / rs[a][nz[0]]
            lam = l if lam is None else lam
            ok &= np.array_equal(o, lam * rs[a])
        else:
            ok &= not np.any(o)
    print(f"  *^-1 Sh_Phi(xi ^ zeta)_a = lam * gamma_a gamma^b gamma^c xi_b zeta_c (b, c, a distinct) exactly, "
          f"lam = {lam:+.0f}: {ok}   (random integer xi, zeta)")

    print("\n=== F1. Chirality ===")
    for name in ('Phi', "Phi'", 'gam'):
        S = SH[name]
        flips, keeps = True, True
        for _ in range(3):
            b, a = rng.choice(n, 2, replace=False)
            for src, dst_same, dst_opp in ((SP, SP, SM), (SM, SM, SP)):
                v = np.zeros(128, dtype=np.int64)
                v[src] = rng.integers(-2, 3, len(src))
                xi = {tuple(sorted((int(a), int(b)))): v}
                res = S(xi)
                flips &= all(not np.any(w[dst_same]) for w in res.values())
                keeps &= all(not np.any(w[dst_opp]) for w in res.values())
        print(f"  Sh_{name:4s}: maps Lambda^2 (x) S^+- into Lambda^13 (x) S^-+ (flips chirality): {flips}; "
              f"preserves chirality: {keeps}")
    print("  d_0 = d + spin connection (Clifford grade 2, even) preserves chirality.  In row zeta-bar_- of (9.16) the entry\n"
          "  *Sh(d_0 + w_{+-}) acts on zeta_- and the entry d_0 + w_{+-} on nu_-: outputs S^+ and S^- respectively, which no\n"
          "  single invariant pairing with zeta-bar_- accepts (B, B' and the dual pairing are chirality-homogeneous).")
    # independence of Sh_Phi and Sh_gam on S^+ (dim Hom = 2, Prop. G6)
    rows = []
    for _ in range(2):
        v = np.zeros(128, dtype=np.int64)
        v[SP] = rng.integers(-2, 3, 64)
        xi = {(0, 1): v, (2, 9): v}
        r1, r2 = SH['Phi'](xi), SH['gam'](xi)
        for K in set(r1) | set(r2):
            rows.append(np.stack([np.asarray(r1.get(K, np.zeros(128, dtype=np.int64))),
                                  np.asarray(r2.get(K, np.zeros(128, dtype=np.int64)))], axis=1))
    Mx = np.vstack(rows)
    import sympy as sp
    print(f"  rank of [Sh_Phi(xi) | Sh_gam(xi)] on sample S^+-valued 2-forms (exact): "
          f"{sp.Matrix(Mx[np.any(Mx != 0, axis=1)][:400].tolist()).rank()}  (2 => independent)")

    print("\n=== F2. Symbol criterion: beta_b(z1,z2) = C(top z1 ^ Sh(e^b ^ z2)) ===")
    z1 = rand_form1(lambda: rng.integers(-2, 3, 128))
    z2 = rand_form1(lambda: rng.integers(-2, 3, 128))
    for Cn, C in (("B (sym)", B), ("B' (antisym)", Bp)):
        for name in ('Phi', "Phi'", 'gam'):
            S = SH[name]
            sym_def, anti_def, mag = 0, 0, 0
            for b in range(n):
                x12 = pair_top(E, z1, S(E.wedge({(b,): 1}, z2, lambda s, v: s * v)), C)
                x21 = pair_top(E, z2, S(E.wedge({(b,): 1}, z1, lambda s, v: s * v)), C)
                sym_def = max(sym_def, abs(x12 - x21))
                anti_def = max(anti_def, abs(x12 + x21))
                mag = max(mag, abs(x12), abs(x21))
            print(f"  C = {Cn:12s} Sh_{name:4s}: max_b|beta_b(z1,z2) - beta_b(z2,z1)| = {sym_def:6d} "
                  f"(Grassmann/Majorana-consistent iff 0); max_b|beta + beta^t| = {anti_def:6d} (bosonic s.a. iff 0); "
                  f"max|beta| = {mag}")
    # null space of the Grassmann defect on span(Sh_Phi, Sh_gam), C = B and B'
    for Cn, C in (("B", B), ("B'", Bp)):
        rows = []
        for _ in range(2):
            y1 = rand_form1(lambda: rng.integers(-2, 3, 128))
            y2 = rand_form1(lambda: rng.integers(-2, 3, 128))
            for b in range(n):
                rows.append([int(pair_top(E, y1, SH[s](E.wedge({(b,): 1}, y2, lambda s_, v: s_ * v)), C)
                                 - pair_top(E, y2, SH[s](E.wedge({(b,): 1}, y1, lambda s_, v: s_ * v)), C))
                             for s in ('Phi', 'gam')])
        M = sp.Matrix(rows)
        print(f"  C = {Cn:2s}: Grassmann-defect matrix on (s Sh_Phi + t Sh_gam): exact rank {M.rank()} of 2; "
              f"null space {M.nullspace()}")

    print("\n=== F3. Symmetry of C gamma_I on S^eps x S^eps (eps = +, -), all monomials of grades 1,3,5,7 ===")
    for Cn, C in (("B", B), ("B'", Bp)):
        line = []
        for k in (1, 3, 5, 7):
            types = set()
            nz = 0
            for I in itertools.combinations(range(n), k):
                M = C @ mono(GAM, I)
                for Sx in (SP, SM):
                    blk = M[np.ix_(Sx, Sx)]
                    if not np.any(blk):
                        types.add('zero')
                        continue
                    nz += 1
                    if np.array_equal(blk, blk.T):
                        types.add('sym')
                    elif np.array_equal(blk, -blk.T):
                        types.add('antisym')
                    else:
                        types.add('neither')
            line.append(f"k={k}: {sorted(types)}")
        print(f"  C = {Cn:2s}: " + "; ".join(line))
    print("  even grades vanish on S^eps x S^eps (C pairs S^eps with S^-eps).  Hence (Schur + dimensions 14+2002 = 2016,\n"
          "  364+1716 = 2080): Lambda^2(S^eps)* = Lambda^1 + Lambda^5 and Sym^2(S^eps)* = Lambda^3 + Lambda^7_eps; in particular\n"
          "  Hom_Spin(7,7)(Lambda^3, Lambda^2 S^eps*) = 0 and Hom(Lambda^3, Lambda^2 S*) = 0 (S^+ (x) S^- = even forms).")

    print("\n=== F4. Cross block (zeta-bar nu, nu-bar zeta) with printed signs (+d_0, -d_0^*) ===")
    for Cn, C in (("B (sym)", B), ("B' (antisym)", Bp)):
        y1, y2 = rand_form1(lambda: rng.integers(-2, 3, 128)), rand_form1(lambda: rng.integers(-2, 3, 128))
        n1, n2 = rng.integers(-2, 3, 128), rng.integers(-2, 3, 128)
        d = 0
        for b in range(n):
            k12 = y1[(b,)] @ C @ n2 + n1 @ C @ y2[(b,)]
            k21 = y2[(b,)] @ C @ n1 + n2 @ C @ y1[(b,)]
            d = max(d, abs(int(k12 - k21)))
        print(f"  C = {Cn:12s}: max_b |k_b(psi1,psi2) - k_b(psi2,psi1)| = {d}  (Grassmann-consistent iff 0)")

    print("\n=== F5. Hermitian (complex) reading, H(psi,chi) = psi^dagger B chi ===")
    cz = lambda: rng.integers(-2, 3, 128) + 1j * rng.integers(-2, 3, 128)
    w1, w2 = rand_form1(cz), rand_form1(cz)
    for lab, ph in (("gamma_abc (real)", 1), ("i gamma_abc (u(64,64))", 1j)):
        S = SH['Phi']
        dh, da = 0, 0
        for b in range(n):
            x12 = ph * pair_top(E, w1, S(E.wedge({(b,): 1}, w2, lambda s, v: s * v)), B, herm=True)
            x21 = ph * pair_top(E, w2, S(E.wedge({(b,): 1}, w1, lambda s, v: s * v)), B, herm=True)
            dh = max(dh, abs(x12 + np.conj(x21)))
            da = max(da, abs(x12 - np.conj(x21)))
        print(f"  Shiab block with {lab:24s}: max|kappa_b(z1,z2) + conj kappa_b(z2,z1)| = {dh:7.0f} "
              f"(real action <=> 0); max|kappa - conj kappa^t| = {da:7.0f}")
    m1, m2 = cz(), cz()
    dh, da = 0, 0
    for b in range(n):
        k12 = np.conj(w1[(b,)]) @ B @ m2 + np.conj(m1) @ B @ w2[(b,)]
        k21 = np.conj(w2[(b,)]) @ B @ m1 + np.conj(m2) @ B @ w1[(b,)]
        dh = max(dh, abs(k12 + np.conj(k21)))
        da = max(da, abs(k12 - np.conj(k21)))
    print(f"  cross block as printed:                    max|k + conj k^t| = {dh:7.0f}; max|k - conj k^t| = {da:7.0f}")

    print("\n=== F6. Principal symbol of D = [[*Sh_Phi d, d], [-d^*, 0]] ===")

    def symbol_matrix(gam, eta, xi, shiab_block=True):
        """rows: (m, i) [13-form missing m] (x) s_i, then ('top', i); cols: (a, j) zeta = e^a (x) s_j, then nu_j."""
        Eg = Ext(eta)
        N128 = gam[0].shape[0]
        dt = gam[0].dtype
        M = np.zeros(((n + 1) * N128, (n + 1) * N128), dtype=dt)
        IdM = np.eye(N128, dtype=dt)
        mm = lambda X, Y: X @ Y
        for a in range(n):
            # zeta = e^a (x) (all basis spinors at once): xi ^ e^a (x) Id
            if shiab_block:
                xa = Eg.wedge({(b,): xi[b] for b in range(n) if xi[b]}, {(a,): IdM}, lambda s, v: s * v)
                res = Eg.wedge(make_phi11(gam, Eg), xa, mm)
                for I, blk in res.items():
                    m = next(i for i in range(n) if i not in I)
                    M[m * N128:(m + 1) * N128, a * N128:(a + 1) * N128] += blk
            # nu-bar row: xi ^ *zeta  (top form)
            top = Eg.top(Eg.wedge({(b,): xi[b] for b in range(n) if xi[b]}, Eg.star({(a,): IdM}), lambda s, v: s * v))
            if not isinstance(top, int):
                M[n * N128:, a * N128:(a + 1) * N128] += top
        # nu column: *(xi (x) nu)
        st = Eg.star({(b,): xi[b] * IdM for b in range(n) if xi[b]})
        for I, blk in st.items():
            m = next(i for i in range(n) if i not in I)
            M[m * N128:(m + 1) * N128, n * N128:] += blk
        return M

    _phi_cache = {}

    def make_phi11(gam, Eg):
        key = id(gam)
        if key not in _phi_cache:
            _phi_cache[key] = Eg.star({abc: mono(gam, abc) for abc in itertools.combinations(range(n), 3)})
        return _phi_cache[key]

    def rank_of(M):
        if np.iscomplexobj(M):
            Re, Im = np.real(M).astype(np.int64), np.imag(M).astype(np.int64)
            assert np.array_equal(Re + 1j * Im, M)
            R = np.block([[Re, -Im], [Im, Re]])
            r, *_ = rank_of(R)
            assert r % 2 == 0
            return (r // 2,)
        rr, cc = np.nonzero(M)
        return sparse_rank(rr, cc, M[rr, cc], M.shape[0], M.shape[1])

    full = (n + 1) * 128
    cases = [("e^0 (eta=+1)", {0: 1}), ("e^7 (eta=-1)", {7: 1}), ("e^0+e^1 (eta=+2)", {0: 1, 1: 1}),
             ("e^0+e^7 (null)", {0: 1, 7: 1}), ("e^0+e^1+e^7+e^8 (null)", {0: 1, 1: 1, 7: 1, 8: 1})]
    rows_S = {'+': SP, '-': SM}
    for lab, xd in cases:
        xi = [xd.get(i, 0) for i in range(n)]
        M = symbol_matrix(GAM, ETA, xi)
        r = rank_of(M)[0]
        M0 = symbol_matrix(GAM, ETA, xi, shiab_block=False)
        r0 = rank_of(M0)[0]
        # Weyl sector: zeta in S^-, nu in S^+  (rows: 13-forms (x) S^+, top (x) S^-)
        cz_ = np.concatenate([a * 128 + SM for a in range(n)] + [n * 128 + SP])
        rz_ = np.concatenate([m * 128 + SP for m in range(n)] + [n * 128 + SM])
        Mw = M[np.ix_(rz_, cz_)]
        off = np.abs(M[np.ix_(np.setdiff1d(np.arange(full), rz_), cz_)]).max()
        rw = rank_of(Mw)[0]
        print(f"  (7,7) xi = {lab:24s}: rank sigma(D) = {r:4d} / {full} (deficit {full - r:4d}); Weyl sector "
              f"(zeta in S^-, nu in S^+): {rw:3d} / {full // 2} (block-diagonal: {off == 0}); "
              f"without the Shiab block (Majorana): {r0:4d}")
    gamE = [GAM[a].astype(complex) * (1 if a < 7 else 1j) for a in range(n)]
    for lab, xd in [("e^0", {0: 1}), ("e^0+e^1", {0: 1, 1: 1}), ("e^0+e^7", {0: 1, 7: 1})]:
        xi = [xd.get(i, 0) for i in range(n)]
        assert all(np.array_equal(gamE[a] @ gamE[b] + gamE[b] @ gamE[a], 2 * (a == b) * np.eye(128)) for a in range(n)
                   for b in range(n))
        _phi_cache.clear()
        M = symbol_matrix(gamE, [1] * 14, xi)
        r = rank_of(M)[0]
        print(f"  (14,0) xi = {lab:8s}: rank sigma(D) = {r:4d} / {full} (deficit {full - r})   "
              f"[gamma'_a = gamma_a (a<7), i gamma_a (a>=7)]")
    print(f"\nruntime {time.time() - t0:.0f} s")


if __name__ == '__main__':
    main()
