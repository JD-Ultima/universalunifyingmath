"""GU round 2 (D): invariant theory of the space of Shiab operators.

Setting (formalization/geometric-unity/D-selfadjoint-shiab.md, D.1):
  K = Spin(7,7) acting on V = R^{7,7} through SO(7,7) and on the gauge algebra g by conjugation
  (spin representation inside Cl(7,7) (x) C).  By the Clifford grade decomposition g is, as a real
  K-representation, a direct sum of exterior powers Lambda^k V over a set of grades:
      u(64,64)  : k = 0..14                      so(64,64): k = 1,2,5,6,9,10,13,14    (draft (8.5)-(8.6))
      so(8,8)   : k = 1,2,13,14 (= span gamma_a, gamma_ab, gamma_a Gamma, Gamma)
      so(7,8)   : k = 1,2       (= span gamma_a, gamma_ab)          so(7,7): k = 2
  All Lambda^k (k != 7) and Lambda^7_+-  are absolutely irreducible real representations, so real
  dimensions of real Hom spaces equal complex dimensions for so(14,C) = D_7.

Quantities (all exact integers):
  dim S    = dim Hom_K(Lambda^2 (x) g, Lambda^13 (x) g) = dim Hom_K(Lambda^2 (x) g, Lambda^1 (x) g)
  dim S_sa = dim Hom_K(Lambda^3, S^2 g)       (Theorem D2: Shiab o d_B formally self-adjoint)
  dim S_asa= dim Hom_K(Lambda^3, Lambda^2 g)  (Shiab o d_B formally skew-adjoint)
  identity check: dim S = dim S_sa + dim S_asa + dim Hom_K(V (+) Gamma_{21}, g (x) g)
  Def. 8.1:  dim [Lambda^i (x) g]^K for i = 0..14
  Einstein channels: multiplicities of trivial, S^2_0 and Weyl in Lambda^1 (x) Lambda^l.
Method: weight multiplicities (exact), irreducible multiplicities n_lambda = sum_w det(w) m(lambda+rho-w rho)
(Weyl character formula), every decomposition checked against the Weyl dimension formula.
Cross-check: for D_2 (n=4) and D_3 (n=6) the same quantities are computed by brute force (explicit
so(2m,C) action on V(x)V(x)V(x)Lambda^j(x)Lambda^l in a weight basis; invariants = zero-weight vectors killed
by the simple raising operators), with the symmetry conditions of Theorem D2 imposed directly as
linear constraints (this tests Theorem D2 itself, not only the character code).
Exact integer arithmetic except the brute-force ranks (numpy SVD; singular-value gap printed).
"""
import itertools
from fractions import Fraction
from math import comb
import numpy as np


def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


# ============================================================ characters of D_m
class Char:
    """Weight multiplicities on the box [-M, M]^m (integer weights in the e_i basis)."""

    def __init__(self, a, M):
        self.a = a
        self.M = M

    @property
    def m(self):
        return self.a.ndim

    def dim(self):
        return int(self.a.sum())

    def pad(self, M):
        if M == self.M:
            return self.a
        d = M - self.M
        return np.pad(self.a, [(d, d)] * self.m)

    def __add__(self, o):
        M = max(self.M, o.M)
        return Char(self.pad(M) + o.pad(M), M)

    def __mul__(self, o):
        M = self.M + o.M
        out = np.zeros((2 * M + 1,) * self.m, dtype=np.int64)
        w = 2 * o.M + 1
        for idx in np.argwhere(self.a != 0):
            sl = tuple(slice(i, i + w) for i in idx)
            out[sl] += self.a[tuple(idx)] * o.a
        return Char(out, M)

    def adams2(self):
        M = 2 * self.M
        out = np.zeros((2 * M + 1,) * self.m, dtype=np.int64)
        for idx in np.argwhere(self.a != 0):
            out[tuple(2 * (idx - self.M) + M)] += self.a[tuple(idx)]
        return Char(out, M)

    def sym2(self):
        sq, ad = (self * self), self.adams2()
        a = sq.a + ad.a
        assert np.all(a % 2 == 0)
        return Char(a // 2, sq.M)

    def alt2(self):
        sq, ad = (self * self), self.adams2()
        a = sq.a - ad.a
        assert np.all(a % 2 == 0)
        return Char(a // 2, sq.M)


class Dm:
    def __init__(self, m):
        self.m = m
        self.rho = np.arange(m - 1, -1, -1)
        perms = list(itertools.permutations(range(m)))
        P = np.array(perms)
        ps = np.array([perm_sign(p) for p in perms])
        S = np.array([s for s in itertools.product((1, -1), repeat=m) if np.prod(s) == 1])
        wrho = (S[None, :, :] * self.rho[P][:, None, :]).reshape(-1, m)
        self.order = len(wrho)
        self.D = self.rho[None, :] - wrho
        self.eps = np.repeat(ps, len(S))
        self._filt = {}
        # exterior powers of the vector representation C^{2m}
        opts = [(0, 0), (1, 1), (1, -1), (2, 0)]
        self.L = {k: np.zeros((3,) * m, dtype=np.int64) for k in range(2 * m + 1)}
        for ch in itertools.product(range(4), repeat=m):
            k = sum(opts[c][0] for c in ch)
            self.L[k][tuple(opts[c][1] + 1 for c in ch)] += 1

    def ext(self, k):
        return Char(self.L[k].copy(), 1)

    def ext_sum(self, grades):
        c = Char(np.zeros((3,) * self.m, dtype=np.int64), 1)
        for k in grades:
            c = c + self.ext(k)
        return c

    def dominant(self, lam):
        return all(lam[i] >= lam[i + 1] for i in range(self.m - 2)) and lam[self.m - 2] >= abs(lam[self.m - 1])

    def weyl_dim(self, lam):
        l = np.array(lam) + self.rho
        r = self.rho
        num, den = 1, 1
        for i in range(self.m):
            for j in range(i + 1, self.m):
                num *= int((l[i] - l[j]) * (l[i] + l[j]))
                den *= int((r[i] - r[j]) * (r[i] + r[j]))
        q = Fraction(num, den)
        assert q.denominator == 1
        return int(q)

    def decompose(self, ch):
        """{highest weight: multiplicity}; asserts sum of mult*dim = dim."""
        M = ch.M
        if M not in self._filt:
            keep = np.abs(self.D).max(axis=1) <= 2 * M
            self._filt[M] = (self.D[keep], self.eps[keep])
        D, eps = self._filt[M]
        out = {}
        for idx in np.argwhere(ch.a != 0):
            lam = tuple(int(x) for x in idx - M)
            if not self.dominant(lam):
                continue
            pts = np.array(lam)[None, :] + D
            ok = np.all(np.abs(pts) <= M, axis=1)
            vals = ch.a[tuple((pts[ok] + M).T)]
            n = int((eps[ok] * vals).sum())
            assert n >= 0
            if n:
                out[lam] = n
        tot = sum(n * self.weyl_dim(l) for l, n in out.items())
        assert tot == ch.dim(), (tot, ch.dim())
        return out


def hom(dA, dB):
    return sum(n * dB.get(l, 0) for l, n in dA.items())


# ============================================================ brute force for so(2m, C), m = 2, 3
class BruteSO:
    """V = C^{2m}, basis f_1..f_m (weight +e_i), f_{-1}..f_{-m} (weight -e_i); form Q(f_i, f_{-i}) = 1."""

    def __init__(self, m):
        self.m = m
        n = 2 * m
        self.n = n
        e = np.eye(m, dtype=np.int64)
        self.wV = np.array([e[i] for i in range(m)] + [-e[i] for i in range(m)])
        roots = [[(i, i + 1, 1), (m + i + 1, m + i, -1)] for i in range(m - 1)]
        roots.append([(m - 2, 2 * m - 1, 1), (m - 1, 2 * m - 2, -1)])
        Q = np.zeros((n, n), dtype=np.int64)
        for i in range(m):
            Q[i, m + i] = Q[m + i, i] = 1
        for X in roots:
            Xm = np.zeros((n, n), dtype=np.int64)
            for p, q, v in X:
                Xm[p, q] += v
            assert np.array_equal(Xm.T @ Q + Q @ Xm, 0 * Q)
        self.roots = roots
        self.basis = {k: list(itertools.combinations(range(n), k)) for k in range(n + 1)}
        self.index = {k: {S: i for i, S in enumerate(B)} for k, B in self.basis.items()}
        self.wt = {k: np.array([self.wV[list(S)].sum(axis=0) if S else np.zeros(m, dtype=np.int64)
                                for S in B]) for k, B in self.basis.items()}
        # sparse action of each simple root vector on each Lambda^k: list of dicts (col -> {row: val})
        self.act = {}
        for k, B in self.basis.items():
            for r, X in enumerate(roots):
                cols = []
                for S in B:
                    out = {}
                    for p, q, v in X:
                        if q in S and (p == q or p not in S):
                            T = list(S)
                            T[S.index(q)] = p
                            Ts = tuple(sorted(T))
                            j = self.index[k][Ts]
                            out[j] = out.get(j, 0) + v * perm_sign(T)
                    cols.append({j: v for j, v in out.items() if v})
                self.act[(k, r)] = cols

    def invariants(self, blocks, perms):
        """blocks: list of factor-degree tuples (each block = a summand Lambda^{k1}(x)...(x)Lambda^{k5});
        perms: list of (pi, sign) acting on positions, mapping block -> block, constraint
        theta(t) - sign * theta(pi t) = 0.  Returns (nullity, smallest kept / largest dropped singular value)."""
        cols, colidx = [], {}
        for b in blocks:
            W = [self.wt[k] for k in b]
            tot = W[0][:, None, None, None, None, :] + W[1][None, :, None, None, None, :] \
                + W[2][None, None, :, None, None, :] + W[3][None, None, None, :, None, :] \
                + W[4][None, None, None, None, :, :]
            for t in np.argwhere(np.all(tot == 0, axis=-1)):
                key = (b, tuple(int(x) for x in t))
                colidx[key] = len(cols)
                cols.append(key)
        rows = {}
        ent = []
        for c, (b, t) in enumerate(cols):
            for r in range(len(self.roots)):
                for f in range(5):
                    for j, v in self.act[(b[f], r)][t[f]].items():
                        t2 = t[:f] + (j,) + t[f + 1:]
                        key = ('E', r, b, t2)
                        if key not in rows:
                            rows[key] = len(rows)
                        ent.append((rows[key], c, v))
        nE = len(rows)
        R = nE
        for pi, sg in perms:
            for c, (b, t) in enumerate(cols):
                b2 = tuple(b[pi[i]] for i in range(5))
                t2 = tuple(t[pi[i]] for i in range(5))
                c2 = colidx[(b2, t2)]
                ent.append((R, c, 1))
                ent.append((R, c2, -sg))
                R += 1
        A = np.zeros((R, len(cols)))
        for i, j, v in ent:
            A[i, j] += v
        if len(cols) == 0:
            return 0, (np.inf, 0.0)
        s = np.linalg.svd(A, compute_uv=False)
        tol = 1e-8 * max(s[0], 1)
        rank = int((s > tol).sum())
        kept = s[rank - 1] if rank else np.inf
        dropped = s[rank] if rank < len(s) else 0.0
        return len(cols) - rank, (kept, dropped)


def brute_check(m):
    n = 2 * m
    Bf = BruteSO(m)
    Dg = Dm(m)
    V, L2, L3 = Dg.ext(1), Dg.ext(2), Dg.ext(3)
    dVL2 = Dg.decompose(V * L2)
    dL3 = Dg.decompose(L3)
    # theta in V (x) V (x) V (x) Lambda^j (x) Lambda^l : positions (a, b, c, X, Y); Theorem D2 constraints:
    sigma = ((0, 2, 1, 3, 4), -1)      # antisymmetry in (b, c): theta depends on v ^ alpha
    tau = ((2, 1, 0, 4, 3), -1)        # self-adjointness: theta_{cba}(Y, X) = -theta_{abc}(X, Y)
    worst = [np.inf, 0.0]
    mism = 0
    totS = totSa = predS = predSa = 0
    for j in range(n + 1):
        for l in range(j, n + 1):
            blocks = [(1, 1, 1, j, l)] + ([(1, 1, 1, l, j)] if l != j else [])
            nS, g1 = Bf.invariants(blocks, [sigma])
            nSa, g2 = Bf.invariants(blocks, [sigma, tau])
            dJL = Dg.decompose(Dg.ext(j) * Dg.ext(l))
            pS = hom(dVL2, dJL) * len(blocks)
            if l != j:
                pSa = hom(dL3, dJL)
            else:
                pSa = hom(dL3, Dg.decompose(Dg.ext(j).sym2()))
            for g in (g1, g2):
                worst[0] = min(worst[0], g[0])
                worst[1] = max(worst[1], g[1])
            mism += (nS != pS) + (nSa != pSa)
            totS += nS
            totSa += nSa
            predS += pS
            predSa += pSa
    print(f"  D_{m} (n = {n}), g = Lambda^0..{n}: brute force dim S = {totS}, dim S_sa = {totSa}; "
          f"characters + Thm D2: {predS}, {predSa}; block mismatches: {mism} "
          f"(of {2 * (n + 1) * (n + 2) // 2}); SVD gap: smallest kept {worst[0]:.2e}, largest dropped {worst[1]:.1e}")
    return mism


def main():
    print("=== D0. Cross-check of the character code and of Theorem D2 by brute force (D_2, D_3) ===")
    bad = brute_check(2) + brute_check(3)
    assert bad == 0

    print("\n=== D1. D_7 = so(14,C): decompositions used below (all pass the Weyl-dimension check) ===")
    G = Dm(7)
    print(f"  |W(D_7)| = {G.order}")
    L = {k: G.ext(k) for k in range(15)}
    dL = {k: G.decompose(L[k]) for k in range(15)}
    print("  Lambda^k irreducible except k = 7: "
          + ", ".join(f"{k}:{list(dL[k].values())}" for k in range(15)))
    dVL2 = G.decompose(L[1] * L[2])
    print(f"  Lambda^1 (x) Lambda^2 = {dVL2}  (= Lambda^3 + Lambda^1 + Gamma_21)")
    triv, S20, Weyl, L4 = (0,) * 7, (2,) + (0,) * 6, (2, 2) + (0,) * 5, (1, 1, 1, 1, 0, 0, 0)
    w3 = (1, 1, 1, 0, 0, 0, 0)
    dS2L2 = G.decompose(L[2].sym2())
    print(f"  S^2(Lambda^2) = {dS2L2}; algebraic curvature tensors Curv = S^2(Lambda^2) - Lambda^4 "
          f"= trivial + S^2_0 + Weyl (dims 1, {G.weyl_dim(S20)}, {G.weyl_dim(Weyl)})")

    print("\n=== D2. Def. 8.1: dim [Lambda^i (x) u(64,64)]^Spin(7,7) ===")
    dU = G.decompose(G.ext_sum(range(15)))
    inv = [hom(dL[i], dU) for i in range(15)]
    print(f"  i = 0..14: {inv}  (total {sum(inv)}; the draft's 'basis {{Phi_i}}_{{i=0}}^{{14}}' has 15 elements)")

    print("\n=== D3. Size of the admissible space and of its self-adjoint part ===")
    dA = {j: G.decompose(L[2] * L[j]) for j in range(15)}
    dB = {l: G.decompose(L[1] * L[l]) for l in range(15)}
    h3 = {}
    for j in range(15):
        for l in range(j, 15):
            if j == l:
                h3[(j, j, 'S')] = hom({w3: 1}, G.decompose(L[j].sym2()))
                h3[(j, j, 'A')] = hom({w3: 1}, G.decompose(L[j].alt2()))
            else:
                h3[(j, l)] = hom({w3: 1}, G.decompose(L[j] * L[l]))
    algebras = [('u(64,64)', list(range(15))), ('so(64,64)', [1, 2, 5, 6, 9, 10, 13, 14]),
                ('so(8,8) [grades 1,2,13,14]', [1, 2, 13, 14]), ('so(7,8) [grades 1,2]', [1, 2]),
                ('so(7,7) [grade 2]', [2])]
    rest = {l: n for l, n in dVL2.items() if l != w3}
    for name, gr in algebras:
        dimS = sum(hom(dA[j], dB[l]) for j in gr for l in gr)
        sa = sum(h3[(j, j, 'S')] for j in gr) + sum(h3[(j, l)] for j in gr for l in gr if j < l)
        asa = sum(h3[(j, j, 'A')] for j in gr) + sum(h3[(j, l)] for j in gr for l in gr if j < l)
        other = sum(hom(rest, G.decompose(L[j] * L[l])) for j in gr for l in gr)
        assert dimS == sa + asa + other
        dimg = sum(comb(14, k) for k in gr)
        print(f"  g = {name:28s} (dim {dimg:5d}): dim S = {dimS:5d}; dim S_sa = {sa:4d}; dim S_asa = {asa:4d}; "
              f"rest (V + Gamma_21 channels) = {other:5d}; identity holds")
    print("  nonzero blocks of Hom(Lambda^3, S^2 u):  (j,l): mult  [j<l: Lambda^j (x) Lambda^l; (j,j): S^2 Lambda^j]")
    print("   " + ", ".join(f"({k[0]},{k[1]}):{v}" for k, v in sorted(h3.items(), key=lambda x: x[0][:2])
                            if v and (len(k) == 2 or k[2] == 'S')))

    print("\n=== D4. Einstein channels: which outputs can an equivariant map Curv -> Lambda^13 (x) g have? ===")
    rows = []
    for l in range(15):
        rows.append((l, dB[l].get(triv, 0), dB[l].get(S20, 0), dB[l].get(Weyl, 0)))
    print("  l: (mult of trivial, S^2_0, Weyl) in Lambda^1 (x) Lambda^l:")
    print("   " + "; ".join(f"{l}:({a},{b},{c})" for l, a, b, c in rows))
    print(f"  Weyl occurs in Lambda^1 (x) g for no l => every K-equivariant Shiab operator annihilates Weyl "
          f"curvature; trivial and S^2_0 occur only for l in {[l for l, a, b, c in rows if a or b]}")
    print(f"  Hom(Lambda^3, Lambda^1 (x) Lambda^2) = {h3[(1, 2)]};  Hom(Lambda^3, Lambda^2 (x) Lambda^13) = {h3[(2, 13)]};"
          f"  (Lambda^3)^K = {dL[3].get(triv, 0)}")
    k2 = h3[(1, 1, 'S')] + h3[(2, 2, 'S')] + h3[(1, 2)]
    k4 = sum(h3[(j, j, 'S')] for j in (1, 2, 13, 14)) + sum(h3[(j, l)] for j in (1, 2, 13, 14) for l in (1, 2, 13, 14) if j < l)
    print(f"  self-adjoint Shiab operators restricted to so(7,8) x so(7,8): dim {k2};  to so(8,8) x so(8,8): dim {k4}")


if __name__ == '__main__':
    main()
