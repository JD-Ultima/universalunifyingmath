"""GU Phase 2(c)-II: the Shiab operator (draft (8.1), (9.2)-(9.3)) in the actual setting
d = 14, metric eta = (+^7, -^7), ad = u(64,64) realised inside Cl(7,7) (x) C = M_128(C) (clifford77.py).
Claims tested: GU-045, GU-047, GU-048, GU-050, GU-051, GU-053 (claims/geometric-unity.md).

Reconstruction family R-S (formalization/geometric-unity/C-gauge-shiab.md):
   Phi1 = sum_m e^m (x) gamma_m,   Phi2 = sum_{m<n} e^m ^ e^n (x) gamma_m gamma_n       (Spin(7,7)-invariant)
   [a ^ b]_C = sum e^I ^ e^J (a_I b_J - b_J a_I),   [a ^ b]_A = i sum e^I ^ e^J (a_I b_J + b_J a_I)
   Shiab_eps(xi) = [Phi1_eps ^ *xi]_{s1}  -  (c/2) * [Phi1_eps ^ *[Phi2_eps ^ *xi]_{s3}]_{s2},  Phi_eps = eps^-1 Phi eps
   Pairing <T, X> = Tr(T ^ X) (top coefficient), T in Omega^1(ad), X in Omega^13(ad).
Tests
  S0  degree bookkeeping: Shiab maps 2-forms to 13-forms; <T, Shiab(.)> is a 14-form.
  S1  on the Levi-Civita-type curvature xi_R = sum_{a<b} e^ab (x) 1/4 R_abcd gamma^c gamma^d (R an algebraic
      curvature tensor): fit Shiab(xi_R) = *( sum_a e^a (x) sum_b (x Ric_ab + y R eta_ab) gamma^b ); Weyl part -> 0?
  S2  tilted equivariance  Shiab_{eps h}(h^-1 xi h) = h^-1 Shiab_eps(xi) h ; untilted control.
  S3  formal self-adjointness of L = Shiab o d_B on Omega^1(ad): necessary and sufficient pointwise condition
      beta_b(T, a) + beta_b(a, T) = 0 for all b, with beta_b(T, a) = <T, Shiab(e^b ^ a)>.  (CS sanity check in d=3.)
  S4  cubic term: tau(T1,T2,T3) = <T1, Shiab([T2 ^ T3]_C)>; test tau(T,T,a) = kappa tau(a,T,T) for a constant kappa.
Floating point (complex128); relative tolerances printed. Fixed RNG seed.
"""
import itertools
import numpy as np
from clifford77 import clifford_pq

rng = np.random.default_rng(20210401)


# ------------------------------------------------------------ exterior algebra with matrix coefficients
def perm_sign(seq):
    seq = list(seq)
    s = 1
    for i in range(len(seq)):
        for j in range(i + 1, len(seq)):
            if seq[i] > seq[j]:
                s = -s
    return s


class Ext:
    def __init__(self, eta):
        self.eta = np.array(eta)
        self.n = len(eta)
        self._sc = {}

    def star_coeff(self, I):
        if I not in self._sc:
            Ic = tuple(i for i in range(self.n) if i not in I)
            etaI = int(np.prod([self.eta[i] for i in I])) if I else 1
            self._sc[I] = (Ic, etaI * perm_sign(I + Ic))
        return self._sc[I]

    def star(self, f):
        out = {}
        for I, v in f.items():
            Ic, c = self.star_coeff(I)
            out[Ic] = out.get(Ic, 0) + c * v
        return out

    def wedge(self, f, g, prod):
        out = {}
        for I, a in f.items():
            for J, b in g.items():
                if set(I) & set(J):
                    continue
                K = tuple(sorted(I + J))
                s = perm_sign(I + J)
                out[K] = out.get(K, 0) + s * prod(a, b)
        return out

    def top(self, f):
        return f.get(tuple(range(self.n)), 0)


def comm(a, b):
    return a @ b - b @ a


def acomm(a, b):
    return 1j * (a @ b + b @ a)


BR = {'C': comm, 'A': acomm}


class Shiab:
    def __init__(self, ext, gam, eps=None):
        self.E = ext
        n = ext.n
        if eps is None:
            g = gam
        else:
            ei = np.linalg.inv(eps)
            g = [ei @ x @ eps for x in gam]
        self.Phi1 = {(m,): g[m] for m in range(n)}
        self.Phi2 = {(m, k): g[m] @ g[k] for m, k in itertools.combinations(range(n), 2)}

    def term1(self, xi, s1):
        return self.E.wedge(self.Phi1, self.E.star(xi), BR[s1])

    def term2(self, xi, s2, s3):
        inner = self.E.star(self.E.wedge(self.Phi2, self.E.star(xi), BR[s3]))   # 0-form
        return self.E.star(self.E.wedge(self.Phi1, inner, BR[s2]))               # 13-form

    def apply(self, xi, s1, s2, s3, c=1.0):
        t1 = self.term1(xi, s1)
        t2 = self.term2(xi, s2, s3)
        out = dict(t1)
        for K, v in t2.items():
            out[K] = out.get(K, 0) - (c / 2) * v
        return out


def pair(E, T, X):
    """Tr of the top coefficient of T ^ X (T a 1-form, X an (n-1)-form)."""
    t = E.top(E.wedge(T, X, lambda a, b: a @ b))
    return np.trace(t) if isinstance(t, np.ndarray) else complex(t)


def norm_form(f):
    return np.sqrt(sum(np.linalg.norm(v) ** 2 for v in f.values())) if f else 0.0


# ------------------------------------------------------------ curvature tensors
def kn(h, k):
    return (np.einsum('ac,bd->abcd', h, k) + np.einsum('bd,ac->abcd', h, k)
            - np.einsum('ad,bc->abcd', h, k) - np.einsum('bc,ad->abcd', h, k))


def ricci(R, etainv):
    Ric = np.einsum('ac,abcd->bd', etainv, R)
    return Ric, np.einsum('bd,bd->', etainv, Ric)


def weyl(R, eta):
    n = eta.shape[0]
    etainv = np.linalg.inv(eta)
    Ric, s = ricci(R, etainv)
    P = (Ric - s / (2 * (n - 1)) * eta) / (n - 2)
    return R - kn(P, eta)


def random_curvature(n, rng, terms=3):
    R = np.zeros((n,) * 4)
    for _ in range(terms):
        h = rng.integers(-3, 4, (n, n)).astype(float)
        k = rng.integers(-3, 4, (n, n)).astype(float)
        R += kn(h + h.T, k + k.T)
    return R


def spin_curvature(R, gam_up, n):
    xi = {}
    for a, b in itertools.combinations(range(n), 2):
        xi[(a, b)] = 0.25 * sum(R[a, b, c, d] * (gam_up[c] @ gam_up[d]) for c in range(n) for d in range(n) if R[a, b, c, d] != 0)
    return xi


def main():
    n = 14
    pos, neg = clifford_pq(7)
    gam = [g.astype(complex) for g in pos + neg]
    eta_d = [1] * 7 + [-1] * 7
    eta = np.diag(eta_d).astype(float)
    gam_up = [eta_d[c] * gam[c] for c in range(n)]
    E = Ext(eta_d)
    Sh = Shiab(E, gam)
    # B: symmetric split form, u(64,64) = {X : X^H B + B X = 0}
    B = np.eye(128, dtype=complex)
    for a in range(7, 14):
        B = B @ gam[a]
    Binv = np.linalg.inv(B)

    def rand_u():
        M = rng.standard_normal((128, 128)) + 1j * rng.standard_normal((128, 128))
        return 0.5 * (M - Binv @ M.conj().T @ B)
    X = rand_u()
    print(f"sanity: random u(64,64) element satisfies X^H B + B X = 0: "
          f"{np.allclose(X.conj().T @ B + B @ X, 0)}; gamma_m in u(64,64): "
          f"{all(np.allclose(g.conj().T @ B + B @ g, 0) for g in gam)}")

    print("\n=== S0. Degrees ===")
    xi = {(0, 1): rand_u()}
    out = Sh.apply(xi, 'C', 'A', 'A')
    print(f"Shiab(2-form) has components of degrees {sorted(set(len(K) for K in out))}; "
          f"<T, Shiab(xi)> degree = 1 + 13 = 14 = dim Y")

    print("\n=== S1. Shiab on Levi-Civita-type spin curvature (random algebraic curvature tensor) ===")
    R = random_curvature(n, rng)
    W = weyl(R, eta)
    etainv = np.linalg.inv(eta)
    Ric, s = ricci(R, etainv)
    RicW, sW = ricci(W, etainv)
    print(f"|R| = {np.linalg.norm(R):.3e}; Bianchi |R_[abc]d| = "
          f"{np.linalg.norm(R + np.einsum('abcd->bcad', R) + np.einsum('abcd->cabd', R)):.1e}; "
          f"|Ric(W)| = {np.linalg.norm(RicW):.1e}")
    xiR = spin_curvature(R, gam_up, n)
    xiW = spin_curvature(W, gam_up, n)
    U = np.array([sum(Ric[a, b] * gam_up[b] for b in range(n)) for a in range(n)])
    V = np.array([s * gam[a] for a in range(n)])

    def as_vector_valued(Xf):
        """13-form X = *( sum_a e^a (x) M_a ) -> M_a."""
        M = np.zeros((n, 128, 128), dtype=complex)
        for a in range(n):
            Ic, c = E.star_coeff((a,))
            M[a] = Xf.get(Ic, 0) / c
        return M
    for s1 in 'CA':
        t1 = as_vector_valued(Sh.term1(xiR, s1))
        t1W = norm_form(Sh.term1(xiW, s1))
        A = np.stack([U.ravel(), V.ravel()], axis=1)
        coef, *_ = np.linalg.lstsq(A, t1.ravel(), rcond=None)
        resid = np.linalg.norm(A @ coef - t1.ravel()) / max(np.linalg.norm(t1), 1e-300)
        print(f"term1 [Phi1 ^ *xi]_{s1}: |term1| = {np.linalg.norm(t1):.3e}; fit x Ric + y R eta: x = {coef[0]:.6f}, "
              f"y = {coef[1]:.6f}, rel. residual {resid:.1e}; |term1(Weyl)| = {t1W:.1e}")
    for s2, s3 in itertools.product('CA', repeat=2):
        t2 = as_vector_valued(Sh.term2(xiR, s2, s3))
        t2W = norm_form(Sh.term2(xiW, s2, s3))
        nt2 = np.linalg.norm(t2)
        if nt2 < 1e-8 * np.linalg.norm(U):
            print(f"term2 [Phi1 ^ *[Phi2 ^ *xi]_{s3}]_{s2}: identically 0 on Riemann-type xi (|term2| = {nt2:.1e})")
            continue
        A = np.stack([U.ravel(), V.ravel()], axis=1)
        coef, *_ = np.linalg.lstsq(A, t2.ravel(), rcond=None)
        resid = np.linalg.norm(A @ coef - t2.ravel()) / nt2
        print(f"term2 [Phi1 ^ *[Phi2 ^ *xi]_{s3}]_{s2}: |term2| = {nt2:.3e}; fit: x = {coef[0]:.6f}, y = {coef[1]:.6f}, "
              f"rel. residual {resid:.1e}; |term2(Weyl)| = {t2W:.1e}")

    print("\n=== S2. Equivariance ===")

    def cayley(Xu):
        Id = np.eye(128)
        return np.linalg.solve(Id - Xu, Id + Xu)
    eps = cayley(0.05 * rand_u())
    h = cayley(0.05 * rand_u())
    print(f"Cayley elements are B-unitary: {np.allclose(h.conj().T @ B @ h, B)}")
    xi = {K: rand_u() for K in [(0, 1), (2, 9), (5, 13)]}
    hi = np.linalg.inv(h)
    xih = {K: hi @ v @ h for K, v in xi.items()}
    for (s1, s2, s3) in [('C', 'A', 'A'), ('C', 'C', 'C')]:
        lhs = Shiab(E, gam, eps @ h).apply(xih, s1, s2, s3)
        rhs = {K: hi @ v @ h for K, v in Shiab(E, gam, eps).apply(xi, s1, s2, s3).items()}
        err = norm_form({K: lhs.get(K, 0) - rhs.get(K, 0) for K in set(lhs) | set(rhs)}) / norm_form(rhs)
        lhs0 = Shiab(E, gam, None).apply(xih, s1, s2, s3)
        rhs0 = {K: hi @ v @ h for K, v in Shiab(E, gam, None).apply(xi, s1, s2, s3).items()}
        err0 = norm_form({K: lhs0.get(K, 0) - rhs0.get(K, 0) for K in set(lhs0) | set(rhs0)}) / norm_form(rhs0)
        print(f"brackets ({s1},{s2},{s3}): tilted rel. error {err:.1e}; untilted (eps fixed = 1) rel. error {err0:.1e}")

    print("\n=== S3. Formal self-adjointness of Shiab o d (pointwise symbol condition) ===")
    # CS sanity: d = 3, Shiab = identity, L = d : Omega^1 -> Omega^2
    E3 = Ext([1, 1, 1])

    def rand_u2():
        M = rng.standard_normal((2, 2)) + 1j * rng.standard_normal((2, 2))
        return 0.5 * (M - M.conj().T)
    T3 = {(a,): rand_u2() for a in range(3)}
    a3 = {(a,): rand_u2() for a in range(3)}
    defect = []
    for b in range(3):
        eb = {(b,): np.eye(2)}
        bTa = pair(E3, T3, E3.wedge(eb, a3, lambda u, v: u @ v))
        baT = pair(E3, a3, E3.wedge(eb, T3, lambda u, v: u @ v))
        defect.append(abs(bTa + baT) / max(abs(bTa), abs(baT)))
    print(f"CS check (d=3, L=d): max_b |beta_b(T,a)+beta_b(a,T)|/|beta| = {max(defect):.1e}  (0 => self-adjoint)")

    Id = np.eye(128)
    variants = [('term1', 'C', None, None), ('term1', 'A', None, None)] + \
               [('term2', None, s2, s3) for s2, s3 in itertools.product('CA', repeat=2)]

    def grade_elt(k):
        """random real combination of grade-k Clifford monomials (k = 1, 2 lie in so(64,64) < u(64,64))"""
        out = np.zeros((128, 128), dtype=complex)
        for K in itertools.combinations(range(n), k):
            M = np.eye(128, dtype=complex)
            for m in K:
                M = M @ gam[m]
            out += rng.standard_normal() * M
        return out
    sectors = [('u(64,64)', rand_u), ('grade 1 (span gamma_a)', lambda: grade_elt(1)),
               ('grade 2 (spin(7,7))', lambda: grade_elt(2)),
               ('grade 1 + grade 2 (so(7,7) (+) R^{7,7} inside so(64,64))', lambda: grade_elt(1) + grade_elt(2))]
    for sname, gen in sectors:
        T = {(a,): gen() for a in range(n)}
        al = {(a,): gen() for a in range(n)}
        print(f"-- T, a in Omega^1 (x) {sname}")
        for (which, s1, s2, s3) in variants:
            anti, sym, mags = [], [], []
            for b in range(n):
                eb = {(b,): Id}
                xa = E.wedge(eb, al, lambda u, v: u @ v)
                xT = E.wedge(eb, T, lambda u, v: u @ v)
                if which == 'term1':
                    Sa, ST = Sh.term1(xa, s1), Sh.term1(xT, s1)
                else:
                    Sa, ST = Sh.term2(xa, s2, s3), Sh.term2(xT, s2, s3)
                bTa, baT = pair(E, T, Sa), pair(E, al, ST)
                m = max(abs(bTa), abs(baT))
                mags.append(m)
                anti.append(abs(bTa + baT) / max(m, 1e-300))
                sym.append(abs(bTa - baT) / max(m, 1e-300))
            lab = f"[Phi1^*xi]_{s1}" if which == 'term1' else f"[Phi1^*[Phi2^*xi]_{s3}]_{s2}"
            if max(mags) < 1e-9:
                print(f"   {lab:28s}: beta_b = 0 for all b (|beta| < 1e-9)")
                continue
            print(f"   {lab:28s}: max|beta_b| = {max(mags):.2e}; antisymmetry defect max/min over b = "
                  f"{max(anti):.2e}/{min(anti):.2e}; symmetry defect max/min = {max(sym):.2e}/{min(sym):.2e}")
    T = {(a,): rand_u() for a in range(n)}
    al = {(a,): rand_u() for a in range(n)}
    # full operator for the bracket choices selected in S1, all c
    print("full Shiab = term1_C - (c/2) term2_(A,A): defect as a function of c")
    b_vals = range(n)
    t1 = {b: (pair(E, T, Sh.term1(E.wedge({(b,): Id}, al, lambda u, v: u @ v), 'C')),
              pair(E, al, Sh.term1(E.wedge({(b,): Id}, T, lambda u, v: u @ v), 'C'))) for b in b_vals}
    t2 = {b: (pair(E, T, Sh.term2(E.wedge({(b,): Id}, al, lambda u, v: u @ v), 'A', 'A')),
              pair(E, al, Sh.term2(E.wedge({(b,): Id}, T, lambda u, v: u @ v), 'A', 'A'))) for b in b_vals}
    # beta_b(T,a)+beta_b(a,T) = P_b - (c/2) Q_b ; vanishing for all b requires a common c
    P = np.array([t1[b][0] + t1[b][1] for b in b_vals])
    Q = np.array([t2[b][0] + t2[b][1] for b in b_vals])
    cbest = 2 * np.vdot(Q, P) / np.vdot(Q, Q) if np.linalg.norm(Q) > 0 else np.nan
    res = np.linalg.norm(P - (cbest / 2) * Q) / np.linalg.norm(P)
    print(f"  least-squares c = {cbest:.4f}; relative defect at best c = {res:.3e} (0 would mean some c works)")

    print("explicit counterexample (proof in C-gauge-shiab.md, Prop. C5): T = gamma_0 gamma_1 e^0, a = gamma_0 e^1, b = 0")
    Tx = {(0,): gam[0] @ gam[1]}
    ax = {(1,): gam[0]}
    for c in (0.0, 1.0, 2.0, -1.0):
        e0 = {(0,): Id}
        bTa = pair(E, Tx, Sh.apply(E.wedge(e0, ax, lambda u, v: u @ v), 'C', 'A', 'A', c))
        baT = pair(E, ax, Sh.apply(E.wedge(e0, Tx, lambda u, v: u @ v), 'C', 'A', 'A', c))
        print(f"  c = {c:+.0f}: beta_0(T,a) = {bTa.real:+.1f}{bTa.imag:+.1f}i, beta_0(a,T) = {baT.real:+.1f}{baT.imag:+.1f}i")

    print("\n=== S4. Cubic term ===")
    for (s1, s2, s3) in [('C', 'A', 'A'), ('C', 'C', 'C')]:
        ratios = []
        for trial in range(3):
            Tt = {(a,): rand_u() for a in range(n)}
            at = {(a,): rand_u() for a in range(n)}
            TT = E.wedge(Tt, Tt, comm)
            Ta = E.wedge(Tt, at, comm)
            tau_aTT = pair(E, at, Sh.apply(TT, s1, s2, s3))
            tau_TTa = pair(E, Tt, Sh.apply(Ta, s1, s2, s3))
            ratios.append(tau_TTa / tau_aTT)
        print(f"brackets ({s1},{s2},{s3}): tau(T,T,a)/tau(a,T,T) over 3 random samples = "
              + ", ".join(f"{r.real:+.4f}{r.imag:+.1e}i" for r in ratios))
    # CS sanity in d=3
    rs = []
    for trial in range(3):
        Tt = {(a,): rand_u2() for a in range(3)}
        at = {(a,): rand_u2() for a in range(3)}
        tau_aTT = pair(E3, at, E3.wedge(Tt, Tt, comm))
        tau_TTa = pair(E3, Tt, E3.wedge(Tt, at, comm))
        rs.append(tau_TTa / tau_aTT)
    print("CS check (d=3, Shiab = identity): ratios = " + ", ".join(f"{r.real:+.4f}" for r in rs))


if __name__ == '__main__':
    main()
