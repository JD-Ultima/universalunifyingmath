"""GU round 14 (R): the even operator L_{T_o} on all of Lambda^1 (x) g_ev, for T_o in the invariant sector.

Formalization: formalization/geometric-unity/R-even-operator.md.  Builds on I (I2: Upsilon_ev = Shiab(d_0 T_o) + kappa *T_e
+ 2B(T_e, T_o); I2.2: slaving) and J (J1: the vector sector V, the source S and the slaving block M on V).  Imports
J_slaved_sector.py (and through it I_elliptic_reduction.py, G_symbols.py, H_ellipticity.py, shiab_checks.py), unmodified.

Setting (as in I, J).  d = 14, eta = (+^14) and (+^7,-^7); g = Cl(eta) (x) C; T_o = t Phi^1 + u Phi^1 Gamma; family
Shiab = x W_A + y W'_A + z W_C.AdG + w W'_C.AdG (code normalisation).  L_{T_o} T_e := Y_ev[kappa T_e + 2B(T_e, T_o)] (the
1-form Y of J.0, Y_a = eta_aa D_a), a pointwise Spin-equivariant endomorphism of Lambda^1 (x) g_ev (dim 14 * 8192).
Constituents (so(14, C)):  Lambda^1 (x) g_ev = 4 Lam1 + 4 Lam3 + 4 Lam5 + 2 Lam7+ + 2 Lam7- + 2 H2 + 2 H4 + 2 H6,
H_k = [2, 1^(k-1)] (hook) inside V (x) Lambda^k.  Copies (omega in Lam^m, X a hook vector):
   c(omega)_a = grade_{m-1}(gamma_a omega),  w(omega)_a = grade_{m+1}(gamma_a omega),  right multiplication by Gamma.
   Lam^m (m = 1, 3, 5): (c, w, w Gamma, c Gamma) in grades (m-1, m+1, 13-m, 15-m);  for m = 1 these are J's P_0, P_2, P_12, P_14.
   Lam7: (c(om), c(om Gamma), w(om), w(om Gamma)) (grades 6, 6, 8, 8);  H_k: (X, X Gamma) (grades k, 14-k).
Sections
  R0  constituents: Weyl dimensions (D_7), multiplicities, dimension count per grade and in total.
  R1  blocks: L on every copy of every constituent, fitted exactly on the copies of the same abstract vector (two test
      vectors each, both maps, all four operators); templates in sigma = tx + eps uy, theta = tw - uz, delta = ty - ux;
      the Lam1 block equals J's M (J's own j1); the Lam7 block splits by R_Gamma into two 2x2 blocks with equal det.
  R1b equivariance spot check: L(delta X) = delta L(X) for two so(eta) generators, every constituent, all operators.
  R1c the source: Y_ev[Shiab(e^c ^ Phi^1)], Y_ev[Shiab(e^c ^ Phi^1 Gamma)] lie in V exactly, all 14 e^c (= J's S).
  R1d traces over all of Lambda^1 (x) g_ev (orbit representatives under signed permutations): tr L_{op,lab} and, in
      (14,0), all 36 tr(L_a L_b); in (7,7) tr(L^2) at a numeric member; equal to the block predictions.
  R2  determinants: every block, det L = prod det_rho^dim rho, det M a factor, coprimality, factorisations at theta = 0
      and delta = 0, a diagonal-dominance bound.
  R3  singular sets: Shiab# and W_A explicitly; the t-axis; exact kernel vectors verified on Clifford elements; (14,0):
      nonempty for every member; (7,7): empty iff (x, y) = 0 or (z, w) = c (y, x) with c^2 > 25/144 (exact certificate).
  R4  consequences: the even equation off and on the singular set; a flat-base configuration of I7(b) in Z_perp.
Exact arithmetic throughout (Python Fractions, sympy rationals and Gaussian rationals).  No randomness.  One core.
"""
import contextlib
import io
import itertools
import os
import sys
import time
from fractions import Fraction as Fr
from math import comb

for _v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_v, '1')
import sympy as sp  # noqa: E402

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import J_slaved_sector as J  # noqa: E402
import I_elliptic_reduction as I  # noqa: E402
from G_symbols import FULL, N, mask  # noqa: E402

T0 = time.time()
pc = lambda v: bin(v).count('1')  # noqa: E731
P_ = J.P_
SIG = I.SIG
OPN = J.OPN
xs, ys, zs, ws = J.xs, J.ys, J.zs, J.ws
t_, u_, kap = J.t_, J.u_, J.kap
sg_, th_, de_ = sp.symbols('sigma theta delta')
LABS = ('t', 'u')
NCHK = [0]
R = sp.Rational


def ok(label, cond):
    NCHK[0] += 1
    assert cond, label
    P_(f"  [{NCHK[0]}] {label}: True")


def q(v):
    return sp.Rational(v.numerator, v.denominator) if isinstance(v, Fr) else sp.sympify(v)


# ====================================================================== R0
def weyl_dim(lam):
    """dimension of the D_7 irrep with highest weight lam (epsilon coordinates)"""
    rho = [6, 5, 4, 3, 2, 1, 0]
    l = [lam[i] + rho[i] for i in range(7)]
    num = den = Fr(1)
    for i, j in itertools.combinations(range(7), 2):
        num *= (l[i] - l[j]) * (l[i] + l[j])
        den *= (rho[i] - rho[j]) * (rho[i] + rho[j])
    v = num / den
    assert v.denominator == 1
    return int(v)


HW = {'Lam1': (1, 0, 0, 0, 0, 0, 0), 'Lam3': (1, 1, 1, 0, 0, 0, 0), 'Lam5': (1, 1, 1, 1, 1, 0, 0),
      'Lam7+': (1, 1, 1, 1, 1, 1, 1), 'Lam7-': (1, 1, 1, 1, 1, 1, -1), 'H2': (2, 1, 0, 0, 0, 0, 0),
      'H4': (2, 1, 1, 1, 0, 0, 0), 'H6': (2, 1, 1, 1, 1, 1, 0)}
# constituents of V (x) Lambda^k, k even (so(14) rule V (x) Lambda^k = Lambda^(k+1) + Lambda^(k-1) + [2,1^(k-1)],
# V (x) Lambda^6 = Lam7+ + Lam7- + Lam5 + H6, Lambda^k = Lambda^(14-k))
GRADE_CONS = {0: ['Lam1'], 2: ['Lam3', 'Lam1', 'H2'], 4: ['Lam5', 'Lam3', 'H4'], 6: ['Lam7+', 'Lam7-', 'Lam5', 'H6']}
for _k in (8, 10, 12, 14):
    GRADE_CONS[_k] = GRADE_CONS[14 - _k]


def r0():
    P_("=== R0. Constituents of Lambda^1 (x) g_ev (so(14, C); both real forms complexify to it) ===")
    dims = {n: weyl_dim(h) for n, h in HW.items()}
    P_(f"  Weyl dimensions (D_7, rho = (6,...,0)): {dims}")
    ok("Lam^m = C(14, m) for m = 1, 3, 5 and Lam7+- = C(14, 7)/2",
       all(dims['Lam%d' % m] == comb(14, m) for m in (1, 3, 5)) and dims['Lam7+'] == dims['Lam7-'] == comb(14, 7) // 2)
    ok("H_k = 14 C(14, k) - C(14, k+1) - C(14, k-1) for k = 2, 4, 6 (the hook is all of ker(wedge) cap ker(contraction))",
       all(dims['H%d' % k] == 14 * comb(14, k) - comb(14, k + 1) - comb(14, k - 1) for k in (2, 4, 6)))
    ok("dim V (x) Lambda^k = sum of its constituents, every even k",
       all(14 * comb(14, k) == sum(dims[c] for c in GRADE_CONS[k]) for k in GRADE_CONS))
    mult = {}
    for k in GRADE_CONS:
        for c in GRADE_CONS[k]:
            mult[c] = mult.get(c, 0) + 1
    P_(f"  multiplicities: {mult}")
    ok("sum mult * dim = 14 * 2^13 = 114688", sum(mult[c] * dims[c] for c in mult) == 14 * 2 ** 13 == 114688)
    return dims, mult


# ====================================================================== patterns
def gpart(X, g):
    return {k: v for k, v in X.items() if pc(k) == g}


def rG(Al, F):
    return [Al.mul(Xa, {FULL: 1}) if Xa else {} for Xa in F]


def cpat(Al, om, m):
    return [gpart(Al.mul(Al.gamma(a), om), m - 1) for a in range(N)]


def wpat(Al, om, m):
    return [gpart(Al.mul(Al.gamma(a), om), m + 1) for a in range(N)]


def hookvec(Al, p, q_, K, k):
    """X_p = gamma_q E_K, X_q = c gamma_p E_K (p, q not in K, |K| = k - 1), c fixed by sum_a gamma^a X_a = 0"""
    Xp = Al.mul({1 << q_: 1}, {mask(K): 1})
    Xq0 = Al.mul({1 << p: 1}, {mask(K): 1})
    (m1, v1), = gpart(Al.mul({1 << p: Al.eta[p]}, Xp), k + 1).items()
    (m2, v2), = gpart(Al.mul({1 << q_: Al.eta[q_]}, Xq0), k + 1).items()
    assert m1 == m2
    X = [dict() for _ in range(N)]
    X[p] = Xp
    X[q_] = {kk: Fr(-v1, v2) * vv for kk, vv in Xq0.items()}
    return X


def contract_all(Al, X):
    S = {}
    for a in range(N):
        if X[a]:
            S = I.add(S, Al.mul({1 << a: Al.eta[a]}, X[a]))
    return S


TESTS = {'Lam1': [(0,), (9,)], 'Lam3': [(0, 1, 2), (2, 7, 11)], 'Lam5': [(0, 1, 2, 3, 4), (1, 4, 8, 9, 13)],
         'Lam7': [(0, 1, 2, 3, 4, 5, 6), (0, 2, 4, 7, 9, 11, 13)],
         'H2': [(0, 1, (2,)), (3, 9, (12,))], 'H4': [(0, 1, (2, 3, 4)), (6, 8, (1, 10, 13))],
         'H6': [(0, 1, (2, 3, 4, 5, 6)), (2, 12, (0, 5, 7, 9, 11))]}
CONS = ('Lam1', 'Lam3', 'Lam5', 'Lam7', 'H2', 'H4', 'H6')
GR = {'Lam1': (0, 2, 12, 14), 'Lam3': (2, 4, 10, 12), 'Lam5': (4, 6, 8, 10), 'Lam7': (6, 6, 8, 8),
      'H2': (2, 12), 'H4': (4, 10), 'H6': (6, 8)}
COPYNAME = {'Lam1': ('c', 'w', 'wG', 'cG'), 'Lam3': ('c', 'w', 'wG', 'cG'), 'Lam5': ('c', 'w', 'wG', 'cG'),
            'Lam7': ('c', 'cR', 'w', 'wR'), 'H2': ('X', 'XG'), 'H4': ('X', 'XG'), 'H6': ('X', 'XG')}


def patterns(Al, name, tv):
    if name.startswith('Lam') and name != 'Lam7':
        m = int(name[3:])
        om = {mask(tv): 1}
        c, w = cpat(Al, om, m), wpat(Al, om, m)
        return [c, w, rG(Al, w), rG(Al, c)]
    if name == 'Lam7':
        om = {mask(tv): 1}
        omG = Al.mul(om, {FULL: 1})
        return [cpat(Al, om, 7), cpat(Al, omG, 7), wpat(Al, om, 7), wpat(Al, omG, 7)]
    k = int(name[1:])
    p, q_, K = tv
    X = hookvec(Al, p, q_, K, k)
    return [X, rG(Al, X)]


def fit(Y, pats):
    """Y = sum_j s_j pats[j] exactly (pivot on a key private to pattern j, then full residual)"""
    s = []
    for j, F in enumerate(pats):
        piv = None
        for a in range(N):
            for kk, v in F[a].items():
                if all(G[a].get(kk, 0) == 0 for i, G in enumerate(pats) if i != j):
                    piv = (a, kk, v)
                    break
            if piv:
                break
        assert piv is not None
        a, kk, v = piv
        s.append(Fr(Y[a].get(kk, 0)) / Fr(v))
    for a in range(N):
        rec = {}
        for j, F in enumerate(pats):
            if s[j] != 0:
                I.acc_into(rec, F[a], s[j])
        if I.add(rec, Y[a], -1) != {}:
            return s, False
    return s, True


def famsum(d):
    return sp.expand(sum(J.UNIT[op] * J.FAM[op] * q(d[op]) for op in OPN if d.get(op, 0) != 0))


def L_raw(C, X, lab, cache):
    """per operator: Y_ev[2B(X, Phi^1)] (lab t) or Y_ev[2B(X, Phi^1 Gamma)] (lab u), unit-free type combinations"""
    To = C.Phi1 if lab == 't' else C.Phi1G
    po = J.per_op(J.q2B(C.Al, X, To, {'y': cache[lab]}))
    return {op: C.hat(po[op], 3) for op in OPN}


def templates(eps):
    I_ = sp.I
    k, s, th, d = kap, sg_, th_, de_
    T = {}
    T['Lam1'] = sp.Matrix([[k, 0, 208 * eps * d, 0], [0, k + 224 * s, 384 * I_ * eps * th, 16 * eps * d],
                           [16 * d, -384 * I_ * th, k + 544 * s, 0], [0, 208 * d, 0, k]])
    T['Lam3'] = sp.Matrix([[k + R(112, 3) * s, 0, R(352, 3) * eps * d, 64 * I_ * eps * th],
                           [0, k + 240 * s, 320 * I_ * eps * th, 32 * eps * d],
                           [32 * d, -320 * I_ * th, k + 400 * s, 0], [-64 * I_ * th, R(352, 3) * d, 0, k + R(272, 3) * s]])
    T['Lam5'] = sp.Matrix([[k + 96 * s, 0, 48 * eps * d, 128 * I_ * eps * th],
                           [0, k + R(704, 3) * s, 256 * I_ * eps * th, R(80, 3) * eps * d],
                           [R(80, 3) * d, -256 * I_ * th, k + R(832, 3) * s, 0], [-128 * I_ * th, 48 * d, 0, k + 160 * s]])
    T['Lam7'] = sp.Matrix([[k + 176 * s, 0, 0, 192 * I_ * eps * th], [0, k + 176 * s, -192 * I_ * th, 0],
                           [0, 192 * I_ * eps * th, k + 208 * s, 0], [-192 * I_ * th, 0, 0, k + 208 * s]])
    for nm, a, b in (('H2', R(-56, 3), R(-136, 3)), ('H4', -24, -40), ('H6', R(-88, 3), R(-104, 3))):
        T[nm] = sp.Matrix([[k + a * s, -32 * I_ * eps * th], [32 * I_ * th, k + b * s]])
    return T


def sub_std(eps):
    return {sg_: t_ * xs + eps * u_ * ys, th_: t_ * ws - u_ * zs, de_: t_ * ys - u_ * xs}


# ====================================================================== R1
def r1(C, d1):
    Al, eta = C.Al, C.eta
    eps = eta[7]
    P_(f"\n=== R1 {SIG[eps]}. The blocks of L_(T_o) on every constituent (exact; Schur) ===")
    cache = {'t': Al.quad3(C.Phi1), 'u': Al.quad3(C.Phi1G)}
    raw, pats0 = {}, {}
    allfit, same, struct = True, True, True
    for name in CONS:
        ref = None
        for tv in TESTS[name]:
            pats = patterns(Al, name, tv)
            struct &= all(any(F[a] for a in range(N)) for F in pats)
            struct &= all(all(pc(kk) == GR[name][j] for a in range(N) for kk in F[a]) for j, F in enumerate(pats))
            if name.startswith('H'):
                struct &= contract_all(Al, pats[0]) == {}
            if name == 'Lam7':          # c(om) Gamma = w(om Gamma), w(om) Gamma = c(om Gamma)
                struct &= rG(Al, pats[0]) == pats[3] and rG(Al, pats[2]) == pats[1]
            if name == 'Lam1':          # J's P_k(e^c) are the copies of omega = gamma^c
                c = tv[0]
                pk = [C.P(k, C.e(c)) for k in J.KS]
                struct &= all([{kk: v * eta[c] for kk, v in F[a].items()} for a in range(N)] == pk[j]
                              for j, F in enumerate(pats))
            mats = {}
            for lab in LABS:
                for i, X in enumerate(pats):
                    outs = L_raw(C, X, lab, cache)
                    for op in OPN:
                        s, okf = fit(outs[op], pats)
                        allfit &= okf
                        for j in range(len(pats)):
                            mats[(lab, op, j, i)] = s[j]
            same &= ref is None or mats == ref
            if ref is None:
                ref, pats0[name] = mats, pats
        raw[name] = ref
    ok("copies: nonzero, in the stated grades; hook vectors in ker(wedge) cap ker(contraction); Lam7: c(om)G = w(omG),"
       " w(om)G = c(omG); Lam1 copies = J's P_0, P_2, P_12, P_14", struct)
    ok("L(copy_i(v)) = sum_j B_ji copy_j(v) exactly (zero residual), every constituent, both test vectors, t and u maps,"
       " four operators", allfit)
    ok("the block matrices B are the same for both test vectors (index sets inside one signature block and across)", same)
    blocks, T = {}, templates(eps)
    for name in CONS:
        n = len(GR[name])
        Bm = kap * sp.eye(n)
        for lab, var in (('t', t_), ('u', u_)):
            Bm += var * sp.Matrix(n, n, lambda j, i: famsum({op: raw[name][(lab, op, j, i)] for op in OPN}))
        blocks[name] = sp.expand(Bm)
    okT = all(sp.expand(T[nm].subs(sub_std(eps)) - blocks[nm]) == sp.zeros(*blocks[nm].shape) for nm in CONS)
    for nm in CONS:
        P_(f"  {nm} (copies {COPYNAME[nm]}, grades {GR[nm]}): B = {T[nm].tolist()}")
    ok(f"every block equals the printed template in sigma = tx {'+' if eps > 0 else '-'} uy, theta = tw - uz,"
       f" delta = ty - ux (eps = {eps:+d})", okT)
    ok("Lam1 block = J's M (J's own j1); J's S recomputed in R1c", sp.expand(blocks['Lam1'] - d1['M']) == sp.zeros(4, 4))
    # Lam7: split by R_Gamma (right composition with omega -> omega Gamma; Gamma^2 = ss(FULL))
    G2 = Al.cl.sign(FULL, FULL)
    Rm = sp.Matrix([[0, G2, 0, 0], [1, 0, 0, 0], [0, 0, 0, G2], [0, 0, 1, 0]])
    B7 = T['Lam7']
    comm = sp.expand(B7 * Rm - Rm * B7) == sp.zeros(4, 4)
    sub7 = {}
    for lam in (sp.sqrt(G2), -sp.sqrt(G2)):
        Pm = sp.Matrix([[lam, 0], [1, 0], [0, lam], [0, 1]])
        Bl = sp.expand((B7 * Pm).extract([1, 3], [0, 1]))       # rows where the eigenvectors have entry 1
        okinv = sp.expand(B7 * Pm - Pm * Bl) == sp.zeros(4, 2)
        sub7[lam] = (sp.expand(Bl), okinv)
    d7 = [sp.expand(v[0].det()) for v in sub7.values()]
    P_(f"  Gamma^2 = {G2}; Lam7 splits on the R_Gamma-eigenspaces (eigenvalues {list(sub7)}):"
       f" {[v[0].tolist() for v in sub7.values()]}")
    ok("B_Lam7 commutes with R_Gamma; both eigenspaces are B-invariant; det B_7+ = det B_7- = "
       f"{sp.factor(d7[0])}; det B_Lam7 = (det B_7+)^2",
       comm and all(v[1] for v in sub7.values()) and sp.expand(d7[0] - d7[1]) == 0
       and sp.expand(B7.det() - d7[0] * d7[1]) == 0
       and sp.expand(d7[0] - ((kap + 176 * sg_) * (kap + 208 * sg_) - eps * 36864 * th_ ** 2)) == 0)
    return dict(raw=raw, blocks=blocks, T=T, pats=pats0, cache=cache, d7=d7[0], sub7=sub7)


# ====================================================================== R1b
def gen_act(Al, X, p, q_):
    """so(eta) generator (p, q) on a g-valued 1-form: ad(gamma_p gamma_q / 2) on values, rotation on the form index"""
    eta = Al.eta
    h = Al.mul({1 << p: Fr(1, 2)}, {1 << q_: 1})
    out = [Al.comm(h, X[a]) if X[a] else {} for a in range(N)]
    if X[p]:
        out[q_] = I.add(out[q_], X[p], -eta[q_])
    if X[q_]:
        out[p] = I.add(out[p], X[q_], eta[p])
    return out


def r1b(C, res):
    Al = C.Al
    P_(f"\n=== R1b {SIG[C.eta[7]]}. Equivariance spot check (two generators, every constituent, all operators) ===")
    okall = True
    inv = all(gen_act(Al, F, p, q_) == [dict() for _ in range(N)] for F in (C.Phi1, C.Phi1G) for p, q_ in ((0, 1), (2, 9)))
    for name in CONS:
        X = res['pats'][name][0]
        for p, q_ in ((0, 1), (2, 9)):
            dX = gen_act(Al, X, p, q_)
            for lab in LABS:
                A = L_raw(C, dX, lab, res['cache'])
                Bv = L_raw(C, X, lab, res['cache'])
                for op in OPN:
                    okall &= A[op] == gen_act(Al, Bv[op], p, q_)
    ok("Phi^1 and Phi^1 Gamma are invariant, and L_op,lab(delta X) = delta L_op,lab(X) for generators (0,1), (2,9)",
       inv and okall)


# ====================================================================== R1c
def r1c(C, d1):
    Al = C.Al
    P_(f"\n=== R1c {SIG[C.eta[7]]}. The source Shiab(d_0 T_o) lies in V (exact, all 14 basis covectors) ===")
    okall, ref = True, None
    for c in range(N):
        pk = [C.P(k, C.e(c)) for k in J.KS]
        fits = {}
        for lab, P1 in (('t', C.Phi1), ('u', C.Phi1G)):
            po = J.per_op(Al.psi_form(C.vwedge(C.e(c), P1)))
            for op in OPN:
                s, okf = fit(C.hat(po[op]), pk)
                okall &= okf
                fits[(lab, op)] = tuple(s)
        okall &= ref is None or fits == ref
        ref = fits
    S = sp.Matrix([[famsum({op: ref[(lab, op)][j] for op in OPN}) for lab in LABS] for j in range(4)])
    P_(f"  S = {S.tolist()}")
    ok("Y_ev[Shiab(e^c ^ Phi^1)], Y_ev[Shiab(e^c ^ Phi^1 Gamma)] = sum_k P_k(S_k e^c) with zero residual (no component"
       " outside V), the same S for all 14 c; S = J's S", okall and sp.expand(S - d1['S']) == sp.zeros(4, 2))
    return S


# ====================================================================== R1d
def orbit_types(eta):
    """(count, a, K): representatives of (a, K) under the signed permutations in SO(eta)"""
    out = []
    if eta[7] == 1:
        for k in range(0, 15, 2):
            if k >= 1:
                out.append((comb(14, k) * k, 0, tuple(range(k))))
            if k <= 13:
                out.append((comb(14, k) * (14 - k), 13, tuple(range(k))))
        return out
    pos, neg = list(range(7)), list(range(7, 14))
    for k in range(0, 15, 2):
        for p in range(max(0, k - 7), min(7, k) + 1):
            qq = k - p
            K = tuple(pos[:p] + neg[:qq])
            base = comb(7, p) * comb(7, qq)
            if p >= 1:
                out.append((base * p, 0, K))
            if p <= 6:
                out.append((base * (7 - p), 6, K))
            if qq >= 1:
                out.append((base * qq, 7, K))
            if qq <= 6:
                out.append((base * (7 - qq), 13, K))
    return out


def block_pred(res, dims, f):
    """sum over constituents of dim * f(block-matrix function); Lam7 counted with dim 1716 and its 4x4 matrix"""
    tot = 0
    for name in CONS:
        d = dims['Lam7+'] if name == 'Lam7' else dims[name]
        tot += d * f(name)
    return tot


def rawmat(res, name, lab, op):
    n = len(GR[name])
    return sp.Matrix(n, n, lambda j, i: q(res['raw'][name][(lab, op, j, i)]))


def r1d(C, res, dims):
    Al, eta = C.Al, C.eta
    P_(f"\n=== R1d {SIG[eta[7]]}. Traces over all of Lambda^1 (x) g_ev (orbit representatives) ===")
    types = orbit_types(eta)
    ok(f"{len(types)} orbit types of basis elements e^a (x) E_K, counts summing to 114688",
       sum(c for c, _, _ in types) == 114688)
    cache = res['cache']
    tr1 = {(lab, op): Fr(0) for lab in LABS for op in OPN}
    tr2 = {}
    numL2 = Fr(0)
    # numeric member for (7,7): unit * coefficient real: (x, y, z, w) = (2, -1, 3i, i), t = 1/7, u = -1/3
    cnum = {'W_A': Fr(-4), "W'_A": Fr(-2), 'W_C.AdG': Fr(-6), "W'_C.AdG": Fr(-2)}
    tn, un = Fr(1, 7), Fr(-1, 3)
    full = eta[7] == 1
    for cnt, a, K in types:
        mK = mask(K)
        X = [dict() for _ in range(N)]
        X[a] = {mK: 1}
        first = {lab: L_raw(C, X, lab, cache) for lab in LABS}
        for lab in LABS:
            for op in OPN:
                tr1[(lab, op)] += cnt * Fr(first[lab][op][a].get(mK, 0))
        if full:
            for lb in LABS:
                for ob in OPN:
                    Yb = first[lb][ob]
                    if not any(Yb):
                        continue
                    sec = {la: L_raw(C, Yb, la, cache) for la in LABS}
                    for la in LABS:
                        for oa in OPN:
                            key = ((la, oa), (lb, ob))
                            tr2[key] = tr2.get(key, Fr(0)) + cnt * Fr(sec[la][oa][a].get(mK, 0))
        else:
            X1 = [dict() for _ in range(N)]
            for lab, sc in (('t', tn), ('u', un)):
                for op in OPN:
                    for b in range(N):
                        I.acc_into(X1[b], first[lab][op][b], cnum[op] * sc)
            v = Fr(0)
            if any(X1):
                sec = {lab: L_raw(C, X1, lab, cache) for lab in LABS}
                for lab, sc in (('t', tn), ('u', un)):
                    for op in OPN:
                        v += cnum[op] * sc * Fr(sec[lab][op][a].get(mK, 0))
            numL2 += cnt * v
    pred1 = {(lab, op): block_pred(res, dims, lambda nm: rawmat(res, nm, lab, op).trace()) for lab in LABS for op in OPN}
    P_(f"  tr L_(op,lab) direct: {[str(v) for v in tr1.values()]}")
    ok("tr L_(op,lab) = sum_rho dim rho tr B_rho,(op,lab) for all 8 (op, lab); all vanish, so tr L_(T_o) = 114688 kappa",
       all(q(tr1[k]) == pred1[k] == 0 for k in tr1))
    if full:
        pred2 = {key: block_pred(res, dims, lambda nm: (rawmat(res, nm, *key[0]) * rawmat(res, nm, *key[1])).trace())
                 for key in tr2}
        allkeys = [((la, oa), (lb, ob)) for la in LABS for oa in OPN for lb in LABS for ob in OPN]
        for key in allkeys:
            tr2.setdefault(key, Fr(0))
            if key not in pred2:
                pred2[key] = block_pred(res, dims,
                                        lambda nm: (rawmat(res, nm, *key[0]) * rawmat(res, nm, *key[1])).trace())
        nz = sum(1 for k in allkeys if tr2[k] != 0)
        ok(f"tr(L_a L_b) direct (second application) = block prediction for all 64 ordered pairs a, b of (op, lab)"
           f" ({nz} nonzero)", all(q(tr2[k]) == pred2[k] for k in allkeys))
    else:
        def bnum(nm):
            n = len(GR[nm])
            Bm = sp.zeros(n, n)
            for lab, sc in (('t', tn), ('u', un)):
                for op in OPN:
                    Bm += q(cnum[op] * sc) * rawmat(res, nm, lab, op)
            return (Bm * Bm).trace()
        pred = block_pred(res, dims, bnum)
        ok(f"tr((L - kappa)^2) at (x, y, z, w) = (2, -1, 3i, i), (t, u) = (1/7, -1/3): direct = {numL2} = block prediction",
           q(numL2) == pred and numL2 != 0)


# ====================================================================== R2
def r2(res_all, dims, mult):
    P_("\n=== R2. Determinants ===")
    expo = {'Lam1': dims['Lam1'], 'Lam3': dims['Lam3'], 'Lam5': dims['Lam5'], 'Lam7': 2 * dims['Lam7+'],
            'H2': dims['H2'], 'H4': dims['H4'], 'H6': dims['H6']}
    D = {}
    for eps in (1, -1):
        T = res_all[eps]['T']
        D[eps] = {}
        for nm in CONS:
            D[eps][nm] = sp.expand(res_all[eps]['d7']) if nm == 'Lam7' else sp.expand(T[nm].det())
        P_(f"  {SIG[eps]}:")
        for nm in CONS:
            P_(f"    det B_{nm} = {sp.factor(D[eps][nm])}" + ("  (each of Lam7+, Lam7-)" if nm == 'Lam7' else ''))
    size = {nm: (2 if nm == 'Lam7' else len(GR[nm])) for nm in CONS}
    ok("det L_(T_o) = det(M)^14 det(B_3)^364 det(B_5)^2002 det(B_7)^3432 det(B_H2)^896 det(B_H4)^11648 det(B_H6)^36608"
       " (Lam7: det B_7+ = det B_7-, each to the power 1716): the degrees add up to 114688",
       sum(expo[nm] * size[nm] for nm in CONS) == 114688)
    okr = all(all(sp.im(cf) == 0 for cf in sp.Poly(D[e][nm], kap, sg_, th_, de_).coeffs()) for e in D for nm in CONS)
    ok("every block determinant is a real polynomial in (kappa, sigma, theta, delta); at T_o = 0 it is kappa^n", okr and all(
        D[e][nm].subs({sg_: 0, th_: 0, de_: 0}) == kap ** (2 if nm == 'Lam7' else len(GR[nm])) for e in D for nm in CONS))
    eps_ = sp.Symbol('epsilon')
    rho = kap ** 2 - 3328 * eps_ * de_ ** 2
    ok("det M = (rho + 224 kappa sigma)(rho + 544 kappa sigma) - eps (384 kappa theta)^2, rho = kappa^2 - 3328 eps delta^2"
       " (J1(b))", all(sp.expand(D[e]['Lam1'] - ((rho + 224 * kap * sg_) * (rho + 544 * kap * sg_)
                                               - eps_ * (384 * kap * th_) ** 2).subs(eps_, e)) == 0 for e in D))
    two = {'Lam7': (176, 208, 192), 'H2': (R(-56, 3), R(-136, 3), 32), 'H4': (-24, -40, 32), 'H6': (R(-88, 3), R(-104, 3), 32)}
    ok("2x2 blocks: det = (kappa + a sigma)(kappa + b sigma) - eps c^2 theta^2 with (a, b, c) = Lam7 (176, 208, 192),"
       " H2 (-56/3, -136/3, 32), H4 (-24, -40, 32), H6 (-88/3, -104/3, 32)",
       all(sp.expand(D[e][nm] - ((kap + a * sg_) * (kap + b * sg_) - e * c * c * th_ ** 2)) == 0
           for e in D for nm, (a, b, c) in two.items()))
    # factorisations at theta = 0 and at delta = 0 (4x4 blocks)
    diag = {'Lam1': (0, 224, 544, 0), 'Lam3': (R(112, 3), 240, 400, R(272, 3)), 'Lam5': (96, R(704, 3), R(832, 3), 160)}
    Pd = {'Lam1': 3328, 'Lam3': R(11264, 3), 'Lam5': 1280}
    qr = {'Lam1': (0, 384), 'Lam3': (64, 320), 'Lam5': (128, 256)}
    okf = True
    for e in D:
        for nm, (ac, aw, awG, acG) in diag.items():
            dc, dw, dwG, dcG = [kap + a * sg_ for a in (ac, aw, awG, acG)]
            f0 = sp.expand((dc * dwG - e * Pd[nm] * de_ ** 2) * (dw * dcG - e * Pd[nm] * de_ ** 2))
            f1 = sp.expand((dc * dcG - e * qr[nm][0] ** 2 * th_ ** 2) * (dw * dwG - e * qr[nm][1] ** 2 * th_ ** 2))
            okf &= sp.expand(D[e][nm].subs(th_, 0) - f0) == 0 and sp.expand(D[e][nm].subs(de_, 0) - f1) == 0
    ok("4x4 blocks (order c, w, wG, cG): at theta = 0, det = [d_c d_wG - eps P delta^2][d_w d_cG - eps P delta^2]"
       " (P = 3328, 11264/3, 1280); at delta = 0, det = [d_c d_cG - eps q^2 theta^2][d_w d_wG - eps r^2 theta^2]"
       " ((q, r) = (0, 384), (64, 320), (128, 256)); d = kappa + (diagonal) sigma", okf)
    cop = True
    for e in D:
        for nm in CONS:
            if nm != 'Lam1':
                cop &= sp.gcd(D[e]['Lam1'], D[e][nm]) == 1
    ok("det M is coprime to every other block determinant (as polynomials in kappa, sigma, theta, delta)", cop)
    fl = {}
    for e in D:
        for nm in CONS:
            f = sp.factor_list(D[e][nm])[1]
            fl[(e, nm)] = len(f)
    P_(f"  number of irreducible factors over Q (generic member): {fl}")
    # diagonal dominance (Gershgorin): row sums bounded by 544|sigma| + 384|theta| + 208|delta|
    okg = True
    for e in D:
        T = res_all[e]['T']
        for nm in CONS:
            M_ = T[nm]
            for i in range(M_.shape[0]):
                cs = abs(sp.Poly(M_[i, i], sg_).coeff_monomial(sg_))
                ct = sum(abs(sp.Poly(M_[i, j], th_).coeff_monomial(th_)) for j in range(M_.shape[0]) if j != i)
                cd = sum(abs(sp.Poly(M_[i, j], de_).coeff_monomial(de_)) for j in range(M_.shape[0]) if j != i)
                okg &= cs <= 544 and ct <= 384 and cd <= 208
    ok("every row of every block: |diagonal sigma-coefficient| <= 544, sum |theta-coefficients| <= 384, sum |delta-"
       "coefficients| <= 208; so L_(T_o) is invertible whenever 544|sigma| + 384|theta| + 208|delta| < |kappa|", okg)
    # real form (14,0): x = -i x', ..., t = i t', u = i u' leaves sigma, theta, delta unchanged in form
    xp, yp, zp, wp, tp, up = sp.symbols("x' y' z' w' t' u'")
    rf = {xs: -sp.I * xp, ys: -sp.I * yp, zs: -sp.I * zp, ws: -sp.I * wp, t_: sp.I * tp, u_: sp.I * up}
    s1 = sub_std(1)
    ok("(14,0) real form: sigma, theta, delta in (x, y, z, w, t, u) = (-i x', ..., i t', i u') are the same expressions in"
       " the primed real variables", all(sp.expand(s1[v].subs(rf, simultaneous=True)
                                                   - s1[v].subs({xs: xp, ys: yp, zs: zp, ws: wp, t_: tp, u_: up})) == 0
                                         for v in (sg_, th_, de_)))
    okj = all(sp.expand(sp.Matrix([[sp.diff(sub_std(e)[v], w_) for w_ in (t_, u_)] for v in (sg_, th_)]).det()
                        + (xs * zs + e * ys * ws)) == 0 for e in (1, -1))
    ok("d(sigma, theta)/d(t, u) = -(xz + eps yw) = -l (H7's factor): l != 0 iff sigma, theta are coordinates on the (t, u)"
       " plane", okj)
    bsh = {k_: v.subs({xs: 1, ys: 0, zs: 1, ws: 0}) for k_, v in sub_std(1).items()}
    bwa = {k_: v.subs({xs: 1, ys: 0, zs: 0, ws: 0}) for k_, v in sub_std(1).items()}
    ok("Shiab#: (sigma, theta, delta) = (t, -u, -u), W_A: (t, 0, -u) (both real forms); the dominance bound reads"
       " 544|t| + 592|u| < |kappa|, resp. 544|t| + 208|u| < |kappa|",
       all(sub_std(e)[v].subs({xs: 1, ys: 0, zs: 1, ws: 0}) == bsh[v] for e in (1, -1) for v in bsh)
       and (bsh[sg_], bsh[th_], bsh[de_]) == (t_, -u_, -u_) and (bwa[sg_], bwa[th_], bwa[de_]) == (t_, 0, -u_))
    return D


# ====================================================================== R3
def apply_L_num(C, res, X, member, tv, uv, kv):
    """L_(T_o) X for a numeric member (sympy numbers), X with Gaussian-rational coefficients: returns 1-form (sympy)"""
    re_ = [{kk: Fr(int(sp.re(v).p), int(sp.re(v).q)) for kk, v in Xa.items() if sp.re(v) != 0} for Xa in X]
    im_ = [{kk: Fr(int(sp.im(v).p), int(sp.im(v).q)) for kk, v in Xa.items() if sp.im(v) != 0} for Xa in X]
    out = [dict() for _ in range(N)]
    for a in range(N):
        for kk, v in X[a].items():
            out[a][kk] = out[a].get(kk, 0) + kv * v
    for part, ph in ((re_, 1), (im_, sp.I)):
        if not any(part):
            continue
        for lab, sc in (('t', tv), ('u', uv)):
            if sc == 0:
                continue
            Y = L_raw(C, part, lab, res['cache'])
            for op in OPN:
                cf = J.UNIT[op] * member[op] * sc * ph
                if cf == 0:
                    continue
                for a in range(N):
                    for kk, v in Y[op][a].items():
                        out[a][kk] = out[a].get(kk, 0) + cf * q(v)
    return [{kk: sp.expand(v) for kk, v in Oa.items() if sp.expand(v) != 0} for Oa in out]


def combo(pats, vec):
    X = [dict() for _ in range(N)]
    for F, c in zip(pats, vec):
        if c != 0:
            for a in range(N):
                for kk, v in F[a].items():
                    X[a][kk] = X[a].get(kk, 0) + c * q(v)
    return [{kk: sp.expand(v) for kk, v in Xa.items() if sp.expand(v) != 0} for Xa in X]


def r3(res_all, D):
    P_("\n=== R3. Singular sets Z = {det L_(T_o) = 0} in real (t, u); kernels ===")
    SH = {xs: 1, ys: 0, zs: 1, ws: 0}
    WA = {xs: 1, ys: 0, zs: 0, ws: 0}
    for e in (1, -1):
        for nmb, mb in (('Shiab#', SH), ('W_A', WA)):
            sb = {k_: v.subs(mb) for k_, v in sub_std(e).items()}
            P_(f"  {SIG[e]} {nmb}: " + ";  ".join(f"{nm}: {sp.factor(D[e][nm].subs(sb))}" for nm in CONS))
    # the t-axis (u = 0): Shiab# and W_A coincide there (sigma = t, theta = delta = 0)
    diagc = {'Lam1': (0, 224, 544, 0), 'Lam3': (R(112, 3), 240, 400, R(272, 3)), 'Lam5': (96, R(704, 3), R(832, 3), 160),
             'Lam7': (176, 176, 208, 208), 'H2': (R(-56, 3), R(-136, 3)), 'H4': (-24, -40), 'H6': (R(-88, 3), R(-104, 3))}
    pts = {}
    for nm, cs in diagc.items():
        for j, a in enumerate(cs):
            if a != 0:
                pts.setdefault(-1 / sp.Rational(a), set()).add(f"{nm}:{COPYNAME[nm][j]}")
    okax = True
    for e in (1, -1):
        T = res_all[e]['T']
        for nm in CONS:
            Bx = T[nm].subs({th_: 0, de_: 0})
            okax &= Bx.is_diagonal() and all(sp.expand(Bx[i, i] - kap - sp.Rational(a) * sg_) == 0
                                             for i, a in enumerate(diagc[nm]))
    neg = sorted(p for p in pts if p < 0)
    posp = sorted(p for p in pts if p > 0)
    P_(f"  t-axis (u = 0, any member with sigma = t there, e.g. Shiab#, W_A; kappa = 1): singular t = -1/a for every"
       f" diagonal entry kappa + a t;  {len(pts)} values:")
    for p in neg + posp:
        P_(f"    t = {p}: kernel = {sorted(pts[p])}")
    ok("on u = 0 every block of Shiab# and W_A is diagonal with entries kappa + a t (a as listed), both real forms; the"
       " component of invertibility containing 0 is -kappa/544 < t < 3 kappa/136 (kappa > 0), while det M vanishes only at"
       " t = -kappa/224, -kappa/544", okax and neg[-1] == R(-1, 544) and posp[0] == R(3, 136)
       and {p for p in pts if any(s.startswith('Lam1') for s in pts[p])} == {R(-1, 224), R(-1, 544)})
    # exact kernel vectors on Clifford elements
    P_("  exact kernel elements T_e (outside V unless stated) with L_(T_o) T_e = 0, checked on the Clifford algebra:")
    I_ = sp.I
    mSH = {'W_A': 1, "W'_A": 0, 'W_C.AdG': 1, "W'_C.AdG": 0}
    mWA = {'W_A': 1, "W'_A": 0, 'W_C.AdG': 0, "W'_C.AdG": 0}
    mm = sp.Symbol('m')
    tq = sp.simplify((R(704, 25) * mm ** 2 - 1200) / (44800 - 11264 * mm ** 2))
    tq1 = tq.subs(mm, 1)
    uq1 = (tq1 + R(1, 400))
    cases = [  # (signature, member name, member, t, u, constituent, kernel vector)
        (1, 'Shiab#', mSH, R(3, 136), 0, 'H2', [0, 1]),
        (1, 'Shiab#', mSH, R(-1, 400), 0, 'Lam3', [0, 0, 1, 0]),
        (1, 'Shiab#', mSH, 0, R(1, 32), 'H2', [1, I_]),
        (1, 'Shiab#', mSH, 0, R(1, 32), 'H4', [1, I_]),
        (1, 'Shiab#', mSH, 0, R(1, 32), 'H6', [1, I_]),
        (-1, 'Shiab#', mSH, R(3, 136), 0, 'H2', [0, 1]),
        (-1, 'W_A', mWA, R(-1, 176), 0, 'Lam7', [1, 0, 0, 0]),
        (-1, 'W_A', mWA, R(-3, 832), 0, 'Lam5', [0, 0, 1, 0]),
        (1, 'W_A', mWA, tq1, uq1, 'Lam3', None),
        (1, 'W_A', mWA, R(-1, 544) + 0, 0, 'Lam1', [0, 0, 1, 0]),
    ]
    okk = True
    for e, nmb, mb, tv, uv, nm, vec in cases:
        res = res_all[e]
        C = res['C']
        Bm = res['blocks'][nm].subs({xs: mb['W_A'], ys: mb["W'_A"], zs: mb['W_C.AdG'], ws: mb["W'_C.AdG"],
                                     t_: tv, u_: uv, kap: 1})
        if vec is None:
            ns = Bm.nullspace()
            assert len(ns) == 1
            vec = list(ns[0])
        v = sp.Matrix(vec)
        okb = sp.expand(Bm * v) == sp.zeros(Bm.shape[0], 1) and sp.expand(Bm.det()) == 0
        X = combo(res['pats'][nm], vec)
        Y = apply_L_num(C, res, X, mb, tv, uv, 1)
        okc = all(Ya == {} for Ya in Y) and any(X)
        # control: a vector outside the kernel is not annihilated
        ctrl = [1] + [0] * (len(vec) - 1) if vec[0] == 0 else [0] * (len(vec) - 1) + [1]
        Yc = apply_L_num(C, res, combo(res['pats'][nm], ctrl), mb, tv, uv, 1)
        okk &= okb and okc and any(Yc)
        P_(f"    {SIG[e]} {nmb}, kappa = 1, (t, u) = ({tv}, {uv}): {nm}, copy vector {[sp.nsimplify(c) for c in vec]}:"
           f" block kernel {okb}, L T_e = 0 on the Clifford element {okc}")
    ok("each listed T_e is a nonzero element of the stated constituent with L_(T_o) T_e = 0 exactly; a non-kernel copy"
       " is not annihilated (control)", okk)
    sb = {k_: v.subs({xs: 1, ys: 0, zs: 1, ws: 0}) for k_, v in sub_std(1).items()}
    at = {nm: D[1][nm].subs(sb).subs({t_: 0, u_: R(1, 32), kap: 1}) for nm in CONS}
    rk = {nm: res_all[1]['blocks'][nm].subs({xs: 1, ys: 0, zs: 1, ws: 0, t_: 0, u_: R(1, 32), kap: 1}).rank()
          for nm in ('H2', 'H4', 'H6')}
    P_(f"    (14,0) Shiab#, (t, u) = (0, kappa/32): block determinants {at}")
    ok("(14,0) Shiab# at (t, u) = (0, kappa/32): the three hook blocks have rank 1 and the others are invertible, so"
       " ker L_(T_o) = H2 + H4 + H6 (one copy each), of dimension 896 + 11648 + 36608 = 49152, and det M != 0 there",
       all(at[nm] == 0 for nm in ('H2', 'H4', 'H6')) and all(at[nm] != 0 for nm in ('Lam1', 'Lam3', 'Lam5', 'Lam7'))
       and all(v == 1 for v in rk.values()) and 896 + 11648 + 36608 == 49152)
    # (14,0): Z_perp nonempty for every nonzero member (2x2 blocks; sigma = lambda theta case by discriminant)
    lam_, a_, b_, c_ = sp.symbols('lambda a b c')
    disc = sp.discriminant(kap ** 2 + (a_ + b_) * lam_ * kap * th_ + (a_ * b_ * lam_ ** 2 - c_ ** 2) * th_ ** 2, th_)
    ok("(14,0), sigma = lambda theta: discriminant of (kappa + a sigma)(kappa + b sigma) - c^2 theta^2 in theta is"
       " kappa^2 ((a - b)^2 lambda^2 + 4 c^2) > 0", sp.expand(disc - kap ** 2 * ((a_ - b_) ** 2 * lam_ ** 2 + 4 * c_ ** 2)) == 0)
    # (7,7) classification
    P_("  (7,7): classification of the members with Z empty")
    m_ = sp.Symbol('m', nonnegative=True)
    thr = {}
    for nm, (a, b, c) in (('H2', (R(-56, 3), R(-136, 3), 32)), ('H4', (-24, -40, 32)), ('H6', (R(-88, 3), R(-104, 3), 32)),
                          ('Lam7', (176, 208, 192))):
        thr[nm] = sp.Rational((a - b) ** 2, 4 * c * c)
    P_(f"    theta = c sigma, m = c^2: 2x2 block positive definite in (kappa, sigma) iff m > {thr}")
    ok("2x2 thresholds (a - b)^2 / (4 c^2): H2 25/144 (the largest), H4 1/16, H6 1/144, Lam7 1/144",
       thr == {'H2': R(25, 144), 'H4': R(1, 16), 'H6': R(1, 144), 'Lam7': R(1, 144)})
    e = -1
    rr = sp.Symbol('r')
    fM = sp.expand((rr + 224 * kap * sg_) * (rr + 544 * kap * sg_) + 384 ** 2 * m_ * kap ** 2 * sg_ ** 2)
    rmin = sp.solve(sp.diff(fM, rr), rr)[0]
    ok("Lam1 at theta = c sigma, as a quadratic in rho = kappa^2 + 3328 delta^2: minimum kappa^2 sigma^2 (384^2 m - 160^2),"
       " positive for m > 25/144, sigma != 0 (and rho^2 > 0 at sigma = 0)",
       sp.expand(D[e]['Lam1'].subs(th_, sp.sqrt(m_) * sg_) - fM.subs(rr, kap ** 2 + 3328 * de_ ** 2)) == 0
       and sp.expand(fM.subs(rr, rmin) - kap ** 2 * sg_ ** 2 * (384 ** 2 * m_ - 160 ** 2)) == 0)
    m0 = R(25, 144)
    okcert = True
    for nm in ('Lam3', 'Lam5'):
        Dt = D[e][nm]
        even = sp.expand(Dt - Dt.subs(th_, -th_)) == 0 and sp.expand(Dt - Dt.subs(de_, -de_)) == 0
        Pm = sp.expand(9 * Dt.subs({kap: 1, th_: sp.sqrt(m_) * sg_, de_: sp.sqrt(rr)}))
        poly_ok = even and Pm.free_symbols <= {m_, rr, sg_} and sp.Poly(Pm, m_, rr, sg_) is not None
        A2, B1, C0 = [sp.expand(sp.Poly(Pm, rr).coeff_monomial(rr ** i)) for i in (2, 1, 0)]
        P_(f"    {nm} (kappa = 1, theta = c sigma, m = c^2, r = delta^2): 9 det = {A2} r^2 + ({sp.factor(B1)}) r"
           f" + {sp.factor(C0)}")
        # (i) d/dm = sigma^2 (alpha m sigma^2 + beta r + Q(sigma)), alpha, beta >= 0, Q positive definite
        dm = sp.expand(sp.diff(Pm, m_))
        dsm, rem = sp.div(sp.Poly(dm, sg_), sp.Poly(sg_ ** 2, sg_))
        dsm = sp.expand(dsm.as_expr())
        pd = sp.Poly(dsm, m_, rr, sg_)
        mons = dict(zip(pd.monoms(), pd.coeffs()))
        Qs = sp.Poly(dsm.subs({m_: 0, rr: 0}), sg_).all_coeffs()
        ok_i = rem.is_zero and set(mons) <= {(1, 0, 2), (0, 1, 0), (0, 0, 2), (0, 0, 1), (0, 0, 0)} and \
            mons.get((1, 0, 2), 0) >= 0 and mons.get((0, 1, 0), 0) >= 0 and len(Qs) == 3 and Qs[0] > 0 and \
            Qs[1] ** 2 - 4 * Qs[0] * Qs[2] < 0
        P_(f"      d/dm (9 det) = sigma^2 ({sp.factor(dsm)})")
        # (ii) at m0: A r^2 + B r + C >= 0 for r >= 0
        C00 = sp.Poly(sp.expand(C0.subs(m_, m0)), sg_)
        B00 = sp.Poly(sp.expand(B1.subs(m_, m0)), sg_)
        dis = sp.Poly(sp.expand(4 * A2 * C00.as_expr() - B00.as_expr() ** 2), sg_)
        c_nonneg = C00.LC() > 0 and all(mlt % 2 == 0 or sp.Poly(f_, sg_).count_roots() == 0
                                        for f_, mlt in sp.factor_list(C00.as_expr())[1])
        iv = B00.intervals(eps=R(1, 10 ** 6))
        if B00.LC() > 0 and not iv:
            b_case, ok_ii = 'B > 0 everywhere', c_nonneg
        else:
            Lo, Hi = iv[0][0][0], iv[-1][0][1]
            ok_ii = c_nonneg and B00.LC() > 0 and dis.count_roots(Lo, Hi) == 0 and dis.eval(Lo) > 0
            b_case = f"B < 0 only inside [{Lo}, {Hi}], where 4AC - B^2 has no root and is > 0"
        ok_iii = A2 > 0 and B1.subs(sg_, 0) > 0 and C0.subs(sg_, 0) > 0
        okcert &= poly_ok and ok_i and ok_ii and ok_iii
        P_(f"      at m = 25/144: C = {sp.factor(C00.as_expr())} (>= 0), B = {sp.factor(B00.as_expr())}, "
           f"4AC - B^2 = {sp.factor(dis.as_expr())}; {b_case}")
    ok("Lam3, Lam5 at theta = c sigma (det even in theta and in delta): d/dm det = sigma^2 (alpha m sigma^2 + beta r +"
       " Q(sigma)) with alpha, beta >= 0 and Q positive definite; at m = 25/144, A r^2 + B r + C >= 0 on r >= 0; at"
       " sigma = 0 all coefficients > 0. Hence det > 0 for all real (sigma, delta) when m > 25/144", okcert)
    # (x, y) = 0: sigma = delta = 0
    okxy = all(all(cf_ > 0 for cf_ in sp.Poly(D[e][nm].subs({sg_: 0, de_: 0}), kap, th_).coeffs()) for nm in CONS)
    ok("(7,7), (x, y) = 0 (sigma = delta = 0): every block determinant is a positive-coefficient polynomial in kappa^2,"
       " theta^2, so Z is empty", okxy)
    # M1 invariance of sigma, theta, delta
    s_ = sp.Symbol('s')
    okM1 = True
    for e2 in (1, -1):
        if e2 == 1:
            cth, sth = (1 - s_ ** 2) / (1 + s_ ** 2), 2 * s_ / (1 + s_ ** 2)
            Rm = sp.Matrix([[cth, -sth], [sth, cth]])
        else:
            ch, sh = (1 + s_ ** 2) / (1 - s_ ** 2), 2 * s_ / (1 - s_ ** 2)
            Rm = sp.Matrix([[ch, sh], [sh, ch]])
        xy = Rm * sp.Matrix([xs, ys])
        zw = Rm * sp.Matrix([zs, ws])
        tu = Rm * sp.Matrix([t_, u_])
        sb = {xs: xy[0], ys: xy[1], zs: zw[0], ws: zw[1], t_: tu[0], u_: tu[1]}
        for v in (sg_, th_, de_):
            okM1 &= sp.simplify(sub_std(e2)[v].subs(sb, simultaneous=True) - sub_std(e2)[v]) == 0
    sbK = {xs: -ys, ys: -xs, zs: -ws, ws: -zs, t_: u_, u_: t_}
    okM1 &= all(sp.expand(sub_std(-1)[v].subs(sbK, simultaneous=True) - sub_std(-1)[v]) == 0 for v in (sg_, th_, de_))
    ok("sigma, theta, delta (hence every block and det L) are invariant under M1(a) (SO(2) in (14,0), SO+(1,1) in (7,7),"
       " acting on (x,y), (z,w), (t,u)) and under M1(b) ((7,7): p -> -Kp, (t,u) -> (u,t))", okM1)
    # threshold attained: at m = 25/144 H2 has a double root in sigma
    h2 = sp.expand(D[-1]['H2'].subs(th_, R(5, 12) * sg_))
    ok("at c^2 = 25/144 the H2 determinant is (kappa - 32 sigma)^2: Z is nonempty for c^2 <= 25/144",
       sp.expand(h2 - (kap - 32 * sg_) ** 2) == 0)


# ====================================================================== R4
def r4(D):
    P_("\n=== R4. Consequences ===")
    s, cV = sp.symbols('s c_V', real=True)
    # (14,0) Shiab#, flat-base family of I7(b)/N12(c): u = 0, 312 s^2 + kappa s + 15/(2 c_V) = 0
    e1 = 312 * s ** 2 + kap * s + R(15, 2) / cV
    sb = {k_: v.subs({xs: 1, ys: 0, zs: 1, ws: 0}) for k_, v in sub_std(1).items()}
    det = {nm: D[1][nm].subs(sb).subs({t_: s, u_: 0}) for nm in CONS}
    hits = {}
    for nm in CONS:
        for f_ in sp.factor_list(sp.expand(det[nm].subs(kap, 1)))[1]:
            rt = sp.solve(f_[0], s)
            for r_ in rt:
                if -R(1, 312) < r_ < 0:
                    hits[r_] = nm
    P_(f"  (14,0) Shiab#, u = 0, kappa = 1: roots of the I7(b) curve 312 s^2 + s + 15/(2 c_V) = 0 lie in (-1/312, 0);"
       f" singular points there: {hits}")
    cv1 = sp.solve(e1.subs({kap: 1, s: R(-1, 400)}), cV)[0]
    cv2 = sp.solve(e1.subs({kap: 1, s: R(-1, 544)}), cV)[0]
    ok(f"the flat-base family meets Z only at s = -kappa/400 (Lam3, copy wG; c_V kappa^2 = {cv1} >= 9360, det M != 0) and"
       f" s = -kappa/544 (det M = 0; c_V kappa^2 = {cv2})",
       hits == {R(-1, 400): 'Lam3', R(-1, 544): 'Lam1'} and cv1 == R(150000, 11) and cv1 >= 9360 and cv2 == R(277440, 29)
       and det['Lam1'].subs({kap: 1, s: R(-1, 400)}) != 0)
    s2 = sp.solve(e1.subs({kap: 1, cV: cv1}), s)
    other = [r_ for r_ in s2 if r_ != R(-1, 400)]
    ok(f"at c_V kappa^2 = 150000/11 the other root s = {other[0]} kappa is regular (every block determinant != 0)",
       len(other) == 1 and other[0] == R(-11, 15600) and all(det[nm].subs({kap: 1, s: other[0]}) != 0 for nm in CONS))


# ====================================================================== main
def main():
    dims, mult = r0()
    res_all = {}
    for eta in (I.EUC, I.E77):
        C = J.Ctx(eta)
        with contextlib.redirect_stdout(io.StringIO()):
            d1 = J.j1(C)
        res = r1(C, d1)
        res['C'] = C
        r1b(C, res)
        r1c(C, d1)
        r1d(C, res, dims)
        res_all[eta[7]] = res
        P_(f"  [{SIG[eta[7]]} done, {time.time() - T0:.0f} s]")
    D = r2(res_all, dims, mult)
    r3(res_all, D)
    r4(D)
    P_(f"\n{NCHK[0]} checks, {NCHK[0]} True (each asserted)")
    P_(f"runtime {time.time() - T0:.0f} s")


if __name__ == '__main__':
    main()
