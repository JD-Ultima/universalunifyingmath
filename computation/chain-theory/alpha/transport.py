"""Gap CT-G01 / CT-G02: transport the qubit letters to the qutrit copies through the read links,
then credit the join.

Write-up: formalization/chain-theory/alpha/transport.md.

operator-origin.md (its Section 5) left one direction untested. Carry each qubit letter at q_i to
the qutrit copy Q_i through the read link |omega_2>_{Q_i q_i}, by the transpose identity
(A (x) 1)|omega> = (1 (x) A^T)|omega>, and then credit the three matched directions of the join
|omega_3>_{Q1 Q2}. Could the Letter's weights sqrt(k_r) (k = 1, 2, 3, 6), in particular k = 3 and 6
on the 12 qubit-only records, come out of such a construction rather than from a clause?

Parts (every printed number is computed here):
  1. letter level (exact, sympy): t(B) = iota B^T iota^dag is the only traceless Hermitian operator
     at Q_i that acts on the read link as B at q_i; it lives on the 2-dim common support S. The
     subgroup H of G_Omega fixing both read links (dim 4). Up to scale, t is the only H-equivariant
     map from qubit letters to qutrit letters, and the only G_Omega-equivariant one is 0.
  2. operator level: six E'-covariant transports on B(H) (identity value b in {0,1} at the third
     direction, three orderings). They are H-equivariant and not G_Omega-equivariant. Their action
     on the 172 records; collisions.
  3. fiber obstruction: colliding records tie diag(T^dag D T) and the Rayleigh quotients of D, for
     every superoperator D (checked with random D)
  4. what transported qubit-only records see of the join: the restricted join, the invariant
     projectors of H and of SU(3), the amplitudes J1 = sqrt3, J2 = 2/sqrt3, J3 = sqrt2; the
     H-orbit credit applied to every record
  5. family F3: frame-local survival x attribution dressings, 338 688 ordered choices, exact
  6. family F3a: words of length <= 3 in the 24 covariant primitives of operator_origin.py and 9
     transport letters, with the tests of operator_origin.py (operator, diagonal, per-isotypic
     traces; also twirled, |C|, sqrt C)
  7. summary

The definitions of operator_origin.py (which imports records_independent.py and symmetry.py) are
loaded by executing its source up to the start of its F1 enumeration. The diagnostic loop of its
section 1 is left out; it draws no random numbers, so the random state is the same. None of the
three files is modified, and the F1/F2 searches are not rerun.
Run: python3 transport.py > transport.out    (~13 min on 4 shared cores; part 6 uses 4 worker processes)
"""
import os
for _v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")         # one BLAS thread per process; part 6 runs 4 processes
import io, contextlib, itertools, sys, time, types, hashlib, multiprocessing
from fractions import Fraction
from collections import Counter, defaultdict
import numpy as np
import sympy as sp

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
NPROC = min(4, os.cpu_count() or 1)


def load_operator_origin():
    p = os.path.join(HERE, "operator_origin.py")
    src = open(p).read()
    s1 = src.index("# ---------------------------------------------------------------- 1. primitives")
    s2 = src.index("# ---------------------------------------------------------------- 2. symmetry structure on R")
    cut = src.index("\nt_words = time.time()")
    sec1 = src[s1:s2]
    code = src[:s1] + sec1[:sec1.index("FACT, DEN = {}, {}")] + src[s2:cut]
    mod = types.ModuleType("operator_origin_defs")
    mod.__file__ = p
    with contextlib.redirect_stdout(io.StringIO()):
        exec(compile(code, p, "exec"), mod.__dict__)
    return mod


oo = load_operator_origin()
rec, sym = oo.rec, oo.sym
nrm = np.linalg.norm
TOL = oo.TOL
B3, B, Bh, NR = oo.B3, oo.B, oo.Bh, oo.NR
cls, rtype, epar = oo.cls, oo.rtype, oo.epar
CLASSES, KOF = oo.CLASSES, oo.KOF
ML, KL = oo.ML, oo.KL
N = 31104
AINV = 137.035999177                       # CODATA 2022, as in operator_origin.py
rng = np.random.default_rng(20261031)
print(f"[setup] operator_origin.py definitions loaded in {time.time() - T_START:.0f} s", file=sys.stderr)


def hdr(s):
    print("\n" + "=" * 78 + "\n" + s + "\n" + "=" * 78)


GMNAME = {0: 1, 1: 2, 6: 3, 2: 4, 3: 5, 4: 6, 5: 7, 7: 8}   # rec.GM index -> standard lambda number
iota = np.zeros((3, 2)); iota[0, 0] = iota[1, 1] = 1
PS = np.diag([1.0, 1.0, 0.0])                                  # projector onto the common support S
E22 = np.diag([0.0, 0.0, 1.0])
t_letter = lambda Bm: iota @ Bm.T @ iota.T

# ================================================================== 1. letter level
hdr("1. Letter level: carrying a qubit letter to the qutrit through a read link (exact)")
asym = sp.symbols('a0:9'); A_ = sp.Matrix(3, 3, asym)
bsym = sp.symbols('b0:4'); B_ = sp.Matrix(2, 2, bsym)
phi = sp.zeros(6, 1); phi[0] = phi[3] = 1                      # |00> + |11> on C^3_Q (x) C^2_q
eqs = list(sp.kronecker_product(A_, sp.eye(2)) * phi - sp.kronecker_product(sp.eye(3), B_) * phi)
sol = sp.solve(eqs, [asym[i] for i in (0, 1, 3, 4, 6, 7)], dict=True)
assert len(sol) == 1
Agen = A_.subs(sol[0])
print("(A (x) 1)|omega_2> = (1 (x) B)|omega_2>  <=>  A =", Agen.tolist())
print("  i.e. A = B^T on S = span{|0>,|1>}, zero from S into |2>; the column of |2> (a2, a5, a8) is not fixed by the link")
x_, y_, z_ = sp.symbols('x y z', real=True)
pr = sp.symbols('p0:6', real=True)
Bherm = sp.Matrix([[z_, x_ - sp.I * y_], [x_ + sp.I * y_, -z_]])
Ah = Agen.subs({bsym[0]: Bherm[0, 0], bsym[1]: Bherm[0, 1], bsym[2]: Bherm[1, 0], bsym[3]: Bherm[1, 1],
                asym[2]: pr[0] + sp.I * pr[1], asym[5]: pr[2] + sp.I * pr[3], asym[8]: pr[4] + sp.I * pr[5]})
conds = [sp.expand(c) for c in list(Ah - Ah.H) + [Ah.trace()]]
solh = sp.solve([sp.re(c) for c in conds] + [sp.im(c) for c in conds], pr, dict=True)
assert len(solh) == 1 and all(v == 0 for v in solh[0].values()) and len(solh[0]) == 6
tB = sp.simplify(Ah.subs(solh[0]))
assert sp.simplify(tB - sp.Matrix(iota.astype(int)) * Bherm.T * sp.Matrix(iota.astype(int)).T) == sp.zeros(3, 3)
print("B traceless Hermitian, A required traceless Hermitian: unique solution A = t(B) = iota B^T iota^dag =", tB.tolist())
A1 = Agen.subs({bsym[0]: 1, bsym[1]: 0, bsym[2]: 0, bsym[3]: 1, asym[2]: 0, asym[5]: 0})
print("B = 1, A Hermitian: A =", A1.tolist(), " -> the value b = a8 at the third direction is free (b = 0: support")
print("  projector, what the link itself gives; b = 1: identity)")
for nm, P in zip(('x', 'y', 'z'), rec.PA):
    tp = t_letter(P)
    hit = [(GMNAME[j], s) for j, G in enumerate(rec.GM) for s in (1, -1) if np.allclose(tp, s * G)]
    print(f"  t(sigma_{nm}) = {'+' if hit[0][1] > 0 else '-'}lambda_{hit[0][0]};  PS t PS = t: {np.allclose(PS @ tp @ PS, tp)};"
          f"  t |2> = 0: {np.allclose(tp[:, 2], 0) and np.allclose(tp[2, :], 0)}")
print("  so the transported qubit letters are the Gell-Mann letters lambda_1, -lambda_2, lambda_3 of the su(2) on S")

# --- the subgroup H of G_Omega fixing both read links
gens = list(sym.K3) + list(sym.K2)          # 11 generators of G_Omega on H (U (x) V (x) Ubar (x) Vbar)
cols = [np.array([(K @ P - P @ K).reshape(-1) for K in gens]).T for P in (oo.PL1, oo.PL2)]
Mc = np.vstack(cols); Mr = np.vstack([Mc.real, Mc.imag])
_, sv, vh = np.linalg.svd(Mr)
HNULL = vh[(sv > 1e-9).sum():]
print("\nLie algebra of the subgroup H of G_Omega = SU(3) x SU(2) fixing both read links P_L1, P_L2: dim", HNULL.shape[0])
assert HNULL.shape[0] == 4
XQ = [sum(c[a] * rec.GM[a] / 2 for a in range(8)) for c in HNULL]           # action at Q1 (U side)
Xq = [sum(c[8 + m] * rec.PA[m] / 2 for m in range(3)) for c in HNULL]       # action at q1 (V side)
q_rank = np.linalg.matrix_rank(np.array([[np.trace(X @ P).real for P in rec.PA] for X in Xq]))
rel = max(nrm(X[2, :2]) + nrm(X[:2, 2]) + nrm((iota.T @ X @ iota + Y.T) - np.trace(iota.T @ X @ iota + Y.T) / 2 * np.eye(2))
          for X, Y in zip(XQ, Xq))
lam8 = rec.GM[7]
Qproj = np.array([[np.trace(X @ P).real for P in rec.PA] for X in Xq])            # (4 basis elements) x (3 Paulis)
_, s_k, v_k = np.linalg.svd(Qproj.T)
dker = v_k[(s_k > 1e-9).sum():]
Xk = sum(d * X for d, X in zip(dker[0], XQ))
k8 = nrm(Xk - np.trace(Xk @ lam8) / 2 * lam8) / nrm(Xk)
print("  projection of h onto su(2)_q has rank", q_rank, "; every element has U|_S = c Vbar and U|2> in C|2>"
      f" (residual {rel:.1e}); the kernel of the projection ({len(dker)}-dim) is lambda_8 at the qutrit (defect {k8:.1e})")
print("  so H = {(U, V): U = c Vbar (+) c^-2} = SU(2) x U(1), the diagonal su(2) on S plus lambda_8 (cf. Supplement S3, Thm S-stabilizer)")
print("  E' exchanges the two read links: ||E'(P_L1) - P_L2|| =", f"{nrm(oo.Ep(oo.PL1[:, :, None])[:, :, 0] - oo.PL2):.1e}")


def commutant_dim(mats, d):
    rows = [np.kron(X, np.eye(d)) - np.kron(np.eye(d), X.T) for X in mats]
    return d * d - np.linalg.matrix_rank(np.vstack(rows), tol=1e-9)


print("  commutant on C^3 of h:", commutant_dim(XQ, 3), "(span{PS, |2><2|}: invariant subspaces 0, S, C|2>, C^3);"
      "  of su(3):", commutant_dim([G / 2 for G in rec.GM], 3), "(invariant subspaces 0, C^3)")


def ad_matrix_basis(X, basis):
    return np.array([[np.trace(Bi.conj().T @ (X @ Bj - Bj @ X)) / np.trace(Bi.conj().T @ Bi) for Bj in basis] for Bi in basis])


def equivariant_maps(pairs):
    """linear maps F: span(Pauli) -> span(Gell-Mann) with ad_XQ F = F ad_Xq for all generator pairs"""
    rows = []
    for XQ_, Xq_ in pairs:
        aQ = ad_matrix_basis(XQ_, rec.GM); aq = ad_matrix_basis(Xq_, rec.PA)
        rows.append(np.kron(np.eye(3), aQ) - np.kron(aq.T, np.eye(8)))      # column-major vec(F)
    M_ = np.vstack(rows)
    _, s_, v_ = np.linalg.svd(M_)
    return v_[(s_ > 1e-9).sum():].conj()


GOm_pairs = [(G / 2, np.zeros((2, 2))) for G in rec.GM] + [(np.zeros((3, 3)), P / 2) for P in rec.PA]
eq_G = equivariant_maps(GOm_pairs)
eq_H = equivariant_maps(list(zip(XQ, Xq)))
Ft = np.array([[np.trace(G @ t_letter(P)).real / 2 for P in rec.PA] for G in rec.GM])
vt = Ft.reshape(-1, order='F'); vt = vt / nrm(vt)
in_span = nrm(vt - eq_H.T @ (eq_H.conj() @ vt)) if len(eq_H) else 1.0
print(f"equivariant maps qubit letters -> qutrit letters: under G_Omega {len(eq_G)} (Schur: (1,3) -> (8,1)),"
      f" under H {len(eq_H)}; t spans the H-equivariant ones: {in_span < 1e-9}")
print("  => the transport is, up to scale, the only H-equivariant one, and no nonzero transport is G_Omega-equivariant")

# ================================================================== 2. operator level
hdr("2. Operator level: transports on B(H) and what they do to the 172 records")
ROWS, COLS = 'ABCD', 'EFGH'


def tb_tensor(b):
    """t_b(E_ac)[x, y], shape (3,3,2,2); t_b(B) = iota B^T iota^dag + b (tr B / 2) |2><2|"""
    t = np.zeros((3, 3, 2, 2))
    for a in range(2):
        for c in range(2):
            t[c, a, a, c] = 1.0
        t[2, 2, a, a] += b / 2.0
    return t


def link_transport(i, order, b):
    """read link i (0, 1): A (x) B at (Q_i, q_i) -> A t_b(B) ('R'), t_b(B) A ('L'), their mean ('S');
    the identity is left at q_i. Returns (map, Hilbert-Schmidt adjoint) acting on (36,36,n) batches."""
    s, u = 2 * i, 2 * i + 1
    tb = tb_tensor(b)
    ins = ROWS + COLS + 'n'

    def one(ordr):
        if ordr == 'R':
            tbs = COLS[s] + 'z' + ROWS[u] + COLS[u]
            outr = ROWS.replace(ROWS[u], 'j'); outc = COLS.replace(COLS[s], 'z').replace(COLS[u], 'k')
        else:
            tbs = 'z' + ROWS[s] + ROWS[u] + COLS[u]
            outr = ROWS.replace(ROWS[s], 'z').replace(ROWS[u], 'j'); outc = COLS.replace(COLS[u], 'k')
        spec = ins + ',' + tbs + ',jk->' + outr + outc + 'n'
        spec_a = outr + outc + 'n,' + tbs + ',jk->' + ins

        def f(Y):
            T = Y.reshape((3, 2, 3, 2, 3, 2, 3, 2, Y.shape[-1]))
            return np.einsum(spec, T, tb, np.eye(2), optimize=True).reshape(36, 36, -1)

        def fa(Y):
            Z = Y.reshape((3, 2, 3, 2, 3, 2, 3, 2, Y.shape[-1]))
            return np.einsum(spec_a, Z, tb, np.eye(2), optimize=True).reshape(36, 36, -1)
        return f, fa
    if order in ('R', 'L'):
        return one(order)
    (fR, aR), (fL, aL) = one('R'), one('L')
    return (lambda Y: 0.5 * (fR(Y) + fL(Y))), (lambda Y: 0.5 * (aR(Y) + aL(Y)))


def full_transport(o1, o2, b):
    f1, a1 = link_transport(0, o1, b)
    f2, a2 = link_transport(1, o2, b)
    return (lambda Y: f1(f2(Y))), (lambda Y: a2(a1(Y)))


Xt = rng.normal(size=(36, 36, 3)) + 1j * rng.normal(size=(36, 36, 3))
Yt = rng.normal(size=(36, 36, 3)) + 1j * rng.normal(size=(36, 36, 3))
for i, (PL, u) in enumerate(((oo.PL1, 1), (oo.PL2, 3))):
    fR, _ = link_transport(i, 'R', 0); fL, _ = link_transport(i, 'L', 0)
    dR = nrm(fR(Xt) - 4 * oo.TAU[u](oo.R_b(PL)(Xt))) / nrm(fR(Xt))
    dL = nrm(fL(Xt) - 4 * oo.TAU[u](oo.L_b(PL)(Xt))) / nrm(fL(Xt))
    print(f"read link {i + 1}: transport with b = 0 equals 2 Tr_q[X P_L] (x) 1 = 4 tau_q(X P_L) (R) and 4 tau_q(P_L X) (L):"
          f" defects {dR:.1e}, {dL:.1e}")
cH = rng.normal(size=HNULL.shape[0]) @ HNULL
Kh = sum(ci * K for ci, K in zip(cH, gens))
wK, VK = np.linalg.eigh(Kh)
gH = oo.Ad_b((VK * np.exp(1j * wK)) @ VK.conj().T)       # a random element of H acting on B(H)
TRANSPORTS = {}
print(f"{'transport':12s} {'adjoint ok':>10s} {'E-cov':>8s} {'H-equiv':>8s} {'G_Om-equiv':>10s}")
for b in (0, 1):
    for o1, o2 in (('R', 'L'), ('L', 'R'), ('S', 'S'), ('R', 'R')):
        f, fa = full_transport(o1, o2, b)
        adj = abs(np.vdot(Yt, f(Xt)) - np.vdot(fa(Yt), Xt)) < 1e-9 * abs(np.vdot(Yt, f(Xt)))
        dE = nrm(oo.Ep(f(Xt)) - f(oo.Ep(Xt))) / nrm(f(Xt))
        dH = nrm(f(gH(Xt)) - gH(f(Xt))) / nrm(f(Xt))
        dG = nrm(f(oo.gW(Xt)) - oo.gW(f(Xt))) / nrm(f(Xt))
        nm = f"T_b{b}({o1},{o2})"
        print(f"{nm:12s} {str(adj):>10s} {dE:8.1e} {dH:8.1e} {dG:10.2f}")
        if dE < 1e-12:
            TRANSPORTS[nm] = (f, fa)
print("E'-covariant transports kept:", list(TRANSPORTS), "((R,R) fails E'-covariance, since T(A t(B)) = t(B)^T A^T turns R into L)")
REF = "T_b1(R,L)"

Q_OCC = {o: oo.occ_proj(o) for o in [(), (0,), (2,), (0, 2)]}
q_occ = oo.lin((1, oo.ID), (-1, oo.compose(oo.TAU[1], oo.TAU[3])))      # [N_q >= 1]
IMG = {}
for nm, (f, fa) in TRANSPORTS.items():
    Y = f(B3)
    IMG[nm] = Y
ztypes = ['QQ', 'qq', 'mix', 'comm']
same_link = oo.lin((1, oo.occ_proj((0, 1))), (1, oo.occ_proj((2, 3))))     # letters at the two ends of one read link
is_same = [nrm(same_link(B3[:, :, r:r + 1]) - B3[:, :, r:r + 1]) < 1e-9 for r in range(NR)]
# mix-same: letters at (Q_i, q_i); mix-cross: letters at (Q_i, q_j), i != j
subtype = np.array([('mix-same' if is_same[r] else 'mix-cross') if rtype[r] == 'mix' else rtype[r] for r in range(NR)])
print("record sub-types:", dict(Counter(subtype)))


def describe(Y):
    nY = nrm(Y)
    if nY < 1e-12:
        return 'zero', None, None
    pat = tuple(o for o in [(), (0,), (2,), (0, 2)] if nrm(Q_OCC[o](Y[:, :, None])) > 1e-9 * nY)
    v = Y.reshape(-1); c = Bh @ v
    leak = nrm(v - B @ c) / nY
    j = int(np.argmax(np.abs(c)))
    match = j if abs(abs(c[j]) - nY) < 1e-9 * nY else None
    return pat, leak, match


OCCN = {(): 'scalar', (0,): 'Q1', (2,): 'Q2', (0, 2): 'Q1Q2'}
for nm, Y in IMG.items():
    qmax = max(nrm(q_occ(Y[:, :, r:r + 1])) for r in range(NR))
    epd = max(nrm(oo.Ep(Y[:, :, r:r + 1])[:, :, 0] - epar[r] * Y[:, :, r]) for r in range(NR))
    tab = Counter()
    for r in range(NR):
        pat, leak, match = describe(Y[:, :, r])
        if pat == 'zero':
            tab[(cls[r], subtype[r], 'zero')] += 1
            continue
        where = ('= record (' + cls[match] + ',' + subtype[match] + ')') if match is not None else ('in R' if leak < 1e-9 else 'leaves R')
        tab[(cls[r], subtype[r], '+'.join(OCCN[o] for o in pat) + ' ' + where)] += 1
    print(f"\n{nm}: max q-occupation of the images {qmax:.1e}; E'-parity kept (max defect {epd:.1e});"
          f" rank of T on R {np.linalg.matrix_rank(Y.reshape(1296, NR), tol=1e-9)} of {NR}")
    if nm in (REF, "T_b0(R,L)", "T_b1(S,S)"):
        for key, n in sorted(tab.items()):
            print(f"    {key[0]:7s} {key[1]:10s} -> {key[2]:55s} x{n}")


def collisions(Y):
    nY = nrm(Y.reshape(1296, NR), axis=0)
    idx = np.flatnonzero(nY > 1e-12)
    Yn = Y.reshape(1296, NR)[:, idx] / nY[idx]
    G = np.abs(Yn.conj().T @ Yn)
    out = []
    for a, bb in zip(*np.nonzero(np.triu(G, 1) > 1 - 1e-9)):
        r, s = idx[a], idx[bb]
        out.append((r, s, (nY[r] / nY[s]) ** 2))
    return out


COLL = {}
for nm, Y in IMG.items():
    cl = collisions(Y)
    COLL[nm] = cl
    summ = Counter()
    for r, s, kap in cl:
        a_, b_ = (r, s) if subtype[r] >= subtype[s] else (s, r)
        kk = (np.linalg.norm(Y[:, :, a_]) / np.linalg.norm(Y[:, :, b_])) ** 2
        summ[(cls[a_] + '/' + subtype[a_], cls[b_] + '/' + subtype[b_], str(Fraction(kk).limit_denominator(100)))] += 1
    print(f"\ncollisions T r || T r' (r != r') under {nm}: {len(cl)} pairs, {len(summ)} kinds"
          + (":" if nm in (REF, "T_b0(R,L)") else " (table printed for T_b1(R,L) and T_b0(R,L) only)"))
    if nm in (REF, "T_b0(R,L)"):
        for key, n in sorted(summ.items()):
            print(f"    {key[0]:18s} with {key[1]:18s} kappa = |T r|^2/|T r'|^2 = {key[2]:6s} x{n}")

# ================================================================== 3. fiber obstruction
hdr("3. Fiber obstruction: a transported frame cannot tell colliding records apart")
Dr = (rng.normal(size=(1296, 1296)) + 1j * rng.normal(size=(1296, 1296))) / 36
for nm in (REF, "T_b0(S,S)"):
    TBm = IMG[nm].reshape(1296, NR)
    CD = TBm.conj().T @ Dr @ TBm                 # <T r, D T r'>
    nT2 = np.real(np.diag(TBm.conj().T @ TBm))
    worst_d = worst_r = 0.0
    for r, s, kap in COLL[nm]:
        worst_d = max(worst_d, abs(CD[r, r] - kap * CD[s, s]) / abs(CD[r, r]))
        worst_r = max(worst_r, abs(CD[r, r] / nT2[r] - CD[s, s] / nT2[s]) / abs(CD[r, r] / nT2[r]))
    print(f"{nm}: random D: max |diag(T^+ D T)[r] - kappa diag[r']| / |diag[r]| = {worst_d:.1e};"
          f" max difference of the Rayleigh quotients of D at T r and T r' = {worst_r:.1e}")
kap_pairs = defaultdict(set)
for r, s, kap in COLL[REF]:
    a_, b_ = (r, s) if subtype[r] >= subtype[s] else (s, r)
    kap_pairs[(cls[a_] + '/' + subtype[a_], cls[b_] + '/' + subtype[b_], KOF[cls[a_]], KOF[cls[b_]])].add(
        Fraction((nrm(IMG[REF][:, :, a_]) / nrm(IMG[REF][:, :, b_])) ** 2).limit_denominator(100))
print(f"colliding records r, r' under {REF}: for every D, diag(T^+ D T)[r] / diag[r'] = kappa and the Rayleigh quotients are equal")
print(f"   {'r':18s} {'r_prime':18s} {'kappa':14s} {'M_L ratio':>9s} {'K_L ratio':>9s}  breaks: diag(T^+DT)=M_L, =K_L; Rayleigh=M_L, =K_L")
vM = vK = vR = False
for (ca, cb, ka, kb), v in sorted(kap_pairs.items()):
    bM = any(abs(float(x) - (ka / kb) ** 0.5) > 1e-9 for x in v)
    bK = any(abs(float(x) - ka / kb) > 1e-9 for x in v)
    bR = ka != kb
    vM |= bM; vK |= bK; vR |= bR
    print(f"   {ca:18s} {cb:18s} {str(sorted(map(str, v))):14s} {(ka / kb) ** 0.5:9.4f} {ka / kb:9.4f}  "
          f"{'x' if bM else '.'} {'x' if bK else '.'} {'x' if bR else '.'} {'x' if bR else '.'}")
print("  => some collision rules out diag(T^+ D T) = M_L:", vM, "; = K_L:", vK,
      "; Rayleigh quotients of D through T = M_L or K_L:", vR, " (for every superoperator D)")
for nm in TRANSPORTS:
    bM = bK = bR = False
    for r, s_, kap in COLL[nm]:
        ka, kb = KOF[cls[r]], KOF[cls[s_]]
        bM |= abs(kap - (ka / kb) ** 0.5) > 1e-9
        bK |= abs(kap - ka / kb) > 1e-9
        bR |= ka != kb
    print(f"   {nm}: collisions rule out diag = M_L: {bM}, diag = K_L: {bK}, Rayleigh = M_L or K_L: {bR}")
# self-consistency of M_L under the transport
YR = IMG[REF].reshape(1296, NR)
inR = [r for r in range(NR) if nrm(YR[:, r]) > 1e-12 and nrm(YR[:, r] - B @ (Bh @ YR[:, r])) < 1e-9 * nrm(YR[:, r])]
pb = {}
for r in inR:
    c = Bh @ YR[:, r]
    pb[r] = (c.conj() @ ML @ c).real / (c.conj() @ c).real
cons = Counter((subtype[r], cls[r], round(pb[r], 4), round(ML[r, r].real, 4)) for r in inR)
print("pull-back of the Letter's M_L along the transport, rho(r) = <T r, M_L T r>/<T r, T r>, on the records whose image lies in R:")
for key, n in sorted(cons.items()):
    print(f"   {key[0]:10s} {key[1]:7s} rho = {key[2]:.4f}   M_L = {key[3]:.4f}   x{n}")
print("   records whose image leaves R (no pull-back):", dict(Counter(subtype[r] for r in range(NR) if r not in inR)))

# ================================================================== 4. the join seen by transported records
hdr("4. The join as seen by the transported qubit-only records")
om3 = np.eye(3).reshape(-1) / np.sqrt(3); om2S = np.zeros(9); om2S[[0, 4]] = 1 / np.sqrt(2)
print("restricted join (PS (x) PS)|omega_3> = sqrt(2/3) |omega_2> on S (x) S: defect",
      f"{nrm(np.kron(PS, PS) @ om3 - np.sqrt(2 / 3) * om2S):.1e}")
Z22 = rec.slot_op({0: E22, 2: E22}).real
P3 = oo.P3
qonly = [r for r in range(NR) if rtype[r] in ('qq', 'comm')]


def join_support(Y):
    """Q1-support of the join-reduced operator X_J = 3 Tr_Q2[Y P_3] (Y |omega_3> = (X_J (x) 1)|omega_3>)"""
    X8 = (Y @ P3).reshape(3, 2, 3, 2, 3, 2, 3, 2)
    XJ = 3 * np.trace(X8, axis1=2, axis2=6)                      # (Q1,q1,q2 ; Q1,q1,q2)
    XJ = XJ.reshape(12, 12)
    G = XJ @ XJ.conj().T + XJ.conj().T @ XJ
    return np.trace(G.reshape(3, 4, 3, 4), axis1=1, axis2=3)


def support_class(G3):
    if nrm(G3) < 1e-12:
        return 'none'
    w = np.linalg.eigvalsh((G3 + G3.conj().T) / 2)
    rank = int((w > 1e-9 * w.max()).sum())
    if abs(G3[2, 2]) < 1e-9 * nrm(G3) and nrm(G3[2, :]) < 1e-9 * nrm(G3):
        return f'in S (rank {rank})'
    if nrm(G3 - G3[2, 2] * E22) < 1e-9 * nrm(G3):
        return 'C|2>'
    return f'reaches |2> (rank {rank})'


def letter_support(Y):
    """3x3 PSD whose range is the union of the supports of Y's letters (traceless parts) at the qutrit copies;
    Q2 is identified with Q1 by the join (a letter B at Q2 acts on |omega_3> as B^T at Q1)"""
    G = np.zeros((3, 3), complex)
    for s_ in (0, 2):
        O = oo.COTAU[s_](Y[:, :, None])[:, :, 0]
        OO = (O @ O.conj().T + O.conj().T @ O).reshape(3, 2, 3, 2, 3, 2, 3, 2)
        rows = list('pqrs'); cols = list('pqrs'); cols[s_] = 'z'
        Gs = np.einsum(''.join(rows) + ''.join(cols) + '->' + rows[s_] + 'z', OO)
        G += Gs.T if s_ == 2 else Gs
    return G


def h_rank(G3):
    """rank of the smallest H-invariant projector (0, PS, |2><2| or 1) containing the range of G3"""
    sc = support_class(G3)
    return 0 if sc == 'none' else (2 if sc.startswith('in S') else (1 if sc == 'C|2>' else 3))


print(f"{'transport':12s} {'max|T r Z22|,|Z22 T r|':>23s} {'T r |omega_3> = 0':>18s}  letter support of the 12 qubit-only records")
for nm, Y in IMG.items():
    z = max(max(nrm(Y[:, :, r] @ Z22), nrm(Z22 @ Y[:, :, r])) / nrm(Y[:, :, r]) for r in qonly)
    ann = sum(nrm(join_support(Y[:, :, r])) < 1e-12 * nrm(Y[:, :, r]) ** 2 for r in qonly)
    reach = Counter(str(rtype[r]) for r in qonly if not support_class(letter_support(Y[:, :, r])).startswith('in S'))
    print(f"{nm:12s} {z:23.1e} {ann:18d}  {dict(Counter(support_class(letter_support(Y[:, :, r])) for r in qonly))}"
          + (f"  reaching |2>: {dict(reach)}" if reach else ""))
print("  (with b = 0 the identity left at the other copy becomes PS, whose traceless part (prop. to lambda_8) reaches |2>;"
      " the letters themselves are t(sigma) on S)")
YREF = IMG[REF]
ann_q = [r for r in qonly if nrm(join_support(YREF[:, :, r])) < 1e-12]
print(f"{REF}: qubit-only records whose image annihilates the join, T r |omega_3> = 0:",
      dict(Counter((cls[r], subtype[r]) for r in ann_q)), "(E'-odd single-slot sums and anticommuting pairs, by the transpose identity)")
natQQ = [r for r in range(NR) if rtype[r] == 'QQ']
print("native QQ records for comparison: letter supports under H:", dict(Counter(support_class(letter_support(B3[:, :, r])) for r in natQQ)),
      "; annihilating the join:", sum(nrm(join_support(B3[:, :, r])) < 1e-12 for r in natQQ),
      "(the Letter gives every E'-even QQ record sqrt3 regardless)")
print("smallest invariant projector containing the support S of a transported letter: under H it is PS (rank 2), under SU(3) it is 1 (rank 3)")
C3 = np.eye(3) / np.sqrt(3)
AMP = {'J1': np.trace(np.eye(3) @ C3), 'J2': np.trace(PS @ C3), 'J3': np.trace(np.eye(2) / np.sqrt(2))}
# J4: the path q1 -> Q1 -> join -> Q2 -> q2 in the Letter's joint state (S7, without the unread link)
om_link = np.zeros(6); om_link[[0, 3]] = 1 / np.sqrt(2)              # |omega_2> on (Q_i, q_i), Q-major
v = np.kron(om_link, om_link)                                       # order Q1 q1 Q2 q2
Om = np.sqrt(6) * (P3 @ v)
rho_q = np.einsum('aibj,akbl->ijkl', Om.reshape(3, 2, 3, 2), Om.reshape(3, 2, 3, 2).conj()).reshape(4, 4)
wq, vq = np.linalg.eigh(rho_q)
psi = vq[:, -1].reshape(2, 2)
AMP['J4'] = np.sum(np.linalg.svd(psi, compute_uv=False))
print(f"joint state: ||Omega|| = {nrm(Om):.12f}; its q1q2 marginal is pure (largest eigenvalue {wq[-1]:.12f}), a qubit pair")
print("amplitude summed over matched directions, a = Tr(P C_omega):")
print(f"   J1 = Tr(1_3 . 1_3/sqrt3) = {AMP['J1'].real:.12f} = sqrt3   (rank-3 projector: the join's three directions)")
print(f"   J2 = Tr(PS  . 1_3/sqrt3) = {AMP['J2'].real:.12f} = 2/sqrt3 (minimal H-invariant support: two of the join's strands)")
print(f"   J3 = Tr C of the restricted join renormalized, |omega_2> on S (x) S = {AMP['J3']:.12f} = sqrt2")
print(f"   J4 = Tr C of the q1q2 pair reached through link-join-link in the joint state = {AMP['J4']:.12f} = sqrt2")
JK = {'J1': Fraction(3), 'J2': Fraction(4, 3), 'J3': Fraction(2)}     # squared credits a^2
for j, kv in JK.items():
    assert abs(float(kv) - abs(AMP[j]) ** 2) < 1e-12

# ================================================================== 5. family F3
hdr("5. Family F3: frame-local survival x attribution dressings (exact)")
EQR, l1 = oo.compress(oo.EpQ); EqR, l2 = oo.compress(oo.Epq); SGR, l3 = oo.compress(oo.SIGMA)
EPR = oo.EPR
TPR = {t: oo.compress(f) for t, f in oo.TYPE_OPS.items()}
assert max(l1, l2, l3, *[v[1] for v in TPR.values()]) < 1e-12
Hj = EPR + np.pi * EQR + np.e * EqR + 7.1 * SGR + sum((i + 1) * 11.3 * TPR[t][0] for i, t in enumerate(oo.TYPE_OPS))
wj, Vj = np.linalg.eigh((Hj + Hj.conj().T) / 2)
labs = []
for j in range(NR):
    vv = Vj[:, j]
    rq = lambda X: int(round((vv.conj() @ X @ vv).real))
    t = [tt for tt in oo.TYPE_OPS if rq(TPR[tt][0]) == 1][0]
    labs.append((t, rq(EPR), rq(EQR), rq(EqR), rq(SGR)))
LABC = Counter(labs)
ULAB = sorted(LABC)
dimv = [LABC[l] for l in ULAB]
print("joint eigen-labels (type, E', E'_Q, E'_q, sigma) and dimensions (as in operator_origin.py Section 5):")
for l in ULAB:
    print("   ", l, LABC[l])
PLAB = {l: Vj[:, [j for j in range(NR) if labs[j] == l]] for l in ULAB}
CONT = np.array([[np.trace(PLAB[l].conj().T @ P @ PLAB[l]).real for P in oo.ISO_P] for l in ULAB])   # dims per (label, isotypic)
SURV = {
    '1': lambda l: 1,
    "Pi+(E')": lambda l: int(l[1] == 1), "Pi+(sigma E')": lambda l: int(l[1] * l[4] == 1),
    "Pi+(E'_Q)": lambda l: int(l[2] == 1), "Pi+(sigma E'_Q)": lambda l: int(l[2] * l[4] == 1),
    "Pi+(E'_q)": lambda l: int(l[3] == 1), "Pi+(sigma E'_q)": lambda l: int(l[3] * l[4] == 1),
}
TYPES = ['QQ', 'qq', 'mix', 'comm']
DOMAINS = [frozenset(c) for r_ in range(5) for c in itertools.combinations(TYPES, r_)]
FRAMES = {
    'orig-Q': {'QQ': 'Q', 'qq': 'q', 'mix': 'Q', 'comm': 'q'},     # where the letters sit; Q when both ends
    'orig-q': {'QQ': 'Q', 'qq': 'q', 'mix': 'q', 'comm': 'q'},     # where the letters sit; q when both ends
    'tr':     {'QQ': 'Q', 'qq': 'Q', 'mix': 'Q', 'comm': 'Q'},     # after transport every letter sits at a copy
}
TRONLY = {'qq', 'comm'}                    # records whose qutrit letters are all transported ones


def letter_class(l):
    if l[0] == 'comm':
        return 'comm'
    return 'odd' if l[1] == -1 else ('evsame' if l[0] in ('QQ', 'qq') else 'evmix')


KL_lab = tuple(Fraction(KOF[letter_class(l)]) for l in ULAB)


def factor_vec(sname, dom, frame, J):
    out = []
    for l in ULAB:
        t = l[0]
        if t not in dom or not SURV[sname](l):
            out.append(Fraction(1)); continue
        end = FRAMES[frame][t]
        if end == 'q':
            out.append(Fraction(2))
        else:
            out.append(JK[J] if (frame == 'tr' and t in TRONLY) else Fraction(3))
    return tuple(out)


def pexp(n, p):
    e = 0
    while n % p == 0:
        n //= p; e += 1
    return e, n


def sqrt_coord(K):
    """sqrt(K), K = 2^a 3^b rational, as rational coordinates on (1, sqrt2, sqrt3, sqrt6)"""
    a1, r1 = pexp(K.numerator, 2); b1, r1 = pexp(r1, 3)
    a2, r2 = pexp(K.denominator, 2); b2, r2 = pexp(r2, 3)
    assert r1 == 1 and r2 == 1
    a, b = a1 - a2, b1 - b2
    coef = Fraction(2) ** (a // 2) * Fraction(3) ** (b // 2)
    v = [Fraction(0)] * 4
    v[{(0, 0): 0, (1, 0): 1, (0, 1): 2, (1, 1): 3}[(a % 2, b % 2)]] = coef
    return v


def numerator(kv):
    c = [Fraction(0)] * 4
    for K, n in zip(kv, dimv):
        for i, x in enumerate(sqrt_coord(K)):
            c[i] += n * x
    return tuple(c)


num_float = lambda c: float(c[0]) + float(c[1]) * 2 ** 0.5 + float(c[2]) * 3 ** 0.5 + float(c[3]) * 6 ** 0.5


def num_str(c):
    parts = []
    for x, s in zip(c, ('', 'r2', 'r3', 'r6')):
        if x != 0:
            parts.append(f"{x}{s}")
    return '+'.join(parts).replace('+-', '-') if parts else '0'


MLn = np.sqrt(319.0)
dist = lambda kv: np.sqrt(sum(n * (float(K) ** 0.5 - float(k0) ** 0.5) ** 2 for K, k0, n in zip(kv, KL_lab, dimv))) / MLn
tvec = lambda kv: CONT.T @ np.array([float(K) ** 0.5 for K in kv])
tM = oo.tM


def show_choice(c):
    s, d, fr = c
    return f"{s} on {{{','.join(t for t in TYPES if t in d) or '-'}}} in {fr}"


F3 = {}
for J in JK:
    singles = defaultdict(list)
    for sname in SURV:
        for dom in DOMAINS:
            for fr in FRAMES:
                singles[factor_vec(sname, dom, fr, J)].append((sname, dom, fr))
    dress = defaultdict(list)                  # K label vector -> list of (factor1, factor2) keys
    for f1 in singles:
        for f2 in singles:
            dress[tuple(a * b for a, b in zip(f1, f2))].append((f1, f2))
    F3[J] = (singles, dress)
n_single = len(SURV) * len(DOMAINS) * len(FRAMES)
print(f"per exchange: {len(SURV)} survival operators x {len(DOMAINS)} domains x {len(FRAMES)} frames = {n_single};"
      f" ordered pairs {n_single ** 2}; x {len(JK)} join credits = {n_single ** 2 * len(JK)} ordered choices")
all_kv = set()
for J, (singles, dress) in F3.items():
    all_kv |= set(dress)
    print(f"  {J}: distinct single-exchange factors {len(singles)}, distinct dressings {len(dress)}")
print("  distinct dressings over all three credits:", len(all_kv))
hits = []
for J, (singles, dress) in F3.items():
    if KL_lab in dress:
        for f1, f2 in dress[KL_lab]:
            for c1 in singles[f1]:
                for c2 in singles[f2]:
                    hits.append((J, c1, c2))
nQocc = oo.lin((1, oo.ID), (-1, oo.compose(oo.TAU[0], oo.TAU[2])))          # [N_Q >= 1]
TBr = IMG[REF].reshape(1296, NR)
PB = TBr.conj().T @ nQocc(IMG[REF]).reshape(1296, NR)
print(f"F3 reads the transported frame record by record. The linear pull-back T^+ [N_Q>=1] T of {REF} on R is not diagonal"
      f" in the record basis: off-diagonal part {nrm(PB - np.diag(np.diag(PB))) / nrm(PB):.3f} of the norm (T has rank"
      f" {np.linalg.matrix_rank(TBr, tol=1e-9)} on R)")
noQ = oo.compose(oo.TAU[0], oo.TAU[2])                                        # [N_Q = 0]
TSel = TRANSPORTS[REF][0](noQ(B3)).reshape(1296, NR)
wS, vS = np.linalg.eigh(TSel.conj().T @ TSel)
Psupp = vS[:, wS > 1e-9] @ vS[:, wS > 1e-9].conj().T
PXR = TPR['QQ'][0] + TPR['qq'][0] + TPR['comm'][0]                            # Pi_X = [N_Q N_q = 0] on R
PQR = TPR['QQ'][0] + TPR['mix'][0]                                            # [N_Q >= 1] on R
print(f"selective transport T' = T [N_Q = 0]: rank {int((wS > 1e-9).sum())} on R, support = qq + comm records:"
      f" {nrm(Psupp - TPR['qq'][0] - TPR['comm'][0]) < 1e-9};  Pi_X = Pi_X [N_Q>=1] + supp(T'^+ T') on R:"
      f" {nrm(PXR - PXR @ PQR - Psupp) < 1e-9}  (the selector [N_Q = 0] is the 12 qubit-only records)")
print("ordered choices giving exactly M_L:", len(hits), "; credits used:", dict(Counter(h[0] for h in hits)),
      "; frame pairs used:", dict(Counter(tuple(sorted((h[1][2], h[2][2]))) for h in hits)))
qqe, cmi = ULAB.index(('qq', 1, 1, 1, 1)), ULAB.index(('comm', -1, 1, -1, -1))
tr_credit = all(factor_vec(*(c1 if c1[2] == 'tr' else c2), J)[qqe] == 3 and factor_vec(*(c1 if c1[2] == 'tr' else c2), J)[cmi] == 3
                for J, c1, c2 in hits)
print("  every hit has one exchange in the transported frame, which credits 3 to the E'-even qq and the commutator records,"
      " and the other in a qubit-credited original frame:",
      all(sorted((h[1][2], h[2][2])) == ['orig-q', 'tr'] for h in hits) and tr_credit)
for h in hits[:3]:
    print("   e.g.", h[0], "|", show_choice(h[1]), "|", show_choice(h[2]))
ok_orig = []
for J, (singles, dress) in F3.items():
    so = {f for f, cs in singles.items() if any(c[2] != 'tr' for c in cs)}
    ok_orig.append(any(KL_lab == tuple(a * b for a, b in zip(f1, f2)) for f1 in so for f2 in so))
print("  M_L from two original-frame exchanges (Prop. 5 of operator-origin.md):", any(ok_orig))
st = {f for f, cs in F3['J1'][0].items() if any(c[2] == 'tr' for c in cs)}
print("  M_L from two transported-frame exchanges:", any(KL_lab == tuple(a * b for a, b in zip(f1, f2)) for f1 in st for f2 in st))

# the 18 Letter-shaped members: Letter's survival and domains, frames and credit varied
print("\nLetter-shaped members (X: Pi+(sigma E') on {QQ,qq,comm}; Y: Pi+(sigma E') on {mix,comm}); frames and credit varied:")
print(f"   {'credit':6s} {'X frame':7s} {'Y frame':7s} {'numerator over 31104':28s} {'1/value':>9s} {'dist':>7s}  = M_L")
LS_rows = {}
for J in JK:
    for fx in ('orig-q', 'tr'):
        for fy in ('orig-q', 'orig-Q', 'tr'):
            kv = tuple(a * b for a, b in zip(factor_vec("Pi+(sigma E')", frozenset({'QQ', 'qq', 'comm'}), fx, J),
                                             factor_vec("Pi+(sigma E')", frozenset({'mix', 'comm'}), fy, J)))
            c = numerator(kv)
            LS_rows[(J, fx, fy)] = (kv, c)
            print(f"   {J:6s} {fx:7s} {fy:7s} {num_str(c):28s} {N / num_float(c):9.3f} {dist(kv):7.4f}  {kv == KL_lab}")
print("   ('orig-q' for X equals 'orig-Q' for X: X's domain has no mixed records.)")
S5 = {
    'computed (Letter)': (79, 48, 42, 3), 'two ends swapped': (79, 42, 48, 3), 'both at qubit end': (85, 90, 0, 0),
    'both at qutrit end': (88, 0, 90, 0), 'X where letters sit': (85, 54, 36, 0), 'qq moved by Y': (79, 48, 36, 9),
    'copies separate': (79, 0, 42, 51), 'no own sign': (82, 48, 42, 0), 'comm survive X only': (79, 48, 45, 0),
    'no weights': (172, 0, 0, 0),
}
num_of = {}
for J, (singles, dress) in F3.items():
    for kv in dress:
        num_of.setdefault(numerator(kv), set()).add(J)
print("Supplement Table S5 weight-law rows inside F3 (credits under which the numerator occurs):")
for nm_, val in S5.items():
    key = tuple(Fraction(x) for x in val)
    print(f"   {nm_:22s} {num_str(key):22s} 1/value {N / num_float(key):9.3f}  {sorted(num_of.get(key, []))}")
print("closest members to M_L that are not M_L (HS distance ||M - M_L|| / ||M_L||), per credit:")
CLOSEST = {}
for J, (singles, dress) in F3.items():
    lst = sorted(((dist(kv), kv) for kv in dress if kv != KL_lab), key=lambda x: x[0])[:4]
    CLOSEST[J] = lst
    for dd, kv in lst:
        f1, f2 = dress[kv][0]
        print(f"   {J}: dist {dd:.4f}  {num_str(numerator(kv)):30s} 1/value {N / num_float(numerator(kv)):9.3f}"
              f"   [{show_choice(singles[f1][0])} | {show_choice(singles[f2][0])}]")
print("per-isotypic traces proportional to the Letter's point (projective residual < 1e-9):")
for J, (singles, dress) in F3.items():
    tp = [kv for kv in dress if oo.proj_resid(tvec(kv).astype(complex), tM) < 1e-9]
    print(f"   {J}: {len(tp)} dressings; equal to M_L: {sum(kv == KL_lab for kv in tp)};"
          f" others: {[num_str(numerator(kv)) for kv in tp if kv != KL_lab][:4]}")
vals = {}
for J, (singles, dress) in F3.items():
    for kv in dress:
        c = numerator(kv)
        vals.setdefault(c, set()).add(J)
near = [c for c in vals if abs(N / num_float(c) / AINV - 1) < 1e-3]
print(f"distinct values in F3: {len(vals)}; within 1e-3 relative of 1/alpha: {len(near)}")
for c in sorted(near, key=lambda c: abs(N / num_float(c) - AINV)):
    print(f"   {num_str(c):30s} 1/value {N / num_float(c):10.5f}  credits {sorted(vals[c])}")
rng_vals = [N / num_float(c) for c in vals]
print(f"1/value range over F3: [{min(rng_vals):.3f}, {max(rng_vals):.3f}]")

# --- the H-orbit credit applied to every record (not only to the transported ones)
print("\nH-orbit credit for every record: X (Pi+(sigma E') on QQ, qq, comm) in the transported frame with squared credit"
      " (rank of the smallest H-invariant projector containing the record's letter supports)^2 / 3; Y (Pi+(sigma E') on mix,"
      " comm) in the original frame, qubit end:")
KX = np.ones(NR); KY = np.ones(NR)
supp_kind = Counter()
for r in range(NR):
    survives = (epar[r] == 1 and rtype[r] in ('QQ', 'qq')) or rtype[r] == 'comm'
    if survives:
        rank_min = h_rank(letter_support(YREF[:, :, r]))
        KX[r] = rank_min ** 2 / 3
        supp_kind[(str(rtype[r]), rank_min)] += 1
    if (epar[r] == 1 and rtype[r] == 'mix') or rtype[r] == 'comm':
        KY[r] = 2
MH = np.diag(np.sqrt(KX * KY)).astype(complex)
print("   X-survivors by (type, rank of minimal H-invariant support):", dict(sorted(supp_kind.items())))
coH = [Fraction(0)] * 4
for r in range(NR):
    for i, x in enumerate(sqrt_coord(Fraction(KX[r] * KY[r]).limit_denominator(1000))):
        coH[i] += x
print(f"   numerator {num_str(tuple(coH))};  1/value {N / np.trace(MH).real:.3f};  dist to M_L {nrm(MH - ML) / nrm(ML):.4f};"
      f"  G_Omega x <E'> twirl defect {nrm(MH - oo.twirl(MH)) / nrm(MH):.3f};"
      f"  per-isotypic residual {oo.proj_resid(oo.iso_traces(MH), tM):.4f}")

# ================================================================== 6. family F3a: words with transport letters
hdr("6. Family F3a: words of length <= 3 in 24 covariant primitives and 9 transport letters")
TL = []
for b in (0, 1):
    for (o1, o2) in (('R', 'L'), ('S', 'S')):
        f, fa = TRANSPORTS[f"T_b{b}({o1},{o2})"]
        TL.append((f"T_b{b}({o1},{o2})", f, 'tr')); TL.append((f"T_b{b}({o1},{o2})^+", fa, 'tr'))
CS_op = rec.slot_op({0: PS, 2: PS}).real
TL.append(("C_S", oo.Ad_b(CS_op), 'tr'))
ALPH = [(n, f, g) for (n, f, g) in oo.LETTERS if g == 'cov'] + TL
NA = len(ALPH)
TRIDX = {i for i, (_, _, g) in enumerate(ALPH) if g == 'tr'}
print(f"alphabet: {NA - len(TL)} covariant letters of operator_origin.py and {len(TL)} transport letters:", [n for n, _, _ in TL])
ID1296 = oo.ID1296.astype(complex)
for nm_, f, _ in TL:
    Mx = f(ID1296).reshape(1296, 1296)
    print(f"   {nm_:16s} denominator {oo.min_denominator(Mx)}; Q|q Schmidt rank {oo.schmidt_rank_Qq(Mx)};"
          f" G_Omega-cov defect {nrm(f(oo.gW(oo.Atest)) - oo.gW(f(oo.Atest))) / nrm(f(oo.Atest)):.2f};"
          f" H-equiv defect {nrm(f(gH(oo.Atest)) - gH(f(oo.Atest))) / nrm(f(oo.Atest)):.1e};"
          f" E'-cov defect {nrm(f(oo.Ep(oo.Atest)) - oo.Ep(f(oo.Atest))) / nrm(f(oo.Atest)):.1e}")
print("   (all transport letters are Q-rational, so Lemma 3 of operator-origin.md applies to every word of F3a;"
      " the eight transports are not Q|q products, so Theorem 4 there does not apply to words containing them)")
del ID1296


def image(word):
    Y = B3
    for li in reversed(word):
        Y = ALPH[li][1](Y)
    return Y


def hkey(Z):
    k_ = oo.canon_key(Z)
    return None if k_ is None else hashlib.sha1(k_).digest()


def w_extend(w2):
    Y2 = image(w2)
    return [(li, hkey(ALPH[li][1](Y2))) for li in range(NA)]


def w_analyse(w):
    Z = image(w)
    Yf = Z.reshape(1296, NR)
    C = Bh @ Yf
    if nrm(C) < 1e-10 * nrm(Yf):
        return w, None
    a = oo.analyse(C)
    a['leak'] = nrm(Yf - B @ C) / nrm(Yf)
    a['trace_ratio'] = np.trace(C)
    a['odd_mean'] = np.diag(C)[oo.CIDX['odd']].mean()
    return w, a


def rat_pattern(v):
    """copied from operator_origin.py (defined there after the point where its source is cut)"""
    j = int(np.argmax(np.abs(v) > 1e-9 * np.abs(v).max()))
    u = v / v[j]
    if np.abs(u.imag).max() > 1e-9:
        return None
    fr = tuple(Fraction(float(x)).limit_denominator(10000) for x in u.real)
    if max(abs(float(f) - x) for f, x in zip(fr, u.real)) > 1e-9:
        return ('irrational',) + tuple(np.round(u.real, 6))
    return fr


def fmt_pat(p):
    return '(' + ', '.join(str(x) for x in p) + ')' if p is not None else 'None'


t6 = time.time()
seen = {}
order = []                                   # first-occurrence words, in enumeration order
nformed = Counter()


def visit(w, key):
    if key is None or key in seen:
        return False
    seen[key] = w; order.append(w)
    return True


visit((), hkey(B3)); nformed[0] = 1
lev1 = []
for li in range(NA):
    nformed[1] += 1
    if visit((li,), hkey(ALPH[li][1](B3))):
        lev1.append((li,))
lev2 = []
for w1 in lev1:
    Y1 = image(w1)
    for li in range(NA):
        nformed[2] += 1
        if visit((li,) + w1, hkey(ALPH[li][1](Y1))):
            lev2.append((li,) + w1)
ctx = multiprocessing.get_context('fork')
with ctx.Pool(NPROC) as pool:
    ext = pool.map(w_extend, lev2, chunksize=4)
    for w2, lst in zip(lev2, ext):
        for li, key in lst:
            nformed[3] += 1
            visit((li,) + w2, key)
    new_tr = [w for w in order if any(i in TRIDX for i in w)]
    res = pool.map(w_analyse, new_tr, chunksize=8)
print(f"[F3a] enumeration and analysis {time.time() - t6:.0f} s ({NPROC} processes)", file=sys.stderr)
nnew = Counter(len(w) for w in order)
for L in range(4):
    print(f"length {L}: {nformed[L]:6d} words formed, {nnew[L]:6d} new distinct nonzero images of R (up to scale)")
total_words = sum(NA ** L for L in range(4))
cov_words = sum((NA - len(TL)) ** L for L in range(4))
print(f"words in F3a: all words of length <= 3 over {NA} letters ({total_words}), of which {total_words - cov_words} contain a"
      f" transport letter; distinct images: {len(order)}, of which {len(new_tr)} first reached by a word with a transport"
      " letter (the others are images of F1 words)")
good = [(w, a) for w, a in res if a is not None]
print("analysed (distinct, nonzero on R, with a transport letter):", len(good), "; zero on R:", sum(a is None for _, a in res))
print("  compression G_Omega-invariant:", sum(a['inv'] < TOL for _, a in good), "| also E'-invariant:",
      sum((a['inv'] < TOL) and (a['epinv'] < TOL) for _, a in good), "| maps R into R:", sum(a['leak'] < TOL for _, a in good),
      "| diagonal in the record basis:", sum(a['offdiag'] < TOL for _, a in good))
wname = lambda w: ' . '.join(ALPH[i][0] for i in w) if w else '1 (identity)'
sqk = np.sqrt(np.array([1, 2, 3, 6.0]))
allpats = []
for tag, cskey, vkey in (('raw C', 'cs', 'csv'), ("twirled over G x <E'>", 'cs_tw', 'csv_tw'), ('|C|', 'abs_cs', 'abs_csv')):
    pats = Counter()
    for w, a in good:
        if a[cskey] < TOL and np.abs(a[vkey]).max() > 1e-9:
            pats[rat_pattern(a[vkey])] += 1
    allpats += list(pats.elements())
    print(f"[{tag}] class-scalar: {sum(pats.values())} words, patterns: "
          f"{[(fmt_pat(p), n) for p, n in pats.most_common(10)]}{' ...' if len(pats) > 10 else ''}")
    print("   prop. to (1, sqrt2, sqrt3, sqrt6):",
          sum(n for p, n in pats.items() if p is not None and p[0] == 'irrational' and np.allclose(np.array(p[1:], float), sqk, atol=1e-6)),
          "  prop. to (1, 2, 3, 6):", sum(n for p, n in pats.items() if p is not None and p[0] != 'irrational' and tuple(p) == (1, 2, 3, 6)))
sqs = [(w, a) for w, a in good if a['psd'] and a['sq_cs'] < TOL and np.abs(a['sq_csv']).max() > 1e-9]
allpats += [rat_pattern(a['sq_csv']) for _, a in sqs]
print(f"[sqrt C, C Hermitian PSD] {sum(a['psd'] for _, a in good)} PSD words; class-scalar sqrt: {len(sqs)};"
      f" patterns {Counter(fmt_pat(rat_pattern(a['sq_csv'])) for _, a in sqs).most_common(6)}")
dconst = [(w, a) for w, a in good if a['dspread'] < TOL and np.abs(a['dmeans']).max() > 1e-9]
print(f"(b') diagonal constant on each class: {len(dconst)} words; patterns "
      f"{Counter(fmt_pat(rat_pattern(a['dmeans'])) for _, a in dconst).most_common(8)}")
allpats += [rat_pattern(a['dmeans']) for _, a in dconst]
ok01 = all(p is not None and p[0] != 'irrational' and all(x in (0, 1, -1) for x in p) for p in allpats)
print("class-scalar or class-constant patterns found in F3a (all tests):", len(allpats),
      "" if ok01 else [fmt_pat(p) for p in set(allpats) if p is None or p[0] == 'irrational' or not all(x in (0, 1, -1) for x in p)][:8])
A_ID = oo.analyse(np.eye(NR).astype(complex))


def value_odd_norm(a):
    """1/value of the word with the scale fixed so that the 79 E'-odd records weigh 1 on average"""
    if abs(a['odd_mean']) < 1e-12:
        return float('nan')
    return N / (a['trace_ratio'] / a['odd_mean']).real


def show(title, key, n=4, dkey='dmeans', filt=lambda a: True):
    print(f"\n{title}   [identity: {A_ID[key]:.4f}]")
    lst = sorted([(a[key], w, a) for w, a in good if key in a and filt(a)], key=lambda x: x[0])[:n]
    for val, w, a in lst:
        dm = a[dkey] / a[dkey][0] if abs(a[dkey][0]) > 1e-12 else a[dkey]
        print(f"   resid {val:.4f}  [{wname(w)}]  G-inv defect {a['inv']:.0e}  class means (odd, evmix, evsame, comm)"
              f" rel. to odd {np.round(dm.real, 4)};  1/value (odd records weigh 1) {value_odd_norm(a):.3f}")
    return lst


show("closest misses, diagonal vs c*sqrt(k_r)", 'rM')
show("closest misses, diagonal vs c*k_r (K-route)", 'rK')
show("closest misses, per-isotypic traces vs the Letter's point", 'rtM')
show("closest misses, per-isotypic traces of |C|", 'abs_rtM', dkey='abs_dmeans')
show("closest misses, per-isotypic traces of sqrt C (C Hermitian PSD)", 'sq_rtM', dkey='sq_dmeans')
print("\nminimum residual over F3a / residual of the identity (0 = exact reproduction):")
for key in ('rM', 'rK', 'rtM', 'abs_rM', 'abs_rtM', 'sq_rM', 'sq_rtM'):
    print(f"   {key:8s} {min(a[key] for _, a in good if key in a):.4f} / {A_ID[key]:.4f}")
print("   F1 (operator_origin.out) minima for comparison: rM 0.2471, rK 0.4572, rtM 0.1008, abs_rM 0.2471, abs_rtM 0.0912,"
      " sq_rM 0.2471, sq_rtM 0.1007")

# ================================================================== 7. summary
hdr("SUMMARY")
print("1. The transport of a qubit letter through its read link is unique (t(B) = iota B^T iota^dag), lives on the 2-dim")
print("   support S, is H-equivariant (H = SU(2) x U(1), dim 4) and has no nonzero G_Omega-equivariant version.")
print("2. Uniform transport (every record) is many-to-one on R: odd cross-mixed records collide with odd QQ records (kappa = 2/3,")
print("   M_L and K_L need 1), even cross-mixed with even QQ, and commutator with odd same-link mixed records (Rayleigh quotients")
print("   forced equal, the Letter has sqrt2 vs sqrt3 and sqrt6 vs 1). So no D in the transported frame, as T^+ D T or as")
print("   Rayleigh quotients, gives M_L or K_L, for any of the six transports.")
nh = len(hits)
print(f"3. F3 ({n_single ** 2 * len(JK)} ordered choices): M_L is reached by {nh} choices, all with the credit J1 = sqrt3 and with X in the")
print("   transported frame and Y in the original frame (qubit end). With the H-invariant credits the Letter-shaped member gives")
for J in ('J2', 'J3'):
    kv, c = LS_rows[(J, 'tr', 'orig-q')]
    print(f"   {J}: {num_str(c)}  1/value {N / num_float(c):.3f}")
print(f"4. F3a ({len(good)} distinct compressions with a transport letter): no class-scalar or class-constant word at all, so none")
print("   reproduces sqrt(k_r) or k_r as operator or diagonal; per-isotypic traces: best residual"
      f" {min(a['rtM'] for _, a in good):.4f} (identity 0.2070, F1 best 0.1008).")
print(f"runtime {time.time() - T_START:.0f} s", file=sys.stderr)
