"""GU Phase 2(b): branching of Spin(14) / Spin(7,7) spinor and Rarita-Schwinger modules under
Spin(4) x Spin(10)  (complexification of Spin(1,3) x Spin(6,4)),
then Spin(10) -> Spin(6) x Spin(4) = SU(4) x SU(2)_L x SU(2)_R  (Pati-Salam)
then -> SU(3) x SU(2)_L x U(1)_Y  with  Y = (B-L)/2 + T3R.

Exact weight-multiset (character) computations; integer arithmetic on doubled weights.
Claims tested: GU-025 ... GU-027, GU-029, GU-058 ... GU-065 (claims/geometric-unity.md).
Write-up: formalization/geometric-unity/B-spinor-branching.md

Complex representations only: reality (Majorana) structure is handled separately in
clifford77.py and in the write-up.
"""
from collections import Counter
from fractions import Fraction
import itertools


# ------------------------------------------------------------ weight systems (doubled weights)
def spinor(m, chir):
    """S^+ (chir=+1): even number of minus signs; S^- (chir=-1): odd. Weights doubled (entries +-1)."""
    out = Counter()
    for s in itertools.product([1, -1], repeat=m):
        neg = s.count(-1)
        if (neg % 2 == 0) == (chir == 1):
            out[s] += 1
    return out


def vector(m):
    out = Counter()
    for i in range(m):
        for sgn in (2, -2):
            w = [0] * m
            w[i] = sgn
            out[tuple(w)] += 1
    return out


def trivial(m):
    return Counter({tuple([0] * m): 1})


def tensor(A, B):
    out = Counter()
    for a, ma in A.items():
        for b, mb in B.items():
            out[tuple(x + y for x, y in zip(a, b))] += ma * mb
    return out


def minus(A, B):
    out = Counter(A)
    for k, v in B.items():
        out[k] -= v
        assert out[k] >= 0, "negative multiplicity: B not contained in A"
        if out[k] == 0:
            del out[k]
    return out


def rarita(m, chir):
    """R^{+-} := V (x) S^{+-}  minus  S^{-+}  (irreducible, highest weight e1 + spinor hw)."""
    return minus(tensor(vector(m), spinor(m, chir)), spinor(m, -chir))


def split(C, k):
    """Restrict a weight multiset of D_m to D_k x D_{m-k} (first k coordinates | rest)."""
    out = Counter()
    for w, mult in C.items():
        out[(w[:k], w[k:])] += mult
    return out


def dim(C):
    return sum(C.values())


# ------------------------------------------------------------ irreducible libraries + peeling
def library_D(m):
    """Named irreps of D_m (m>=2) used here, keyed by doubled highest weight."""
    lib = {}
    lib[tuple([0] * m)] = ('1', trivial(m))
    lib[tuple([2] + [0] * (m - 1))] = ('V', vector(m))
    for chir, nm in [(1, 'S+'), (-1, 'S-')]:
        S = spinor(m, chir)
        lib[max(S, key=lambda w: fD(w))] = (nm, S)
        R = rarita(m, chir)
        lib[max(R, key=lambda w: fD(w))] = ('R' + nm[1], R)
    return lib


def fD(w):
    """Linear functional strictly positive on positive roots e_i +- e_j (i<j) of D_m."""
    m = len(w)
    return sum((m - i + 1) * x for i, x in enumerate(w))


def decompose_product(C, libs):
    """Decompose a character of D_{m1} x D_{m2} (keys (w1,w2)) into products of library irreps."""
    C = Counter(C)
    result = Counter()
    lib1, lib2 = libs
    while C:
        hw = max(C, key=lambda k: (fD(k[0]) + fD(k[1]), fD(k[0])))
        n1, ch1 = lib1[hw[0]]
        n2, ch2 = lib2[hw[1]]
        mult = C[hw]
        prod = Counter({(a, b): ma * mb for a, ma in ch1.items() for b, mb in ch2.items()})
        for k, v in prod.items():
            C[k] -= v * mult
            assert C[k] >= 0
            if C[k] == 0:
                del C[k]
        result[(n1, n2)] += mult
    return result


# ------------------------------------------------------------ SM quantum numbers
def su3_weights(a, b):
    """Weights of SU(3) irrep with Dynkin labels (a,b) as Dynkin-coordinate pairs, via SSYT."""
    lam1, lam2 = a + b, b
    out = Counter()
    letters = (1, 2, 3)
    for row1 in itertools.combinations_with_replacement(letters, lam1):
        for row2 in itertools.combinations_with_replacement(letters, lam2):
            if all(row2[i] > row1[i] for i in range(lam2)):
                c = Counter(row1 + row2)
                out[(c[1] - c[2], c[2] - c[3])] += 1
    return out


SU3_NAMES = {(0, 0): '1', (1, 0): '3', (0, 1): '3b', (2, 0): '6', (0, 2): '6b', (1, 1): '8',
             (3, 0): '10', (0, 3): '10b', (2, 1): '15', (1, 2): '15b'}


def sm_charges(w10):
    """Doubled D5 weight (w1..w5) -> (SU3 Dynkin pair, 2*T3L, 6*Y, 3*Q).
    SU(4)=Spin(6) on w1..w3, Spin(4)=SU(2)_L x SU(2)_R on w4,w5 with
    T3L=(w4+w5)/2, T3R=(w4-w5)/2, B-L = -(2/3)(w1+w2+w3), Y=(B-L)/2+T3R (actual weights)."""
    w = [Fraction(x, 2) for x in w10]
    dyn = (int(w[0] - w[1]), int(w[1] - w[2]))
    T3L = (w[3] + w[4]) / 2
    T3R = (w[3] - w[4]) / 2
    BL = -Fraction(2, 3) * (w[0] + w[1] + w[2])
    Y = BL / 2 + T3R
    Q = T3L + Y
    return dyn, int(2 * T3L), int(6 * Y), Q


def sm_decompose(C10):
    """Decompose a D5 weight multiset into SU(3) x SU(2)_L x U(1)_Y irreps.
    Returns Counter of (su3 name, su2 dim, n=6Y)."""
    pts = Counter()
    for w, m in C10.items():
        dyn, t3, n6, _ = sm_charges(w)
        pts[(dyn, t3, n6)] += m
    res = Counter()
    while pts:
        # highest weight: maximise a functional positive on the simple roots of SU(3)xSU(2)
        # (Dynkin coords: alpha1=(2,-1), alpha2=(-1,2); f=(1,1) is positive on both)
        hw = max(pts, key=lambda k: (k[2], 10 * (k[0][0] + k[0][1]) + k[1]))
        # restrict to same Y and take the highest there
        same = [k for k in pts if k[2] == hw[2]]
        hw = max(same, key=lambda k: 10 * (k[0][0] + k[0][1]) + k[1])
        (a, b), t3, n6 = hw
        assert a >= 0 and b >= 0 and t3 >= 0
        mult = pts[hw]
        for (wt, mw) in su3_weights(a, b).items():
            for t in range(-t3, t3 + 1, 2):
                key = (wt, t, n6)
                pts[key] -= mw * mult
                assert pts[key] >= 0, (hw, key)
                if pts[key] == 0:
                    del pts[key]
        res[(SU3_NAMES.get((a, b), str((a, b))), t3 + 1, n6)] += mult
    return res


def ps_decompose(C10):
    """Spin(10) -> Spin(6) x Spin(4) with SU(4)/SU(2)xSU(2) names."""
    names6 = {'1': '1', 'S+': '4', 'S-': '4b', 'V': '6', 'R+': '20', 'R-': '20b'}
    names4 = {'1': '(1,1)', 'S+': '(2,1)', 'S-': '(1,2)', 'V': '(2,2)', 'R+': '(3,2)', 'R-': '(2,3)'}
    d = decompose_product(split(C10, 3), (library_D(3), library_D(2)))
    return Counter({(names6[a], names4[b]): m for (a, b), m in d.items()})


def fmt_sm(res):
    rows = sorted(res.items(), key=lambda kv: (-kv[0][2], kv[0][0], kv[0][1]))
    return ", ".join(f"[{s}x{d}]^{n:+d}" + (f" x{m}" if m > 1 else "") for (s, d, n), m in rows)


# ------------------------------------------------------------ main
def main():
    print("=== B0. Dimensions ===")
    for m in (2, 5, 7):
        print(f"D{m}: dim S+ = {dim(spinor(m, 1))}, dim S- = {dim(spinor(m, -1))}, "
              f"dim V = {dim(vector(m))}, dim R+ = {dim(rarita(m, 1))}, dim R- = {dim(rarita(m, -1))}")
    # irreducibility sanity: R^+ highest weight has multiplicity 1 and V(x)S^+ = S^- + R^+
    assert tensor(vector(7), spinor(7, 1)) == spinor(7, -1) + rarita(7, 1)

    lib4, lib10 = library_D(2), library_D(5)
    nm4 = {'1': '1', 'S+': '2+=(2,1)', 'S-': '2-=(1,2)', 'V': '4=(2,2)', 'R+': '6+=(3,2)', 'R-': '6-=(2,3)'}
    nm10 = {'1': '1', 'S+': '16+', 'S-': '16-', 'V': '10', 'R+': '144+', 'R-': '144-'}

    def show(label, C14):
        d = decompose_product(split(C14, 2), (lib4, lib10))
        parts = " + ".join(f"{m}*({nm4[a]} x {nm10[b]})" if m > 1 else f"({nm4[a]} x {nm10[b]})"
                           for (a, b), m in sorted(d.items()))
        print(f"{label} [dim {dim(C14)}] -> {parts}")
        return d

    print("\n=== B1. Spin(14) -> Spin(4) x Spin(10) (complexified Spin(7,7) -> Spin(1,3) x Spin(6,4)) ===")
    dSp = show("S+(14)", spinor(7, 1))
    dSm = show("S-(14)", spinor(7, -1))
    show("V(14)", vector(7))
    dRp = show("R+(14)", rarita(7, 1))
    dRm = show("R-(14)", rarita(7, -1))

    print("\n=== B2. Field content nu in Omega^0(S^a), zeta in Omega^1(S^b): count of spin-1/2 x 16 pieces ===")
    for a, b in [(1, -1), (1, 1)]:
        nu = spinor(7, a)
        zeta = tensor(vector(7), spinor(7, b))        # T*Y (x) S^b  (T* = V as Spin rep)
        zdec = show(f"  zeta = V(x)S{'+' if b == 1 else '-'}", zeta)
        tot = decompose_product(split(nu, 2), (lib4, lib10)) + zdec
        F_same = tot[('S+', 'S+')] + tot[('S-', 'S-')]
        F_mirror = tot[('S+', 'S-')] + tot[('S-', 'S+')]
        print(f"  nu in S{'+' if a == 1 else '-'}, zeta in V(x)S{'+' if b == 1 else '-'}: "
              f"#(2+x16+ or 2-x16-) = {F_same}, #(2+x16- or 2-x16+) = {F_mirror}")

    print("\n=== B3. Spin(10) -> Pati-Salam SU(4) x SU(2)_L x SU(2)_R ===")
    for nm, C in [('16+', spinor(5, 1)), ('16-', spinor(5, -1)), ('10', vector(5)),
                  ('144+', rarita(5, 1)), ('144-', rarita(5, -1))]:
        d = ps_decompose(C)
        print(f"{nm:>5} -> " + " + ".join(f"({a},{b[1:-1]})" + (f"x{m}" if m > 1 else "")
                                         for (a, b), m in sorted(d.items())))

    print("\n=== B4. Spin(10) -> SU(3) x SU(2)_L x U(1)_Y, entries [SU3 x SU2]^{n = 6Y} ===")
    sm = {}
    for nm, C in [('16+', spinor(5, 1)), ('16-', spinor(5, -1)), ('10', vector(5)),
                  ('144+', rarita(5, 1)), ('144-', rarita(5, -1))]:
        sm[nm] = sm_decompose(C)
        assert sum(int(SU3dim(s)) * d * m for (s, d, n), m in sm[nm].items()) == dim(C)
        print(f"{nm:>5}: {fmt_sm(sm[nm])}")

    print("\n=== B5. Comparison with the draft's tables (Sec. 11.3, pp. 52-53) ===")
    gen_src = Counter({('3', 2, 1): 1, ('3b', 1, 2): 1, ('3b', 1, -4): 1, ('1', 2, -3): 1,
                       ('1', 1, 6): 1, ('1', 1, 0): 1})
    print("source 'Left' generation table == SM content of 16+ :", gen_src == sm['16+'])
    q32_src = Counter({('3b', 2, -1): 1, ('3', 1, -2): 1, ('3', 1, 4): 1, ('1', 2, 3): 1,
                       ('1', 1, -6): 1, ('1', 1, 0): 1})
    print("source spin-3/2 'cousin' table == SM content of 16- :", q32_src == sm['16-'])
    # source 144 table (p.53), transcribed from a 220-dpi render (bars legible), entries
    # (SU3, SU2 dim, n=6Y): total multiplicity (rows 20 and 21 occur twice in the table and are summed)
    src144 = Counter({('8', 2, -3): 1, ('8', 1, 0): 1, ('8', 1, 6): 1, ('6b', 2, 1): 1, ('6', 1, 2): 1,
                      ('6', 1, -4): 1, ('3b', 3, 2): 1, ('3b', 3, -4): 1, ('3', 2, 7): 1, ('3', 2, -5): 1,
                      ('3', 2, 1): 3, ('3b', 2, 5): 1, ('3', 1, -2): 1, ('3', 1, -8): 1, ('3b', 1, 2): 2,
                      ('3b', 1, -4): 2, ('1', 3, 6): 1, ('1', 2, -9): 1, ('1', 1, 0): 1, ('1', 3, 0): 1,
                      ('1', 2, -3): 2, ('1', 1, 6): 1, ('1', 2, 3): 1})
    print("source 144 table total dimension:", sum(SU3dim(a) * b * m for (a, b, n), m in src144.items()))
    for nm in ('144+', '144-'):
        same = sm[nm] == src144
        print(f"source 144 table == SM content of {nm}: {same}")
        if not same:
            print("   in source not ours:", dict(src144 - sm[nm]))
            print("   in ours not source:", dict(sm[nm] - src144))
    # electric charges Q = T3L + Y of each 144- block, compare with the source's charge column
    print("electric charges in 144- blocks (Q = T3 + n/6):")
    for (s, d, n), m in sorted(sm['144-'].items(), key=lambda kv: (kv[0][0], -kv[0][1], kv[0][2])):
        qs = [Fraction(t, 2) + Fraction(n, 6) for t in range(d - 1, -d, -2)]
        print(f"   [{s}x{d}]^{n:+d} x{m}: Q = " + ", ".join(str(q) for q in qs))

    print("\n=== B6. Does 144 contain a copy of the SM content of 16+ ? (multiset containment) ===")
    for nm in ('144+', '144-'):
        k = min(sm[nm][x] // c for x, c in sm['16+'].items())
        k2 = min(sm[nm][x] // c for x, c in sm['16-'].items())
        print(f"{nm}: max copies of SM(16+) contained = {k}; of SM(16-) = {k2}")

    print("\n=== B8. Majorana-Weyl reading: left-handed spin-1/2 content of nu (+) zeta, nu in S+, zeta in T*(x)S- ===")
    # For a real field whose complexification is (2,1)(x)X + (1,2)(x)conj(X), the independent
    # left-handed Weyl content is X.  Read off the (2+ x .) spin-1/2 pieces of nu and zeta.
    nu = decompose_product(split(spinor(7, 1), 2), (lib4, lib10))
    ze = decompose_product(split(tensor(vector(7), spinor(7, -1)), 2), (lib4, lib10))
    lh = Counter()
    for d in (nu, ze):
        for (a, b), m in d.items():
            if a == 'S+':
                lh[b] += m
    print("left-handed spin-1/2 Spin(10) content:", dict(lh))
    tot = Counter()
    rep = {'S+': '16+', 'S-': '16-', 'R+': '144+', 'R-': '144-'}
    for b, m in lh.items():
        for k, v in sm[rep[b]].items():
            tot[k] += v * m
    ncopies = min(tot[x] // c for x, c in sm['16+'].items())
    print(f"max number of copies of the SM content of one generation (16+) contained: {ncopies}")
    print(f"remaining spin-1/2 states after removing them: {sum(SU3dim(s)*d*m for (s,d,n),m in tot.items()) - 16*ncopies}")
    lh32 = Counter()
    for (a, b), m in ze.items():
        if a == 'R+':
            lh32[b] += m
    print("left-handed spin-3/2 [(3,2) x .] Spin(10) content of zeta:", dict(lh32))

    print("\n=== B9. Reality readings: left-handed spin-1/2 Spin(10) content, net chirality #16+ - #16- ===")
    # complex reading: LH content of (2,1)(x)X is X and of (1,2)(x)X is conj(X); conj swaps +/- for D5.
    conj = {'S+': 'S-', 'S-': 'S+', 'R+': 'R-', 'R-': 'R+', 'V': 'V', '1': '1'}

    def lh_content(C14, real):
        d = decompose_product(split(C14, 2), (lib4, lib10))
        out = Counter()
        for (a, b), m in d.items():
            if a == 'S+':
                out[b] += m
            elif a == 'S-' and not real:
                out[conj[b]] += m
        return out
    readings = [('Majorana-Weyl (real)', True, [1]), ('complex Weyl', False, [1]), ('Dirac (S+ + S-)', False, [1, -1])]
    for label, real, chis in readings:
        for zchis in ([-1], [1]) if len(chis) == 1 else ([1, -1],):
            C = Counter()
            for c in chis:
                C += spinor(7, c)
            Z = Counter()
            for c in zchis:
                Z += tensor(vector(7), spinor(7, c))
            lh = lh_content(C, real) + lh_content(Z, real)
            net = lh['S+'] - lh['S-'] + lh['R+'] * 0
            zlab = '+'.join('S+' if c == 1 else 'S-' for c in zchis)
            print(f"  {label:22s} zeta in T*(x)({zlab}): LH spin-1/2 = {dict(sorted(lh.items()))}; "
                  f"#16+ - #16- = {lh['S+'] - lh['S-']}; #144- - #144+ = {lh['R-'] - lh['R+']}")

    print("\n=== B7. Source formula (4.2): internal spinor dimension for normal bundle of rank n(n+1)/2 ===")
    for n in range(2, 9):
        r = n * (n + 1) // 2
        if r % 2 == 0:
            print(f"n={n}: rank N = {r}, Dirac 2^(r/2) = {2 ** (r // 2)}, Weyl = {2 ** (r // 2 - 1)}; "
                  f"total dim Y = {n * (n + 3) // 2}")
        else:
            print(f"n={n}: rank N = {r} odd -> no Weyl splitting")


def SU3dim(name):
    return {'1': 1, '3': 3, '3b': 3, '6': 6, '6b': 6, '8': 8, '10': 10, '10b': 10, '15': 15, '15b': 15}[name]


if __name__ == '__main__':
    main()
