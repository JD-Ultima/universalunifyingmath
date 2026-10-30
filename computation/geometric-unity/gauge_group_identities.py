"""GU Phase 2(c)-I: exact symbolic check of the inhomogeneous gauge group G = H x| N (draft Secs. 5-7).

Model (local, structure-group independent): base R^2 with coordinates (x, y); structure group GL(2)
realised by matrix-valued polynomial functions; an ad-valued 1-form is a pair (A_x, A_y) of 2x2
matrices; a connection is represented by its potential relative to the trivialisation.
Conventions (reconstruction R-C1, formalization/geometric-unity/C-gauge-shiab.md):
   A . h        = h^-1 A h + h^-1 dh                        (right gauge action)
   d_{A0} h     = dh + [A0, h]                              (h as a section of Ad P)
   (e1,w1)(e2,w2) = (e1 e2, e2^-1 w1 e2 + w2)                (draft (5.9), Aut(e^-1, w) := e^-1 w e)
   A . (e, w)   = A . e + w                                 (draft (5.10))
   (e, w) . A   = (A + w) . e^-1                            (draft (5.11))
   tau(h)       = (h, A0 - A0.h) = (h, -h^-1 d_{A0} h)      (draft (6.4))
   pi(e, w)     = e w e^-1 + (d_{A0} e) e^-1                (draft (6.11))
   g . gamma    = e (w + gamma) e^-1 + (d_{A0} e) e^-1      (draft (6.15))
   mu(e, w)     = (A0 + w, A0 . e)                          (draft (6.18))
   T(e, w)      = w - e^-1 d_{A0} e                         (draft (7.3))
All identities are polynomial; each difference is expanded and must vanish identically.
"""
import sympy as sp

x, y = sp.symbols('x y')
X = (x, y)


def M(a, b, c, d):
    return sp.Matrix([[a, b], [c, d]])


def unipotent_product(p, q):
    """Invertible with polynomial inverse: [[1,p],[0,1]] [[1,0],[q,1]]."""
    return M(1, p, 0, 1) * M(1, 0, q, 1)


def inv(g):
    r = sp.simplify(g.inv())
    return r.applyfunc(sp.expand)


# generic polynomial test data
e1 = unipotent_product(x**2 + 3*x*y - y, 2*x - y**2 + 1) * M(1, 0, x*y, 1)
e2 = unipotent_product(y**3 - x + 2, x*y + 3*y) * M(1, x - 2*y**2, 0, 1)
e3 = unipotent_product(1 + x - y, x**2 * y - 1)
h1 = unipotent_product(x - y**2, 3 + x*y)
h2 = unipotent_product(2*x*y - 1, y**2 + x)


def oneform(seed):
    a = sp.Integer(seed)
    return (M(x + a, y**2 - a*x, a*x*y, 1 - y), M(x**2 - a, a + y, x - a*y, x*y + a))


A0 = oneform(1)
w1, w2, w3, gam, Aany = oneform(2), oneform(3), oneform(5), oneform(7), oneform(11)


def add(a, b):
    return tuple(ai + bi for ai, bi in zip(a, b))


def neg(a):
    return tuple(-ai for ai in a)


def conj(g, a):           # g^-1 a g
    gi = inv(g)
    return tuple(gi * ai * g for ai in a)


def conj_inv(g, a):       # g a g^-1
    gi = inv(g)
    return tuple(g * ai * gi for ai in a)


def d(g):
    return tuple(g.diff(v) for v in X)


def gauge(A, h):          # A . h
    hi = inv(h)
    return tuple(hi * Ai * h + hi * dh for Ai, dh in zip(A, d(h)))


def dA0(h):               # d_{A0} h
    return tuple(dh + A0i * h - h * A0i for dh, A0i in zip(d(h), A0))


def lmul(g, a):           # g * a  (matrix times each component)
    return tuple(g * ai for ai in a)


def rmul(a, g):
    return tuple(ai * g for ai in a)


def mul(g1, g2):
    (a, u), (b, v) = g1, g2
    return (a * b, add(conj(b, u), v))


def ginv(g):
    a, u = g
    return (inv(a), neg(conj_inv(a, u)))


def act_right(A, g):
    e, w = g
    return add(gauge(A, e), w)


def act_left(g, A):
    e, w = g
    return gauge(add(A, w), inv(e))


def tau(h):
    return (h, neg(lmul(inv(h), dA0(h))))


def pi(g):
    e, w = g
    return add(conj_inv(e, w), rmul(dA0(e), inv(e)))


def act_N(g, gm):
    e, w = g
    return add(conj_inv(e, add(w, gm)), rmul(dA0(e), inv(e)))


def mu(g):
    e, w = g
    return (add(A0, w), gauge(A0, e))


def T(g):
    e, w = g
    return add(w, neg(lmul(inv(e), dA0(e))))


def is_zero(obj):
    if isinstance(obj, sp.MatrixBase):
        return all(sp.expand(v) == 0 for v in obj)
    return all(is_zero(o) for o in obj)


def diff(a, b):
    if isinstance(a, sp.MatrixBase):
        return a - b
    return tuple(diff(ai, bi) for ai, bi in zip(a, b))


def check(name, lhs, rhs):
    ok = is_zero(diff(lhs, rhs))
    print(f"{name:78s} {'holds' if ok else 'FAILS'}")
    return ok


def main():
    g1, g2, g3 = (e1, w1), (e2, w2), (e3, w3)
    results = []
    results.append(check("G1  (5.9) associativity (g1 g2) g3 = g1 (g2 g3)", mul(mul(g1, g2), g3), mul(g1, mul(g2, g3))))
    results.append(check("G1' inverse: g g^-1 = (1, 0)", mul(g1, ginv(g1)), (sp.eye(2), (sp.zeros(2), sp.zeros(2)))))
    results.append(check("G2  (5.10) right action: (A.g1).g2 = A.(g1 g2)", act_right(act_right(Aany, g1), g2),
                         act_right(Aany, mul(g1, g2))))
    results.append(check("G3  (5.11) left action: (g1 g2).A = g1.(g2.A)", act_left(mul(g1, g2), Aany),
                         act_left(g1, act_left(g2, Aany))))
    results.append(check("G4  (6.4) tau(h) = (h, A0 - A0.h) agrees with (h, -h^-1 d_A0 h)",
                         tau(h1)[1], add(A0, neg(gauge(A0, h1)))))
    results.append(check("G4' Lemma 6.1: tau(h1 h2) = tau(h1) tau(h2)", tau(h1 * h2), mul(tau(h1), tau(h2))))
    results.append(check("G5  (6.6)-(6.9) A0 . tau(h) = A0 (tau(H) = Stab(A0))", act_right(A0, tau(h1)), A0))
    results.append(check("G6  Lemma 6.3: pi(g tau(h)) = pi(g)", pi(mul(g1, tau(h1))), pi(g1)))
    results.append(check("G6' pi(g) = A0 - A0.g^-1 = -(g.A0 - A0)  (pi is the orbit map of the left action)",
                         pi(g1), add(A0, neg(act_right(A0, ginv(g1))))))
    results.append(check("G7  Lemma 6.4: (g1 g2).gamma = g1.(g2.gamma) on N", act_N(mul(g1, g2), gam),
                         act_N(g1, act_N(g2, gam))))
    results.append(check("G7' pi is G-equivariant: pi(g1 g2) = g1 . pi(g2)", pi(mul(g1, g2)), act_N(g1, pi(g2))))
    mg = mu(mul(g1, tau(h1)))
    m0 = mu(g1)
    results.append(check("G8  Lemma 6.5: mu(g tau(h)) = mu(g) . h (both components)", mg,
                         (gauge(m0[0], h1), gauge(m0[1], h1))))
    results.append(check("G9  (7.2) T = delta(mu(g)) = (A0 + w) - A0.e", T(g1), add(m0[0], neg(m0[1]))))
    results.append(check("G9' Lemma 7.2: T(g tau(h)) = h^-1 T(g) h", T(mul(g1, tau(h1))), conj(h1, T(g1))))
    # negative control: untilted subgroup (h, 0) does NOT act covariantly on T
    ok = is_zero(diff(T(mul(g1, (h1, (sp.zeros(2), sp.zeros(2))))), conj(h1, T(g1))))
    print(f"{'N1  control: T(g (h,0)) = h^-1 T(g) h  (expected to FAIL: untilted H)':78s} {'holds' if ok else 'fails (as expected)'}")
    assert not ok
    # residual of the control equals -h^-1 d_{A0} h, i.e. the 'common disease' of (2.8)
    res = diff(T(mul(g1, (h1, (sp.zeros(2), sp.zeros(2))))), conj(h1, T(g1)))
    print(f"{'N1'' residual equals -h^-1 d_A0 h':78s} {'holds' if is_zero(diff(res, neg(lmul(inv(h1), dA0(h1))))) else 'FAILS'}")
    print(f"\nall positive identities hold: {all(results)}")
    assert all(results)


if __name__ == '__main__':
    main()
