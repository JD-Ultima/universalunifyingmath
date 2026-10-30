"""GU, Prop. C3' / referee issue GU-4 (reports/reviews/REVIEW-2026-10-30.md): what the tilt-sign repair changes.
Imports the exact GL(2)-over-R^2 model of gauge_group_identities.py unchanged (conventions R-C1 are listed there).

R1  With tau' (repaired (6.4)) and the printed pi (6.11), Lemma 6.3 fails with the exact residual
       pi(g tau'(h)) - pi(g) = 2 e (d_A0 h) h^-1 e^-1 .
R2  Repair A = tau' together with pi'(e, w) = e w e^-1 - (d_A0 e) e^-1: Lemmas 6.1, 6.3, 6.5, 7.2 all hold.
R3  Repair B = printed tau together with mu''(e, w) = (A0 - w, A0.e) and T''(e, w) = -w - e^-1 d_A0 e (flipped w in
    (6.18) and (7.3)): Lemmas 6.1, 6.3, 6.5, 7.2 all hold; T'' is still delta(mu'').
R4  No single printed formula among (6.4), (6.11), (6.18), (7.3) can be changed alone: with tau', (6.11) fails (R1);
    with tau, (6.18) and (7.3) fail (G8, G9' of gauge_group_identities.out).
R5  The left action of Lemma 6.4. The draft obtains (6.15) "as a benefit of this homogeneous principal fibration",
    i.e. as the action induced on N = G/tau(H) through pi: pi(g1 g2) = g1 . pi(g2) (G7'). Under repair B this holds
    as printed. Under repair A the induced action is g * gamma = e (w + gamma) e^-1 - (d_A0 e) e^-1, so (6.15) must
    change too if it is to stay the induced action (it remains a left action either way, G7).
R6  Other printed formulas that restate tau or mu (added after review 2026-10-30, N1).
    (6.13): g . tau(h) = (e h, h^-1 w h - h^-1 d_A0 h). It holds for tau (repair B) and fails for tau'
    (repair A), so repair A also changes (6.13), the expansion of (6.4).
    (6.22): A_w = nabla0 + w, B_w = nabla0 + e^-1 (nabla0 e), i.e. (A, B) = mu(g). It holds for the printed mu
    (repair A) and fails for mu'' (repair B), so repair B also changes (6.22).
R7  (10.6): d/ds g . tau(exp(s gamma)) at s = 0 is printed as (e gamma, d_A0 gamma - [gamma, w]) = (e gamma, d_A gamma),
    computed there with the tau' sign (exp(s gamma)^-1 d_A0 exp(s gamma)). Under repair A (tau', A = A0 + w) this holds.
    Under repair B (tau, A'' = A0 - w) the derivative is (e gamma, -d_A0 gamma - [gamma, w]) = (e gamma, -d_A'' gamma),
    so (10.6) conflicts with repair B. Exact: h_s = 1 + s gamma and h_s^-1 = 1 - s gamma agree with exp(+-s gamma)
    to first order, and the s-derivative at s = 0 sees only first order.
All identities are exact polynomial identities (sympy).
"""
import sympy as sp
import gauge_group_identities as gg
from gauge_group_identities import (add, neg, mul, tau, tau_p, pi, mu, T, gauge, conj, conj_inv, rmul, lmul,
                                    inv, dA0, d, check, diff, A0, e1, e2, h1, h2, w1, w2, M, x, y)

s = sp.Symbol('s')
gamma = M(x * y - 2, x + y ** 2, 3 * x - y, x ** 2 + 1)   # generic section of ad P (polynomial 2x2)


def pi_p(g):
    e, w = g
    return add(conj_inv(e, w), neg(rmul(dA0(e), inv(e))))


def act_N_p(g, gm):       # (6.15) with the sign of (d_A0 e) e^-1 flipped
    e, w = g
    return add(conj_inv(e, add(w, gm)), neg(rmul(dA0(e), inv(e))))


def mu_pp(g):
    e, w = g
    return (add(A0, neg(w)), gauge(A0, e))


def T_pp(g):
    e, w = g
    return add(neg(w), neg(lmul(inv(e), dA0(e))))


def right_tilt_first_order(g, sign):
    """g . (h_s, sign h_s^-1 d_A0 h_s) with h_s = 1 + s gamma, h_s^-1 := 1 - s gamma (exact to first order in s)."""
    e, w = g
    hs = sp.eye(2) + s * gamma
    hsi = sp.eye(2) - s * gamma
    return (e * hs, tuple(hsi * wi * hs + sign * hsi * c for wi, c in zip(w, dA0(hs))))


def dds0(obj):
    if isinstance(obj, sp.MatrixBase):
        return obj.diff(s).subs(s, 0).applyfunc(sp.expand)
    return tuple(dds0(o) for o in obj)


def d_conn(A, gm):        # d_A gamma = d gamma + [A, gamma], A a connection potential
    return tuple(dg + Ai * gm - gm * Ai for dg, Ai in zip(d(gm), A))


def bracket(gm, a):       # [gamma, a] componentwise
    return tuple(gm * ai - ai * gm for ai in a)


def main():
    g1, g2 = (e1, w1), (e2, w2)
    print("=== R1. tau' with the printed pi (6.11) ===")
    res = diff(pi(mul(g1, tau_p(h1))), pi(g1))
    check("R1  pi(g tau'(h)) - pi(g) = 2 e (d_A0 h) h^-1 e^-1", res,
          tuple(2 * x for x in conj_inv(e1, rmul(dA0(h1), inv(h1)))))
    print("\n=== R2. Repair A: tau' with pi' ===")
    check("R2a Lemma 6.1: tau'(h1 h2) = tau'(h1) tau'(h2)", tau_p(h1 * h2), mul(tau_p(h1), tau_p(h2)))
    check("R2b Lemma 6.3: pi'(g tau'(h)) = pi'(g)", pi_p(mul(g1, tau_p(h1))), pi_p(g1))
    check("R2c Lemma 6.5: mu(g tau'(h)) = mu(g).h", mu(mul(g1, tau_p(h1))),
          tuple(gauge(m, h1) for m in mu(g1)))
    check("R2d Lemma 7.2: T(g tau'(h)) = h^-1 T(g) h", T(mul(g1, tau_p(h1))), conj(h1, T(g1)))
    print("\n=== R3. Repair B: printed tau with mu'', T'' (w flipped in (6.18) and (7.3)) ===")
    check("R3a Lemma 6.3 (printed pi): pi(g tau(h)) = pi(g)", pi(mul(g1, tau(h1))), pi(g1))
    check("R3b Lemma 6.5: mu''(g tau(h)) = mu''(g).h", mu_pp(mul(g1, tau(h1))),
          tuple(gauge(m, h1) for m in mu_pp(g1)))
    check("R3c Lemma 7.2: T''(g tau(h)) = h^-1 T''(g) h", T_pp(mul(g1, tau(h1))), conj(h1, T_pp(g1)))
    check("R3d T'' = delta(mu'') = first - second component", T_pp(g1), add(mu_pp(g1)[0], neg(mu_pp(g1)[1])))
    print("\n=== R4. One formula alone does not suffice ===")
    check("R4a tau' with printed pi: Lemma 6.3 (expected false)", pi(mul(g1, tau_p(h1))), pi(g1))
    check("R4b tau with printed mu: Lemma 6.5 (expected false)", mu(mul(g1, tau(h1))),
          tuple(gauge(m, h1) for m in mu(g1)))
    check("R4c tau with printed T: Lemma 7.2 (expected false)", T(mul(g1, tau(h1))), conj(h1, T(g1)))
    print("\n=== R5. The left action (6.15) induced through the projection ===")
    check("R5a repair B: pi(g1 g2) = g1 . pi(g2) with printed (6.15)", pi(mul(g1, g2)), gg.act_N(g1, pi(g2)))
    check("R5b repair A: pi'(g1 g2) = g1 . pi'(g2) with printed (6.15) (expected false)",
          pi_p(mul(g1, g2)), gg.act_N(g1, pi_p(g2)))
    check("R5c repair A: pi'(g1 g2) = g1 * pi'(g2) with flipped (6.15)", pi_p(mul(g1, g2)), act_N_p(g1, pi_p(g2)))
    check("R5d flipped (6.15) is a left action: (g1 g2) * gamma = g1 * (g2 * gamma)",
          act_N_p(mul(g1, g2), gg.gam), act_N_p(g1, act_N_p(g2, gg.gam)))
    print("\n=== R6. Printed restatements of tau and mu: (6.13) and (6.22) ===")
    rhs613 = (e1 * h1, add(conj(h1, w1), neg(lmul(inv(h1), dA0(h1)))))
    check("R6a (6.13) as printed, with tau (repair B)", mul(g1, tau(h1)), rhs613)
    check("R6b (6.13) as printed, with tau' (repair A) (expected false)", mul(g1, tau_p(h1)), rhs613)
    rhs622 = (add(A0, w1), add(A0, lmul(inv(e1), dA0(e1))))
    check("R6c (6.22) as printed, (A_w, B_w) = mu(g) (repair A)", mu(g1), rhs622)
    check("R6d (6.22) as printed, (A_w, B_w) = mu''(g) (repair B) (expected false)", mu_pp(g1), rhs622)
    check("R6e B_w = nabla0 + e^-1 (nabla0 e) under mu'' (second entry unchanged)", mu_pp(g1)[1], rhs622[1])
    print("\n=== R7. The infinitesimal tilted action (10.6) ===")
    rhs106 = (e1 * gamma, add(dA0(gamma), neg(bracket(gamma, w1))))
    dA = dds0(right_tilt_first_order(g1, +1))
    dB = dds0(right_tilt_first_order(g1, -1))
    check("R7a repair A (tau'): d/ds g.tau'(e^{s gamma}) = (e gamma, d_A0 gamma - [gamma, w])", dA, rhs106)
    check("R7b   and = (e gamma, d_A gamma) with A = mu(g)[0] = A0 + w", dA, (e1 * gamma, d_conn(mu(g1)[0], gamma)))
    check("R7c repair B (tau): d/ds g.tau(e^{s gamma}) = printed (10.6) (expected false)", dB, rhs106)
    check("R7d   it is (e gamma, -d_A'' gamma) with A'' = mu''(g)[0] = A0 - w", dB,
          (e1 * gamma, neg(d_conn(mu_pp(g1)[0], gamma))))


if __name__ == '__main__':
    main()
