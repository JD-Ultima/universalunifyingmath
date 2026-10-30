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
All identities are exact polynomial identities (sympy).
"""
import gauge_group_identities as gg
from gauge_group_identities import (add, neg, mul, tau, tau_p, pi, mu, T, gauge, conj, conj_inv, rmul, lmul,
                                    inv, dA0, check, diff, A0, e1, e2, h1, h2, w1, w2)


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


if __name__ == '__main__':
    main()
