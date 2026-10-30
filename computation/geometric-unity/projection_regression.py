"""Exact regression for the (P) solution-set claim in I6 (counterexample I6.1).

Run: python3 computation/geometric-unity/projection_regression.py
Uses only the Python standard library; no randomness or repository imports.

Definitions: formalization/geometric-unity/F-zorro-reduction.md R-Z3 defines
(P) by Hodge dualising Upsilon and pulling back the form index along s(X).
Clifford indices are retained (see I7(c)'s 4*91 projected components).
I6 gives c1 = kappa*t + 312*t^2 - 104*u^2 and
c13 = u*(kappa+208*t) for W_A in (7,7). F5 gives Lambda_P=9*c_H/(2*c_V)
in Z2. These rational identities supply the algebraic part of the example.

The derivative term is evaluated from the definitions of D.1 and I.0:
  Phi^1 = sum_a e^a tensor gamma_a; Gamma = gamma_0 ... gamma_13;
  Psi(alpha ^ zeta) = sum_a sum_{b<c} Tr(alpha_a M_abc(zeta_bc));
  M^A_abc(Y) = -{gamma^a gamma^b gamma^c, Y};
  M^twist_abc(Y) = i[gamma^a gamma^b gamma^c, Ad_Gamma(Y)].
This is the same normalisation as xcheck_round5.py, implemented here directly.
The printed twisted coefficients omit their common, nonzero phase i.

A local extension u=1/2+rho(y), with rho zero on the section and d rho normal
there, exists in a tubular chart. The test is for its first jet along s(X).
It does not claim a full 14D solution on a neighbourhood.
"""
from fractions import Fraction
from itertools import combinations

N = 14
FULL = (1 << N) - 1
ETA = (1,) * 7 + (-1,) * 7


def multiply(left, right):
    """Integer Clifford multiplication: gamma_A gamma_B = sign gamma_(A xor B)."""
    result = {}
    for a, av in left.items():
        for b, bv in right.items():
            inversions = sum((a >> (j + 1)).bit_count() for j in range(N) if b & (1 << j))
            sign = -1 if inversions % 2 else 1
            sign *= (-1) ** ((a & b & sum(1 << j for j in range(7, N))).bit_count())
            key = a ^ b
            result[key] = result.get(key, 0) + sign * av * bv
    return {key: value for key, value in result.items() if value}


def add(left, right, factor=1):
    result = dict(left)
    for key, value in right.items():
        result[key] = result.get(key, 0) + factor * value
    return {key: value for key, value in result.items() if value}


def scale(value, factor):
    return {key: factor * coefficient for key, coefficient in value.items() if factor * coefficient}


def derivative_output(normal, twisted):
    """Dual coefficients at each form slot for Shiab(du ^ Phi^1 Gamma), du=e^normal."""
    phi1gamma = [multiply({1 << a: 1}, {FULL: 1}) for a in range(N)]
    output = [dict() for _ in range(N)]
    for b, c in combinations(range(N), 2):
        if normal not in (b, c):
            continue
        value = phi1gamma[c] if b == normal else scale(phi1gamma[b], -1)
        if twisted:
            value = {key: (-coefficient if key.bit_count() % 2 else coefficient)
                     for key, coefficient in value.items()}
        for a in range(N):
            if a in (b, c):
                continue
            g3 = multiply(multiply({1 << a: ETA[a]}, {1 << b: ETA[b]}), {1 << c: ETA[c]})
            image = add(multiply(g3, value), multiply(value, g3), -1 if twisted else 1)
            output[a] = add(output[a], image, 1 if twisted else -1)
    return output


def main():
    # Flat Lorentzian base, Z2; positive c_H,c_V and lambda>1/4 give (7,7).
    c_h, c_v = Fraction(1), Fraction(3, 52)
    kappa, t, u = Fraction(208), Fraction(-1), Fraction(1, 2)
    c1 = kappa * t + 312 * t * t - 104 * u * u
    c13 = u * (kappa + 208 * t)
    lambda_p = 9 * c_h / (2 * c_v)
    print("Z2, flat Lorentzian base, lambda=1/2, c_H=1, c_V=3/52, kappa=208")
    print("t=-1, u|s=1/2; dt=0; du=e^1 normal to the section")
    print(f"c1={c1}; c13={c13}; Lambda_P={lambda_p}; Lambda_4={lambda_p - c_h * c1}")
    assert c1 == lambda_p and c13 == 0

    # H has signature (1,3); index 1 belongs to its nondegenerate vertical complement.
    tangents = (0, 7, 8, 9)
    wa = derivative_output(normal=1, twisted=False)
    twist = derivative_output(normal=1, twisted=True)
    print(f"W_A derivative output vanishes at every form slot: {not any(wa)}")
    print("Twist derivative output at retained tangential form slots (common phase i omitted):")
    for a in tangents:
        print(f"  a={a}: {dict(sorted(twist[a].items()))}")
    assert not any(wa)
    assert all(twist[a] for a in tangents)
    assert all(key.bit_count() == 12 for a in tangents for key in twist[a])
    print("Exact result: this first jet satisfies (P) for W_A and fails (P) for W_A+W_C.Ad_Gamma.")
    print("The 4D algebraic equation agrees; the projected normal-derivative constraint differs.")


if __name__ == "__main__":
    main()
