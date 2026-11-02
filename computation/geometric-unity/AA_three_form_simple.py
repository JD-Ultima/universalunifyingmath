"""Exact exterior-algebra controls for AA7's real Euclidean 14-dimensional result.

For a real three-form omega, Z_a=(i_a omega) wedge omega and
F=sum_a (i_a omega) wedge (i_a omega), the condition H5(Z)=Z-w(F)/10=0
is equivalent to omega being zero or decomposable. The universal proof is
the annihilator/support argument in AA7, not a classification by sampling.

This helper uses only standard-library exact rational exterior algebra.
It checks the elementary symplectic commutator on every basis form of the
degrees and support dimensions used in the proof, projector conventions,
decomposable and nondecomposable controls, and the dimension-seven boundary.
No Clifford engine, optimization, random inputs or off-repository files.
"""
import time
from fractions import Fraction as F
from itertools import combinations

N = 14
CHECKS = []


def add(left, right, factor=1):
    out = dict(left)
    for mask, coefficient in right.items():
        out[mask] = out.get(mask, F(0)) + factor * coefficient
    return {mask: coefficient for mask, coefficient in out.items() if coefficient}


def scale(value, factor):
    return {mask: coefficient * factor for mask, coefficient in value.items() if coefficient * factor}


def form(indices, coefficient=1):
    return {sum(1 << index for index in indices): F(coefficient)}


def wedge(left, right):
    out = {}
    for a, ca in left.items():
        for b, cb in right.items():
            if a & b:
                continue
            inversions = sum((a >> (j + 1)).bit_count() for j in range(N) if b & (1 << j))
            out[a | b] = out.get(a | b, F(0)) + (-1)**inversions * ca * cb
    return {mask: coefficient for mask, coefficient in out.items() if coefficient}


def contract(index, value):
    return {mask ^ (1 << index): (-1)**((mask & ((1 << index) - 1)).bit_count()) * coefficient
            for mask, coefficient in value.items() if mask & (1 << index)}


def sum_forms(values):
    out = {}
    for value in values:
        out = add(out, value)
    return out


def lifts(value, dimension, exterior):
    return [wedge(form((a,)), value) if exterior else contract(a, value) for a in range(dimension)]


def traces(field):
    return (sum_forms(wedge(form((a,)), value) for a, value in enumerate(field)),
            sum_forms(contract(a, value) for a, value in enumerate(field)))


def hook(field):
    dimension = len(field)
    alt, ctr = traces(field)
    return [add(add(value, contract(a, alt), F(-1, 6)),
                wedge(form((a,)), ctr), -F(1, dimension - 4))
            for a, value in enumerate(field)]


def tensors(omega, dimension=N):
    contractions = [contract(a, omega) for a in range(dimension)]
    z = [wedge(value, omega) for value in contractions]
    f4 = sum_forms(wedge(value, value) for value in contractions)
    return z, f4


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def symplectic_checks():
    checked = 0
    for r in range(2, 7):
        theta = sum_forms(form((2 * i, 2 * i + 1)) for i in range(r))

        def adjoint(value):
            return sum_forms(contract(2 * i + 1, contract(2 * i, value)) for i in range(r))

        count = 0
        for degree in range(4):
            for indices in combinations(range(2 * r), degree):
                eta = form(indices)
                commutator = add(adjoint(wedge(theta, eta)), wedge(theta, adjoint(eta)), -1)
                assert commutator == scale(eta, r - degree)
                count += 1
        checked += count
        check(f"symplectic support dimension {2*r}: all {count} basis commutators have coefficient r-j", count > 0)
    check("the symplectic commutator was checked on all 625 relevant basis forms", checked == 625)

    # On S4, Lambda2 has dimension six and L_theta maps to the one-dimensional
    # Lambda4. Its nonzero image therefore has rank one and kernel dimension five.
    theta = add(form((0, 1)), form((2, 3)))
    images = [wedge(theta, form(indices)) for indices in combinations(range(4), 2)]
    check("rank-four contraction: primitive Lambda2(S4) has dimension five",
          len(images) == 6 and any(images) and all(set(value) <= {15} for value in images))
    check("in dimension 14 the r=2 annihilator dimension bound is positive", (14 - 1 - 4) - 5 == 4)


def exterior_checks():
    examples = (
        ("zero", {}, True),
        ("simple coordinate form", form((0, 1, 2)), True),
        ("simple rational change of factors",
         wedge(add(form((0,)), form((5,)), F(2, 3)),
               wedge(add(form((1,)), form((7,)), F(-3, 5)),
                     add(form((2,)), form((9,)), F(4, 7)))), True),
        ("disjoint triples with F=0", add(form((0, 1, 2)), form((3, 4, 5))), False),
        ("shared-index nondecomposable form", add(form((0, 1, 2)), form((0, 3, 4))), False),
    )
    for name, omega, expected_zero in examples:
        z, f4 = tensors(omega)
        alt, ctr = traces(z)
        hz = hook(z)
        check(name + ": traces, projector and expected hook status",
              not alt and ctr == f4 and hook(hz) == hz and not any(traces(hz))
              and (not any(hz)) == expected_zero)
        residual = sum_forms(wedge(contract(a, omega), hz[a]) for a in range(N))
        check(name + ": exact (n-7) wedge identity",
              residual == scale(wedge(f4, omega), F(N - 7, N - 4)))

    omega = add(form((0, 1, 2)), form((0, 3, 4)))
    z, f4 = tensors(omega)
    check("the explicit AA6 excluded form has F=2e1234 and H5(Z)_5=-e12345/5",
          f4 == form((1, 2, 3, 4), 2) and hook(z)[5] == form((1, 2, 3, 4, 5), F(-1, 5)))
    check("an unused direction annihilates omega and F, but detects a nonzero hook component",
          not contract(13, omega) and not contract(13, f4)
          and hook(z)[13] == scale(wedge(form((13,)), f4), F(-1, 10))
          and bool(hook(z)[13]))

    # Same explicit seven-dimensional form in two ambient dimensions. This
    # is a boundary control, not a dimension-independent classification.
    phi = sum_forms((form((0, 1, 2)), form((0, 3, 4)), form((0, 5, 6)),
                     form((1, 3, 5)), form((1, 4, 6), -1),
                     form((2, 3, 6), -1), form((2, 4, 5), -1)))
    z7, f7 = tensors(phi, 7)
    z14, f14 = tensors(phi, 14)
    check("the seven-dimensional comparison form has zero seven-dimensional hook but F!=0",
          not any(hook(z7)) and bool(f7))
    check("embedding that form in dimension 14 has nonzero hook", f14 == f7 and bool(any(hook(z14))))
    residual = sum_forms(wedge(contract(a, phi), hook(z14)[a]) for a in range(N))
    expected = scale(wedge(f14, phi), F(N - 7, N - 4))
    check("the (n-7) wedge identity also holds when its seven-form side is nonzero",
          residual == expected and bool(expected))

    # Checks the isolated x-containing term in the final decomposability
    # step: x=e0,u=e1,v=e2,sigma=e34,lambda=1.
    perturbed = add(form((0, 1, 2)), form((1, 3, 4)))
    square = wedge(contract(1, perturbed), contract(1, perturbed))
    x_part = {mask: coefficient for mask, coefficient in square.items() if mask & 1}
    expected = scale(wedge(form((2,)), wedge(form((0,)), form((3, 4)))), 2)
    check("the rank-two contraction proof has the stated unique x-containing square term",
          x_part == expected and bool(expected))


def main():
    start = time.time()
    print("AA7. Exterior-algebra controls for the real Euclidean dimension-14 hook zero locus", flush=True)
    symplectic_checks()
    exterior_checks()
    print("The proof classifies the hook zero locus as zero or decomposable real three-forms in positive dimension 14.", flush=True)
    print("This remains a necessary pointwise condition for AA6; no field-equation or PDE sufficiency is asserted.", flush=True)
    print(f"\n{len(CHECKS)} checks, all True; runtime {time.time()-start:.1f} s", flush=True)


if __name__ == "__main__":
    main()
