"""Portable exact spectral-orbit identity verifier for a future main driver.

Accepts the freshly reconstructed rational quartic and its four variables.
The identity is checked by Fraction polynomial arithmetic, without an engine,
Groebner calculation, stored source polynomial or off-repository path.
"""
from fractions import Fraction
from itertools import permutations, product
import json
from pathlib import Path

DATA = Path(__file__).with_name('HH_orbit_certificate.json')
ZERO = (0, 0, 0, 0)


def add(left, right):
    result = dict(left)
    for exponent, coefficient in right.items():
        result[exponent] = result.get(exponent, Fraction(0)) + coefficient
        if not result[exponent]:
            del result[exponent]
    return result


def multiply(left, right):
    result = {}
    for first, a in left.items():
        for second, b in right.items():
            exponent = tuple(x+y for x, y in zip(first, second))
            result[exponent] = result.get(exponent, Fraction(0)) + a*b
    return {exponent: coefficient for exponent, coefficient in result.items() if coefficient}


def power(polynomial, exponent):
    result = {ZERO: Fraction(1)}
    for _ in range(exponent):
        result = multiply(result, polynomial)
    return result


def substitute(polynomial, images):
    """Simultaneous substitution, by four independent monomial expansions."""
    result = {}
    for exponent, coefficient in polynomial.items():
        term = {ZERO: coefficient}
        for image, count in zip(images, exponent):
            term = multiply(term, power(image, count))
        result = add(result, term)
    return result


def verify_orbit_identity(polynomial, variables, emit=print):
    """Verify the fixed certificate against the main driver's actual quartic.

    `polynomial` may be a SymPy expression or Poly; `variables` gives the two
    leading eigenvalues of W+ followed by those of W-. Their names are free.
    """
    count = 0

    def check(label, condition):
        nonlocal count
        count += 1
        assert bool(condition), 'spectral-orbit identity: '+label

    variables = tuple(variables)
    assert len(variables) == len(set(variables)) == 4
    rational = polynomial.as_poly(*variables, domain='QQ')
    assert rational is not None
    source = {tuple(exponent): Fraction(str(coefficient)) for exponent, coefficient in rational.as_dict().items() if coefficient}
    check('the fresh input is an exact homogeneous quartic', source and all(sum(e) == 4 for e in source))

    data = json.loads(DATA.read_text(encoding='utf-8'))
    assert data['schema_version'] == 1 and data['variable_order'] == ['p', 'q', 'r', 's']
    seed_data = data['seed']
    prefactor = Fraction(seed_data['numerator_factor'], seed_data['denominator'])
    seed = {}
    for term in seed_data['terms']:
        exponent = tuple(term['exponents'])
        assert len(exponent) == 4 and all(isinstance(e, int) and e >= 0 for e in exponent)
        assert sum(exponent) == 2 and exponent not in seed
        seed[exponent] = prefactor*Fraction(term['coefficient'])
    check('the recorded fixed seed is a nonzero exact quadratic', seed and all(seed.values()))

    expected_frames = tuple((a, b) for a in permutations(range(3)) for b in permutations(range(3)))
    frames = tuple((tuple(item['plus']), tuple(item['minus'])) for item in data['frames'])
    check('all 36 independent spectral permutations occur once', frames == expected_frames and len(set(frames)) == 36)

    coordinates = [{tuple(int(i == j) for i in range(4)): Fraction(1)} for j in range(4)]
    p, q, r, s = coordinates
    negative_plus = {e: -c for e, c in add(p, q).items()}
    negative_minus = {e: -c for e, c in add(r, s).items()}
    plus, minus = (p, q, negative_plus), (r, s, negative_minus)
    X = add(add(power(p, 2), power(q, 2)), power(negative_plus, 2))
    Y = add(add(power(r, 2), power(s, 2)), power(negative_minus, 2))
    total = {}
    preserve_norms = True
    for first, second in frames:
        images = (plus[first[0]], plus[first[1]], minus[second[0]], minus[second[1]])
        preserve_norms &= substitute(X, images) == X and substitute(Y, images) == Y
        total = add(total, multiply(substitute(source, images), substitute(seed, images)))
    check('the simultaneous permutations preserve X and Y exactly', preserve_norms)

    target = power(add(X, Y), 3)
    exponents6 = tuple(e for e in product(range(7), repeat=4) if sum(e) == 6)
    check('all output terms have degree six', all(sum(e) == 6 for e in total))
    check('all 84 degree-six slots, including zero target slots, agree', len(exponents6) == 84
          and all(total.get(e, Fraction(0)) == target.get(e, Fraction(0)) for e in exponents6))
    emit(f'  orbit certificate: {count} exact checks; 36 frames, all 84 degree-six coefficient slots agree.')
    return dict(checks=count, frames=36, degree_six_slots=84,
                zero_target_slots=sum(e not in target for e in exponents6),
                source_terms=len(source), seed_terms=len(seed))
