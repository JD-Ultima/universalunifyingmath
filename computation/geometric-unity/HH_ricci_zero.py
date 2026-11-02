"""Fresh exact certificates at every nonflat Ricci-zero four-metric point.

The hand degree bounds are A<=2 and RHS<=3 before interpolation. The required
union of 45 systems is assembled once in memory: 35 compatibility nodes and
ten additional p=-1 normalized rank nodes. There is no disk cache, raw stored
node matrix or scratch-path dependency. Both-singular chiral blocks use GG's
separately reviewed theorem; new rank data address the nonsingular case.
"""
import itertools
import os
import pathlib
import sys
import time
from fractions import Fraction as Fr

for _key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ.setdefault(_key, '1')

import sympy as sp

HERE = pathlib.Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import U_rank_one_frontier as U
import HH_rank_blocks as Rank
from HH_orbit_identity import verify_orbit_identity

p, q, r, s = sp.symbols('p q r s', real=True)
VARS = (p, q, r, s)
LAYOUT = (2, 3, 5, 7)
ROWS = (168, 337, 492)
CHECKS = []


def compositions(total, n):
    if n == 1:
        yield (total,)
    else:
        for first in range(total+1):
            for tail in compositions(total-first, n-1):
                yield (first,)+tail


COMPAT_NODES = tuple(node for total in range(4) for node in compositions(total, 4))
ALL_NODES = tuple(dict.fromkeys(COMPAT_NODES+tuple((sign,)+node for sign in (1, -1) for node in Rank.NODES)))
MON3 = tuple(sp.prod(variable**degree for variable, degree in zip(VARS, node)) for node in COMPAT_NODES)
MON2 = MON3[:15]
V2 = sp.Matrix([[mon.subs(dict(zip(VARS, node))) for mon in MON2] for node in COMPAT_NODES[:15]])
V3 = sp.Matrix([[mon.subs(dict(zip(VARS, node))) for mon in MON3] for node in COMPAT_NODES])
INV2, INV3 = V2.inv(), V3.inv()


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f'[{len(CHECKS):02d}] {label}: {ok}', flush=True)
    assert ok, label


def two_form(i, j):
    return [[Fr(int((a, b) == (i, j))-int((a, b) == (j, i))) for b in range(4)] for a in range(4)]


THETA = (two_form(0, 1), two_form(0, 2), two_form(0, 3))
TAU = (two_form(2, 3), two_form(3, 1), two_form(1, 2))
KEYS = tuple(itertools.product(range(4), repeat=4))


def curvature(node):
    pp, qq, rr, ss = map(Fr, node)
    electric = ((pp+rr)/2, (qq+ss)/2, -(pp+qq+rr+ss)/2)
    magnetic = ((pp-rr)/2, (qq-ss)/2, -(pp+qq-rr-ss)/2)
    return {key: sum(electric[i]*(THETA[i][key[0]][key[1]]*THETA[i][key[2]][key[3]]
                                  +TAU[i][key[0]][key[1]]*TAU[i][key[2]][key[3]])
                     +magnetic[i]*(THETA[i][key[0]][key[1]]*TAU[i][key[2]][key[3]]
                                   +TAU[i][key[0]][key[1]]*THETA[i][key[2]][key[3]]) for i in range(3)) for key in KEYS}


def system(node):
    R = curvature(node)
    base = U.X.weyl_base(R, (1, 1, 1, 1))
    data = U.X.point_data(base, U.diagonal(LAYOUT), Fr(1, 5), Fr(2), Fr(3), riem=True)
    return U.relaxed(data)


def interpolate(values, degree):
    if all(value == values[0] for value in values):
        return sp.Rational(values[0])
    monomials, inverse = (MON2, INV2) if degree == 2 else (MON3, INV3)
    coefficients = inverse*sp.Matrix([sp.Rational(value) for value in values])
    return sp.expand(sum(coefficient*mon for coefficient, mon in zip(coefficients, monomials)))


def compatibility(systems):
    node_systems = [systems[node] for node in COMPAT_NODES]
    coeff = [[interpolate([A[row][col] for A, _ in node_systems[:15]], 2) for col in range(42)] for row in ROWS]
    rhs = [interpolate([b[row] for _, b in node_systems], 3) for row in ROWS]
    allowed = {27, 32, 36, 39, 41}
    check('selected rows have only xi_13 and the same curvature-independent n term',
          all(poly == 0 for row in coeff for col, poly in enumerate(row) if col not in allowed)
          and all(row[col] == coeff[0][col] and not row[col].has(*VARS) for row in coeff for col in (32, 36, 39, 41)))
    check('all 30 unused coefficient nodes check the reconstructed selected rows',
          all(poly.subs(dict(zip(VARS, node))) == A[row][col] for node, (A, _) in systems.items()
              if node not in COMPAT_NODES[:15] for row, polys in zip(ROWS, coeff) for col, poly in enumerate(polys)))
    check('the ten additional normalized nodes check every selected RHS polynomial',
          all(poly.subs(dict(zip(VARS, node))) == b[row] for node, (_, b) in systems.items()
              if node not in COMPAT_NODES for row, poly in zip(ROWS, rhs)))
    P = sp.Poly(sp.expand(sp.det(sp.Matrix([[coeff[i][27], 1, rhs[i]] for i in range(3)]))), *VARS, domain=sp.QQ)
    check('the actual compatibility determinant is a homogeneous quartic with 35 terms',
          P.total_degree() == 4 and len(P.monoms()) == 35 and all(sum(mon) == 4 for mon in P.monoms()))
    print('  actual compatibility determinant P =', P.as_expr(), flush=True)
    orbit = verify_orbit_identity(P, VARS, emit=lambda text: print(text, flush=True))
    rejected = False
    try:
        verify_orbit_identity(2*P.as_expr(), VARS, emit=lambda _: None)
    except AssertionError:
        rejected = True
    check('the orbit certificate rejects an incorrectly doubled determinant', rejected)
    return orbit


def main():
    start = time.time()
    print('HH. Every positive Ricci-zero point with nonzero Riemann curvature', flush=True)
    print('Fresh portable path: 45 required node systems once in memory; no raw stored matrices or disk caches.', flush=True)
    basis = [curvature(tuple(int(i == j) for i in range(4))) for j in range(4)]
    check('all four curvature basis tensors satisfy Bianchi, pair symmetries and zero Ricci',
          all(W[(i, j, k, l)] == -W[(j, i, k, l)] == -W[(i, j, l, k)] == W[(k, l, i, j)]
              and W[(i, j, k, l)]+W[(i, k, l, j)]+W[(i, l, j, k)] == 0
              for W in basis for i, j, k, l in KEYS)
          and all(sum(W[(i, j, i, k)] for i in range(4)) == 0 for W in basis for j, k in itertools.product(range(4), repeat=2)))
    check('compatibility nodes and normalized-rank nodes match the prior total-degree spaces',
          len(COMPAT_NODES) == 35 and len(ALL_NODES) == 45 and Rank.LAYOUT == LAYOUT
          and abs(V2.det()) == 16 and abs(V3.det()) == 84934656)
    systems = {node: system(node) for node in ALL_NODES}
    orbit = compatibility(systems)
    rank_checks = 0
    for sign in (1, -1):
        matrices = [systems[(sign,)+node][0] for node in Rank.NODES]
        rank_checks += Rank.verify_rank_case(sign, matrices, emit=lambda text: print(text, flush=True))['checks']
    flat, flat_b = systems[(0, 0, 0, 0)]
    check('the complete flat quadratic control retains ranks 28/28',
          U.qmatrix(flat).rank() == 28 and not any(flat_b)
          and U.qmatrix([row+[rhs] for row, rhs in zip(flat, flat_b)]).rank() == 28)
    check('all 45 exact systems were assembled once', len(systems) == len(ALL_NODES) == 45)
    print('Both-singular chiral blocks use reviewed GG1-GG6. Otherwise choose a nonsingular W+, then a nonzero orbit term, THEN normalize its first eigenvalue p.', flush=True)
    print('The positive covariance transports rank42 and incompatibility TOGETHER from that selected normalized witness. No common original numeric fibre is asserted.', flush=True)
    print('Arbitrary higher jets use DD3/GG5, and the local fixed-minor step uses nonzero Riemann curvature. The flat origin and indefinite real forms remain outside this theorem.', flush=True)
    print(f'\n{len(CHECKS)} main + {rank_checks} rank + {orbit["checks"]} orbit checks, all True; runtime {time.time()-start:.1f} s', flush=True)


if __name__ == '__main__':
    main()
