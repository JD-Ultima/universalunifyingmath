"""Exact Banach certificate for a real mixed scalar-stress point.

Only canonical integer rows exported by BL_symmetric_mixed are inputs.
The center and radius are rational, and all derivative and self-map bounds
are checked exactly. No numerical solver or earlier main is run; this
helper writes stdout only and does not read recorded output.
"""
import hashlib
import itertools
import json
from pathlib import Path
from time import perf_counter
import sys
import sympy as s

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE))
import BL_symmetric_mixed as REDUCTION

ROW_JSON_SHA = '8d365144715aac0307b26b0051d5c59610bc6ead63a903f66d727d005788d6bb'
CENTER_NUMERATORS = (4780996876601263938643, -17430528795097021097941,
                     -33861696489766133108980, 289635016937541269176233)
CENTER_DENOMINATOR = 10**24


def main():
    start = perf_counter()
    rows = REDUCTION.QUADRATIC_ROWS
    pairs = REDUCTION.AMPLITUDE_PAIRS
    assert tuple(pairs) == tuple(itertools.combinations_with_replacement(range(6), 2))
    assert len(rows) == 4 and all(len(row) == 21 and all(isinstance(v, int) for v in row) for row in rows)
    row_bytes = json.dumps(rows,separators=(',',':')).encode('ascii')
    assert row_bytes == REDUCTION.canonical_rows_bytes()
    assert hashlib.sha256(row_bytes).hexdigest() == ROW_JSON_SHA
    x = s.symbols('a0:4')
    amplitudes = (*x, -s.Rational(1, 100), s.Integer(1))
    f = s.Matrix([sum(coefficient*amplitudes[i]*amplitudes[j]
                     for coefficient, (i, j) in zip(row, pairs)) for row in rows])
    J = f.jacobian(x)
    center = s.Matrix([s.Rational(v, CENTER_DENOMINATOR) for v in CENTER_NUMERATORS])
    substitutions = dict(zip(x, center))
    Jc = J.subs(substitutions)
    assert Jc.det() != 0
    R = Jc.inv()
    assert R*Jc == s.eye(4)
    correction = R*f.subs(substitutions)
    eta = max(abs(value) for value in correction)
    H = [s.hessian(value, x) for value in f]
    assert all(not value.free_symbols for hessian in H for value in hessian)
    # I-RJ(center+h)=-sum_k (R H[:, :, k])*h_k exactly, since f is quadratic.
    L = max(sum(abs(sum(R[i, ell]*H[ell][j, k] for ell in range(4)))
                for j in range(4) for k in range(4)) for i in range(4))
    radius = s.Rational(1, 10**10)
    q = radius*L
    image_radius = eta+q*radius
    assert 0 < eta < s.Rational(1, 10**23)
    assert 0 < L < 10**9
    assert 0 < q < s.Rational(1, 10) < 1
    assert image_radius < radius
    assert all(abs(value) > radius for value in center)
    assert not any(value.atoms(s.Float) for matrix in (f, Jc, R, correction, *H) for value in matrix)
    print('BL: exact real-zero certificate for the four canonical mixed-scalar-stress quadrics.')
    print('Canonical integer-row JSON SHA256:', hashlib.sha256(row_bytes).hexdigest())
    print('Fixed projective coordinates: a4=-1/100, a5=1.')
    print('Rational center:', tuple(center))
    print('Closed infinity-norm box radius:', radius)
    print('R=J(center)^(-1) verified exactly; determinant J(center):', Jc.det())
    print('eta=||R*f(center)||_infinity:', eta, '<', s.Rational(1, 10**23))
    print('L=max_i sum_jk |sum_l R_il H_ljk|:', L, '<', 10**9)
    print('q=radius*L:', q, '<', s.Rational(1, 10), '<1')
    print('eta+q*radius:', image_radius, '<radius:', radius)
    print('Exact positive self-map margin:', radius-image_radius)
    print('All four variable-coordinate intervals exclude zero.')
    print('Banach contraction gives one and only one real common zero in this closed rational box.')
    print('This is a real pointwise scalar-stress zero, not an explicit rational tensor or a field solution.')
    print('All inequalities use rational arithmetic; no recorded output or numerical solver is an input.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
