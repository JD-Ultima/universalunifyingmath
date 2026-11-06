"""Exact sphere constants for the H6 orbit-average profile.

The group action, integral construction, dominated convergence and analytic
continuation are separate hand arguments. This script checks normalization
and constants only; it executes no previous helper.
"""
from time import perf_counter
import sympy as s


def main():
    start = perf_counter()
    t = s.Symbol('t', real=True)
    density = s.Rational(3, 4)*(1-t*t)
    zonal = ((7*t**3-3*t)/4, (33*t**5-30*t**3+5*t)/8)
    dimensions = (30, 91)
    degrees = (3, 5)
    assert s.integrate(density, (t, -1, 1)) == 1
    coefficients = []
    for degree, dimension, z in zip(degrees, dimensions, zonal):
        assert z.subs(t, 1) == 1 and z.subs(t, -t) == -z
        assert s.expand((1-t*t)*s.diff(z, t, 2)-4*t*s.diff(z, t)
                        +degree*(degree+3)*z) == 0
        assert s.integrate(z*z*density, (t, -1, 1)) == s.Rational(1, dimension)
        moment = 2*s.integrate(z*density, (t, 0, 1))
        coefficients.append(s.factor(dimension*moment))
        print('Degree', degree, ': zonal=', z, '; norm squared=', s.Rational(1, dimension))
        print('Signed hemisphere moment=', moment, '; normalized field limit=', coefficients[-1])
    assert s.integrate(zonal[0]*zonal[1]*density, (t, -1, 1)) == 0
    assert coefficients == [-s.Rational(15, 8), s.Rational(273, 128)]
    assert sum(coefficients) == s.Rational(33, 128)
    print('Orthogonality, sphere Casimir eigenvalues18/40 and both boundary constants verified.')
    print('The averaged lift of Bstar has boundary multiplier33/128; the affine seed has multiplier1.')
    print('No field profile, dominated-convergence argument or nonlinear equation is verified by these scalar checks.')
    print(f'Runtime {perf_counter()-start:.3f} seconds.')


if __name__ == '__main__':
    main()
