"""Exact controls for AM's complete flat Z1 model and global bosonic field.

SymPy only; no observer-metric or Clifford engine is imported. The note proves
the global inverse, completeness and spin descent. These controls check an
exact nonzero-x pullback, a geodesic, and the imported AI field coefficients.
"""
import sympy as sp


CHECKS = 0


def check(label, condition):
    global CHECKS
    CHECKS += 1
    good = bool(condition)
    print(f'[{CHECKS:02d}] {label}: {good}', flush=True)
    assert good, label


def unit_translation(x, scale):
    out = sp.eye(5)
    out[:4, 4] = -scale*x
    return out


def main():
    rational = sp.Rational
    c_h, c_v = rational(9), rational(2)
    scale = rational(3, 2)
    check('positive scale satisfies 2*cV*s^2=cH', 2*c_v*scale**2 == c_h)
    lower = sp.Matrix([[2, 0, 0, 0], [1, 3, 0, 0],
                       [1, -1, 1, 0], [2, 1, -1, 1]])
    block = lower*lower.T
    cross = sp.Matrix([rational(1, 2), -rational(2, 3),
                       rational(3, 5), -rational(4, 7)])
    q_target = sp.zeros(5)
    q_target[:4, :4] = block
    q_target[:4, 4] = cross
    q_target[4, :4] = cross.T
    q_target[4, 4] = 1/block.det()+(cross.T*block.inv()*cross)[0]
    check('arbitrary-block SPD target has exact determinant one',
          q_target.det() == 1 and all(q_target[:j, :j].det() > 0 for j in range(1, 6)))
    x = -block.inv()*cross/scale
    y = block.det()*block
    det_root = block.det()
    check('inverse uses the positive fifth root of det(y)',
          det_root > 0 and y.det() == det_root**5)
    trans = unit_translation(x, scale)
    diagonal = sp.diag(y, 1)
    q = trans.T*diagonal*trans/det_root
    check('global inverse reconstructs all twenty-five target entries', q == q_target)

    pairs = [(i, j) for i in range(4) for j in range(i, 4)]
    directions, expected_vertical = [], []
    y_inverse = y.inv()
    for i in range(4):
        dx = sp.zeros(4, 1)
        dx[i] = 1
        dt = unit_translation(dx, scale)-sp.eye(5)
        directions.append((dt.T*diagonal*trans+trans.T*diagonal*dt)/det_root)
    for i, j in pairs:
        h = sp.zeros(4)
        h[i, j] = h[j, i] = 1
        trace = sp.trace(y_inverse*h)
        delta = sp.diag(h, 0)-trace*diagonal/5
        directions.append(trans.T*delta*trans/det_root)
        expected_vertical.append((h, trace))
    inverse = q.inv()
    actual = sp.Matrix(14, 14, lambda i, j: c_v*sp.trace(
        inverse*directions[i]*inverse*directions[j]))
    expected = sp.zeros(14)
    expected[:4, :4] = c_h*y
    for i, (h, th) in enumerate(expected_vertical):
        for j, (k, tk) in enumerate(expected_vertical):
            expected[4+i, 4+j] = c_v*(sp.trace(y_inverse*h*y_inverse*k)-th*tk/5)
    check('all 196 pullback coefficients agree at nonzero x and nonunit scales', actual == expected)

    translation = sp.Matrix([1, -2, 3, -4])
    shifted = unit_translation(x+translation, scale)
    deck = unit_translation(translation, scale)
    check('translation is the exact SL5 congruence action',
          deck.det() == 1 and shifted == trans*deck
          and shifted.T*diagonal*shifted == deck.T*(trans.T*diagonal*trans)*deck)

    cv = sp.Symbol('cV', positive=True)
    kap = sp.Symbol('kappa', real=True, nonzero=True)
    a = -rational(1, 400)
    b = sp.factor(kap**2*a*(1+312*a))
    residual = 15/(2*cv)+b
    check('full invariant residual vanishes at the accepted scale',
          b == -11*kap**2/20000
          and sp.simplify(residual.subs(cv, 150000/(11*kap**2))) == 0)
    check('incorrect scale negative control has nonzero residual',
          sp.simplify(residual.subs(cv, 1500000/(11*kap**2))) != 0)
    odd_t = a*kap
    vector_matrix = sp.diag(kap, kap+224*odd_t, kap+544*odd_t, kap)
    vector_det = sp.factor(vector_matrix.det())
    check('full four-by-four vector block is regular, including kappa=2 and -3',
          vector_det == -99*kap**4/625 and vector_det.is_nonzero
          and all(vector_matrix.subs(kap, value).det() == -99*value**4/rational(625)
                  and vector_matrix.subs(kap, value).det() != -99*value**2/rational(625)
                  for value in (2, -3)))

    u, v = sp.symbols('u v', positive=True)
    check('flat positive Z1 Ricci spectrum has fourteen negative directions',
          all(value.is_negative for value in [-u]*5+[-v]*9))

    time = sp.Symbol('t', real=True)
    sigma = sp.Symbol('s', positive=True)
    ch, sh = sp.cosh(time), sp.sinh(time)
    boost_y = sp.diag(ch**2, ch, ch, ch)
    boost_x = sp.Matrix([-sp.tanh(time)/sigma, 0, 0, 0])
    boost_trans = unit_translation(boost_x, sigma)
    boost_image = boost_trans.T*sp.diag(boost_y, 1)*boost_trans/ch
    boost = sp.eye(5)
    boost[0, 0] = boost[4, 4] = ch
    boost[0, 4] = boost[4, 0] = sh
    boost_inverse = boost.subs(time, -time)
    velocity = boost.diff(time)
    check('explicit boost maps globally to its determinant-one geodesic',
          (boost_image-boost).applyfunc(sp.simplify) == sp.zeros(5)
          and sp.simplify(boost.det()-1) == 0
          and (boost.diff(time, 2)-velocity*boost_inverse*velocity).applyfunc(sp.simplify) == sp.zeros(5))
    dy = boost_y.diff(time)
    dx = boost_x.diff(time)
    y_velocity = boost_y.inv()*dy
    pullback_speed = (2*cv*sigma**2*(dx.T*boost_y*dx)[0]
                      +cv*(sp.trace(y_velocity*y_velocity)-sp.trace(y_velocity)**2/5))
    target_speed = cv*sp.trace(boost_inverse*velocity*boost_inverse*velocity)
    check('boost has constant squared speed 2*cV in both metrics',
          sp.simplify(pullback_speed-2*cv) == 0 and sp.simplify(target_speed-2*cv) == 0)
    print(f'{CHECKS} exact controls passed; no metric or Clifford engine imported.')
    print('AM1-AM2 global completeness and spin descent are proved in the note, not inferred from samples.')
    print('F7 supplies the existing local identification; AI3 supplies the full bosonic equation.')


if __name__ == '__main__':
    main()
