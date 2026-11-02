"""Exact algebraic controls for AJ's curved-base Lambda3 endpoint theorem.

The geometric inputs are F4 on radial fibres, N1's Euler Ricci identity,
S4 rationality and S5(b) Einstein-open flatness. This helper verifies their
algebraic combination; it does not independently reconstruct the metric.
"""
import sympy as sp

CHECKS = []


def check(label, condition):
    ok = bool(condition)
    CHECKS.append(ok)
    print(f"[{len(CHECKS):02d}] {label}: {ok}", flush=True)
    assert ok, label


def polynomial(B):
    tau = sp.trace(B)
    return (3300 * B * B - 380 * tau * B
            + 9 * tau**2 * sp.eye(B.rows)).applyfunc(sp.expand)


def main():
    print("AJ: first Lambda3 endpoint t/kappa=-1/400, u=0")
    print(f"SymPy {sp.__version__}; exact symbolic controls; geometry imported.")
    z, delta, tau = sp.symbols("z delta tau", real=True)
    parallel, perpendicular = -sp.Rational(220, 3)*z, -180*z
    B = sp.diag(*([parallel]*3 + [perpendicular]*11))
    check("simple-form eigenvalues give trace -2200*z", sp.trace(B) == -2200*z)
    check("the full spectral polynomial vanishes identically",
          polynomial(B) == sp.zeros(14))
    check("trace-square elimination retains the zero-form case",
          sp.expand(1650*sp.trace(B*B) - 127*sp.trace(B)**2) == 0)
    check("scalar limit forces delta=0",
          polynomial(delta*sp.eye(14)) == -256*delta**2*sp.eye(14))

    q, cH, cV, h = sp.symbols("q cH cV h", positive=True)
    ricci = sp.symbols("rho0:4", real=True)
    R = sum(ricci)
    radial = [sp.Rational(15, 2)/cV-h + (rho-R/2)/(q*cH)
              for rho in ricci]
    radial += [sp.Rational(15, 2)/cV-h-R/(2*q*cH)]*10
    limits = [sp.limit(value, q, sp.oo) for value in radial]
    check("F4 radial eigenvalues have the claimed scalar limit",
          all(sp.simplify(value - (sp.Rational(15, 2)/cV-h)) == 0
              for value in limits))
    section = sp.diag(*([rho-R/2 for rho in ricci] + [-R/2]*10))
    check("independent section control forces scalar curvature zero",
          sp.expand(polynomial(section)[4, 4] - 9*R**2) == 0)

    scalar, c = sp.symbols("scalar c", real=True, nonzero=True)
    h_euler = 6/c
    trace_B = -6*scalar - 14*h_euler
    euler_eigenvalue = -1/c - scalar/2 - h_euler
    check("N1 Euler eigenvalue equals trace(B)/12 at the forced scale",
          sp.expand(euler_eigenvalue - trace_B/12) == 0)
    euler_coefficient = (3300*(tau/12)**2
                         - 380*tau*(tau/12) + 9*tau**2)
    check("spectral polynomial on Euler field gives tau^2/4",
          sp.expand(euler_coefficient - tau**2/4) == 0)
    kappa = sp.symbols("kappa", real=True, nonzero=True)
    scale = sp.Rational(150000, 11)/kappa**2
    check("baseline and scale cancel exactly for either kappa sign",
          sp.simplify(sp.Rational(15, 2)/scale
                      - sp.Rational(11, 20000)*kappa**2) == 0
          and scale.is_positive is True)

    # The spectral condition alone allows nonzero B. The Euler geometry is essential.
    check("nonzero spectral control prevents an unconditional B=0 inference",
          B.subs(z, 1) != sp.zeros(14)
          and polynomial(B.subs(z, 1)) == sp.zeros(14))
    # Nilpotence need not imply zero for a self-adjoint indefinite endomorphism.
    eta = sp.diag(1, -1)
    nilpotent = sp.Matrix([[1, 1], [-1, -1]])
    check("indefinite control keeps positivity explicit",
          nilpotent != sp.zeros(2) and nilpotent**2 == sp.zeros(2)
          and nilpotent.T*eta == eta*nilpotent)
    print(f"{len(CHECKS)} exact checks passed.")
    print("AK subsequently classifies both Lambda3 branches on all positive Z1 metrics.")


if __name__ == "__main__":
    main()
