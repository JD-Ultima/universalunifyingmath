"""TH-005..TH-012: Howard's numerical statements, tested independently.

All exact statements use sympy; floating statements state their precision.
Sources: JRE #2152 transcript (timestamps mm:ss from the transcript file; see sources/howard.md).
"""
import math
from fractions import Fraction as Fr
import sympy as sp
import mpmath as mp

print("sympy", sp.__version__, "mpmath", mp.__version__)
x = sp.symbols('x', positive=True)
r2 = sp.sqrt(2)

# ---------------------------------------------------------------- TH-005: the sqrt(2) "loop"
print("\n== TH-005 sqrt(2) loop (JRE #2152, 124:37-125:38) ==")
print("exact (sqrt2)^3            =", sp.simplify(r2**3), "  == 2*sqrt2 :", sp.simplify(r2**3 - 2 * r2) == 0)
print("exact (sqrt2)^3 / 2        =", sp.simplify(r2**3 / 2))
print("fixed points of x^3/2      :", sp.solve(sp.Eq(x**3 / 2, x), x), "(positive); 0 and -sqrt2 also")
c = sp.symbols('c', positive=True)
print("general: sqrt(c) fixed by x^3/c :", sp.simplify((sp.sqrt(c))**3 / c - sp.sqrt(c)) == 0)
print("derivative of x^3/2 at sqrt2 :", sp.diff(x**3 / 2, x).subs(x, r2), "(> 1: repelling)")
print("x^3 = 2x  <=>  x^2 = 2 (x != 0):", sp.solve(sp.Eq(x**3, 2 * x), x))
print("quoted digits  sqrt2      : 1.414213562373095   true:", mp.nstr(mp.sqrt(2), 20))
print("quoted digits  (sqrt2)^3  : 2.82842712174619    true:", mp.nstr(2 * mp.sqrt(2), 20),
      " (transcript differs in the 9th decimal: ...712[1]74... vs ...712[4]74...)")
with mp.workdps(40):
    print("float64 sqrt(2) =", repr(math.sqrt(2)), " error vs true:",
          mp.nstr(mp.mpf(math.sqrt(2)) - mp.sqrt(2), 5))

def escape(step_fn, x0, target, tol=1e-3, nmax=400):
    """Step at which |x_n - target| > tol; 'stuck' if the rounded map fixes x exactly."""
    xn = x0
    for n in range(1, nmax + 1):
        xp, xn = xn, step_fn(xn)
        if abs(xn - target) > tol:
            return n
        if xn == xp:
            return f"stuck at step {n} (rounded map fixes the rounded sqrt2 exactly)"
    return None

n64 = escape(lambda t: t**3 / 2.0, math.sqrt(2), math.sqrt(2))
print(f"float64 iteration x -> x^3/2 from fl(sqrt2): leaves |x - sqrt2| <= 1e-3 at step {n64}")
for dps in (15, 30, 50, 100):
    with mp.workdps(dps):
        s2 = mp.sqrt(2)
        n = escape(lambda t: t**3 / 2, s2, s2)
        eps = abs(mp.sqrt(2) - s2) if dps < 120 else None
    with mp.workdps(dps + 60):
        err0 = abs(mp.mpf(s2) - mp.sqrt(2))
        pred = mp.log(mp.mpf('1e-3') / err0) / mp.log(3) if err0 != 0 else None
    print(f"mpmath dps={dps:3d}: escape step {n}; first-order prediction log(1e-3/eps0)/log 3 = "
          f"{mp.nstr(pred, 4) if pred else 'n/a'} (eps0 = {mp.nstr(err0, 3)})")
tally = {"escape": 0, "stuck": 0}
for dps in range(10, 121):
    with mp.workdps(dps):
        s2 = mp.sqrt(2)
        n = escape(lambda t: t**3 / 2, s2, s2)
    tally["stuck" if isinstance(n, str) else "escape"] += 1
print(f"mpmath dps 10..120: escapes in {tally['escape']} precisions, exact floating fixed point in {tally['stuck']}")
print("calculator loop alternative (cube, then divide by 2, start 1.414 exactly as typed):")
t = Fr(1414, 1000)
seq = []
for _ in range(4):
    t = t**3 / 2
    seq.append(float(t))
print("   ", [f"{v:.6f}" for v in seq])

# ---------------------------------------------------------------- TH-006: 'x + x > x^2 for most numbers'
print("\n== TH-006 'square root of a number added to itself ... more than that number' (39:50) ==")
n_ = sp.symbols('n', positive=True)
print("reading R1: 2*sqrt(n) > n   <=> n in", sp.solve_univariate_inequality(2 * sp.sqrt(n_) > n_, n_, relational=False))
print("reading R2: x + x > x^2     <=> x in", sp.solve_univariate_inequality(2 * x > x**2, x, relational=False))
for N in (10, 100, 10**6):
    k1 = sum(1 for m in range(1, N + 1) if 4 * m > m * m)          # 2 sqrt m > m  <=> 4m > m^2
    k2 = sum(1 for m in range(1, N + 1) if 2 * m > m * m)
    print(f"   integers 1..{N}: R1 holds for {k1}, R2 holds for {k2}")

# ---------------------------------------------------------------- TH-008: speed of light figures
print("\n== TH-008 speed-of-light figures (33:25 and 126:59) ==")
c_SI = 299792458                       # m/s, exact by SI definition (BIPM 2019)
mile = Fr(1609344, 1000)               # international mile, exact (1959)
c_mi = Fr(c_SI) / mile
print("c (m/s, exact)            :", c_SI, "  quoted '299,752,400 something'")
print("   difference             :", c_SI - 299752400, "m/s; relative", f"{(c_SI-299752400)/c_SI:.3e}")
print("c (mi/s, exact rational)  :", c_mi, "=", f"{float(c_mi):.6f}", "  quoted 186,624")
print("   difference             :", f"{186624 - float(c_mi):.6f}", "mi/s; relative", f"{(186624-float(c_mi))/float(c_mi):.3e}")
print("186624 = 432^2 ?", 432**2 == 186624, "; factorization", sp.factorint(186624))
print("'300 million km/second squared': c in km/s =", c_SI / 1000, "(a speed, not an acceleration)")

# ---------------------------------------------------------------- TH-009: 'element tones' in Hz
print("\n== TH-009 element 'tones' (10:22-10:28): H 40.5, C 81, Si 162, Co 324 Hz ==")
tones = [Fr(81, 2), 81, 162, 324, 648]
print("factorizations:", [(str(t), sp.factorrat(sp.Rational(t.numerator, t.denominator))) for t in map(Fr, tones)])
print("successive ratios:", [str(Fr(tones[i + 1]) / Fr(tones[i])) for i in range(len(tones) - 1)])
print("Pythagorean E over C=256 Hz: 256*81/64 =", Fr(256) * Fr(81, 64), "; /8 =", Fr(256) * Fr(81, 64) / 8)
E1_ET = 440 * 2 ** (-41 / 12)         # E1 (MIDI 28) is 41 semitones below A4 (MIDI 69)
print(f"12-TET (A4=440) E1 = {E1_ET:.4f} Hz; 40.5 Hz is {1200*math.log2(40.5/E1_ET):+.2f} cents from it")
print(f"Pythagorean third 81/64 vs 12-TET third: {1200*math.log2(81/64) - 400:+.2f} cents")
try:
    import mendeleev
    from mendeleev import element
    print("mendeleev", mendeleev.__version__, "(data bundled with the package)")
    for group in (["H", "C", "Si", "Co", "Rh"], ["Li", "Na", "K"]):
        els = [element(s) for s in group]
        Z = [e.atomic_number for e in els]
        W = [e.atomic_weight for e in els]
        I = [e.ionenergies[1] for e in els]
        rat = lambda v: [round(v[i + 1] / v[i], 4) for i in range(len(v) - 1)]
        print(f"  {group}: Z={Z} ratios {rat(Z)}")
        print(f"  {'':{len(str(group))}}  atomic weight={W} ratios {rat(W)}")
        print(f"  {'':{len(str(group))}}  1st ionization energy (eV)={I} ratios {rat(I)}")
except Exception as ex:  # pragma: no cover
    print("mendeleev unavailable:", ex)

# ---------------------------------------------------------------- TH-010: "64 circles" of the flower of life
print("\n== TH-010 'these 64 circles overlapping' (29:01) ==")
hexn = [1 + 3 * k * (k + 1) for k in range(0, 8)]
print("centered hexagonal numbers 1+3k(k+1):", hexn, "; 64 among them:", 64 in hexn)
# circles centred on triangular-lattice points within Euclidean radius R (unit spacing)
pts = [(i + j / 2, j * math.sqrt(3) / 2) for i in range(-12, 13) for j in range(-12, 13)]
d = sorted(set(round(math.hypot(px, py), 9) for px, py in pts))
counts = sorted(set(sum(1 for px, py in pts if math.hypot(px, py) <= R + 1e-9) for R in d[:40]))
print("circle counts achievable by a disc truncation of the lattice (first values):", counts[:20],
      "; 64 achievable:", 64 in counts)
print("64 = 4^3 = 2^6:", 4**3, 2**6)

# ---------------------------------------------------------------- P-TH-1: 3-smoothness of the quoted constants
print("\n== P-TH-1 3-smooth (2^a 3^b) constants ==")
smooth = lambda q: set(sp.factorrat(sp.Rational(q)).keys()) <= {2, 3}
for q in ["81/2", 81, 162, 324, 648, 186624, 432, 299752400, 300000000]:
    print(f"  {q}: 3-smooth {smooth(q)}  {sp.factorrat(sp.Rational(q))}")
print("  cross-framework remark: Chain Theory readout dimension 31104 =", sp.factorint(31104),
      "(3-smooth because built from qubit/qutrit dimensions 2 and 3; no structural link claimed)")
