"""Algebraic data of alpha_cand = (79 + 48 r2 + 42 r3 + 3 r6)/31104 and the
general record-count formula f(dQ, dq) of the Letter's construction.
Run: python3 algebraic.py"""
import sympy as sp
r2, r3, r6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
a_c = (79 + 48*r2 + 42*r3 + 3*r6) / 31104
t = sp.symbols('t')
mp_ = sp.minimal_polynomial(a_c, t)
print("minimal polynomial of alpha_cand:", sp.factor(mp_))
mpi = sp.minimal_polynomial(1/a_c, t)
print("minimal polynomial of 1/alpha_cand:", mpi)
print("Galois conjugates of 1/alpha_cand:")
for s2 in (1, -1):
    for s3 in (1, -1):
        v = 31104/(79 + 48*s2*r2 + 42*s3*r3 + 3*s2*s3*r6)
        print(f"   signs (r2,r3)=({s2:+d},{s3:+d}):", sp.N(v, 20))
# general formula: dQ = read-end (center) dim, dq = far-end dim
def f(dQ, dq, comm=True):
    a, b = dQ**2 - 1, dq**2 - 1
    n_odd = sp.binomial(a, 2) + sp.binomial(b, 2) + 2*a*b
    n_em = 2*a*b
    n_es = a*(a+1)//2 + b*(b+1)//2
    n_c = sp.binomial(b, 2) if (comm and dq % 2 == 0) else 0
    N = (dQ**2 - 1)*(dq**2 - 1)*(dQ*dq)**4
    num = n_odd + n_em*sp.sqrt(dq) + n_es*sp.sqrt(dQ) + n_c*sp.sqrt(dQ*dq)
    return num, N
num, N = f(3, 2)
assert sp.simplify(num/N - a_c) == 0
print("general formula reproduces alpha_cand at (dQ,dq)=(3,2)")
print("same construction at other (dQ,dq) [records, 1/value]:")
for dQ in (2, 3, 4, 5):
    for dq in (2, 3, 4):
        num, N = f(dQ, dq)
        print(f"   ({dQ},{dq}): numerator {sp.radsimp(num)}, N={N}, 1/value = {sp.N(N/num, 12)}")
# deviations from measurements (value, standard uncertainty) of 1/alpha
import mpmath as mp
mp.mp.dps = 40
inv = mp.mpf(31104) / (79 + 48*mp.sqrt(2) + 42*mp.sqrt(3) + 3*mp.sqrt(6))
meas = {"CODATA 2022": ("137.035999177", "0.000000021"), "e g-2 (Fan 2023)": ("137.035999166", "0.000000015"),
        "Rb (Morel 2020)": ("137.035999206", "0.000000011"), "Cs (Parker 2018)": ("137.035999046", "0.000000027")}
for k, (v, u) in meas.items():
    print(f"   {k:18s} (cand - meas)/u = {mp.nstr((inv - mp.mpf(v))/mp.mpf(u), 4)}")
print("   Rb vs g-2 tension:", mp.nstr((mp.mpf('206') - 166)/mp.sqrt(11**2 + 15**2), 4), "sigma")
