"""TH-020..TH-024: the Lynchpin / "where n bubbles meet" constructions, under explicit reconstructions.

Reconstruction L-T (formalization/howard/geometry.md): the six "pentangular areas" are planar walls
lying on the cone over the 1-skeleton of a regular tetrahedron centred at O (Taylor's T-cone).
Everything below is exact (sympy) unless a float is printed with stated digits.
"""
import itertools
import numpy as np
import sympy as sp

print("sympy", sp.__version__, "numpy", np.__version__)
V = [sp.Matrix(v) for v in [(1, 1, 1), (1, -1, -1), (-1, 1, -1), (-1, -1, 1)]]
dot = lambda p, q: (p.T * q)[0]
thetaT = sp.acos(sp.Rational(-1, 3))
deg = lambda t: sp.N(t * 180 / sp.pi, 12)

print("\n== L-T.1  rays O->V_i ==")
angles = {sp.simplify(sp.acos(dot(V[i], V[j]) / sp.sqrt(dot(V[i], V[i]) * dot(V[j], V[j]))))
          for i, j in itertools.combinations(range(4), 2)}
print("pairwise ray angles:", angles, "=", [deg(t) for t in angles], "deg")
top = V[0]
perp = [V[j] - dot(V[j], top) / dot(top, top) * top for j in (1, 2, 3)]
ang120 = {sp.simplify(sp.acos(dot(p, q) / sp.sqrt(dot(p, p) * dot(q, q)))) for p, q in itertools.combinations(perp, 2)}
print("'one on the top and three ... 120 degrees apart': azimuths of the other three rays differ by",
      [deg(t) for t in ang120], "deg")

print("\n== L-T.2  the six wall planes span(V_i, V_j) are exactly the mirror planes of T_d ==")
vset = {tuple(v) for v in V}
for i, j in itertools.combinations(range(4), 2):
    n = V[i].cross(V[j])
    refl = lambda x: x - 2 * dot(x, n) / dot(n, n) * n
    image = {tuple(refl(v)) for v in V}
    print(f"  plane(V{i+1},V{j+1}) normal {tuple(n)}: reflection permutes vertices: {image == vset}")

print("\n== L-T.3  dihedral angles between the three walls meeting along a ray ==")
for i in range(4):
    others = [j for j in range(4) if j != i]
    w = [V[j] - dot(V[j], V[i]) / dot(V[i], V[i]) * V[i] for j in others]
    a = sorted({deg(sp.acos(sp.simplify(dot(p, q) / sp.sqrt(dot(p, p) * dot(q, q))))) for p, q in itertools.combinations(w, 2)})
    print(f"  along ray V{i+1}: {a} deg")

print("\n== L-T.4  complement of the six walls = four open trihedral cones ('four pyramidal inner spaces') ==")
print("  sum of V_i =", tuple(sum(V, sp.zeros(3, 1))), "(positive spanning set => the four face-cones cover R^3)")
rng = np.random.default_rng(0)
Vn = np.array([[float(c) for c in v] for v in V])
counts = np.zeros(5, dtype=int)
for d in rng.normal(size=(200000, 3)):
    k = 0
    for l in range(4):
        M = np.array([Vn[m] for m in range(4) if m != l]).T
        lam = np.linalg.solve(M, d)
        k += bool(np.all(lam > 0))
    counts[k] += 1
print("  200000 random directions (seed 0), number of open face-cones containing each:", dict(enumerate(counts)))

print("\n== L-T.5  regular pentagons cannot be the walls ==")
pent = 3 * sp.pi / 5
print("  regular pentagon interior angle:", deg(pent), "deg; wall apex angle arccos(-1/3):", deg(thetaT), "deg")
print("  mismatch arccos(-1/3) - 108 deg =", deg(thetaT - pent), "deg")

def sym_equilateral_pentagon(theta):
    """Unit-edge pentagon O,A,B,B',A' with mirror symmetry and apex angle theta at O (convex branch)."""
    A = sp.Matrix([sp.cos(theta / 2), sp.sin(theta / 2)])
    yB = sp.Rational(1, 2)
    xB = A[0] + sp.sqrt(1 - (yB - A[1])**2)
    B = sp.Matrix([xB, yB]); O = sp.Matrix([0, 0])
    ang = lambda P, Q, R: sp.acos(sp.simplify(dot(Q - P, R - P) / sp.sqrt(dot(Q - P, Q - P) * dot(R - P, R - P))))
    return [deg(theta), deg(ang(A, O, B)), deg(ang(B, A, sp.Matrix([xB, -yB])))]
print("  symmetric equilateral pentagon, apex 108 deg  -> angles (O, A, B):", sym_equilateral_pentagon(pent))
print("  symmetric equilateral pentagon, apex arccos(-1/3) -> angles (O, A, B):", sym_equilateral_pentagon(thetaT))

print("\n== L-T.6  an explicit well-posed 'six pentangular areas' object ==")
h, g = 3, 2   # truncate the T-cone by K = {x : V_k . x <= h for all k} ∩ [-g, g]^3 ; pentagons iff h/3 < g < h
s, t = sp.symbols('s t', nonnegative=True)
P = sp.solve([3 * s - t - h, s + t - g], [s, t])
pts2 = [(0, 0), (sp.Rational(h, 3), 0), (P[s], P[t]), (P[t], P[s]), (0, sp.Rational(h, 3))]
emb = lambda st: st[0] * V[0] + st[1] * V[1]
pts3 = [emb(q) for q in pts2]
print("  wall in plane(V1,V2): vertices", [tuple(p) for p in pts3])
L = [sp.sqrt(dot(pts3[(k + 1) % 5] - pts3[k], pts3[(k + 1) % 5] - pts3[k])) for k in range(5)]
print("  edge lengths:", [sp.nsimplify(l) for l in L])
ang = []
for k in range(5):
    p0, p1, p2 = pts3[k - 1], pts3[k], pts3[(k + 1) % 5]
    u, v = p0 - p1, p2 - p1
    ang.append(deg(sp.acos(sp.simplify(dot(u, v) / sp.sqrt(dot(u, u) * dot(v, v))))))
print("  interior angles (deg):", ang, " sum:", sp.N(sum(ang), 12))

print("\n== L-T.7  sphere-truncated walls ('circular internal lynchpin', D1002750 title) ==")
r = sp.symbols('r', positive=True)
area = 6 * thetaT / 2 * r**2
print("  total wall area inside a ball of radius r: 3*arccos(-1/3)*r^2 =", sp.N(area / r**2, 12), "* r^2")
print("  area density relative to a flat disc (pi r^2):", sp.N(3 * thetaT / sp.pi, 12), " (Y-cone: 3/2)")

print("\n== L-T.8  the 108 deg vs Plateau mismatch for a 12-neighbour (dodecahedral) cell ==")
print("  dodecahedron face angle 108 deg vs Plateau vertex angle", deg(thetaT))
print("  dodecahedron dihedral angle arccos(-1/sqrt5) =", deg(sp.acos(-1 / sp.sqrt(5))), "deg vs Plateau 120 deg")

print("\n== L-B.1  'negative space where four bubbles meet' (4 tangent unit balls at tetrahedron vertices) ==")
Om = 3 * sp.acos(sp.Rational(1, 3)) - sp.pi
print("  vertex solid angle 3*arccos(1/3) - pi =", sp.N(Om, 12), "; equals arccos(23/27):",
      sp.N(Om - sp.acos(sp.Rational(23, 27)), 15) == 0 or abs(sp.N(Om - sp.acos(sp.Rational(23, 27)), 30)) < 1e-25)
Vtet = (2)**3 / (6 * sp.sqrt(2))
inter = sp.simplify(Vtet - 4 * Om / (4 * sp.pi) * sp.Rational(4, 3) * sp.pi)
print("  tetrahedron (edge 2) volume:", sp.nsimplify(Vtet), "=", sp.N(Vtet, 12))
print("  interstice volume (r=1):", inter, "=", sp.N(inter, 12))
print("  fraction of tetrahedron covered by the balls sqrt2*Omega =", sp.N(sp.sqrt(2) * Om, 12), "(Rogers 1958)")

print("\n== L-B.2  eight tangent unit balls at cube vertices ==")
print("  interstice inside cube (edge 2):", 8 - sp.Rational(4, 3) * sp.pi, "=", sp.N(8 - sp.Rational(4, 3) * sp.pi, 12),
      "; central gap radius sqrt3-1 =", sp.N(sp.sqrt(3) - 1, 12))

print("\n== L-B.3  vertex/face/edge counts of candidate n-ball clusters vs stated 'pole' counts ==")
table = [("tetrahedron", 4, 6, 4), ("octahedron", 6, 12, 8), ("cube", 8, 12, 6), ("square antiprism", 8, 16, 10),
         ("icosahedron", 12, 30, 20), ("cuboctahedron", 12, 24, 14), ("truncated tetrahedron", 12, 18, 8),
         ("truncated octahedron", 24, 36, 14), ("truncated cube", 24, 36, 14), ("rhombicuboctahedron", 24, 48, 26),
         ("snub cube", 24, 60, 38)]
stated = {4: "4 contractive poles", 6: "30 poles; 6 unaccounted electrical", 8: "8 contractive, 6 magnetic",
          12: "4 unaccounted electrical", 24: "all accounted (0 unaccounted)"}
for name, Vv, E, F in table:
    assert Vv - E + F == 2
    print(f"  {name:22s} V={Vv:2d} E={E:2d} F={F:2d}  V-F={Vv-F:+3d}  | Howard (n={Vv}): {stated[Vv]}")

print("\n== G1  closed sphere-type polyhedra with all faces pentagons need F >= 12 (Lean: pentagon_faces_ge_12) ==")
for F in range(1, 13):
    E = sp.Rational(5 * F, 2); Vv = 2 + E - F
    print(f"  F={F:2d}: E={E}, V={Vv}, 2E>=3V? {2*E >= 3*Vv}") if F in (6, 11, 12) else None
