#!/usr/bin/env python3
"""Round 3: the mass cells of the compendium -- derivation or selection?

Usage: python3 -I masscells.py <compendium_dir> [--pdg-mcd mass_width_2024.mcd] [--tsv OUT.tsv]

Object: m_p/m_e ~ (2871/67)^2 (+1.84e-5) and m_mu/m_e ~ (963/67)^2, m_tau/m_e ~ (3969/67)^2 (scorecard SC-35/36/38).
Blocks (each prints its own header):
  ANCHORS   every curated string (masscells_data.py) is checked verbatim against its statement unit
  RECORDS   the (2,3) record census re-derived combinatorially: 172 records, closed count 79 = 48+28+3, Born count 319
  WEBS      the electron web parsed from the 152.1 figure; muon, tau and proton webs built from the printed
            descriptions; measure, strike, rank, hub in-legs, dressing (198.3) checked against the printed values
  PATHS     the measuring read (20.1, 186.1) by breadth-first search from the entered site: in-legs passed / unpassed
  READINGS  the uniform-reading grid: records typing x strike rule x dressing rule (the axes 20.1, 21.1, 32.1, 186.1,
            155.1, 198.3 name), validated against the six printed tuples; which readings meet the cells / the ratios
  FACES     the 155.1 faces (dock, hub, mixed) per body and the ratio at each
  MASSLAW   the agreeing cells under every printed mass law (S, S^2, S^2 + closure of 195.1, T1(1+G2^2) of 176.2)
  SPACES    distribution of predicted ratios and look-elsewhere estimates: printed cells; per-body choice over the
            reading grid; the round-2 declared space; the narrow family 319*3^u/(67*3^v); any-integer baselines
  LEPTON    the face-free content (tau/mu) and the factor nine
  CROSS     the same segment rule at the quark sinks of 181.5
  DIGITS    digits put in against digits matched
  SWEEP     mechanical sweep of integers in proton first-power contexts, classified (audit of the printed space)
Exact arithmetic with sympy; 30-digit mpmath for logs. References: CODATA 2022 from scipy.constants 1.18.1 (asserted
against the transcription); PDG 2024 lepton and quark masses from the MC table mass_width_2024.mcd (checked when
--pdg-mcd is given). Read-only: nothing in the source is executed.
"""
import collections, itertools, math, pathlib, re, sys
import sympy as sp
import mpmath as mp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))   # our own directory only
import scorecard_extract as X
import masscells_data as D

mp.mp.dps = 30
R = sp.Rational

# ------------------------------------------------------------------ references
import scipy, scipy.constants as sc
assert scipy.__version__.startswith("1.18"), scipy.__version__
PC = sc.physical_constants
CODATA = {"proton-electron mass ratio": (1836.152673426, 3.2e-8), "muon-electron mass ratio": (206.7682827, 4.6e-6),
          "tau-electron mass ratio": (3477.23, 0.23)}
for k, (v, e) in CODATA.items():
    assert abs(PC[k][0] - v) <= 1e-12 * v and abs(PC[k][2] - e) <= 1e-6 * e, (k, PC[k])
MP, sMP = CODATA["proton-electron mass ratio"]
MMU, sMMU = CODATA["muon-electron mass ratio"]
PDG = {11: (5.1099895000e-4, 1.5e-13), 13: (1.056583755e-1, 2.3e-9), 15: (1.77693, 9.0e-5),   # GeV, mass_width_2024.mcd
       2: (2.16e-3, 7e-5), 1: (4.70e-3, 7e-5), 3: (9.35e-2, 8e-4), 4: (1.273, 5e-3), 5: (4.183, 7e-3), 6: (1.7257e2, 2.9e-1)}
MTAU = PDG[15][0] / PDG[11][0]
sMTAU = MTAU * math.hypot(PDG[15][1] / PDG[15][0], PDG[11][1] / PDG[11][0])
TAUMU = PDG[15][0] / PDG[13][0]
sTAUMU = TAUMU * math.hypot(PDG[15][1] / PDG[15][0], PDG[13][1] / PDG[13][0])
_m = [PDG[11][0], PDG[13][0], PDG[15][0]]
KOIDE = sum(_m) / sum(math.sqrt(x) for x in _m) ** 2
_kq = lambda m: sum(m) / sum(math.sqrt(x) for x in m) ** 2
sKOIDE = math.sqrt(sum(((_kq([x + (PDG[i][1] if j == k else 0) for k, (x, i) in enumerate(zip(_m, (11, 13, 15)))]) - KOIDE) / 1.0) ** 2
                       for j in range(3)))


def check_mcd(path):
    rows = {}
    for line in open(path, encoding="ascii", errors="replace"):
        if line.startswith("*") or len(line) < 60 or not line[33:51].strip(): continue
        for i in [int(x) for x in line[:32].split()]:
            rows[i] = (float(line[33:51]), float(line[52:60]))
    for k, (m, e) in PDG.items():
        assert abs(rows[k][0] - m) <= 1e-9 * m and abs(rows[k][1] - e) <= 1e-9 * e, (k, rows[k])
    return len(rows)


def ff(x, n=10):
    return sp.N(x, n) if isinstance(x, sp.Basic) else x


def rel(v, t):
    return float(v) / t - 1


def field(t):          # (a, b, c, d) rational strings -> a + b sqrt2 + c sqrt3 + d sqrt6
    a, b, c, d = (R(x) for x in t)
    return a + b * sp.sqrt(2) + c * sp.sqrt(3) + d * sp.sqrt(6)


# ------------------------------------------------------------------ ANCHORS
def anchors(src):
    U = {u["num"]: u for u in X.units(src / "tex")}
    UN = {k: X.norm(u["body"] + " " + u["rest"]) for k, u in U.items()}
    groups = {"trace": D.TRACE, "closures": [a for _, a in D.CLOSURES.values()] + [D.G2[1]],
              "printed cells": [c[2] for L in (D.PROTON_CELLS, D.ELECTRON_UNITS, D.MUON_CELLS, D.TAU_CELLS) for c in L]}
    print("[ANCHORS] verbatim check of the curated strings (whitespace-normalised) against their statement units")
    tot = miss = 0
    for g, L in groups.items():
        bad = [(n, s) for n, s in L if X.norm(s) not in UN.get(n, "")]
        tot += len(L); miss += len(bad)
        print(f"  {g:14s} {len(L):3d} strings, {len(bad)} missing" + "".join(f"\n    MISSING {n}: {s[:90]}" for n, s in bad))
    fb = [(f, s) for f, s in D.FILE_ANCHORS if X.norm(s) not in X.norm((src / "tex" / f).read_text(encoding="utf-8"))]
    print(f"  file-level     {len(D.FILE_ANCHORS):3d} strings, {len(fb)} missing")
    grades = {k: U[k]["grade"] for k in sorted({n for n, _ in D.TRACE}, key=lambda s: [int(x) if x.isdigit() else x for x in re.split(r'[.]', s)])}
    print("  grades of the statements cited: " + ", ".join(f"{k} {v[0]}" for k, v in grades.items()))
    print(f"  total {tot + len(D.FILE_ANCHORS)} strings, {miss + len(fb)} missing\n")
    assert miss == 0 and not fb
    return U


# ------------------------------------------------------------------ RECORDS: the (2,3) link census
def records():
    print("[RECORDS] the (2,3) record census re-derived (letters: 8 traceless at each qutrit slot Q1,Q2; 3 at each qubit slot q1,q2;")
    print("          two-slot products; exchange E: Q1<->Q2, q1<->q2; records p +/- E(p), E-fixed products, 3 commutator records)")
    slots = {"Q1": 8, "q1": 3, "Q2": 8, "q2": 3}
    swap = {"Q1": "Q2", "Q2": "Q1", "q1": "q2", "q2": "q1"}
    prods = []
    for s1, s2 in itertools.combinations(["Q1", "q1", "Q2", "q2"], 2):
        for a in range(slots[s1]):
            for b in range(slots[s2]):
                prods.append(frozenset({(s1, a), (s2, b)}))
    E = lambda p: frozenset((swap[s], a) for s, a in p)
    cls = lambda p: "same" if {s[0] for s, _ in p} in ({"Q"}, {"q"}) else "mixed"   # both qutrit or both qubit slots = same-end
    kind = lambda p: "QQ" if {s[0] for s, _ in p} == {"Q"} else ("qq" if {s[0] for s, _ in p} == {"q"} else "Qq")
    fixed = [p for p in prods if E(p) == p]
    pairs = {frozenset({p, E(p)}) for p in prods if E(p) != p}
    odd = collections.Counter(kind(next(iter(pr))) for pr in pairs)
    even = collections.Counter(cls(next(iter(pr))) for pr in pairs) + collections.Counter(cls(p) for p in fixed)
    n_comm = 3
    k = {"odd": 1, "mixed": 2, "same": 3, "comm": 6}     # the alpha Letter's squared weights (F-CT-2 class function)
    n = {"odd": sum(odd.values()), "mixed": even["mixed"], "same": even["same"], "comm": n_comm}
    born = sum(n[c] * k[c] for c in n)
    print(f"  {len(prods)} products, {len(fixed)} E-fixed, {len(pairs)} pairs; records {sum(n.values())} = odd {n['odd']} + even mixed "
          f"{n['mixed']} + even same-end {n['same']} + commutator {n['comm']}")
    print(f"  closed (odd, k = 1) records {n['odd']} = mixed {odd['Qq']} + two qutrit slots {odd['QQ']} + two qubit slots {odd['qq']} "
          f"(176.2 prints 79 = 48 + 28 + 3)")
    print(f"  Born count with the Letter's weights k = 1, 2, 3, 6: {n['odd']}*1 + {n['mixed']}*2 + {n['same']}*3 + {n['comm']}*6 = {born} "
          f"(176.2: 319, 'the cap' of the (2,3) genre)")
    # dependence of 319 on the weight law: the alpha audit's F-CT-4 squared-factor space (k in {1,2,3,6} per product class;
    # commutator factor {1,2,3,4,6,9}, count {0,3,6}); classes odd 79, even mixed 48, even QQ 36, even qq 6
    bc = collections.Counter(79 * a + 48 * b + 36 * c + 6 * d + nc * kc for a, b, c, d in itertools.product((1, 2, 3, 6), repeat=4)
                             for kc in (1, 2, 3, 4, 6, 9) for nc in (0, 3, 6))
    print(f"  over the alpha audit's 4^4 x 6 x 3 = {sum(bc.values())} weight rules the Born count takes {len(bc)} values "
          f"({min(bc)} to {max(bc)}); 319 occurs at {bc[319]} rules; so 319 is fixed once the alpha weight law is fixed\n")
    assert (n["odd"], n["mixed"], n["same"], n["comm"], born) == (79, 48, 42, 3, 319) and (odd["Qq"], odd["QQ"], odd["qq"]) == (48, 28, 3)
    return born


# ------------------------------------------------------------------ WEBS
class Web:
    def __init__(self, name, dims, links, struck, entered, capped_entry=False):
        self.name, self.dims, self.links, self.struck, self.entered, self.capped_entry = name, dims, links, set(struck), entered, capped_entry
        for a, b in links: assert (a, b) in links
    def measure(self): return sum(self.dims[a] * self.dims[b] for a, b in self.links)
    def strike(self): return sum(self.dims[a] * self.dims[b] for a, b in self.links if (a, b) in self.struck)
    def rank(self): return self.measure() - self.strike()
    def hubs(self): return [s for s, d in self.dims.items() if d == 4]
    def inlegs(self):     # (3,4) links: qutrit -> ququart hub
        return [(a, b) if self.dims[b] == 4 else (b, a) for a, b in self.links if {self.dims[a], self.dims[b]} == {3, 4}]
    def dressing(self):   # 198.3: product over the distinct hub in-legs of the far-end dimension
        return math.prod(self.dims[q] for q, h in self.inlegs())
    def adj(self):
        A = collections.defaultdict(set)
        for a, b in self.links: A[a].add(b); A[b].add(a)
        return A


def parse_figure(unit):
    """Sites and links of the picture environment in a statement unit: circles with a dimension and a label, qbezier
    links with a genre label; each bezier end is matched to the nearest circle."""
    t = unit["rest"]
    pic = t[t.index("\\begin{picture}"): t.index("\\end{picture}")]
    sites = {}
    for m in re.finditer(r"\\put\(([\d.]+),([\d.]+)\)\{\\circle\{14\}\}\s*\\put\(\1,\2\)\{\\makebox\(0,0\)\{\\footnotesize \$(\d)\$\}\}\s*"
                         r"\\put\([\d.]+,[\d.]+\)\{\\makebox\(0,0\)\[b\]\{\\scriptsize \$([^$]+)\$\}\}", pic):
        sites[m.group(4).replace("_", "")] = (float(m.group(1)), float(m.group(2)), int(m.group(3)))
    links, genres = [], []
    for m in re.finditer(r"\\qbezier\(([\d.]+),([\d.]+)\)\([\d.]+,[\d.]+\)\(([\d.]+),([\d.]+)\)\s*\\put\([\d.]+,[\d.]+\)\{\\makebox\(0,0\)"
                         r"\{\\scriptsize \$\((\d),(\d)\)\$\}\}", pic):
        ends = []
        for x, y in ((float(m.group(1)), float(m.group(2))), (float(m.group(3)), float(m.group(4)))):
            ends.append(min(sites, key=lambda s: (sites[s][0] - x) ** 2 + (sites[s][1] - y) ** 2))
        links.append(tuple(ends)); genres.append((int(m.group(5)), int(m.group(6))))
    return {s: v[2] for s, v in sites.items()}, links, genres


def assembly(prefix, core_site, dims, links):
    h = "H" if prefix == "a" else "H" + {"b": "2", "c": "3"}[prefix]
    dims[h] = 4; links.append((core_site, h))
    for i in (1, 2, 3):
        dims[f"{prefix}{i}"] = 2; links.append((h, f"{prefix}{i}"))
    links.append((f"{prefix}1", f"{prefix}2"))           # the qubit pair of the assembly
    return h


def build_webs(U):
    print("[WEBS] electron parsed from the 152.1 figure; muon, tau, proton built from the printed descriptions (masscells_data.WEB_NOTES)")
    dims, links, genres = parse_figure(U["152.1"])
    ok = all(tuple(sorted((dims[a], dims[b]))) == g for (a, b), g in zip(links, genres))
    struck_e = [(a, b) for a, b in links if {dims[a], dims[b]} == {2, 3}]
    W = {"e": Web("electron", dims, links, struck_e, "w1")}
    print(f"  152.1 figure: {len(dims)} sites {sorted(dims.items())}, {len(links)} links, every printed genre = (dim, dim) of its ends: {ok}")
    core = {"A": 3, "B": 3, "C": 3}
    tri = [("A", "B"), ("B", "C"), ("C", "A")]
    d, L = dict(core), list(tri)
    assembly("a", "A", d, L); assembly("b", "B", d, L)
    d["w2"] = 2; L.append(("w2", "C")); d["w1"] = 2; L.append(("w1", "b3"))
    W["mu"] = Web("muon", d, L, [("w2", "C"), ("w1", "b3")], "w2")
    d, L = dict(core), list(tri)
    assembly("a", "A", d, L); assembly("b", "B", d, L); assembly("c", "C", d, L)
    d["w1"] = 2; L.append(("w1", "b3")); d["w2"] = 2; L.append(("w2", "c3"))
    W["tau"] = Web("tau", d, L, [("w1", "b3"), ("w2", "c3")], "w1")
    d = {"J": 3, "Q1": 3, "Q2": 3, "Q3": 3, "H": 4, "f1": 2, "f2": 2, "n1": 2, "n2": 2, "n3": 2}
    L = [("J", f"Q{i}") for i in (1, 2, 3)] + [(f"n{i}", f"Q{i}") for i in (1, 2, 3)] + [(f"Q{i}", "H") for i in (1, 2, 3)] + \
        [("H", "f1"), ("H", "f2"), ("f1", "f2")]
    W["p_seed"] = Web("proton seed", d, L, [(f"n{i}", f"Q{i}") for i in (1, 2, 3)], "n1", capped_entry=True)
    d = {"s0": 2, "s1": 2, "s2": 2, "s3": 3, "s4": 3, "s5": 3, "s6": 3, "s7": 4}
    L = [("s0", "s3"), ("s1", "s4"), ("s2", "s5"), ("s3", "s7"), ("s4", "s7"), ("s5", "s7"), ("s6", "s3"), ("s6", "s4"), ("s6", "s5")]
    W["p_sink"] = Web("proton sink of record", d, L, [("s0", "s3"), ("s1", "s4"), ("s2", "s5")], "s0", capped_entry=True)
    printed = {"e": (9, 10, 79, 12, 67, 3), "mu": (13, 15, 117, 10, 107, 9), "tau": (17, 20, 155, 8, 147, 27),
               "p_seed": (10, 12, 101, 18, 83, 27), "p_sink": (8, 9, 81, 18, 63, 27)}
    print(f"  {'web':22s} sites links measure strike rank in-legs dressing(198.3) | printed (sites, links, measure, strike, rank, dressing)")
    for k, w in W.items():
        got = (len(w.dims), len(w.links), w.measure(), w.strike(), w.rank(), w.dressing())
        print(f"  {w.name:22s} {got[0]:5d} {got[1]:5d} {got[2]:7d} {got[3]:6d} {got[4]:4d} {len(w.inlegs()):7d} {got[5]:8d}      | "
              f"{printed[k]} {'equal' if got == printed[k] else 'DIFFERS'}")
        assert got == printed[k], (k, got)
    print("  printed sources: 152.1 (79 = 27+40+12, ranks 67/107/147), 186.1/195.1 (117, 10; 155, 8; links 10/15/20), 195.1 (sites 9/13/17;")
    print("  sink 8 sites, 9 links, measure 81, dock set 18, dressing 27 at the mass-register row), 175.3 (seed 101 = 27+18+36+16+4, ten")
    print("  sites, twelve links, three in-legs). Ranks 27 + 40k: " + str([W[k].rank() for k in ("e", "mu", "tau")]) + "\n")
    return W


# ------------------------------------------------------------------ PATHS
def paths(W):
    print("[PATHS] the measuring read (20.1, 186.1): breadth-first from the entered site to the first hub reached;")
    print("        an in-leg on a geodesic to that hub is 'passed', the other hub in-legs are 'not passed'")
    out = {}
    for k in ("e", "mu", "tau", "p_seed", "p_sink"):
        w = W[k]; A = w.adj()
        dist = {w.entered: 0}; q = collections.deque([w.entered]); par = collections.defaultdict(set)
        while q:
            x = q.popleft()
            for y in A[x]:
                if y not in dist: dist[y] = dist[x] + 1; q.append(y)
                if dist[y] == dist[x] + 1: par[y].add(x)
        step = min(dist[h] for h in w.hubs())
        first = [h for h in w.hubs() if dist[h] == step]
        passed = {(q_, h) for q_, h in w.inlegs() if h in first and q_ in par[h]}
        tie = len(first) > 1
        out[k] = dict(n=len(w.inlegs()), passed=len(passed), tie=tie, step=step)
        print(f"  {w.name:22s} entered at {w.entered}: first hub(s) {first} at step {step}; in-legs {len(w.inlegs())}, on the geodesic "
              f"{len(passed)}{' (a tie: two hubs at once)' if tie else ''}")
    print("  printed: electron 'reaches the hub at the third step with the one in-leg on the geodesic passed'; muon 'reaches two hubs at")
    print("  the third step at once through two in-legs, a tie, both passed' (20.1) / 'with both in-legs unpassed' (186.1); tau 'at the")
    print("  second step through a (2,4) leg with no in-leg on the geodesic'; proton 'at the second step through one in-leg passed and")
    print("  leaves two in-legs not passed' (20.1)\n")
    assert [out[k]["step"] for k in ("e", "mu", "tau", "p_seed")] == [3, 3, 2, 2] and out["mu"]["tie"]
    assert [out[k]["passed"] for k in ("e", "tau", "p_seed", "p_sink")] == [1, 0, 1, 1]
    return out


# ------------------------------------------------------------------ READINGS
BODIES = ("e", "mu", "tau", "p")
TGT = {"mu": (MMU, sMMU), "tau": (MTAU, sMTAU), "p": (MP, sMP)}
DERIVED = {"e": 67, "mu": 963, "tau": 3969, "p": 2871}


def reading_axes(W, P, born):
    web = {"e": W["e"], "mu": W["mu"], "tau": W["tau"], "p": W["p_seed"]}
    reg = {"e": 79, "mu": 79, "tau": 24, "p": born}        # 20.1 registers (closed count of the entered link; cap if capped)
    capped = {"e": False, "mu": False, "tau": False, "p": True}   # 32.1 / 175.4 adopted typing
    n = {b: P["p_seed" if b == "p" else b]["n"] for b in BODIES}
    passed = {b: P["p_seed" if b == "p" else b]["passed"] for b in BODIES}
    tie = {b: P["p_seed" if b == "p" else b]["tie"] for b in BODIES}
    rec = {
        "reach (32.1; 186.1)": ({b: (born if capped[b] else web[b].measure()) for b in BODIES}, capped),
        "register (20.1)": (reg, capped),
        "measure, all uncapped, seed (179.1)": ({b: web[b].measure() for b in BODIES}, {b: False for b in BODIES}),
        "measure, all uncapped, sink (195.1)": ({**{b: web[b].measure() for b in BODIES}, "p": W["p_sink"].measure()}, {b: False for b in BODIES}),
    }
    strike = {"every entry (20.1)": lambda b, cap: web[b].strike() if b != "p" else W["p_seed"].strike(),
              "uncapped entry alone (186.1)": lambda b, cap: 0 if cap[b] else (web[b].strike() if b != "p" else W["p_seed"].strike()),
              "none ('nothing removed', 20.1)": lambda b, cap: 0}
    dress = {"path, tie passed (20.1)": lambda b, cap: n[b] - passed[b],
             "path, tie unpassed (186.1)": lambda b, cap: n[b] - (0 if tie[b] else passed[b]),
             "all in-legs (198.3, 185.6; hub face)": lambda b, cap: n[b],
             "all but one (155.1 dock face)": lambda b, cap: n[b] - 1,
             "stop clause (21.1)": lambda b, cap: n[b] if cap[b] else n[b] - (0 if tie[b] else passed[b]),
             "none": lambda b, cap: 0}
    return rec, strike, dress


def readings(W, P, born):
    print("[READINGS] uniform readings: one records typing, one strike rule, one dressing rule for every body; S = (records - strike) 3^u")
    rec, strike, dress = reading_axes(W, P, born)
    grid = {}
    for (rn, (rv, cap)), (sn, sf), (dn, df) in itertools.product(rec.items(), strike.items(), dress.items()):
        grid[(rn, sn, dn)] = {b: (rv[b] - sf(b, cap)) * 3 ** df(b, cap) for b in BODIES}
    printed = {
        "20.1 adopted reading": (("register (20.1)", "every entry (20.1)", "path, tie passed (20.1)"), (67, 69, 432, 2709)),
        "21.1 stop clause": (("register (20.1)", "uncapped entry alone (186.1)", "stop clause (21.1)"), (67, 621, 432, 8613)),
        "186.1 segment rule": (("reach (32.1; 186.1)", "uncapped entry alone (186.1)", "path, tie unpassed (186.1)"), (67, 963, 3969, 2871)),
        "185.12 / 152.1 hub (leptons)": (("reach (32.1; 186.1)", "every entry (20.1)", "all in-legs (198.3, 185.6; hub face)"), (201, 963, 3969, None)),
        "179.1 typing uncapped (proton)": (("measure, all uncapped, seed (179.1)", "every entry (20.1)", "path, tie passed (20.1)"), (None, None, None, 747)),
        "195.1 mass-register row": (("measure, all uncapped, sink (195.1)", "every entry (20.1)", "all in-legs (198.3, 185.6; hub face)"), (201, 963, 3969, 1701)),
    }
    print(f"  grid: {len(rec)} records typings x {len(strike)} strike rules x {len(dress)} dressing rules = {len(grid)} readings "
          f"(20.1 names a grid of 'four of the register, three of the records and five of the path' = 60)")
    print("  validation against the printed tuples (e, mu, tau, p):")
    for lab, (key, want) in printed.items():
        got = tuple(grid[key][b] for b in BODIES)
        okk = all(w is None or w == g for w, g in zip(want, got))
        print(f"    {lab:32s} {got}  printed {tuple('-' if w is None else w for w in want)}  {'equal' if okk else 'DIFFERS'}")
        assert okk, lab
    # which readings meet the derived cells
    meet = {b: [k for k, v in grid.items() if v[b] == DERIVED[b]] for b in BODIES}
    allfour = [k for k, v in grid.items() if all(v[b] == DERIVED[b] for b in BODIES)]
    no186 = [k for k in grid if k[1] != "uncapped entry alone (186.1)"]
    print(f"  readings meeting each derived cell: e 67 at {len(meet['e'])}, mu 963 at {len(meet['mu'])}, tau 3969 at {len(meet['tau'])}, "
          f"p 2871 at {len(meet['p'])} of {len(grid)}")
    p_nostr = [k for k in meet["p"] if k[1] == "none ('nothing removed', 20.1)"]
    print(f"  p 2871 with the strike rules of 20.1 only (every entry, none): {len([k for k in meet['p'] if k in no186])} readings, all at "
          f"'nothing removed' and a capped register ({len(p_nostr)}; 20.1: 'The proton is met at six readings, every one at the capped register of its entered link with nothing removed')")
    leptons3 = [k for k in no186 if all(grid[k][b] == DERIVED[b] for b in ("e", "mu", "tau"))]
    print(f"  the three lepton cells at once, 20.1 strike rules: {len(leptons3)} reading(s): {[k for k in leptons3]}")
    print(f"  all four derived cells at once: {len(allfour)} of {len(grid)} readings: {allfour}; without the 186.1 strike clause: "
          f"{len([k for k in allfour if k in no186])} (20.1: 'no one reading meets the four derived cells at once')")
    # which readings put the ratios within 1% (exponent 2 and exponent 1)
    print("  ratios within 1% of the measured value, by exponent (mass = S^p):")
    for p in (2, 1):
        ok = {b: {k for k, v in grid.items() if v["e"] > 0 and v[b] > 0 and abs(rel((R(v[b], v["e"])) ** p, TGT[b][0])) <= 0.01} for b in ("mu", "tau", "p")}
        both = ok["mu"] & ok["tau"]; three = both & ok["p"]
        print(f"    p = {p}: mu/e {len(ok['mu'])}, tau/e {len(ok['tau'])}, p/e {len(ok['p'])}, mu and tau {len(both)}, all three {len(three)} of {len(grid)}")
        if p == 2:
            for lab, s in (("mu and tau", both), ("p/e", ok["p"]), ("all three", three)):
                for k in sorted(s): print(f"      {lab:10s} {k}  -> {tuple(grid[k][b] for b in BODIES)}")
    print()
    return grid


# ------------------------------------------------------------------ FACES (155.1)
def faces(src):
    print("[FACES] the 155.1 faces: dock = every body with one in-leg absorbed into the unit (M 3^k / 67), hub = no leg absorbed")
    print("        (M 3^(k+1) / 201), mixed = 'the hub numerator over the dock unit'; ranks M and in-legs as in WEBS, proton records 319")
    hub = {"e": 67 * 3, "mu": 107 * 9, "tau": 147 * 27, "p": 319 * 27}
    dock = {b: v // 3 for b, v in hub.items()}
    path = DERIVED
    rows = [("dock / dock", dock, dock["e"]), ("hub / hub", hub, hub["e"]), ("mixed: hub / dock unit", hub, dock["e"]),
            ("186.1 path reading", path, path["e"])]
    print(f"  numerators: dock {dock}; hub {hub}; 186.1 {path}")
    print(f"  {'face':24s} {'m_mu/m_e':>12s} {'m_tau/m_e':>12s} {'m_p/m_e':>12s}   (mass = S^2; measured {MMU:.4f}, {MTAU:.3f}, {MP:.4f})")
    for lab, num, den in rows:
        r = [R(num[b], den) ** 2 for b in ("mu", "tau", "p")]
        flags = ["ok" if abs(rel(x, TGT[b][0])) <= 0.01 else "x%.2f" % (float(x) / TGT[b][0]) for x, b in zip(r, ("mu", "tau", "p"))]
        print(f"  {lab:24s} {float(r[0]):12.4f} {float(r[1]):12.3f} {float(r[2]):12.4f}   {flags}")
    print("  in-legs absorbed per body: dock face 1, 1, 1, 1; hub face 0, 0, 0, 0; 186.1 path reading e 1 of 1, mu 0 of 2, tau 0 of 3,")
    print("  p 1 of 3. The proton ratio is the same at the dock and hub faces (2871/67 = 8613/201) and at 186.1; the lepton ratios")
    print("  meet only where the electron absorbs a leg and the muon and tau do not (the mixed face). No single face gives all three.")
    tex = "\n".join(p_.read_text(encoding="utf-8") for p_ in sorted((src / "tex").rglob("*.tex")))
    cnt = {v: len(re.findall(r"(?<![\d.,])" + str(v) + r"(?![\d])", tex)) for v in (67, 201, 321, 963, 1323, 3969, 2871, 8613)}
    print(f"  occurrences of the numerators as whole numbers in the tex sources: {cnt} (the muon and tau dock-face numerators are never printed)\n")


# ------------------------------------------------------------------ MASSLAW
def masslaw():
    print("[MASSLAW] the agreeing cells under each printed mass law (exact; float to 12 digits)")
    ce = field(D.CLOSURES["e"][0]); cmu = field(D.CLOSURES["mu"][0]); ctau = field(D.CLOSURES["tau"][0])
    cpm = field(D.CLOSURES["p_massrow"][0]); g2 = field(D.G2[0])
    S = {"e": 67, "mu": 963, "tau": 3969, "p": 2871}
    laws = [
        ("m = S^2 (185.9 leading term; 176.2 first order)", {b: S[b] ** 2 for b in S}),
        ("m = office = S^2 + closure (195.1 T; = 185.9 with its body term)", {"e": 67 ** 2 + ce, "mu": 963 ** 2 + cmu, "tau": 3969 ** 2 + ctau, "p": 2871 ** 2}),
        ("m = office, proton closure at the mass-row weight (O.1)", {"e": 67 ** 2 + ce, "mu": 963 ** 2 + cmu, "tau": 3969 ** 2 + ctau, "p": 2871 ** 2 + cpm}),
        ("m = T1 (1 + G2^2) (176.2), one G2 at e, mu and the (2,3) proton entry", {"e": 67 ** 2 * (1 + g2 ** 2), "mu": 963 ** 2 * (1 + g2 ** 2), "tau": 3969 ** 2, "p": 2871 ** 2 * (1 + g2 ** 2)}),
        ("m = T1 (1 + G2^2) (176.2), G2 at e, mu only (capped proton entry, tau pendant)", {"e": 67 ** 2 * (1 + g2 ** 2), "mu": 963 ** 2 * (1 + g2 ** 2), "tau": 3969 ** 2, "p": 2871 ** 2}),
        ("m = S (office = m^2: 181.5, 200.6 bridges; 200.3 exponent one)", {b: S[b] for b in S}),
    ]
    vals = {}
    print(f"  {'law':78s} {'m_p/m_e':>16s} {'rel':>10s} {'m_mu/m_e':>11s} {'rel':>9s} {'m_tau/m_e':>10s} {'rel':>8s} {'Koide':>9s}")
    for lab, m in laws:
        rp, rm, rt = (sp.N(m[b] / m["e"], 15) for b in ("p", "mu", "tau"))
        sq = [sp.sqrt(m[b]) for b in ("e", "mu", "tau")]
        kq = sp.N(sum(m[b] for b in ("e", "mu", "tau")) / sum(sq) ** 2, 10)
        vals[lab] = float(rp)
        print(f"  {lab:78s} {float(rp):16.6f} {rel(rp, MP):+10.2e} {float(rm):11.4f} {rel(rm, MMU):+9.2e} {float(rt):10.2f} {rel(rt, MTAU):+8.2e} {float(kq):9.6f}")
    for lab, m in laws[:2]:
        ns = [(float(sp.N(m[b] / m["e"], 20)) - TGT[b][0]) / TGT[b][1] for b in ("p", "mu", "tau")]
        print(f"  as exact statements, {lab.split(' (')[0]}: n_sigma p {ns[0]:+.2e}, mu {ns[1]:+.2e}, tau {ns[2]:+.0f}")
    lo, hi = min(v for k, v in vals.items() if "m = S (" not in k), max(v for k, v in vals.items() if "m = S (" not in k)
    print(f"  electron closure / 67^2 = {float(ce / 67**2):.3e}; G2^2 = {float(g2**2):.3e}; proton mass-row closure / 2871^2 = {float(cpm / 2871**2):.1e}")
    print(f"  spread of m_p/m_e over the printed second-order laws (S^2 family): {lo:.6f} to {hi:.6f}, relative width {hi/lo-1:.2e};")
    print(f"  CODATA 2022 {MP} lies {'inside' if lo <= MP <= hi else 'outside'} that band; deviations {lo/MP-1:+.2e} to {hi/MP-1:+.2e};")
    print(f"  in sigma (3.2e-8): {(lo-MP)/sMP:+.2e} to {(hi-MP)/sMP:+.2e}. Measured Koide (PDG 2024) {KOIDE:.7f}.\n")
    return hi / lo - 1


# ------------------------------------------------------------------ SPACES
def lnd(v, t): return abs(math.log(float(v) / t))


def density(vals, t, w):
    vs = {float(sp.N(v, 14)) if isinstance(v, sp.Basic) else round(float(v), 10) for v in vals if float(v) > 0}
    return sum(1 for v in vs if abs(math.log(v / t)) < w) / (2 * w), len(vs)


def lee(name, vals, t, delta, extra=""):
    vals = [v for v in vals if float(v) > 0]
    ds = sorted({round(lnd(v, t), 14) for v in vals})
    best = ds[0]
    out = [f"  {name}: {len(vals)} variants, {len({round(float(v), 10) for v in vals})} distinct; closest |ln dev| {best:.2e}"
           f" (rank 1); next {ds[1]:.2e}" if len(ds) > 1 else ""]
    for w in (0.5, 1.0):
        rho, _ = density(vals, t, w)
        lam = 2 * delta * rho
        out.append(f"    window |ln| < {w}: {rho:.1f} distinct values per unit ln -> expected chance hits within {delta:.2e}: {lam:.2e}, "
                   f"P(>=1) = {1 - math.exp(-lam):.2e}")
    print("\n".join(out) + extra)


def spaces(grid, W, spread):
    print("[SPACES] distributions of predicted ratios and look-elsewhere estimates. A 'variant' is one (numerator cell, electron cell,")
    print("         exponent) choice; lambda = 2 delta x (distinct values per unit ln near the target); P = 1 - exp(-lambda)")
    dp = lnd(R(2871, 67) ** 2, MP); dmu = lnd(R(963, 67) ** 2, MMU); dtau = lnd(R(3969, 67) ** 2, MTAU)
    print(f"  achieved |ln dev|: p {dp:.3e}, mu {dmu:.3e}, tau {dtau:.3e}; source's own convention spread for p (MASSLAW) {spread:.2e}")
    # (1) printed cells
    pc = [(sp.Rational(460867, 3072) + sp.Rational(25, 1536) * sp.sqrt(3)) if v == "214.6" else v for v, _, _ in D.PROTON_CELLS]
    eu = [v for v, _, _ in D.ELECTRON_UNITS]
    rows = []
    for (c, lab, _), cv in zip(D.PROTON_CELLS, pc):
        for e in eu:
            for p in (1, 2):
                x = (cv / sp.Integer(e)) ** p
                rows.append((float(sp.N(x, 15)), c, e, p, lab))
    rows.sort(key=lambda r: lnd(r[0], MP))
    print(f"  (1) printed space for m_p/m_e: {len(D.PROTON_CELLS)} printed proton cells x {len(eu)} printed electron cells {eu} x p in {{1,2}}")
    print("      ten closest (value, cell, e, p, rel. dev.):")
    for r in rows[:10]: print(f"        {r[0]:14.6f}  {str(r[1]):>6s} / {r[2]:3d}  p={r[3]}  {r[0]/MP-1:+.3e}   {r[4]}")
    q = sorted(abs(math.log(r[0] / MP)) for r in rows)
    print(f"      |ln dev| quantiles: median {q[len(q)//2]:.2f}, 10% {q[len(q)//10]:.3f}; within 1%: {sum(1 for x in q if x < 0.01)}, "
          f"within 10%: {sum(1 for x in q if x < 0.0953)}, within a factor 2: {sum(1 for x in q if x < 0.693)}")
    lee("printed space", [r[0] for r in rows], MP, dp)
    lee("printed space at the convention spread", [r[0] for r in rows], MP, spread)
    for b, cells, tg, dd in (("mu", D.MUON_CELLS, MMU, dmu), ("tau", D.TAU_CELLS, MTAU, dtau)):
        vv = [(R(c, e)) ** p for c, _, _ in cells for e in eu for p in (1, 2)]
        lee(f"printed space for m_{b}/m_e ({len(cells)} cells x {len(eu)} x 2)", vv, tg, dd)
    jt = [(e, m, t) for e in eu for m, _, _ in D.MUON_CELLS for t, _, _ in D.TAU_CELLS
          if abs(rel(R(m, e) ** 2, MMU)) <= 0.01 and abs(rel(R(t, e) ** 2, MTAU)) <= 0.01]
    print(f"    joint, printed: mu/e and tau/e both within 1% (p = 2) at {len(jt)} of {len(eu) * len(D.MUON_CELLS) * len(D.TAU_CELLS)} (e, mu, tau): {jt}")
    # (2) per-body choice over the reading grid
    cells = {b: sorted({v[b] for v in grid.values() if v[b] > 0}) for b in BODIES}
    print(f"  (2) reading grid, each body's cell chosen independently: cells e {cells['e']}; mu {len(cells['mu'])}; tau {len(cells['tau'])}; p {cells['p']}")
    for b, dd in (("p", dp), ("mu", dmu), ("tau", dtau)):
        vv = [(R(c, e)) ** p for c in cells[b] for e in cells["e"] for p in (1, 2)]
        lee(f"per-body choice, m_{b}/m_e", vv, TGT[b][0], dd)
    vv = [(R(v["p"], v["e"])) ** p for v in grid.values() for p in (1, 2) if v["e"] > 0 and v["p"] > 0]
    lee("uniform readings only, m_p/m_e", vv, MP, dp)
    # (3) the round-2 declared space (scorecard_freedom.py part P), recomputed here on the same footing
    Rr = [96, 319, 906, 1036, 2053, 5100, 81, 85, 101, 125, 150]
    dress = sorted({3 ** a * 4 ** b for a in range(5) for b in range(5) if a + b <= 4})
    c3 = sorted({(r - s) * d for r in Rr for s in (0, 18) for d in dress})
    vv = [(R(c, e)) ** p for c in c3 for e in (67, 201) for p in (1, 2)]
    lee(f"(3) round-2 declared space ({len(c3)} cells x {{67, 201}} x 2)", vv, MP, dp)
    # (4) the narrow family with 319 and 67 fixed upstream
    fam = [(R(319 * 3 ** u, 67 * 3 ** v)) ** p for u in range(5) for v in range(2) for p in (1, 2)]
    hits = sorted({(u, v, p) for u in range(5) for v in range(2) for p in (1, 2) if lnd(R(319 * 3 ** u, 67 * 3 ** v) ** p, MP) <= dp * 1.0000001})
    print(f"  (4) narrow family (319 3^u / (67 3^v))^p, u 0..4, v 0..1, p 1..2: hits at the achieved precision {hits}; spacing in ln is p ln 3,")
    lamf = sum(2 * dp / (p * math.log(3)) for p in (1, 2))
    print(f"      so a fixed (u - v, p) rule is hit by chance with probability {2*dp/(2*math.log(3)):.1e} (p = 2) and some p with {lamf:.1e}")
    # (5) the source's own lepton censuses: only the number of distinct first powers and their range are printed
    print("  (5) the saturation censuses as candidate spaces (154.1: 300 three-hub sinks, 63 distinct first powers 369..4617; 189.1:")
    print("      231 two-hub sinks, 92 values 95..3321); values not printed, so log-uniform spacing is assumed; e = 67, p = 2:")
    for lab, nd, lo_, hi_, tg, dd in (("tau", 63, 369, 4617, MTAU, dtau), ("mu", 92, 95, 3321, MMU, dmu)):
        rho = nd / (2 * math.log(hi_ / lo_))
        print(f"      m_{lab}/m_e: target cell 67 sqrt(ratio) = {67 * math.sqrt(tg):.1f} inside [{lo_}, {hi_}]; {rho:.1f} values per unit ln -> "
              f"expected chance hits within {dd:.2e}: {2 * dd * rho:.3f}")
    print("      the proton saturation (175.1, 175.3: 16946 and 44341 sinks) prints no first powers, so no density can be formed;")
    print("      the source reports only that none equals the candidate integer")
    # (6) any integer
    st = 67 * math.sqrt(MP)
    print(f"  (6) any-integer baselines over e = 67, p = 2: target numerator {st:.4f}; nearest integer 2871 (off {2871 - st:+.4f});")
    for lab, dl in (("achieved", dp), ("convention spread", spread)):
        print(f"      {lab:18s} delta {dl:.2e}: any integer {min(1, dl * st):.3f}; multiples of 9 {min(1, dl * st / 9):.4f}; "
              f"with the electron cell free in {eu}: any integer {1 - math.prod(1 - min(1, dl * e * math.sqrt(MP)) for e in eu):.3f}")
    print()


# ------------------------------------------------------------------ LEPTON
def lepton():
    print("[LEPTON] face-free content and the factor nine")
    tm = R(3969, 963) ** 2
    print(f"  m_tau/m_mu = (3969/963)^2 = (441/107)^2 = {float(tm):.4f} at every face (both numerators hub); PDG 2024 {TAUMU:.4f}({sTAUMU*1e4:.0f}): "
          f"{rel(tm, TAUMU):+.3%}, {(float(tm)-TAUMU)/sTAUMU:+.0f} sigma")
    h_mu, h_tau = float(R(963, 201) ** 2), float(R(3969, 201) ** 2)
    print(f"  face-consistent ratios (176.2 prints them): m_mu/m_e = {h_mu:.4f}, m_tau/m_e = {h_tau:.3f}; measured/predicted = "
          f"{MMU/h_mu:.4f}, {MTAU/h_tau:.4f}; nearest power of 3: 9 ({math.log(MMU/h_mu/9):+.2e}, {math.log(MTAU/h_tau/9):+.2e} in ln)")
    for lab, base in (("3^n", 3), ("9^n", 9)):
        print(f"  chance that a fixed prediction lands within {abs(math.log(MMU/h_mu/9)):.2e} of a factor {lab} away from the measured m_mu/m_e: "
              f"{2*abs(math.log(MMU/h_mu/9))/math.log(base):.1e}")
    print(f"  given mu/e, tau/e follows from tau/mu; its deviation {rel(R(3969,67)**2, MTAU):+.3%} is tau/mu's {rel(tm, TAUMU):+.3%} plus mu/e's "
          f"{rel(R(963,67)**2, MMU):+.3%}")
    for lab, cell, tg, step in (("mu", 963, MMU, 9), ("tau", 3969, MTAU, 27)):
        t_ = 67 * math.sqrt(tg); near = round(t_ / step) * step
        print(f"  integer level, m_{lab}/m_e at e = 67: target cell {t_:.2f}, printed {cell} (off {cell - t_:+.2f}); nearest integer "
              f"{round(t_)}, nearest multiple of {step} {near}; P(any integer this close) = {min(1, 2 * abs(cell - t_)):.2f}, "
              f"P(multiple of {step}) = {min(1, 2 * abs(cell - t_) / step):.3f}")
    for lab, e in (("mixed (e 67)", 67), ("hub (e 201)", 201)):
        k = R(e ** 2 + 963 ** 2 + 3969 ** 2, (e + 963 + 3969) ** 2)
        print(f"  Koide at the {lab} face: {k} = {float(k):.6f} ({rel(k, KOIDE):+.3%} against {KOIDE:.7f}, {(float(k)-KOIDE)/sKOIDE:+.0f} sigma)")
    print()


# ------------------------------------------------------------------ CROSS
def cross():
    print("[CROSS] the same segment walk at the quark sinks (181.5, THEOREM: 'entered at the dock and stopped at the hub', capped entry)")
    me = PDG[11][0]
    q = [("up arm", 319, 2), ("down arm", 319, 1), ("strange sink", 25839, 3), ("charm sink", 8613, 4), ("bottom sink", 8613, 5), ("top sink", 2871, 6)]
    print(f"  {'body':14s} {'S':>6s} {'(S/67)^2':>12s} {'S/67':>8s} {'measured m/m_e':>15s}  (PDG 2024: u, d, s MS-bar 2 GeV; c, b MS-bar at m; t direct)")
    for lab, s_, pid in q:
        print(f"  {lab:14s} {s_:6d} {float(R(s_,67)**2):12.2f} {s_/67:8.2f} {PDG[pid][0]/me:15.1f}")
    print(f"  the top and the proton share the cell 319*9: at either exponent m_t = m_p, against m_t/m_p = {PDG[6][0]/(MP*me):.1f}; ordering s > c = b > t\n")


# ------------------------------------------------------------------ DIGITS
def digits():
    print("[DIGITS] information put in against agreement obtained")
    for lab, num, den, t in (("m_p/m_e", 2871, 67, MP), ("m_mu/m_e", 963, 67, MMU), ("m_tau/m_e", 3969, 67, MTAU)):
        d = -math.log10(abs(rel(R(num, den) ** 2, t)))
        print(f"  {lab:9s} ({num}/{den})^2: integers carry log10({num}) + log10({den}) = {math.log10(num) + math.log10(den):.2f} digits; "
              f"agreement {d:.2f} digits; numerator alone {math.log10(num):.2f}")
    print("  For a nearest-integer numerator the expected agreement is log10(S) + log10(2) digits: 3.76 for S ~ 2871; obtained 4.74.\n")


# ------------------------------------------------------------------ SWEEP
KNOWN = {**{v: "P" for v, _, _ in D.PROTON_CELLS if isinstance(v, int)},
         **{v: "L" for L in (D.ELECTRON_UNITS, D.MUON_CELLS, D.TAU_CELLS) for v, _, _ in L},
         # registers, measures, ranks, counts, squares
         **{v: "N" for v in (64, 79, 80, 81, 83, 85, 101, 107, 113, 117, 125, 133, 144, 155, 219, 234, 243, 276, 300, 319, 1296, 2304,
                             4489, 7225, 6561, 20736, 31104, 101761, 8242641)},
         # other bodies: neutron (8127, 10773, 8046, 9315, 1125, 1596, 3591), deuteron (1647, 1971, 2873, 4941, 5913), carrier-kind
         # sinks of 69.1 (2565, 2889, 3213), tau-census range of 189.1 (95, 3321) and 154.1 (369, 4617), quark bodies (120, 1044),
         # 116.1 insertion readings (71, 75, 999, 1035, 1071, 4077, 4293) and metric increment 542, 229.1 joined network (1746),
        # 189.1 lepton-census ranges and counts (123, 231, 630), 179.1 sink id 94, 177.1 strange-docked body 720; census sizes
         **{v: "O" for v in (8127, 10773, 8046, 9315, 1125, 1596, 3591, 1647, 1971, 2873, 4941, 5913, 2565, 2889, 3213, 95, 3321, 369,
                             4617, 120, 1044, 71, 75, 999, 1035, 1071, 4077, 4293, 1746, 957, 3857, 1236, 516, 1972, 25558, 16946, 44341,
                             94, 123, 231, 542, 630, 720)}}


def sweep(U):
    print("[SWEEP] integers 50..99999 within 120 characters after a first-power keyword, in statements whose box names the proton or")
    print("        nucleon (fractions, \\texttt office strings, picture coordinates and surd coefficients removed first). Classes: P printed")
    print("        proton cell (in the space); L lepton cell; N register, measure, rank, count or square; O another body or a census size;")
    print("        ? unclassified")
    key = re.compile(r"first powers?|S_\{?1\}?\s*=|S_\{?p\}?\s*=|candidate cell|numerator|spelling", re.I)
    num = re.compile(r"(?<![\d.])(\d{1,3}(?:\\,\d{3})+|\d+)(?![\d.])")
    strip = [r"\\t?frac\{[^{}]*\}\{[^{}]*\}", r"\\texttt\{(?:[^{}]|\{\})*\}", r"\d+\s+of\s+\d+", r"\\(?:qbezier|put)\([^)]*\)(?:\([^)]*\))*", r"\d+\s*/\s*\d+",
             r"\d[\d\\, ]*\s*over\s*\d[\d\\, ]*", r"\d+\s*\\sqrt", r"\\label\{[^}]*\}"]
    seen = collections.defaultdict(set)
    for n_, u in U.items():
        if not re.search(r"proton|nucleon", u["body"], re.I): continue
        t = re.sub(r"\s+", " ", u["body"] + " " + u["rest"])
        for s_ in strip: t = re.sub(s_, " ", t)
        for m in key.finditer(t):
            for k in num.finditer(t[m.start(): m.end() + 120]):
                v = int(k.group(1).replace("\\,", ""))
                if 50 <= v < 100000: seen[v].add(n_)
    byc = collections.defaultdict(list)
    for v in sorted(seen): byc[KNOWN.get(v, "?")].append(v)
    for c in "PLNO": print(f"  {c}: {len(byc[c])} {byc[c]}")
    print(f"  ?: {len(byc['?'])} " + ", ".join(f"{v} {sorted(seen[v])}" for v in byc["?"]))
    miss = [v for v, _, _ in D.PROTON_CELLS if isinstance(v, int) and v not in seen]
    print(f"  printed proton cells not caught by the sweep (anchored by hand): {miss}")
    print("  excluded by reason (masscells_data.EXCLUDED):")
    for s, r in D.EXCLUDED: print(f"    {s}: {r}")
    print()


def main():
    src = pathlib.Path(sys.argv[1])
    if "--pdg-mcd" in sys.argv:
        n = check_mcd(sys.argv[sys.argv.index("--pdg-mcd") + 1])
        print(f"[PDG] transcription checked against mass_width_2024.mcd ({n} particle ids)")
    print(f"[TARGETS] CODATA 2022: m_p/m_e = {MP}({sMP*1e9:.0f}), m_mu/m_e = {MMU}({sMMU*1e7:.0f}); PDG 2024: m_tau/m_e = {MTAU:.3f}({sMTAU*1e3:.0f}) "
          f"(CODATA 2022 {CODATA['tau-electron mass ratio'][0]}(23)), m_tau/m_mu = {TAUMU:.5f}({sTAUMU*1e5:.0f}), Koide {KOIDE:.7f}({sKOIDE*1e7:.0f})\n")
    U = anchors(src)
    born = records()
    W = build_webs(U)
    P = paths(W)
    grid = readings(W, P, born)
    faces(src)
    spread = masslaw()
    spaces(grid, W, spread)
    lepton()
    cross()
    digits()
    sweep(U)
    if "--tsv" in sys.argv:
        with open(sys.argv[sys.argv.index("--tsv") + 1], "w", encoding="utf-8") as f:
            f.write("records\tstrike\tdressing\tS_e\tS_mu\tS_tau\tS_p\tmu_e_p2\ttau_e_p2\tp_e_p2\tp_e_p1\n")
            for (rn, sn, dn), v in grid.items():
                r2 = [f"{float(R(v[b], v['e'])**2):.8g}" if v["e"] > 0 else "" for b in ("mu", "tau", "p")]
                f.write(f"{rn}\t{sn}\t{dn}\t{v['e']}\t{v['mu']}\t{v['tau']}\t{v['p']}\t" + "\t".join(r2) +
                        f"\t{float(R(v['p'], v['e'])):.8g}\n" if v["e"] > 0 else "\n")


if __name__ == "__main__":
    main()
