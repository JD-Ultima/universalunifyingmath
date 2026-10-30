#!/usr/bin/env python3
"""Round 5: the hydrogen body registers of 143.1, the ladder beyond n = 8, and the comparison with experiment.

Usage: python3 -I hydrogen_registers.py <compendium_dir>

Blocks (each prints its own header):
  ANCHORS   every quoted string is checked verbatim (whitespace-normalised) against its statement unit or file
  FACES     the per-link faces from the end dimensions (16.1, 103.1): whole (fc)^2, cut (f^2-1)(c^2-1), rank fc,
            support min(f,c)^4; checked against the 16 printed genre faces, 122.1, 164.4 and the count registers
  BODIES    the product rule "the register of a body at a face the product over its links of the per-link face"
            (103.1; 147.1) at the three bodies whose registers the source prints (W word 147.1, Z body and gluon
            word 149.1); cut at the port only against the cut at every link
  ELECTRON  the electron web parsed from the 152.1 figure; the whole and cut body registers under each reading;
            the dock and hub units 67 and 201
  CARRIER   the photon word parsed from the 143.1 figure; the printed link counts and the expressions equal to 36
  E1        the level-one binding 36 u / R at every printed unit and every register reading; the printed equalities
  LADDER    (b) every printed ladder sentence (range of n); reconstruction R exactly at n = 1..200 (commutant of the
            sl2 action on Sym^(n-1)(C^2)), the support dimension m in {2,3,4,6}, order and per-record readings
  COMPARE   (c) E1 against R_inf hc / m_e c^2 and the measured ionization energy, a declared conversion family, and
            the unit-free ratio E_2/E_1 against the measured 2S/1S binding ratio
Exact arithmetic with fractions and sympy; floats only in COMPARE. Read-only: nothing in the source is executed.
References (COMPARE): CODATA 2022 via scipy.constants 1.18.1 (alpha, R_inf hc, m_e c^2, h, e, c), asserted against the
transcription; H ionization energy 109678.77174307(10) cm^-1 (NIST ASD, Kramida, At. Data Nucl. Data Tables 96, 586
(2010)) = 13.598434599702 eV; 1S-2S hyperfine centroid 2466061413187035(10) Hz (Parthey et al., PRL 107, 203001
(2011); arXiv:1107.3101, where the value and the 'hyperfine centroid' label are stated).
"""
import itertools, math, pathlib, re, sys
from fractions import Fraction as Fr
import sympy as sp

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))   # our own directory only
import scorecard_extract as X
import masscells as MC          # parse_figure (round 3); importing it also asserts scipy 1.18 and the CODATA transcription

import scipy.constants as sc

FX = lambda q: "*".join(f"{p}^{e}" if e > 1 else f"{p}" for p, e in sorted(sp.factorint(q).items())) if q > 1 else "1"


def frac_str(q):
    q = Fr(q)
    if q.denominator == 1: return FX(q.numerator)
    d = FX(q.denominator)
    return f"{FX(q.numerator)}/" + (f"({d})" if "*" in d or "^" in d else d)


# ------------------------------------------------------------------ ANCHORS
ANCHORS = [
    # 143.1: the binding, the registers, the read rule, the mass face, the ladder
    ("143.1", "by the read through the electron body at its dock port $B$-$w_{1}$ (the mass face multiplied by the unit over the register of the body at that face)"),
    ("143.1", "the whole face with the hub unit two hundred and one over the body register two to the thirtieth multiplied by three to the eighteenth"),
    ("143.1", "the cut face with the dock unit sixty-seven over the body register two to the thirty-first multiplied by three to the seventeenth (the port at its cut twenty-four)"),
    ("143.1", "both faces printed, and the port-register readings beside, the mass face multiplied by sixty-seven over the port registers thirty-six and twenty-four"),
    ("143.1", "E_{1}^{\\mathrm{whole}} = 36 \\cdot \\frac{201}{2^{30}\\cdot 3^{18}} = \\frac{67}{2^{28}\\cdot 3^{15}}"),
    ("143.1", "E_{1}^{\\mathrm{cut}} = 36 \\cdot \\frac{67}{2^{31}\\cdot 3^{17}} = \\frac{67}{2^{29}\\cdot 3^{15}}"),
    ("143.1", "The two body registers of the level-one binding are cells of the read tables of this result, printed in their prime factorizations and re-derived on the page not"),
    ("143.1", "by the docking act (the distilling pair, one link made at each end and one contraction at one link)"),
    ("143.1", "the mass face of the held chain thirty-six and the exchanged mass face six root six"),
    ("143.1", "\\text{the held chain: } 6 \\ \\text{sites}, \\ 5 \\ \\text{links of type } (2,3); \\qquad M_{L} = 36"),
    ("143.1", "The fine structure is the dark cut of the bond at every one of the eight links of the held chain"),
    ("143.1", "By exhaustive enumeration of the eight ranks of the held chain and of the eight links of the bond"),
    ("143.1", "Docking at one end alone moves the hold cell from a theorem to none with the mass face undefined"),
    ("143.1", "the level-$n$ binding as the level-one binding over $n$ squared at every rank from one to eight, the register of the ladder equal to the dimension of the register of the pair"),
    ("143.1", "E_{n} = \\frac{E_{1}}{n^{2}}, \\qquad n = 1, \\dots, 8"),
    ("143.1", "the per-record ladder is exact at every rank and no pure power of the rank, and the two-vertex ladder equals one over the rank squared at all eight ranks with no step posed"),
    ("143.1", "exact at every rank and a pure power of the rank at none"),
    ("143.1", "its level-one binding through the electron body is an exact element of the field in the electron unit"),
    ("143.1", "the level-one binding is the two-vertex term at the bond, the Rydberg constant in the electron unit the ladder cell of the held carrier"),
    ("143.1", "The levels of hydrogen fall as the inverse square of the principal quantum number"),
    ("143.1", "The correspondence is by name"),
    # the per-link faces and the product rule
    ("16.1", "From the whole register $f^{2} c^{2}$ and the cut register $(f^{2} - 1)(c^{2} - 1)$, by the register of the count at order one, the cut multiplied by the whole multiplied by the whole"),
    ("16.1", "two thousand three hundred four at (2,2) and fourteen million seven hundred forty-five thousand six hundred at (4,4)"),
    ("103.1", "the register of a body at a face the product over its links of the per-link face"),
    ("103.1", "its four faces are the whole, the square of the product of the ends, the cut, the product of the squared ends less one each, the rank, the product of the ends, and the support, the whole again at equal ends"),
    ("103.1", "the faces of the four genres from the dimensions read 81, 64, 9, 81 at (3,3), 36, 24, 6, 16 at (2,3), 256, 225, 16, 256 at (4,4) and 64, 45, 8, 16 at (2,4)"),
    ("122.1", "the whole register $36=4\\cdot 9$ and the cut register $24=3\\cdot 8$"),
    ("164.4", "At the qutrit-qubit genre the whole 36, the cut 24, the common-support register 16 and the two links."),
    ("214.6", "The register of the common support is the smaller end of the port genre to the fourth power."),
    ("147.1", "The product register of the body at the whole port face is the product over its four links of the squared measures, the registers sixty-four, thirty-six, thirty-six and sixty-four"),
    ("147.1", "R = (8 \\cdot 6 \\cdot 6 \\cdot 8)^{2} = 5308416"),
    ("147.1", "over the star of the one link, thirty-one thousand one hundred and four"),
    ("147.1", "its measure twenty-eight, its first power twenty-eight and its cut chain sum twenty-two, its port the second link of the qutrit-qubit type at the common support two with the register sixteen at the support and thirty-six at the whole port"),
    ("147.1", "multiplied by the first power of the body over its product register"),
    ("149.1", "the links $s_{1}s_{2}$ at (2,4), $s_{2}s_{3}$ at (2,3), $s_{3}s_{4}$ at (2,3) and $s_{4}s_{5}$ at (2,4), the measure 28 its rest energy, and the port $s_{2}s_{3}$ at (2,3)"),
    ("149.1", "At the Z body the registers are 5308416 at the whole face and 3538944 at the cut face, the chain sums 28 and 22, and the factors $7/1327104$ and $11/1769472$"),
    ("149.1", "the links $s_{1}s_{2}$, $s_{2}s_{3}$ and $s_{3}s_{4}$ every one of genre (3,3), and the measure 27; its read link is the first (3,3) link in the order of the library, $s_{1}s_{2}$"),
    ("149.1", "Striking the scale body to the gluon word itself moves the factors to $1/19683$ and $1/23328$"),
    # units and the mass unit
    ("152.1", "so the unit of the mass readout stands at 201 at the hub reading and 67 at the dock reading"),
    ("152.1", "27 + 40 + 12 = 79, \\qquad \\operatorname{rank}(e) = 79 - 12 = 67"),
    ("192.1", "the tick at the one tick as one, the whole rung register as thirty-six, the common-support register of the electron dock as sixteen, the dock-struck unit as sixty-seven, and the first power at the probe hub leg as sixty-seven"),
    ("192.1", "so no one of the five spellings is the office and the mass unit stays a question of which cell is meant"),
    ("229.1", "the hub unit is reached by no joined read, so the hub-or-dock question of the previous releases resolves to the dock wherever the joined read opens at all"),
    ("229.1", "so the unit is the dock and the register is"),
    ("229.1", "16 = 2 \\cdot 2 \\cdot 2 \\cdot 2"),
    ("176.1", "The proton in the electron unit stands outside the comparison by orders"),
    ("9.4", "The open (2,3) link carries no mass face"),
    ("200.4", "The photon word open, its ends of dimensions 2 and 3."),
    # the ladder elsewhere
    ("45.1", "at the held pair the order-$k$ term at rank $n$ the term at rank one over $n$ squared at all eight ranks"),
    ("45.1", "the Rydberg ladder falls as the inverse square of the rank"),
]
FILE_ANCHORS = [
    ("front.tex", "The bond is the link that joins the electron body to the proton body in the hydrogen web, the link at which the two-vertex term is read."),
    ("back.tex", "The Rydberg constant as the two vertex term at the bond, a conditional theorem, treated at chapter 16-7 and 7-4, rests on Dark link cut admissibility (an open result);"),
    ("back.tex", "Rydberg ladder from the per record read, a conditional theorem, treated at chapter 16-7 and 7-4"),
]
SAME_LINE = ("back.tex", "The Rydberg constant as the two vertex term at the bond", "Mass unit as the rung register (an open question)")


def anchors(src):
    U = {u["num"]: u for u in X.units(src / "tex")}
    UN = {k: X.norm(u["body"] + " " + u["rest"]) for k, u in U.items()}
    print("[ANCHORS] verbatim check (whitespace-normalised) of the quoted strings against their statement units / files")
    miss = [(n, s) for n, s in ANCHORS if X.norm(s) not in UN.get(n, "")]
    fmiss = [(f, s) for f, s in FILE_ANCHORS if X.norm(s) not in X.norm((src / "tex" / f).read_text(encoding="utf-8"))]
    lines = (src / "tex" / SAME_LINE[0]).read_text(encoding="utf-8").splitlines()
    same = [i + 1 for i, l in enumerate(lines) if SAME_LINE[1] in l and SAME_LINE[2] in l]
    print(f"  statement units: {len(ANCHORS)} strings, {len(miss)} missing" + "".join(f"\n    MISSING {n}: {s[:90]}" for n, s in miss))
    print(f"  file level:      {len(FILE_ANCHORS)} strings, {len(fmiss)} missing" + "".join(f"\n    MISSING {f}: {s[:90]}" for f, s in fmiss))
    print(f"  back.tex: 'The Rydberg constant as the two vertex term at the bond' and 'Mass unit as the rung register (an open "
          f"question)' on one item line: lines {same}")
    cited = sorted({n for n, _ in ANCHORS}, key=lambda s: [int(x) for x in s.split(".")])
    G = dict(THEOREM="T", CONDITIONAL="C", OPEN="O", TABLE="Tab", DEFINITION="D")
    print("  grades of the units cited: " + ", ".join(f"{k} {G[U[k]['grade']]}" for k in cited))
    print(f"  locations: " + ", ".join(f"{k} {U[k]['file']}:{U[k]['line']}" for k in ("143.1", "16.1", "103.1", "147.1", "149.1", "152.1", "192.1", "229.1", "45.1")))
    assert not miss and not fmiss and same
    print()
    return U


# ------------------------------------------------------------------ FACES
def faces(d1, d2):
    f, c = sorted((d1, d2))
    return dict(whole=(f * c) ** 2, cut=(f * f - 1) * (c * c - 1), rank=f * c, support=f ** 4)


def faces_block():
    print("[FACES] per-link faces from the end dimensions (16.1, 103.1): whole (fc)^2, cut (f^2-1)(c^2-1), rank fc, support f^4 (f <= c)")
    printed = {(3, 3): (81, 64, 9, 81), (2, 3): (36, 24, 6, 16), (4, 4): (256, 225, 16, 256), (2, 4): (64, 45, 8, 16)}
    ok = 0
    for g, p in printed.items():
        F = faces(*g); got = (F["whole"], F["cut"], F["rank"], F["support"])
        ok += sum(a == b for a, b in zip(got, p))
        print(f"  {g}: whole {got[0]:4d} cut {got[1]:4d} rank {got[2]:3d} support {got[3]:4d}   printed (103.1) {p}  {'equal' if got == p else 'DIFFERS'}")
    print(f"  {ok} of 16 printed faces reproduced; (2,3) whole 36 = 4*9, cut 24 = 3*8 (122.1), support 16 (164.4, 214.6: smaller end^4)")
    cnt = {g: faces(*g)["cut"] * faces(*g)["whole"] ** 2 for g in ((2, 2), (2, 3), (3, 3), (4, 4))}
    print(f"  count register cut*whole^2 (16.1): (2,2) {cnt[(2,2)]}, (2,3) {cnt[(2,3)]} (147.1 star 31104), (3,3) {cnt[(3,3)]}, "
          f"(4,4) {cnt[(4,4)]} (16.1 prints 2304 and 14745600)")
    assert ok == 16 and cnt[(2, 2)] == 2304 and cnt[(4, 4)] == 14745600 and cnt[(2, 3)] == 31104
    print()


# ------------------------------------------------------------------ BODIES
def register(dims, links, face_at):
    """Product over the links of the per-link face; face_at(link) names the face read at that link."""
    return math.prod(faces(dims[a], dims[b])[face_at((a, b))] for a, b in links)


def measure(dims, links):
    return sum(dims[a] * dims[b] for a, b in links)


CUT = {"C1": "the port at its cut, every other link whole",
       "C2": "every link at its cut (the 103.1 rule read literally at the cut face)",
       "C3": "every (2,3) link at its cut (at the electron: both docks)"}


def cut_readings(dims, links, port, docks):
    """The register readings at the cut face that are tried (code -> face rule)."""
    return {"C1": lambda l: "cut" if l == port else "whole", "C2": lambda l: "cut",
            "C3": lambda l: "cut" if l in docks else "whole"}


def bodies_block():
    print("[BODIES] the product rule at the bodies whose registers the source prints; which cut-face reading reproduces them")
    Z = dict(s1=4, s2=2, s3=3, s4=2, s5=4); ZL = [("s1", "s2"), ("s2", "s3"), ("s3", "s4"), ("s4", "s5")]; Zp = ("s2", "s3")
    G = dict(s1=3, s2=3, s3=3, s4=3); GL = [("s1", "s2"), ("s2", "s3"), ("s3", "s4")]; Gp = ("s1", "s2")
    out = {}
    for name, d, L, p, S1, printed in (("W word / Z body (147.1, 149.1)", Z, ZL, Zp, 28, (5308416, 3538944, Fr(7, 1327104), Fr(11, 1769472))),
                                       ("gluon word (149.1 counter-test)", G, GL, Gp, 27, (None, None, Fr(1, 19683), Fr(1, 23328)))):
        assert measure(d, L) == S1
        Rw = register(d, L, lambda l: "whole")
        Scut = S1 - d[p[0]] * d[p[1]]         # cut chain sum: the port's measure removed (147.1: 28 -> 22)
        print(f"  {name}: measure {S1}, per-link whole {[faces(d[a], d[b])['whole'] for a, b in L]}, R_whole = {Rw} = {FX(Rw)}; "
              f"cut chain sum S1 - m(port) = {Scut}")
        if printed[0] is not None:
            print(f"    printed R_whole {printed[0]}: {'equal' if Rw == printed[0] else 'DIFFERS'}")
        print(f"    whole factor S1/R_whole = {Fr(S1, Rw)}  printed {printed[2]}: {'equal' if Fr(S1, Rw) == printed[2] else 'DIFFERS'}")
        for rname, rule in cut_readings(d, L, p, [l for l in L if {d[l[0]], d[l[1]]} == {2, 3}]).items():
            Rc = register(d, L, rule)
            fac = Fr(Scut, Rc)
            hit_R = printed[1] is not None and Rc == printed[1]
            hit_f = fac == printed[3]
            out[(name, rname)] = (hit_R or printed[1] is None) and hit_f
            print(f"    {rname} ({CUT[rname]}): R_cut = {Rc} = {FX(Rc)}; S_cut/R_cut = {fac}"
                  f"{'; R printed ' + str(printed[1]) + (' EQUAL' if hit_R else ' differs') if printed[1] else ''}; factor printed {printed[3]} "
                  f"{'EQUAL' if hit_f else 'differs'}")
    port_only = [k for k, v in out.items() if v and k[1] == "C1"]
    pairs = [(S, Fr(23328 * S, 81 * 81)) for S in range(1, 28) if Fr(23328 * S, 81 * 81).denominator == 1 and 23328 * S // (81 * 81) <= 81]
    print(f"  gluon word: (cut chain sum S, port factor c) with S/(81*81*c) = 1/23328, S <= 27, c <= 81: "
          f"{[(S, int(c)) for S, c in pairs]}; c = 64 is the (3,3) cut face, 32 is no face of 103.1")
    print(f"  the port-only cut reproduces the printed cut values at {len(port_only)} of 2 bodies (the gluon word's cut factor with the cut "
          f"chain sum 27 - 9 = 18, the 147.1 pattern); every-link and every-(2,3)-link cuts reproduce none")
    assert [(S, int(c)) for S, c in pairs] == [(9, 32), (18, 64)]
    assert len(port_only) == 2 and not any(v for k, v in out.items() if k[1] != "C1")
    print()


# ------------------------------------------------------------------ ELECTRON
def electron_block(U):
    print("[ELECTRON] the electron web parsed from the 152.1 figure; its body registers under each reading")
    dims, links, genres = MC.parse_figure(U["152.1"])
    assert all(tuple(sorted((dims[a], dims[b]))) == g for (a, b), g in zip(links, genres))
    docks = [l for l in links if {dims[l[0]], dims[l[1]]} == {2, 3}]
    port = next(l for l in links if set(l) == {"B", "w1"})
    deg = {s: sum(s in l for l in links) for s in dims}
    m = measure(dims, links); strike = sum(dims[a] * dims[b] for a, b in docks); rank = m - strike
    inlegs = [l for l in links if {dims[l[0]], dims[l[1]]} == {3, 4}]
    hub = math.prod(dims[a] if dims[a] == 3 else dims[b] for a, b in inlegs)
    print(f"  sites {len(dims)} {dims}; links {len(links)}; genres {genres}")
    print(f"  dock links (2,3) {docks}; port B-w1 = {port}, genre {tuple(sorted((dims[port[0]], dims[port[1]])))}; degrees {deg}")
    print(f"  measure {m}, dock strike {strike}, rank {rank} (dock unit), hub product {hub}, hub unit {rank * hub}  (152.1: 79, 12, 67; 201)")
    assert (len(dims), len(links), m, strike, rank, rank * hub) == (9, 10, 79, 12, 67, 201)
    Rw = register(dims, links, lambda l: "whole")
    via_deg = math.prod(dims[s] ** (2 * deg[s]) for s in dims)
    print(f"  per-link whole faces {[faces(dims[a], dims[b])['whole'] for a, b in links]}")
    print(f"  R_whole = prod over the 10 links of (d_a d_b)^2 = {Rw} = {FX(Rw)}  (= prod over sites of d^(2 deg) = {FX(via_deg)}); "
          f"printed 2^30*3^18: {'EQUAL' if Rw == 2**30 * 3**18 else 'differs'}")
    assert Rw == via_deg == 2**30 * 3**18
    R, NAME = {}, {}
    NAME["W0"] = "whole, all 10 links (the struck docks counted)"; R["W0"] = Rw
    NAME["W1"] = "whole, the 8 unstruck links only"; R["W1"] = register(dims, [l for l in links if l not in docks], lambda l: "whole")
    for code, rule in cut_readings(dims, links, port, docks).items():
        NAME[code] = CUT[code]; R[code] = register(dims, links, rule)
    NAME["C4"] = "the port at its support face (16), every other link whole"
    R["C4"] = register(dims, links, lambda l: "support" if l == port else "whole")
    NAME["C5"] = "the port at its rank face (6), every other link whole"
    R["C5"] = register(dims, links, lambda l: "rank" if l == port else "whole")
    NAME["F1"] = "every link at its rank face"; R["F1"] = register(dims, links, lambda l: "rank")
    NAME["F2"] = "every link at its support face"; R["F2"] = register(dims, links, lambda l: "support")
    print("  register readings (printed: whole 2^30*3^18, cut 2^31*3^17):")
    for k, v in R.items():
        tag = "= printed whole" if v == 2**30 * 3**18 else ("= printed cut" if v == 2**31 * 3**17 else "")
        print(f"    {k} {NAME[k]:68s} {FX(v):16s} {tag}")
    hits_w = [k for k, v in R.items() if v == 2**30 * 3**18]
    hits_c = [k for k, v in R.items() if v == 2**31 * 3**17]
    print(f"  readings that give the printed whole register: {hits_w}; the printed cut register: {hits_c}")
    print(f"  coincidence: C3 and C4 give the same number {FX(R['C4'])}, since (24/36)^2 = 16/36 = 4/9")
    assert hits_w == ["W0"] and hits_c == ["C1"] and R["C3"] == R["C4"]
    print()
    return dict(dims=dims, links=links, port=port, docks=docks, R=R, rank=rank, hub=rank * hub)


# ------------------------------------------------------------------ CARRIER
def carrier_block(U):
    print("[CARRIER] the photon word parsed from the 143.1 figure; the mass face M_L = 36 (printed without a formula)")
    dims, links, genres = MC.parse_figure(U["143.1"])
    m = measure(dims, links)
    print(f"  photon word: {len(dims)} sites {list(dims.values())}, {len(links)} links, genres {set(genres)}, measure {m}")
    assert (len(dims), len(links), m, set(genres)) == (6, 5, 30, {(2, 3)})
    w23 = faces(2, 3)["whole"]
    rows = [("measure of the held chain as displayed (6 sites, 5 links)", 5 * 6),
            ("measure after the docking act: 5 chain + 2 made - 1 contracted = 6 links of (2,3)", (5 + 2 - 1) * 6),
            ("measure of 'the eight links of the held chain' (488, 511)", 8 * 6),
            ("whole register of one (2,3) link (16.1, 103.1; 192.1 'the whole rung register as thirty-six')", w23)]
    for t, v in rows:
        print(f"    {t:96s} {v:3d} {'= 36' if v == 36 else ''}")
    print(f"  exchanged mass face 6 sqrt6: (6 sqrt6)^2 = {6**3} = 6^3, M_exch / M_L = 1/sqrt6; no formula printed for either")
    print(f"  link counts the text gives the held chain or the bond: 5 (display), 8 (fine-structure bullet, census), 1 (glossary: "
          f"'the link that joins the electron body to the proton body')\n")


# ------------------------------------------------------------------ E1
def e1_block(E):
    print("[E1] the level-one binding E1 = M_L * u / R ('the mass face multiplied by the unit over the register of the body at that face')")
    M = 36
    R = E["R"]
    regs = [("W0 whole", R["W0"]), ("C1 port cut", R["C1"]), ("C2 every link cut", R["C2"]),
            ("C3/C4 docks cut / port support", R["C3"]), ("C5 port rank", R["C5"]), ("W1 unstruck whole", R["W1"]),
            ("port register 36", 36), ("port register 24", 24), ("port register 16 (229.1)", 16)]
    printed = {("W0 whole", 201): Fr(67, 2**28 * 3**15), ("C1 port cut", 67): Fr(67, 2**29 * 3**15),
               ("port register 36", 67): Fr(67), ("port register 24", 67): Fr(201, 2)}
    print(f"  {'register reading':30s} | u = 67 (dock unit)                 | u = 201 (hub unit)")
    found = 0
    for rn, rv in regs:
        cells = []
        for u in (67, 201):
            v = Fr(M * u, rv)
            p = printed.get((rn, u))
            if p is not None:
                assert v == p; found += 1
            cells.append(f"{frac_str(v):26s}{' PRINTED' if p is not None else '        '}")
        print(f"  {rn:30s} | {cells[0]} | {cells[1]}")
    print(f"  printed equalities reproduced: {found} of 4 (two faces; two port-register readings)")
    w = Fr(M * 201, regs[0][1]); c = Fr(M * 67, regs[1][1])
    print(f"  whole / cut = {w / c} = (201/67) * (24/36) = {Fr(201, 67)} * {Fr(24, 36)}: the hub dressing 3 times the (2,3) cut/whole ratio 2/3")
    print(f"  pairing: the text pairs (whole register, hub unit 201) and (cut register, dock unit 67); the cross pairings give "
          f"{frac_str(Fr(M * 67, regs[0][1]))} and {frac_str(Fr(M * 201, regs[1][1]))}")
    print(f"  229.1 resolves the unit to the dock at the register 16: 36*67/16 = {Fr(M * 67, 16)}; with the port at its support face in the body "
          f"register: {frac_str(Fr(M * 67, R['C4']))}\n")
    assert found == 4 and w / c == 2
    return {"whole": w, "cut": c, "port36": Fr(67), "port24": Fr(201, 2)}


# ------------------------------------------------------------------ LADDER
LADDER_RX = re.compile(r"held pendant pair|held pair|two-vertex ladder|per-record ladder|per record read|level-\$n\$|E_\{n\}|Rydberg ladder|"
                       r"over \$n\$ squared|rank squared|eight ranks", re.I)
NUMWORD = {"nine": 9, "ten": 10, "eleven": 11, "twelve": 12, "sixteen": 16, "twenty": 20}


def sl2_commutant_is_scalar(n):
    """Sym^(n-1)(C^2) with basis x^a y^(n-1-a), a = n-1..0: H = diag(a-b), E = x d/dy (integer entries b at (a+1, a)).
    A matrix commuting with H is diagonal iff the H eigenvalues are distinct; a diagonal matrix commuting with E is scalar iff
    the support graph of E connects all basis vectors. Both are checked exactly; then the commutant is C*I."""
    k = n - 1
    basis = [(a, k - a) for a in range(k, -1, -1)]
    hdiag = [a - b for a, b in basis]
    distinct = len(set(hdiag)) == len(hdiag)
    idx = {v: i for i, v in enumerate(basis)}
    edges = [(idx[(a + 1, b - 1)], idx[(a, b)]) for a, b in basis if b > 0]      # E x^a y^b = b x^(a+1) y^(b-1), b != 0
    parent = list(range(len(basis)))
    def find(x):
        while parent[x] != x: parent[x] = parent[parent[x]]; x = parent[x]
        return x
    for i, j in edges: parent[find(i)] = find(j)
    connected = len({find(i) for i in range(len(basis))}) == 1
    return distinct and connected


def sl2_commutant_bruteforce(n):
    k = n - 1
    basis = [(a, k - a) for a in range(k, -1, -1)]
    idx = {v: i for i, v in enumerate(basis)}
    H = sp.diag(*[a - b for a, b in basis]) if n > 1 else sp.Matrix([[0]])
    Em = sp.zeros(n, n); Fm = sp.zeros(n, n)
    for a, b in basis:
        if b > 0: Em[idx[(a + 1, b - 1)], idx[(a, b)]] = b
        if a > 0: Fm[idx[(a - 1, b + 1)], idx[(a, b)]] = a
    xs = sp.symbols(f"x0:{n*n}")
    Xm = sp.Matrix(n, n, xs)
    eqs = [e for A in (H, Em, Fm) for e in (A * Xm - Xm * A)]
    Mx = sp.Matrix([[sp.diff(e, x) for x in xs] for e in eqs])
    return n * n - Mx.rank()


def is_pure_power(q, n):
    q = Fr(q)
    if n == 1: return q == 1
    if q <= 0: return False
    k0 = round(math.log(q.numerator / q.denominator) / math.log(n))
    return any(Fr(n) ** k == q for k in (k0 - 1, k0, k0 + 1))


def ladder_block(src):
    print("[LADDER] (b) the source's ladder beyond n = 8, and reconstruction R")
    print("  (1) every sentence of the tex that states the hydrogen ladder (mechanical scan, regex below), with the range it states")
    print(f"      regex: {LADDER_RX.pattern}")
    rows = []
    for f in sorted((src / "tex").rglob("*.tex")):
        for ln, line in enumerate(f.read_text(encoding="utf-8").splitlines(), 1):
            if not LADDER_RX.search(line): continue
            for s in re.split(r"(?<=[.;])\s+", line):
                if LADDER_RX.search(s):
                    eight = bool(re.search(r"\beight\b|1, \\dots, 8|one to eight", s))
                    beyond = [w for w in NUMWORD if re.search(rf"\b{w}\b", s)] + re.findall(r"n\s*=\s*(?:9|1\d)", s)
                    rows.append((f"{f.relative_to(src)}:{ln}", eight, beyond, s))
    hyd = [r for r in rows if "pendant pair" in r[3] or "held pair" in r[3] or "ladder" in r[3].lower() or "level-$n$" in r[3] or "E_{n}" in r[3]]
    for loc, eight, beyond, s in hyd:
        print(f"      {loc:44s} range 'eight'/'1..8': {'yes' if eight else 'no '}  other numerals: {beyond or '-'}  | {X.norm(s)[:150]}")
    n8 = sum(1 for r in hyd if r[1]); nb = sum(1 for r in hyd if r[2])
    print(f"      {len(hyd)} sentences; {n8} state the range eight / 1..8; {nb} name any rank beyond eight; no ladder value or formula "
          f"beyond n = 8 is printed; the per-record ladder has no printed value at any n")
    lean_hits = []
    for p in sorted((src / "lean").rglob("*.lean")):
        for ln, line in enumerate(p.read_text(encoding="utf-8", errors="replace").splitlines(), 1):
            if re.search(r"RYDBERG|Rydberg|HYDROGEN|hydrogen", line): lean_hits.append(f"{p.relative_to(src)}:{ln}")
    idx_hits = [ln for ln, line in enumerate((src / "index" / "statements.tsv").read_text(encoding="utf-8").splitlines(), 1)
                if re.search(r"Rydberg ladder", line)]
    print(f"      Lean tree: {len(lean_hits)} lines name the Rydberg ladder or hydrogen ({', '.join(lean_hits)}), as record names only; "
          f"index statements.tsv: {len(idx_hits)} rows with 'Rydberg ladder' (section titles)")
    assert nb == 0 and n8 >= 5

    N = 200
    print(f"  (2) R exactly: twirled state on Sym^(n-1)(C^2) = I/n, single read 1/n, two-vertex read (1/n)^2, at n = 1..{N}")
    ok = all(sl2_commutant_is_scalar(n) for n in range(1, N + 1))
    bf = [sl2_commutant_bruteforce(n) for n in range(1, 6)]
    print(f"      commutant of the sl2 action is C*I at {sum(sl2_commutant_is_scalar(n) for n in range(1, N + 1))} of {N} ranks "
          f"(H eigenvalues distinct, E support connected); brute-force commutant dimension over Q at n = 1..5: {bf}")
    print(f"      so twirl(rho) = tr(rho) I/n (Schur) and the two-vertex ladder is exactly 1/n^2 at every n <= {N}; the argument holds at every n")
    assert ok and bf == [1] * 5
    print("  (3) the support dimension m: dim Sym^(n-1)(C^m) = C(n+m-2, m-1); two-vertex ladder 1/dim^2 against 1/n^2 at n = 1..200")
    nn = sp.symbols("n", positive=True, integer=True)
    for mm in (2, 3, 4, 6):
        eq = sum(1 for n in range(1, N + 1) if sp.binomial(n + mm - 2, mm - 1) == n)
        dim = sp.factor(sp.expand_func(sp.binomial(nn + mm - 2, mm - 1)))
        sol = sp.solve(sp.Eq(sp.expand_func(sp.binomial(nn + mm - 2, mm - 1)), nn), nn)
        print(f"      m = {mm}: dim = {dim}; equal to 1/n^2 at {eq} of {N} ranks; dim = n solved over positive integers: "
              f"{'all n' if mm == 2 else sorted(sol)}")
        assert eq == (N if mm == 2 else 1)
    print("  (4) order readings: under R a term that reads the twirled state r times is n^(-r) exactly (r reads of I/n at a unit vector)")
    for r in (0, 1, 2):
        vals = [Fr(1, n ** r) for n in range(1, 9)]
        print(f"      r = {r}: n = 1..8 -> {[str(v) for v in vals]}; equal to 1/n^2 at {sum(v == Fr(1, n*n) for n, v in zip(range(1, 9), vals))} of 8")
    print("      45.1 (prose): 'the order-k term at rank n the term at rank one over n squared at all eight ranks' -- under R this holds for")
    print("      every k only if every order-k term at the held pair reads the state twice; the source does not define these terms")
    print("  (5) per-record readings under R (143.1: 'the per-record ladder is exact at every rank and no pure power of the rank')")
    P1 = [Fr(1, n) for n in range(1, N + 1)]
    P2 = []
    for n in range(1, N + 1):
        j2 = Fr(n - 1, 2)
        ms = [j2 - i for i in range(n)]
        P2.append(sum(mv * mv for mv in ms) / n)          # <J_z^2> at I/n
    assert all(P2[n - 1] == Fr(n * n - 1, 12) for n in range(1, N + 1))
    p1 = sum(is_pure_power(v, n) for n, v in enumerate(P1, 1)); p2 = sum(is_pure_power(v, n) for n, v in enumerate(P2, 1))
    p2n = [n for n, v in enumerate(P2, 1) if is_pure_power(v, n)]
    print(f"      P1, one read of the twirled state, 1/n: a pure power of n at {p1} of {N} ranks (parts from 'at none' at every rank)")
    print(f"      P2, the second moment <J_z^2> of one traceless letter at I/n, (n^2-1)/12 (checked exactly): a pure power of n at {p2} of {N}, "
          f"n = {p2n} ((4-1)/12 = 1/4 = 2^-2), so it parts from 'at none' at that rank only")
    print("      P0, the first moment of a traceless letter, is 0 at every n: never a power (agrees, trivially)")
    print("      so the per-record sentence decides R only once 'record' is defined; it is not defined for the held pair\n")
    assert p1 == N and p2n == [2]


# ------------------------------------------------------------------ COMPARE
def compare_block(E1):
    print("[COMPARE] (c) the source's hydrogen numbers against experiment")
    PC = sc.physical_constants
    alpha = PC["fine-structure constant"][0]; Ryd_eV = PC["Rydberg constant times hc in eV"][0]
    mec2_eV = PC["electron mass energy equivalent in MeV"][0] * 1e6
    assert abs(alpha - 7.2973525643e-3) < 1e-16 and abs(Ryd_eV - 13.605693122990) < 1e-11 and abs(mec2_eV - 510998.95069) < 1e-4
    IE_cm = 109678.77174307
    IE_eV = IE_cm * 100 * sc.h * sc.c / sc.e
    assert abs(IE_eV - 13.598434599702) < 1e-11
    t_inf = Ryd_eV / mec2_eV; t_H = IE_eV / mec2_eV
    print(f"  identification as printed: E1 'in the electron unit' (box), 'the Rydberg constant in the electron unit the ladder cell of the")
    print(f"  held carrier' (bullet), 'The correspondence is by name' (bridge); 192.1: 'the mass unit stays a question of which cell is")
    print(f"  meant'; back.tex: 'The Rydberg constant as the two vertex term at the bond' rests on 'Mass unit as the rung register (an open")
    print(f"  question)'. The electron unit is mapped to no physical quantity; 176.1 uses 'the proton in the electron unit' for m_p/m_e.")
    print(f"  reading 'X in the electron unit' = X / m_e c^2 (as 176.1 uses it): targets R_inf hc/m_e c^2 = alpha^2/2 = {t_inf:.7e} "
          f"(alpha^2/2 = {alpha**2/2:.7e}); measured IE(H)/m_e c^2 = {t_H:.7e}")
    gaps = {}
    for k in ("whole", "cut"):
        v = float(E1[k]); gaps[k] = math.log10(t_H / v)
        print(f"    E1 {k:5s} = {frac_str(E1[k]):18s} = {v:.6e}: / (R_inf/m_e) = {v/t_inf:.4e}, / (IE/m_e) = {v/t_H:.4e}, "
              f"log10 gap {math.log10(t_H/v):.3f}")
    print("  declared conversion family (the electron unit as a cell is open in the source): the four printed binding numbers")
    print("  {E1 whole, E1 cut, 67, 201/2} divided by u^p, u in {1, 16, 36, 67, 201} (the five 192.1 spellings and the hub unit), p in {1, 2},")
    print("  or undivided: 4 x 9 = 36 variants, compared with IE/m_e c^2")
    conv = [(1, 0)] + [(u, p) for u in (16, 36, 67, 201) for p in (1, 2)]
    var = [(name, u, p, float(Fr(v) / Fr(u) ** p)) for name, v in E1.items() for u, p in conv]
    var.sort(key=lambda t: abs(math.log(t[3] / t_H)))
    print(f"    variants {len(var)}; within x2: {sum(abs(math.log10(t[3]/t_H)) < math.log10(2) for t in var)}, within x10: "
          f"{sum(abs(math.log10(t[3]/t_H)) < 1 for t in var)}, within x100: {sum(abs(math.log10(t[3]/t_H)) < 2 for t in var)}")
    for name, u, p, v in var[:3]:
        print(f"    nearest: {name} / {u}^{p} = {v:.4e}, ratio to IE/m_e c^2 = {v/t_H:.3g}")
    print(f"    log10 range of the family: {min(math.log10(t[3]) for t in var):.2f} to {max(math.log10(t[3]) for t in var):.2f}")
    # the unit-free ladder ratio
    nu12 = 2466061413187035.0                # Hz, 1S-2S hyperfine centroid, Parthey et al. 2011
    B1 = IE_eV * sc.e / sc.h                 # Hz, binding of 1S (NIST ionization energy)
    B2 = B1 - nu12
    d = 4 * B2 / B1 - 1
    print(f"  E_2/E_1 = 1/4 (143.1; the electron unit cancels) against the measured B(2S)/B(1S) = 1 - nu(1S-2S)/IE = {B2/B1:.12f}")
    print(f"    4 B(2S)/B(1S) - 1 = {d:+.4e}; alpha^2 = {alpha**2:.4e}; the Dirac energies alone give alpha^2/16 = {alpha**2/16:.4e}, "
          f"the remainder {d - alpha**2/16:+.2e} comes from the QED and recoil terms beyond them (not evaluated here)")
    sIE, snu = 1.0e-7 / IE_cm, 10.0 / nu12          # relative: IE (10) in the last digits of the cm^-1 value; nu 10 Hz
    sd = 4 * (nu12 / B1) * math.hypot(sIE, snu)
    print(f"    measurement uncertainty of 4 B(2S)/B(1S) - 1: {sd:.1e} (IE {sIE:.1e} relative, nu {snu:.1e} relative)")
    print("  the residual is of the order of the fine structure and Lamb shift, which 143.1 states as further terms without values")
    print("  (Lamb sign and hyperfine 0 vs 1420.405751768 MHz: see hydrogen_audit.py)\n")
    within10 = sum(abs(math.log10(t[3] / t_H)) < 1 for t in var)
    return dict(gaps=gaps, within10=within10, nvar=len(var), nearest=var[0][3] / t_H, delta=d)


def main():
    src = pathlib.Path(sys.argv[1])
    U = anchors(src)
    faces_block()
    bodies_block()
    E = electron_block(U)
    carrier_block(U)
    E1 = e1_block(E)
    ladder_block(src)
    C = compare_block(E1)
    print("[SUMMARY]")
    print("  (a) 2^30*3^18 = prod over the 10 electron links of (d_a d_b)^2 (103.1 product rule, 16.1/103.1 whole face), given the 152.1")
    print("      web and the struck docks counted; 2^31*3^17 = the same with the port B-w1 at its cut 24 = (2^2-1)(3^2-1), the convention")
    print("      that also reproduces the printed Z-body and gluon-word cut values. Units 67, 201 from the web; mass face 36 entered.")
    print("  (b) no ladder value or formula beyond n = 8 is printed; R reproduces the 8 printed values and gives 1/n^2 exactly at every n;")
    print("      it parts from the per-record sentence at reading P1 (every n) and P2 (n = 2), and from the all-orders sentence of 45.1 for")
    print("      terms that read the state other than twice; the source defines neither 'record' nor 'order' at the held pair.")
    print(f"  (c) the unit map is open in the source (192.1); under X/m_e c^2 E1 is {C['gaps']['whole']:.2f} / {C['gaps']['cut']:.2f} orders low; "
          f"{C['within10']} of {C['nvar']} declared")
    print(f"      conversions within x10 (nearest x{C['nearest']:.0f}); E_2/E_1 = 1/4 agrees with the measured ratio to {C['delta']:.1e}, the size of the")
    print("      terms the source does not evaluate.")


if __name__ == "__main__":
    main()
