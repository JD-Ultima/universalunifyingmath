"""TH-001/TH-002/TH-004: which subsets of {H1,H2,LD,RD,AS,CM,INC} are jointly satisfiable over R?

The seven rules (INC = "multiplication makes more"; EFF-strong and BAL are NOT among them, see Part 2b):
  H1  1(x)1 = 2            H2  2(x)2 = 4
  LD  a(x)(b+c) = a(x)b + a(x)c          RD  (a+b)(x)c = a(x)c + b(x)c
  AS  associativity        CM  commutativity
  INC a(x)b > max(a,b) for all positive a, b

Part 1. Ground certificates. For an UNINTERPRETED operation op : R x R -> R, z3 is given finitely many
        ground instances of each rule. UNSAT of the ground instances of a subset is a proof that the
        subset (with its rules universally quantified) is inconsistent. Nine cores are certified this
        way: the five H1+H2 cores, and H1+LD+INC, H1+RD+INC, H2+LD+INC, H2+RD+INC (H1 or H2 fixes a
        value, e.g. 1(x)1/4 = 1/2 or 2(x)1/4 = 1/2, that INC forbids). The two Archimedean cores LD+INC
        and RD+INC have NO such certificate: {LD,RD,AS,CM,INC} has a model over the non-Archimedean
        real closed field R((t^Q)) (formalization/howard/arithmetic.md, A6), and a finite set of ground
        instances true there is satisfiable over R. They are certified by Theorem A6 instead (Lean: archimedean_ld_generic / archimedean_rd_generic, instantiated for Q
        as core_ld_inc_rat / core_rd_inc_rat in proofs/lean/howard/HowardArith.lean; for R the
        generic hypotheses are the standard ordered-field and Archimedean properties of R).
Part 2. Every subset of the 7 rules is classified: inconsistent if it contains a core, consistent if
        one of the explicit witnesses (checked universally in arith_candidates.py / HowardArith.lean,
        2a below) satisfies it, otherwise 'open'. Counts are reported per certificate type.
Part 2b. The same with BAL (a(x)b >= a+b for positive a, b) added as an eighth rule.
Part 3. Finite models: all bi-additive operations on Z/n (n <= 12) satisfying H1, H2.
        Brute force over ALL n^(n^2) operations for n <= 3 confirms that bi-additive ops on Z/n
        are exactly a(x)b = k*a*b.
"""
import itertools
import z3

op = z3.Function('op', z3.RealSort(), z3.RealSort(), z3.RealSort())
R = z3.RealVal

def LD(p, q, r): return op(p, q + r) == op(p, q) + op(p, r)
def RD(p, q, r): return op(p + q, r) == op(p, r) + op(q, r)

GROUND = {
    "H1": [op(1, 1) == 2],
    "H2": [op(2, 2) == 4],
    "LD": [LD(1, 1, 1), LD(2, 1, 1), LD(1, R(1) / 2, R(1) / 2), LD(1, R(1) / 4, R(1) / 4),
           LD(2, R(1) / 2, R(1) / 2), LD(2, R(1) / 4, R(1) / 4)],
    "RD": [RD(1, 1, 1), RD(1, 1, 2), RD(R(1) / 2, R(1) / 2, 1), RD(R(1) / 4, R(1) / 4, 1),
           RD(R(1) / 2, R(1) / 2, 2), RD(R(1) / 4, R(1) / 4, 2)],
    "AS": [op(op(1, 1), 1) == op(1, op(1, 1))],
    "CM": [op(2, 1) == op(1, 2)],
    # INC instances used: 1 (x) 1/4 > max(1, 1/4) = 1, 2 (x) 1/4 > 2, and their mirrors
    "INC": [op(1, R(1) / 4) > 1, op(R(1) / 4, 1) > 1, op(2, R(1) / 4) > 2, op(R(1) / 4, 2) > 2],
    # BAL instances (Part 2b only): a (x) b >= a + b at the same four points
    "BAL": [op(1, R(1) / 4) >= R(5) / 4, op(R(1) / 4, 1) >= R(5) / 4,
            op(2, R(1) / 4) >= R(9) / 4, op(R(1) / 4, 2) >= R(9) / 4],
}

CORES = [("H1", "H2", "LD", "RD"), ("H1", "H2", "LD", "AS"), ("H1", "H2", "LD", "CM"),
         ("H1", "H2", "RD", "AS"), ("H1", "H2", "RD", "CM"), ("H1", "LD", "INC"), ("H1", "RD", "INC"),
         ("H2", "LD", "INC"), ("H2", "RD", "INC")]

def check_core(core):
    s = z3.Solver()
    for k in core:
        s.add(*GROUND[k])
    return s.check()

# Witness operations over R (properties certified universally by arith_candidates.py, z3 nlsat)
WITNESSES = {
    "a+b":       {"H1", "H2", "AS", "CM", "INC", "BAL"},
    "2ab":       {"H1", "LD", "RD", "AS", "CM"},
    "2b":        {"H1", "H2", "LD"},
    "a(b+1)":    {"H1", "RD"},
    "ab":        {"H2", "LD", "RD", "AS", "CM"},
    "a+b+2":     set(),  # placeholder, not used
}
# a mirror of "2b" (a(x)b = 2a) gives {H1,H2,RD}; add it explicitly (z3-verified below)
WITNESSES["2a"] = {"H1", "H2", "RD"}
del WITNESSES["a+b+2"]

def verify_2a():
    A, B, C = z3.Reals('a b c')
    f = lambda p, q: 2 * p
    s = z3.Solver(); s.add(z3.Not(f(A + B, C) == f(A, C) + f(B, C)))
    ok_rd = s.check() == z3.unsat
    return ok_rd and f(1, 1) == 2 and f(2, 2) == 4

AX = ["H1", "H2", "LD", "RD", "AS", "CM", "INC"]

# Archimedean cores: Theorem A6 (formalization/howard/arithmetic.md). Lean: archimedean_ld_generic,
# archimedean_rd_generic (no axioms), instances core_ld_inc_rat, core_rd_inc_rat over Q. Not z3: no finite
# set of ground instances is UNSAT ({LD,RD,AS,CM,INC} has a model over the real closed field R((t^Q))).
ARCH_CORES = [("LD", "INC"), ("RD", "INC")]

def classify(ax, z3_cores, arch_cores, witnesses):
    """Return {'z3': [...], 'A6': [...], 'consistent': [...], 'open': [...]} over all subsets of ax.
    A subset counts as 'z3' if it contains a z3-certified core, else 'A6' if it contains an Archimedean core."""
    stats = {"z3": [], "A6": [], "consistent": [], "open": []}
    for r in range(len(ax) + 1):
        for sub in itertools.combinations(ax, r):
            S = set(sub)
            if any(set(c) <= S for c in z3_cores):
                stats["z3"].append(sub)
            elif any(set(c) <= S for c in arch_cores):
                stats["A6"].append(sub)
            elif any(S <= w for w in witnesses):
                stats["consistent"].append(sub)
            else:
                stats["open"].append(sub)
    return stats

def allsubs7(ax):
    return [c for r in range(len(ax) + 1) for c in itertools.combinations(ax, r)]

def minimal(core, witnesses):
    return all(any(set(core) - {k} <= w for w in witnesses) for k in core)

def zn_models(nmax=12):
    out = []
    for n in range(1, nmax + 1):
        for k in range(n):
            f = lambda p, q: (k * p * q) % n
            if f(1, 1) == 2 % n and f(2, 2) == 4 % n:
                out.append((n, k, (2 % n) != 0))
    return out

def brute_biadditive(n):
    """Enumerate all n^(n*n) tables on Z/n; return the bi-additive ones as k-values (or None)."""
    elems = range(n)
    found = []
    for table in itertools.product(elems, repeat=n * n):
        f = lambda p, q: table[p * n + q]
        if all(f(p, (q + r) % n) == (f(p, q) + f(p, r)) % n and f((p + q) % n, r) == (f(p, r) + f(q, r)) % n
               for p in elems for q in elems for r in elems):
            k = f(1, 1)
            found.append(k if all(f(p, q) == (k * p * q) % n for p in elems for q in elems) else None)
    return found

if __name__ == "__main__":
    print("z3", z3.get_version_string())
    W7 = [w - {"BAL"} for w in WITNESSES.values()]
    print("\nPart 1: z3 ground certificates (UNSAT of the subset's own ground instances)")
    for core in CORES:
        print("  ", "+".join(core).ljust(18), check_core(core))
    print("  minimality (every proper subset satisfied by a witness):")
    for core in CORES + ARCH_CORES:
        tag = "minimal" if minimal(core, W7) else "not minimal"
        if len(core) == 3:
            tag += (f": contains the Archimedean core {'+'.join(core[1:])} (Thm A6);"
                    f" {core[0]} makes it finitely certifiable")
        src = "z3" if core in CORES else "A6 (Lean generic; Q instance)"
        print("  ", "+".join(core).ljust(18), src.ljust(30), tag)
    print("  witness a(x)b = 2a satisfies RD, H1, H2 (z3):", verify_2a())
    print("\nPart 2: all 2^7 = 128 subsets of", AX)
    st = classify(AX, CORES, ARCH_CORES, W7)
    print(f"  inconsistent: {len(st['z3']) + len(st['A6'])}"
          f"  (z3 ground certificate: {len(st['z3'])}; Archimedean Thm A6 only: {len(st['A6'])})")
    print(f"  consistent: {len(st['consistent'])}")
    print(f"  open: {len(st['open'])}")
    print("  subsets certified by A6 only:")
    print("    ", st["A6"])
    expect = [s for s in allsubs7(AX) if set(s) <= {"LD", "RD", "AS", "CM", "INC"} and "INC" in s
              and ({"LD", "RD"} & set(s))]
    print("  check: A6-only subsets = subsets of {LD,RD,AS,CM,INC} containing INC and LD or RD:",
          sorted(map(sorted, st["A6"])) == sorted(map(sorted, expect)))
    print("  maximal consistent subsets:")
    cons = [set(s) for s in st["consistent"]]
    maxi = [s for s in cons if not any(s < t for t in cons)]
    for s in maxi:
        print("    ", sorted(s, key=AX.index))
    allsubs = [set(c) for r in range(len(AX) + 1) for c in itertools.combinations(AX, r)]
    free = [x for x in allsubs if not any(set(k) <= x for k in CORES + ARCH_CORES)]
    maxfree = [x for x in free if not any(x < y for y in free)]
    print("  maximal consistent = maximal core-free:", sorted(map(sorted, maxi)) == sorted(map(sorted, maxfree)))
    print("  sanity: ground instances of each maximal consistent set are SAT (z3):",
          [str(check_core(sorted(m, key=AX.index))) for m in maxi])
    print("\nPart 2b: BAL added as an eighth rule (2^8 = 256 subsets)")
    AX8 = AX + ["BAL"]
    EXTRA8 = [("H1", "LD", "BAL"), ("H1", "RD", "BAL"), ("H2", "LD", "BAL"), ("H2", "RD", "BAL")]
    CORES8 = CORES + EXTRA8
    for c in EXTRA8:
        print("  ", "+".join(c).ljust(18), check_core(c))
    A_BAL = ARCH_CORES + [("LD", "BAL"), ("RD", "BAL")]   # BAL implies a(x)b >= a and >= b: A6 weak form
    st8 = classify(AX8, CORES8, A_BAL, list(WITNESSES.values()))
    print(f"  inconsistent: {len(st8['z3']) + len(st8['A6'])} (z3: {len(st8['z3'])}; A6 only: {len(st8['A6'])});"
          f" consistent: {len(st8['consistent'])}; open: {len(st8['open'])}")
    same = all(((s in st8["consistent"]) ==
                (tuple(sorted((set(s) - {"BAL"}) | ({"INC"} if "BAL" in s else set()), key=AX.index))
                 in st["consistent"])) for s in [x for r in range(9) for x in itertools.combinations(AX8, r)])
    print("  every subset has the verdict of the subset with BAL replaced by INC:", same)
    print("\nPart 3: bi-additive a(x)b = k*a*b on Z/n with H1 and H2, n <= 12  -> (n, k, 2 != 0 in Z/n)")
    print("  ", zn_models())
    for n in (2, 3):
        res = brute_biadditive(n)
        print(f"   brute force Z/{n}: {n**(n*n)} tables, {len(res)} bi-additive, k-values {sorted(res)}"
              f" (None would mean a non-k*a*b solution)")
