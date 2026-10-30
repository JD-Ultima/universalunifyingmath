"""TH-001/TH-002: which subsets of {H1,H2,LD,RD,AS,CM,INC} are jointly satisfiable?

Part 1. For an UNINTERPRETED operation op : R x R -> R, z3 is given only finitely many ground
        instances of each axiom (enough to derive the contradiction). UNSAT on ground instances
        is a proof that the full (universal) axiom set is inconsistent.
Part 2. Every subset of the 7 axioms is classified: inconsistent if it contains one of the
        certified cores, consistent if one of the explicit witnesses (checked universally in
        arith_candidates.py / HowardArith.lean) satisfies it, otherwise 'open'.
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
    "LD": [LD(1, 1, 1), LD(2, 1, 1), LD(1, R(1) / 2, R(1) / 2), LD(1, R(1) / 4, R(1) / 4)],
    "RD": [RD(1, 1, 1), RD(1, 1, 2)],
    "AS": [op(op(1, 1), 1) == op(1, op(1, 1))],
    "CM": [op(2, 1) == op(1, 2)],
    # INC instance used: 1 (x) 1/4 > max(1, 1/4) = 1   (only meaningful on a carrier containing 1/4)
    "INC": [op(1, R(1) / 4) > 1],
}

CORES = [("H1", "H2", "LD", "RD"), ("H1", "H2", "LD", "AS"), ("H1", "H2", "LD", "CM"),
         ("H1", "H2", "RD", "AS"), ("H1", "H2", "RD", "CM"), ("H1", "LD", "INC")]

def check_core(core):
    s = z3.Solver()
    for k in core:
        s.add(*GROUND[k])
    return s.check()

# Witness operations over R (properties certified universally by arith_candidates.py, z3 nlsat)
WITNESSES = {
    "a+b":       {"H1", "H2", "AS", "CM", "INC"},
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

# Cores proved on paper only (formalization/howard/arithmetic.md, Thm A6: an additive map that is
# positive on positives is linear, then Archimedes). No finite set of ground instances suffices.
PAPER_CORES = [("LD", "INC"), ("RD", "INC")]

def classify(use_paper=False):
    stats = {"inconsistent": [], "consistent": [], "open": []}
    cores = CORES + (PAPER_CORES if use_paper else [])
    for r in range(len(AX) + 1):
        for sub in itertools.combinations(AX, r):
            S = set(sub)
            if any(set(c) <= S for c in cores):
                stats["inconsistent"].append(sub)
            elif any(S <= w for w in WITNESSES.values()):
                stats["consistent"].append(sub)
            else:
                stats["open"].append(sub)
    return stats

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
    print("\nPart 1: certified inconsistent cores (ground instances only)")
    for core in CORES:
        print("  ", "+".join(core).ljust(18), check_core(core))
    print("  minimality: every proper subset of each core is satisfied by a witness:")
    for core in CORES:
        subs = [set(core) - {k} for k in core]
        ok = all(any(sb <= w for w in WITNESSES.values()) for sb in subs)
        print("  ", "+".join(core).ljust(18), "minimal" if ok else
              "not minimal: LD+INC alone is inconsistent (Thm A6); H1 only shortens the certificate")
    print("  witness a(x)b = 2a satisfies RD, H1, H2 (z3):", verify_2a())
    print("\nPart 2: all 2^7 = 128 subsets of", AX)
    st = classify()
    for k, v in st.items():
        print(f"  {k}: {len(v)}")
    print("  maximal consistent subsets found:")
    cons = [set(s) for s in st["consistent"]]
    maxi = [s for s in cons if not any(s < t for t in cons)]
    for s in maxi:
        print("    ", sorted(s, key=AX.index))
    print("  open subsets:", st["open"])
    st2 = classify(use_paper=True)
    print("  adding the paper cores LD+INC, RD+INC (Thm A6):",
          {k: len(v) for k, v in st2.items()})
    print("  every previously open subset contains LD+INC or RD+INC:",
          all(({"LD", "INC"} <= set(o)) or ({"RD", "INC"} <= set(o)) for o in st["open"]))
    print("\nPart 3: bi-additive a(x)b = k*a*b on Z/n with H1 and H2, n <= 12  -> (n, k, 2 != 0 in Z/n)")
    print("  ", zn_models())
    for n in (2, 3):
        res = brute_biadditive(n)
        print(f"   brute force Z/{n}: {n**(n*n)} tables, {len(res)} bi-additive, k-values {sorted(res)}"
              f" (None would mean a non-k*a*b solution)")
