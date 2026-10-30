"""TH-001/TH-002/TH-003: candidate operations (x) consistent with Howard's stated data.

For each candidate a(x)b we decide, over R (and over Z / Z_{>0} where it differs):
  H1: 1(x)1 = 2            (Rolling Stone 2015; JRE #2152 39:03)
  H2: 2(x)2 = 4            (inferred from "the square root of four is two", RS 2015)
  SQ: x(x)x = 2 has root x = 1 ("what's the square root of two? Should be one", RS 2015)
  LD, RD: left/right distributivity over ordinary +
  AS, CM: associativity, commutativity
  ID: a two-sided identity exists (and which)
  INC_Z: a(x)b > max(a,b) for all integers a,b >= 1   ("multiplication ... means to make more")
  INC_R: a(x)b > max(a,b) for all reals a,b > 0
  BAL_Z: a(x)b >= a+b for integers a,b >= 1            ("what happened to the other one")
Universal statements are decided with z3 (nonlinear real/integer arithmetic): "yes" means the
negation is UNSAT (a proof); "no" comes with a concrete counterexample from z3's model.
"""
import z3
import sympy as sp

a, b, c, x = sp.symbols('a b c x', real=True)

def zmin(p, q):
    return z3.If(p <= q, p, q)

# name -> (sympy expression or None, z3 lambda, motivation)
CANDS = [
    ("ordinary a*b",      lambda p, q: p * q,            "baseline (standard multiplication)"),
    ("a+b",               lambda p, q: p + q,            "'balance': both inputs accounted for; fits RS sqrt data"),
    ("2ab",               lambda p, q: 2 * p * q,        "forced by LD+RD+H1 (Thm A3); 'every unit is a pair'"),
    ("ab+1",              lambda p, q: p * q + 1,        "'one times itself must have an effect' (+1)"),
    ("ab+a = a(b+1)",     lambda p, q: p * q + p,        "'b times more than a' (a plus b further copies)"),
    ("ab+min(a,b)",       lambda p, q: p * q + zmin(p, q), "coordinator-suggested min-type correction"),
    ("2b (left-only)",    lambda p, q: 2 * q + 0 * p,    "LD witness keeping H1,H2"),
    ("patch: ab, but 1(x)1:=2",
     lambda p, q: z3.If(z3.And(p == 1, q == 1), p + q, p * q),
     "minimal edit of the times table"),
]

def decide(formula_vars, body, sort):
    """Return ('yes', None) if forall vars. body is valid, else ('no', model-dict)."""
    s = z3.Solver()
    s.set("timeout", 20000)
    s.add(z3.Not(body))
    r = s.check()
    if r == z3.unsat:
        return "yes", None
    if r == z3.sat:
        m = s.model()
        return "no", {str(v): m.eval(v, model_completion=True) for v in formula_vars}
    return "unknown", None

def run(sort):
    mk = z3.Reals if sort == "R" else z3.Ints
    A, B, C, E = mk('a b c e')
    one = z3.RealVal(1) if sort == "R" else z3.IntVal(1)
    rows = []
    for name, f, why in CANDS:
        res = {}
        res["H1"] = decide([], f(one, one) == 2, sort)
        res["H2"] = decide([], f(2 * one, 2 * one) == 4, sort)
        res["LD"] = decide([A, B, C], f(A, B + C) == f(A, B) + f(A, C), sort)
        res["RD"] = decide([A, B, C], f(A + B, C) == f(A, C) + f(B, C), sort)
        res["AS"] = decide([A, B, C], f(f(A, B), C) == f(A, f(B, C)), sort)
        res["CM"] = decide([A, B], f(A, B) == f(B, A), sort)
        # identity: exists e forall b: f(e,b)=b and f(b,e)=b   (z3 quantifier query)
        s = z3.Solver(); s.set("timeout", 20000)
        s.add(z3.ForAll([B], z3.And(f(E, B) == B, f(B, E) == B)))
        r = s.check()
        res["ID"] = ("yes", {"e": s.model().eval(E)}) if r == z3.sat else (("no", None) if r == z3.unsat else ("unknown", None))
        if sort == "Z":
            pos = z3.And(A >= 1, B >= 1)
            res["INC"] = decide([A, B], z3.Implies(pos, z3.And(f(A, B) > A, f(A, B) > B)), sort)
            res["BAL"] = decide([A, B], z3.Implies(pos, f(A, B) >= A + B), sort)
        else:
            pos = z3.And(A > 0, B > 0)
            res["INC"] = decide([A, B], z3.Implies(pos, z3.And(f(A, B) > A, f(A, B) > B)), sort)
            res["BAL"] = decide([A, B], z3.Implies(pos, f(A, B) >= A + B), sort)
        rows.append((name, why, res))
    return rows

def show(rows, sort):
    keys = ["H1", "H2", "LD", "RD", "AS", "CM", "ID", "INC", "BAL"]
    print(f"\n=== carrier {'R' if sort=='R' else 'Z'} (INC/BAL over {'reals > 0' if sort=='R' else 'integers >= 1'}) ===")
    print("candidate".ljust(26) + "".join(k.ljust(5) for k in keys))
    for name, why, res in rows:
        print(name[:25].ljust(26) + "".join(res[k][0][0].upper().ljust(5) if res[k][0] != "unknown" else "?".ljust(5) for k in keys))
    print("counterexamples / witnesses:")
    for name, why, res in rows:
        det = "; ".join(f"{k}:{res[k][1]}" for k in keys if res[k][1])
        print(f"  {name}: {det}")

def sqrt_table():
    """Solve x(x)x = y for y in {2,4}, x > 0, for the polynomial candidates (sympy, exact)."""
    print("\n=== (x)-square roots: positive real solutions of x(x)x = y (exact) ===")
    polys = [("ordinary a*b", x * x), ("a+b", 2 * x), ("2ab", 2 * x**2), ("ab+1", x**2 + 1),
             ("ab+a", x**2 + x), ("ab+min(a,b)", x**2 + x), ("2b (left-only)", 2 * x)]
    for name, expr in polys:
        r2 = [s for s in sp.solve(sp.Eq(expr, 2), x) if s.is_positive]
        r4 = [s for s in sp.solve(sp.Eq(expr, 4), x) if s.is_positive]
        print(f"  {name:16s} sqrt(2) -> {r2}   sqrt(4) -> {r4}")
    print("  Howard (RS 2015) requires sqrt(4) = 2 and sqrt(2) = 1.")

if __name__ == "__main__":
    print("z3", z3.get_version_string(), "sympy", sp.__version__)
    for sort in ("R", "Z"):
        show(run(sort), sort)
    sqrt_table()
