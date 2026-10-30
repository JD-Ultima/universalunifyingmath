"""Python transcription of ChainAlpha.lean, Parts 0-1: the definitions on which the readout
`readout objectWorld objectTwoPoint` depends (audit F-CT-4, next step 4). Used by mechanical_mutations.py.

Each Lean definition is transcribed line by line with Lean's semantics: `getD` defaults for out-of-range indices,
truncated subtraction, x / 0 = 0 and x % 0 = x on Nat, Euclidean division and remainder on Int (Int.ediv, Int.emod;
checked with Lean 4.35.0: (-7)/2 = -4, 7/(-2) = -3, (-7)/(-2) = 4), truncating zip, Option as None. All functions
are pure, so per-evaluation memoisation does not change any value.

Every mutation site found by mechanical_mutations.py has one hook here, named '<definition>:<operator>:<k>' (the k-th
site of that operator in that definition, in source order); `H` holds the hooks of the mutant being evaluated.
Swapping the branches of an if is evaluated as negating its condition (the same function). Hooks named 'A:<name>'
are the author's rewordings from mutations.py that are not single or paired mutants; they are used only to reproduce
upstream_mutations_lean.out.

The only departure from Lean is a resource bound: a readout list longer than MAX_TERMS = 4000 terms raises
Infeasible. One single mutant duplicates the list at every record (about 2^172 terms; Lean aborts with "deep
recursion"); the longest list of any feasible mutant is printed by mechanical_mutations.py.
"""

# Literal sites of type Int (all other literal sites are Nat); a Nat literal cannot become negative.
INT_LITS = {"pairKet:LIT:0", "pairKet:LIT:1", "coeffOf:LIT:0", "orderSign:LIT:0", "orderSign:LIT:1", "wY:LIT:0",
            "kOf:LIT:0", "kOf:LIT:2", "sigX:LIT:0", "sigX:LIT:1", "sigY:LIT:0", "sigY:LIT:1"} \
    | {f"tracelessBasis:LIT:{k}" for k in (0, 1, 2, 3, 4, 6, 7)} | {f"parity:LIT:{k}" for k in range(3)} \
    | {f"recs:LIT:{k}" for k in range(6)} | {f"addRoot:LIT:{k}" for k in range(3)}
# Projection sites on a pair of two different types (Int x Assign): exchanging .1 and .2 there is ill-typed.
HETERO_PROJ = {"coeffOf:VAR:0", "coeffOf:VAR:1", "commutatorCandidates:VAR:12", "commutatorCandidates:VAR:13",
               "recs:VAR:2", "recs:VAR:3", "survives:VAR:0", "survives:VAR:1", "sigY:VAR:0", "sigY:VAR:1",
               "wY:VAR:0"} | {f"parity:VAR:{k}" for k in range(4)} | {f"sameForm:VAR:{k}" for k in (2, 3, 6, 7)}

H = {}            # active hooks of the mutant being evaluated: site id -> variant
CACHE = {}        # per-evaluation memo (all functions are pure)
NAT, INT = "N", "I"

# ---------------- hooks ----------------

def L(s, v):                       # literal
    return H.get(s, v)


def V(s, a, b):                    # label replaced by its partner
    return b if s in H else a


def P(s, tup, idx):                # projection .1 / .2 (idx 0 / 1), replaced by the other component
    return tup[1 - idx] if s in H else tup[idx]


def IFC(s, c):                     # if-then-else with swapped branches  ==  negated condition
    return (not c) if s in H else c


def EQ(s, a, b):
    return (a != b) if s in H else (a == b)


def NE(s, a, b):
    return (a == b) if s in H else (a != b)


def LT(s, a, b):
    m = H.get(s)
    return a < b if m is None else (a <= b if m == "strict" else a >= b)


def LE(s, a, b):
    m = H.get(s)
    return a <= b if m is None else (a < b if m == "strict" else a > b)


def NOTB(s, x):                    # negation of an eqZ / neZ / isSome test
    return (not x) if s in H else x


def AND(s, a, b):
    m = H.get(s)
    return (a and b) if m is None else (a or b) if m == "swap" else b if m == "dropL" else a


def OR(s, a, b):
    m = H.get(s)
    return (a or b) if m is None else (a and b) if m == "swap" else b if m == "dropL" else a


def QALL(s, xs, p):
    return any(p(x) for x in xs) if s in H else all(p(x) for x in xs)


def NEGI(s, x):                    # unary minus, deleted
    return x if s in H else -x


def NOT(s, b):                     # Boolean not, deleted
    return b if s in H else (not b)


def ediv(a, b):
    if b == 0:
        return 0
    r = a % abs(b)
    return (a - r) // b


def emod(a, b):
    return a if b == 0 else a % abs(b)


def binop(op, a, b, ty):
    if op == "+":
        return a + b
    if op == "-":
        return max(a - b, 0) if ty == NAT else a - b
    if op == "*":
        return a * b
    if op == "/":
        return (a // b if b else 0) if ty == NAT else ediv(a, b)
    if op == "%":
        return (a % b if b else a) if ty == NAT else emod(a, b)
    raise ValueError(op)


def AR(s, op, a, b, ty):
    return binop(H.get(s, op), a, b, ty)


def MN(s, a, b):
    return max(a, b) if s in H else min(a, b)


def _seq(s, items):
    m = H.get(s)
    if m is None:
        return list(items)
    tag, x = m
    items = list(items)
    return items[:x] + items[x + 1:] if tag == "drop" else items[:x + 1] + [items[x]] + items[x + 1:]


def LST(s, entries):               # list literal: drop / duplicate an entry
    return _seq(s, entries)


def CAT(s, parts):                 # ++ chain: drop / duplicate a summand
    out = []
    for p in _seq(s, parts):
        out += p
    return out


def CONS(s, elems, tail):          # e1 :: e2 :: tail: drop / duplicate an element
    return _seq(s, elems) + list(tail)


def memo(f):
    name = f.__name__

    def g(*a):
        k = (name,) + a
        try:
            return CACHE[k]
        except KeyError:
            v = CACHE[k] = f(*a)
            return v
    g.__name__ = name
    return g


# ---------------- Part 0 (held fixed) ----------------

def rankAux(fuel, rows):
    if fuel == 0:
        return 0
    rows = [r for r in rows if any(x != 0 for x in r)]
    if not rows:
        return 0
    r, rest = rows[0], rows[1:]
    c = next((i for i, x in enumerate(r) if x != 0), 0)
    a = r[c] if c < len(r) else 0
    return 1 + rankAux(fuel - 1, [[a * x - (s[c] if c < len(s) else 0) * y for x, y in zip(s, r)] for s in rest])


def rank(rows):
    return rankAux(len(rows), rows)


def mget(z, m, i, j):
    row = m[i] if i < len(m) else []
    return row[j] if j < len(row) else z


def tabulate(n, k, f):
    return [[f(i, j) for j in range(k)] for i in range(n)]


def mulZ(a, b):
    k = len(b[0]) if b else 0
    return [[sum(x * mget(0, b, i, j) for i, x in enumerate(row)) for j in range(k)] for row in a]


def traceZ(a):
    return sum(mget(0, a, i, i) for i in range(len(a)))


def transposeZ(a):
    k = len(a[0]) if a else 0
    return tabulate(k, len(a), lambda i, j: mget(0, a, j, i))


def subZ(a, b):
    return [[x - y for x, y in zip(r, s)] for r, s in zip(a, b)]


def pairZ(a, b):
    return traceZ(mulZ(transposeZ(a), b))


# ---------------- Part 1 ----------------

class World:
    __slots__ = ("qudits", "dim", "links")

    def __init__(self, qudits, dim, links):
        self.qudits, self.dim, self.links = qudits, dim, links


class TwoPoint:
    __slots__ = ("center", "read1", "read2", "state")

    def __init__(self, center, read1, read2, state):
        self.center, self.read1, self.read2, self.state = center, read1, read2, state


def pairKet(dA, dB, i, j):
    c = AND("pairKet:LOG:0",
            EQ("pairKet:REL:0", V("pairKet:VAR:0", i, j), V("pairKet:VAR:1", j, i)),
            LT("pairKet:REL:1", V("pairKet:VAR:2", i, j),
               MN("pairKet:AOR:0", V("pairKet:VAR:3", dA, dB), V("pairKet:VAR:4", dB, dA))))
    return L("pairKet:LIT:0", 1) if IFC("pairKet:IF:0", c) else L("pairKet:LIT:1", 0)


def steered(dA, dB, i):
    return [pairKet(V("steered:VAR:1", dA, dB), V("steered:VAR:2", dB, dA), i, j)
            for j in range(V("steered:VAR:0", dB, dA))]


@memo
def fixesFar(dA, dB):
    if "A:fixesFar_always_true" in H:
        return True
    rows = [steered(V("fixesFar:VAR:1", dA, dB), V("fixesFar:VAR:2", dB, dA), i)
            for i in range(V("fixesFar:VAR:0", dA, dB))]
    return EQ("fixesFar:REL:0", rank(rows), V("fixesFar:VAR:3", dB, dA))


def readEnd(dA, dB):
    c = fixesFar(V("readEnd:VAR:0", dA, dB), V("readEnd:VAR:1", dB, dA))
    return V("readEnd:VAR:2", dA, dB) if IFC("readEnd:IF:0", c) else V("readEnd:VAR:3", dB, dA)


def farEnd(dA, dB):
    c = fixesFar(V("farEnd:VAR:0", dA, dB), V("farEnd:VAR:1", dB, dA))
    return V("farEnd:VAR:2", dB, dA) if IFC("farEnd:IF:0", c) else V("farEnd:VAR:3", dA, dB)


def readoutDim(w, t):
    acc = L("readoutDim:LIT:2", 1)
    for l in w.links:
        d1 = w.dim(P("readoutDim:VAR:0", l, 0))
        d2 = w.dim(P("readoutDim:VAR:1", l, 1))
        c = OR("readoutDim:LOG:0",
               EQ("readoutDim:REL:0", l, V("readoutDim:VAR:2", t.read1, t.read2)),
               EQ("readoutDim:REL:1", l, V("readoutDim:VAR:3", t.read2, t.read1)))
        if IFC("readoutDim:IF:0", c):
            x = AR("readoutDim:AOR:1", "*",
                   AR("readoutDim:AOR:0", "*", V("readoutDim:VAR:4", d1, d2), V("readoutDim:VAR:5", d2, d1), NAT),
                   AR("readoutDim:AOR:2", "*", V("readoutDim:VAR:6", d1, d2), V("readoutDim:VAR:7", d2, d1), NAT), NAT)
        elif "A:register_full_unread" in H:
            x = (d1 * d2) * (d1 * d2)
        else:
            x = AR("readoutDim:AOR:5", "*",
                   AR("readoutDim:AOR:4", "-",
                      AR("readoutDim:AOR:3", "*", V("readoutDim:VAR:8", d1, d2), V("readoutDim:VAR:9", d1, d2), NAT),
                      L("readoutDim:LIT:0", 1), NAT),
                   AR("readoutDim:AOR:7", "-",
                      AR("readoutDim:AOR:6", "*", V("readoutDim:VAR:10", d2, d1), V("readoutDim:VAR:11", d2, d1), NAT),
                      L("readoutDim:LIT:1", 1), NAT), NAT)
        acc = AR("readoutDim:AOR:8", "*", acc, x, NAT)
    return acc


@memo
def tracelessBasis(d):
    pairs = []
    for i in range(d):
        js = [j for j in range(d) if LT("tracelessBasis:REL:0", V("tracelessBasis:VAR:0", i, j),
                                        V("tracelessBasis:VAR:1", j, i))]
        pairs += [(V("tracelessBasis:VAR:2", i, j), V("tracelessBasis:VAR:3", j, i)) for j in js]

    def symf(p):
        def f(a, b):
            c = OR("tracelessBasis:LOG:1",
                   AND("tracelessBasis:LOG:0",
                       EQ("tracelessBasis:REL:1", V("tracelessBasis:VAR:4", a, b), P("tracelessBasis:VAR:5", p, 0)),
                       EQ("tracelessBasis:REL:2", V("tracelessBasis:VAR:6", b, a), P("tracelessBasis:VAR:7", p, 1))),
                   AND("tracelessBasis:LOG:2",
                       EQ("tracelessBasis:REL:3", V("tracelessBasis:VAR:8", a, b), P("tracelessBasis:VAR:9", p, 1)),
                       EQ("tracelessBasis:REL:4", V("tracelessBasis:VAR:10", b, a), P("tracelessBasis:VAR:11", p, 0))))
            return L("tracelessBasis:LIT:0", 1) if IFC("tracelessBasis:IF:0", c) else L("tracelessBasis:LIT:1", 0)
        return tabulate(d, d, f)

    def antif(p):
        def f(a, b):
            c1 = AND("tracelessBasis:LOG:3",
                     EQ("tracelessBasis:REL:5", V("tracelessBasis:VAR:12", a, b), P("tracelessBasis:VAR:13", p, 0)),
                     EQ("tracelessBasis:REL:6", V("tracelessBasis:VAR:14", b, a), P("tracelessBasis:VAR:15", p, 1)))
            if IFC("tracelessBasis:IF:1", c1):
                return L("tracelessBasis:LIT:2", 1)
            c2 = AND("tracelessBasis:LOG:4",
                     EQ("tracelessBasis:REL:7", V("tracelessBasis:VAR:16", a, b), P("tracelessBasis:VAR:17", p, 1)),
                     EQ("tracelessBasis:REL:8", V("tracelessBasis:VAR:18", b, a), P("tracelessBasis:VAR:19", p, 0)))
            return L("tracelessBasis:LIT:3", -1) if IFC("tracelessBasis:IF:2", c2) else L("tracelessBasis:LIT:4", 0)
        return tabulate(d, d, f)

    def diagf(k):
        def f(a, b):
            c1 = AND("tracelessBasis:LOG:5",
                     EQ("tracelessBasis:REL:10", V("tracelessBasis:VAR:20", a, b), V("tracelessBasis:VAR:21", b, a)),
                     LT("tracelessBasis:REL:11", V("tracelessBasis:VAR:22", a, b), k))
            if IFC("tracelessBasis:IF:3", c1):
                return L("tracelessBasis:LIT:6", 1)
            c2 = AND("tracelessBasis:LOG:6",
                     EQ("tracelessBasis:REL:12", V("tracelessBasis:VAR:23", a, b), V("tracelessBasis:VAR:24", b, a)),
                     EQ("tracelessBasis:REL:13", V("tracelessBasis:VAR:25", a, b), k))
            return NEGI("tracelessBasis:NEG:0", k) if IFC("tracelessBasis:IF:4", c2) else L("tracelessBasis:LIT:7", 0)
        return tabulate(d, d, f)

    sym = [symf(p) for p in pairs]
    anti = [antif(p) for p in pairs]
    ks = [k for k in range(d) if LE("tracelessBasis:REL:9", L("tracelessBasis:LIT:5", 1), k)]
    diag = [diagf(k) for k in ks]
    return CAT("tracelessBasis:CAT:0", [sym, anti, diag])


@memo
def slotEnd(s):
    return AR("slotEnd:AOR:0", "%", s, L("slotEnd:LIT:0", 2), NAT)


@memo
def slotDim(dA, dB, s):
    c = EQ("slotDim:REL:0", slotEnd(s), L("slotDim:LIT:0", 0))
    return V("slotDim:VAR:0", dA, dB) if IFC("slotDim:IF:0", c) else V("slotDim:VAR:1", dB, dA)


@memo
def slotBasis(dA, dB, s):
    return tracelessBasis(slotDim(V("slotBasis:VAR:0", dA, dB), V("slotBasis:VAR:1", dB, dA), s))


@memo
def exch(s):
    return AR("exch:AOR:1", "%", AR("exch:AOR:0", "+", s, L("exch:LIT:0", 2), NAT), L("exch:LIT:1", 4), NAT)


def letterAt(a, s):
    return a[s] if s < len(a) else None


@memo
def letterSlots(a):
    return [s for s in range(L("letterSlots:LIT:0", 4)) if NOTB("letterSlots:REL:0", letterAt(a, s) is not None)]


@memo
def exchAssign(a):
    return tuple(letterAt(a, exch(s)) for s in range(L("exchAssign:LIT:0", 4)))


@memo
def letterCode(x):
    if x is None:
        return L("letterCode:LIT:0", 0)
    return AR("letterCode:AOR:0", "+", x, L("letterCode:LIT:1", 1), NAT)


@memo
def key(a):
    acc = L("key:LIT:2", 0)
    for s in range(L("key:LIT:0", 4)):
        acc = AR("key:AOR:1", "+", AR("key:AOR:0", "*", acc, L("key:LIT:1", 100), NAT), letterCode(letterAt(a, s)), NAT)
    return acc


@memo
def products(dA, dB):
    out = []
    for s1 in range(L("products:LIT:0", 4)):
        s2s = [s2 for s2 in range(L("products:LIT:1", 4))
               if LT("products:REL:0", V("products:VAR:0", s1, s2), V("products:VAR:1", s2, s1))]
        for s2 in s2s:
            n1 = len(slotBasis(V("products:VAR:2", dA, dB), V("products:VAR:3", dB, dA), V("products:VAR:4", s1, s2)))
            n2 = len(slotBasis(V("products:VAR:5", dA, dB), V("products:VAR:6", dB, dA), V("products:VAR:7", s2, s1)))
            ts = range(L("products:LIT:2", 4))
            for i1 in range(n1):
                for i2 in range(n2):
                    asg = []
                    for t in ts:
                        if IFC("products:IF:0", EQ("products:REL:1", t, V("products:VAR:8", s1, s2))):
                            asg.append(V("products:VAR:9", i1, i2))
                        elif IFC("products:IF:1", EQ("products:REL:2", t, V("products:VAR:10", s2, s1))):
                            asg.append(V("products:VAR:11", i2, i1))
                        else:
                            asg.append(None)
                    out.append(tuple(asg))
    return out


def coeffOf(f, a):
    acc = L("coeffOf:LIT:0", 0)
    for t in f:
        if EQ("coeffOf:REL:0", t[1], a):
            acc = AR("coeffOf:AOR:0", "+", acc, t[0], INT)
    return acc


def parity(f):
    c0 = QALL("parity:QNT:0", f, lambda t: NOTB("parity:REL:0", coeffOf(f, exchAssign(t[1])) == t[0]))
    if IFC("parity:IF:0", c0):
        return L("parity:LIT:0", 1)
    c1 = QALL("parity:QNT:1", f, lambda t: NOTB("parity:REL:1", coeffOf(f, exchAssign(t[1])) == NEGI("parity:NEG:0", t[0])))
    return L("parity:LIT:1", -1) if IFC("parity:IF:1", c1) else L("parity:LIT:2", 0)


@memo
def weylPhaseExp(d, j):
    A = AR
    x = A("weylPhaseExp:AOR:1", "+", A("weylPhaseExp:AOR:0", "*", L("weylPhaseExp:LIT:0", 2), j, NAT),
          A("weylPhaseExp:AOR:2", "*", L("weylPhaseExp:LIT:1", 2), d, NAT), NAT)
    y = A("weylPhaseExp:AOR:4", "*", L("weylPhaseExp:LIT:2", 2),
          A("weylPhaseExp:AOR:6", "%", A("weylPhaseExp:AOR:5", "+", j, L("weylPhaseExp:LIT:3", 1), NAT), d, NAT), NAT)
    return A("weylPhaseExp:AOR:7", "%", A("weylPhaseExp:AOR:3", "-", x, y, NAT),
             A("weylPhaseExp:AOR:8", "*", L("weylPhaseExp:LIT:4", 2), d, NAT), NAT)


@memo
def weylIsSign(d):
    return QALL("weylIsSign:QNT:0", range(d), lambda j: EQ("weylIsSign:REL:0", weylPhaseExp(d, j), d))


def coeffOn(basis, m, k):
    bk = basis[k] if k < len(basis) else []
    return AR("coeffOn:AOR:0", "/", pairZ(bk, m), pairZ(bk, bk), INT)


@memo
def commAt(dA, dB, i, j, s):
    b = slotBasis(V("commAt:VAR:0", dA, dB), V("commAt:VAR:1", dB, dA), s)

    def g(x):
        return b[x] if x < len(b) else []
    c = subZ(mulZ(g(V("commAt:VAR:2", i, j)), g(V("commAt:VAR:3", j, i))),
             mulZ(g(V("commAt:VAR:4", j, i)), g(V("commAt:VAR:5", i, j))))
    out = []
    for k in range(len(b)):
        if NOTB("commAt:REL:0", coeffOn(b, c, k) != 0):
            asg = tuple((k if IFC("commAt:IF:0", EQ("commAt:REL:1", t, s)) else None)
                        for t in range(L("commAt:LIT:0", 4)))
            out.append((coeffOn(b, c, k), asg))
    return out


@memo
def orderSign(d):
    return L("orderSign:LIT:0", -1) if IFC("orderSign:IF:0", weylIsSign(d)) else L("orderSign:LIT:1", 1)


@memo
def commutatorCandidates(dA, dB, i, j):
    c1 = commAt(V("commutatorCandidates:VAR:0", dA, dB), V("commutatorCandidates:VAR:1", dB, dA),
                V("commutatorCandidates:VAR:2", i, j), V("commutatorCandidates:VAR:3", j, i),
                L("commutatorCandidates:LIT:0", 1))
    c3 = commAt(V("commutatorCandidates:VAR:4", dA, dB), V("commutatorCandidates:VAR:5", dB, dA),
                V("commutatorCandidates:VAR:6", i, j), V("commutatorCandidates:VAR:7", j, i),
                L("commutatorCandidates:LIT:1", 3))
    e0 = CAT("commutatorCandidates:CAT:0", [V("commutatorCandidates:VAR:8", c1, c3),
                                            V("commutatorCandidates:VAR:9", c3, c1)])
    e1 = CAT("commutatorCandidates:CAT:1", [V("commutatorCandidates:VAR:10", c1, c3),
                                            [(NEGI("commutatorCandidates:NEG:0", t[0]), t[1])
                                             for t in V("commutatorCandidates:VAR:11", c3, c1)]])
    return LST("commutatorCandidates:LIST:0", [e0, e1])


def sameEnd(f):
    if f:
        ls = letterSlots(f[0][1])
        if len(ls) == 2:
            s1, s2 = ls
            return EQ("sameEnd:REL:0", slotEnd(V("sameEnd:VAR:0", s1, s2)), slotEnd(V("sameEnd:VAR:1", s2, s1)))
        return L("sameEnd:BLIT:0", False)
    return L("sameEnd:BLIT:1", False)


@memo
def commutatorRecs(dA, dB):
    cond = weylIsSign(V("commutatorRecs:VAR:0", dB, dA))
    if "A:no_commutator_records" in H:
        cond = False
    if not IFC("commutatorRecs:IF:0", cond):
        return []
    nb = len(slotBasis(V("commutatorRecs:VAR:1", dA, dB), V("commutatorRecs:VAR:2", dB, dA), L("commutatorRecs:LIT:0", 1)))
    out = []
    for i in range(nb):
        js = [j for j in range(nb) if LT("commutatorRecs:REL:0", V("commutatorRecs:VAR:3", i, j), V("commutatorRecs:VAR:4", j, i))]
        for j in js:
            A = commutatorCandidates(V("commutatorRecs:VAR:5", dA, dB), V("commutatorRecs:VAR:6", dB, dA),
                                     V("commutatorRecs:VAR:7", i, j), V("commutatorRecs:VAR:8", j, i))
            B = commutatorCandidates(V("commutatorRecs:VAR:9", dA, dB), V("commutatorRecs:VAR:10", dB, dA),
                                     V("commutatorRecs:VAR:11", j, i), V("commutatorRecs:VAR:12", i, j))
            for fg in zip(A, B):
                os_ = orderSign(V("commutatorRecs:VAR:14", dB, dA))
                if "A:commutators_keep_sums" in H:
                    os_ = -os_
                c = NOTB("commutatorRecs:REL:1", parity(P("commutatorRecs:VAR:13", fg, 0)) == os_)
                if "A:commutators_keep_both" in H:
                    c = True
                if IFC("commutatorRecs:IF:1", c):
                    reads = None if "A:commutator_records_carry_no_reads" in H else P("commutatorRecs:VAR:17", fg, 1)
                    out.append((P("commutatorRecs:VAR:15", fg, 0), orderSign(V("commutatorRecs:VAR:16", dB, dA)), reads))
    return out


@memo
def recs(dA, dB):
    ps = products(V("recs:VAR:0", dA, dB), V("recs:VAR:1", dB, dA))

    def mk(f):
        reads = [(u[0], exchAssign(u[1])) for u in f] if IFC("recs:IF:0", sameEnd(f)) else None
        return (f, L("recs:LIT:0", 1), reads)
    fixedP = [mk(LST("recs:LIST:0", [(L("recs:LIT:1", 1), p)])) for p in ps if EQ("recs:REL:0", exchAssign(p), p)]
    reps = [p for p in ps if LT("recs:REL:1", key(p), key(exchAssign(p)))]
    odd = [mk(LST("recs:LIST:1", [(L("recs:LIT:2", 1), p), (L("recs:LIT:3", -1), exchAssign(p))])) for p in reps]
    even = [mk(LST("recs:LIST:2", [(L("recs:LIT:4", 1), p), (L("recs:LIT:5", 1), exchAssign(p))])) for p in reps]
    return CAT("recs:CAT:0", [V("recs:VAR:4", odd, even), V("recs:VAR:5", even, odd), fixedP,
                              commutatorRecs(V("recs:VAR:6", dA, dB), V("recs:VAR:7", dB, dA))])


@memo
def slotQudit(t, s):
    if "A:copies_separate" in H:
        if slotEnd(s) == 0:
            return 100 + s // 2
    elif IFC("slotQudit:IF:0", EQ("slotQudit:REL:0", slotEnd(s), L("slotQudit:LIT:0", 0))):
        return t.center
    if IFC("slotQudit:IF:1", EQ("slotQudit:REL:1", AR("slotQudit:AOR:0", "/", s, L("slotQudit:LIT:1", 2), NAT),
                                 L("slotQudit:LIT:2", 0))):
        return P("slotQudit:VAR:1", V("slotQudit:VAR:0", t.read1, t.read2), 1)
    return P("slotQudit:VAR:3", V("slotQudit:VAR:2", t.read2, t.read1), 1)


@memo
def linkSwap(t, q):
    if "A:copies_separate" in H:
        if q == 100:
            return 101
        if q == 101:
            return 100
    r1, r2 = t.read1, t.read2
    if IFC("linkSwap:IF:0", EQ("linkSwap:REL:0", q, P("linkSwap:VAR:1", V("linkSwap:VAR:0", r1, r2), 1))):
        return P("linkSwap:VAR:3", V("linkSwap:VAR:2", r2, r1), 1)
    if IFC("linkSwap:IF:1", EQ("linkSwap:REL:1", q, P("linkSwap:VAR:5", V("linkSwap:VAR:4", r2, r1), 1))):
        return P("linkSwap:VAR:7", V("linkSwap:VAR:6", r1, r2), 1)
    if IFC("linkSwap:IF:2", EQ("linkSwap:REL:2", q, P("linkSwap:VAR:9", V("linkSwap:VAR:8", r1, r2), 0))):
        return P("linkSwap:VAR:11", V("linkSwap:VAR:10", r2, r1), 0)
    if IFC("linkSwap:IF:3", EQ("linkSwap:REL:3", q, P("linkSwap:VAR:13", V("linkSwap:VAR:12", r2, r1), 0))):
        return P("linkSwap:VAR:15", V("linkSwap:VAR:14", r1, r2), 0)
    return q


def occupies(t, a):
    return [slotQudit(t, s) for s in letterSlots(a)]


def sameSet(xs, ys):
    a = QALL("sameSet:QNT:0", V("sameSet:VAR:0", xs, ys), lambda x: x in V("sameSet:VAR:1", ys, xs))
    b = QALL("sameSet:QNT:1", V("sameSet:VAR:2", ys, xs), lambda y: y in V("sameSet:VAR:3", xs, ys))
    return AND("sameSet:LOG:0", a, b)


def moved(t, f):
    if f:
        a = f[0][1]
        if "A:moved_per_letter" in H:
            return any(linkSwap(t, q) != q for q in occupies(t, a))
        o = occupies(t, a)
        return NOT("moved:NEG:0", sameSet([linkSwap(t, q) for q in o], o))
    return L("moved:BLIT:0", False)


def sameForm(f, g):
    F, G = V("sameForm:VAR:0", f, g), V("sameForm:VAR:1", g, f)
    a = QALL("sameForm:QNT:0", F, lambda u: NOTB("sameForm:REL:0", coeffOf(G, u[1]) == u[0]))
    G2, F2 = V("sameForm:VAR:4", g, f), V("sameForm:VAR:5", f, g)
    b = QALL("sameForm:QNT:1", G2, lambda u: NOTB("sameForm:REL:1", coeffOf(F2, u[1]) == u[0]))
    return AND("sameForm:LOG:0", a, b)


def survives(r, image):
    if image is not None:
        if "A:survival_ignores_own_sign" in H:
            return sameForm(image, r[0])
        return sameForm(image, [(AR("survives:AOR:0", "*", r[1], u[0], INT), u[1]) for u in r[0]])
    return L("survives:BLIT:0", False)


def sigX(r):
    return L("sigX:LIT:0", -1) if IFC("sigX:IF:0", survives(r, r[2])) else L("sigX:LIT:1", 1)


def sigY(t, r):
    img = [(u[0], exchAssign(u[1])) for u in r[0]] if IFC("sigY:IF:1", moved(t, r[0])) else None
    return L("sigY:LIT:0", -1) if IFC("sigY:IF:0", survives(r, img)) else L("sigY:LIT:1", 1)


def wX(dA, dB):
    if "A:wX_dA_equiv" in H:
        return dA
    return (farEnd if "wX:VAR:0" in H else readEnd)(V("wX:VAR:1", dA, dB), V("wX:VAR:2", dB, dA))


def wY(w, t, r):
    hd = r[0][0] if r[0] else (L("wY:LIT:0", 0), ())
    a = hd[1]

    def at_(e):
        ss = [s for s in letterSlots(a)
              if AND("wY:LOG:0", EQ("wY:REL:0", slotEnd(s), e),
                     NE("wY:REL:1", linkSwap(t, slotQudit(t, s)), slotQudit(t, s)))]
        if ss:
            return w.dim(t.center) if "A:Y_at_center" in H else w.dim(slotQudit(t, ss[0]))
        return L("wY:LIT:1", 1)
    return AR("wY:AOR:0", "*", at_(L("wY:LIT:2", 0)), at_(L("wY:LIT:3", 1)), NAT)


def kOf(x, w, t, r):
    sx = IFC("kOf:IF:0", NOTB("kOf:REL:0", sigX(r) == L("kOf:LIT:0", -1)))
    sy = IFC("kOf:IF:1", NOTB("kOf:REL:1", sigY(t, r) == L("kOf:LIT:2", -1)))
    if "A:X_where_letters_sit" in H:
        hd = r[0][0] if r[0] else (0, ())
        a = (2 if all(slotEnd(s) == 1 for s in letterSlots(hd[1])) else x) if sx else 1
        return a * (wY(w, t, r) if sy else 1)
    if "A:X_Y_factors_swapped" in H:
        return (2 if sx else 1) * (x if sy else 1)
    a = x if sx else L("kOf:LIT:1", 1)
    b = wY(w, t, r) if sy else L("kOf:LIT:3", 1)
    return AR("kOf:AOR:0", "*", a, b, NAT)


class Infeasible(Exception):
    """the readout list outgrows MAX_TERMS (a mutant whose fold grows exponentially)"""


MAX_TERMS = 4000


def addRoot(k, lst):
    # the structural recursion, unrolled: addRoot k (x :: rest) = seq(x) ++ addRoot k rest in the third case
    head = []
    for x in range(len(lst)):
        r, c = lst[x]
        rest = lst[x + 1:]
        if IFC("addRoot:IF:0", EQ("addRoot:REL:0", r, k)):
            out = head + CONS("addRoot:CONS:0", [(r, AR("addRoot:AOR:0", "+", c, L("addRoot:LIT:1", 1), INT))], rest)
            break
        if IFC("addRoot:IF:1", LT("addRoot:REL:1", k, r)):
            out = head + CONS("addRoot:CONS:1", [(k, L("addRoot:LIT:2", 1)), (r, c)], rest)
            break
        head += CONS("addRoot:CONS:2", [(r, c)], [])
        if len(head) > MAX_TERMS:
            raise Infeasible
    else:
        out = head + LST("addRoot:LIST:0", [(k, L("addRoot:LIT:0", 1))])
    if len(out) > MAX_TERMS:
        raise Infeasible
    return out


def readout(w, t):
    dc = w.dim(t.center)
    df = w.dim(P("readout:VAR:1", V("readout:VAR:0", t.read1, t.read2), 1))
    dA = (farEnd if "readout:VAR:2" in H else readEnd)(V("readout:VAR:3", dc, df), V("readout:VAR:4", df, dc))
    dB = (readEnd if "readout:VAR:5" in H else farEnd)(V("readout:VAR:6", dc, df), V("readout:VAR:7", df, dc))
    x = wX(V("readout:VAR:10", dc, df), V("readout:VAR:11", df, dc))
    acc = []
    for r in recs(V("readout:VAR:8", dA, dB), V("readout:VAR:9", dB, dA)):
        k = kOf(x, w, t, r)
        if "A:born_weights" in H:
            k = k * k
        acc = addRoot(k, acc)
    return acc, readoutDim(w, t)


def objectWorld():
    def dim(q):
        if "A:copies_separate" in H:
            return 3 if (q == 0 or 100 <= q) else 2
        return L("objectWorld:LIT:2", 3) if IFC("objectWorld:IF:0", EQ("objectWorld:REL:0", q, L("objectWorld:LIT:1", 0))) \
            else L("objectWorld:LIT:3", 2)
    links = LST("objectWorld:LIST:0", [(L("objectWorld:LIT:4", 0), L("objectWorld:LIT:5", 1)),
                                      (L("objectWorld:LIT:6", 0), L("objectWorld:LIT:7", 2)),
                                      (L("objectWorld:LIT:8", 0), L("objectWorld:LIT:9", 3))])
    return World(L("objectWorld:LIT:0", 4), dim, links)


def objectTwoPoint():
    return TwoPoint(L("objectTwoPoint:LIT:0", 0), (L("objectTwoPoint:LIT:1", 0), L("objectTwoPoint:LIT:2", 1)),
                    (L("objectTwoPoint:LIT:3", 0), L("objectTwoPoint:LIT:4", 2)),
                    (L("objectTwoPoint:LIT:5", 0), L("objectTwoPoint:LIT:6", 3)))


# ---- the other conjuncts of `main` (outside the readout; used only to reproduce the author's log) ----

def wellFormed(w):
    ok = all(2 <= w.dim(q) for q in range(w.qudits))
    ok = ok and all(l[0] != l[1] and l[0] < w.qudits and l[1] < w.qudits for l in w.links)
    ls = w.links
    ok = ok and all(not (ls[x] == ls[y] or (ls[x][0] == ls[y][1] and ls[x][1] == ls[y][0]))
                    for x in range(len(ls)) for y in range(x + 1, len(ls)))
    return ok


def valid(w, t):
    ls = [t.read1, t.read2, t.state]
    typeOf = lambda l: (w.dim(l[0]), w.dim(l[1]))
    return (wellFormed(w) and all(l in w.links and l[0] == t.center for l in ls) and t.read1 != t.read2
            and t.read1 != t.state and t.read2 != t.state and all(typeOf(l) == typeOf(t.read1) for l in ls)
            and len(w.links) == 3)


BASE_READOUT = ([(1, 79), (2, 48), (3, 42), (6, 3)], 31104)


def evaluate(hooks):
    """readout of objectWorld objectTwoPoint under the given hooks"""
    global H
    H = hooks
    CACHE.clear()
    try:
        return readout(objectWorld(), objectTwoPoint())
    except Infeasible:
        return "INFEASIBLE"


def main_and_comm(hooks):
    """(main proved?, commutators proved?, readout) as the author's mutations.py reports them"""
    global H
    H = hooks
    CACHE.clear()
    w, t = objectWorld(), objectTwoPoint()
    ro = readout(w, t)
    dc, df = w.dim(t.center), w.dim(t.read1[1])
    m = (valid(w, t) and fixesFar(dc, df) is True and fixesFar(df, dc) is False
         and all(parity(r[0]) != 0 for r in recs(readEnd(dc, df), farEnd(dc, df))) and ro == BASE_READOUT)
    c = all(parity(r[0]) == -1 for r in commutatorRecs(3, 2))
    return m, c, ro
