"""Mechanically generated rule family for the Chain Theory alpha readout: rule-family sensitivity (audit F-CT-4 E;
gap CT-G10).

F-CT-4 B measures sensitivity over a hand-declared rule space B (46 080 rules, freedom.py). This script adds a second
declared family, generated mechanically from the Lean definitions themselves (the operator set is also ours):

  1. It reads the pinned `ChainAlpha.lean` (sha256 below) as text, masks comments, splits Part 1 into definitions and
     computes the dependency cone of the readout `readout objectWorld objectTwoPoint` (the last conjunct of `main`):
     44 definitions. Part 0 (exact-arithmetic utilities) is held fixed.
  2. A token scanner finds every site of 15 declared operators in the 44 bodies (literals, tests, connectives,
     quantifiers, branches, match arms, arithmetic operators, end/role labels, lists and chains; OP_DESC below) and
     builds every single mutant and every pair of mutants with disjoint sites. Each mutant is a set of text edits of
     the Lean file and the matching hook settings of the Python transcription.
  3. The readout of every mutant is evaluated exactly by `chainalpha_transcription.py` (a Python transcription of
     Parts 0-1 with one hook per site). The transcription is validated here on the unmutated file and on the
     author's 19 rewordings (it reproduces upstream_mutations_lean.out line for line) and, against Lean itself, on
     every single mutant and a seeded sample of pairs by proofs/lean/chain-theory/mutation_spotcheck.sh.
  4. It prints the F-CT-4 metrics for the mechanical family (distinct exact values, density near 1/alpha, closest
     values) and compares them with B. As in freedom.py, an illustrative random-target model is reported separately:
     the target U is uniform on (120,155) while the rule outputs remain fixed, and the expected number of distinct
     outputs within 1 sigma of U is the total length of their +-sigma intervals, clipped to the band, divided by 35.
     This reference count is not a calibrated p-value for the observed alpha, and neither declared family identifies
     the construction's search history. The density depends on the mutation depth (here at most 2).

Values are compared exactly: sqrt k is reduced to s sqrt m with m squarefree, and distinct squarefree radicals are
linearly independent over Q. 1/value is evaluated with mpmath (40 digits).

The author's sources are only read (never executed). Run, from this directory:
  UPSTREAM=/path/to/chain-theory python3 mechanical_mutations.py > mechanical_mutations.out
  (or --lean /path/to/ChainAlpha.lean).  2 worker processes; see README.md for the runtime.
Modes for the Lean check: --emit-lean DIR (write the Lean files and the Python expectations), --compare DIR.
"""
import hashlib
import itertools
import re
import sys
from collections import Counter, OrderedDict, defaultdict
PIN = "522844e911cf5f6091b603db0d1d042f50544d2ec8d10ed9336b722945c331da"

# ----------------------------------------------------------------------------------------------
# lexing
# ----------------------------------------------------------------------------------------------

def mask_comments(src):
    """Replace comments (--, nested /- -/ incl. docstrings) by spaces, keeping newlines and offsets."""
    out = list(src)
    i, n, depth = 0, len(src), 0
    while i < n:
        if depth == 0 and src.startswith("--", i):
            j = src.find("\n", i)
            j = n if j < 0 else j
            for k in range(i, j):
                out[k] = " "
            i = j
        elif src.startswith("/-", i):
            depth += 1
            out[i] = out[i + 1] = " "
            i += 2
        elif depth > 0 and src.startswith("-/", i):
            depth -= 1
            out[i] = out[i + 1] = " "
            i += 2
        else:
            if depth > 0 and src[i] != "\n":
                out[i] = " "
            i += 1
    return "".join(out)


TOK = re.compile(r"(?P<nl>\n)|(?P<ws>[ \t\r]+)|(?P<num>\d+)|(?P<ident>[A-Za-z_][A-Za-z0-9_']*\??)"
                 r"|(?P<op>:=|=>|==|!=|&&|\|\||\+\+|::|<;>|<\||\|>|->|<=|>=|.)")
KEYWORDS = {"if", "then", "else", "fun", "let", "match", "with", "do", "where", "def", "in"}
OPEN, CLOSE = "([{⟨", ")]}⟩"


class Tok:
    __slots__ = ("kind", "text", "s", "e", "i", "col", "line_indent")

    def __init__(self, kind, text, s, e):
        self.kind, self.text, self.s, self.e = kind, text, s, e

    def __repr__(self):
        return f"{self.text!r}@{self.s}"


def tokenize(ms, s0, e0):
    """Tokens of masked source ms[s0:e0] (offsets are absolute). Newlines are kept as tokens."""
    toks = []
    for m in TOK.finditer(ms, s0, e0):
        if m.lastgroup == "ws":
            continue
        toks.append(Tok(m.lastgroup, m.group(), m.start(), m.end()))
    for i, t in enumerate(toks):
        t.i = i
        ls = ms.rfind("\n", 0, t.s) + 1
        t.col = t.s - ls
        k = ls
        while k < len(ms) and ms[k] == " ":
            k += 1
        t.line_indent = k - ls
    return toks


def nonl(toks):
    return [t for t in toks if t.kind != "nl"]


# ----------------------------------------------------------------------------------------------
# items, bodies, cone
# ----------------------------------------------------------------------------------------------

class Body:
    def __init__(self, owner, label, toks, params):
        self.owner, self.label, self.toks, self.params = owner, label, toks, params
        self.name = owner if label is None else f"{owner}.{label}"
        for i, t in enumerate(toks):
            t.i = i


def depth_at(toks):
    """bracket depth before each token"""
    d, out = 0, []
    for t in toks:
        if t.text in CLOSE:
            d -= 1
        out.append(d)
        if t.text in OPEN:
            d += 1
    return out


def parse_items(src, ms):
    a = src.index("/-! # Part 1")
    b = src.index("/-! ## 1.9")
    starts = [0] + [m.end() for m in re.finditer(r"\n", src)]
    tops = [p for p in starts if a <= p < b and ms[p] not in " \n"]
    items = []
    ns = None
    for x, p in enumerate(tops):
        q = tops[x + 1] if x + 1 < len(tops) else b
        toks = tokenize(ms, p, q)
        t0 = toks[0].text
        if t0 == "namespace":
            ns = toks[1].text
            continue
        if t0 == "end":
            ns = None
            continue
        if t0 not in ("def", "abbrev", "structure"):
            raise ValueError(f"unexpected item {t0} at {p}")
        k = 1
        name = toks[k].text
        while toks[k + 1].text == "." and toks[k + 1].e == toks[k + 2].s:
            name += "." + toks[k + 2].text
            k += 2
        if ns and t0 == "def":
            name = ns + "." + name
        items.append((t0, name, toks, p, q))
    return items


def header_split(toks):
    """index of the body start token and the body form ('expr', 'where', 'arms')"""
    d = 0
    for k, t in enumerate(toks):
        if t.text in OPEN:
            d += 1
        elif t.text in CLOSE:
            d -= 1
        elif d == 0 and t.text == ":=":
            return k + 1, "expr"
        elif d == 0 and t.text == "where":
            return k + 1, "where"
        elif d == 0 and t.text == "|" and k > 0 and toks[k - 1].kind == "nl":
            return k, "arms"
    raise ValueError("no body")


def header_params(toks):
    ps, d, grp = [], 0, []
    for t in toks:
        if t.text in OPEN:
            d += 1
            if d == 1:
                grp = []
                continue
        elif t.text in CLOSE:
            d -= 1
            if d == 0:
                txt = [g.text for g in grp]
                if ":" in txt:
                    ps += [g.text for g in grp[:txt.index(":")] if g.kind == "ident"]
                continue
        if d >= 1:
            grp.append(t)
    return ps


def bodies_of(items):
    defs = OrderedDict()
    for kind, name, toks, p, q in items:
        if kind != "def":
            defs[name] = (kind, [], [])
            continue
        k, form = header_split(toks)
        params = header_params(toks[:k])
        if form in ("expr", "arms"):
            bl = [Body(name, None, toks[k:], params)]
        else:
            rest = toks[k:]
            fstarts = [j for j in range(len(rest) - 1)
                       if rest[j].kind == "ident" and rest[j + 1].text == ":=" and (j == 0 or rest[j - 1].kind == "nl")]
            bl = []
            for x, j in enumerate(fstarts):
                jj = fstarts[x + 1] if x + 1 < len(fstarts) else len(rest)
                bl.append(Body(name, rest[j].text, rest[j + 2:jj], params))
        defs[name] = (kind, bl, params)
    return defs


def refs(body, names):
    out = set()
    ts = body.toks
    for k, t in enumerate(ts):
        if t.kind != "ident":
            continue
        if t.text in names:
            out.add(t.text)
        if k > 0 and ts[k - 1].text == "." and ts[k - 1].e == t.s:
            for n in names:
                if n.endswith("." + t.text):
                    out.add(n)
    return out


def cone(defs, roots):
    names = set(defs)
    seen, todo = [], list(roots)
    while todo:
        n = todo.pop(0)
        if n in seen:
            continue
        seen.append(n)
        for b in defs[n][1]:
            for r in sorted(refs(b, names)):
                if r not in seen and defs[r][0] == "def":
                    todo.append(r)
    return seen


# ----------------------------------------------------------------------------------------------
# expression structure (enough of Lean 4's term grammar for these bodies)
# ----------------------------------------------------------------------------------------------

BINPREC = {"||": (30, "l"), "&&": (35, "l"), "∨": (30, "r"), "∧": (35, "r"),
           "==": (50, "n"), "!=": (50, "n"), "=": (50, "n"), "≠": (50, "n"), "<": (50, "n"),
           "≤": (50, "n"), ">": (50, "n"), "≥": (50, "n"), "++": (65, "l"), "+": (65, "l"),
           "-": (65, "l"), "*": (70, "l"), "/": (70, "l"), "%": (70, "l"), "::": (67, "r"),
           "<|": (10, "r"), "|>": (10, "l"), "$": (10, "r")}
STOP_L = {",", ":=", "=>", "then", "else", "if", "|", "with", "fun", "let", "match", "do"}
STOP_R = {",", ":=", "=>", "then", "else", "|", "with", "do"}
UNARY_CTX = {"(", ",", "then", "else", ":=", "=>", "[", "⟨", "if"}


def prev_tok(toks, k):
    j = k - 1
    while j >= 0 and toks[j].kind == "nl":
        j -= 1
    return toks[j] if j >= 0 else None


def next_tok_idx(toks, k):
    j = k + 1
    while j < len(toks) and toks[j].kind == "nl":
        j += 1
    return j if j < len(toks) else None


def operand_end(t):
    return t.kind == "num" or (t.kind == "ident" and t.text not in KEYWORDS) or t.text in CLOSE or t.text == "·"


def is_binary(toks, k):
    t = toks[k]
    if t.text not in BINPREC:
        return False
    if t.text == "-":
        p = prev_tok(toks, k)
        return p is not None and operand_end(p)
    return True


def trim(toks, a, b):
    """shrink the token index range [a, b) so it starts and ends with non-newline tokens"""
    while a < b and toks[a].kind == "nl":
        a += 1
    while b > a and toks[b - 1].kind == "nl":
        b -= 1
    return a, b


def left_operand(toks, k):
    p, a = BINPREC[toks[k].text]
    j, d = k - 1, 0
    while j >= 0:
        t = toks[j]
        if t.kind == "nl":
            if d == 0:
                break
        elif t.text in CLOSE:
            d += 1
        elif t.text in OPEN:
            if d == 0:
                break
            d -= 1
        elif d == 0:
            if t.text in STOP_L:
                break
            if is_binary(toks, j):
                q, b = BINPREC[t.text]
                if q < p or (q == p and not (a == "l" and b == "l")):
                    break
        j -= 1
    return trim(toks, j + 1, k)


def right_operand(toks, k):
    p, a = BINPREC[toks[k].text]
    j, d = k + 1, 0
    while j < len(toks):
        t = toks[j]
        if t.kind == "nl":
            if d == 0:
                break
        elif t.text in OPEN:
            d += 1
        elif t.text in CLOSE:
            if d == 0:
                break
            d -= 1
        elif d == 0:
            if t.text in STOP_R:
                break
            if is_binary(toks, j):
                q, b = BINPREC[t.text]
                if q < p or (q == p and not (a == "r" and b == "r")):
                    break
        j += 1
    return trim(toks, k + 1, j)


def parse_if(toks, k):
    """(cond range, then range, else range) as token index ranges, for the `if` at index k"""
    j, d, pend = k + 1, 0, 0
    while True:
        t = toks[j]
        if t.text in OPEN:
            d += 1
        elif t.text in CLOSE:
            d -= 1
        elif d == 0 and t.text == "if":
            pend += 1
        elif d == 0 and t.text == "else" and pend:
            pend -= 1
        elif d == 0 and t.text == "then" and pend == 0:
            break
        j += 1
    tk = j
    j, d, pend = tk + 1, 0, 0
    while True:
        t = toks[j]
        if t.text in OPEN:
            d += 1
        elif t.text in CLOSE:
            d -= 1
        elif d == 0 and t.text == "if":
            pend += 1
        elif d == 0 and t.text == "else":
            if pend == 0:
                break
            pend -= 1
        j += 1
    ek = j
    base = toks[k].line_indent
    j, d, pend = ek + 1, 0, 0
    while j < len(toks):
        t = toks[j]
        if t.kind == "nl":
            n = next_tok_idx(toks, j)
            if d == 0 and (n is None or toks[n].line_indent < base):
                break
        elif t.text in OPEN:
            d += 1
        elif t.text in CLOSE:
            if d == 0:
                break
            d -= 1
        elif d == 0:
            if t.text in (",", "|", ":="):
                break
            if t.text == "if":
                pend += 1
            elif t.text == "else":
                if pend == 0:
                    break
                pend -= 1
            elif t.text == "then" and pend == 0:
                break
        j += 1
    return trim(toks, k + 1, tk), trim(toks, tk + 1, ek), trim(toks, ek + 1, j)


def parse_arms(toks, j):
    """arms starting at token j (a `|`): list of (bar, arrow, rhs range); and the end index"""
    arms = []
    n = len(toks)
    while j is not None and j < n and toks[j].text == "|":
        bar = j
        d = 0
        j += 1
        while not (d == 0 and toks[j].text == "=>"):
            if toks[j].text in OPEN:
                d += 1
            elif toks[j].text in CLOSE:
                d -= 1
            j += 1
        arrow = j
        j += 1
        d = 0
        nxt = None
        while j < n:
            t = toks[j]
            if t.kind == "nl" and d == 0:
                m = next_tok_idx(toks, j)
                if m is None or toks[m].text == "|" or toks[m].line_indent < toks[bar].col:
                    nxt = m if (m is not None and toks[m].text == "|") else None
                    break
            elif t.text in OPEN:
                d += 1
            elif t.text in CLOSE:
                if d == 0:
                    break
                d -= 1
            elif d == 0 and t.text == "|":
                nxt = j
                break
            j += 1
        arms.append((bar, arrow, trim(toks, arrow + 1, j)))
        if nxt is None:
            break
        j = nxt
    return arms, j


def pattern_ranges(toks):
    """token index ranges of match patterns (from `|` to `=>`)"""
    out = []
    for k, t in enumerate(toks):
        if t.text == "|":
            j, d = k + 1, 0
            while j < len(toks) and not (d == 0 and toks[j].text == "=>"):
                if toks[j].text in OPEN:
                    d += 1
                elif toks[j].text in CLOSE:
                    d -= 1
                j += 1
            out.append((k, j))
    return out


def group_end(toks, k):
    """index of the closer of the innermost bracket group containing token k (or len)"""
    j, d = k + 1, 0
    while j < len(toks):
        t = toks[j]
        if t.text in OPEN:
            d += 1
        elif t.text in CLOSE:
            if d == 0:
                return j
            d -= 1
        j += 1
    return len(toks)


def binders(body):
    """binder token indices and scopes {name: [(a, b)]} (token index ranges)"""
    toks = body.toks
    btok, scopes = set(), defaultdict(list)
    for p in body.params:
        scopes[p].append((0, len(toks)))
    for k, t in enumerate(toks):
        if t.text == "fun":
            j = k + 1
            while toks[j].text != "=>":
                if toks[j].kind == "ident":
                    btok.add(j)
                    scopes[toks[j].text].append((j, group_end(toks, k)))
                j += 1
        elif t.text == "let":
            btok.add(k + 1)
            scopes[toks[k + 1].text].append((k + 1, group_end(toks, k)))
            j = k + 2
            if toks[j].text == "(":
                while toks[j].text != ":":
                    if toks[j].kind == "ident":
                        btok.add(j)
                        scopes[toks[j].text].append((j, group_end(toks, k)))
                    j += 1
    for a, b in pattern_ranges(toks):
        arms, _ = parse_arms(toks, a)
        rhs_end = arms[0][2][1] if arms else b
        for j in range(a + 1, b):
            if toks[j].kind == "ident" and toks[j].text not in ("none", "some", "true", "false", "_"):
                btok.add(j)
                scopes[toks[j].text].append((b, rhs_end))
    return btok, scopes


def in_scope(scopes, name, k):
    return any(a <= k < b for a, b in scopes.get(name, []))


# ----------------------------------------------------------------------------------------------
# sites
# ----------------------------------------------------------------------------------------------

_VP = [("dA", "dB"), ("d1", "d2"), ("dc", "df"), ("i", "j"), ("a", "b"), ("s1", "s2"), ("i1", "i2"),
       ("c1", "c3"), ("f", "g"), ("xs", "ys"), ("odd", "even")]
VARPAIRS = {}
for _a, _b in _VP:
    VARPAIRS[_a], VARPAIRS[_b] = _b, _a
FIELDPAIRS = {"read1": "read2", "read2": "read1"}
GLOBALPAIRS = {"readEnd": "farEnd", "farEnd": "readEnd"}
PROJPAIRS = {"1": "2", "2": "1"}
PAIRGROUP = {}
for _d, _g in ((VARPAIRS, None), (FIELDPAIRS, None), (GLOBALPAIRS, None), (PROJPAIRS, None)):
    for _a, _b in _d.items():
        PAIRGROUP[_a] = "/".join(sorted((_a, _b))) if _d is not PROJPAIRS else ".1/.2"
REL_VARIANTS = {"=": [("neg", "≠")], "≠": [("neg", "=")], "==": [("neg", "!=")], "!=": [("neg", "==")],
                "<": [("strict", "≤"), ("neg", "≥")], "≤": [("strict", "<"), ("neg", ">")],
                "eqZ": [("neg", "!eqZ")], "neZ": [("neg", "!neZ")], "isSome": [("neg", "isNone")],
                "isNone": [("neg", "isSome")]}
AOPS = ["+", "-", "*", "/", "%"]


class Site:
    def __init__(self, sid, op, defname, body, k, span, **info):
        self.id, self.op, self.defname, self.body, self.k, self.span = sid, op, defname, body, k, span
        self.info = info

    def ctx(self, ms, w=28):
        s, e = self.span
        a = max(0, s - w)
        return (ms[a:s] + "«" + ms[s:e] + "»" + ms[e:e + w]).replace("\n", "↵")


def find_sites(defs, names, ms):
    sites = []
    count = defaultdict(int)

    def add(op, dn, body, k, span, **info):
        sid = f"{dn}:{op}:{count[(dn, op)]}"
        count[(dn, op)] += 1
        s = Site(sid, op, dn, body, k, span, **info)
        sites.append(s)
        return s

    for dn in names:
        for body in defs[dn][1]:
            toks = body.toks
            pats = pattern_ranges(toks)
            inpat = set()
            for a, b in pats:
                inpat.update(range(a, b + 1))
            btok, scopes = binders(body)
            events = []          # (token index, order, op, payload) to keep textual order per op
            # arms body: the def-level match
            if toks and toks[0].text == "|":
                arms, _ = parse_arms(toks, 0)
                events.append((0, "MATCH", ("arms", arms)))
            for k, t in enumerate(toks):
                if k in inpat:
                    continue
                p = prev_tok(toks, k)
                adj_dot = k > 0 and toks[k - 1].text == "." and toks[k - 1].e == t.s
                if t.kind == "num" and not adj_dot:
                    neg = (p is not None and p.text == "-" and p.e == t.s and
                           (prev_tok(toks, p.i) is None or prev_tok(toks, p.i).text in UNARY_CTX))
                    span = (p.s, t.e) if neg else (t.s, t.e)
                    events.append((k, "LIT", (span, -int(t.text) if neg else int(t.text))))
                elif t.kind == "num" and adj_dot and t.text in PROJPAIRS:
                    events.append((k, "VAR", ((t.s, t.e), t.text, PROJPAIRS[t.text])))
                elif t.kind == "ident" and t.text in ("true", "false"):
                    events.append((k, "BLIT", ((t.s, t.e), t.text)))
                elif t.text in REL_VARIANTS:
                    events.append((k, "REL", ((t.s, t.e), t.text)))
                elif t.kind == "ident" and t.text in ("eqZ", "neZ") and not adj_dot:
                    events.append((k, "REL", ((t.s, t.e), t.text)))
                elif t.kind == "ident" and t.text in ("isSome", "isNone") and adj_dot:
                    events.append((k, "REL", ((t.s, t.e), t.text)))
                elif t.kind == "ident" and t.text in ("all", "any") and adj_dot:
                    events.append((k, "QNT", ((t.s, t.e), t.text)))
                elif t.text in ("&&", "||", "∧", "∨"):
                    events.append((k, "LOG", None))
                elif t.text == "!" or (t.text == "-" and not is_binary(toks, k)):
                    nx = toks[k + 1]
                    if t.text == "-" and nx.kind == "num" and nx.s == t.e:
                        continue                                   # sign of a negative literal (LIT site)
                    events.append((k, "NEG", ((t.s, t.e), t.text)))
                elif t.text in AOPS and is_binary(toks, k):
                    events.append((k, "AOR", None))
                elif t.kind == "ident" and t.text in ("min", "max") and not adj_dot:
                    events.append((k, "AOR", None))
                elif t.kind == "ident" and t.text == "if":
                    events.append((k, "IF", None))
                elif t.kind == "ident" and t.text == "match":
                    j = k + 1
                    d = 0
                    while not (d == 0 and toks[j].text == "with"):
                        d += toks[j].text in OPEN
                        d -= toks[j].text in CLOSE
                        j += 1
                    arms, _ = parse_arms(toks, next_tok_idx(toks, j))
                    events.append((k, "MATCH", ("match", arms)))
                elif t.text == "[" and toks[k + 1].text != "]":
                    events.append((k, "LIST", None))
                elif t.kind == "ident" and k not in btok:
                    if t.text in VARPAIRS and not adj_dot and in_scope(scopes, VARPAIRS[t.text], k) \
                            and in_scope(scopes, t.text, k):
                        events.append((k, "VAR", ((t.s, t.e), t.text, VARPAIRS[t.text])))
                    elif t.text in FIELDPAIRS and adj_dot:
                        events.append((k, "VAR", ((t.s, t.e), t.text, FIELDPAIRS[t.text])))
                    elif t.text in GLOBALPAIRS and not adj_dot:
                        events.append((k, "VAR", ((t.s, t.e), t.text, GLOBALPAIRS[t.text])))
            # chains of ++ and ::
            for opx, grp_key in (("++", "CAT"), ("::", "CONS")):
                chains = OrderedDict()
                for k, t in enumerate(toks):
                    if t.text == opx and k not in inpat:
                        L, R = left_operand(toks, k), right_operand(toks, k)
                        key = L[0] if opx == "++" else R[1]
                        chains.setdefault(key, []).append((k, L, R))
                for key, ops in chains.items():
                    ops.sort()
                    events.append((ops[0][1][0], grp_key, ops))
            events.sort(key=lambda e: e[0])
            for k, op, pay in events:
                t = toks[k]
                if op == "LIT":
                    span, v = pay
                    add("LIT", dn, body, k, span, value=v)
                elif op == "BLIT":
                    add("BLIT", dn, body, k, pay[0], value=pay[1])
                elif op == "REL":
                    add("REL", dn, body, k, pay[0], rel=pay[1])
                elif op == "LOG":
                    L, R = left_operand(toks, k), right_operand(toks, k)
                    add("LOG", dn, body, k, (t.s, t.e), conn=t.text, L=L, R=R,
                        node=(toks[L[0]].s, toks[R[1] - 1].e))
                elif op == "QNT":
                    add("QNT", dn, body, k, pay[0], q=pay[1])
                elif op == "NEG":
                    add("NEG", dn, body, k, pay[0], sym=pay[1])
                elif op == "AOR":
                    if t.kind == "ident":
                        add("AOR", dn, body, k, (t.s, t.e), aop=t.text)
                    else:
                        L, R = left_operand(toks, k), right_operand(toks, k)
                        add("AOR", dn, body, k, (t.s, t.e), aop=t.text, L=L, R=R)
                elif op == "IF":
                    c, th, el = parse_if(toks, k)
                    add("IF", dn, body, k, (t.s, toks[el[1] - 1].e), then=th, els=el)
                elif op == "MATCH":
                    kind, arms = pay
                    last = arms[-1][2]
                    add("MATCH", dn, body, k, (t.s, toks[last[1] - 1].e), arms=arms, kind=kind)
                elif op == "LIST":
                    j, d, ents, st = k + 1, 0, [], k + 1
                    while True:
                        u = toks[j]
                        if u.text in OPEN:
                            d += 1
                        elif u.text in CLOSE:
                            if d == 0:
                                break
                            d -= 1
                        elif d == 0 and u.text == ",":
                            ents.append(trim(toks, st, j))
                            st = j + 1
                        j += 1
                    ents.append(trim(toks, st, j))
                    add("LIST", dn, body, k, (t.s, toks[j].e), entries=ents)
                elif op == "VAR":
                    span, frm, to = pay
                    add("VAR", dn, body, k, span, frm=frm, to=to, group=PAIRGROUP[frm])
                elif op in ("CAT", "CONS"):
                    ops = pay
                    if op == "CAT":
                        start, end = ops[0][1][0], max(R[1] for _, _, R in ops)
                    else:
                        start, end = min(L[0] for _, L, _ in ops), ops[-1][2][1]
                    cuts = [kk for kk, _, _ in ops]
                    segs, st = [], start
                    for kk in cuts:
                        segs.append(trim(toks, st, kk))
                        st = kk + 1
                    segs.append(trim(toks, st, end))
                    add(op, dn, body, ops[0][0], (toks[start].s, toks[end - 1].e), segs=segs)
    return sites


# ----------------------------------------------------------------------------------------------
# mutations
# ----------------------------------------------------------------------------------------------

class Mut:
    __slots__ = ("id", "op", "defname", "hooks", "spans", "edits")

    def __init__(self, mid, op, defname, hooks, spans, edits):
        self.id, self.op, self.defname = mid, op, defname
        self.hooks, self.spans, self.edits = hooks, tuple(spans), edits


def lit_text(v):
    return str(v) if v >= 0 else f"({v})"


def build_mutations(sites, ms, lit_type):
    muts = []

    def tx(body, r):
        return ms[body.toks[r[0]].s:body.toks[r[1] - 1].e]

    def wrap_edits(body, L, k, R, newop):
        toks = body.toks
        lm, rm = (L[1] - L[0]) > 1, (R[1] - R[0]) > 1
        ed = []
        if lm:
            ed.append((toks[L[0]].s, toks[L[0]].s, "("))
        ed.append((toks[k].s, toks[k].e, (") " if lm else "") + newop + (" (" if rm else "")))
        if rm:
            ed.append((toks[R[1] - 1].e, toks[R[1] - 1].e, ")"))
        return ed

    by_def = defaultdict(list)
    for s in sites:
        by_def[s.defname].append(s)
        sp = [s.span]
        if s.op == "LIT":
            v = s.info["value"]
            for nv in (v - 1, v + 1) + ((-v,) if v < 0 else ()):
                tag = "NEG" if nv == -v and v < 0 else "LIT"
                muts.append(Mut(f"{s.id}={nv}", tag, s.defname, {s.id: nv}, sp, [(s.span[0], s.span[1], lit_text(nv))]))
        elif s.op == "BLIT":
            nv = {"true": "false", "false": "true"}[s.info["value"]]
            muts.append(Mut(f"{s.id}={nv}", "BLIT", s.defname, {s.id: nv == "true"}, sp, [(s.span[0], s.span[1], nv)]))
        elif s.op == "REL":
            for tag, new in REL_VARIANTS[s.info["rel"]]:
                muts.append(Mut(f"{s.id}={tag}", "REL", s.defname, {s.id: tag}, sp, [(s.span[0], s.span[1], new)]))
        elif s.op == "QNT":
            new = {"all": "any", "any": "all"}[s.info["q"]]
            muts.append(Mut(f"{s.id}={new}", "QNT", s.defname, {s.id: new}, sp, [(s.span[0], s.span[1], new)]))
        elif s.op == "LOG":
            b, L, R, k = s.body, s.info["L"], s.info["R"], s.k
            dual = {"&&": "||", "||": "&&", "∧": "∨", "∨": "∧"}[s.info["conn"]]
            muts.append(Mut(f"{s.id}=swap", "LOG", s.defname, {s.id: "swap"}, sp, wrap_edits(b, L, k, R, dual)))
            node = [s.info["node"]]
            for tag, keep in (("dropL", R), ("dropR", L)):
                t = tx(b, keep)
                t = t if keep[1] - keep[0] == 1 else f"({t})"
                muts.append(Mut(f"{s.id}={tag}", "LOG", s.defname, {s.id: tag}, node,
                                [(node[0][0], node[0][1], t)]))
        elif s.op == "NEG":
            muts.append(Mut(f"{s.id}=del", "NEG", s.defname, {s.id: "del"}, sp, [(s.span[0], s.span[1], "")]))
        elif s.op == "AOR":
            a = s.info["aop"]
            if a in ("min", "max"):
                new = {"min": "max", "max": "min"}[a]
                muts.append(Mut(f"{s.id}={new}", "AOR", s.defname, {s.id: new}, sp, [(s.span[0], s.span[1], new)]))
            else:
                for new in AOPS:
                    if new != a:
                        muts.append(Mut(f"{s.id}={new}", "AOR", s.defname, {s.id: new}, sp,
                                        wrap_edits(s.body, s.info["L"], s.k, s.info["R"], new)))
        elif s.op == "IF":
            b, th, el = s.body, s.info["then"], s.info["els"]
            st, en = b.toks[th[0]].s, b.toks[el[1] - 1].e
            muts.append(Mut(f"{s.id}=swap", "IF", s.defname, {s.id: "swap"}, sp,
                            [(st, en, f"({tx(b, el)}) else ({tx(b, th)})")]))
        elif s.op == "MATCH":
            arms = s.info["arms"]
            for x, y in itertools.combinations(range(len(arms)), 2):
                rx, ry = arms[x][2], arms[y][2]
                b = s.body
                ed = [(b.toks[rx[0]].s, b.toks[rx[1] - 1].e, tx(b, ry)),
                      (b.toks[ry[0]].s, b.toks[ry[1] - 1].e, tx(b, rx))]
                muts.append(Mut(f"{s.id}=swap{x}{y}", "MATCH", s.defname, {s.id: f"swap{x}{y}"}, sp, ed))
        elif s.op in ("LIST", "CAT", "CONS"):
            b = s.body
            parts = s.info["entries"] if s.op == "LIST" else s.info["segs"]
            texts = [tx(b, r) for r in parts]
            nel = len(texts) if s.op != "CONS" else len(texts) - 1
            for x in range(nel):
                for tag in ("drop", "dup"):
                    if tag == "drop":
                        new = texts[:x] + texts[x + 1:]
                    else:
                        new = texts[:x + 1] + [texts[x]] + texts[x + 1:]
                    if s.op == "LIST":
                        t = "[" + ", ".join(new) + "]"
                    elif s.op == "CAT":
                        t = " ++ ".join(new)
                    else:
                        t = " :: ".join(new)
                    muts.append(Mut(f"{s.id}={tag}{x}", s.op, s.defname, {s.id: (tag, x)}, sp,
                                    [(s.span[0], s.span[1], t)]))
        elif s.op == "VAR":
            muts.append(Mut(f"{s.id}={s.info['to']}", "VAR", s.defname, {s.id: s.info["to"]}, sp,
                            [(s.span[0], s.span[1], s.info["to"])]))
    # literal swaps within a definition (same type, different values)
    for dn, ss in by_def.items():
        lits = [s for s in ss if s.op == "LIT"]
        for a, b in itertools.combinations(lits, 2):
            va, vb = a.info["value"], b.info["value"]
            if va != vb and lit_type[a.id] == lit_type[b.id]:
                muts.append(Mut(f"{dn}:LSW:{a.id.split(':')[-1]},{b.id.split(':')[-1]}", "LSW", dn,
                                {a.id: vb, b.id: va}, [a.span, b.span],
                                [(a.span[0], a.span[1], lit_text(vb)), (b.span[0], b.span[1], lit_text(va))]))
        # end-label swap: every occurrence of a declared pair in one definition
        groups = defaultdict(list)
        for s in ss:
            if s.op == "VAR":
                groups[s.info["group"]].append(s)
        for g, vs in groups.items():
            if len(vs) >= 2:
                muts.append(Mut(f"{dn}:ENDSWAP:{g}", "ENDSWAP", dn, {s.id: s.info["to"] for s in vs},
                                [s.span for s in vs], [(s.span[0], s.span[1], s.info["to"]) for s in vs]))
    return muts


def apply_edits(src, edits):
    out = src
    for s, e, t in sorted(edits, key=lambda x: (-x[0], -x[1])):
        out = out[:s] + t + out[e:]
    return out


def compatible(m1, m2):
    for a, b in m1.spans:
        for c, d in m2.spans:
            if a < d and c < b:
                return False
            if a == b == c == d:
                return False
    return True


def match_swap_unbound(site):
    """for a MATCH site: does swapping arm right-hand sides move a pattern variable out of its arm?"""
    toks = site.body.toks
    arms = site.info["arms"]
    pv = []
    for bar, arrow, rhs in arms:
        pv.append({toks[j].text for j in range(bar + 1, arrow)
                   if toks[j].kind == "ident" and toks[j].text not in ("none", "some", "true", "false", "_")})
    used = [{toks[j].text for j in range(r[0], r[1]) if toks[j].kind == "ident"} for _, _, r in arms]
    return any((used[y] & pv[y]) - pv[x] or (used[x] & pv[x]) - pv[y]
               for x, y in itertools.combinations(range(len(arms)), 2))


def ill_typed(m, sites_by_id, int_lits, hetero):
    """predicted elaboration failure of a single mutation (checked against Lean on every single mutant)"""
    for sid, v in m.hooks.items():
        s = sites_by_id[sid]
        if s.op == "LIT" and sid not in int_lits and v < 0:
            return "negative Nat literal"
        if s.op == "VAR" and sid in hetero:
            return "projection swap on Int x Assign"
        if s.op == "MATCH" and match_swap_unbound(s):
            return "match arm uses a variable of the other pattern"
    return None








# ----------------------------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------------------------
import argparse
import ast
import math
import multiprocessing as mpc
import os
import random
import mpmath as mp

import chainalpha_transcription as T

HERE = os.path.dirname(os.path.abspath(__file__))
mp.mp.dps = 40
INV_CODATA = mp.mpf("137.035999177")          # CODATA 2022, 1/alpha = 137.035 999 177(21)
SIG = mp.mpf("2.1e-8")
SEED = 20261030
N_PAIR_SAMPLE = 300
ROOTS = ["readout", "objectWorld", "objectTwoPoint"]
OP_ORDER = ["LIT", "NEG", "LSW", "BLIT", "REL", "QNT", "LOG", "IF", "MATCH", "AOR", "VAR", "ENDSWAP", "LIST", "CAT",
            "CONS"]
OP_DESC = {
    "LIT": "integer literal -> value - 1, value + 1",
    "NEG": "delete a unary minus or Boolean not; sign of a negative literal",
    "LSW": "swap two literals of one type and different value within one definition",
    "BLIT": "Boolean literal true <-> false",
    "REL": "negate a test (= / ≠, == / !=, eqZ, neZ, isSome / isNone); < -> ≤, ≥ and ≤ -> <, >",
    "QNT": "list quantifier all <-> any",
    "LOG": "&& <-> ||, ∧ <-> ∨; drop the left or the right operand",
    "IF": "swap the two branches of an if-then-else",
    "MATCH": "swap the right-hand sides of two match arms",
    "AOR": "binary + - * / % -> each of the other four; min <-> max",
    "VAR": "one occurrence of a declared label -> its partner (end/role pairs, .1 <-> .2, readEnd <-> farEnd)",
    "ENDSWAP": "every occurrence of one declared label pair in one definition exchanged",
    "LIST": "list literal: drop or duplicate an entry",
    "CAT": "++ chain: drop or duplicate a summand",
    "CONS": ":: chain: drop or duplicate an element",
}


def load_lean(path):
    with open(path, "rb") as f:
        raw = f.read()
    h = hashlib.sha256(raw).hexdigest()
    assert h == PIN, f"sha256 mismatch: {h}"
    return raw.decode("utf-8")


def setup(path):
    src = load_lean(path)
    ms = mask_comments(src)
    assert len(ms) == len(src)
    defs = bodies_of(parse_items(src, ms))
    names = cone(defs, ROOTS)
    order = [n for n in defs if n in names]
    sites = find_sites(defs, order, ms)
    byid = {s.id: s for s in sites}
    lit_type = {s.id: ("I" if s.id in T.INT_LITS else "N") for s in sites if s.op == "LIT"}
    muts = build_mutations(sites, ms, lit_type)
    ill = [ill_typed(m, byid, T.INT_LITS, T.HETERO_PROJ) for m in muts]
    return src, ms, defs, order, sites, byid, muts, ill


# ---------- exact values ----------

def sqfree(k):
    s, m, p = 1, k, 2
    while p * p <= m:
        while m % (p * p) == 0:
            m //= p * p
            s *= p
        p += 1
    return s, m


def canonical(ro):
    """exact class of sum(c sqrt k)/N: ('N=0',), ('zero',) or (N', ((m, a'), ...)) reduced, m squarefree"""
    if isinstance(ro, str):
        return (ro,)
    lst, N = ro
    coef = {}
    for k, c in lst:
        if k == 0 or c == 0:
            continue
        s, m = sqfree(k)
        coef[m] = coef.get(m, 0) + c * s
    items = sorted((m, a) for m, a in coef.items() if a != 0)
    if N == 0:
        return ("N=0",)
    if not items:
        return ("zero",)
    g = math.gcd(N, *[a for _, a in items])
    return (N // g, tuple((m, a // g) for m, a in items))


def inv_value(key):
    N, items = key
    return mp.mpf(N) / mp.fsum(a * mp.sqrt(m) for m, a in items)


def fmt_key(key):
    if len(key) == 1:
        return key[0]
    N, items = key
    parts = []
    for m, a in items:
        t = f"{a}" if m == 1 else f"{a}√{m}"
        parts.append(t)
    return "(" + " + ".join(parts).replace("+ -", "- ") + f")/{N}"


BASE_KEY = canonical(T.BASE_READOUT)


# ---------- evaluation (worker side) ----------

_W = {}


def _init_worker(hook_list, ill_list):
    _W["hooks"], _W["ill"] = hook_list, ill_list


def _eval_pairs(chunk):
    hooks, ill = _W["hooks"], _W["ill"]
    out, longest = [], 0
    for a, b in chunk:
        if ill[a] or ill[b]:
            out.append(("ill-typed",))
        else:
            h = dict(hooks[a])
            h.update(hooks[b])
            ro = T.evaluate(h)
            if not isinstance(ro, str):
                longest = max(longest, len(ro[0]))
            out.append(canonical(ro))
    return out, longest


# ---------- statistics ----------

def uniform_target_count(values, lo, hi):
    """Illustrative model of freedom.py: U uniform on (lo,hi), the outputs `values` fixed. E[#{v : |v-U| <= SIG}] is
    the total length of the intervals [v-SIG, v+SIG] clipped to (lo,hi), divided by hi-lo. No independence between
    the intervals is needed."""
    lo, hi = mp.mpf(lo), mp.mpf(hi)
    return mp.fsum(max(mp.mpf(0), min(hi, v + SIG) - max(lo, v - SIG)) for v in values) / (hi - lo)


def stats(keys, label):
    """keys: list of canonical keys (one per mutant). Prints the F-CT-4 metrics for the well-defined ones."""
    good = [k for k in keys if len(k) == 2]
    distinct = sorted(set(good))
    inv = {k: inv_value(k) for k in distinct}
    nbase = sum(1 for k in good if k == BASE_KEY)
    neg = sum(1 for k in distinct if inv[k] < 0)
    print(f"  {label}: {len(keys)} mutants, {len(good)} well-defined, {nbase} equal to the Letter's value, "
          f"{len(distinct)} distinct values ({neg} negative)")
    band = sorted(inv[k] for k in distinct if 120 < inv[k] < 155)
    dens = mp.mpf(len(band)) / 35
    ref = uniform_target_count(inv.values(), 120, 155)
    print(f"    distinct values with 1/value in (120,155): {len(band)}; density {mp.nstr(dens, 4)} per unit")
    print("    illustrative model: target U uniform on (120,155); rule outputs fixed")
    print(f"    expected match count in that model within +-1 sigma: {mp.nstr(ref, 4)}")
    print("    this reference count is not a calibrated p-value for the observed alpha or a search-history estimate")
    for lo, hi in ((130, 145), (136, 138)):
        nb = sum(1 for k in distinct if lo < inv[k] < hi)
        print(f"    distinct values in ({lo},{hi}): {nb}; density {mp.nstr(mp.mpf(nb) / (hi - lo), 4)} per unit; "
              f"same model with U uniform on ({lo},{hi}): {mp.nstr(uniform_target_count(inv.values(), lo, hi), 4)}")
    for w in ("1e-2", "1e-3", "1e-4", "1e-6"):
        n = sum(1 for k in distinct if abs(inv[k] - INV_CODATA) / INV_CODATA < mp.mpf(w))
        print(f"    distinct values within relative {w} of 1/alpha: {n}")
    if len(band) > 1:
        gaps = sorted(b - a for a, b in zip(band, band[1:]))
        print(f"    median gap between adjacent distinct values in (120,155): {mp.nstr(gaps[len(gaps) // 2], 4)}")
    return distinct, inv, dens, ref


def closest(distinct, inv, who, n=10):
    rows = sorted(((abs(inv[k] - INV_CODATA) / SIG, k) for k in distinct), key=lambda t: t[0])[:n]
    for d, k in rows:
        ex = who[k]
        print(f"    {mp.nstr(d, 6):>12s} sigma  1/value {mp.nstr(inv[k], 15):>18s}  {fmt_key(k):40s} "
              f"{len(ex)} mutant(s), e.g. {ex[0]}")


def author_rows():
    pat = re.compile(r"^(\S+)\s+main\s+(proved|FAILS)\s+commutators\s+(proved|FAILS)\s+readout\s+(\(.*\))\s*$")
    rows = []
    with open(os.path.join(HERE, "upstream_mutations_lean.out")) as f:
        for line in f:
            m = pat.match(line.rstrip("\n"))
            rows.append((m.group(1), line.rstrip("\n"), ast.literal_eval(m.group(4))))
    return rows


# the author's 19 rewordings (mutations.py), as hooks of the transcription. Five are elements of the mechanical space
# (their hooks are those of one single or one pair); the others are written as dedicated variants ('A:...').
AUTHOR = {
    "BASE": {},
    "fixesFar_always_true": {"A:fixesFar_always_true": 1},
    "reads_at_qubit": {"readEnd:IF:0": "swap", "farEnd:IF:0": "swap"},
    "X_at_far_end": {"wX:VAR:0": "farEnd"},
    "X_where_letters_sit": {"A:X_where_letters_sit": 1},
    "Y_at_center": {"A:Y_at_center": 1},
    "X_Y_factors_swapped": {"A:X_Y_factors_swapped": 1},
    "moved_per_letter": {"A:moved_per_letter": 1},
    "copies_separate": {"A:copies_separate": 1},
    "survival_ignores_own_sign": {"A:survival_ignores_own_sign": 1},
    "commutator_records_carry_no_reads": {"A:commutator_records_carry_no_reads": 1},
    "register_full_unread": {"A:register_full_unread": 1},
    "no_commutator_records": {"A:no_commutator_records": 1},
    "commutators_keep_sums": {"A:commutators_keep_sums": 1},
    "commutators_keep_both": {"A:commutators_keep_both": 1},
    "exch_shift1": {"exch:LIT:0": 1},
    "drop_fixed": {"recs:CAT:0": ("drop", 2)},
    "drop_even": {"recs:CAT:0": ("drop", 1)},
    "born_weights": {"A:born_weights": 1},
    "wX_dA_equiv": {"A:wX_dA_equiv": 1},
}


def pair_sample(muts, pairs):
    rng = random.Random(SEED)
    return sorted(rng.sample(range(len(pairs)), N_PAIR_SAMPLE))


def all_pairs(muts):
    n = len(muts)
    return [(a, b) for a in range(n) for b in range(a + 1, n) if compatible(muts[a], muts[b])]


def lean_text(src, edits):
    cut = src.index("/-! ## 1.9")
    assert all(e[1] <= cut for e in edits)
    return apply_edits(src[:cut], edits) + "\n#eval readout objectWorld objectTwoPoint\nend ChainAlpha\n"


def expected_str(key_ro):
    return key_ro if isinstance(key_ro, str) else repr(key_ro)


def emit(src, muts, ill, pairs, outdir):
    os.makedirs(outdir, exist_ok=True)
    sample = pair_sample(muts, pairs)
    with open(os.path.join(outdir, "expected.tsv"), "w") as f:
        for i, m in enumerate(muts):
            with open(os.path.join(outdir, f"s{i:04d}.lean"), "w") as g:
                g.write(lean_text(src, m.edits))
            exp = "ILL" if ill[i] else expected_str(T.evaluate(dict(m.hooks)))
            f.write(f"s{i:04d}\tsingle\t{m.id}\t{exp}\n")
        for x, p in enumerate(sample):
            a, b = pairs[p]
            with open(os.path.join(outdir, f"p{x:04d}.lean"), "w") as g:
                g.write(lean_text(src, muts[a].edits + muts[b].edits))
            if ill[a] or ill[b]:
                exp = "ILL"
            else:
                h = dict(muts[a].hooks)
                h.update(muts[b].hooks)
                exp = expected_str(T.evaluate(h))
            f.write(f"p{x:04d}\tpair\t{muts[a].id} & {muts[b].id}\t{exp}\n")
    print(f"wrote {len(muts)} single and {len(sample)} pair mutants to {outdir}")


def compare(outdir):
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(outdir, "expected.tsv"))]
    tally = defaultdict(lambda: [0, 0])
    bad, beyond = [], []
    resource = Counter()
    for f, kind, mid, exp in rows:
        out = open(os.path.join(outdir, f + ".out")).read()
        ex = re.search(r"exit (\d+)\s*$", out)
        code = int(ex.group(1)) if ex else None
        lines = out.split("\n")
        st = next((x for x, l in enumerate(lines) if l.startswith("([")), None)
        if code == 124 or "deep recursion" in out or "out of memory" in out.lower():
            got = "RESOURCE"
            resource[("timeout" if code == 124 else "deep recursion" if "deep recursion" in out else "out of memory")] += 1
        elif "error" in out:
            got = "ILL"
        elif code == 0 and st is not None:
            got = "".join(l for l in lines[st:] if not l.startswith("exit"))   # #eval output, possibly wrapped
        else:
            got = "?"
        cls = exp if exp in ("ILL", "INFEASIBLE") else "value"
        nows = lambda x: re.sub(r"\s+", "", x)
        if exp == "INFEASIBLE" and got.startswith("(["):
            beyond.append((f, mid, len(ast.literal_eval(got)[0])))
            continue
        ok = (got == "RESOURCE") if exp == "INFEASIBLE" else (nows(got) == nows(exp))
        tally[(kind, cls)][0 if ok else 1] += 1
        if not ok:
            bad.append((f, mid, exp, got))
    print("Lean 4.35.0 against the Python transcription (mechanical_mutations.py --emit-lean):")
    print("  Lean result counted as: timeout (120 s) or resource abort (deep recursion, out of memory) -> infeasible; "
          "else an error message -> fails to elaborate; else the #eval output (compared without whitespace)")
    for (kind, cls), (y, n) in sorted(tally.items()):
        lab = {"value": "Python value", "ILL": "Python: fails to elaborate", "INFEASIBLE": "Python: infeasible"}[cls]
        print(f"  {kind:6s} {lab:30s} agree {y:5d}  disagree {n}")
    tot = sum(y for y, _ in tally.values()), sum(n for _, n in tally.values())
    print(f"  total: {tot[0]} agree, {tot[1]} disagree")
    print("  Lean resource failures by kind: " + (", ".join(f"{v} {k}" for k, v in sorted(resource.items())) or "none"))
    print(f"  beyond the transcription's bound (Python: more than {T.MAX_TERMS} terms, not evaluated) but evaluated by "
          f"Lean: {len(beyond)}" + "".join(f"; {f} {mid} ({n} terms)" for f, mid, n in beyond))
    for f, mid, exp, got in bad:
        print(f"  DISAGREE {f} {mid}: python {exp[:120]} lean {got[:120]}")


def main_mech_names(src, defs):
    a = src.index("theorem main :")
    b = src.index(":= by", a)
    toks = tokenize(mask_comments(src), a, b)
    return sorted({t.text for t in toks if t.kind == "ident" and t.text in defs} |
                  {"TwoPoint.valid" for t in toks if t.text == "valid"})


def part0_defs(src):
    """Part 0 definitions: qualified name (without the file's namespace ChainAlpha) -> identifiers of its body"""
    a = src.index("/-! # Part 1")
    ms = mask_comments(src[:a])
    out, stack, cur = {}, [], None
    for line in ms.split("\n"):
        m = re.match(r"(namespace|end|def) ([\w.]+)", line)
        if m and m.group(1) == "namespace":
            stack.append(m.group(2))
            cur = None
            continue
        if m and m.group(1) == "end":
            stack.pop()
            cur = None
            continue
        if m and m.group(1) == "def":
            cur = ".".join([x for x in stack if x != "ChainAlpha"] + [m.group(2)])
            out[cur] = set()
            line = line[m.end():]
        if cur:
            ns = ".".join(x for x in stack if x != "ChainAlpha")
            for t in re.findall(r"[A-Za-z_][\w.]*", line):
                out[cur].add(t)
                if ns:
                    out[cur].add(ns + "." + t)
    return out


def report(src, ms, defs, order, sites, byid, muts, ill, pairs):
    print("Mechanically generated rule family for the readout of ChainAlpha.lean: rule-family sensitivity (audit F-CT-4 E)")
    print(f"Source: ChainAlpha.lean, sha256 {PIN} (pinned upstream file; read as data, never run here).")
    print("\n== 1. Inventory and declared space (printed in this run before any value is computed) ==")
    p1 = [n for n, v in defs.items() if v[0] == "def"]
    print(f"Part 1 has {len(p1)} definitions and {len(defs) - len(p1)} structures/abbreviations.")
    print(f"The readout `readout objectWorld objectTwoPoint` (the last conjunct of `main`) depends on {len(order)} of "
          f"them (identifier references, transitively, from {', '.join(ROOTS)}):")
    for n in order:
        c = Counter(s.op for s in sites if s.defname == n)
        print(f"  {n:22s} " + ", ".join(f"{op} {c[op]}" for op in OP_ORDER if c[op]))
    c = Counter(s.op for s in sites)
    print(f"  total: {len(sites)} sites = " + " + ".join(f"{c[op]} {op}" for op in OP_ORDER if c[op]))
    p0d = part0_defs(src)
    p0 = set(p0d)
    used0 = set()
    for n in order:
        for b in defs[n][1]:
            for k, t in enumerate(b.toks):
                if t.kind != "ident":
                    continue
                nm = t.text
                if k >= 2 and b.toks[k - 1].text == "." and b.toks[k - 2].kind == "ident":
                    nm2 = b.toks[k - 2].text + "." + nm
                    if nm2 in p0:
                        used0.add(nm2)
                if nm in p0:
                    used0.add(nm)
    todo, clo = sorted(used0), set()
    while todo:
        x = todo.pop()
        if x in clo:
            continue
        clo.add(x)
        todo += [y for y in p0d[x] if y in p0 and y not in clo]
    print(f"Held fixed: the Part 0 utilities the cone calls (exact arithmetic): {', '.join(sorted(used0))}; "
          f"with what they call: {', '.join(sorted(clo - used0))}.")
    extra = [n for n in main_mech_names(src, defs) if n not in order]
    xc = [n for n in cone(defs, extra) if n not in order and defs[n][0] == "def"]
    print(f"In `main` but outside the readout (cannot change it): {', '.join(xc)}.")
    print("Types not mutated: " + ", ".join(n for n, v in defs.items() if v[0] != "def") + ".")
    print("\nDeclared operators (each applied at every site the scanner finds in the 44 bodies):")
    for op in OP_ORDER:
        print(f"  {op:8s} {OP_DESC[op]}")
    print("  label pairs (VAR, ENDSWAP): " + ", ".join(f"{a}/{b}" for a, b in _VP) +
          ", read1/read2, readEnd/farEnd, .1/.2 (a label is replaced only where its partner is in scope)")
    print("  pairs of mutations: every unordered pair whose sites are disjoint (a site is the token it rewrites, "
          "or the whole if / match / list / chain / connective node it moves or drops)")
    c = Counter(m.op for m in muts)
    ci = Counter(m.op for m, i in zip(muts, ill) if i)
    print(f"\nSingle mutants: {len(muts)} = " + " + ".join(f"{c[op]} {op}" for op in OP_ORDER if c[op]))
    print(f"  predicted to fail to elaborate: {sum(1 for i in ill if i)} "
          f"({', '.join(f'{v} {k}' for k, v in sorted(Counter(i for i in ill if i).items()))}); "
          f"by operator: " + ", ".join(f"{op} {ci[op]}" for op in OP_ORDER if ci[op]))
    n = len(muts)
    npw = sum(1 for a, b in pairs if not ill[a] and not ill[b])
    print(f"Pairs: {n * (n - 1) // 2} unordered pairs, {len(pairs)} with disjoint sites; "
          f"{len(pairs) - npw} of these contain a mutant that fails to elaborate; {npw} are evaluated.")
    print(f"Space: {n} + {len(pairs)} = {n + len(pairs)} mutants.")

    print("\n== 2. Validation of the Python transcription (chainalpha_transcription.py) ==")
    tsrc = open(T.__file__).read()
    hooks_used = set(re.findall(r'"([A-Za-z]+:[A-Z]+:\d+)"', tsrc)) | T.INT_LITS | T.HETERO_PROJ
    need = {s.id for s in sites if s.op != "MATCH"}
    print(f"  hooks named in the transcription (with its tables INT_LITS, HETERO_PROJ) = sites of the generator other "
          f"than MATCH: {hooks_used == need} ({len(need)} sites; the {sum(1 for s in sites if s.op == 'MATCH')} MATCH "
          f"sites need no hook: every match swap moves a pattern variable out of its arm and fails to elaborate)")
    assert hooks_used == need
    assert T.INT_LITS <= {s.id for s in sites if s.op == "LIT"} and T.HETERO_PROJ <= {s.id for s in sites if s.op == "VAR"}
    assert all(match_swap_unbound(s) for s in sites if s.op == "MATCH")
    m = re.search(r"readout objectWorld objectTwoPoint = (\(\[.*?\], \d+\)) := by", src)
    lean_main = ast.literal_eval(m.group(1))
    base = T.evaluate({})
    print(f"  unmutated readout: {base}; the statement of `dressing_of_the_qutrit_qubit_link`: {lean_main}; "
          f"equal: {base == lean_main}")
    assert base == lean_main == T.BASE_READOUT
    rows = author_rows()
    same = 0
    for name, line, ro in rows:
        mm, cc, ro2 = T.main_and_comm(AUTHOR[name])
        mine = f"{name:24s} main {'proved' if mm else 'FAILS '}  commutators {'proved' if cc else 'FAILS '}  readout {ro2!r}"
        same += mine == line
        if mine != line:
            print(f"  MISMATCH {name}: {mine}")
    print(f"  upstream_mutations_lean.out (Lean 4.35.0 on the author's 19 rewordings + BASE): {same} of {len(rows)} "
          f"lines reproduced character for character (main / commutators verdicts and readout)")
    assert same == len(rows) == 20
    hk = {tuple(sorted(mu.hooks.items())): mu.id for mu in muts}
    for name in AUTHOR:
        h = tuple(sorted(AUTHOR[name].items()))
        if name == "BASE" or any(k.startswith("A:") for k, _ in h):
            continue
        if h in hk:
            print(f"  {name} is the single mutant {hk[h]}")
        else:
            parts = [mu.id for mu in muts if set(mu.hooks.items()) <= set(h)]
            print(f"  {name} is the pair {' & '.join(parts)}")
    print("  The other 14 rewordings change several tokens or add code; they are not single or paired mutants.")
    print("  Lean itself on the generated mutants: proofs/lean/chain-theory/mutation_spotcheck.sh -> "
          "mutation_spotcheck.out (all singles and a seeded sample of pairs).")

    print("\n== 3. Evaluation ==")
    single_keys, longest = [], 0
    for mu, i in zip(muts, ill):
        if i:
            single_keys.append(("ill-typed",))
            continue
        ro = T.evaluate(dict(mu.hooks))
        if not isinstance(ro, str):
            longest = max(longest, len(ro[0]))
        single_keys.append(canonical(ro))
    hooks = [tuple(mu.hooks.items()) for mu in muts]
    chunks = [pairs[x:x + 2000] for x in range(0, len(pairs), 2000)]
    pair_keys = []
    with mpc.get_context("fork").Pool(2, initializer=_init_worker, initargs=(hooks, ill)) as pool:
        for res, lg in pool.imap(_eval_pairs, chunks):
            pair_keys += res
            longest = max(longest, lg)
    cats = ("ill-typed", "INFEASIBLE", "N=0", "zero")
    names = {"ill-typed": "fail to elaborate", "INFEASIBLE": f"infeasible (readout list > {T.MAX_TERMS} terms)",
             "N=0": "readout dimension 0 (value undefined)", "zero": "numerator 0 (1/value undefined)"}
    print(f"  longest readout list of a feasible mutant: {longest} terms (infeasible: more than {T.MAX_TERMS})")
    dup = [x for x, mu in enumerate(muts) if mu.id == "addRoot:CONS:2=dup0"][0]
    ninf = sum(1 for (a, b), k in zip(pairs, pair_keys) if k == ("INFEASIBLE",))
    ninf_dup = sum(1 for (a, b), k in zip(pairs, pair_keys) if k == ("INFEASIBLE",) and dup in (a, b))
    print(f"  infeasible pairs: {ninf}, of which {ninf_dup} contain addRoot:CONS:2=dup0 (it duplicates every list "
          f"element it passes, so the list doubles at each record whose radicand exceeds the head's)")
    for label, ks in (("singles", single_keys), ("pairs", pair_keys)):
        cc = Counter(k[0] if len(k) == 1 else "ok" for k in ks)
        print(f"  {label}: {len(ks)} mutants: " + ", ".join(f"{cc[x]} {names[x]}" for x in cats) +
              f", {cc['ok']} well-defined")

    print("\n== 4. Results ==")
    print("Single mutants by operator (well-defined / distinct values / equal to the Letter's value):")
    for op in OP_ORDER:
        ks = [k for mu, k in zip(muts, single_keys) if mu.op == op]
        good = [k for k in ks if len(k) == 2]
        print(f"  {op:8s} {len(ks):4d} mutants, {len(good):4d} well-defined, {len(set(good)):4d} distinct, "
              f"{sum(1 for k in good if k == BASE_KEY):4d} unchanged")
    who = defaultdict(list)
    for mu, k in zip(muts, single_keys):
        who[k].append(mu.id)
    for (a, b), k in zip(pairs, pair_keys):
        if len(who[k]) < 3:
            who[k].append(f"{muts[a].id} & {muts[b].id}")
        else:
            who[k].append(None)
    print("\nF-CT-4 metrics (distinct exact values; CODATA 2022 1/alpha = 137.035999177, sigma = 2.1e-8):")
    stats(single_keys, "singles")
    stats(pair_keys, "pairs")
    dist, inv, dens, eh = stats(single_keys + pair_keys, "singles + pairs")
    print("  closest distinct values to CODATA 2022 (singles + pairs):")
    closest(dist, inv, who, 12)

    print("\nComparison with the hand-declared space B of F-CT-4 (freedom.py; rebuilt from mutations_classified.py):")
    import mutations_classified as MC
    bkeys = []
    for ks, kc, nc, power, nN, num in MC.B.values():
        bkeys.append(canonical(([(1, num[0]), (2, num[1]), (3, num[2]), (6, num[3])], MC.NS[nN])))
    assert len(set(bkeys)) == 15639
    bd, binv, bdens, beh = stats(bkeys, "B (one entry per distinct value)")
    sb, sm = set(bd), set(dist)
    both = sb & sm
    band = lambda S, I: {k for k in S if 120 < I[k] < 155}
    print(f"  distinct values in both spaces: {len(both)}; with 1/value in (120,155): "
          f"{len(band(both, inv))} (mechanical {len(band(sm, inv))}, B {len(band(sb, binv))})")
    uinv = dict(inv)
    uinv.update(binv)
    print(f"  union of the two spaces: density in (120,155) {mp.nstr(mp.mpf(len(band(sm, inv) | band(sb, binv))) / 35, 4)}"
          f" per unit; expected match count in the illustrative model above: "
          f"{mp.nstr(uniform_target_count(uinv.values(), 120, 155), 4)}")
    print("  the author's 19 rewordings (upstream_mutations_lean.out), by value:")
    sk, pk = set(single_keys), set(pair_keys)
    for name, line, ro in rows[1:]:
        k = canonical(ro)
        where = "single" if k in sk else "pair" if k in pk else "not reached"
        print(f"    {name:36s} 1/value {mp.nstr(inv_value(k), 10):>12s}  in B: {'yes' if k in sb else 'no '}  "
              f"mechanical: {where}" + (f" (e.g. {who[k][0]})" if k in who and where != "not reached" else ""))
    print("  the 15 closed-form alternatives of Table S5 (mutations_classified.S5), by value:")
    for name, vec, N, printed in MC.S5[1:]:
        k = canonical(([(1, vec[0]), (2, vec[1]), (3, vec[2]), (6, vec[3])], N))
        where = "single" if k in sk else "pair" if k in pk else "not reached"
        print(f"    {name:42s} 1/value {printed:>8.3f}  in B: {'yes' if k in sb else 'no '}  mechanical: {where}"
              + (f" (e.g. {who[k][0]})" if where != "not reached" else ""))
    # every mutant with the value of the Born-weight row (review R4-CT-3: the rows above print one example per value)
    vec, N = [(v, N) for name, v, N, printed in MC.S5[1:] if name == "Born weights k_r"][0]
    born = canonical(([(1, vec[0]), (2, vec[1]), (3, vec[2]), (6, vec[3])], N))
    bw = [mu.id for mu, k in zip(muts, single_keys) if k == born]
    nbs = len(bw)
    bw += [f"{muts[a].id} & {muts[b].id}" for (a, b), k in zip(pairs, pair_keys) if k == born]
    print(f"  every mutant with the value of 'Born weights k_r' ({fmt_key(born)}): {len(bw)} ({nbs} singles, "
          f"{len(bw) - nbs} pairs)")
    for x in bw:
        print(f"    {x}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lean", help="path of ChainAlpha.lean (default $UPSTREAM/papers/alpha/lean/ChainAlpha.lean)")
    ap.add_argument("--emit-lean", metavar="DIR", help="write every single and a seeded sample of pairs as Lean files")
    ap.add_argument("--compare", metavar="DIR", help="compare Lean outputs DIR/*.out with DIR/expected.tsv")
    args = ap.parse_args()
    if args.compare:
        compare(args.compare)
        return
    path = args.lean or os.path.join(os.environ["UPSTREAM"], "papers", "alpha", "lean", "ChainAlpha.lean")
    src, ms, defs, order, sites, byid, muts, ill = setup(path)
    pairs = all_pairs(muts)
    if args.emit_lean:
        emit(src, muts, ill, pairs, args.emit_lean)
        return
    report(src, ms, defs, order, sites, byid, muts, ill, pairs)


if __name__ == "__main__":
    main()
