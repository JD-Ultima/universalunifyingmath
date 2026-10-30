#!/usr/bin/env python3
"""Statement-level premises of the Chain Theory compendium (release 0.1), round 4.

Usage:
  python3 -I premises.py <compendium_dir> [--tsv out.tsv]

<compendium_dir> is the `compendium/` directory of the source repository
(commit a28db8998bce5bd5ef11195c05f13c86cd9dac5c). Everything is read as text; nothing in
the source is executed.

Why. The index premise lists (`premises` column; the "From ... Used by ..." line under every
box) are, by the source's own description, chapter-level: "one line names the chapters it rests
on and the chapters that rest on it, generated mechanically from the ancestry of the chapters"
(tex/parts/01-framework.tex:109). The statement-level premises are in the derivations
(`chproof` environments), each of which opens with a premise sentence ("The premises are ...",
"The premises of this derivation are ...", "\\textit{The premises.} ...").

Extraction rule (layer M, mechanical; every step below is implemented in this file):
  M1  Segmentation. Sections, boxes and derivations are numbered as in the index (same rule as
      boxes.py / depgraph.py). A box's derivation is the next `chproof` of its section; boxes
      between two derivations share the later one ("derivation group"). A `chproof` with no box
      before it continues the previous group's derivation.
  M2  Premise sentence. The first sentence (paren/brace-balanced) after "The premises are",
      "The premises of this derivation are"; or, after "\\textit{The premises.}", the text up to
      "\\textit{The steps", "Nothing outside", "The chain" or a blank line.
  M3  Items. Every chapter reference in the premise sentence is an item:
        REF  chapter~\\ref{mod:L} / Section~\\ref{mod:L}  -> the one section L (printed as its number)
        DES  "chapter P-c" (designator)                 -> reading MOD: the sections whose label
             begins "<roman P>-c-" (module numbering; "1-0" -> FRAMEWORK); reading FIF: the
             designators.tsv first_section..last_section range (the index's own mapping).
        SELF "this chapter", "the present chapter", "the ... chapter above" -> the boxes of the own
             module chapter other than the group's ("chapter above": only the boxes before the group);
             "the <grade> stated above" -> the boxes of the own module chapter that stand above the
             derivation, the group's own boxes included (13-thermodynamics.tex:45, "the parent theorem
             and the open result stated above", shared by 217.1-217.3: 217.3 is the open result)
        FRAME "the framework" -> section 1
        NOWHERE "(... outside the chapters of this compendium)", "(... placed at no chapter ...)"
      The item's grade word is (i) the grade phrase inside the item's own parenthesis, or
      (ii) "the theorem/open result/table of chapter~\\ref", or "the theorem/open result/table
      stated above" for a SELF back-reference, or (iii) ", <grade>" / " (<grade>)" right after
      the reference, or (iv) a plural grade ("definitions", "theorems", "open results") closing a
      run of ungraded references. Otherwise the item is ungraded. A conjunctive grade phrase that
      names an open result ("a theorem and an open result", "an open result and two theorems")
      is class TO: its candidates are THEOREM, CONDITIONAL and OPEN, and it declares an OPEN
      premise (the item names all the boxes it lists).
  M4  Candidates. The statements of the resolved sections whose index grade is compatible with
      the grade word: "theorem(s)" -> THEOREM or CONDITIONAL (the source writes "theorem" for
      conditional theorems, e.g. 9.1 is "a theorem" at 01-framework.tex:582 and "a conditional
      theorem" at :691); "conditional theorem" -> CONDITIONAL; "open result(s)" -> OPEN;
      "table"/"instrument" -> TABLE; "definition(s)"/"the text" -> DEFINITION (the axiom and
      definitions 1.1-1.6 when the description names the axiom, definitions or the text);
      ungraded -> every statement of the section.
  M5  Edges. candidate -> box for every box of the derivation group (self-loops dropped).
      An item with one candidate is "unique"; with several, an "any-of" item: the box rests on
      at least one of them (the text does not say which).

Layer J (human judgment, separate and labelled): for the derivation groups of the statements
under audit, each item is resolved to the statement(s) its description names, by the box
title or box text quoted in JUDGMENT below. Layer J never adds a premise the derivation does
not cite; it only narrows M's candidate sets.

Analyses printed (outputs/premises.txt):
  [0] the chain rule and the premise-line description, verbatim with file:line;
  [1] coverage of layer M, verbatim check of every item quote, evidence for the designator reading;
  [2] the index graph as in depgraph.py (8367 edges, one SCC of 22);
  [3]-[4] layer J and the statement-level premises of the reconstructed statements;
  [5] re-grade: a THEOREM violates the rule if one item's candidate set (self excluded) is
      non-empty and all CONDITIONAL/OPEN, or the item's own grade word is "open result" or a
      conjunction naming one (class TO) ("declared OPEN"); undetermined if an item mixes proven
      and unproven candidates;
  [6] cycles: the index two-cycles through 124.1 checked edge by edge in the derivations; SCCs of
      the any-of graphs (an upper bound); the strict graph (J choices, single-candidate items);
      cycles forced under every resolution of the any-of items (fixpoint below);
  [7] load-bearing premises: upper bound = reachability in the any-of graph; lower bound =
      "definite" dependents (fixpoint: a box is added when one item of its derivation has a
      non-empty candidate set inside the current set, so the dependence holds whichever
      candidate the text means);
  [8] all 100 THEOREMs at derivation level, with reading sensitivity.

Interpretive choices (see README.md, round 4): reading MOD for designators in derivation text
(evidence printed by this script); "theorem" read as THEOREM or CONDITIONAL; TABLE and
instruments treated as neutral (the front matter: "A presentation, such as a table, ... makes
no claim of its own").
"""
import collections, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))   # our own directory only
from boxes import balanced                                  # noqa: E402
from depgraph import load_tsv, parse_refs                     # noqa: E402

ENV = "chtheorem|chconditional|chwall|chdefinition|chtable"
GRADE = dict(chtheorem="THEOREM", chconditional="CONDITIONAL", chwall="OPEN",
             chdefinition="DEFINITION", chtable="TABLE")
G1 = dict(THEOREM="T", CONDITIONAL="C", OPEN="O", TABLE="Tab", DEFINITION="D")
ROMAN = {1: "I", 2: "II", 3: "III", 4: "IV", 5: "V", 6: "VI", 7: "VII", 8: "VIII", 9: "IX",
         10: "X", 11: "XI", 12: "XII", 13: "XIII", 14: "XIV", 15: "XV", 16: "XVI"}

# ----------------------------------------------------------------------------- M1 segmentation

def segment(texdir):
    """Sections (index numbering), boxes, derivations and derivation groups."""
    texdir = pathlib.Path(texdir)
    order = ["front"] + [p.stem for p in sorted((texdir/"parts").glob("*.tex"))] + ["back", "instruments"]
    tok = re.compile(r"\\part\*\{|\\appendix\b|\\section\{|\\begin\{(" + ENV + r"|chproof)\}")
    secnum, appendix, instr = 0, False, False
    sections = collections.OrderedDict()     # number -> dict(label, title, file, line)
    boxes = collections.OrderedDict()        # number -> dict(env, grade, title, file, line, section)
    groups = []                              # dict(section, boxes, proof, file, line, lines)
    files = {}
    for f in order:
        p = texdir/(f + ".tex") if f in ("front", "back", "instruments") else texdir/"parts"/(f + ".tex")
        s = p.read_text(encoding="utf-8"); files[p.name] = s
        pos, cur, k, pending = 0, None, 0, []
        while True:
            m = tok.search(s, pos)
            if not m: break
            t = m.group(0)
            if t.startswith("\\part*"):
                title, e = balanced(s, m.end()-1)
                if title.strip() == "Instruments": instr, appendix, secnum = True, False, 0
                pos = e
            elif t.startswith("\\appendix"):
                appendix, secnum, pos = True, 0, m.end()
            elif t.startswith("\\section{"):
                title, e = balanced(s, m.end()-1); secnum += 1; k = 0; pending = []
                cur = f"I.{secnum}" if instr else (chr(64+secnum) if appendix else str(secnum))
                lab = re.match(r"\s*\\label\{mod:([^}]*)\}", s[e:e+400])
                sections[cur] = dict(label=lab.group(1) if lab else None, title=title, file=p.name,
                                     line=s.count("\n", 0, m.start()) + 1)
                pos = e
            elif m.group(1) == "chproof":
                end = s.find("\\end{chproof}", m.end())
                body = s[m.end():end]
                g = dict(section=cur, boxes=list(pending), proof=body, file=p.name,
                         offset=m.end(), line=s.count("\n", 0, m.start()) + 1)
                if pending or not groups or groups[-1]["section"] != cur:
                    groups.append(g)
                else:                                   # continuation of the previous derivation
                    groups[-1].setdefault("parts", []).append(g)
                pending = []; pos = end + 1
            else:
                env = m.group(1); k += 1; e = m.end(); opt = ""
                if s[e:e+1] == "[": opt, e = balanced(s, e, "[", "]")
                end = s.find("\\end{" + env + "}", e)
                num = f"{cur}.{k}"
                boxes[num] = dict(env=env, grade=GRADE[env], title=opt, body=s[e:end], file=p.name,
                                  line=s.count("\n", 0, m.start()) + 1, section=cur)
                pending.append(num); pos = end + 1
    for i, g in enumerate(groups):
        g["id"] = i
    return sections, boxes, groups, files

# ----------------------------------------------------------------------------- M2 premise sentence

MARK = re.compile(r"The premises (?:of this derivation )?are\b|\\textit\{The premises\.\}")
C_END = re.compile(r"\\textit\{The steps|Nothing outside|The chain|\n\s*\n|\\textit\{The ")

def sentence_end(s, i):
    """Index just after the first full stop at paren/brace depth 0 that ends a sentence."""
    d = 0
    for j in range(i, len(s)):
        ch = s[j]
        if ch in "({[": d += 1
        elif ch in ")}]": d = max(0, d-1)
        elif ch == "." and d == 0:
            nxt = s[j+1:j+3]
            if nxt == "" or nxt[:1] == "\n" or (nxt[:1] == " " and (nxt[1:2].isupper() or nxt[1:2] in "\\")):
                return j + 1
    return len(s)

def premise_sentence(g, files):
    """(text, absolute offset in file, file) of the premise sentence of a derivation group."""
    for part in [g] + g.get("parts", []):
        m = MARK.search(part["proof"])
        if not m: continue
        s = files[part["file"]]; a = part["offset"] + m.start()
        if m.group(0).startswith("\\textit"):
            e = C_END.search(s, a + len(m.group(0)))
            b = e.start() if e else a + 2000
        else:
            b = sentence_end(s, a)
        return s[a:b], a, part["file"]
    return None, None, g["file"]

# ----------------------------------------------------------------------------- M3 items

GRADE_PHRASES = [   # (phrase, class); alternatives tried longest first
    ("definitions and a theorem", "DT"),
    ("a theorem and an open result", "TO"),
    ("an open result and two theorems", "TO"),
    ("two open results", "O"), ("open results", "O"), ("an open result", "O"),
    ("conditional theorems", "C"), ("a conditional theorem", "C"),
    ("theorems", "T"), ("a theorem", "T"),
    ("definitions", "D"), ("a definition", "D"), ("the text", "D"),
    ("an instrument", "Tab"), ("among the instruments", "Tab"), ("instruments", "Tab"), ("a table", "Tab"),
]
GALT = "|".join(re.escape(p) for p, _ in GRADE_PHRASES)
GRADE_RE = re.compile(r"\b(" + GALT + r")\b")
GCLASS = dict(GRADE_PHRASES)
PRE_GRADE = re.compile(r"\bthe (theorem|open result|table|conditional theorem) of (?:the chapter titled [^()]*?\()?\s*$")
SELF_PRE = re.compile(r"\bthe (theorem|open result|table|conditional theorem)\s+$")   # "the open result stated above"
POST_GRADE = re.compile(r"\s*(?:,\s*|\()(" + GALT + r")\b")
PLURAL = re.compile(r"(?:,\s*|\()(definitions|theorems|open results)\b")
REF = re.compile(r"(?:chapter|Section|section)~?\\ref\{mod:([^}]+)\}")
DES1 = r"\d+-(?:\d+[a-z]?|L)"
DES = re.compile(r"\bchapters?\s+(" + DES1 + r"(?:\s*(?:,|and|to)\s*(?:chapter\s+)?" + DES1 + r")*)")
SELF = re.compile(r"\b(?:this chapter|the present chapter|chapter above|stated above)\b")
FRAME = re.compile(r"\bthe framework\b")
NOWHERE = re.compile(r"outside the chapters of this compendium|placed at no chapter of this compendium")
ALLOWED = {"T": {"THEOREM", "CONDITIONAL"}, "C": {"CONDITIONAL"}, "O": {"OPEN"},
           "D": {"DEFINITION"}, "Tab": {"TABLE"}, "DT": {"DEFINITION", "THEOREM", "CONDITIONAL"},
           "TO": {"THEOREM", "CONDITIONAL", "OPEN"},
           None: {"THEOREM", "CONDITIONAL", "OPEN", "TABLE", "DEFINITION"}}

def paren_span(s, i):
    """(start, end) of the innermost parenthesis enclosing position i, or None."""
    d, st = 0, None
    for j in range(i - 1, -1, -1):
        if s[j] == ")": d += 1
        elif s[j] == "(":
            if d == 0: st = j; break
            d -= 1
    if st is None: return None
    d = 0
    for j in range(st, len(s)):
        if s[j] == "(": d += 1
        elif s[j] == ")":
            d -= 1
            if d == 0: return (st, j + 1)
    return (st, len(s))

def last_sentence_start(s, a, b):
    """Start of the last top-level sentence (after '. ' or '; ') inside s[a:b]."""
    d, k = 0, a
    for j in range(a, b):
        ch = s[j]
        if ch in "({[": d += 1
        elif ch in ")}]": d = max(0, d - 1)
        elif ch in ".;" and d == 0 and s[j+1:j+2] == " ": k = j + 2
    return k

def items_of(sent):
    """Parse a premise sentence into items (rule M3)."""
    anchors = []
    for m in REF.finditer(sent): anchors.append(dict(kind="REF", target=[m.group(1)], a=m.start(), b=m.end()))
    for m in DES.finditer(sent):
        ds = re.findall(DES1, m.group(1))
        if re.search(r"\bto\b", m.group(1)) and len(ds) == 2: ds = [ds[0] + ".." + ds[1]]
        anchors.append(dict(kind="DES", target=ds, a=m.start(), b=m.end()))
    for m in SELF.finditer(sent): anchors.append(dict(kind="SELF", target=[], a=m.start(), b=m.end()))
    for m in FRAME.finditer(sent): anchors.append(dict(kind="FRAME", target=["FRAMEWORK"], a=m.start(), b=m.end()))
    for m in NOWHERE.finditer(sent): anchors.append(dict(kind="NOWHERE", target=[], a=m.start(), b=m.end()))
    anchors.sort(key=lambda x: (x["a"], -x["b"]))
    out = []
    for x in anchors:                      # drop anchors nested in an earlier one
        if out and x["a"] < out[-1]["b"]: continue
        out.append(x)
    anchors = out
    for i, x in enumerate(anchors):
        x["grade"] = x["gword"] = None
        nxt = anchors[i+1]["a"] if i + 1 < len(anchors) else len(sent)
        x["desc_end"], x["item_end"] = x["a"], x["b"]
        ps = paren_span(sent, x["a"])
        if ps and ps[1] > x["b"]:                                   # (i) inside a parenthesis
            x["desc_end"], x["item_end"] = ps[0], ps[1]
            gm = GRADE_RE.search(sent, ps[0], ps[1])
            if gm: x["grade"], x["gword"] = GCLASS[gm.group(1)], gm.group(1); continue
            nxt = anchors[i+1]["a"] if i + 1 < len(anchors) else len(sent)
            gm = POST_GRADE.match(sent, ps[1], nxt)
            if gm: x["grade"], x["gword"], x["item_end"] = GCLASS[gm.group(1)], gm.group(1), gm.end()
            continue
        pre = PRE_GRADE.search(sent, max(0, x["a"] - 200), x["a"])   # (ii) "the theorem of chapter~\ref"
        if x["kind"] == "SELF" and sent[x["a"]:x["b"]] == "stated above":   # (ii) "the open result stated above"
            sp = SELF_PRE.search(sent, max(0, x["a"] - 40), x["a"])
            if sp:
                w = sp.group(1)
                x["grade"] = {"theorem": "T", "open result": "O", "table": "Tab", "conditional theorem": "C"}[w]
                x["gword"] = "the " + w + " stated above"
                continue
        if pre:
            w = pre.group(1)
            x["grade"] = {"theorem": "T", "open result": "O", "table": "Tab", "conditional theorem": "C"}[w]
            x["gword"] = "the " + w; x["desc_end"] = pre.start()
            if w == "table":
                ps2 = paren_span(sent, x["a"])
                if ps2: x["item_end"] = ps2[1]
            continue
        gm = POST_GRADE.match(sent, x["b"], nxt)                     # (iii) ", a theorem" / " (theorems)"
        if gm:
            x["grade"], x["gword"], x["item_end"] = GCLASS[gm.group(1)], gm.group(1), gm.end()
            if gm.group(0).lstrip().startswith("("):
                ps2 = paren_span(sent, gm.end() - 1)
                if ps2: x["item_end"] = ps2[1]
    for i, x in enumerate(anchors):                                  # (iv) plural grade closing a run
        if x["grade"] is not None or x["kind"] == "NOWHERE": continue
        j, ok = i, None
        while True:
            end = anchors[j+1]["a"] if j + 1 < len(anchors) else len(sent)
            gap_ok = "(" not in sent[anchors[j]["b"]:end] and last_sentence_start(sent, anchors[j]["b"], end) == anchors[j]["b"]
            gm = PLURAL.search(sent, anchors[j]["b"], end)
            if gm and "(" not in sent[anchors[j]["b"]:gm.start()] and \
                    last_sentence_start(sent, anchors[j]["b"], gm.start()) == anchors[j]["b"]:
                ok = (gm.group(1), j, gm.end()); break
            if j + 1 >= len(anchors) or not gap_ok: break
            y = anchors[j+1]
            if y["grade"] is not None:
                if y["gword"] in ("definitions", "theorems", "open results"): ok = (y["gword"], j, None)
                break
            j += 1
        if ok:
            w, j, e = ok
            for y in anchors[i:j+1]:
                y["grade"], y["gword"] = GCLASS[w], w + " (shared)"
            if e: anchors[j]["item_end"] = max(anchors[j]["item_end"], e)
    prev = 0                                                         # descriptions and quotes
    for x in anchors:
        q0 = last_sentence_start(sent, prev, x["desc_end"])
        d = sent[q0:x["desc_end"]]
        d = re.sub(r"^\s*(?:\\textit\{The premises\.\}|The premises (?:of this derivation )?are)", "", d)
        d = re.sub(r"^[\s,;.]*(?:and\s+)?", "", d).strip()
        x["desc"] = d
        q = sent[q0:max(x["item_end"], x["b"])].strip(" ,;")
        x["quote"] = re.sub(r"^(?:and|,)\s+", "", re.sub(r"^\s*(?:\\textit\{The premises\.\}|The premises (?:of this derivation )?are)\s*", "", q))
        prev = max(prev, x["item_end"])
        if x["kind"] == "REF" and x["target"] == ["FRAMEWORK"] and x["grade"] is None and re.search(r"axiom|definition", d, re.I):
            x["grade"], x["gword"] = "D", "(the axiom and the definitions)"
        if x["kind"] == "FRAME" and x["grade"] is None:
            x["grade"], x["gword"] = "D", "(the framework)"
    return anchors

# ----------------------------------------------------------------------------- resolution

def module_sections(sections, des):
    """Reading MOD: sections whose label starts with '<roman P>-c-' (or 'FRAMEWORK' for 1-0)."""
    if des == "1-0": return [n for n, s in sections.items() if s["label"] == "FRAMEWORK"]
    p, c = des.split("-", 1)
    pre = f"{ROMAN.get(int(p), '?')}-{c}-"
    return [n for n, s in sections.items() if s["label"] and (s["label"] + "-").startswith(pre)]

def fifteen_sections(sections, desig, des):
    """Reading FIF: designators.tsv first_section..last_section (else module_numbering range)."""
    d = desig.get(des)
    if not d: return []
    a, b = d["first_section"], d["last_section"]
    if not a: a, b = d["module_numbering_first"], d["module_numbering_last"]
    labs = [s["label"] for s in sections.values()]
    nums = list(sections.keys())
    if a not in labs or b not in labs: return []
    i, j = labs.index(a), labs.index(b)
    return nums[min(i, j):max(i, j)+1]

def chapter_prefix(label):
    """Module chapter of a section label: 'VI-10-higgs-...' -> 'VI-10'; 'FRAMEWORK' -> 'FRAMEWORK'."""
    if not label: return None
    if label == "FRAMEWORK": return label
    m = re.match(r"([IVX]+-(?:\d+[a-z]?|L|C))-", label + "-")
    return m.group(1) if m else None

def resolve(item, g, sections, boxes, desig, reading):
    """Sections an item points to under a reading ('MOD', 'FIF' or 'UNION')."""
    k = item["kind"]
    if k == "REF":
        return [n for n, s in sections.items() if s["label"] == item["target"][0]]
    if k == "FRAME":
        return [n for n, s in sections.items() if s["label"] == "FRAMEWORK"]
    if k == "NOWHERE":
        return []
    if k == "SELF":
        own = chapter_prefix(sections[g["section"]]["label"])
        return [n for n, s in sections.items() if chapter_prefix(s["label"]) == own]
    secs = []
    for d in item["target"]:
        ds = [d]
        if ".." in d:                      # "chapters a to b": every designator between them
            a, b = d.split("..")
            keys = list(desig.keys())
            ds = keys[keys.index(a):keys.index(b)+1] if a in keys and b in keys else [a, b]
        for x in ds:
            if reading in ("MOD", "UNION"): secs += module_sections(sections, x)
            if reading in ("FIF", "UNION"): secs += fifteen_sections(sections, desig, x)
    return list(dict.fromkeys(secs))

def candidates(item, g, sections, boxes, desig, reading):
    """Rule M4: statements of the resolved sections with a grade compatible with the grade word."""
    secs = set(resolve(item, g, sections, boxes, desig, reading))
    allowed = ALLOWED[item["grade"]]
    c = [n for n, b in boxes.items() if b["section"] in secs and b["grade"] in allowed]
    stated_above = item["kind"] == "SELF" and re.search(r"\bstated above\b", item["quote"])
    if stated_above:                       # the boxes standing above the derivation, the group's own included
        order = list(boxes.keys()); last = order.index(g["boxes"][-1])
        return [n for n in c if order.index(n) <= last]
    if item["kind"] == "SELF" and "above" in item["quote"]:
        order = list(boxes.keys()); first = order.index(g["boxes"][0])
        c = [n for n in c if order.index(n) < first]
    if item["grade"] == "D" and not c and (re.search(r"axiom|definition", item["desc"], re.I) or item["gword"] in ("the text", "(the framework)")):
        c = [n for n, b in boxes.items() if b["grade"] == "DEFINITION" and b["section"] == "1"]
    if item["kind"] == "SELF":
        c = [n for n in c if n not in g["boxes"]] or c
    return c

# ----------------------------------------------------------------------------- graphs

def build(sections, boxes, groups, files, desig, reading):
    """Rule M5. Returns (items per group, edges {(p, s): [item refs]}, box -> group)."""
    gitems, edges, gof = {}, collections.defaultdict(list), {}
    for g in groups:
        for b in g["boxes"]: gof[b] = g["id"]
        sent, off, f = premise_sentence(g, files)
        g["sentence"], g["s_file"] = sent, f
        g["s_line"] = files[f].count("\n", 0, off) + 1 if off is not None else None
        its = items_of(sent) if sent else []
        for i, x in enumerate(its):
            x["cand"] = candidates(x, g, sections, boxes, desig, reading)
            x["id"] = f"g{g['id']}.{i}"
            for c in x["cand"]:
                for b in g["boxes"]:
                    if c != b: edges[(c, b)].append(x["id"])
        gitems[g["id"]] = its
    return gitems, edges, gof

def sccs(nodes, edges):
    ch = collections.defaultdict(set)
    for (p, s) in edges: ch[p].add(s)
    idx, low, onst, st, out, c = {}, {}, set(), [], [], [0]
    sys.setrecursionlimit(20000)
    def strong(v):
        idx[v] = low[v] = c[0]; c[0] += 1; st.append(v); onst.add(v)
        for w in ch[v]:
            if w not in idx: strong(w); low[v] = min(low[v], low[w])
            elif w in onst: low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = st.pop(); onst.discard(w); comp.append(w)
                if w == v: break
            if len(comp) > 1: out.append(comp)
    for v in nodes:
        if v not in idx: strong(v)
    return out

def shortest_cycle(x, edges):
    """Shortest premise cycle through x (list of statements, x first), or None."""
    ch = collections.defaultdict(set)
    for (p, s) in edges: ch[p].add(s)
    prev, frontier, seen = {}, [x], {x}
    while frontier:
        nxt = []
        for v in frontier:
            for w in sorted(ch[v]):
                if w == x:
                    path = [v]
                    while path[-1] != x: path.append(prev[path[-1]])
                    return [x] + path[::-1][1:]
                if w not in seen: seen.add(w); prev[w] = v; nxt.append(w)
        frontier = nxt
    return None

# ----------------------------------------------------------------------------- layer J (judgment)
# Key: first box of the derivation group -> list of (start of the item's quote, choice, evidence).
# choice: list of statements (must be a subset of layer M's candidates, MOD reading), or
#         []   no boxed statement is meant (stated in the note), or
#         None the description names no box: layer M's candidate set is kept.
# evidence: (box, verbatim words from that box's title or text) supporting the choice, or a note.
# Layer J only narrows M; it never adds a premise the derivation does not cite.
DEFS = ["1.1", "1.2", "1.3", "1.4", "1.5", "1.6"]
JUDGMENT = {
    "2.1": [
        ("the axiom and its six definitions", DEFS, "the text = the axiom and Definitions 1.1-1.6"),
        ("the one-act clause as a theorem", ["2.3"], ("2.3", "The one-act clause is a theorem of the derived act")),
        ("the listed webs with their end lines", None, "no THEOREM/CONDITIONAL box of chapter 1-5 names the listed webs; the one box of 1-5 whose title names the list is H.2 (OPEN)"),
    ],
    "9.1": [
        ("the axiom with its definition of a read", DEFS, "the text"),
        ("the carrier as the closed read of its type", ["66.4"], ("66.4", "The carriers as the closed reads of their type faces")),
        ("the read engine", [], "an instrument outside the chapters (no box)"),
    ],
    "18.1": [
        ("the axiom with its six definitions", DEFS, "the text"),
        ("the link with its six types and the readability order", ["9.1"], ("9.1", "Readability order at the six genres")),
        ("the explicit two-point function with its records", ["124.1"], ("124.1", "the records of its two reads are 172")),
        ("the sector-to-ends map at the six types", ["I.1"], ("I.1", "The sector-to-ends map at the parting web")),
    ],
    "124.1": [
        ("the axiom and the six definitions", DEFS, "Definitions 1.1-1.6"),
        ("the open result of chapter~\\ref{mod:I-1", None, "section 2 holds four OPEN boxes (2.4, 2.5, 2.7, 2.8); the derivation does not say which"),
        ("the theorem of chapter~\\ref{mod:I-2", ["9.1", "9.2", "9.4"], ("124.1", "the $(2,3)$ link alone carries a central generator acting on its pair, abelian, conserved and massless")),
        ("the theorem of chapter~\\ref{mod:I-3", ["18.2", "18.3"], ("124.1", "with the weights one, root two, root three and root six by the gain law")),
        ("the theorem of chapter~\\ref{mod:I-6", None, "any of 28.1, 28.2, 28.5 (all THEOREM); the derivation does not say which"),
        ("the theorem of chapter~\\ref{mod:II-1", ["35.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:II-2", ["40.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:III-1", ["50.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:III-2", ["53.1"], "unique"),
        ("the open result of chapter~\\ref{mod:III-5", ["57.2"], "unique"),
    ],
    "16.1": [
        ("the readability order at the six genres", ["9.1"], ("9.1", "Readability order at the six genres")),
        ("the count law at the six genres", None, "no THEOREM/CONDITIONAL box of chapter 8-1 names the count law (193.1, OPEN, mentions it)"),
        ("the Higgs and carrier offices with the body-alone rule", ["200.6"], ("200.6", "The Higgs and carrier offices at their sinks")),
        ("the stabiliser of the pair ket", ["131.1"], ("131.1", "The stabilizer of the pair ket and the unbroken charge")),
    ],
    "58.1": [
        ("the stabiliser of the pair ket with its identification algebra", ["131.1"], ("131.1", "The stabilizer of the pair ket and the unbroken charge")),
        ("the Lorentz group at the qubit pair", ["54.1"], ("54.1", "The Lorentz group as the isometries of the Hermitian split")),
        ("the boost open result", [], "declared OPEN; chapter 3-3 holds no OPEN box (56.2 'The boost at the one unit' is in chapter 3-4)"),
        ("the dock open result and the dock-register open result", None, "declared OPEN; 57.2 is the only OPEN box of chapter 3-5 and its title does not name the dock"),
        ("the photon row at the centre", None, "declared OPEN; no OPEN box of chapter 6-7 names a photon row"),
    ],
    "200.1": [
        ("the link and its six genres of chapter 1-2", [], "chapters 1-2, 1-3, 1-4 cited as 'definitions'; they hold no DEFINITION box"),
        ("the reads and the Born count of chapter 1-3", [], "as above"),
        ("the moves of chapter 1-4", [], "as above"),
        ("the read at a register of chapter 1-6", None, "chapter-level citation without a grade word"),
        ("the gauge chapters 4-2 and 4-3", None, "chapter-level citation without a grade word"),
        ("the doublets of chapter 6-3", None, "chapter-level citation without a grade word"),
        ("the mass unit of chapter 8-1", ["185.1"], ("185.1", "The mass unit and the first power")),
    ],
    "200.5": [
        ("the Higgs chapter above", None, "the boxes above in the chapter: 200.1-200.4"),
        ("the framework", DEFS, "the framework chapter"),
    ],
    "200.6": [
        ("the Higgs chapter above", None, "the boxes above in the chapter: 200.1-200.5"),
        ("the body-alone office rule of chapter 8-2", ["195.1"], ("195.1", "The body-alone office rule")),
    ],
    "8.1": [
        ("the three readout rules at the three levels of the wash", ["24.1"], ("24.1", "The three readout rules as one read at three levels")),
        ("the wash move at the cap", [], "chapter 11-6 holds no THEOREM/CONDITIONAL box"),
        ("the body-alone office of every word by the one rule", None, "no box title of chapter 6-10 names it"),
        ("the body-alone square at the five bodies", ["195.1"], ("195.1", "The body-alone office rule")),
        ("the two definitions of open at every end", None, "declared OPEN; any of 2.4, 2.5, 2.7, 2.8, 4.1"),
    ],
    "27.1": [
        ("the confinement as the readability of the dark genre", ["171.1"], ("171.1", "Confinement as the dark-genre readability")),
        ("the gluon pair ket through the scale body", ["149.1"], ("149.1", "The strong coupling as the gluon pair-ket office through the scale body")),
        ("the quark mixing block", ["166.6"], ("166.6", "The quark mixing block at the read sites")),
        ("the gain lemma", ["18.2"], ("18.2", "The gain lemma at the spectrum of records")),
        ("the conjugate-letter sorting", None, "no box title of chapter 7-7 names it"),
        ("the body-alone office of every body", ["195.1"], ("195.1", "The body-alone office rule")),
    ],
    "39.1": [
        ("the Higgs and carrier offices with the body-alone rule", ["200.6"], ("200.6", "The Higgs and carrier offices at their sinks")),
        ("the unpassed in-legs at the five bodies", ["191.1"], ("191.1", "The dressing as the unpassed (3,4) in-legs")),
        ("the Higgs leg unit", None, "declared OPEN; 200.4 or 200.5"),
        ("the readability order at the six genres", ["9.1"], ("9.1", "Readability order at the six genres")),
        ("the forces as the move classes", ["E.1"], ("E.1", "The forces as the move classes at the genres")),
        ("the twelve-slot table", None, "instrument tables I.9.1, I.10.1"),
    ],
    "69.1": [
        ("the terminated saturation at the census", ["H.3", "U.1"], ("H.3", "The termination of seeded saturation")),
        ("the body-alone office rule at the measure-struck reading", ["195.1"], ("195.1", "The body-alone office rule")),
        ("the Higgs and carrier offices with the word faces", ["200.6"], ("200.6", "The Higgs and carrier offices at their sinks")),
        ("the Higgs attachment at the (4,4) count", ["200.5"], ("200.5", "The Higgs attachment at the saturated web")),
        ("the count law at the six genres", None, "as for 16.1"),
    ],
    "70.1": [
        ("the Green-function row", None, "declared OPEN; no box title of chapter 10-2 names it"),
        ("the two-two two-point function at its register", None, "declared OPEN; no box title of chapter 10-2 names it"),
        ("the Ward-functional row", [], "declared OPEN; chapter 4-8 holds no OPEN box"),
        ("the serial insertion row at the chains", None, "declared OPEN; 116.1, 117.1 or 119.1"),
    ],
    "74.1": [
        ("the moves in one table", None, "D.3 or D.4 (both THEOREM)"),
        ("the content exponents of every web", None, "no box of chapter 1-5 names content exponents"),
        ("the move matrix of the parent result", None, "D.3 or D.4 (both THEOREM)"),
    ],
    "78.1": [
        ("the one-act clause as a theorem", ["2.3"], ("2.3", "The one-act clause is a theorem of the derived act")),
        ("conservation as one theorem of one link", ["74.3"], ("74.3", "conservation is one theorem of one link")),
        ("the explicit two-point function with its one hundred seventy-two records", ["124.1"], ("124.1", "the records of its two reads are 172")),
        ("the dressing at the listed face", None, "no box title of chapter 1-3 names it"),
    ],
    "146.1": [
        ("the axiom and the six definitions", DEFS, "Definitions 1.1-1.6"),
        ("the theorem of chapter~\\ref{mod:I-3", None, "any of 18.1-18.4"),
        ("the theorem of chapter~\\ref{mod:I-6", None, "any of 28.1, 28.2, 28.5"),
        ("the theorem of chapter~\\ref{mod:II-6", ["44.1"], "unique"),
        ("the open result of chapter~\\ref{mod:II-7", ["45.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:V-1", ["124.1"], "unique"),
    ],
    "147.1": [
        ("the axiom and the six definitions", DEFS, "Definitions 1.1-1.6"),
        ("the theorem of chapter~\\ref{mod:I-2", None, "any of 9.1-9.7"),
        ("the theorem of chapter~\\ref{mod:I-3", None, "any of 18.1-18.4"),
        ("the theorem of chapter~\\ref{mod:I-6", None, "any of 28.1, 28.2, 28.5"),
        ("the theorem of chapter~\\ref{mod:III-1", ["50.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:III-7", ["62.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:V-1", ["124.1"], "unique"),
        ("the theorem of chapter~\\ref{mod:V-6", ["146.1"], "unique"),
    ],
    # groups on the derivation-level paths between 124.1 and its premises
    "64.1": [
        ("the axiom with the definition", DEFS, "the text"),
        ("the admissible matchings as the whole unitary group", None, "59.1 or 61.1"),
        ("the identification algebra at the (2,3) tether", [], "a theorem placed at no chapter (no box)"),
        ("the explicit two-point function with its records", ["124.1"], ("124.1", "the records of its two reads are 172")),
        ("the Born count as the readout of a record", ["40.1"], "unique"),
    ],
    "66.1": [
        ("the listed webs with their end lines and the four bodies", None, "as for 2.1-2.9"),
        ("the Cartan weights at the open ends", ["72.1"], ("72.1", "The Cartan weights and the charge")),
        ("the gauge group with its sixteen rows as a comparison target and no input", [], "the derivation states 'as a comparison target and no input'"),
        ("the dart charge at the dock face", [], "declared OPEN, placed at no chapter (no box)"),
    ],
    "131.1": [
        ("the axiom and the six definitions", DEFS, "Definitions 1.1-1.6"),
        ("the link and its six genres", None, "any of 9.1-9.7"),
        ("the records, reads and the Born count", None, "any of 18.1-18.4"),
        ("the gauge group of chapter", ["64.2"], ("64.2", "The gauge group as the face group of the tether")),
        ("the fine-structure constant at the one link", ["124.1"], ("124.1", "The fine-structure constant at the one link")),
    ],
}
# What the shared derivation of a group concludes (verbatim; judgment: the procedure derives that box).
GROUP_NOTES = {
    "2.1": ("2.7", "so the hub step is a third letter, which is the open result"),
    "9.1": ("9.1", "so the readability order is a typing of the type, which is the theorem"),
    "18.1": ("18.1", "So the Born count is a measure and the readout its root, which is the theorem"),
}
REQUESTED = (["2.2", "2.3", "2.4", "2.5", "2.6", "2.7", "2.8", "9.1", "9.7", "9.8"] +
             [f"18.{k}" for k in range(1, 7)] + ["124.1", "16.1", "58.1"] + [f"200.{k}" for k in range(1, 7)])

def norm(t):
    return re.sub(r"\s+", " ", t).strip()

def apply_judgment(groups, gitems, boxes, gof):
    """Attach item['J'] and item['Jnote'] for judged groups; return list of problems (must be empty)."""
    probs = []
    for key, rules in JUDGMENT.items():
        g = groups[gof[key]]
        its = gitems[g["id"]]
        if g["boxes"][0] != key: probs.append(f"{key}: not the first box of its group")
        for pre, choice, ev in rules:
            hit = [x for x in its if x["quote"].startswith(pre)]
            if len(hit) != 1:
                probs.append(f"{key}: rule '{pre}' matches {len(hit)} items"); continue
            x = hit[0]
            if choice is not None and not set(choice) <= set(x["cand"]):
                probs.append(f"{key}: rule '{pre}' chooses {choice} outside M {x['cand']}")
            if isinstance(ev, tuple):
                b = boxes[ev[0]]
                if norm(ev[1]) not in norm(b["title"] + " " + b["body"]):
                    probs.append(f"{key}: evidence for '{pre}' not verbatim in {ev[0]}")
                note = f"{ev[0]}: \"{ev[1]}\""
            else:
                note = ev
            x["J"], x["Jnote"] = choice, note
        for x in its:
            x.setdefault("J", "unjudged")
    return probs

def layer_cands(x, layer):
    """Candidate set of an item in a layer: 'M' (mechanical) or 'J' (judgment where given, else M)."""
    if layer == "M" or "J" not in x or x["J"] in (None, "unjudged"): return x["cand"]
    return x["J"]

def layer_edges(groups, gitems, layer):
    E = collections.defaultdict(list)
    for g in groups:
        for x in gitems[g["id"]]:
            for c in layer_cands(x, layer):
                for b in g["boxes"]:
                    if c != b: E[(c, b)].append(x["id"])
    return E

# ----------------------------------------------------------------------------- chain rule

UNPROVEN = {"CONDITIONAL", "OPEN"}

def item_status(x, b, boxes, layer):
    """'open-declared', 'unproven', 'mixed', 'proven' or 'none' for item x as a premise of box b."""
    c = [n for n in layer_cands(x, layer) if n != b]
    gs = {boxes[n]["grade"] for n in c}
    if x["grade"] in ("O", "TO"): return "open-declared"         # "open result", or a conjunction naming one
    if x["grade"] == "C" and not c: return "unproven"
    if not c: return "none"
    if gs <= UNPROVEN: return "unproven"
    if gs & UNPROVEN: return "mixed"
    return "proven"

def verdict(b, its, boxes, layer):
    st = [item_status(x, b, boxes, layer) for x in its]
    if "open-declared" in st: return "violates (declared OPEN premise)", st
    if "unproven" in st: return "violates (CONDITIONAL/OPEN premise)", st
    if "mixed" in st: return "undetermined", st
    return "consistent", st

# ----------------------------------------------------------------------------- index graph (as depgraph.py, C3)

def index_graph(rows, sections, desig, module=()):
    """Index premise graph as depgraph.py (C3); designators in `module` read by module numbering."""
    by_num = {r["number"]: r for r in rows}
    by_section = collections.defaultdict(list)
    for r in rows: by_section[r["number"].rsplit(".", 1)[0]].append(r)
    dsec = collections.OrderedDict((d, module_sections(sections, d) if (module == "all" or d in module)
                                    else fifteen_sections(sections, desig, d)) for d in desig)
    prem = {}
    for r in rows:
        p = parse_refs(r["premises"], by_section, dsec)[0]
        prem[r["number"]] = {x for x in p if x in by_num and x != r["number"]}
    return {(p, s) for s, ps in prem.items() for p in ps}, prem

# ----------------------------------------------------------------------------- reachability

def closure(x, E, down=True):
    nb = collections.defaultdict(set)
    for (p, s) in E:
        if down: nb[p].add(s)
        else: nb[s].add(p)
    seen, st = set(), [x]
    while st:
        y = st.pop()
        for z in nb[y]:
            if z not in seen: seen.add(z); st.append(z)
    seen.discard(x)
    return seen

def definite_dependents(P, groups, gitems, layer):
    """Lower bound: boxes that rest on some member of P whatever the any-of items resolve to.
    Fixpoint: a box is added when one item of its derivation has a non-empty candidate set
    contained in the current set."""
    D = set(P); changed = True
    while changed:
        changed = False
        for g in groups:
            for x in gitems[g["id"]]:
                c = set(layer_cands(x, layer))
                if c and c <= D:
                    for b in g["boxes"]:
                        if b not in D: D.add(b); changed = True
    return D - set(P)

def shortest_path(a, b, E):
    ch = collections.defaultdict(set)
    for (p, s) in E: ch[p].add(s)
    prev, fr = {a: None}, [a]
    while fr:
        nx = []
        for v in fr:
            for w in sorted(ch[v]):
                if w in prev: continue
                prev[w] = v; nx.append(w)
                if w == b:
                    out = [b]
                    while prev[out[-1]] is not None: out.append(prev[out[-1]])
                    return out[::-1]
        fr = nx
    return None

def forced_cycles(groups, gitems, gof, boxes, layer):
    """Boxes X that rest on X under every resolution of the any-of items (X definitely depends on X).
    A SELF item never names the box itself (rule M3), so X is dropped from its candidates
    ("the open result stated above" is 217.3 for 217.1 and 217.2, not a premise of 217.3)."""
    out = []
    for X in boxes:
        if X not in gof: continue
        D = {X} | definite_dependents([X], groups, gitems, layer)
        for x in gitems[gof[X]]:
            c = set(layer_cands(x, layer)) - ({X} if x["kind"] == "SELF" else set())
            if c and c <= D:
                out.append(X); break
    return out

def overlap(desc, titles):
    stop = {"with", "their", "every", "from", "that", "this", "which", "the", "and", "chapter", "chapters", "theorem", "result"}
    w = lambda t: {x for x in re.findall(r"[a-z]{4,}", t.lower()) if x not in stop}
    d = w(desc)
    return max([len(d & w(t)) for t in titles] + [0])

# ----------------------------------------------------------------------------- report

def fmt(c, boxes, n=12):
    s = ", ".join(f"{x} {G1[boxes[x]['grade']]}" for x in c[:n])
    return s + (f", ... (+{len(c)-n})" if len(c) > n else "") if c else "-"

def main(cdir, tsv=None):
    cdir = pathlib.Path(cdir)
    rows = load_tsv(cdir/"index"/"statements.tsv")
    by_num = {r["number"]: r for r in rows}
    desig = collections.OrderedDict((d["number"], d) for d in load_tsv(cdir/"index"/"designators.tsv"))
    sections, boxes, groups, files = segment(cdir/"tex")
    P = print
    P("Statement-level premises of the Chain Theory compendium 0.1 (round 4)")
    P("=" * 72)
    mism = [r["number"] for r in rows if r["number"] not in boxes or boxes[r["number"]]["grade"] != r["grade"]]
    P(f"segmentation: {len(sections)} sections, {len(boxes)} boxes ({len(rows)} index rows; numbering/grade mismatches: {len(mism)})")
    P(f"derivation groups: {len(groups)} (box counts {dict(sorted(collections.Counter(len(g['boxes']) for g in groups).items()))}); "
      f"continuation parts merged: {sum(len(g.get('parts', [])) for g in groups)}")
    nod = [b for b in boxes if not any(b in g["boxes"] for g in groups)]
    P(f"boxes with no derivation in their section: {len(nod)} {nod}")

    # --- the source's own rules, verbatim
    P("\n[0] The rule and the premise line, verbatim")
    for f, needle in [("front.tex", "under the chain rule of this compendium every theorem prints at its statement the union of the conditions of the results it rests on, a composed claim carrying the lowest status in its chain"),
                      ("front.tex", "A theorem is a statement derived from the axiom, the definitions and the results proved ahead of it"),
                      ("01-framework.tex", "one line names the chapters it rests on and the chapters that rest on it, generated mechanically from the ancestry of the chapters")]:
        s = files[f]; i = s.find(needle)
        P(f"  {f}:{s.count(chr(10), 0, i) + 1 if i >= 0 else 'NOT FOUND'}  \"{needle}\"")

    # --- layer M over every derivation, readings MOD / FIF / UNION
    layers = {}
    for rd in ("MOD", "FIF", "UNION"):
        layers[rd] = build(sections, boxes, groups, files, desig, rd)
    gitems, EM, gof = layers["MOD"]
    nosent = [g for g in groups if g["sentence"] is None]
    allit = [x for g in groups for x in gitems[g["id"]]]
    P("\n[1] Layer M coverage (all derivations)")
    P(f"  groups with a premise sentence: {len(groups) - len(nosent)} of {len(groups)}; without: {len(nosent)} "
      f"(boxes {sum(len(g['boxes']) for g in nosent)}, grades {dict(sorted(collections.Counter(boxes[b]['grade'] for g in nosent for b in g['boxes']).items()))})")
    P(f"  items: {len(allit)}; by kind {dict(sorted(collections.Counter(x['kind'] for x in allit).items()))}")
    P(f"  by grade word {dict(sorted(collections.Counter(str(x['grade']) for x in allit).items()))}")
    nc = collections.Counter(min(len(x['cand']), 2) for x in allit)
    P(f"  candidates per item (MOD): none {nc[0]}, one {nc[1]}, several {nc[2]}")
    # quotes verbatim
    bad = [x["id"] for g in groups for x in gitems[g["id"]] if x["quote"] not in files[g["s_file"]]]
    P(f"  item quotes found verbatim in the source: {len(allit) - len(bad)} of {len(allit)}")
    numcit = sum(len(re.findall(r"(?:Theorems?|Open results?|Tables?|Definitions?)~?\s*(?:\d+|[A-Z]|I\.\d+)\.\d+", part["proof"]))
                 for g in groups for part in [g] + g.get("parts", []))
    boxref = sum(part["proof"].count("\\ref{box:") for g in groups for part in [g] + g.get("parts", []))
    P(f"  statement-number citations inside derivations ('Theorem 9.1', 'Open result 2.4', ...): {numcit}; \\ref{{box:...}} inside derivations: {boxref}")
    # reading evidence
    des = [(g, x) for g in groups for x in gitems[g["id"]] if x["kind"] == "DES"]
    rs = {rd: [candidates(x, g, sections, boxes, desig, rd) for g, x in des] for rd in ("MOD", "FIF")}
    nz = {rd: sum(1 for c in rs[rd] if c) for rd in rs}
    differ = [i for i in range(len(des)) if set(rs["MOD"][i]) != set(rs["FIF"][i])]
    win = collections.Counter()
    for i in differ:
        g, x = des[i]
        sm = overlap(x["desc"], [boxes[c]["title"] for c in rs["MOD"][i]] + [sections[boxes[c]["section"]]["title"] for c in rs["MOD"][i]])
        sf = overlap(x["desc"], [boxes[c]["title"] for c in rs["FIF"][i]] + [sections[boxes[c]["section"]]["title"] for c in rs["FIF"][i]])
        win["MOD" if sm > sf else "FIF" if sf > sm else "tie"] += 1
    P(f"  designator items (chapter P-c): {len(des)}; with a grade-compatible box under MOD {nz['MOD']}, under FIF {nz['FIF']}; "
      f"readings differ at {len(differ)}")
    P(f"  where they differ, description-title word overlap is larger under MOD at {win['MOD']}, under FIF at {win['FIF']}, tie {win['tie']}")

    # --- index graph
    EI, prem_idx = index_graph(rows, sections, desig)
    SI = sccs(list(by_num), EI)
    P(f"\n[2] Index graph (premises column, C3 as depgraph.py): {len(EI)} edges; SCCs > 1: {sorted(len(c) for c in SI)}")
    for nm, mod in (("4-5 read by module numbering", ("4-5",)), ("every chapter-form designator read by module numbering", "all")):
        E2, _ = index_graph(rows, sections, desig, mod)
        S2 = sccs(list(by_num), E2)
        two = sorted({tuple(sorted(e)) for e in E2 if (e[1], e[0]) in E2 and "124.1" in e})
        P(f"  sensitivity, {nm}: {len(E2)} edges; SCCs > 1: {sorted(len(c) for c in S2)}; two-cycles through 124.1: {two}")
    d45 = desig["4-5"]
    P(f"  designators.tsv row 4-5: first_section {d45['first_section']}, mapped_by {d45['mapped_by']}, "
      f"module_numbering_first {d45['module_numbering_first']} (two readings)")

    # --- judgment
    probs = apply_judgment(groups, gitems, boxes, gof)
    P(f"\n[3] Layer J: {len(JUDGMENT)} derivation groups judged; consistency problems: {len(probs)} {probs}")
    judged_items = [x for k in JUDGMENT for x in gitems[gof[k]]]
    jc = collections.Counter("narrowed" if isinstance(x["J"], list) and x["J"] and set(x["J"]) != set(x["cand"]) else
                             "same as M" if isinstance(x["J"], list) and set(x["J"]) == set(x["cand"]) else
                             "no box" if x["J"] == [] else "unresolved (M kept)" for x in judged_items)
    P(f"  items in judged groups: {len(judged_items)}; {dict(sorted(jc.items()))}")
    EJ = layer_edges(groups, gitems, "J")
    R = [b for k in JUDGMENT for b in groups[gof[k]]["boxes"]]
    scc22 = sorted(max(SI, key=len)) if SI else []
    REQ = list(dict.fromkeys(REQUESTED + scc22))
    P(f"  reconstructed statements R (all boxes of judged groups): {len(R)}; requested (targets + index SCC): {len(REQ)}; "
      f"requested not in R: {[b for b in REQ if b not in R]}")

    # --- per-statement table
    P("\n[4] Statement-level premises of the reconstructed statements")
    for k in JUDGMENT:
        g = groups[gof[k]]
        bl = ", ".join(f"{b} {G1[boxes[b]['grade']]}" for b in g["boxes"])
        P(f"\n  derivation of {bl}"
          f"   ({g['s_file']}:{g['s_line']})")
        for x in gitems[g["id"]]:
            j = x["J"]
            js = ("= M" if isinstance(j, list) and set(j) == set(x["cand"]) else fmt(j, boxes) if j else
                  "no box" if j == [] else "unresolved, M kept")
            P(f"    {x['id']:>8} [{x['kind']} {','.join(x['target'])[:28]}; {x['gword']}] \"{x['quote'][:150]}{'...' if len(x['quote']) > 150 else ''}\"")
            P(f"             M: {fmt(x['cand'], boxes)}  |  J: {js}  ({x['Jnote']})")

    # --- regrade
    P("\n[5] Re-grade under the chain rule (front.tex:34): a THEOREM rests only on the axiom, definitions and THEOREMs")
    for k, (bx, q) in GROUP_NOTES.items():
        g = groups[gof[k]]
        i = files[g["file"]].find(q)
        P(f"  shared derivation of {g['boxes'][0]}-{g['boxes'][-1]} (begins {g['file']}:{g['line']}) concludes at "
          f"{g['file']}:{files[g['file']].count(chr(10), 0, i) + 1 if i >= 0 else '?'}: \"{q}\" -> {bx}; verbatim: {norm(q) in norm(g['proof'])}")
    gitF = layers["FIF"][0]
    P("  | box | grade | index: O / C premises | M verdict (MOD) | M verdict (FIF) | J verdict | unproven premises (J) |")
    P("  |---|---|---|---|---|---|---|")
    vt = collections.Counter(); out_rows = []
    for b in R:
        g = groups[gof[b]]; its = gitems[g["id"]]
        io = sum(1 for p in prem_idx[b] if boxes[p]["grade"] == "OPEN"); ic = sum(1 for p in prem_idx[b] if boxes[p]["grade"] == "CONDITIONAL")
        vm, _ = verdict(b, its, boxes, "M"); vj, sj = verdict(b, its, boxes, "J")
        vf, _ = verdict(b, gitF[g["id"]], boxes, "M")
        out_rows.append((b, vm, vj))
        if boxes[b]["grade"] == "THEOREM":
            up = []
            for x, st in zip(its, sj):
                if st in ("open-declared", "unproven"):
                    c = [n for n in layer_cands(x, "J") if n != b]
                    up.append(f"{x['id']}: " + (fmt(c, boxes, 4) if c else "unboxed") + (" [declared O]" if st == "open-declared" else ""))
            P(f"  | {b} | T | {io} / {ic} | {vm} | {vf} | {vj} | {'; '.join(up) or '-'} |")
            vt[("R", vj.split(" (")[0], "declared" in vj)] += 1
            if b in REQ: vt[("REQ", vj.split(" (")[0], "declared" in vj)] += 1
    thR = [b for b in R if boxes[b]["grade"] == "THEOREM"]; thQ = [b for b in thR if b in REQ]
    for nm, S in (("requested", thQ), ("R", thR)):
        io = sum(1 for b in S if any(boxes[p]["grade"] == "OPEN" for p in prem_idx[b]))
        ioc = sum(1 for b in S if any(boxes[p]["grade"] in UNPROVEN for p in prem_idx[b]))
        jv = [verdict(b, gitems[gof[b]], boxes, "J")[0] for b in S]
        mv = [verdict(b, gitems[gof[b]], boxes, "M")[0] for b in S]
        cnt = lambda v, s: sum(1 for z in v if z.startswith(s))
        P(f"  THEOREMs in {nm}: {len(S)} {S}")
        P(f"    index (F3 basis): list an OPEN premise {io}; OPEN or CONDITIONAL {ioc}")
        P(f"    layer M: declared OPEN {cnt(mv, 'violates (declared')}; violates {cnt(mv, 'violates')}; undetermined {cnt(mv, 'undet')}; consistent {cnt(mv, 'consistent')}")
        P(f"    layer J: declared OPEN {cnt(jv, 'violates (declared')}; violates {cnt(jv, 'violates')}; undetermined {cnt(jv, 'undet')}; consistent {cnt(jv, 'consistent')}")
    # CONDITIONALs resting on OPEN (literal 'lowest status')
    coR = [b for b in R if boxes[b]["grade"] == "CONDITIONAL"]
    co_open = [b for b in coR if any(item_status(x, b, boxes, "J") == "open-declared" or
                                     (layer_cands(x, "J") and all(boxes[c]["grade"] == "OPEN" for c in layer_cands(x, "J") if c != b) and [c for c in layer_cands(x, "J") if c != b])
                                     for x in gitems[gof[b]])]
    P(f"  CONDITIONALs in R: {len(coR)}; whose derivation rests on an OPEN result (J): {len(co_open)} {co_open}; index: "
      f"{sum(1 for b in coR if any(boxes[p]['grade'] == 'OPEN' for p in prem_idx[b]))}")

    # --- cycles
    P("\n[6] Cycles")
    big = max(SI, key=len)
    P(f"  index SCC containing 124.1: {len(big)} statements")
    two = sorted({tuple(sorted(e)) for e in EI if (e[1], e[0]) in EI and "124.1" in e})
    for a, c in two:
        for (p, s) in ((a, c), (c, a)):
            g = groups[gof[s]]
            its = [x for x in gitems[g["id"]] if p in x["cand"]]
            jits = [x for x in gitems[g["id"]] if p in layer_cands(x, "J")]
            P(f"    index edge {p} -> {s}: derivation of {s} cites {p}? M: {'yes ' + str([x['id'] for x in its]) if its else 'no'}; J: {'yes' if jits else 'no'}")
    for rd in ("MOD", "FIF", "UNION"):
        E = layers[rd][1]
        S = sccs(list(boxes), E)
        s124 = [c for c in S if "124.1" in c]
        P(f"  layer M ({rd}): {len(E)} edges; SCCs > 1: {sorted((len(c) for c in S), reverse=True)}; SCC of 124.1: {len(s124[0]) if s124 else 0}; "
          f"shortest cycle through 124.1: {shortest_cycle('124.1', E)}")
    S = sccs(list(boxes), EJ)
    s124 = [c for c in S if "124.1" in c]
    P(f"  layer J (J where judged, M elsewhere): {len(EJ)} edges; SCC of 124.1: {len(s124[0]) if s124 else 0}")
    for nm, E in (("M MOD", layers["MOD"][1]), ("M FIF", layers["FIF"][1]), ("M UNION", layers["UNION"][1]), ("J", EJ)):
        anc = closure("124.1", E, down=False)
        P(f"  transitive premises of 124.1 ({nm}): {len(anc)}; 16.1 among them: {'16.1' in anc}; 58.1 among them: {'58.1' in anc}")
    pth = shortest_path("16.1", "124.1", EJ)
    P(f"  shortest J path 16.1 => 124.1: {pth}")
    for u, v in zip(pth or [], (pth or [])[1:]):
        x = next(x for x in gitems[gof[v]] if u in layer_cands(x, "J"))
        P(f"     {u} -> {v}: item {x['id']} (any-of, {len(layer_cands(x, 'J'))} candidates) \"{x['quote'][:110]}\"")
    EU = collections.defaultdict(list)
    for g in groups:
        for x in gitems[g["id"]]:
            if len(x["cand"]) == 1:
                for b in g["boxes"]:
                    if x["cand"][0] != b: EU[(x["cand"][0], b)].append(x["id"])
    SU = sccs(list(boxes), EU)
    P(f"  strict graph without J (single-candidate items only): {len(EU)} edges; SCCs > 1: {[sorted(c) for c in SU]}")
    # J-only cycles: edges from judged items with a narrowed/unique choice, or unique M items elsewhere
    EJs = collections.defaultdict(list)
    for g in groups:
        for x in gitems[g["id"]]:
            c = x["J"] if isinstance(x.get("J"), list) else (x["cand"] if len(x["cand"]) == 1 else [])
            for p in c:
                for b in g["boxes"]:
                    if p != b: EJs[(p, b)].append(x["id"])
    S2 = sccs(list(boxes), EJs)
    s2 = [sorted(c) for c in S2 if "124.1" in c]
    P(f"  strict graph (J choices where judged, single-candidate items elsewhere): {len(EJs)} edges; SCCs > 1: "
      f"{sorted((len(c) for c in S2), reverse=True)}; SCC of 124.1: {s2[0] if s2 else '-'}")
    for nm, lay, gi2 in (("M MOD", "M", gitems), ("M FIF", "M", layers["FIF"][0]), ("J", "J", gitems)):
        fc = forced_cycles(groups, gi2, gof, boxes, lay)
        P(f"  boxes on a cycle under every resolution of the any-of items ({nm}): {len(fc)} {fc[:40]}{' ...' if len(fc) > 40 else ''}")
    x = next(x for x in gitems[gof["124.1"]] if x["quote"].startswith("the theorem of chapter~\\ref{mod:I-3"))
    y = next(y for y in gitems[gof["18.1"]] if y["quote"].startswith("the explicit two-point function"))
    d124 = definite_dependents(["124.1"], groups, gitems, "M")
    P(f"  witness (M MOD): 124.1 rests on one of {x['cand']} ({x['id']}, 07-cell-3-2-photon.tex:{groups[gof['124.1']]['s_line']}); all of them share the "
      f"derivation of 18.1-18.6, which rests on one of {y['cand']} ({y['id']}, 01-framework.tex:{groups[gof['18.1']]['s_line']}); "
      + ", ".join(f"{c} rests on 124.1: {c == '124.1' or c in d124}" for c in y["cand"]))
    for nm, lay, gi2 in (("M MOD", "M", gitems), ("M FIF", "M", layers["FIF"][0]), ("J", "J", gitems)):
        comps = {}
        for X in forced_cycles(groups, gi2, gof, boxes, lay):
            dX = definite_dependents([X], groups, gi2, lay)
            comp = sorted({X} | {Y for Y in dX if X in definite_dependents([Y], groups, gi2, lay)}, key=lambda n: list(boxes).index(n))
            comps[tuple(comp)] = 1
        P(f"  forced cycle components ({nm}): {[list(c) for c in comps]}")
    cyc = shortest_cycle("124.1", EJs)
    P(f"  shortest strict cycle through 124.1: {cyc}")
    if cyc:
        for p, s in zip(cyc, cyc[1:] + cyc[:1]):
            x = next(x for x in gitems[gof[s]] if p in (x["J"] if isinstance(x.get("J"), list) else x["cand"]))
            g = groups[gof[s]]
            P(f"     {p} -> {s}: item {x['id']} ({g['s_file']}:{g['s_line']}) \"{x['quote'][:120]}\"  [{x['Jnote'] if 'Jnote' in x else 'unique in M'}]")
    for z in ("16.1", "58.1", "18.2", "9.1", "64.1", "66.4", "200.6", "69.1"):
        P(f"     shortest strict cycle through {z}: {shortest_cycle(z, EJs)}")

    # --- F4 at statement level (global)
    P("\n[7] Load-bearing premises (F4) at statement level; T+C = THEOREM and CONDITIONAL boxes (209)")
    claims = {b for b in boxes if boxes[b]["grade"] in ("THEOREM", "CONDITIONAL")}
    sets = [("2.2", ["2.2"]), ("2.4,2.5,2.7,2.8", ["2.4", "2.5", "2.7", "2.8"]), ("9.1,9.7", ["9.1", "9.7"]), ("9.8", ["9.8"]),
            ("18.1,18.3,18.4", ["18.1", "18.3", "18.4"]), ("18.5,18.6", ["18.5", "18.6"]), ("124.1", ["124.1"])]
    P("  | premise set | index transitive T+C | derivation M (upper) | derivation J (upper) | definite M (lower) | definite J (lower) | direct users index / J-definite |")
    P("  |---|---|---|---|---|---|---|---|")
    f4 = {}
    for nm, S0 in sets:
        idx = set().union(*[closure(s, EI) for s in S0]) - set(S0)
        up = set().union(*[closure(s, EM) for s in S0]) - set(S0)
        upj = set().union(*[closure(s, EJ) for s in S0]) - set(S0)
        lo = definite_dependents(S0, groups, gitems, "J")
        lom = definite_dependents(S0, groups, gitems, "M")
        di = {s for (p, s) in EI if p in S0} - set(S0)
        dj = {b for g in groups for x in gitems[g["id"]] if layer_cands(x, "J") and set(layer_cands(x, "J")) <= set(S0) for b in g["boxes"]} - set(S0)
        dm = {b for g in groups for x in gitems[g["id"]] if x["cand"] and set(x["cand"]) <= set(S0) for b in g["boxes"]} - set(S0)
        f4[nm] = dict(index=len(idx & claims), M_upper=len(up & claims), J_upper=len(upj & claims), J_lower=len(lo & claims),
                      direct_index=len(di), direct_J=len(dj))
        P(f"  | {nm} | {len(idx & claims)} | {len(up & claims)} | {len(upj & claims)} | {len(lom & claims)} | {len(lo & claims)} | {len(di)} / {len(dj)} |")
        if nm == "124.1":
            P(f"  124.1, all grades: index direct {len(di)}, transitive {len(idx)}; M-definite direct {len(dm)}, J-definite direct {len(dj)}, definite transitive "
              f"J {len(lo)}, M {len(lom)}; upper J {len(upj)}")
    d124 = sorted({b for g in groups for x in gitems[g['id']] if "124.1" in layer_cands(x, "J") for b in g["boxes"]} - {"124.1"},
                  key=lambda n: list(boxes).index(n))
    P(f"  boxes whose derivation cites 124.1 (J, incl. any-of): {len(d124)} {d124}")

    # --- global layer-M chain-rule tally (all 100 THEOREMs)
    P("\n[8] All 100 THEOREMs, derivation level (J where judged, M elsewhere)")
    th = [b for b in boxes if boxes[b]["grade"] == "THEOREM"]
    tally = collections.Counter(); per = {}
    for b in th:
        g = groups[gof[b]] if b in gof else None
        if g is None or g["sentence"] is None:
            per[b] = "no premise sentence"; tally["no premise sentence"] += 1; continue
        v, _ = verdict(b, gitems[g["id"]], boxes, "J"); per[b] = v; tally[v] += 1
    for k, v in sorted(tally.items()): P(f"  {k}: {v}")
    io = sum(1 for b in th if any(boxes[p]["grade"] == "OPEN" for p in prem_idx[b]))
    P(f"  (index: {io} of 100 list an OPEN premise)")
    for rd in ("FIF", "UNION"):
        gF = layers[rd][0]; tF = collections.Counter()
        for b in th:
            g = groups[gof[b]]
            tF["no premise sentence" if g["sentence"] is None else verdict(b, gF[g["id"]], boxes, "M")[0]] += 1
        P(f"  sensitivity, layer M with reading {rd} (no J): {dict(sorted(tF.items()))}")
    tM = collections.Counter()
    for b in th:
        g = groups[gof[b]]
        tM["no premise sentence" if g["sentence"] is None else verdict(b, gitems[g["id"]], boxes, "M")[0]] += 1
    P(f"  layer M with reading MOD (no J): {dict(sorted(tM.items()))}")
    for k in ("consistent", "no premise sentence", "violates (declared OPEN premise)"):
        P(f"  {k}: {[b for b in th if per[b] == k]}")
    both = sum(1 for b in th if per[b].startswith("violates (declared") and any(boxes[p]["grade"] == "OPEN" for p in prem_idx[b]))
    P(f"  declared OPEN in the derivation and OPEN in the index list: {both}; declared only: "
      f"{sum(1 for b in th if per[b].startswith('violates (declared')) - both}; index only: {io - both}")

    if tsv:
        with open(tsv, "w", encoding="utf-8") as f:
            f.write("group_first_box\tboxes\tsentence_at\titem\tkind\ttarget\tgrade_word\tM_candidates_MOD\tM_candidates_FIF\tJ\tJ_note\tquote\n")
            for g in groups:
                fifi = layers["FIF"][0][g["id"]]
                for x, y in zip(gitems[g["id"]], fifi):
                    j = x.get("J", "unjudged")
                    js = "unjudged" if j == "unjudged" else "unresolved" if j is None else ",".join(j) if j else "none"
                    f.write("\t".join([g["boxes"][0] if g["boxes"] else "", ",".join(g["boxes"]), f"{g['s_file']}:{g['s_line']}", x["id"],
                                       x["kind"], ",".join(x["target"]), str(x["gword"]), ",".join(x["cand"]), ",".join(y["cand"]), js,
                                       x.get("Jnote", ""), norm(x["quote"])]) + "\n")
    return 0

if __name__ == "__main__":
    a = sys.argv[1:]
    t = a[a.index("--tsv")+1] if "--tsv" in a else None
    sys.exit(main(a[0], t))
