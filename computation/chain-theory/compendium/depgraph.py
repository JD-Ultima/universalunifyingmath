#!/usr/bin/env python3
"""Premise / used_by dependency graph of the Chain Theory compendium (release 0.1).

Usage:
  python3 -I depgraph.py <compendium_dir> [--json out.json]

<compendium_dir> is the `compendium/` directory of the source repository
(commit a28db8998bce5bd5ef11195c05f13c86cd9dac5c). Everything is read as data.

Inputs
  index/statements.tsv   one row per boxed statement (number, grade, premises, used_by, ...)
  index/designators.tsv  chapter designators ("1-2") -> section labels
  tex/                   LaTeX source, read only to map section labels -> section numbers

Interpretive choices (documented in README.md):
  C1  A premise string is parsed by extracting every statement number n.k; "a to b"
      (same section) expands to every statement a..b of that section. Kind words
      ("Theorem", "Open result") are ignored because the source labels conditional
      statements as "Theorem" in these lists (e.g. "Theorem 124.1" is CONDITIONAL).
  C2  Appendix sections are lettered A..W in statements.tsv but numbered 238..260 in
      the premise text; letter L <-> 237 + (ord(L) - 64). Verified by the grade
      pattern of D <-> 241 and H <-> 245 (see output).
  C3  "the theorems of Chapters X-Y ..." is resolved with designators.tsv: the
      first_section..last_section range when present, else the module_numbering range,
      taking every THEOREM or CONDITIONAL statement in those sections.
  C4  An edge p -> s means "p is listed among the premises of s".
"""
import csv, re, sys, json, pathlib, collections

def load_tsv(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    return [{k: (v or "").strip() for k, v in r.items()} for r in rows]

def balanced(s, i):
    depth = 0
    for j in range(i, len(s)):
        if s[j] == "{": depth += 1
        elif s[j] == "}":
            depth -= 1
            if depth == 0: return s[i+1:j], j+1
    return s[i+1:], len(s)

def section_numbers(texdir):
    """Reproduce LaTeX section numbering: label -> (number, order index)."""
    order = ["front"] + [p.stem for p in sorted((texdir/"parts").glob("*.tex"))] + ["back", "instruments"]
    tok = re.compile(r"\\part\*\{|\\appendix\b|\\section\{")
    secnum = 0; appendix = False; instr = False
    out = collections.OrderedDict()
    for f in order:
        path = texdir/(f + ".tex") if f in ("front", "back", "instruments") else texdir/"parts"/(f + ".tex")
        s = path.read_text(encoding="utf-8")
        pos = 0
        while True:
            m = tok.search(s, pos)
            if not m: break
            t = m.group(0)
            if t.startswith("\\part*"):
                title, e = balanced(s, m.end()-1)
                if title.strip() == "Instruments":
                    instr = True; appendix = False; secnum = 0
                pos = e
            elif t.startswith("\\appendix"):
                appendix = True; secnum = 0; pos = m.end()
            else:
                title, e = balanced(s, m.end()-1)
                secnum += 1
                lab = re.match(r"\s*\\label\{mod:([^}]*)\}", s[e:e+400])
                num = f"I.{secnum}" if instr else (chr(64+secnum) if appendix else str(secnum))
                if lab: out[lab.group(1)] = (num, len(out))
                pos = e
    return out

NUM = re.compile(r"(?<![\w.-])(\d+)\.(\d+)(?:\s+to\s+(\d+)\.(\d+))?")

def to_tsv_section(sec):
    n = int(sec)
    return chr(64 + n - 237) if 238 <= n <= 262 else str(n)

def parse_refs(text, by_section, desig_sections):
    """Return (set of statement numbers, list of unresolved tokens, used chapter form?)"""
    refs, unresolved = set(), []
    chap_form = False
    m = re.search(r"the theorems of Chapters? (.*)", text)
    if m:
        chap_form = True
        body = m.group(1)
        text = text[:m.start()]
        items = [x.strip() for x in re.split(r",|\band\b", body) if x.strip()]
        dlist = list(desig_sections.keys())
        for it in items:
            r = re.match(r"([0-9A-Za-z]+-[0-9A-Za-z]+)(?:\s+to\s+([0-9A-Za-z]+-[0-9A-Za-z]+))?$", it)
            if not r: unresolved.append(it); continue
            a, b = r.group(1), r.group(2) or r.group(1)
            if a not in desig_sections or b not in desig_sections:
                unresolved.append(it); continue
            for d in dlist[dlist.index(a): dlist.index(b)+1]:
                for sec in desig_sections[d]:
                    for row in by_section.get(sec, []):
                        if row["grade"] in ("THEOREM", "CONDITIONAL"):
                            refs.add(row["number"])
    for m in NUM.finditer(text):
        s1, k1, s2, k2 = m.groups()
        sec = to_tsv_section(s1)
        if s2 is None:
            refs.add(f"{sec}.{k1}")
        else:
            if s2 != s1: unresolved.append(m.group(0)); continue
            for k in range(int(k1), int(k2)+1):
                refs.add(f"{sec}.{k}")
    return refs, unresolved, chap_form

def main(cdir, json_out=None):
    cdir = pathlib.Path(cdir)
    rows = load_tsv(cdir/"index"/"statements.tsv")
    by_num = {r["number"]: r for r in rows}
    by_section = collections.defaultdict(list)
    for r in rows:
        by_section[r["number"].rsplit(".", 1)[0]].append(r)
    secmap = section_numbers(cdir/"tex")
    labels = list(secmap.keys())
    # check numbering reproduces statements.tsv section numbers
    sec_in_tsv = set(by_section)
    sec_from_tex = {v[0] for v in secmap.values()}
    print(f"sections with statements: {len(sec_in_tsv)}; of them found in tex numbering: {len(sec_in_tsv & sec_from_tex)}")
    # appendix letter check (C2)
    for L in "DH":
        print(f"C2 check: appendix {L} <-> {237 + ord(L) - 64}: grades {[r['grade'][0] for r in by_section[L]]}")
    # designators -> section numbers
    desig = collections.OrderedDict()
    amb = 0
    for d in load_tsv(cdir/"index"/"designators.tsv"):
        a, b = d["first_section"], d["last_section"]
        if d["first_section"] and d["module_numbering_first"]: amb += 1
        if not a: a, b = d["module_numbering_first"], d["module_numbering_last"]
        secs = []
        if a in secmap and b in secmap:
            i, j = secmap[a][1], secmap[b][1]
            secs = [secmap[labels[x]][0] for x in range(min(i, j), max(i, j)+1)]
        desig[d["number"]] = secs
    print(f"designators: {len(desig)}; with two readings (both numberings populated): {amb}")

    prem, used, unres = {}, {}, []
    chapform = 0
    for r in rows:
        p, u1, cf = parse_refs(r["premises"], by_section, desig)
        q, u2, _ = parse_refs(r["used_by"], by_section, desig)
        chapform += cf
        prem[r["number"]] = p; used[r["number"]] = q
        unres += [(r["number"], x) for x in u1 + u2]
    missing = sorted({x for s in prem.values() for x in s if x not in by_num} |
                     {x for s in used.values() for x in s if x not in by_num})
    print(f"rows: {len(rows)}; rows with chapter-form premises: {chapform}; unresolved tokens: {len(unres)}; "
          f"referenced numbers not in statements.tsv: {len(missing)} {missing[:10]}")
    for r in rows:
        prem[r["number"]] = {x for x in prem[r["number"]] if x in by_num and x != r["number"]}
        used[r["number"]] = {x for x in used[r["number"]] if x in by_num and x != r["number"]}

    # graph
    edges = {(p, s) for s, ps in prem.items() for p in ps}
    print(f"edges (premise -> statement): {len(edges)}")
    # reciprocity: p in prem[s]  <=>  s in used[p]
    recip = sum(1 for (p, s) in edges if s in used[p])
    print(f"reciprocity: {recip}/{len(edges)} premise edges are echoed in the premise's used_by "
          f"({100*recip/max(1,len(edges)):.1f}%)")
    uedges = {(p, s) for p, ss in used.items() for s in ss}
    recip2 = sum(1 for (p, s) in uedges if p in prem[s])
    print(f"reverse reciprocity: {recip2}/{len(uedges)} used_by edges are echoed in the user's premises")

    children = collections.defaultdict(set)
    for (p, s) in edges: children[p].add(s)
    def descendants(x):
        seen, st = set(), [x]
        while st:
            y = st.pop()
            for z in children[y]:
                if z not in seen: seen.add(z); st.append(z)
        return seen
    def ancestors(x):
        seen, st = set(), [x]
        while st:
            y = st.pop()
            for z in prem[y]:
                if z not in seen: seen.add(z); st.append(z)
        return seen

    # cycles (Tarjan)
    idx, low, onst, st, sccs, c = {}, {}, set(), [], [], [0]
    sys.setrecursionlimit(10000)
    def strong(v):
        idx[v] = low[v] = c[0]; c[0] += 1; st.append(v); onst.add(v)
        for w in children[v]:
            if w not in idx: strong(w); low[v] = min(low[v], low[w])
            elif w in onst: low[v] = min(low[v], idx[w])
        if low[v] == idx[v]:
            comp = []
            while True:
                w = st.pop(); onst.discard(w); comp.append(w)
                if w == v: break
            if len(comp) > 1: sccs.append(sorted(comp))
    for v in by_num:
        if v not in idx: strong(v)
    print(f"strongly connected components of size > 1 (circular premise chains): {len(sccs)}; "
          f"sizes {sorted((len(s) for s in sccs), reverse=True)[:10]}")
    for s in sorted(sccs, key=len, reverse=True)[:3]:
        print("   e.g.", s[:12], "..." if len(s) > 12 else "")

    G = lambda x: by_num[x]["grade"]
    # chain-rule audit
    th = [r["number"] for r in rows if r["grade"] == "THEOREM"]
    th_direct = [t for t in th if any(G(p) in ("OPEN", "CONDITIONAL") for p in prem[t])]
    th_direct_open = [t for t in th if any(G(p) == "OPEN" for p in prem[t])]
    th_trans_open = [t for t in th if any(G(p) == "OPEN" for p in ancestors(t))]
    th_nopremise = [t for t in th if not prem[t]]
    print(f"\nTHEOREMs: {len(th)}; with no parsed premise: {len(th_nopremise)}")
    print(f"THEOREMs listing >=1 OPEN or CONDITIONAL premise directly: {len(th_direct)}")
    print(f"THEOREMs listing >=1 OPEN premise directly: {len(th_direct_open)}")
    print(f"THEOREMs with an OPEN statement among their transitive premises: {len(th_trans_open)}")
    co = [r["number"] for r in rows if r["grade"] == "CONDITIONAL"]
    co_open = [t for t in co if any(G(p) == "OPEN" for p in prem[t])]
    print(f"CONDITIONALs listing >=1 OPEN premise directly: {len(co_open)}/{len(co)}")

    # load-bearing unproven premises
    claims = [r["number"] for r in rows if r["grade"] in ("THEOREM", "CONDITIONAL")]
    def stat(x):
        d = descendants(x)
        return dict(number=x, grade=G(x), title=by_num[x]["title"], part=by_num[x]["part"],
                    direct=len(children[x]),
                    direct_claims=sum(1 for y in children[x] if G(y) in ("THEOREM", "CONDITIONAL")),
                    trans_claims=sum(1 for y in d if G(y) in ("THEOREM", "CONDITIONAL")),
                    trans_theorems=sum(1 for y in d if G(y) == "THEOREM"))
    unproven = [stat(x) for x in by_num if G(x) in ("OPEN", "CONDITIONAL")]
    unproven.sort(key=lambda d: (-d["trans_claims"], -d["direct_claims"], d["number"]))
    print(f"\nTop unproven (OPEN/CONDITIONAL) premises by number of THEOREM+CONDITIONAL statements "
          f"resting on them transitively (of {len(claims)}):")
    print("| No. | grade | Part | direct users (T+C) | transitive T+C | transitive T | title |")
    print("|---|---|---|---|---|---|---|")
    for d in unproven[:25]:
        print(f"| {d['number']} | {d['grade']} | {d['part']} | {d['direct_claims']} | {d['trans_claims']} | "
              f"{d['trans_theorems']} | {d['title']} |")
    allst = [stat(x) for x in by_num]
    allst.sort(key=lambda d: (-d["direct"], d["number"]))
    print("\nTop statements by direct in-degree (times listed as a premise):")
    for d in allst[:12]:
        print(f"  {d['number']:>6} {d['grade']:<11} in-degree {d['direct']:>3}  {d['title']}")
    roots = [x for x in by_num if not prem[x] and children[x]]
    print(f"\nroots (no premise, >=1 user): {len(roots)} -> grades {dict(collections.Counter(G(x) for x in roots))}")

    # physics claims of the abstract
    focus = ["35.1", "40.1", "41.1", "53.1", "54.1", "59.1", "61.1", "64.1", "64.2", "65.1", "74.1",
             "83.1", "124.1", "141.1", "143.1", "146.1", "147.1", "163.1", "171.1", "172.1", "173.1",
             "185.1", "197.1", "197.3", "197.4", "197.5", "198.3", "206.3", "90.2", "115.1", "217.1",
             "220.3", "221.1", "225.1"]
    print("\nPhysics claims: direct premises by grade, and OPEN premises listed")
    out_focus = {}
    for x in focus:
        if x not in by_num: continue
        ps = sorted(prem[x], key=lambda s: [int(t) if t.isdigit() else ord(t[0]) * 1000 for t in re.split(r"\.", s)])
        gc = collections.Counter(G(p) for p in ps)
        anc = ancestors(x)
        ga = collections.Counter(G(p) for p in anc)
        op = [p for p in ps if G(p) == "OPEN"]
        cd = [p for p in ps if G(p) == "CONDITIONAL"]
        out_focus[x] = dict(grade=G(x), direct=dict(gc), transitive=dict(ga), open=op, conditional=cd)
        print(f"  {x} [{G(x)}] direct {dict(gc)}; transitive {dict(ga)}\n      CONDITIONAL premises: {cd}\n      OPEN premises: {op}")
    if json_out:
        json.dump(dict(edges=sorted(edges), unproven_top=unproven[:40], focus=out_focus,
                       sccs=sccs, theorems_with_open_or_conditional_premise=th_direct),
                  open(json_out, "w"), indent=1)
    return 0

if __name__ == "__main__":
    a = sys.argv[1:]
    j = a[a.index("--json")+1] if "--json" in a else None
    sys.exit(main(a[0], j))
