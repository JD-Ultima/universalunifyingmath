"""Referee cross-check of compendium round 7 (10-lean-table-rows.md; REVIEW4-CT).

Recomputes from the source tree, with parsers of its own (it imports nothing from lean_tables.py or its helpers):
  [1] the 48 rows' tied names from index/statements.tsv and their declaration kinds in the Lean tree;
  [2] for the 37 table checks: the data module's Row field types, entered rows, check_i proofs and the `checks`
      conjuncts, classified by an independent rule (LEN / KV / CAT / REL);
  [3] the 13 opaque hypotheses: declared opaque, binder unused, and mentioned by no other theorem;
  [4] the 94 opaque Props, records and aliases: mentioned by no theorem of the library;
  [5] the 15 printed numbers said to occur nowhere in the raw text of the Lean tree;
  [6] entered-table facts quoted in the note (SLBC, NCRE, WAMR, ENT2, CASM, ANDB) and, for comparison, the entered
      step table behind round 6's 9.4 (`LemmaB.zeroTransfer`);
  [7] the box text of the rows audited by hand in REVIEW4-CT, verbatim, with their tied names and kinds.
The source is read as text only. Writes nothing in the repository.
Run from this directory (single process, about 5 s):
  python3 -I -B xcheck_round7.py $SRC > outputs/xcheck_round7.txt      # SRC = .../chain-theory/compendium
"""
import os
import re
import sys
from collections import Counter, defaultdict

SRC = sys.argv[1]
LEAN = os.path.join(SRC, "lean", "ChainTheory")
SKIP_FILES = ("AxiomsCheck.lean", "AlphaAxiomsCheck.lean", "Certificate.lean", "Roots.lean")


def lean_files():
    out = []
    for d, _, fs in os.walk(os.path.join(SRC, "lean")):
        out += [os.path.join(d, f) for f in fs if f.endswith(".lean")]
    return sorted(out)


def parse_decls(path):
    """top-level declarations of a file: qualified name -> (keyword, block text)"""
    text = open(path, encoding="utf-8").read()
    ns, decls, cur = [], {}, None
    for line in text.split("\n"):
        m = re.match(r"(namespace|end)\s+([\w.]+)\s*$", line)
        if m:
            if m.group(1) == "namespace":
                ns.append(m.group(2))
            elif ns:
                ns.pop()
            cur = None
            continue
        m = re.match(r"(?:@\[[^\]]*\]\s*)?(?:private\s+|noncomputable\s+)?(theorem|lemma|def|abbrev|opaque|structure|"
                     r"inductive|instance|example)\s+([^\s:(]+)", line)
        if m:
            q = ".".join(ns + [m.group(2)])
            cur = q
            decls[q] = [m.group(1), line]
            continue
        if line[:1] not in ("", " ", "\t") and not line.startswith("--"):
            cur = None
            continue
        if cur:
            decls[cur][1] += "\n" + line
    return decls


def split_top(s, sep):
    parts, d, cur, i = [], 0, "", 0
    while i < len(s):
        if s.startswith(sep, i) and d == 0:
            parts.append(cur)
            cur, i = "", i + len(sep)
            continue
        ch = s[i]
        d += ch in "([{"
        d -= ch in ")]}"
        cur += ch
        i += 1
    parts.append(cur)
    return [p.strip() for p in parts]


def unparen(s):
    s = s.strip()
    while s.startswith("(") and s.endswith(")"):
        d = 0
        for i, ch in enumerate(s):
            d += ch == "("
            d -= ch == ")"
            if d == 0 and i < len(s) - 1:
                return s
        s = s[1:-1].strip()
    return s


def main():
    files = lean_files()
    print("Referee cross-check of compendium round 7 (REVIEW4-CT)")
    print(f"  .lean files in the tree: {len(files)}")
    decls = {}
    file_of = {}
    for f in files:
        for q, v in parse_decls(f).items():
            decls[q] = v
            file_of[q] = f

    # ---------------------------------------------------------------- [1]
    print("\n[1] The 48 rows and their tied names (index/statements.tsv, lean column)")
    rows = [l.rstrip("\n").split("\t") for l in open(os.path.join(SRC, "index", "statements.tsv"), encoding="utf-8")]
    hdr, rows = rows[0], rows[1:]
    col = {h: i for i, h in enumerate(hdr)}
    thm = [r for r in rows if r[col["grade"]] == "THEOREM" and r[col["lean"]].strip()]
    print(f"  THEOREM rows with a lean entry: {len(thm)}")
    names_of = {r[col["number"]]: r[col["lean"]].split() for r in thm}
    restated = set()
    for num, ns in names_of.items():
        for q in ns:
            kw, body = decls.get(q, ("?", ""))
            if kw == "theorem" and ("AlphaDerivation" in body or "Kernel" in body or "LGIS" in body or "ESRL" in body):
                restated.add(num)
    # the 16 are those whose tied set holds a restated AlphaDerivation theorem or a kernel module (round 6)
    r6 = sorted(restated, key=lambda x: [int(t) if t.isdigit() else t for t in re.split(r"[.]", x)])
    print(f"  rows tying a theorem whose statement mentions AlphaDerivation / a kernel module: {len(r6)}: {', '.join(r6)}")
    the48 = [r[col["number"]] for r in thm if r[col["number"]] not in restated]
    sections = {r[col["section"]] for r in thm if r[col["number"]] in the48}
    mult = sum(len(names_of[n]) for n in the48)
    distinct = sorted({q for n in the48 for q in names_of[n]})
    print(f"  the other rows: {len(the48)}; sections {len(sections)}; names with multiplicity {mult}; distinct {len(distinct)}")

    def kind(q):
        kw, body = decls.get(q, ("missing", ""))
        if kw == "theorem" and "type_of%" in body and ".checks" in body:
            head = body.split(":=")[0]
            return "cond_table_check" if "->" in head.split("type_of%")[0] else "table_check"
        if kw == "opaque" and re.search(r":\s*Prop\b", body):
            return "opaque_prop"
        if kw == "def" and "OpenResult" in body.split(":=")[0]:
            return "open_record"
        if kw == "abbrev" and re.search(r":\s*Prop\s*:=", body):
            return "alias_prop"
        return kw
    kinds = {q: kind(q) for q in distinct}
    print(f"  kinds of the distinct names: {dict(sorted(Counter(kinds.values()).items()))}")
    tables = [q for q in distinct if kinds[q] in ("table_check", "cond_table_check")]
    others = [q for q in distinct if kinds[q] in ("opaque_prop", "open_record", "alias_prop")]

    # ---------------------------------------------------------------- [2]
    print("\n[2] The 37 table checks: modules, fields, rows, proofs and conjuncts (independent classifier)")
    forms = Counter()
    by_mod = {}
    str_only = kvtab = 0
    int_mods, sixty, proofs = [], [], Counter()
    rel_mods = set()
    kv_not_one = 0
    for q in tables:
        body = decls[q][1]
        mod = re.search(r"type_of%\s+ChainTheory\.Data\.(D\w+)\.checks", body).group(1)
        path = os.path.join(LEAN, "ChainTheory", "Data", mod + ".lean")
        t = open(path, encoding="utf-8").read()
        st = t[t.index("structure Row where"):t.index("def rows")]
        fields = dict(re.findall(r"^\s+(\w+)\s*:\s*(\w+)", st, re.M))
        ints = {f for f, ty in fields.items() if ty == "Int"}
        nrows = len(re.findall(r"^\s*\{ ", t[t.index("def rows"):], re.M))
        for m in re.finditer(r"theorem check_\d+ :.*?:= (by \w+)", t, re.S):
            proofs[m.group(1)] += 1
        stmt = t[t.index("theorem checks :") + len("theorem checks :"):]
        stmt = stmt[:stmt.index(":=\n")]
        conj = [unparen(c) for c in split_top(stmt.replace("\n", " "), "/\\")]
        fm = []
        for c in conj:
            c = re.sub(r"\s+", " ", c)
            if re.fullmatch(r"rows\.length = \d+", c):
                fm.append("LEN")
            elif re.fullmatch(r'\(rows\.filter \(fun r => r\.key == "[^"]*" && r\.value == "[^"]*"\)\)\.length = 1', c):
                fm.append("KV")
            elif any(re.search(r"\br\." + f + r"\b", c) for f in ints) or "decide (" in c:
                fm.append("REL")
                rel_mods.add(mod)
            else:
                refs = set(re.findall(r"\br\.(\w+)", c))
                assert refs <= {f for f, ty in fields.items() if ty == "String"}, (mod, c)
                if re.search(r'r\.key == "[^"]*" && r\.value == "[^"]*"', c):
                    kv_not_one += 1
                fm.append("CAT")
        forms.update(fm)
        by_mod[mod] = (nrows, Counter(fm), fields)
        if not ints:
            str_only += 1
            kvtab += set(fields) == {"key", "value"}
        else:
            int_mods.append(mod)
        if nrows == 60:
            keys = re.findall(r'key := "([^"]*)"', t)
            detail = [k for k in keys if k.startswith(("cas.", "web."))]
            sixty.append((mod, "last key " + (keys[-1] if keys else "-"),
                          "last cas./web. key " + (detail[-1] if detail else "-")))
    print(f"  table-check declarations: {len(tables)} ({Counter(kinds[q] for q in tables)}); modules: {len(by_mod)}")
    print(f"  conjuncts by form: {dict(forms)}; total {sum(forms.values())}")
    print(f"  check_i proofs: {dict(proofs)}")
    print(f"  modules with String fields only: {str_only} (key/value tables: {kvtab}); with Int fields: {sorted(int_mods)}")
    print(f"  modules with a REL conjunct: {sorted(rel_mods)}")
    print(f"  CAT conjuncts that are key/value filters with a count other than 1: {kv_not_one}")
    print(f"  modules with exactly 60 entered rows: {sixty}")

    # ---------------------------------------------------------------- [3]
    print("\n[3] The 13 opaque hypotheses")
    thm_blocks = {q: v[1] for q, v in decls.items() if v[0] in ("theorem", "lemma", "example")
                  and not file_of[q].endswith(SKIP_FILES)}
    for q in tables:
        if kinds[q] != "cond_table_check":
            continue
        body = decls[q][1]
        h = re.search(r":\s*(\S+)\s*->", body).group(1)
        hq = [d for d in decls if d.split(".")[-1] == h.split(".")[-1]]
        hk = sorted({decls[d][0] for d in hq})
        unused = bool(re.search(r":=\s*fun _ =>", body))
        users = [t for t, b in thm_blocks.items() if t != q and re.search(r"\b" + re.escape(h.split(".")[-1]) + r"\b", b)]
        print(f"  {q.split('.')[-1]:22s} hypothesis {h:22s} declared as {hk} binder unused {unused}; other theorems "
              f"mentioning it: {len(users)}")

    # ---------------------------------------------------------------- [4]
    print("\n[4] The 94 opaque Props, records and aliases: theorems of the library that mention them")
    users = defaultdict(list)
    for q in others:
        s = q.split(".")[-1]
        for t, b in thm_blocks.items():
            if re.search(r"(?<![\w])" + re.escape(s) + r"(?![\w])", b):
                users[q].append(t)
    print(f"  names: {len(others)} {dict(Counter(kinds[q] for q in others))}; mentioned by a theorem: {len(users)}")
    for q, ts in users.items():
        print(f"    {q}: {ts[:3]}")

    # ---------------------------------------------------------------- [5]
    print("\n[5] Printed numbers said to occur nowhere in the raw text of the Lean tree")
    raw = "\n".join(open(f, encoding="utf-8").read() for f in files)
    for n in (17689, 101761, 185496, 322516, 138129, 570349, 101973, 44341, 218340, 93481, 93312, 9270, 46412, 46563,
              270584):
        print(f"  {n:7d}: occurrences {len(re.findall(r'(?<![0-9])' + str(n) + r'(?![0-9])', raw))}")

    # ---------------------------------------------------------------- [6]
    print("\n[6] Entered-table facts")

    def recs_of(mod):
        t = open(os.path.join(LEAN, "ChainTheory", "Data", mod + ".lean"), encoding="utf-8").read()
        out = []
        for m in re.finditer(r"^\s*\{ (.*) \},?\s*$", t, re.M):
            d = {}
            for k, v in re.findall(r'(\w+) := ("[^"]*"|-?\d+)', m.group(1)):
                d[k] = v.strip('"') if v.startswith('"') else int(v)
            out.append(d)
        return out
    sl = recs_of("D525335_SLBC")
    print(f"  SLBC: {len(sl)} rows; W_diff_site = 0 in {sum(r['W_diff_site'] == 0 for r in sl)}; "
          f"W_diff_free = 0 in {sum(r.get('W_diff_free') == 0 for r in sl)}; R_diff_site > 0 in "
          f"{sum(r['R_diff_site'] > 0 for r in sl)}; R_target_site - W_source_site = R_diff_site in "
          f"{sum(r['R_target_site'] - r['W_source_site'] == r['R_diff_site'] for r in sl)}; c_open "
          f"{dict(Counter(r.get('c_open') for r in sl))}")
    nc = recs_of("D525516_NCRE")
    print(f"  NCRE: {len(nc)} rows; order3 = 0 in {sum(r['order3'] == 0 for r in nc)}; SLOT rows outside KIT by group "
          f"{dict(Counter(r['group'] for r in nc if r['polarity'] == 'SLOT' and r['group'] != 'KIT'))}; networks "
          f"{len({(r['group'], r['web']) for r in nc})}; polarities {dict(Counter(r['polarity'] for r in nc))}")
    wa = recs_of("D526019_WAMR")
    print(f"  WAMR: {len(wa)} rows; wash_moved_acts = 0 in {sum(r['wash_moved_acts'] == 0 for r in wa)}; "
          f"sum of acts {sum(r['acts'] for r in wa)}")
    ent2 = open(os.path.join(LEAN, "ChainTheory", "Data", "D530584_ENT2.lean"), encoding="utf-8").read()
    ent2_hit = re.search(r'theorem check_\d+ : \(rows\.filter \(fun r => r\.key == "onsg\.asym\.pairs" && '
                         r'r\.value == "24"\)\)\.length = 1 := by decide', ent2)
    print(f"  ENT2 enters and checks onsg.asym.pairs = \"24\": {bool(ent2_hit)}")
    thermo = open(os.path.join(SRC, "tex", "parts", "13-thermodynamics.tex"), encoding="utf-8").read()
    print(f"  13-thermodynamics.tex contains 'at every one of the 43': {'at every one of the 43' in thermo}")
    alpha = os.path.join(LEAN, "AlphaDerivation", "Data", "Records.lean")
    st = re.findall(r'\("R\d+", \((-?\d+), (-?\d+), (-?\d+)\), \((-?\d+), (-?\d+), (-?\d+)\), \((-?\d+), (-?\d+), (-?\d+)\)\)',
                    open(alpha, encoding="utf-8").read())
    st = [tuple(map(int, x)) for x in st]
    print(f"  round 6, 9.4 (LemmaB.zeroTransfer over AlphaDerivation/Data/Records.lean steps): {len(st)} entered step "
          f"triples; third triple (0,0,0) in {sum(x[6:] == (0, 0, 0) for x in st)}; first + second = 0 in "
          f"{sum(all(x[i] + x[i + 3] == 0 for i in range(3)) for x in st)}; first triples "
          f"{dict(Counter(x[:3] for x in st))}")

    # ---------------------------------------------------------------- [7]
    print("\n[7] Boxes audited by hand in REVIEW4-CT (verbatim) and their tied names")
    tex = {}
    for d, _, fs in os.walk(os.path.join(SRC, "tex")):
        for f in fs:
            if f.endswith(".tex"):
                tex[os.path.relpath(os.path.join(d, f), os.path.join(SRC, "tex"))] = open(os.path.join(d, f),
                                                                                          encoding="utf-8").read()
    title = {r[col["number"]]: r[col["title"]] for r in rows}
    audited = ["9.4", "2.1", "18.2", "217.1", "234.1", "28.1", "28.2", "223.2", "145.1", "3.1", "57.1", "181.5", "H.3",
               "235.4", "170.1", "64.1"]
    for num in audited:
        hits = []
        for fn, t in sorted(tex.items()):
            for m in re.finditer(r"\\begin\{chtheorem\}\[" + re.escape(title[num]) + r"\][^\n]*\n(.*?)\\end\{chtheorem\}",
                                 t, re.S):
                hits.append((fn, t[:m.start()].count("\n") + 1, re.sub(r"\s+", " ", m.group(1)).strip()))
        print(f"  {num} [{title[num]}]  boxes with this title: {len(hits)}")
        for fn, ln, body in hits[:1]:
            print(f"    {fn}:{ln}: {body}")
        print(f"    tied: " + "; ".join(f"{q.split('.')[-1]} ({kind(q)})" for q in names_of.get(num, [])))


if __name__ == "__main__":
    main()
