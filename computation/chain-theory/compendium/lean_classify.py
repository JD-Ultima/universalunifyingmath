#!/usr/bin/env python3
"""Classify the Lean declarations that statements.tsv ties to each boxed statement.

Usage: python3 -I lean_classify.py <compendium_dir> [<number> ...]

Reads compendium/lean/**/*.lean as text (nothing is built or executed) and classifies
every declaration named in the `lean` column of index/statements.tsv:

  opaque_prop      `opaque X : Prop`               (an uninterpreted proposition; no content)
  alias_prop       `abbrev X : Prop := Y`          (renames another proposition)
  open_record      `def X : ...OpenResult := {...}` (a record describing an open result)
  table_check      `theorem X : type_of% <Data module>.checks`  (re-states a data module's checks)
  cond_table_check `theorem X : H -> type_of% <...>.checks` with H an opaque hypothesis
  kernel_check     as table_check but the module lives under ChainTheory/Kernel
  other_theorem    any other theorem
  other / missing

For data modules, it also records whether every field of the `Row` structure is a String
(the checks are then counts of string labels such as "YES"/"EQUAL").
With statement numbers as extra arguments it prints the per-declaration detail for those.
"""
import re, sys, csv, pathlib, collections

def load_decls(leandir):
    decls = {}
    data_string_only = {}
    for p in leandir.rglob("*.lean"):
        if ".lake" in p.parts: continue
        txt = p.read_text(encoding="utf-8", errors="ignore")
        # data modules: Row fields
        m = re.search(r"^namespace\s+(ChainTheory\.(?:Data|Kernel)\.\S+)", txt, re.M)
        if m:
            sm = re.search(r"structure Row where(.*?)deriving", txt, re.S)
            if sm:
                types = re.findall(r"^\s*\w+\s*:\s*(\S+)", sm.group(1), re.M)
                data_string_only[m.group(1)] = bool(types) and all(t == "String" for t in types)
            else:
                data_string_only[m.group(1)] = None
        # walk namespaces
        ns_stack = []
        lines = txt.split("\n")
        for i, line in enumerate(lines):
            mm = re.match(r"^namespace\s+(\S+)", line)
            if mm: ns_stack.append(mm.group(1)); continue
            if re.match(r"^end\s+\S+", line) and ns_stack: ns_stack.pop(); continue
            dm = re.match(r"^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*"
                          r"(theorem|def|abbrev|opaque|lemma|structure|inductive|axiom)\s+([^\s:({\[]+)(.*)", line)
            if not dm: continue
            kind, name, rest = dm.groups()
            full = (ns_stack[-1] + "." + name) if ns_stack else name
            # capture the declaration head up to ':=' (may span 3 lines)
            head = rest + " " + " ".join(lines[i+1:i+4])
            head = head.split(":=")[0]
            decls[full] = (kind, head.strip(), str(p.relative_to(leandir)))
            decls.setdefault(name, (kind, head.strip(), str(p.relative_to(leandir))))
    return decls, data_string_only

def classify(name, decls, dso):
    if name not in decls:
        return "missing", ""
    kind, head, path = decls[name]
    if kind == "opaque" and re.match(r":\s*Prop\b", head): return "opaque_prop", path
    if kind == "abbrev" and re.match(r":\s*Prop\b", head): return "alias_prop", path
    if kind == "def" and "OpenResult" in head: return "open_record", path
    if kind in ("theorem", "lemma"):
        m = re.search(r"type_of%\s+(\S+)\.checks", head)
        if m:
            mod = m.group(1)
            # resolve module alias chains like ChainTheory.Fields.Ch09.row_x
            if ".Kernel." in mod: c = "kernel_check"
            elif "->" in head.split("type_of%")[0]: c = "cond_table_check"
            else: c = "table_check"
            return c + ("" if dso.get(mod) is not True else "[strings]"), path
        if re.search(r"type_of%\s+(\S+)", head):
            return "restated_theorem", path
        return "other_theorem", path
    return "other_" + kind, path

def main(cdir, focus):
    cdir = pathlib.Path(cdir)
    decls, dso = load_decls(cdir/"lean")
    n_data = len(dso); n_str = sum(1 for v in dso.values() if v)
    print(f"Lean declarations indexed: {len(decls)}; data/kernel modules with a Row structure: {n_data}, "
          f"of which every Row field is a String: {n_str}")
    rows = list(csv.DictReader(open(cdir/"index"/"statements.tsv", newline="", encoding="utf-8"), delimiter="\t"))
    rows = [{k: (v or "").strip() for k, v in r.items()} for r in rows]
    by_grade = collections.defaultdict(collections.Counter)
    th_any_substantive = 0
    th_lean = 0
    for r in rows:
        names = r["lean"].split()
        if r["grade"] == "THEOREM" and names:
            th_lean += 1
        cs = [classify(n, decls, dso)[0] for n in names]
        for c in cs: by_grade[r["grade"]][c] += 1
        if r["grade"] == "THEOREM" and any(c.startswith(("other_theorem", "restated_theorem", "kernel_check")) for c in cs):
            th_any_substantive += 1
    print(f"THEOREM rows naming Lean declarations: {th_lean}")
    tiers = collections.Counter()
    for r in rows:
        if r["grade"] != "THEOREM" or not r["lean"]: continue
        cs = {classify(n, decls, dso)[0] for n in r["lean"].split()}
        if any(c.startswith(("restated_theorem", "kernel_check", "other_theorem")) for c in cs):
            tiers["cites a restated (AlphaDerivation) theorem or a kernel module"] += 1
        elif any(c.startswith(("table_check", "cond_table_check")) for c in cs):
            tiers["best evidence: data-table check(s)"] += 1
        else:
            tiers["only opaque Props / aliases / open-result records (no checkable content)"] += 1
    print(f"  the 64 Lean-tied THEOREM rows by strongest declaration kind cited: {dict(tiers)}")
    nosec = collections.Counter(r["grade"] for r in rows if not r["lean"])
    print(f"  rows with no Lean name at all, by grade: {dict(nosec)}")
    # does the Lean column vary within a section? (tie is per chapter, not per statement)
    persec = collections.defaultdict(set)
    for r in rows: persec[r["number"].rsplit(".", 1)[0]].add(r["lean"])
    print(f"  sections whose statements carry different Lean lists: {sum(1 for v in persec.values() if len(v) > 1)} of {len(persec)}"
          " (the tie is per section: every statement of a section inherits the section's names)")
    print(f"THEOREM rows naming at least one declaration that is not opaque/alias/open-record/table-check: {th_any_substantive}")
    for g in ("THEOREM", "CONDITIONAL", "OPEN"):
        print(f"  declaration kinds cited by {g} rows (with multiplicity): {dict(by_grade[g].most_common())}")
    # distinct names cited by theorem rows
    distinct = collections.Counter()
    seen = set()
    for r in rows:
        if r["grade"] != "THEOREM": continue
        for n in r["lean"].split():
            if n in seen: continue
            seen.add(n); distinct[classify(n, decls, dso)[0]] += 1
    print(f"  distinct declarations cited by THEOREM rows: {len(seen)} -> {dict(distinct.most_common())}")
    for num in focus:
        r = next((x for x in rows if x["number"] == num), None)
        if not r: continue
        print(f"\n{num} [{r['grade']}] {r['title']}")
        for n in r["lean"].split():
            c, p = classify(n, decls, dso)
            print(f"   {c:<24} {n}   ({p})")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1], sys.argv[2:]))
