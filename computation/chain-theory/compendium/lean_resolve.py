#!/usr/bin/env python3
"""Round 6 (proposed further work, item 6): resolve the restated AlphaDerivation theorems and the kernel modules
that compendium theorems are tied to, and re-state the 64-theorem Lean coverage table.

Usage:
  python3 -I lean_resolve.py <compendium_dir> --emit-axioms-file FILE
  python3 -I lean_resolve.py <compendium_dir> [--axioms AX] [--checks CHK] [--build-record REC_B]
                             [--build-record-a REC_A] [--tsv OUT.tsv]

Reads compendium/lean/**/*.lean, compendium/tex and compendium/index as text; nothing of the source is executed.
  --emit-axioms-file FILE  writes a Lean file of ours: `import ChainTheory`, `import AlphaDerivation`, then one
                           `#print axioms` line per name: every declaration named in index/statements.tsv (all grades),
                           every name of the source's own ChainTheory/AxiomsCheck.lean and AlphaAxiomsCheck.lean (read as
                           data, each validated as an identifier), and the Statements conjuncts and Proofs-level theorems
                           resolved below. Then exits.
  --axioms AX              output of `lean FILE` against Build B (lean_audit.sh), parsed per declaration
  --checks CHK             output of lean_audit/ResolveChecks.lean (our own checks)
  --build-record REC_B     lean_build.py record of Build B (all 493 modules, two-line compatibility patch)
  --build-record-a REC_A   lean_build.py record of Build A (published sources, AxiomDeclaration library)
Reuses lean_classify.py (imported, unchanged) for the population and the round-1 classes.
"""
import argparse, collections, csv, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lean_classify as lc          # noqa: E402
import boxes as bx                  # noqa: E402
import lean_resolve_data as D       # noqa: E402
import lean_build as LB             # noqa: E402

DECL_RE = re.compile(r"^\s*(?:@\[[^\]]*\]\s*)?(?:private\s+|protected\s+|noncomputable\s+)*"
                     r"(theorem|def|abbrev|opaque|lemma|structure|inductive|axiom|instance)\s+([^\s:({\[]+)")


def norm(s):
    return re.sub(r"\s+", " ", s).strip()


# ---------------------------------------------------------------- indexing with full namespaces
def index_tree(leandir):
    """full name -> (kind, path, line number (1-based), block text); block = declaration up to the next blank line."""
    idx = {}
    for p in sorted(leandir.rglob("*.lean")):
        if ".lake" in p.parts:
            continue
        lines = p.read_text(encoding="utf-8", errors="ignore").split("\n")
        ns = []
        for i, line in enumerate(lines):
            m = re.match(r"^namespace\s+(\S+)", line)
            if m:
                ns.append(m.group(1)); continue
            m = re.match(r"^end\s+(\S+)", line)
            if m and ns:
                ns.pop(); continue
            dm = DECL_RE.match(line)
            if not dm:
                continue
            j = i
            while j + 1 < len(lines) and lines[j + 1].strip() != "":
                j += 1
            full = ".".join(ns + [dm.group(2)])
            idx.setdefault(full, (dm.group(1), p, i + 1, "\n".join(lines[i:j + 1])))
    return idx


def split_type(block):
    """return (binders text, type text, proof text) of a theorem block: binders before the type colon,
    the type up to the top-level ':=', the proof after it."""
    head = DECL_RE.match(block.split("\n")[0])
    s = block[head.end():]
    depth, colon, assign = 0, None, None
    for k, ch in enumerate(s):
        if ch in "([{⟨":
            depth += 1
        elif ch in ")]}⟩":
            depth -= 1
        elif depth == 0 and ch == ":" and colon is None and s[k:k + 2] != ":=":
            colon = k
        elif depth == 0 and s[k:k + 2] == ":=" and colon is not None:
            assign = k; break
    if colon is None or assign is None:
        return "", "", s
    return s[:colon].strip(), s[colon + 1:assign].strip(), s[assign + 2:].strip()


def proof_method(proof):
    p = norm(proof)
    if p.startswith("by decide +kernel") or p == "by decide +kernel":
        return "decide +kernel"
    if re.fullmatch(r"by decide", p):
        return "decide"
    if p.startswith("by"):
        tacs = sorted(set(re.findall(r"\b(rw|simp only|simp|unfold|constructor|intro|exact|show|decide|omega|rfl|rcases|cases|induction)\b", p)))
        return "tactic proof (" + ", ".join(tacs) + ")"
    if p.startswith("⟨"):
        return "anonymous constructor (" + p[:120] + ")"
    if p == "rfl":
        return "rfl"
    return "term (" + p[:120] + ")"


# ---------------------------------------------------------------- numbers printed in a box
from number_words import numbers_in  # noqa: E402; shared cardinal/ordinal parser


def code_in_title(code, title):
    """the code's letters occur in the title in order, the first one as the initial of a title word"""
    letters = [ch for ch in code.lower() if ch.isalpha()]
    t = title.lower()
    for m in re.finditer(r"\b" + re.escape(letters[0]), t):
        if re.search(".*".join(re.escape(ch) for ch in letters[1:]), t[m.end():]) or len(letters) == 1:
            return True
    return False


def has_number(n, text):
    return re.search(r"(?<![\d])" + str(n) + r"(?![\d])", text) is not None


# ---------------------------------------------------------------- axioms output
def parse_axioms(path):
    out = {}
    txt = pathlib.Path(path).read_text(encoding="utf-8", errors="replace")
    txt = re.sub(r"\n\s+", " ", txt)  # join wrapped continuation lines
    for line in txt.split("\n"):
        m = re.fullmatch(r"'(.+)' depends on axioms: \[([^\]]*)\]\s*", line)
        if m:
            out[m.group(1)] = tuple(sorted(a.strip() for a in m.group(2).split(",") if a.strip()))
            continue
        m = re.fullmatch(r"'(.+)' does not depend on any axioms\s*", line)
        if m:
            out[m.group(1)] = ()
    errors = [l for l in txt.split("\n") if "error" in l.lower()]
    return out, errors


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdir")
    ap.add_argument("--emit-axioms-file")
    ap.add_argument("--axioms")
    ap.add_argument("--checks")
    ap.add_argument("--build-record")
    ap.add_argument("--build-record-a")
    ap.add_argument("--tsv")
    a = ap.parse_args()
    C = pathlib.Path(a.cdir)
    L = C / "lean"
    decls, dso = lc.load_decls(L)
    IDX = index_tree(L)
    rows = list(csv.DictReader(open(C / "index" / "statements.tsv", newline="", encoding="utf-8"), delimiter="\t"))
    rows = [{k: (v or "").strip() for k, v in r.items()} for r in rows]
    byno = {r["number"]: r for r in rows}
    B = bx.boxes(C / "tex")
    cls = lambda n: lc.classify(n, decls, dso)[0]
    SUB = ("restated_theorem", "kernel_check", "other_theorem")

    # population
    th = [r for r in rows if r["grade"] == "THEOREM" and r["lean"]]
    sixteen = [r for r in th if any(cls(n).startswith(SUB) for n in r["lean"].split())]
    tied = sorted({n for r in sixteen for n in r["lean"].split() if cls(n).startswith(SUB)})

    # resolve the conjuncts
    def conjuncts(name):
        blk = IDX[name][3]
        _, ty, _ = split_type(blk)
        cs = re.findall(r"type_of%\s+([A-Za-z0-9_.]+)", ty)
        return ty, cs

    def resolve_stmt(full):
        """AlphaDerivation.Statements.X -> (Statements block, target full name, target block)"""
        k, p, ln, blk = IDX[full]
        _, _, proof = split_type(blk)
        tok = norm(proof).split(" ")[0]
        tgt = "AlphaDerivation." + tok
        return blk, tgt, IDX.get(tgt)

    opaque_props = {n for n, v in IDX.items() if v[0] == "opaque" and re.search(r":\s*Prop\s*$", v[3].split("\n")[0])}
    opaque_short = collections.defaultdict(set)
    for n in opaque_props:
        opaque_short[n.split(".")[-1]].add(n)
    expected = {}
    for n, v in IDX.items():
        if n.startswith("AlphaDerivation.Expected.") and v[0] == "def":
            expected[n.split(".")[-1]] = norm(v[3].split(":=", 1)[1])

    if a.emit_axioms_file:
        names = []
        seen = set()
        def add(n):
            if n not in seen and re.fullmatch(r"[A-Za-z_][A-Za-z0-9_.']*", n):
                seen.add(n); names.append(n)
        for f in ("ChainTheory/AxiomsCheck.lean", "ChainTheory/AlphaAxiomsCheck.lean"):
            for line in (L / f).read_text(encoding="utf-8").split("\n"):
                m = re.fullmatch(r"#print axioms (\S+)", line.strip())
                if m:
                    add(m.group(1))
        for r in rows:
            for n in r["lean"].split():
                if n in IDX:
                    add(n)
        for n in tied:
            add(n)
            if n in IDX and cls(n).startswith("restated"):
                for c in conjuncts(n)[1]:
                    add(c); add(resolve_stmt(c)[1])
        for n in IDX:
            if n.startswith("ChainTheory.Kernel."):
                add(n)
        with open(a.emit_axioms_file, "w", encoding="utf-8") as fh:
            fh.write("import ChainTheory\nimport AlphaDerivation\n\n")
            for n in names:
                fh.write(f"#print axioms {n}\n")
        print(f"wrote {len(names)} #print axioms lines to {a.emit_axioms_file}")
        return 0

    out_rows = []
    P = print
    P("Round 6: restated AlphaDerivation theorems and kernel modules tied to compendium theorems")
    P("=" * 100)
    # ------------------------------------------------------------ [1]
    tiers = collections.Counter()
    for r in th:
        cs = {cls(n) for n in r["lean"].split()}
        if any(c.startswith(SUB) for c in cs): tiers["restated/kernel"] += 1
        elif any(c.startswith(("table_check", "cond_table_check")) for c in cs): tiers["table"] += 1
        else: tiers["opaque/records"] += 1
    P(f"\n[1] Population (lean_classify.py, unchanged): THEOREM rows tied to Lean {len(th)}; round-1 tiers "
      f"restated/kernel {tiers['restated/kernel']}, table checks {tiers['table']}, opaque/records only {tiers['opaque/records']}")
    P(f"    the {len(sixteen)} rows: {', '.join(r['number'] for r in sixteen)}")
    P(f"    distinct restated/kernel declarations they cite: {len(tied)}")
    for n in tied:
        P(f"      {cls(n):<22} {n}  ({IDX[n][1].relative_to(L).as_posix()}:{IDX[n][2]})  cited by "
          + ", ".join(r["number"] for r in sixteen if n in r["lean"].split()))
    allg = collections.defaultdict(list)
    for r in rows:
        for n in r["lean"].split():
            if cls(n).startswith(SUB):
                allg[n].append(f"{r['number']}{r['grade'][0]}")
    P(f"    over all grades, {len(allg)} distinct restated/kernel/other-theorem declarations are cited; those of the 16 rows "
      "are cited by:")
    for n in tied:
        P(f"      {n.split('.')[-1]:<16} {' '.join(allg[n])}")

    # ------------------------------------------------------------ [2]
    P("\n[2] What each declaration states (verbatim type), and what each conjunct resolves to")
    anchor_fail = []
    alpha_txt = norm("\n".join(p.read_text(encoding="utf-8") for p in sorted((L / "ChainTheory" / "AlphaDerivation").rglob("*.lean"))))
    kinds = collections.Counter()
    conj_rows = []
    for n in tied:
        k, p, ln, blk = IDX[n]
        P(f"\n  {n}   ({p.relative_to(L).as_posix()}:{ln})")
        P("    " + blk.replace("\n", "\n    "))
        if cls(n).startswith("kernel"):
            continue
        for c in conjuncts(n)[1]:
            sblk, tgt, tinfo = resolve_stmt(c)
            short = c.split(".")[-1]
            kind, note, anchors = D.CONJ[short]
            kinds[kind] += 1
            tk, tp, tln, tblk = tinfo
            binders, ty, proof = split_type(tblk)
            meth = proof_method(proof)
            bnames = re.findall(r"\((\w+)\s*:", binders)
            btypes = re.findall(r":\s*([A-Za-z0-9_.]+)\)", binders)
            opq = [bt for bt in btypes if bt.split(".")[-1] in opaque_short]
            unused = [b for b in bnames if b.startswith("_") and not re.search(r"(?<![\w])" + re.escape(b) + r"(?![\w])", proof)]
            exps = sorted(set(re.findall(r"Expected\.(\w+)", ty)))
            for an in anchors:
                if norm(an) not in alpha_txt:
                    anchor_fail.append((short, an))
            P(f"\n    conjunct {short}  ->  {tgt}  ({tp.relative_to(L).as_posix()}:{tln})")
            P(f"      statement: {norm(ty)}")
            if binders:
                P(f"      binders: {norm(binders)}" + (f"  [opaque Prop hypotheses: {', '.join(opq)}; unused in the proof: {', '.join(unused)}]" if opq else ""))
            P(f"      proof: {meth}")
            if exps:
                P("      entered constants (Data/Expected.lean): " + "; ".join(f"{e} := {expected.get(e, '?')}" for e in exps))
            P(f"      kind: {kind} -- {note}")
            conj_rows.append((n, short, tgt, kind, meth, len(opq), len(unused)))
    P(f"\n  conjuncts resolved: {len(conj_rows)} (distinct {len({c[1] for c in conj_rows})}); by kind: {dict(kinds)}")
    P(f"  proof methods: {dict(collections.Counter(c[4].split(' (')[0] for c in conj_rows))}")
    P(f"  conjuncts with opaque Prop hypotheses: {sum(1 for c in conj_rows if c[5])} (all such hypotheses unused: "
      f"{all(c[5] == c[6] for c in conj_rows if c[5])})")
    P(f"  Lean anchors checked: {sum(len(v[2]) for v in D.CONJ.values())}, not found: {len(anchor_fail)} {anchor_fail}")
    # entered step table behind LemmaB_zero_transfer
    rec_txt = (L / "ChainTheory/AlphaDerivation/Data/Records.lean").read_text(encoding="utf-8")
    steps = re.findall(r'\("(R\d+)", \(([-\d]+), ([-\d]+), ([-\d]+)\), \(([-\d]+), ([-\d]+), ([-\d]+)\), \(([-\d]+), ([-\d]+), ([-\d]+)\)\)', rec_txt)
    third0 = sum(1 for s in steps if s[7:10] == ("0", "0", "0"))
    neg = sum(1 for s in steps if all(int(s[1 + i]) == -int(s[4 + i]) for i in range(3)))
    P(f"  entered step table (Data/Records.lean): {len(steps)} entries; third triple (0,0,0) as entered in {third0}; "
      f"second = -first in {neg}; first triples {dict(collections.Counter('(' + ','.join(s[1:4]) + ')' for s in steps))}")
    used_cd = [v[1] for nm, v in IDX.items() if v[0] == "theorem" and "centreDims" in v[3]]
    P(f"  Expected.centreDims (the (2,2) and (2,3) centre dimensions, entered [1, 1]) is used by {len(used_cd)} theorems")
    restated_by = collections.Counter()
    for nm, v in IDX.items():
        for c in re.findall(r"type_of%\s+AlphaDerivation\.Statements\.(\w+)", v[3]):
            if nm.startswith("ChainTheory."):
                restated_by[c.split("_")[0]] += 1
    stm_mods = collections.Counter(nm.split(".")[2].split("_")[0] for nm in IDX if nm.startswith("AlphaDerivation.Statements."))
    P(f"  AlphaDerivation.Statements groups restated by chapter declarations (group: restatements / statements): "
      + ", ".join(f"{g}: {restated_by.get(g, 0)}/{stm_mods[g]}" for g in sorted(stm_mods)))

    heavy = [r for r in csv.DictReader(open(L / "build_record.tsv", encoding="utf-8"), delimiter="\t") if float(r["peak_memory_MB"]) >= 4000]
    grp = lambda m: {"IdentificationAlgebra": "IdAlg"}.get(m.split(".")[-1], m.split(".")[-1])
    P(f"  modules at >= 4 GB in the source's build record: {len(heavy)}; of these, AlphaDerivation modules whose Statements "
      "group no chapter declaration restates: " + ", ".join(r["module"].split(".")[-1] for r in heavy
      if r["module"].startswith("AlphaDerivation.Proofs") and restated_by.get(grp(r["module"]), 0) == 0)
      + "; the others: " + ", ".join(r["module"] for r in heavy if not (r["module"].startswith("AlphaDerivation.Proofs")
      and restated_by.get(grp(r["module"]), 0) == 0)))

    # ------------------------------------------------------------ [3] kernel modules
    P("\n[3] The three kernel modules (ChainTheory/Kernel)")
    kmods = sorted({".".join(nm.split(".")[:3]) for nm in IDX if nm.startswith("ChainTheory.Kernel.")})
    for km in kmods:
        rb = IDX.get(km + ".rows")
        rdef = norm(rb[3].split(":=", 1)[1])[:110] if rb else "-"
        lit = rb is not None and rdef.startswith("[")
        chk = IDX[km + ".checks"]
        over = "the rewritten target list (rw [rows_eq_record]) after a computed equation" if "rw [rows_eq_record]" in norm(
            IDX.get(km + ".check_1", ("", "", 0, ""))[3]) else "the computed rows (decide +kernel)"
        cit = collections.Counter()
        citing = []
        for r in rows:
            for n in r["lean"].split():
                if n in decls and decls[n][1].replace(" ", "").find(km.replace(" ", "") + ".checks") >= 0:
                    cit[r["grade"]] += 1; citing.append(r["number"])
        P(f"  {km}  ({IDX[km + '.checks'][1].relative_to(L).as_posix()})")
        P(f"    rows := {rdef}{' ...' if len(rdef) >= 110 else ''}   [{'literal' if lit else 'computed'}]")
        P(f"    checks decided over {over}; checks conjuncts: {norm(split_type(chk[3])[1]).count('/\\') + 1}")
        P(f"    cited by statement rows: {dict(cit)} {citing}")
    P("  inputs entered in KernelDerivationC.Inputs: " + str(sum(1 for nm in IDX if nm.startswith("KernelDerivationC.Inputs.") and IDX[nm][0] == "def")) + " definitions")

    # ------------------------------------------------------------ [4] per-row verdicts
    P("\n[4] Does the tied proposition state the box? (clauses quoted verbatim from the box; counterparts by name)")
    box_fail = []
    verdicts = {}
    tree_txt = "\n".join(v[3] for v in IDX.values())
    for r in sixteen:
        no = r["number"]
        body = norm(B[no]["body"])
        cl = D.ROWS[no]
        names = [n for n in r["lean"].split() if cls(n).startswith(SUB)]
        allowed = set()
        closure = []
        for n in names:
            closure.append(IDX[n][3])
            if cls(n).startswith("restated"):
                for c in conjuncts(n)[1]:
                    allowed.add(c.split(".")[-1])
                    sb, tgt, ti = resolve_stmt(c)
                    closure += [sb, ti[3]]
                    closure += [f"{e} := {expected.get(e, '')}" for e in re.findall(r"Expected\.(\w+)", ti[3])]
            else:
                km = n and decls[n][1].split("type_of%")[1].strip().rsplit(".", 1)[0]
                allowed |= {nm for nm in IDX if nm.startswith(km + ".check_")}
                closure += [v[3] for nm, v in IDX.items() if nm.startswith(km + ".")]
        ctext = "\n".join(closure)
        hit = 0
        for clause, cps, note in cl:
            if norm(clause) not in body:
                box_fail.append((no, clause))
            bad = [c for c in cps if c not in allowed]
            if bad:
                box_fail.append((no, "counterpart not among the tied conjuncts: " + ",".join(bad)))
            hit += bool(cps)
        v = "MATCH" if hit == len(cl) else ("PARTIAL" if hit else "NONE")
        verdicts[no] = v
        nums = sorted(x for x in numbers_in(B[no]["body"]) if x >= 5)
        inside = [x for x in nums if has_number(x, ctext)]
        nowhere = [x for x in nums if not has_number(x, tree_txt)]
        P(f"\n  {no} [{r['grade']}] {r['title']}  --  {v} ({hit} of {len(cl)} clauses)   tied: "
          + ", ".join(n.split(".")[-1] for n in names) + f"  (named for box {D.OWN_BOX[names[0]][0]})")
        for clause, cps, note in cl:
            P(f"     {'+' if cps else '-'} \"{clause}\"" + (f"  <- {', '.join(c.split('.')[-1] for c in cps)}" if cps else "")
              + (f"  [{note}]" if note else ""))
        P(f"     numbers >= 5 printed in the box: {nums}; in the tied statements: {inside}; nowhere in the Lean tree: {nowhere}")
    vc = collections.Counter(verdicts.values())
    P(f"\n  verdicts: {dict(vc)}; box anchors checked: {sum(len(v) for v in D.ROWS.values())}, failures: {box_fail}")

    # ------------------------------------------------------------ [5] attribution
    P("\n[5] Which box each restated/kernel declaration was made for (content, with the 4-letter code as a second sign)")
    for n, (bno, g, ev) in D.OWN_BOX.items():
        ok = "" if bno == "-" or norm(D.OWN_BOX_ANCHORS.get(bno, "")) in norm(B[bno]["body"]) else "  [ANCHOR NOT FOUND]"
        if bno != "-":
            sec = bno.rsplit(".", 1)[0]
            subseq = [b for b in B if b.rsplit(".", 1)[0] == sec and code_in_title(n.split("_")[-1], B[b]["title"])]
        else:
            subseq = []
        gr = byno[bno]["grade"] if bno in byno else "-"
        P(f"  {n.split('.')[-1]:<16} -> {bno:<6} {gr:<12} {ev}{ok}   [boxes of the section whose title contains the code letters in order: {sorted(set(subseq))}]")
    own_c = [n for n, v in D.OWN_BOX.items() if v[1] == "C"]
    P(f"  {len(own_c)} of {len(D.OWN_BOX)} declarations carry the content of a CONDITIONAL box of their section; "
      "they reach THEOREM rows only because the tie is per section")

    # ------------------------------------------------------------ [6] data module provenance
    P("\n[6] Data and kernel modules with a Row structure: are the rows entered or computed?")
    prov = collections.Counter()
    comp_mods = []
    for mod, so in dso.items():
        if so is None:
            continue
        rb = IDX.get(mod + ".rows")
        if not rb:
            prov[("no rows def", so)] += 1; continue
        body = norm(rb[3].split(":=", 1)[1])
        lit = body.startswith("[")
        prov[("entered" if lit else "computed", so)] += 1
        if not lit:
            comp_mods.append((mod, body[:90]))
    P("  (provenance, all fields String) -> modules: " + str(dict(sorted(prov.items()))))
    for m, b in sorted(comp_mods):
        P(f"    computed rows: {m}: rows := {b}")

    # ------------------------------------------------------------ [7] build records
    P("\n[7] Build (Lean 4.35.0, core only; lean_build.py)")
    mods, imps = LB.discover(L)
    need = LB.closure(["ChainTheory", "AlphaDerivation"], imps)
    users = collections.defaultdict(set)
    for m in need:
        for i in imps[m]:
            users[i].add(m)
    def upward(x):
        seen, st = set(), [x]
        while st:
            for z in users[st.pop()]:
                if z not in seen:
                    seen.add(z); st.append(z)
        return seen
    patch_rows = list(csv.DictReader(open(pathlib.Path(__file__).resolve().parent / "lean_build_patch.tsv", encoding="utf-8"), delimiter="\t"))
    pmods = sorted({r["module"] for r in patch_rows})
    pthms = []
    for r in patch_rows:
        line = mods[r["module"]].read_text(encoding="utf-8").split("\n")[int(r["line"]) - 1]
        pthms.append(re.match(r"theorem\s+(\S+)", line).group(1))
    blocked = set().union(*[upward(m) for m in pmods]) - set(pmods)
    first_only = upward("AxiomDeclaration.Theorems")
    lake_files = sorted(q.relative_to(L).as_posix() for q in L.glob("*/lakefile.toml"))
    pins = sorted({(L / q).read_text(encoding="utf-8").strip() for q in [str(x.relative_to(L)) for x in L.glob("*/lean-toolchain")]})
    req = sorted(set(re.findall(r'path\s*=\s*"([^"]+)"', " ".join((L / f).read_text(encoding="utf-8") for f in lake_files))))
    git = sum(("git" in (L / f).read_text(encoding="utf-8")) for f in lake_files)
    roots = collections.Counter(i.split(".")[0] for m in mods for i in imps[m])
    P(f"  lakefiles {len(lake_files)}; toolchain pins {pins}; requirements by path only {req}; git requirements {git}")
    P(f"  import roots of the {len(mods)} library modules: {dict(roots)} (no Mathlib, no Std; `import Lean` only outside the libraries)")
    P(f"  modules of the source's `lake build` (closure of ChainTheory + AlphaDerivation): {len(need)}")
    P(f"  patched modules {pmods}; patched theorems ({len(pthms)}): {pthms}")
    P(f"  modules importing AxiomDeclaration.Theorems transitively: {len(first_only)}; importing a patched module: "
      f"{len(blocked)}; importing neither: {len(need - blocked - set(pmods))} "
      f"(AlphaDerivation modules among them: {sum(1 for m in need - blocked if m.startswith('AlphaDerivation'))})")
    refs = sorted({q.relative_to(L).as_posix() for q in L.rglob("*.lean") for t in pthms
                   if re.search(r"(?<![\w])" + re.escape(t) + r"(?![\w])", q.read_text(encoding="utf-8"))})
    lib_refs = [f for f in refs if f.count("/") >= 2]
    P(f"  files mentioning a patched theorem: {refs}; library modules among them: {lib_refs} "
      "(the others are the source's AxiomsCheck/Roots scripts, outside the libraries)")
    modof = {}
    for nm, v in IDX.items():
        rel = v[1].relative_to(L).parts
        if rel[-1].endswith(".lean") and len(rel) >= 2:
            modof[nm] = ".".join(list(rel[1:-1]) + [rel[-1][:-5]])
    th_names = sorted({n for r in th for n in r["lean"].split() if n in IDX})
    dep = [n for n in th_names if modof.get(n) in blocked]
    P(f"  of the {len(th_names)} declarations cited by the 64 THEOREM rows, {len(dep)} live in modules that import a patched module; "
      f"the 8 restated/kernel declarations: " + ", ".join(f"{n.split('.')[-1]}:{'B only' if modof.get(n) in blocked else 'A and B'}" for n in tied))
    if a.build_record_a:
        ra = list(csv.DictReader(open(a.build_record_a, encoding="utf-8"), delimiter="\t"))
        P(f"  Build A (published sources, unchanged): {sum(1 for r in ra if r['exit'] == '0')} of {len(ra)} modules of the "
          f"AxiomDeclaration library built; failures:")
        for r in ra:
            if r["exit"] not in ("0",):
                P(f"    {r['module']}: exit {r['exit']}: {r['errors'][:300]}")
    if a.build_record:
        rb_ = list(csv.DictReader(open(a.build_record, encoding="utf-8"), delimiter="\t"))
        ok = [r for r in rb_ if r["exit"] == "0"]
        P(f"  Build B (five proofs re-tacticked, lean_build_patch.tsv): {len(ok)} of {len(rb_)} modules built; "
          f"sum of module times {sum(float(r['seconds']) for r in rb_):.0f} s; max peak {max(int(r['peak_rss_MB']) for r in rb_)} MB; "
          f"warnings {sum(int(r['warnings']) for r in rb_)}; 'sorry' warnings {sum(int(r['sorry_warnings']) for r in rb_)}")
        for r in rb_:
            if r["exit"] != "0":
                P(f"    FAILED {r['module']}: exit {r['exit']}: {r['errors'][:300]}")
        P("  patched: " + "; ".join(f"{r['module']} lines {r['patched_lines']}" for r in rb_ if r["patched_lines"]))
        # compare with the source's record
        src = {r["module"]: r for r in csv.DictReader(open(L / "build_record.tsv", encoding="utf-8"), delimiter="\t")}
        common = [r for r in rb_ if r["module"] in src]
        P(f"  source record (4.32.1): {len(src)} modules, {sum(float(r['seconds']) for r in src.values()):.0f} s, "
          f"max {max(float(r['peak_memory_MB']) for r in src.values()):.0f} MB; modules in both records: {len(common)}")

    # ------------------------------------------------------------ [8] axioms
    ax = {}
    if a.axioms:
        ax, errs = parse_axioms(a.axioms)
        P(f"\n[8] #print axioms (our file, against Build B): {len(ax)} declarations reported; error lines {len(errs)}")
        for e in errs[:5]:
            P("    " + e[:200])
        prof = collections.Counter(ax.values())
        P("  axiom sets over all reported declarations: " + ", ".join(f"{'[' + ', '.join(k) + ']' if k else 'none'}: {v}" for k, v in prof.most_common()))
        allax = {x for v in ax.values() for x in v}
        P(f"  axioms occurring: {sorted(allax)}; sorryAx: {sum('sorryAx' in v for v in ax.values())}; "
          f"outside propext/Quot.sound/Classical.choice: {sorted(allax - {'propext', 'Quot.sound', 'Classical.choice'})}")
        th_names = sorted({n for r in th for n in r["lean"].split()})
        miss = [n for n in th_names if n not in ax]
        P(f"  the {len(th_names)} names cited by the 64 THEOREM rows: reported {len(th_names) - len(miss)}, not reported {len(miss)} {miss[:5]}")
        pr = collections.Counter((cls(n).split("[")[0], ax.get(n)) for n in th_names if n in ax)
        for (c, s), v in sorted(pr.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
            P(f"    {c:<18} {('[' + ', '.join(s) + ']') if s else 'none':<42} {v}")
        P("  the restated/kernel declarations and their conjuncts:")
        for n in tied:
            P(f"    {n.split('.')[-1]:<16} {ax.get(n)}")
        for c in conj_rows:
            P(f"      {c[1]:<44} {ax.get('AlphaDerivation.Statements.' + c[1])}  (Proofs: {ax.get(c[2])})")
        # compare with the source's 4.32.1 report
        srcax, _ = parse_axioms(L / "ChainTheory" / "axioms_report.txt")
        both = [n for n in srcax if n in ax]
        same = sum(1 for n in both if srcax[n] == ax[n])
        trans = collections.Counter((srcax[n], ax[n]) for n in both if srcax[n] != ax[n])
        P(f"  against the source's axioms_report.txt (Lean 4.32.1): {len(srcax)} declarations there, {len(both)} also here; "
          f"same axiom set {same}, different {len(both) - same}")
        for (s0, s1), v in trans.most_common(8):
            P(f"    {list(s0)} -> {list(s1)}: {v}")
        P(f"  at 4.32.1, names with sorryAx {sum('sorryAx' in v for v in srcax.values())}, with Classical.choice "
          f"{sum('Classical.choice' in v for v in srcax.values())}")
        P("  the 195 THEOREM-tied names at 4.32.1 (source report): " + str(dict(collections.Counter(
            "[" + ", ".join(srcax[n]) + "]" if srcax.get(n) else ("none" if n in srcax else "absent") for n in th_names))))
        for lab, A in (("4.35.0 (ours)", ax), ("4.32.1 (source)", srcax)):
            u = collections.Counter("[" + ", ".join(sorted({x for n in r["lean"].split() for x in A.get(n, ())})) + "]" for r in th)
            P(f"  per THEOREM row, union of the axioms of its cited names, {lab}: {dict(u)}")
        P(f"  restated/kernel declarations at 4.32.1: " + ", ".join(f"{n.split('.')[-1]} {list(srcax.get(n, ('?',)))}" for n in tied))
    # ------------------------------------------------------------ [9] own checks
    if a.checks:
        ck, errs = parse_axioms(a.checks)
        txt = pathlib.Path(a.checks).read_text(encoding="utf-8", errors="replace")
        P(f"\n[9] Own checks (lean_audit/ResolveChecks.lean): {len(ck)} theorems reported; error lines: {len(errs)}")
        for n, v in ck.items():
            P(f"    {n:<52} axioms {list(v)}")
        for e in errs[:5]:
            P("    " + e[:200])

    # ------------------------------------------------------------ [10] coverage table
    P("\n[10] The 64-theorem coverage table, re-stated")
    newt = collections.Counter()
    for r in th:
        names = r["lean"].split()
        cs = [cls(n) for n in names]
        if any(c.startswith(SUB) for c in cs):
            old = "restated/kernel"
        elif any(c.startswith(("table_check", "cond_table_check")) for c in cs):
            old = "table"
        else:
            old = "opaque"
        no = r["number"]
        if old == "restated/kernel":
            if verdicts[no] in ("MATCH", "PARTIAL"):
                new = "computation stating part of the box" if verdicts[no] == "PARTIAL" else "computation stating the box"
            else:
                rest = [c for c in cs if not c.startswith(SUB)]
                new = ("unrelated computation; best other: table checks" if any(c.startswith(("table_check", "cond_table_check")) for c in rest)
                       else "unrelated computation; otherwise opaque/records only")
        else:
            new = {"table": "table checks (not resolved here)", "opaque": "opaque/records only"}[old]
        newt[new] += 1
        names_ax = [ax.get(n) for n in names if n in IDX]
        axu = sorted({x for v in names_ax if v for x in v})
        n_op = sum(1 for c in cs if c.startswith("opaque_prop"))
        n_cond = sum(1 for c in cs if c.startswith("cond_table_check"))
        n_rec = sum(1 for c in cs if c.startswith("open_record"))
        out_rows.append([no, r["title"], old, new, verdicts.get(no, ""), len(names), n_op, n_cond, n_rec,
                         ",".join(axu) if ax else "", ("yes" if any("sorryAx" in (v or ()) for v in names_ax) else "no") if ax else "",
                         sum(1 for v in names_ax if v is not None) if ax else ""])
    for k, v in newt.most_common():
        P(f"  {k:<58} {v}")
    P(f"  (round-1 tiers: restated/kernel {tiers['restated/kernel']}, table {tiers['table']}, opaque {tiers['opaque/records']})")
    if a.tsv:
        with open(a.tsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["number", "title", "round1_tier", "resolved_tier", "verdict_16", "lean_names", "opaque_props",
                        "opaque_hypothesis_checks", "open_records", "axioms_union_4.35.0", "sorryAx", "names_with_axiom_report"])
            for o in out_rows:
                w.writerow(o)
    return 0


if __name__ == "__main__":
    sys.exit(main())
