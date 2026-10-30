#!/usr/bin/env python3
"""Round 7: the 37 table-check rows and the 11 opaque-only rows of the 64 Lean-tied THEOREMs, resolved clause by clause
as round 6 resolved the 16 restated/kernel rows; then the 64-row content table re-stated.

Usage:
  python3 -I lean_tables.py <compendium_dir> [--axioms AX] [--checks CHK] [--round6-tsv T6]
                            [--tsv OUT.tsv] [--clauses-tsv OUT2.tsv]

Reads compendium/lean/**/*.lean, compendium/tex and compendium/index as text; nothing of the source is executed.
  --axioms AX       round 6's `#print axioms` output against Build B (outputs/lean_axioms.txt), parsed per declaration
  --checks CHK      output of lean_audit/TableChecks.lean (our own checks on the entered rows)
  --round6-tsv T6   round 6's outputs/lean_resolve.tsv, to cross-check the 16 verdicts recomputed here
Imports lean_classify.py, lean_resolve.py, lean_resolve_data.py and boxes.py unchanged; the curated layer of this
round is lean_tables_data.py (every quote checked here).
"""
import argparse, collections, csv, pathlib, re, sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import lean_classify as lc          # noqa: E402
import lean_resolve as lr           # noqa: E402
import lean_resolve_data as D6      # noqa: E402
import lean_tables_data as D        # noqa: E402
import boxes as bx                  # noqa: E402

norm = lr.norm
SUB = ("restated_theorem", "kernel_check", "other_theorem")
TAB = ("table_check", "cond_table_check")
EXCL = ("AxiomsCheck", "Certificate", "Roots", "AlphaAxiomsCheck")   # the source's own audit scripts, not library

from number_words import printed_numbers, legacy_round6_numbers_in  # noqa: E402


# ---------------------------------------------------------------- parsing the source's data modules
def module_info(IDX, mod):
    """mod = 'ChainTheory.Data.D..._XXXX' -> dict(fields, rows, checks_type, conjuncts, check proofs, path)"""
    chk = IDX[mod + ".checks"]
    path = chk[1]
    txt = path.read_text(encoding="utf-8")
    sm = re.search(r"structure Row where(.*?)deriving", txt, re.S)
    fields = dict(re.findall(r"^\s*(\w+)\s*:\s*(\S+)", sm.group(1), re.M)) if sm else {}
    body = txt.split("def rows : List Row := [", 1)[1].split("\n]", 1)[0]
    rows = []
    for line in body.split("\n"):
        line = line.strip().rstrip(",")
        if not line.startswith("{"):
            continue
        d = {}
        for fm in re.finditer(r'(\w+) := ("(?:[^"\\]|\\.)*"|-?\d+)', line):
            v = fm.group(2)
            d[fm.group(1)] = v[1:-1] if v.startswith('"') else int(v)
        rows.append(d)
    _, ty, proof = lr.split_type(chk[3])
    conj = split_conj(ty)
    proofs = re.findall(r"^theorem (check_\d+) : .*?:= (by decide|by decide \+kernel|rfl|\S+)\s*$", txt, re.M)
    return dict(path=path, txt=txt, fields=fields, rows=rows, checks_type=ty, checks_proof=proof, conj=conj,
                check_proofs=proofs, rows_text=body)


def split_conj(ty):
    out, depth, cur = [], 0, ""
    i = 0
    while i < len(ty):
        ch = ty[i]
        if ch in "([{⟨":
            depth += 1
        elif ch in ")]}⟩":
            depth -= 1
        if depth == 0 and ty[i:i + 2] == "/\\":
            out.append(cur); cur = ""; i += 2; continue
        cur += ch; i += 1
    out.append(cur)
    out = [norm(c) for c in out if c.strip()]
    return [c[1:-1].strip() if c.startswith("(") and c.endswith(")") and balanced(c[1:-1]) else c for c in out]


def balanced(s):
    d = 0
    for ch in s:
        d += ch in "([{"
        d -= ch in ")]}"
        if d < 0:
            return False
    return d == 0


def check_kind(c, fields):
    """LEN: rows.length = N; KV: a key/value pair occurs; CAT: filter/all over String fields compared with literals;
    REL: a relation or value test over Int fields"""
    if re.fullmatch(r"rows\.length = \d+", c):
        return "LEN"
    if "r.key ==" in c and "r.value ==" in c:
        return "KV"
    used = set(re.findall(r"\br\.(\w+)", c))
    if any(fields.get(f) in ("Int", "Nat") for f in used):
        return "REL"
    return "CAT"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cdir")
    ap.add_argument("--axioms")
    ap.add_argument("--checks")
    ap.add_argument("--round6-tsv")
    ap.add_argument("--tsv")
    ap.add_argument("--clauses-tsv")
    a = ap.parse_args()
    C = pathlib.Path(a.cdir)
    L = C / "lean"
    decls, dso = lc.load_decls(L)
    IDX = lr.index_tree(L)
    rows = list(csv.DictReader(open(C / "index" / "statements.tsv", newline="", encoding="utf-8"), delimiter="\t"))
    rows = [{k: (v or "").strip() for k, v in r.items()} for r in rows]
    byno = {r["number"]: r for r in rows}
    B = bx.boxes(C / "tex")
    cls = lambda n: lc.classify(n, decls, dso)[0]
    P = print
    th = [r for r in rows if r["grade"] == "THEOREM" and r["lean"]]

    def tier(r):
        cs = {cls(n) for n in r["lean"].split()}
        if any(c.startswith(SUB) for c in cs):
            return "restated/kernel"
        if any(c.startswith(TAB) for c in cs):
            return "table"
        return "opaque"
    tiers = {r["number"]: tier(r) for r in th}
    pop = [r for r in th if tiers[r["number"]] != "restated/kernel"]
    sixteen = [r for r in th if tiers[r["number"]] == "restated/kernel"]
    raw_tree = "\n".join(p.read_text(encoding="utf-8", errors="ignore") for p in sorted(L.rglob("*.lean")) if ".lake" not in p.parts)
    lib_blocks = {nm: v for nm, v in IDX.items() if not any(e in v[1].name for e in EXCL)}
    fails = []

    P("Round 7: the table-check and opaque-only rows of the 64 Lean-tied THEOREMs")
    P("=" * 100)
    # ------------------------------------------------------------ [1] population
    tc = collections.Counter(tiers.values())
    P(f"\n[1] Population (lean_classify.py, unchanged): THEOREM rows tied to Lean {len(th)}; round-1 tiers "
      f"restated/kernel {tc['restated/kernel']}, table {tc['table']}, opaque {tc['opaque']}")
    P(f"    table rows ({tc['table']}): " + ", ".join(r["number"] for r in pop if tiers[r["number"]] == "table"))
    P(f"    opaque rows ({tc['opaque']}): " + ", ".join(r["number"] for r in pop if tiers[r["number"]] == "opaque"))
    names = collections.OrderedDict()
    for r in pop:
        for n in r["lean"].split():
            names.setdefault(n, []).append(r["number"])
    nc = collections.Counter(cls(n) for n in names)
    P(f"    distinct names cited by the {len(pop)} rows: {len(names)} -> {dict(nc.most_common())}")
    P(f"    names cited with multiplicity: {sum(len(r['lean'].split()) for r in pop)}; sections: "
      f"{len({r['number'].rsplit('.', 1)[0] for r in pop})}")

    # ------------------------------------------------------------ [2] the table-check declarations
    P("\n[2] The table-check declarations: verbatim theorem, its module's Row fields, entered rows, the `checks` statement "
      "(verbatim, one conjunct per line) and the kind of each conjunct")
    tabnames = [n for n in names if cls(n).startswith(TAB)]
    MOD = {}
    kinds_all = collections.Counter()
    cond_info = {}
    for n in tabnames:
        k, p, ln, blk = IDX[n]
        mod = re.search(r"type_of%\s+(\S+)\.checks", blk).group(1)
        mi = module_info(IDX, mod)
        MOD[n] = (mod, mi)
        ks = [check_kind(c, mi["fields"]) for c in mi["conj"]]
        kinds_all.update(ks)
        types = collections.Counter(mi["fields"].values())
        P(f"\n  {n}  ({p.relative_to(L).as_posix()}:{ln})  cited by {', '.join(names[n])}")
        P("    " + blk.replace("\n", "\n    "))
        if cls(n).startswith("cond"):
            hyp = re.match(r"(\w+)\s*->", lr.split_type(blk)[1]).group(1)
            hfull = n.rsplit(".", 1)[0] + "." + hyp
            hk = IDX.get(hfull)
            unused = norm(lr.split_type(blk)[2]).startswith("fun _ =>")
            uses = [nm for nm, v in lib_blocks.items() if v[0] == "theorem" and nm != n and re.search(r"(?<![\w])" + re.escape(hyp) + r"(?![\w])", v[3])]
            cond_info[n] = (hyp, hk[0] if hk else "?", unused, uses)
            P(f"    hypothesis: {hk[3].strip() if hk else '?'}  [{hk[0] if hk else '?'}]; proof `fun _ => ...` (binder unused): {unused}; "
              f"other theorems mentioning it: {len(uses)}")
        P(f"    module {mod.split('.')[-1]} ({mi['path'].relative_to(L).as_posix()}): Row fields {dict(types)}; entered rows {len(mi['rows'])}; "
          f"check_i proofs {dict(collections.Counter(x[1] for x in mi['check_proofs']))}; checks := {norm(mi['checks_proof'])[:60]}...")
        P(f"    checks ({len(mi['conj'])} conjuncts; kinds {dict(collections.Counter(ks))}):")
        for c, kk in zip(mi["conj"], ks):
            P(f"      [{kk}] {c}")
    P(f"\n  {len(tabnames)} declarations, {len({MOD[n][0] for n in tabnames})} modules; conjuncts by kind: {dict(kinds_all)}")
    rel_mods = sorted({MOD[n][0].split('.')[-1] for n in tabnames if any(check_kind(c, MOD[n][1]['fields']) == 'REL' for c in MOD[n][1]['conj'])})
    P(f"  modules with at least one REL conjunct: {rel_mods}")
    P(f"  conditional declarations {len(cond_info)}: hypotheses all opaque {all(v[1] == 'opaque' for v in cond_info.values())}, "
      f"all unused {all(v[2] for v in cond_info.values())}, mentioned by other theorems {sum(len(v[3]) for v in cond_info.values())}")
    allp = collections.Counter(x[1] for n in tabnames for x in MOD[n][1]["check_proofs"])
    P(f"  proofs of the check_i lemmas: {dict(allp)}; every `checks` is an And.intro of them: "
      f"{all(norm(MOD[n][1]['checks_proof']).startswith('And.intro check_1') or norm(MOD[n][1]['checks_proof']) == 'check_1' for n in tabnames)}")
    lit = sum(1 for n in tabnames if re.search(r"^def rows : List Row := \[", MOD[n][1]["txt"], re.M))
    so = sum(1 for n in tabnames if set(MOD[n][1]["fields"].values()) == {"String"})
    kvm = sum(1 for n in tabnames if set(MOD[n][1]["fields"]) == {"key", "value"})
    P(f"  rows entered as a literal list: {lit} of {len(tabnames)}; Row fields all String: {so} (key/value tables {kvm}); "
      f"with Int fields: {[MOD[n][0].split('.')[-1] for n in tabnames if 'Int' in MOD[n][1]['fields'].values()]}")
    sixty = [MOD[n][0].split('.')[-1] for n in tabnames if len(MOD[n][1]["rows"]) == 60]
    P(f"  modules with exactly 60 entered rows: {sixty}")

    # ------------------------------------------------------------ [3] opaque Props, aliases, records
    P("\n[3] Opaque Props, aliases and open-result records cited by the 48 rows: what the Lean asserts about them")
    other = [n for n in names if not cls(n).startswith(TAB)]
    oc = collections.Counter()
    for n in other:
        k, p, ln, blk = IDX[n]
        c = cls(n)
        short = n.split(".")[-1]
        th_uses = [nm for nm, v in lib_blocks.items() if v[0] == "theorem" and re.search(r"(?<![\w])" + re.escape(short) + r"(?![\w])", v[3])]
        target = ""
        if c == "alias_prop":
            target = norm(blk.split(":=", 1)[1])
            parts = []
            for t in re.findall(r"[A-Za-z_][\w.]*", target):
                full = t if t in IDX else n.rsplit(".", 1)[0] + "." + t
                parts.append(f"{t} [{IDX[full][0] if full in IDX else '?'}]")
            target = "; ".join(parts)
        if c == "open_record":
            m = re.search(r"statement := ([\w.]+)", blk)
            if m:
                st = m.group(1)
                full = st if st in IDX else n.rsplit(".", 1)[0] + "." + st
                target = f"statement := {st} [{IDX[full][0] if full in IDX else '?'}]"
            else:
                target = "record alias: " + norm(blk.split(":=", 1)[1])[:80]
        oc[(c, bool(th_uses))] += 1
        P(f"  {c:<12} {n}  ({p.relative_to(L).as_posix()}:{ln})  [{', '.join(names[n])}]")
        P(f"      {norm(blk)[:240]}")
        if target:
            P(f"      -> {target}")
        if th_uses:
            P(f"      theorems mentioning it: {th_uses[:4]}")
    P(f"  by (class, mentioned by a library theorem): {dict(oc)}")
    P("  An `opaque X : Prop` has no definition and no proof; a record is a value of type OpenResult whose `statement` field "
      "is such a constant; an alias names one. None of these is a theorem.")

    # ------------------------------------------------------------ [4] per-row verdicts
    P("\n[4] Clause-level verdicts (round 6's rule; clauses quoted verbatim from the box; counterparts by check anchor)")
    P("    rule: a clause has a counterpart when a tied declaration states its finite content under the library's own")
    P("    definitions of its words (kind ENTERED: a relation evaluated over entered per-instance values). A check that an")
    P("    entered string records the clause's verdict, count or letter is RECORDED: not a counterpart; it gives the")
    P("    label-inclusive upper bound. MATCH = every clause has a counterpart, PARTIAL = some, NONE = none.")
    verdict, verdict_l, verdict_a, rowinfo, clause_rows, rel_recorded = {}, {}, {}, {}, [], []

    def tied_text(r):
        chk, rws, oth = [], [], []
        for n in r["lean"].split():
            if n not in IDX:
                continue
            oth.append(IDX[n][3])
            if n in MOD:
                chk.append(MOD[n][1]["checks_type"]); rws.append(MOD[n][1]["rows_text"])
                if n in cond_info:
                    hf = n.rsplit(".", 1)[0] + "." + cond_info[n][0]
                    if hf in IDX:
                        oth.append(IDX[hf][3])
            else:
                blk = IDX[n][3]
                for t in re.findall(r"(?:statement := |:= )([\w.]+)", blk):
                    full = t if t in IDX else n.rsplit(".", 1)[0] + "." + t
                    if full in IDX:
                        oth.append(IDX[full][3])
        return "\n".join(chk), "\n".join(rws), "\n".join(oth)

    def vword(hit, tot):
        return "MATCH" if hit == tot else ("PARTIAL" if hit else "NONE")

    for r in pop:
        no = r["number"]
        body = norm(B[no]["body"])
        cl = D.ROWS.get(no)
        if cl is None:
            fails.append((no, "no clauses")); continue
        tiedn = {n.split(".")[-1]: n for n in r["lean"].split()}
        hit = hitl = hita = 0
        for i, (clause, cps, note) in enumerate(cl):
            if norm(clause) not in body:
                fails.append((no, "clause not in box: " + clause[:60]))
            ok_stated = ok_rec = ok_alt = False
            for dn, anchor, kind in cps:
                full = tiedn.get(dn)
                if full is None or full not in MOD:
                    fails.append((no, f"{dn} is not a tied table check")); continue
                if norm(anchor) not in norm(MOD[full][1]["checks_type"]):
                    fails.append((no, f"anchor not in {dn} checks: {anchor[:60]}"))
                ck = [check_kind(c, MOD[full][1]["fields"]) for c in MOD[full][1]["conj"] if norm(anchor) in c]
                if kind == D.E:
                    ok_stated = True
                    if "REL" not in ck and "LEN" not in ck:
                        fails.append((no, f"ENTERED counterpart is not a REL/LEN check: {anchor[:50]}"))
                elif kind == D.R:
                    ok_rec = True
                    if "REL" in ck:
                        ok_alt = True
                        rel_recorded.append((no, i + 1, dn, anchor))
                else:
                    fails.append((no, "bad kind " + kind))
            hit += ok_stated
            hitl += ok_stated or ok_rec
            hita += ok_stated or ok_alt
            clause_rows.append([no, i + 1, "stated" if ok_stated else ("recorded" if ok_rec else "-"),
                                "; ".join(f"{dn}:{kind}" for dn, _, kind in cps), clause, note])
        verdict[no], verdict_l[no] = vword(hit, len(cl)), vword(hitl, len(cl))
        verdict_a[no] = vword(hita, len(cl))
        chk, rws, oth = tied_text(r)
        nums, ordadd = printed_numbers(B[no]["body"])
        nums5 = sorted(x for x in nums if x >= 5)
        where = {}
        for x in nums5:
            if lr.has_number(x, chk):
                where[x] = "checks"
            elif lr.has_number(x, rws):
                where[x] = "rows"
            elif lr.has_number(x, oth):
                where[x] = "other"
            else:
                where[x] = "-"
        nowhere = [x for x in nums5 if not lr.has_number(x, raw_tree)]
        ccount = collections.Counter(cls(n).split("[")[0] for n in r["lean"].split())
        whole = []
        for n in r["lean"].split():
            if n in MOD:
                for c in MOD[n][1]["conj"]:
                    m = re.search(r'r\.key == "(cell\.C-[A-Z0-9]+)" && r\.value == "([^"]+)"', c)
                    if m and m.group(1).split("C-")[1].startswith("EQ"):
                        whole.append(f"{m.group(1)}={m.group(2)}")
        rowinfo[no] = dict(clauses=len(cl), hit=hit, hitl=hitl, nums5=nums5, where=where, nowhere=nowhere,
                           counts=ccount, whole=whole, ordadd=sorted(ordadd))
        P(f"\n  {no} [{r['grade']}] {r['title']}  ({B[no]['file']}:{B[no]['line']})  --  {verdict[no]} ({hit} of {len(cl)} "
          f"clauses stated); label-inclusive {verdict_l[no]} ({hitl} of {len(cl)})")
        P(f"     tied: {dict(ccount)}; table checks: " + (", ".join(n.split('.')[-1] for n in r['lean'].split() if n in MOD) or "none")
          + (f"; whole-box verdict cells: {whole}" if whole else ""))
        for (clause, cps, note) in cl:
            tag = "+" if any(k == D.E for _, _, k in cps) else ("~" if cps else "-")
            P(f"     {tag} \"{clause}\"" + (f"  <- {', '.join(dn.split('_')[-1] + ':' + k for dn, _, k in cps)}" if cps else "")
              + (f"  [{note}]" if note else ""))
        P(f"     numbers >= 5 printed: {nums5}" + (f" (ordinal words: {sorted(ordadd)})" if ordadd else ""))
        P(f"       in the tied checks: {[x for x in nums5 if where[x] == 'checks']}; in entered rows only: "
          f"{[x for x in nums5 if where[x] == 'rows']}; in other tied declarations only: {[x for x in nums5 if where[x] == 'other']}; "
          f"not in the tied declarations: {[x for x in nums5 if where[x] == '-']}")
        P(f"       nowhere in the Lean tree (raw text of all {len(list(L.rglob('*.lean')))} files): {nowhere}")
    vc, vcl = collections.Counter(verdict.values()), collections.Counter(verdict_l.values())
    ncl = sum(len(D.ROWS[r["number"]]) for r in pop)
    P(f"\n  verdicts (48 rows): {dict(vc)}; label-inclusive: {dict(vcl)}")
    for t in ("table", "opaque"):
        sel = [r["number"] for r in pop if tiers[r["number"]] == t]
        P(f"    {t} rows ({len(sel)}): {dict(collections.Counter(verdict[x] for x in sel))}; label-inclusive "
          f"{dict(collections.Counter(verdict_l[x] for x in sel))}")
    P(f"  clauses {ncl}: stated {sum(1 for c in clause_rows if c[2] == 'stated')}, recorded only "
      f"{sum(1 for c in clause_rows if c[2] == 'recorded')}, none {sum(1 for c in clause_rows if c[2] == '-')}")
    P(f"  rows with a whole-box verdict cell (cell.C-EQ.. in a checked pair): {sum(1 for v in rowinfo.values() if v['whole'])}")
    alln = sum(len(v["nums5"]) for v in rowinfo.values())
    P(f"  printed numbers >= 5 (row, number): {alln}; in tied checks {sum(1 for v in rowinfo.values() for x in v['nums5'] if v['where'][x] == 'checks')}, "
      f"entered rows only {sum(1 for v in rowinfo.values() for x in v['nums5'] if v['where'][x] == 'rows')}, other tied only "
      f"{sum(1 for v in rowinfo.values() for x in v['nums5'] if v['where'][x] == 'other')}, not in tied "
      f"{sum(1 for v in rowinfo.values() for x in v['nums5'] if v['where'][x] == '-')}; nowhere in the tree "
      f"{sum(len(v['nowhere']) for v in rowinfo.values())}")
    P(f"  REL-form checks graded RECORDED (an entered per-instance value that is the clause's conclusion): {rel_recorded}")
    P(f"  sensitivity, counting them as ENTERED (round 6's treatment of an entered conclusion): verdicts "
      f"{dict(collections.Counter(verdict_a.values()))}; rows that change: {[x for x in verdict if verdict_a[x] != verdict[x]]}")
    diffs = []
    for no in [r["number"] for r in pop + sixteen]:
        old_ = sorted(x for x in legacy_round6_numbers_in(B[no]["body"]) if x >= 5)
        new_ = sorted(x for x in printed_numbers(B[no]["body"])[0] if x >= 5)
        if old_ != new_:
            diffs.append((no, old_, new_))
    P("  number extraction, round 6's numbers_in -> this round's parser (numbers >= 5), boxes where they differ:")
    for no, o_, n_ in diffs:
        P(f"    {no:<6} {o_} -> {n_}")
    P(f"  quote checks: {ncl} clauses; failures: {fails}")

    # ------------------------------------------------------------ [5] own derivations from the entered rows
    P("\n[5] Printed numbers against the entered rows (our computation from the rows' text; the source's checks do not state these)")
    gm = lambda code: next(MOD[n][1] for n in MOD if n.endswith(code))
    kvd = lambda code: {x.get("key"): x.get("value") for x in gm(code)["rows"]}
    W = gm("WAMR")["rows"]; S = gm("SLBC")["rows"]; N_ = gm("NCRE")["rows"]
    H = gm("HTOR")["rows"]; RL = gm("RLRO")["rows"]
    casm = sorted({tuple(k.split(".")[4:7]) for k in kvd("CASM") if k.startswith("cas.")}, key=lambda t: (t[0], t[1], int(t[2][1:])))
    andb = sorted({tuple(k.split(".")[1:5]) for k in kvd("ANDB") if k.startswith("web.")})
    pchs = sorted({tuple(k.split(".")[1:5]) for k in kvd("PCHS") if k.startswith("web.")})
    dst2 = sorted({k[len("eqdj."):].rsplit(".", 2)[0] for k in kvd("DST2")})
    ermt = collections.Counter(v for k, v in kvd("ERMT").items() if k.endswith(".result"))
    DER = [
        ("28.1", 228, "WAMR", "sum of acts over the 14 entered bodies", sum(x["acts"] for x in W)),
        ("217.1", 114, "WAMR (tied to 28)", "sum of links", sum(x["links"] for x in W)),
        ("217.1", 78, "WAMR (tied to 28)", "sum of c_open", sum(x["c_open"] for x in W)),
        ("217.1", 150, "WAMR (tied to 28)", "sum of washed", sum(x["washed"] for x in W)),
        ("217.1", 228, "SLBC", "entered reads (rows.length, checked)", len(S)),
        ("217.1", 114, "SLBC", "distinct (body, link)", len({(x["body"], x["link"]) for x in S})),
        ("217.1", 78, "SLBC", "reads with c_open = OPEN", sum(x["c_open"] == "OPEN" for x in S)),
        ("217.1", 150, "SLBC", "reads with c_open = WASHED", sum(x["c_open"] == "WASHED" for x in S)),
        ("217.1", 14, "SLBC", "distinct bodies (13 + U0:NEUTRAL:U)", len({x["body"] for x in S})),
        ("217.1", None, "SLBC", "reads with W unchanged (W_diff_site = W_diff_free = 0)", sum(x["W_diff_site"] == 0 and x["W_diff_free"] == 0 for x in S)),
        ("217.1", None, "SLBC", "reads with R_diff_site > 0", sum(x["R_diff_site"] > 0 for x in S)),
        ("217.1", None, "SLBC", "reads with R_target_site - W_source_site = R_diff_site", sum(x["R_target_site"] - x["W_source_site"] == x["R_diff_site"] for x in S)),
        ("234.1", 64, "NCRE", "KIT SLOT rows defined YES (checked)", sum(x["group"] == "KIT" and x["polarity"] == "SLOT" and x["defined"] == "YES" for x in N_)),
        ("234.1", 18, "NCRE", "SLOT rows outside KIT (BODY 13, NEUTRAL 3, STAR 2)", sum(x["group"] != "KIT" and x["polarity"] == "SLOT" for x in N_)),
        ("234.1", 16, "NCRE", "SLOT rows in BODY and NEUTRAL", sum(x["group"] in ("BODY", "NEUTRAL") and x["polarity"] == "SLOT" for x in N_)),
        ("234.1", None, "NCRE", "rows with order3 = 0", sum(x["order3"] == 0 for x in N_)),
        ("145.1", 44, "CASM", "(chain, separation) pairs in the detail keys: " + " ".join(".".join(t) for t in casm), len(casm)),
        ("225.1", 20, "ANDB", "rings in the per-ring keys", len(andb)),
        ("226.1", 28, "PCHS", "rings in the per-ring keys", len(pchs)),
        ("220.3", 12, "DST2", "saturations in the keys: " + " ".join(dst2), len(dst2)),
        ("D.3", 16, "RLRO", "entered webs (BODY 13, BESIDE 3; checked)", len(RL)),
        ("223.2", 43, "ENT2", "value entered and checked under onsg.asym.pairs", kvd("ENT2")["onsg.asym.pairs"]),
        ("H.3", None, "ERMT", "demonstration results (unchecked meaning)", dict(ermt)),
        ("3.1", None, "HTOR (tied to 57.1)", "OPEN-READ-OF-RECORD rows (cap read, fixed set, fixed, unfixed)",
         [(x["cap_read"], x["fixed_set"], x["fixed_count"], x["unfixed_count"]) for x in H if x["reading"] == "OPEN-READ-OF-RECORD"]),
    ]
    for no, pr, mod, desc, val in DER:
        flag = "" if pr is None else ("  = printed" if str(val) == str(pr) else f"  != printed {pr}")
        P(f"  {no:<6} {mod:<20} {desc}: {val}{flag}")
    andb_last = [k for k in kvd("ANDB")][-1]
    P(f"  last entered key of ANDB: {andb_last}; of CASM detail: {[k for k in kvd('CASM') if k.startswith('cas.')][-1]}")

    # ------------------------------------------------------------ [6] the 16 rows of round 6: their table checks
    P("\n[6] The 16 restated/kernel rows of round 6: their table-check declarations (round 6 graded the restated/kernel ones only)")
    v6, v6l, t6 = {}, {}, collections.Counter()
    for r in sixteen:
        no = r["number"]
        cl = D6.ROWS[no]
        hit6 = sum(1 for _, cps, _ in cl if cps)
        v6[no] = vword(hit6, len(cl))
        tabs = [n for n in r["lean"].split() if cls(n).startswith(TAB)]
        ks = collections.Counter()
        for n in tabs:
            mod = re.search(r"type_of%\s+(\S+)\.checks", IDX[n][3]).group(1)
            mi = module_info(IDX, mod)
            MOD.setdefault(n, (mod, mi))
            ks.update(check_kind(c, mi["fields"]) for c in mi["conj"])
            t6[n] = tuple(sorted(collections.Counter(check_kind(c, mi["fields"]) for c in mi["conj"]).items()))
        extra = D.LABELS16.get(no, {})
        tiedn = {n.split(".")[-1]: n for n in r["lean"].split()}
        hitl = 0
        for i, (clause, cps, note) in enumerate(cl):
            rec = extra.get(i, [])
            for dn, anchor, kind in rec:
                full = tiedn.get(dn)
                if full is None or full not in MOD or norm(anchor) not in norm(MOD[full][1]["checks_type"]):
                    fails.append((no, f"LABELS16 anchor {dn}: {anchor[:50]}"))
            if norm(clause) not in norm(B[no]["body"]):
                fails.append((no, "round-6 clause not in box"))
            hitl += bool(cps) or bool(rec)
        v6l[no] = vword(hitl, len(cl))
        rowinfo[no] = dict(clauses=len(cl), hit=hit6, hitl=hitl)
        P(f"  {no:<6} round 6 {v6[no]:<8} table checks {[n.split('.')[-1] for n in tabs]} kinds {dict(ks)}; "
          f"label-inclusive {v6l[no]} ({hitl} of {len(cl)})"
          + ("; recorded: " + "; ".join(f"clause {i + 1} <- " + ", ".join(dn.split('_')[-1] for dn, _, _ in v) for i, v in extra.items()) if extra else ""))
    k9 = collections.Counter()
    for k in t6.values():
        for x, c in k:
            k9[x] += c
    P(f"  conjunct kinds over these distinct declarations: {dict(k9)}")
    P(f"  distinct table declarations of the 16 rows: {len(t6)}; with a REL conjunct: "
      f"{[n.split('.')[-1] for n, k in t6.items() if any(x == 'REL' for x, _ in k)]}")
    P("  -> no check of these can be a counterpart under the rule: the round-6 verdicts hold over all tied declarations")
    if a.round6_tsv:
        prev = {x["number"]: x["verdict_16"] for x in csv.DictReader(open(a.round6_tsv, encoding="utf-8"), delimiter="\t")}
        diff = [no for no in v6 if prev.get(no) != v6[no]]
        P(f"  recomputed round-6 verdicts against {pathlib.Path(a.round6_tsv).name}: {len(v6) - len(diff)} equal, differing {diff}")
    P(f"  label-inclusive for the 16: {dict(collections.Counter(v6l.values()))} (round 6: {dict(collections.Counter(v6.values()))})")

    # ------------------------------------------------------------ [7] box labels naming a data module
    P("\n[7] Boxes whose \\label names a data module (box:...-D<address>-<code>-...): is that module's declaration tied?")
    tex_all = {f: (C / "tex" / ("parts/" + f if not (C / "tex" / f).exists() else f)) for f in {B[r['number']]['file'] for r in pop}}
    lab_rows = []
    for r in pop:
        no = r["number"]
        f = tex_all[B[no]["file"]]
        lines = f.read_text(encoding="utf-8").split("\n")
        seg = "\n".join(lines[B[no]["line"] - 1:B[no]["line"] + 2])
        m = re.search(r"\\label\{box:[^}]*-D(\d{6})-([a-z0-9]{4})-", seg)
        if not m:
            continue
        addr, code = m.group(1), m.group(2).upper()
        mod = f"ChainTheory.Data.D{addr}_{code}"
        exists = mod + ".checks" in IDX
        citing = sorted(nm for nm, v in lib_blocks.items() if v[0] == "theorem" and (mod + ".checks") in v[3] and not nm.startswith(mod))
        tied = [c for c in citing if c in r["lean"].split()]
        lab_rows.append((no, addr, code, exists, citing, tied))
        P(f"  {no:<6} label D{addr}-{code.lower()}: module {'present' if exists else 'absent'}; declarations citing its checks: "
          f"{[c.split('.', 1)[1] for c in citing]}; tied to the row: {bool(tied)}")
        if exists and not tied:
            mi = module_info(IDX, mod)
            P(f"         its checks: {mi['conj']}")
    # sections whose own-labelled box is elsewhere (I.7.1)
    for mod, (bno, how, ev) in D.OWN_BOX.items():
        ok = ""
        if mod in D.OWN_BOX_ANCHORS and norm(D.OWN_BOX_ANCHORS[mod]) not in norm(B[bno]["body"]):
            ok = "  [ANCHOR NOT FOUND]"; fails.append((mod, "own-box anchor"))
        tied_to = sorted({x for n in MOD if n.endswith(mod.split("_")[1]) and n.split(".")[-1].startswith("row_" + mod[1:7]) for x in names.get(n, [])})
        P(f"  own box of {mod}: {bno} by {how} ({ev}){ok}; the module is tied to {tied_to or 'none of the 48'}")

    # ------------------------------------------------------------ [8] axioms
    ax = {}
    if a.axioms:
        ax, errs = lr.parse_axioms(a.axioms)
        prof = collections.Counter((cls(n).split("[")[0], ax.get(n, ("absent",))) for n in names)
        P(f"\n[8] #print axioms (round 6's file against Build B): the {len(names)} names of the 48 rows")
        for (c, s), v in sorted(prof.items(), key=lambda kv: (kv[0][0], str(kv[0][1]))):
            P(f"    {c:<18} {('[' + ', '.join(s) + ']') if s else 'none':<42} {v}")
        P(f"  sorryAx among them: {sum('sorryAx' in ax.get(n, ()) for n in names)}")
    if a.checks:
        ck, errs = lr.parse_axioms(a.checks)
        P(f"\n[9] Own checks (lean_audit/TableChecks.lean): {len(ck)} theorems reported; error lines: {len(errs)}")
        for n, v in ck.items():
            P(f"    {n:<28} axioms {list(v)}")

    # ------------------------------------------------------------ [10] the 64-row content table
    P("\n[10] The 64 Lean-tied THEOREMs by content (verdicts: round 6 for the 16, this round for the 48)")
    allv = dict(v6); allv.update(verdict)
    allvl = dict(v6l); allvl.update(verdict_l)
    out, ct = [], collections.Counter()
    for r in th:
        no = r["number"]
        t = tiers[no]
        v = allv[no]
        if v == "MATCH":
            content = "states every clause"
        elif v == "PARTIAL":
            content = "states some clause"
        elif t == "restated/kernel":
            content = "states no clause (computation of another box)"
        elif t == "table":
            content = "states no clause (table checks record labels)" if allvl[no] != "NONE" else "states no clause (table checks unrelated)"
        else:
            content = "states no clause (opaque/records only)"
        ct[(t, content)] += 1
        ri = rowinfo[no]
        cs = collections.Counter(cls(n).split("[")[0] for n in r["lean"].split())
        ks = collections.Counter()
        for n in r["lean"].split():
            if n in MOD:
                ks.update(check_kind(c, MOD[n][1]["fields"]) for c in MOD[n][1]["conj"])
        names_ax = [ax.get(n) for n in r["lean"].split() if n in IDX]
        out.append([no, r["title"], t, v, allvl[no], ri["clauses"], ri["hit"], ri["hitl"] - ri["hit"],
                    content, "7" if no in verdict else "6",
                    ",".join(map(str, ri.get("nums5", []))) if no in verdict else "",
                    ",".join(map(str, [x for x in ri.get("nums5", []) if ri["where"][x] != "-"])) if no in verdict else "",
                    ",".join(map(str, ri.get("nowhere", []))) if no in verdict else "",
                    len(r["lean"].split()), cs.get("table_check", 0) + cs.get("cond_table_check", 0), cs.get("cond_table_check", 0),
                    cs.get("opaque_prop", 0), cs.get("alias_prop", 0), cs.get("open_record", 0),
                    " ".join(f"{k}:{ks[k]}" for k in ("LEN", "KV", "CAT", "REL") if ks[k]),
                    ",".join(sorted({x for s in names_ax if s for x in s})) if ax else ""])
    for (t, c), v in sorted(ct.items()):
        P(f"  {t:<16} {c:<52} {v}")
    tot = collections.Counter(allv.values()); totl = collections.Counter(allvl.values())
    P(f"  all 64: {dict(tot)}; label-inclusive: {dict(totl)}")
    P(f"  rows stating at least one clause: {tot['MATCH'] + tot['PARTIAL']} of 64 "
      f"(round 6: {sum(1 for x in v6.values() if x != 'NONE')} of the 16, the 48 not resolved)")
    # sensitivity of the 64-row table to how entered conclusions and partly stated clauses are counted (REVIEW4-CT,
    # R4-CT-8 and R4-CT-9); D.PART_ONLY lists the credited clauses that are stated only in part
    for (no, i), why in D.PART_ONLY.items():
        credited = (any(c[0] == no and c[1] == i and c[2] == "stated" for c in clause_rows) if no in verdict
                    else no in D6.ROWS and bool(D6.ROWS[no][i - 1][1]))
        if not credited:
            fails.append((no, f"PART_ONLY clause {i} is not a credited clause"))
    strict = {no: vword(rowinfo[no]["hit"] - sum(1 for n_, _ in D.PART_ONLY if n_ == no), rowinfo[no]["clauses"])
              for no in allv}
    r6v = dict(allv); r6v.update(verdict_a)
    fmt3 = lambda c: f"{c['MATCH']} / {c['PARTIAL']} / {c['NONE']}"
    P("  sensitivity of the 64-row table (MATCH / PARTIAL / NONE; rows stating at least one clause):")
    for label, vv in (("round 6's treatment: an entered value that is the clause's conclusion counts (the REL checks "
                       "graded RECORDED in [4] count)", r6v),
                      ("headline: round 7's ENTERED/RECORDED split; a clause counts when any part of its finite content "
                       "is stated", allv),
                      ("strict: round 7's split; a clause counts only when all its finite content is stated", strict)):
        c = collections.Counter(vv.values())
        P(f"    {fmt3(c)}; {c['MATCH'] + c['PARTIAL']} rows  {label}"
          + (f"; rows that change: {[x for x in vv if vv[x] != allv[x]]}" if vv is not allv else ""))
    P("    credited clauses stated only in part (not counted under the strict reading):")
    for (no, i), why in D.PART_ONLY.items():
        P(f"      {no:<6} clause {i}: {why}")
    s48 = [no for no in verdict if strict[no] != "NONE"]
    P(f"    the 48 under the strict reading: rows stating a clause {len(s48)} {s48}; clauses stated "
      f"{sum(1 for c in clause_rows if c[2] == 'stated') - sum(1 for n_, _ in D.PART_ONLY if n_ in verdict)}")
    if a.tsv:
        with open(a.tsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["number", "title", "round1_tier", "verdict", "verdict_label_inclusive", "clauses", "clauses_stated",
                        "clauses_recorded_only", "content", "round", "numbers_ge5", "numbers_in_tied", "numbers_nowhere_in_tree",
                        "lean_names", "table_checks", "conditional", "opaque_props", "aliases", "open_records",
                        "check_kinds", "axioms_union_4.35.0"])
            for o in out:
                w.writerow(o)
    if a.clauses_tsv:
        with open(a.clauses_tsv, "w", newline="", encoding="utf-8") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["number", "clause", "status", "counterparts", "text", "note"])
            for c in clause_rows:
                w.writerow(c)
    P(f"\n  quote and anchor failures (all sections): {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
