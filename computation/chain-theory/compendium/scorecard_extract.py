#!/usr/bin/env python3
"""Mechanical screen for the numeric scorecard (round 2, task 1).

Usage: python3 -I scorecard_extract.py <compendium_dir> [--tsv OUT.tsv]

1. Segments compendium/tex into statement units: a boxed statement plus the prose, bridge
   ("In the language of the literature") and derivation that follow it up to the next box or
   (sub)section. Numbering follows the index (same rule as boxes.py); grades are checked against
   index/statements.tsv. A section's bridge is attached to every box of that (sub)section.
2. Flags a unit as a candidate when its box or its bridge (with citations, cross-reference lists and
   the boilerplate "the record supplies ..." clauses removed) names a physical quantity from a fixed
   lexicon (LEXICON below), and counts the numerals in its box and display equations.
3. Checks every curated scorecard row (scorecard_data.py): each anchor string (a number as printed,
   the identification sentence, the source's own verdict) must occur verbatim in the named unit
   (whitespace-normalised). Reports coverage: which candidate units carry a scorecard row and the
   reason code for each candidate that does not.
Read-only; nothing in the source is executed.
"""
import csv, pathlib, re, sys
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))   # our own directory only

ENV = "chtheorem|chconditional|chwall|chdefinition|chtable"
GRADE = dict(chtheorem="THEOREM", chconditional="CONDITIONAL", chwall="OPEN",
             chdefinition="DEFINITION", chtable="TABLE")


def balanced(s, i, o="{", c="}"):
    d = 0
    for j in range(i, len(s)):
        if s[j] == o: d += 1
        elif s[j] == c:
            d -= 1
            if d == 0: return s[i+1:j], j+1
    return s[i+1:], len(s)


def norm(t):
    return re.sub(r"\s+", " ", t).strip()


def units(texdir):
    """Return the list of statement units in index order."""
    texdir = pathlib.Path(texdir)
    order = ["front"] + [p.stem for p in sorted((texdir/"parts").glob("*.tex"))] + ["back", "instruments"]
    tok = re.compile(r"\\part\*\{|\\appendix\b|\\section\{|\\subsection\{|\\begin\{(" + ENV + r")\}")
    out = []
    secnum, appendix, instr, cur, k = 0, False, False, None, 0
    for f in order:
        p = texdir/(f + ".tex") if f in ("front", "back", "instruments") else texdir/"parts"/(f + ".tex")
        s = p.read_text(encoding="utf-8"); pos = 0; marks = []; group = 0
        while True:
            m = tok.search(s, pos)
            if not m: break
            t = m.group(0)
            if t.startswith("\\part*"):
                title, e = balanced(s, m.end()-1)
                if title.strip() == "Instruments": instr, appendix, secnum = True, False, 0
                marks.append((m.start(), "part", None)); group += 1; pos = e
            elif t.startswith("\\appendix"):
                appendix, secnum, pos = True, 0, m.end()
            elif t.startswith("\\section{") or t.startswith("\\subsection{"):
                title, e = balanced(s, m.end()-1)
                if t.startswith("\\section{"):
                    secnum += 1; k = 0
                    cur = f"I.{secnum}" if instr else (chr(64+secnum) if appendix else str(secnum))
                marks.append((m.start(), "sec", title)); group += 1; pos = e
            else:
                env = m.group(1); k += 1; e = m.end(); opt = ""
                if s[e:e+1] == "[": opt, e = balanced(s, e, "[", "]")
                end = s.find("\\end{" + env + "}", e)
                u = dict(num=f"{cur}.{k}", env=env, grade=GRADE[env], title=opt, body=s[e:end],
                         file=p.name, line=s.count("\n", 0, m.start()) + 1, bend=end, group=(f, group))
                marks.append((m.start(), "box", u)); out.append(u); pos = end + 1
        boxes = [(i, pay) for i, (st, kind, pay) in enumerate(marks) if kind == "box"]
        for i, u in boxes:
            nxt = marks[i+1][0] if i + 1 < len(marks) else len(s)
            u["rest"] = s[u["bend"]:nxt]
            u["last_line"] = s.count("\n", 0, nxt) + 1
        # bridges: per (sub)section group, attached to every box of the group
        bygroup = {}
        for _, u in boxes:
            bygroup.setdefault(u["group"], []).append(u)
        for g, us in bygroup.items():
            br = []
            for u in us:
                br += re.findall(r"\\begin\{chbridge\}(.*?)\\end\{chbridge\}", u["rest"], re.S)
            for u in us:
                u["bridge"] = " ".join(br)
                u["text"] = norm(u["body"] + " " + u["rest"])
    return out


def clean_bridge(b):
    b = re.sub(r"\\cite\{[^}]*\}", "", b)
    b = re.sub(r"(chapter|Section)~\\ref\{[^}]*\}", "", b)
    b = re.sub(r"[Tt]he (record|compendium|construction|derivation) supplies .*?(?=\. |$)", "", b, flags=re.S)
    b = re.sub(r"the chapter titled .*?\)", "", b)
    return norm(b)


LEXICON = [
    ("alpha(0)", r"fine-structure constant"),
    ("a_e, a_mu", r"anomalous (magnetic )?dipole|dipole factor|anomalous magnetic"),
    ("hydrogen levels", r"levels of hydrogen|Rydberg|hydrogen spectrum"),
    ("Lamb shift", r"Lamb shift"), ("hyperfine", r"hyperfine"), ("Casimir", r"Casimir"),
    ("running of alpha", r"coupling of quantum electrodynamics runs|effective coupling"),
    ("alpha_s", r"strong coupling"), ("weak angle", r"weak mixing angle"),
    ("W, Z masses", r"W mass|Z mass|W and the Z|ratio of the W"),
    ("CKM", r"quark mixing matrix|Cabibbo|Jarlskog"), ("PMNS", r"lepton mixing"),
    ("widths", r"width"), ("lifetimes, branching", r"lifetime|branching"),
    ("Fermi constant", r"Fermi constant"),
    ("Higgs mass, quartic, vev", r"Higgs mass|mass of the boson|quartic|vacuum expectation|self-coupling"),
    ("Yukawa", r"Yukawa"), ("lepton masses", r"lepton mass|ratio of two lepton|charged-lepton|Koide"),
    ("quark masses", r"quark mass|mass of the heavy quark|top and W|squared mass of a light quark|squared masses"),
    ("nucleon masses", r"proton.{0,40}mass|mass difference|nucleon"),
    ("hadron masses", r"mass of a hadron|hadron spectroscopy|mesons and baryons"),
    ("binding", r"binding energy|mass defect"), ("generations", r"three generations|generation count"),
    ("G_N", r"gravitational constant|gravitational coupling"),
    ("Lambda", r"cosmological constant|cosmological term"), ("Hubble", r"Hubble"),
    ("BH entropy", r"one quarter"), ("baryon asymmetry", r"baryon asymmetry"),
    ("dark matter", r"dark matter"), ("abundances", r"abundances"),
    ("Bell, Tsirelson", r"Tsirelson|Bell inequality"), ("Chern, Hall", r"Chern|Hall conductance"),
    ("theta", r"theta angle|strong CP"), ("spin-statistics", r"spin-statistics|exchange of two identical"),
    ("signature", r"signature"), ("gauge group", r"gauge group|SU\(3\)"),
    ("photon", r"photon.{0,40}(massless|transverse)"),
    ("decay rates", r"decay rate|beta decay|golden rule"), ("nuclei", r"nucleus|nuclei|deuteron"),
    ("phase transition", r"phase transition"), ("bands", r"energy bands"), ("periodic table", r"periodic table"),
]
LEX = [(k, re.compile(v, re.I)) for k, v in LEXICON]
NUMERAL = re.compile(r"\d+\s*/\s*\d+|\\t?frac|\\sqrt|\d\^|\^\{?-?\d|\b\d+\b")


def screen(U):
    rows = []
    for u in U:
        if u["grade"] == "DEFINITION": continue
        txt = clean_bridge(u["bridge"]) + " " + norm(u["body"])
        keys = [k for k, r in LEX if r.search(txt)]
        if not keys: continue
        eqs = re.findall(r"\\begin\{equation\*?\}(.*?)\\end\{equation\*?\}", u["rest"], re.S)
        n = len(NUMERAL.findall(" ".join(eqs))) + len(NUMERAL.findall(re.sub(r"\\label\{[^}]*\}", "", u["body"])))
        rows.append((u, keys, n))
    return rows


def main():
    src = pathlib.Path(sys.argv[1]); tsv = None
    if "--tsv" in sys.argv: tsv = sys.argv[sys.argv.index("--tsv") + 1]
    U = units(src/"tex"); byn = {u["num"]: u for u in U}
    idx = {r["number"]: r for r in csv.DictReader(open(src/"index"/"statements.tsv", encoding="utf-8", newline=""), delimiter="\t")}
    bad = [u["num"] for u in U if u["num"] not in idx or idx[u["num"]]["grade"] != u["grade"]]
    print(f"[segmentation] {len(U)} boxed statements; grades {dict(Counter(u['grade'] for u in U))}")
    print(f"  index rows {len(idx)}; numbering/grade mismatches against statements.tsv: {len(bad)} {bad[:5]}")
    tex_all = " ".join(norm(u["body"] + " " + u["rest"]) for u in U)
    ncmp = tex_all.count("comparison table")
    ctitles = [u["num"] + " " + u["title"] for u in U if "compar" in u["title"].lower()]
    print(f"  mentions of 'the comparison table' in the statement units: {ncmp}; boxed statements whose title contains 'compar': {ctitles}")
    S = screen(U)
    withnum = [r for r in S if r[2] > 0]
    print(f"\n[screen] units naming a physical quantity (box or cleaned bridge): {len(S)}")
    print(f"  of which with >= 1 numeral in box or display equations: {len(withnum)}; numerals counted: {sum(r[2] for r in S)}")
    print(f"  by grade: {dict(Counter(r[0]['grade'] for r in S))}")

    import scorecard_data as D
    # anchor verification
    nanch, fails = 0, []
    for row in D.ROWS:
        quotes = [row["ident_quote"]] + ([row["verdict"]] if row.get("verdict") else [])
        for stmt, a in row["anchors"] + quotes:
            nanch += 1
            u = byn.get(stmt)
            if u is None or norm(a) not in (u["text"] + " " + norm(u["bridge"])):
                fails.append((row["id"], stmt, a[:70]))
    print(f"\n[anchors] scorecard rows {len(D.ROWS)}, faces {sum(len(r['faces']) for r in D.ROWS)}, "
          f"verbatim anchors checked {nanch} (numbers, identification sentences, verdicts), missing {len(fails)}")
    for f in fails: print("  MISSING", f)
    # grades of scorecard statements
    for row in D.ROWS:
        for stmt in row["stmts"]:
            g = idx[stmt]["grade"]
            if stmt in row.get("grades", {}) and row["grades"][stmt] != g:
                print("  GRADE MISMATCH", row["id"], stmt, row["grades"][stmt], g)
    covered = {s for row in D.ROWS for s in row["stmts"]}
    print(f"\n[coverage] candidate units with numerals: {len(withnum)}")
    print(f"  carrying a scorecard row: {sum(1 for r in withnum if r[0]['num'] in covered)}")
    codes = Counter(); unassigned = []
    for u, keys, n in withnum:
        if u["num"] in covered: continue
        c = D.EXCLUDED.get(u["num"])
        if c is None: unassigned.append(u["num"])
        else: codes[c[0]] += 1
    print(f"  excluded with a reason code: {sum(codes.values())} {dict(sorted(codes.items()))}")
    print(f"  unassigned: {len(unassigned)} {unassigned}")
    extra = sorted(covered - {r[0]['num'] for r in S})
    print(f"  scorecard statements outside the screen (added by reading): {extra}")
    print("\n  reason codes:")
    for k, v in D.CODES.items(): print(f"    {k}: {v}")
    print("\n[candidates] num grade numerals | title | lexicon keys | scorecard row or exclusion code")
    for u, keys, n in S:
        sc = [r["id"] for r in D.ROWS if u["num"] in r["stmts"]]
        tag = ",".join(sc) if sc else ("-" if n == 0 else D.EXCLUDED.get(u["num"], ("?", ""))[0])
        print(f"  {u['num']:7s} {u['grade'][:4]} {n:4d} | {u['title'][:58]:58s} | {','.join(keys)[:48]:48s} | {tag}")
    if tsv:
        with open(tsv, "w", encoding="utf-8", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["number", "grade", "file", "line", "numerals", "lexicon", "scorecard", "excluded"])
            for u, keys, n in S:
                sc = [r["id"] for r in D.ROWS if u["num"] in r["stmts"]]
                w.writerow([u["num"], u["grade"], u["file"], u["line"], n, ";".join(keys), ";".join(sc),
                            "" if sc or n == 0 else D.EXCLUDED.get(u["num"], ("?",))[0]])


if __name__ == "__main__":
    main()
