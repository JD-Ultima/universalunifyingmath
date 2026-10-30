#!/usr/bin/env python3
"""Census of compendium/index/statements.tsv (Chain Theory compendium 0.1).

Usage: python3 -I census.py <path/to/statements.tsv>
Reads the TSV as data only (CRLF tolerated). Prints:
  - counts by grade
  - counts by Part x grade, with 'theorems tied to Lean' (non-empty lean column)
  - comparison with the headline numbers 100 / 109 / 198 / 64
"""
import csv, sys, collections

PART_ORDER = ["I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII",
              "XIII","XIV","XV","XVI","App","Instr"]

def load(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    for r in rows:
        for k in list(r):
            r[k] = (r[k] or "").strip()
    return rows

def main(path):
    rows = load(path)
    print(f"rows: {len(rows)}")
    g = collections.Counter(r["grade"] for r in rows)
    print("by grade:", dict(sorted(g.items())))
    tab = collections.defaultdict(collections.Counter)
    for r in rows:
        tab[r["part"]][r["grade"]] += 1
        if r["grade"] == "THEOREM" and r["lean"]:
            tab[r["part"]]["THEOREM_LEAN"] += 1
        if r["grade"] == "CONDITIONAL" and r["lean"]:
            tab[r["part"]]["COND_LEAN"] += 1
    cols = ["THEOREM","THEOREM_LEAN","CONDITIONAL","COND_LEAN","OPEN","TABLE","DEFINITION"]
    print("\n| Part | " + " | ".join(cols) + " |")
    print("|---|" + "---|"*len(cols))
    tot = collections.Counter()
    for p in PART_ORDER:
        c = tab.get(p, collections.Counter())
        tot.update(c)
        print(f"| {p} | " + " | ".join(str(c[k]) for k in cols) + " |")
    print("| all | " + " | ".join(str(tot[k]) for k in cols) + " |")
    parts_only = collections.Counter()
    for p in PART_ORDER[:16]:
        parts_only.update(tab.get(p, collections.Counter()))
    print("\nParts I-XVI only:", {k: parts_only[k] for k in ["THEOREM","CONDITIONAL","OPEN"]})
    cells = collections.Counter()
    for p in ["VI","VII","VIII","IX","X","XI"]:
        cells.update(tab.get(p, collections.Counter()))
    print("Parts VI-XI (the six cells):", {k: cells[k] for k in ["THEOREM","CONDITIONAL","OPEN"]})
    head = {"THEOREM":100, "CONDITIONAL":109, "OPEN":198, "THEOREM_LEAN":64}
    print("\nheadline check:")
    for k, v in head.items():
        print(f"  {k}: headline {v}, computed {tot[k]}, {'MATCH' if v == tot[k] else 'MISMATCH'}")
    # Lean-tie granularity: lean column is shared by all statements of a section?
    sec = collections.defaultdict(set)
    for r in rows:
        ch = r["number"].split(".")[0]
        sec[(r["part"], ch)].add(r["lean"])
    multi = sum(1 for v in sec.values() if len(v) > 1)
    print(f"\nchapters: {len(sec)}; chapters whose statements carry differing lean lists: {multi}")
    return 0

if __name__ == "__main__":
    sys.exit(main(sys.argv[1]))
