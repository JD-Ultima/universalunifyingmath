"""Referee cross-check of the mechanical trial-factor space (F-CT-4 E; REVIEW4-CT).

Checks `mechanical_mutations.py` / `mechanical_mutations.out` / `mechanical-trial-factor.md` without rerunning the
full evaluation:
  [1] the readout cone, recomputed by an independent regex closure over Part 1 of the pinned ChainAlpha.lean;
  [2] coarse site counts per operator class, recounted by regular expressions over the cone bodies;
  [3] every example mutant printed in out section 4 (closest values, the author's rewordings, Table S5), re-evaluated
      with the transcription; value and 1/value recomputed (mpmath, 50 digits, own code);
  [4] a seeded uniform sample of 2000 of the 332 308 mutants: outcome proportions against the counts of out section 3;
  [5] the Born-weight value 319/31104: every pair through `recs:REL:0=neg` or `sameForm:REL:0=neg`;
  [6] depth dependence: seeded samples of 12 000 pairs and 12 000 triples (pairwise-disjoint sites, same rule as the
      generator's pairs): distinct values with 1/value in (120,155), singleton/doubleton counts and the Chao1
      lower-bound estimate of the number of such values in the whole depth-k space;
  [7] optional (LEAN set): Lean 4.35.0 on 6 pairs outside the 300-pair spot-check sample and 4 sampled triples.
  [8] the pairs excluded by the disjoint-site rule, by reason (the note says: 'inside a branch it moves').
The arithmetic of the F-CT-4 metrics from the printed counts is redone in REVIEW4-CT, not here.

Imports mechanical_mutations.py (setup, compatible, canonical, lean_text) and chainalpha_transcription.py unchanged;
reads the pinned upstream file as text only. Writes nothing in the repository; Lean files go to a temporary directory.
Single process. Run from this directory:
  python3 -I -B xcheck_mechanical.py /path/to/ChainAlpha.lean > xcheck_mechanical.out
  LEAN=/path/to/lean-4.35.0/bin/lean python3 -I -B xcheck_mechanical.py /path/to/ChainAlpha.lean > xcheck_mechanical.out
About 3 min (Lean part about 15 s).
"""
import os
import random
import re
import subprocess
import sys
import tempfile
import time
from collections import Counter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True
import mpmath as mp                       # noqa: E402
import mechanical_mutations as M          # noqa: E402
import chainalpha_transcription as T      # noqa: E402

mp.mp.dps = 50
INV_CODATA = mp.mpf("137.035999177")
SIG = mp.mpf("2.1e-8")
SEED = 20261030


def inv_exact(key):
    """1/value of a canonical key (N, ((m, a), ...)), own code"""
    N, items = key
    return mp.mpf(N) / mp.fsum(a * mp.sqrt(m) for m, a in items)


# ---------------------------------------------------------------- [1] cone, independently
def strip_comments(s):
    out, i, depth = [], 0, 0
    while i < len(s):
        if depth == 0 and s.startswith("--", i):
            j = s.find("\n", i)
            i = len(s) if j < 0 else j
        elif s.startswith("/-", i):
            depth, i = depth + 1, i + 2
        elif depth and s.startswith("-/", i):
            depth, i = depth - 1, i + 2
        else:
            if depth == 0:
                out.append(s[i])
            elif s[i] == "\n":
                out.append("\n")
            i += 1
    return "".join(out)


def part1_items(src):
    a, b = src.index("/-! # Part 1"), src.index("/-! ## 1.9")
    text = strip_comments(src[a:b])
    items, ns, cur = {}, None, None
    for line in text.split("\n"):
        m = re.match(r"(namespace|end|def|structure|abbrev)\s+([\w.]+)", line)
        if m:
            kw, nm = m.groups()
            if kw == "namespace":
                ns, cur = nm, None
                continue
            if kw == "end":
                ns, cur = None, None
                continue
            full = (ns + "." + nm) if (ns and kw == "def") else nm
            cur = full
            items[cur] = [kw, line[m.end():]]
            continue
        if cur:
            items[cur][1] += "\n" + line
    return items


def regex_cone(items, roots):
    defs = {n for n, (kw, _) in items.items() if kw == "def"}
    short = {}
    for n in defs:
        short.setdefault(n.split(".")[-1], set()).add(n)
    seen, todo = set(), list(roots)
    while todo:
        n = todo.pop()
        if n in seen:
            continue
        seen.add(n)
        for ident in re.findall(r"[A-Za-z_][\w']*(?:\.[A-Za-z_][\w']*)*", items[n][1]):
            parts = ident.split(".")
            for k in range(len(parts)):
                tail = ".".join(parts[k:])
                if tail in defs:
                    todo.append(tail)
            for p in parts[1:]:                      # dot notation w.typeOf, t.valid (over-approximation)
                todo += [d for d in short.get(p, ()) if "." in d]
    return seen


# ---------------------------------------------------------------- [2] coarse counts
def coarse_counts(items, cone):
    c = Counter()
    for n in cone:
        body = items[n][1]
        body = body[body.index(":=") + 2:] if (":=" in body and not re.match(r"[^\n]*\n\s*\|", body)) else body
        body = re.sub(r"(?<!\|)\|(?!\|)[^=|\n]*=>", " ", body)          # drop match patterns `| ... =>`
        c["IF"] += len(re.findall(r"\bif\b", body))
        c["MATCH"] += len(re.findall(r"\bmatch\b", body)) + (1 if re.match(r"[^\n]*\n\s*\|", items[n][1]) and
                                                                 ":=" not in items[n][1].split("\n")[0] else 0)
        c["LOG"] += len(re.findall(r"&&|\|\||∧|∨", body))
        c["QNT"] += len(re.findall(r"\.(?:all|any)\b", body))
        c["BLIT"] += len(re.findall(r"\b(?:true|false)\b", body))
        c["LIT"] += len(re.findall(r"(?<![\w.])\d+", body))
        c["negative LIT"] += len(re.findall(r"(?:[(\[,]|then|else|:=)\s*-\d", body))
        c["REL"] += len(re.findall(r"==|!=|≠|≤|(?<![:=!<>])=(?![=>])|<(?![|;=])|\beqZ\b|\bneZ\b|\.isSome\b", body))
        c["++"] += len(re.findall(r"\+\+", body))
        c["::"] += len(re.findall(r"::", body))
        c["list literals"] += len(re.findall(r"\[(?!\])", body))
        c["AOR (+ * / % min max)"] += len(re.findall(r"(?<!\+)\+(?!\+)|\*|/|%|\bmin\b|\bmax\b", body))
        c["binary -"] += len(re.findall(r"[\w)]\s+-\s", body))
        c["unary - / !"] += len(re.findall(r"!(?!=)|(?:[(\[,]|then|else|:=)\s*-(?!\d)", body))
    return c


def main():
    t_start = time.time()
    lean_path = sys.argv[1]
    quick = "--quick" in sys.argv[2:]
    src, ms, defs, order, sites, byid, muts, ill = M.setup(lean_path)
    print("Referee cross-check of mechanical_mutations.py (REVIEW4-CT)")
    print(f"source sha256 pinned and checked by mechanical_mutations.setup: {M.PIN}")

    print("\n[1] Readout cone, independent regex closure over Part 1")
    items = part1_items(src)
    p1defs = [n for n, (kw, _) in items.items() if kw == "def"]
    cone = regex_cone(items, M.ROOTS)
    print(f"  Part 1 items: {len(items)} ({len(p1defs)} def, {len(items) - len(p1defs)} structure/abbrev)")
    print(f"  cone: {len(cone)} definitions; equal to the generator's 44: {cone == set(order)}")
    outside = sorted(set(p1defs) - cone)
    main_txt = src[src.index("theorem main :"):src.index(":= by", src.index("theorem main :"))]
    used_by_main = set()
    for n in outside:
        if re.search(r"\b" + re.escape(n.split(".")[-1]) + r"\b", main_txt):
            used_by_main.add(n)
    main_cone = regex_cone(items, sorted(used_by_main))
    print(f"  outside the cone ({len(outside)}): {', '.join(outside)}")
    print(f"    named in `main` or reached from those names: {', '.join(sorted(n for n in outside if n in main_cone))}")
    print(f"    in neither: {', '.join(sorted(n for n in outside if n not in main_cone))}")
    print(f"  any cone body mentions an outside definition: "
          f"{any(re.search(r'(?<![A-Za-z_])' + re.escape(o.split('.')[-1]) + r'\b', items[n][1]) for n in cone for o in outside if o.split('.')[-1] not in ('valid',))}")

    print("\n[2] Coarse site counts over the 44 cone bodies (regex), against the generator")
    gen = Counter(s.op for s in sites)
    cc = coarse_counts(items, cone)
    nneg = sum(1 for s in sites if s.op == "LIT" and s.info["value"] < 0)
    rows = [("IF", cc["IF"], gen["IF"]), ("MATCH", cc["MATCH"], gen["MATCH"]), ("LOG", cc["LOG"], gen["LOG"]),
            ("QNT", cc["QNT"], gen["QNT"]), ("BLIT", cc["BLIT"], gen["BLIT"]),
            ("LIT (integer literals not after '.')", cc["LIT"], gen["LIT"]),
            ("negative literals", cc["negative LIT"], nneg),
            ("REL", cc["REL"], gen["REL"]),
            ("AOR (+ * / % min max, binary -)", cc["AOR (+ * / % min max)"] + cc["binary -"], gen["AOR"]),
            ("NEG (unary -, ! on a non-literal)", cc["unary - / !"], gen["NEG"]),
            ("list literals", cc["list literals"], gen["LIST"])]
    for lab, a, b in rows:
        print(f"  {lab:40s} regex {a:4d}   generator {b:4d}   {'equal' if a == b else 'DIFFER'}")
    print(f"  `++` tokens {cc['++']} in {gen['CAT']} chains; `::` tokens outside patterns {cc['::']} in {gen['CONS']} chains "
          f"(generator)")

    print("\n[3] Example mutants printed in out section 4, re-evaluated")
    idx = {m.id: i for i, m in enumerate(muts)}
    out = open(os.path.join(HERE, "mechanical_mutations.out")).read()
    ok = bad = 0

    def ev(mid):
        hs = {}
        for part in mid.split(" & "):
            hs.update(muts[idx[part]].hooks)
        if any(ill[idx[p]] for p in mid.split(" & ")):
            return ("ill-typed",)
        return M.canonical(T.evaluate(hs))
    for m in re.finditer(r"^\s+([\d.e+]+) sigma\s+1/value\s+([\d.]+)\s+(\(.*?\)/\d+)\s+(\d+) mutant\(s\), e\.g\. (.+)$",
                         out, re.M):
        dev, invp, form, cnt, ex = m.groups()
        k = ev(ex.strip())
        inv = inv_exact(k)
        devm = abs(inv - INV_CODATA) / SIG
        good = (M.fmt_key(k) == form and abs(inv - mp.mpf(invp)) < mp.mpf("1e-9") * 1000 and
                abs(devm - mp.mpf(dev)) / mp.mpf(dev) < mp.mpf("1e-4"))
        ok, bad = ok + good, bad + (not good)
        print(f"  {'ok ' if good else 'BAD'} {ex.strip():50s} {M.fmt_key(k):38s} 1/value {mp.nstr(inv, 15)}  "
              f"{mp.nstr(devm, 6)} sigma")
    for m in re.finditer(r"^\s{4}(\S.*?)\s+1/value\s+([\d.]+)\s+in B: (?:yes|no )\s+mechanical: (single|pair) "
                         r"\(e\.g\. (.+)\)$", out, re.M):
        name, invp, kind, ex = m.groups()
        k = ev(ex.strip())
        inv = inv_exact(k)
        nd = len(invp.split(".")[1])
        good = (mp.nstr(inv, 30)[:len(invp)] == invp or abs(inv - mp.mpf(invp)) <= mp.mpf(10) ** (-nd)) and \
            (kind == "single") == (" & " not in ex)
        ok, bad = ok + good, bad + (not good)
        print(f"  {'ok ' if good else 'BAD'} {name.strip()[:40]:40s} {ex.strip():48s} 1/value {mp.nstr(inv, 12)} "
              f"(printed {invp})")
    print(f"  {ok} agree, {bad} disagree")

    if quick:
        return
    print("\n[4] Seeded uniform sample of 2000 of the 332 308 mutants (seed 20261030)")
    pairs = M.all_pairs(muts)
    n1 = len(muts)
    space = n1 + len(pairs)
    rng = random.Random(SEED)
    pick = rng.sample(range(space), 2000)
    cat = Counter()
    for x in pick:
        if x < n1:
            mids = [x]
        else:
            mids = list(pairs[x - n1])
        if any(ill[i] for i in mids):
            cat["fail to elaborate"] += 1
            continue
        hs = {}
        for i in mids:
            hs.update(muts[i].hooks)
        k = M.canonical(T.evaluate(hs))
        lab = {"INFEASIBLE": "infeasible", "N=0": "dimension 0", "zero": "numerator 0"}.get(k[0], "well-defined") \
            if len(k) == 1 else "well-defined"
        cat[lab] += 1
        if k == M.BASE_KEY:
            cat["  of which the Letter's value"] += 1
    full = {"fail to elaborate": 47 + 37105, "infeasible": 1 + 715, "dimension 0": 24 + 17924, "numerator 0": 4 + 3500,
            "well-defined": 742 + 272246, "  of which the Letter's value": 29292}
    for lab, tot in full.items():
        p = tot / space
        e = 2000 * p
        sd = (2000 * p * (1 - p)) ** 0.5
        print(f"  {lab:30s} sample {cat[lab]:5d}   expected from the .out {e:8.1f} +- {sd:5.1f}   "
              f"z = {(cat[lab] - e) / sd:+.2f}")

    print("\n[5] The Born-weight value 319/31104 (Table S5 'Born weights k_r')")
    born = M.canonical(([(1, 319)], 31104))
    hits = set()
    for anchor in ("recs:REL:0=neg", "sameForm:REL:0=neg"):
        a = idx[anchor]
        for b in range(n1):
            if b == a or ill[b] or not M.compatible(muts[a], muts[b]):
                continue
            hs = dict(muts[a].hooks)
            hs.update(muts[b].hooks)
            if M.canonical(T.evaluate(hs)) == born:
                hits.add(" & ".join(sorted((muts[a].id, muts[b].id))))
    print(f"  pairs through either anchor that give 319/31104: {len(hits)}")
    for h in sorted(hits):
        print(f"    {h}")
    print(f"  singles: recs:REL:0=neg -> {M.fmt_key(M.canonical(T.evaluate(dict(muts[idx['recs:REL:0=neg']].hooks))))}; "
          f"sameForm:REL:0=neg -> {M.fmt_key(M.canonical(T.evaluate(dict(muts[idx['sameForm:REL:0=neg']].hooks))))}")

    print("\n[6] Depth dependence: seeded samples of 12 000 pairs and 12 000 triples (pairwise-disjoint sites)")
    lean_triples = []
    for depth, seed in ((2, SEED + 2), (3, SEED + 3)):
        t0 = time.time()
        r = random.Random(seed)
        keys, combos = [], []
        while len(keys) < 12000:
            c = sorted(r.sample(range(n1), depth))
            if not all(M.compatible(muts[x], muts[y]) for i, x in enumerate(c) for y in c[i + 1:]):
                continue
            combos.append(c)
            if any(ill[x] for x in c):
                keys.append(("ill-typed",))
                continue
            hs = {}
            for x in c:
                hs.update(muts[x].hooks)
            keys.append(M.canonical(T.evaluate(hs)))
        good = [k for k in keys if len(k) == 2]
        inband = Counter(k for k in good if 120 < inv_exact(k) < 155)
        f = Counter(inband.values())
        S, f1, f2 = len(inband), f[1], f[2]
        chao = S + (f1 * f1 / (2 * f2) if f2 else f1 * (f1 - 1) / 2)
        near = sorted((abs(inv_exact(k) - INV_CODATA) / SIG, k) for k in set(good) if k != M.BASE_KEY)[:1]
        print(f"  depth {depth}: {len(keys)} sampled, {len(good)} well-defined, {sum(1 for k in good if k == M.BASE_KEY)} "
              f"with the Letter's value, {len(set(good))} distinct values, {S} distinct with 1/value in (120,155) "
              f"(seen once {f1}, twice {f2}); Chao1 estimate of the in-band values of the whole depth-{depth} space: "
              f"{chao:.0f} (density {chao / 35:.1f} per unit); closest other value {mp.nstr(near[0][0], 4)} sigma "
              f"[{time.time() - t0:.0f} s]")
        if depth == 3:
            for c, k in zip(combos, keys):
                if len(k) == 2 and k != M.BASE_KEY and 120 < inv_exact(k) < 155 and len(lean_triples) < 4:
                    lean_triples.append(c)
    print("  reference: the full depth-2 space has 1105 in-band values (out section 4). The Chao1 estimate from the "
          "pair sample above is to be read against that number.")
    tri = sum(1 for a in range(n1) for b in range(a + 1, n1) if M.compatible(muts[a], muts[b]))
    print(f"  space sizes: singles {n1}, pairs with disjoint sites {tri}; C(818,3) = {n1 * (n1 - 1) * (n1 - 2) // 6}, "
          f"C(818,4) = {n1 * (n1 - 1) * (n1 - 2) * (n1 - 3) // 24}, C(818,5) = "
          f"{n1 * (n1 - 1) * (n1 - 2) * (n1 - 3) * (n1 - 4) // 120} (before the disjointness rule)")

    print("\n[7] Lean 4.35.0 on mutants outside the spot-check sample")
    lean = os.environ.get("LEAN")
    if not lean:
        print("  skipped (LEAN not set)")
    else:
        sample = {pairs[p] for p in M.pair_sample(muts, pairs)}
        chosen = ["recs:REL:1=strict & linkSwap:VAR:0=read2", "recs:REL:0=neg & sameForm:REL:0=neg",
                  "readoutDim:AOR:5=/ & objectWorld:LIT:2=2", "wX:VAR:0=farEnd & wY:AOR:0=+",
                  "commutatorRecs:VAR:16=dA & slotQudit:AOR:0=+", "recs:REL:0=neg & addRoot:IF:0=swap"]
        combos = [sorted(idx[p] for p in c.split(" & ")) for c in chosen] + lean_triples
        with tempfile.TemporaryDirectory() as d:
            agree = 0
            for j, c in enumerate(combos):
                assert len(c) != 2 or tuple(c) not in sample
                f = os.path.join(d, f"m{j}.lean")
                with open(f, "w") as g:
                    g.write(M.lean_text(src, sum((muts[x].edits for x in c), [])))
                res = subprocess.run([lean, f], capture_output=True, text=True, timeout=300)
                lines = res.stdout.split("\n")          # skip linter warnings, as mutation_spotcheck --compare does
                st = next((x for x, l in enumerate(lines) if l.startswith("([")), None)
                got = re.sub(r"\s+", "", "".join(lines[st:])) if st is not None else "?"
                warn = sum(1 for l in lines if ": warning:" in l)
                hs = {}
                for x in c:
                    hs.update(muts[x].hooks)
                exp = re.sub(r"\s+", "", repr(T.evaluate(hs)))
                agree += got == exp and res.returncode == 0 and "error" not in res.stdout
                print(f"  {'agree' if got == exp else 'DIFFER'}  {' & '.join(muts[x].id for x in c):75s} "
                      f"lean {got[:48]}{'...' if len(got) > 48 else ''}" + (f" ({warn} linter warning)" if warn else ""))
            print(f"  {agree} of {len(combos)} agree (in the 300-pair sample: none of the 6 pairs)")
    print("\n[8] Pairs excluded by the disjoint-site rule, by reason")
    why = Counter()
    for a in range(n1):
        for b in range(a + 1, n1):
            A, B = muts[a], muts[b]
            if M.compatible(A, B):
                continue
            if set(A.hooks) & set(B.hooks):
                why["the two share a site"] += 1
                continue

            def reason(X, Y):
                if X.op == "IF":
                    st = byid[next(iter(X.hooks))]
                    t0 = st.body.toks[st.info["then"][0]].s
                    return ("IF swap x mutation in its own condition" if all(e <= t0 for _, e in Y.spans)
                            else "IF swap x mutation in a branch")
                if X.op == "LOG" and X.id.endswith(("dropL", "dropR")):
                    st = byid[next(iter(X.hooks))]
                    kp = st.info["R"] if X.id.endswith("dropL") else st.info["L"]
                    ks, ke = st.body.toks[kp[0]].s, st.body.toks[kp[1] - 1].e
                    return ("LOG drop x mutation in the operand it keeps" if all(ks <= x and e <= ke for x, e in Y.spans)
                            else "LOG drop x mutation in the operand it drops")
                if X.op in ("LIST", "CAT", "CONS"):
                    return f"{X.op} x mutation inside it"
                return None
            why[reason(A, B) or reason(B, A) or "other overlap"] += 1
    for k, v in why.most_common():
        print(f"  {v:5d}  {k}")
    print(f"  total {sum(why.values())} = {n1 * (n1 - 1) // 2} - {len(pairs)}")
    print(f"\ntotal time {time.time() - t_start:.0f} s")


if __name__ == "__main__":
    main()
