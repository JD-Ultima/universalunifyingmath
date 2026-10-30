"""Classification of the 19 Lean mutation readouts (audit F-CT-4 D; referee issue CT-2, second-round review 2026-10-30).

Input: `upstream_mutations_lean.out` (this directory), the log of the author's `mutations.py` run with Lean 4.35.0 on
the pinned ChainAlpha.lean. Each line gives the readout Lean evaluated, as a list of (k, count) with value
sum(count * sqrt(k)) / N.

For each mutation the script prints
  - the numerator, reduced exactly to a + b sqrt2 + c sqrt3 + d sqrt6 (integers), and N;
  - 1/value (mpmath, 30 digits; 12 shown);
  - whether the value equals a closed-form row of Table S5 (supplement.tex, Sec. "What else the count could have
    given"; rows transcribed below), and which;
  - whether the value lies in the declared rule space B of freedom.py (46 080 rules), with one witness rule
    (k_odd, k_mixed, k_QQ, k_qq; commutator factor, commutator count; amp/Born; N) whose numerator equals the Lean
    numerator exactly.
Equality is exact: 1, sqrt2, sqrt3, sqrt6 are linearly independent over Q, so two values are equal iff their
rational coefficient vectors (a, b, c, d) / N are equal.

A second table checks the 16 closed-form rows of Table S5 against their printed inverses (3 decimals), gives the
relative change against the computed row, of the value and of 1/value, and whether each row lies in B.
Run: python3 mutations_classified.py
"""
import ast
import itertools
import os
import re
from fractions import Fraction
import mpmath as mp

mp.mp.dps = 30
HERE = os.path.dirname(os.path.abspath(__file__))

# exact reduction of sqrt(k) to the basis (1, sqrt2, sqrt3, sqrt6)
SQRT = {1: (1, 0, 0, 0), 2: (0, 1, 0, 0), 3: (0, 0, 1, 0), 4: (2, 0, 0, 0), 6: (0, 0, 0, 1), 9: (3, 0, 0, 0),
        36: (6, 0, 0, 0)}
BORN = {k: (k, 0, 0, 0) for k in (1, 2, 3, 4, 6, 9, 36)}
BASIS = (mp.mpf(1), mp.sqrt(2), mp.sqrt(3), mp.sqrt(6))


def vadd(u, v):
    return tuple(a + b for a, b in zip(u, v))


def vscale(c, v):
    return tuple(c * a for a in v)


def key(vec, N):
    return tuple(Fraction(a, N) for a in vec)


def inv_value(vec, N):
    return mp.mpf(N) / sum(c * b for c, b in zip(vec, BASIS))


def fmt_num(vec):
    names = ("", "√2", "√3", "√6")
    parts = [f"{c}{n}" if n else f"{c}" for c, n in zip(vec, names) if c]
    return " + ".join(parts) if parts else "0"


# ---------------- declared rule space B (as in freedom.py, Part B) ----------------
CLASSES = (79, 48, 36, 6)          # odd; even mixed; even same-end QQ; even same-end qq
NS = {"24*36^2": 24 * 36 * 36, "36^3": 36 ** 3, "24^3": 24 ** 3, "24*35^2": 24 * 35 * 35, "35^3": 35 ** 3}
B = {}
nrules = 0
for ks in itertools.product([1, 2, 3, 6], repeat=4):
    for kc in [1, 2, 3, 4, 6, 9]:
        for nc in [0, 3, 6]:
            for power, w in (("amp", SQRT), ("born", BORN)):
                num = (0, 0, 0, 0)
                for n, k in zip(CLASSES, ks):
                    num = vadd(num, vscale(n, w[k]))
                num = vadd(num, vscale(nc, w[kc]))
                for nN, N in NS.items():
                    nrules += 1
                    B.setdefault(key(num, N), (ks, kc, nc, power, nN, num))
assert nrules == 46080

# ---------------- Table S5 closed-form rows (supplement.tex), with printed inverses ----------------
N0 = 24 * 36 * 36
S5 = [  # name, numerator (a, b, c, d), N, printed inverse
    ("computed", (79, 48, 42, 3), N0, 137.036),
    ("the two ends swapped", (79, 42, 48, 3), N0, 135.894),
    ("both exchanges at the qubit end", (85, 90, 0, 0), N0, 146.524),
    ("both exchanges at the qutrit end", (88, 0, 90, 0), N0, 127.536),
    ("exchange of reads where the letters sit", (85, 54, 36, 0), N0, 139.030),
    ("qubit-qubit records moved by Y", (79, 48, 36, 9), N0, 134.485),
    ("copies of the qutrit as separate qudits", (79, 0, 42, 51), N0, 112.423),
    ("survival without the record's own sign", (82, 48, 42, 0), N0, 139.713),
    ("commutator records survive X only", (79, 48, 45, 0), N0, 138.348),
    ("sums and differences of commutators kept", (79, 48, 42, 6), N0, 132.739),
    ("no commutator records", (79, 48, 42, 0), N0, 141.621),
    ("Born weights k_r", (319, 0, 0, 0), N0, 97.505),
    ("count of read ends", (158, 96, 82, 3), N0, 70.190),
    ("normalised by the 172 records", (79, 48, 42, 3), 172, 0.758),
    ("no weights", (172, 0, 0, 0), N0, 180.837),
    ("smallest readable configuration", (30, 30, 0, 0), 256, 3.535),
]
S5key = {key(v, N): name for name, v, N, _ in S5}
base_key = key(S5[0][1], S5[0][2])
v0 = 1 / inv_value(S5[0][1], S5[0][2])


def main():
    print(f"Declared space B: {nrules} rules, {len(B)} distinct values (freedom.py prints 15 639).")
    print("\nTable S5 closed-form rows: printed vs recomputed inverse; relative change of the value and of 1/value; in B?")
    n_alt = n_altB = 0
    changes = []
    for name, vec, N, printed in S5:
        iv = inv_value(vec, N)
        assert abs(float(iv) - printed) < 6e-4, (name, iv, printed)
        inB = key(vec, N) in B
        rel = abs((1 / iv) / v0 - 1)              # change of the value
        reli = abs(iv * v0 - 1)                   # change of 1/value
        if name != "computed":
            n_alt += 1
            n_altB += inB
            changes.append((reli, rel, name))
        print(f"  {name:42s} {'(' + fmt_num(vec) + ')/' + str(N):34s} 1/value {mp.nstr(iv, 9):>11s} (printed {printed:>7.3f})"
              f"  change: value {mp.nstr(100 * rel, 4):>9s}%, 1/value {mp.nstr(100 * reli, 4):>7s}%  in B: {'yes' if inB else 'no'}")
    changes.sort()
    print(f"  => {n_alt} closed-form alternatives; {n_altB} in B, {n_alt - n_altB} outside B; all printed inverses match.")
    print(f"  => inverses of the alternatives span [{mp.nstr(min(inv_value(v, N) for n, v, N, _ in S5[1:]), 5)}, "
          f"{mp.nstr(max(inv_value(v, N) for n, v, N, _ in S5[1:]), 5)}]")
    for reli, rel, name in changes[:2]:
        print(f"  => within 1%: {name}: 1/value {mp.nstr(100 * reli, 4)}%, value {mp.nstr(100 * rel, 4)}%")
    print(f"  => all other alternatives change both by more than 1%: "
          f"{all(min(a, b) > mp.mpf('0.01') for a, b, _ in changes[2:])}")

    pat = re.compile(r"^(\S+)\s+main\s+(proved|FAILS)\s+commutators\s+(proved|FAILS)\s+readout\s+\((\[.*\]),\s*(\d+)\)\s*$")
    rows = []
    with open(os.path.join(HERE, "upstream_mutations_lean.out")) as f:
        for line in f:
            m = pat.match(line.strip())
            if m and m.group(1) != "BASE":
                rows.append((m.group(1), ast.literal_eval(m.group(4)), int(m.group(5))))
    assert len(rows) == 19, len(rows)
    print("\nThe 19 Lean mutations (upstream_mutations_lean.out)")
    print(f"  {'mutation':36s} {'numerator / N':34s} {'1/value':>13s}  {'in Table S5?':42s} {'in B?':5s} witness in B / class")
    tally = {"unchanged": 0, "in B, Table S5 row": 0, "in B only": 0, "outside B": 0}
    for name, readout, N in rows:
        vec = (0, 0, 0, 0)
        for k, c in readout:
            vec = vadd(vec, vscale(c, SQRT[k]))
        kk = key(vec, N)
        s5 = S5key.get(kk, "-")
        wit = B.get(kk)
        if kk == base_key:
            cls = "unchanged"
        elif wit and s5 != "-":
            cls = "in B, Table S5 row"
        elif wit:
            cls = "in B only"
        else:
            cls = "outside B"
        tally[cls] += 1
        if wit:
            ks, kc, nc, power, nN, num = wit
            assert key(num, NS[nN]) == kk
            wtxt = f"k={ks} comm k={kc} x{nc} {power} N={nN}"
        else:
            wtxt = f"no rule of B ({sum(c for _, c in readout)} weighted records)"
        print(f"  {name:36s} {'(' + fmt_num(vec) + ')/' + str(N):34s} {mp.nstr(inv_value(vec, N), 12):>13s}  "
              f"{s5:42s} {'yes' if wit else 'no':5s} {wtxt} [{cls}]")
    print(f"  => {tally}, total {sum(tally.values())}")
    assert tally == {"unchanged": 2, "in B, Table S5 row": 10, "in B only": 3, "outside B": 4}
    print("  => 13 of 19 lie in B (10 of them Table S5 rows, 3 in B only); 4 lie outside B; 2 leave the value unchanged.")
    print("     For each of the 13, the Lean numerator equals the numerator of the witness rule of B exactly.")


if __name__ == "__main__":
    main()
