"""Independent enumeration of the configurations of Supplement Sec. S5
(Theorem S-config) and of the uniqueness criteria of Theorem 1 under the
charge conventions of F-CT-3.

Configuration: centre c in {2,3}; multiset of arms (far dim e in {2,3}, R/U);
links between 2 and c+1; at least one read.  Read cap: read arm needs c*e <= 6.
'Readable': >= 2 reads, all read arms readable.  (The record-space condition of
the paper is never binding here; the table S2 count is reproduced without it.)
Uniqueness criteria (Sec. S5): (a) end rule: a read arm is read at the larger
end, i.e. needs c >= e; (b) two-point: exactly three arms of one type, two read;
(c) acting charge on every link, under reading A (c != e) or A* (always).
Run: python3 config_enumeration.py
"""
from itertools import combinations_with_replacement as cwr
ARMS = [(2, 'R'), (2, 'U'), (3, 'R'), (3, 'U')]
configs = []
for c in (2, 3):
    for L in range(2, c + 2):
        for arms in cwr(ARMS, L):
            if any(a[1] == 'R' for a in arms):
                configs.append((c, arms))
admitted = [k for k in configs if all(k[0] * e <= 6 for e, t in k[1] if t == 'R')]
readable = [k for k in admitted if sum(t == 'R' for _, t in k[1]) >= 2]
print("configurations:", len(configs), " by centre:", {c: sum(k[0] == c for k in configs) for c in (2, 3)})
print("admitted by read cap:", len(admitted), "  readable:", len(readable))
assert (len(configs), len(admitted), len(readable)) == (76, 42, 23)

def end_rule(k): return all(k[0] >= e for e, t in k[1] if t == 'R')
def two_point(k):
    arms = k[1]
    return len(arms) == 3 and sum(t == 'R' for _, t in arms) == 2 and len({e for e, _ in arms}) == 1
for name, acting in [("A  (Letter)", lambda c, e: c != e), ("A* (conjugate)", lambda c, e: True),
                     ("B  (matched phase)", lambda c, e: False)]:
    sel = [k for k in readable if end_rule(k) and two_point(k) and all(acting(k[0], e) for e, _ in k[1])]
    fmt = lambda k: f"centre {k[0]}: " + " ".join(f"{e}{t}" for e, t in k[1])
    print(f"reading {name}: {len(sel)} configuration(s) satisfy (a)+(b)+(c):", [fmt(k) for k in sel])
