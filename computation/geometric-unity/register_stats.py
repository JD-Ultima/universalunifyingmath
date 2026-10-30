"""Counts quoted in reports/geometric-unity-phase1.md, computed from the register files.
A re-graded claim is counted by its current status, the last bolded status letter in its status cell.
Run from this directory: python3 register_stats.py > register_stats.out"""
import collections
import re
from pathlib import Path

root = Path(__file__).resolve().parents[2]
rows = [l for l in (root / 'claims/geometric-unity.md').read_text().splitlines() if l.startswith('| GU-0')]
primary = collections.Counter()
status = collections.Counter()
names = {'V': 'verified', 'C': 'conditional/repaired', 'O': 'open/not well-posed', 'F': 'false as stated',
         'T': 'trivially true', '—': 'not mathematical'}
falses = []
splits = []
for line in rows:
    cells = [c.strip() for c in line.replace('\\|', '¦').strip().strip('|').split('|')]
    assert len(cells) == 6, cells[0]
    first = re.split(r'[ /+(]', cells[2])[0]
    primary[first] += 1
    # current status: the last bolded status letter (re-graded rows read "O ... -> **C** ..."), else the leading one
    bold = re.findall(r'\*\*([VFCTO—])\*\*', cells[5])
    s = bold[-1] if bold else re.match(r'\**([VFCTO—])', cells[5]).group(1)
    if len(set(bold)) > 1:
        splits.append(f"{cells[0]} ({'/'.join(dict.fromkeys(bold))}, counted as {s})")
    status[s] += 1
    if s == 'F':
        falses.append(cells[0])
gaps = [l for l in (root / 'gaps/geometric-unity.md').read_text().splitlines() if l.startswith('| GU-G')]
print(f"claims: {len(rows)}")
print("primary type (first type listed):", dict(primary.most_common()))
print("status:", {names[k]: v for k, v in status.most_common()})
print("false as stated:", falses)
print("re-graded with more than one status (counted by the last):", splits)
print(f"gaps: {len(gaps)}")
