"""Counts quoted in reports/geometric-unity-phase1.md, computed from the register files.
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
for line in rows:
    cells = [c.strip() for c in line.replace('\\|', '¦').strip().strip('|').split('|')]
    assert len(cells) == 6, cells[0]
    first = re.split(r'[ /+(]', cells[2])[0]
    primary[first] += 1
    s = re.match(r'\**([VFCTO—])', cells[5]).group(1)
    status[s] += 1
    if s == 'F':
        falses.append(cells[0])
gaps = [l for l in (root / 'gaps/geometric-unity.md').read_text().splitlines() if l.startswith('| GU-G')]
print(f"claims: {len(rows)}")
print("primary type (first type listed):", dict(primary.most_common()))
print("status:", {names[k]: v for k, v in status.most_common()})
print("false as stated:", falses)
print(f"gaps: {len(gaps)}")
