"""Evaluate 1/value for every readout printed by the author's mutations.py
(run with Lean 4.35.0 on the pinned upstream file; log: upstream_mutations_lean.out).
Each readout is ([(radicand, coefficient), ...], N).  Run: python3 mutations_values.py"""
import ast, mpmath as mp
mp.mp.dps = 30
INV = mp.mpf('137.035999177'); SIG = mp.mpf('2.1e-8')
for line in open('upstream_mutations_lean.out'):
    name = line.split()[0]
    lst, N = ast.literal_eval(line[line.index('readout') + 8:].strip())
    num = sum(c * mp.sqrt(r) for r, c in lst)
    inv = N / num
    print(f"{name:34s} 1/value = {mp.nstr(inv, 12):>16s}   ({mp.nstr((inv - INV) / SIG, 4)} sigma)")
