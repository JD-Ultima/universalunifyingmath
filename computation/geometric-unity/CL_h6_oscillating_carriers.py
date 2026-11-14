"""Verify the CL H6 Cartan carrier fields and 52-dimensional carrier space."""

import argparse
from fractions import Fraction as F
from itertools import combinations
import json
from math import prod
from pathlib import Path
from time import perf_counter

import numpy as np

HERE = Path(__file__).resolve().parent

import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BO_h6_full_kernel as BO
import BR_h6_an_closed as BR

AXIS = 13
G = tuple(BC.METRIC)
V = range(14)
EVEN = tuple(i for i in range(32) if i.bit_count() % 2 == 0)


def metric_dot(a,b):
    return sum((v*b.get((i,m),0)/(G[i]*prod(G[j] for j in BC.bit_indices(m))) for (i,m),v in a.items()),F(0))


def verify_character_differential():
    # The SPD tangent is Phi=5Phi/5, and the upper triangular infinitesimal
    # lift A satisfies A+A^T=Phi when using the congruence realization.
    Q = BG.integer_phi()
    lifts = tuple(tuple(tuple(F(m[i][j],5) if i<j else F(m[i][i],10) if i==j else F(0)
                               for j in range(5)) for i in range(5)) for m in Q)
    omega = [F(-1,4)] * 4 + [F(1)]
    assert sum(omega) == 0
    derivative = [sum(w * A[i][i] for i,w in enumerate(omega)) for A in lifts]
    assert derivative == [F(0)] * AXIS + [F(1)]
    assert sum(value * value / gram for value, gram in zip(derivative, G)) == F(5, 16)
    assert all(not any(BG.bracket(A,B)[i][i] for i in range(5)) for A,B in combinations(lifts,2))
    print("Cartan additive character chi(s)=-(1/4)sum_{j<4} log(s_jj)+log(s_44) has dchi=e^13.")


def verify_section(entry, R):
    ch = entry["character"]
    data = entry["entries"]
    b = {(a,m):F(c) for a,m,c in data}
    assert b and len(b)==len(data)
    assert all(a in V and m.bit_count()==6 and m < 1<<14 and m>>AXIS & 1 and c and BR.chars[a]^BR.char(m)==ch
               for (a,m),c in b.items())
    assert BC.hook_traces(b)==({}, {}) and not BR.check_sparse(b)
    assert not BG.full_sparse_curvature(R,b)
    return b


def verify_rank(ch, expected):
    columns, rows = BR.stack(ch)
    indexes = [j for j,(a,m) in enumerate(columns) if m>>AXIS & 1]
    lookup = {a:j for j,a in enumerate(indexes)}
    local = [{lookup[a]:v for a,v in row.items() if a in lookup} for row in rows]
    local = [z for z in local if z]
    A=np.zeros((len(local),len(indexes)),dtype=np.int64)
    for j,row in enumerate(local):
        for i,c in row.items():
            A[j,i]=c%BO.PRIME
    rank=BO.modular_rank(A)
    rank_exact,basis=BR.exact_kernel(local,len(indexes))
    count=len(basis)
    assert rank==rank_exact==len(indexes)-count
    assert [len(indexes),len(local),rank,count]==expected
    # BR's kernel routine checks every exact candidate by multiplication.
    print(f"Sign character {ch:02d}: {len(indexes)} columns, rank {rank}, exact nullity {count}.",flush=True)
    return count


def run(full):
    t=perf_counter()
    verify_character_differential()
    data=json.loads((HERE/"CL_h6_carrier_witness.json").read_text(encoding="utf-8"))
    assert data["axis"]==AXIS and data["prime"]==BO.PRIME
    R=BG.curvature_model()
    m=[verify_section(x,R) for x in data["fields"]]
    assert len(m)==2
    gram=[[metric_dot(a,b) for b in m] for a in m]
    assert gram==[[F(1026505,524288),F(0)],[F(0),F(1710915,32768)]]
    assert gram==[list(map(F,row)) for row in data["gram"]]
    assert [len(z) for z in m]==[222,225]
    print("Two exact AN carrier fields pass full hooks, covariant exterior closure, and coefficient curvature. Gram:",gram)
    expected={int(ch):moves for ch,moves in data["table"].items()}
    assert tuple(sorted(expected))==EVEN
    if full:
        dims={ch:verify_rank(ch,expected[ch]) for ch in EVEN}
        assert sorted(dims.values())==[2]+[3]*10+[4]*5 and sum(dims.values())==52
        print("All 16 ranks give exactly 52 real carrier directions.")
    else:
        assert verify_rank(0,expected[0])==2
    print("Unit pair yields all-frequency bounded constant-norm closed profiles. Runtime",round(perf_counter()-t,3),"seconds.")


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--full",action="store_true",help="rebuild all 16 exact rank blocks")
    run(p.parse_args().full)
