"""Check exact failure of the BV receiver at one genuine imaginary AN character."""

from fractions import Fraction as F
import gzip
import hashlib
import json
from math import isqrt, prod
from pathlib import Path
from time import perf_counter

import sympy as s

HERE = Path(__file__).resolve().parent

import BC_affine_hooks as BC
import BC_quadratic_xcheck as X
import BG_h6_curvature as BG
import BN_h6_invariant_grade9 as BN
import BR_h6_an_closed as BR

G = tuple(BC.METRIC)
OMEGA = 4
AXIS = 13
N = len(G)
FULL = (1 << N) - 1


def metric(a, m):
    return G[a] * prod(G[j] for j in BC.bit_indices(m))


def input_tensor(rows, expected=2638):
    t = {(a,m):F(v) for a,m,v in rows}
    assert len(t)==len(rows)==expected
    assert all(0<=a<N and m.bit_count()==6 and m<=FULL and v and BR.chars[a]^BR.char(m)==0 for (a,m),v in t.items())
    return t


def verify_geometry(t):
    R = BG.curvature_model()
    for z in t:
        assert BC.hook_traces(z)==({}, {}) and not BG.full_sparse_curvature(R,z)
    Ds=list(map(BR.check_sparse,t))
    for source,target,sign in ((1,0,-1),(0,1,1)):
        for (a,m),c in t[source].items():
            w=BR.wedge(AXIS,m)
            if w:
                BR.add(Ds[target],(a,w[0]),sign*OMEGA*BR.den*w[1]*c)
    assert not any(Ds)
    print("Both rational real components pass hooks, curvature, and the exact twisted closure equation.")


def square_root_classes(g):
    out = {}
    for q in g:
        magnitude=q.numerator*q.denominator
        powers=s.factorint(magnitude)
        base=prod(int(k) for k,v in powers.items() if v%2)
        factor=F(prod(int(k)**(v//2) for k,v in powers.items()), q.denominator)
        assert factor**2*base==q
        out[q]=(base,factor)
    return out


def evaluate(t,receiver):
    den = receiver["denominator"]
    roots=square_root_classes({metric(a,m) for a,m,c in receiver["rows"]})
    algebra = X.BladeAlgebra(G)
    measure={}
    denom = F(0)
    for z in t:
        norm=sum((x*x/metric(a,m) for (a,m),x in z.items()),F(0))
        denom += norm
        field=[{} for _ in G]
        for (a,m),x in z.items():
            h,scale=algebra.multiply(m,FULL)
            BC.plus(field[a],h,x*scale)
        assert algebra.multiply(FULL,FULL)==(0,-1/prod(G))
        q=BN.PhysicalQuadratic(field)
        for j,(a,m,v) in enumerate(receiver["rows"]):
            original=q.coefficient(a,m)
            if j in (0,10,30,1000,2000,3000):
                assert original==prod(G)*algebra.quadratic_coefficient(field,a,m)
            physical=original+int(m==1<<a)*F(26,21)*norm*G[a]
            r,scale=roots[metric(a,m)]
            measure[r]=measure.get(r,F(0))+physical*F(v,den)/(r*scale)
    measure={a:v for a,v in measure.items() if v}
    return measure, denom


def enclose(coeff, norm):
    u=10**32
    center=uncertainty=F(0)
    r=F(1,u)
    for n,c in coeff.items():
        b=F(isqrt(n*u*u),u)
        assert b*b<=n<(b+r)**2
        center+=b*c
        uncertainty+=abs(c)*r
    return (center-uncertainty)/norm,(center+uncertainty)/norm


def main():
    t=perf_counter()
    raw=gzip.decompress((HERE/"CO_h6_twist_witness.json.gz").read_bytes())
    record=json.loads((HERE/"CO_h6_twist_certificate.json").read_text(encoding="utf-8"))
    assert hashlib.sha256(raw).hexdigest()==record["witness_sha256"]
    data=json.loads(raw)
    assert data["frequency"]==OMEGA and data["character"]==0
    negative=[input_tensor(data[k]) for k in ("real","imaginary")]
    positive=[input_tensor(data["positive"],762),{}]
    verify_geometry(negative)
    verify_geometry(positive)
    assert all(m >> AXIS & 1 for _, m in positive[0])
    assert F(OMEGA**2, G[AXIS]) == 5
    source=Path(BR.__file__).with_name("BV_receiver.json")
    receiver_bytes=source.read_bytes()
    assert hashlib.sha256(receiver_bytes).hexdigest()==record["BV_file_sha256"]
    receiver=json.loads(receiver_bytes)
    coeff,norm=evaluate(negative,receiver)
    assert coeff=={int(k):F(v) for k,v in record["radical_coefficients"].items()} and norm>0
    low,high=enclose(coeff,norm)
    assert (low,high)==(F(record["normalized_lower"]),F(record["normalized_upper"]))
    assert F(-314427,10**6)<low<=high<F(-314426,10**6)
    print("Exact Hermitian BV moment per squared norm is enclosed between -0.314427 and -0.314426.")
    coeff,norm=evaluate(positive,receiver)
    low,high=enclose(coeff,norm)
    assert F(0)<low<=high
    print(f"Independent positive witness: certified normalized interval ({float(low):.12g}, {float(high):.12g}).")
    print(f"Runtime {perf_counter()-t:.3f} seconds.")

if __name__=="__main__":
    main()
