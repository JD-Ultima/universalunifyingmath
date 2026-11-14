"""Verify CS's exact high-pass threshold on the imaginary H6 Cartan line."""

from fractions import Fraction as F
import json
from math import isqrt, prod
from pathlib import Path
from time import perf_counter

HERE=Path(__file__).resolve().parent

import BC_affine_hooks as BC
import BG_h6_curvature as BG
import BR_h6_an_closed as BR

V=range(14)
G=tuple(BC.METRIC)
AXIS=13


def norm_data(project):
    # Row/column Schur bound in the unit physical Clifford coframe.
    blade=[prod(G[j] for j in BC.bit_indices(m)) for m in range(1<<14)]
    scale=10**6
    sqrt_upper={}
    d=[]
    for ch in range(32):
        if ch.bit_count()%2:
            continue
        cols,_=BR.stack(ch)
        receiver={}
        upper=F(0)
        for a,m in cols:
            v=F(0)
            for (b,k),c in BR.closure_column(a,m).items():
                if project and not k>>AXIS&1:
                    continue
                q=G[a]*blade[m]/(G[b]*blade[k])
                if q not in sqrt_upper:
                    h=isqrt(q.numerator*scale*scale//q.denominator)+1
                    square=F(h,scale)
                    assert square*square>q
                    sqrt_upper[q]=square
                x=abs(F(c,BR.den))*sqrt_upper[q]
                receiver[b,k]=receiver.get((b,k),F(0))+x
                v+=x
            upper=max(v,upper)
        R=max(receiver.values())
        d.append(dict(character=ch,cols=len(cols),rows=len(receiver),product=str(upper*R)))
    assert len(d)==16 and sum(x["cols"] for x in d)==42042
    return d


def main():
    t=perf_counter()
    Q=BG.integer_phi()
    g=[F(-1,4)]*4+[F(1)]
    tangent=[sum(F(q[i][i],10)*v for i,v in enumerate(g)) for q in Q]
    assert tangent==[F(0)]*13+[F(1)]
    assert 1/G[AXIS]==F(5,16)
    ref=json.loads((HERE/"CS_h6_highpass_certificate.json").read_text(encoding="utf-8"))
    data=norm_data(True)
    assert data==ref["projected_rows"]
    c=max(F(x["product"]) for x in data)
    assert c==F(9233556899517,160000000000)<58
    a=norm_data(False)
    assert a==ref["full_rows"]
    m=max(F(x["product"]) for x in a)
    assert m==F(ref["ambient_bound"])<64
    rho=F(3,7)-c/(F(5,16)*21**2)
    assert rho==F(ref["axial_margin"])>F(9,1000)
    # Fixed cone example: unit direction within 1/20 of e13 or its negative
    # and physical character norm >=14. Then δ+8/14 < sqrt(3/7).
    eta=F(1,20)+F(8,14)
    assert F(3,7)-eta*eta>0
    print("All 16 exact sign block Schur bounds checked in the unit metric.")
    print("Axis-13 projected closure norm squared <=",c,"< 58.")
    print("Ambient closure norm squared <=",m,"< 64.")
    print("For |omega|>=21, the fixed axial Hermitian stress exceeds",float(rho),"times squared norm.")
    print("Uniform cone of radius 1/20 around +/- unit e13 and |alpha|>=14 also has positive axial stress.")
    print("Runtime",round(perf_counter()-t,3),"seconds.")


if __name__=="__main__":
    main()
