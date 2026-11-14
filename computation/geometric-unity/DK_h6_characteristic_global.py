"""Exact calibrations for the globally regular balanced characteristic lift."""
import argparse
import json
from pathlib import Path
import time

import sympy as s

from DK_h6_characteristic_field import verify as verify_local_field


HERE = Path(__file__).resolve().parent
parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--write-receipt", action="store_true")
args = parser.parse_args()
local_result, p, metric, v, leaf, stabilizer = verify_local_field()
tick = time.time()
order = [1, 3, 4, 0, 2]
G=metric.extract(order,order);Gi=G.inv();Q=s.Matrix([[0,1],[1,0]])
# Explicit Lie algebra description of the new leaf subgroup.
for X in leaf + stabilizer:
    X=X.extract(order,order)
    assert X[2:,:2]==s.zeros(3,2) and X[3:,:3]==s.zeros(2,3)
    assert X[3:,3:]==Q*X[:2,:2]*Q
assert len(leaf) + len(stabilizer) == 12
assert local_result["origin_entries"] == [[1, 4875, "1"], [1, 8971, "-1"]]
# Arbitrary upper block group matrix in this subgroup.
a,b,c,d,ell0,ell1,n00,n01,n10,n11=s.symbols('a b c d ell0 ell1 n00 n01 n10 n11')
L=s.Matrix([[a,b],[c,d]]);D=Q*L*Q;middle=1/L.det()**2
h=s.zeros(5);h[:2,:2]=L;h[:2,2]=s.Matrix([ell0,ell1]);h[:2,3:]=s.Matrix([[n00,n01],[n10,n11]]);h[2,2]=middle;h[3:,3:]=D
assert s.simplify(h.det())==1
y=h*Gi*h.T;R=y[:2,3:]
assert s.simplify(R*Q-Q*R.T)==s.zeros(2)
Cv=v.extract(order,order)*Gi
assert Cv[:2,:2]==Q/8 and Cv[2:,:]==s.zeros(3,5) and Cv[:,2:]==s.zeros(5,3)
assert s.simplify((h*Cv*h.T)[:2,:2]-R*Q/4)==s.zeros(2)
# Branch cancellation: r=sqrt(z)=exp(u), du=dr/r.
r,rr12,rr21,rr22=s.symbols('r rr12 rr21 rr22',nonzero=True)
R=s.Matrix([[r*r*rr22,rr12],[rr21,rr22]]);w=r*r/(1+r*r)
H=s.Matrix([[rr12*w,R[0,0]*(1-w)],[R[0,0]*(1-w),rr21*(1-w)]])
coefficient=s.Matrix([[rr12*r,R[0,0]/r],[R[0,0]/r,rr21/r]])/4
multiplier=r**3/(1+r*r)**3
assert s.simplify(multiplier*coefficient/r-H*s.diff(w,r)/8)==s.zeros(2)
# Base differential normalization in the exact BO frame.
dt=[];dw=[]
for P in p:
    dy=P.extract(order,order)*Gi;F=dy[3:,3:];e=dy[2,3:];C=dy[:2,3:]
    dt.append([F[0,0],F[0,1],e[0,0],F[1,1],e[0,1]])
    dw.append((C[0,0]-C[1,1])/2)
assert s.Matrix(dt).T==s.eye(14).extract([0,1,3,8,9],list(range(14)))
assert dw==[0]*12+[s.Rational(1,2),s.Rational(-1,2)]
Hbase=Q/4
assert Hbase/2==Cv[:2,:2]
result=dict(global_denominator='tr(R), strictly positive on real Y by CU',extra_coordinate='w=R11/tr(R)',regular_coefficient='n_W H n_W^T',H='[[R12*w,R11*(1-w)],[R11*(1-w),R21*(1-w)]]',local_multiplier_before_normalizing='exp(3u)/(1+exp(2u))^3',normalizing_factor=8,base_dw='(e12-e13)/2',origin_matches_DK_local=True,global_holomorphic_branch_required=False,scope='Exact group/branch-cancellation/base calibrations; global real analyticity and closure follow by CU spectral regularity and the analytic identity theorem as detailed in the proof.',seconds=time.time()-tick)
receipt=HERE/'DK_h6_characteristic_global.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in known.items() if k!='seconds'}=={k:v for k,v in result.items() if k!='seconds'}
print(json.dumps(result,indent=2))
