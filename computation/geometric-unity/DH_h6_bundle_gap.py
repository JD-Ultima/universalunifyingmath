"""Exact trace, weight and inducing-parameter calibrations for DH."""
from pathlib import Path
from collections import Counter
import argparse,json
import sympy as s
weights=Counter()
# Sym^10 of the four-dimensional Spin5 module has highest55 and
# dimension286, hence equals V55 by the Weyl dimension formula.
for a in range(11):
    for b in range(11-a):
        for c in range(11-a-b):
            d=10-a-b-c;weights[(a+b-c-d)//2,(a-b+c-d)//2]+=1
assert sum(weights.values())==286 and weights[5,5]==1
assert all(abs(a)<=5 and abs(b)<=5 for a,b in weights)
rho=[2,1,0,-1,-2];assert sum(x*x for x in rho)==10
H=s.diag(1,-1);assert s.trace(H*H)==2
records=[]
for cuts in range(1,16):
    bb=[0]+[i+1 for i in range(4) if cuts>>i&1]+[5];blocks=[b-a for a,b in zip(bb,bb[1:])]
    if any(n>=3 for n in blocks):
        assert any(n-1>n//2 for n in blocks if n>=3)
        continue
    nu=s.symbols('v0:'+str(len(blocks)-1),real=True)
    last=-sum(n*v for n,v in zip(blocks[:-1],nu))/blocks[-1];vs=[*nu,last]
    ell=s.symbols('l0:'+str(blocks.count(2)),real=True);j=0;parameter=[]
    for n,v in zip(blocks,vs):
        if n==2:parameter.extend([ell[j]/2+s.I*v,-ell[j]/2+s.I*v]);j+=1
        else:parameter.append(s.I*v)
    omega=s.expand(sum(x*x for x in parameter)-10)
    expected=sum(x*x for x in ell)/2-10-sum(n*v*v for n,v in zip(blocks,vs))
    assert s.expand(omega-expected)==0
    bound=s.simplify(omega.subs({**{v:0 for v in nu},**{x:4 for x in ell}}))
    assert bound<=6
    records.append(dict(blocks=blocks,split_dimension=len(nu),number_SL2_blocks=len(ell),casimir_upper_bound=int(bound)))
assert len(records)==8 and max(x['casimir_upper_bound'] for x in records)==6
assert s.Rational(35-6,4)==s.Rational(29,4) and s.Rational(29,4)-7==s.Rational(1,4)
result=dict(weight_dimension=286,SO2_frequency_bound=5,SL2_lowest_K_weight_bound=5,SL2_Harish_Chandra_parameter_bound=4,trace_rho_squared=10,tempered_Casimir_upper_bound=6,rough_Laplacian_lower_bound='29/4',closed_field_eigenvalue=7,spectral_gap='1/4',families=records,scope='Exact finite calibrations; transfer to the full L2 bundle spectrum uses the representation-theoretic proof in DH_h6_bundle_gap_proof.txt.')
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args();receipt=Path(__file__).with_suffix('.json')
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:assert json.loads(receipt.read_text())==result
print(json.dumps(result,indent=2))
