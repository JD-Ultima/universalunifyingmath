"""Exact normalization and proper-parabolic checks for the DF spectral proof."""
from pathlib import Path
import argparse,json
from fractions import Fraction as F
import sympy as s
HERE=Path(__file__).resolve().parent
de=json.loads((HERE/'DE_h6_elliptic_periodic.json').read_text())
assert de['automatic_coclosure'] and de['weitzenbock_curvature_scalar']=='-7' and de['rough_laplacian_eigenvalue']=='7'
compact_casimir=5*(5+3)+5*(5+1);assert compact_casimir==70
group_casimir=F(compact_casimir,2)-4*F(de['rough_laplacian_eigenvalue']);assert group_casimir==7
# All proper real standard parabolics of split SL5 are ordered block
# decompositions. This includes any possible cuspidal Levi, without
# assuming connectedness of its M factor.
records=[]
for cuts in range(1,16):
    boundaries=[0]+[i+1 for i in range(4) if cuts>>i&1]+[5]
    blocks=[b-a for a,b in zip(boundaries,boundaries[1:])];n=len(blocks)-1
    gram=s.Matrix(n,n,lambda i,j:(blocks[i] if i==j else 0)+s.Rational(blocks[i]*blocks[j],blocks[-1]))
    minors=[gram[:j,:j].det() for j in range(1,n+1)]
    assert n>0 and all(x>0 for x in minors)
    records.append(dict(blocks=blocks,split_dimension=n,positive_gram_leading_minors=list(map(str,minors))))
assert len(records)==15
result=dict(compact_casimir=70,rough_laplacian_eigenvalue=7,trace_group_casimir=7,rank_complex_G=4,rank_K=2,proper_parabolic_checks=records,scope='Checks constants and positive continuous split parameter spaces. L2 nonexistence additionally invokes the full Harish-Chandra Plancherel and compact-Cartan discrete-series theorems, precisely as stated in DF_h6_l2_vanishing_proof.txt.')
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args();receipt=HERE/'DF_h6_l2_vanishing.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:assert json.loads(receipt.read_text())==result
print(json.dumps({k:v for k,v in result.items() if k!='proper_parabolic_checks'},indent=2))
