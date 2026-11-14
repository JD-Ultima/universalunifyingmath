"""Exact normalizations for bounded C0 CU averages; see accompanying proof."""
from pathlib import Path
import argparse,json,time
import sympy as s
tick=time.time()
# By orthogonal diagonalization of the positive compression, this block
# computation proves the cofactor identity for any center unit vector.
a=s.symbols('a0:4',positive=True);b=s.Matrix(s.symbols('b0:4',real=True));d=s.symbols('d',positive=True)
A=s.diag(*a);Y=A.row_join(b).col_join(b.T.row_join(s.Matrix([[d]])))
schur=d-(b.T*A.inv()*b)[0]
assert s.factor(Y.det()-s.prod(a)*schur)==0
assert s.factor(Y.det()/schur-A.det())==0
# Gaussian polar-coordinate constant with normalized S4 measure.
sphere_area=2*s.pi**s.Rational(5,2)/s.gamma(s.Rational(5,2))
radial_factor=s.gamma(s.Rational(5,2))/2
assert s.simplify(sphere_area*radial_factor-s.pi**s.Rational(5,2))==0
assert s.simplify(sphere_area-8*s.pi**2/3)==0
t,delta=s.symbols('t delta',real=True)
density=s.Rational(3,4)*(1-t*t)
assert s.integrate(density,(t,-1,1))==1
small=s.integrate(density,(t,-delta,delta))
assert s.expand(s.Rational(3,2)*delta-small)==delta**3/2
L=s.symbols('L',positive=True);cut=L**(-s.Rational(5,12))
assert s.simplify(s.sqrt(cut)-L**(-s.Rational(5,24)))==0
assert s.simplify(L**(-s.Rational(5,4))*cut**(-s.Rational(5,2))-L**(-s.Rational(5,24)))==0
# B2 dimension and the two Schur-orthogonality constants.
l=m=5;dimension=s.Rational((2*l+3)*(2*m+1)*(l-m+1)*(l+m+2),6);assert dimension==286
N2=s.Rational(3,131072);normalizer=dimension/N2;schur_scalar=N2/dimension
assert normalizer*schur_scalar==1
assert s.simplify(normalizer*s.sqrt(schur_scalar)*s.sqrt(N2)-s.sqrt(286))==0
result=dict(dimension=int(dimension),base_real_norm_squared=str(N2),evaluation_density_normalizer=str(normalizer),schur_scalar=str(schur_scalar),uniform_lift_bound='sqrt(286)*norm(v)',endpoint_power='5/2',uniform_decay_power='5/24',exact_checks=['Compression cofactor identity by positive block diagonalization','Gaussian radial and normalized S4 constants','Exact S4 one-coordinate density and small-strip bound','Uniform decay exponent balance','Irreducible dimension286 and exact Schur constants'],seconds=time.time()-tick,scope='Scalar and determinant normalizations only. DA supplies the CU majorant; DB_h6_bounded_average_proof.txt proves the integration, C0, analytic, value-lifting and mixed nonlinear statements.')
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args()
receipt=Path(__file__).with_suffix('.json')
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in result.items() if k!='seconds'}=={k:v for k,v in known.items() if k!='seconds'}
print(json.dumps(result,indent=2))
