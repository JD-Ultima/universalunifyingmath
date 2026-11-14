"""Exact algebraic checks for finite CU-plus-bounded-field degeneration."""
from pathlib import Path
from fractions import Fraction as F
from itertools import combinations
import argparse,json,sys,time
parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--source-dir',type=Path);parser.add_argument('--write-receipt',action='store_true');args=parser.parse_args()
HERE=Path(__file__).resolve().parent
SOURCE=args.source_dir or (HERE if (HERE/'BO_h6_full_kernel.py').exists() else HERE.parents[2]/'universalunifyingmath/computation/geometric-unity')
sys.path.insert(0,str(SOURCE));import BO_h6_full_kernel as BO
s=BO.s;tick=time.time();weights,inverse,curvature,roots=BO.build_model()
metric=s.zeros(5);metric[0,1]=metric[1,0]=metric[2,3]=metric[3,2]=2;metric[4,4]=1
p=[]
for i in range(5):
    for j in range(i,5):
        if (i,j) in ((0,1),(2,3),(4,4)):continue
        e=BO.matrix_units(i,j);p.append((e if i==j else e+e.T)*metric)
for i,j in [(0,1),(2,3)]:p.append((BO.matrix_units(i,j)+BO.matrix_units(j,i)-4*BO.matrix_units(4,4))*metric)
dual=[sum((s.Rational(inverse[a][j])*p[j] for j in range(14)),s.zeros(5)) for a in range(14)]
k=[(BO.matrix_units(i,j)-BO.matrix_units(j,i))*metric for i,j in combinations(range(5),2)]
W=[1,3];outside=[0,2,4];conditions=s.Matrix([[a[i,j] for a in k] for i in outside for j in W])
stab=[sum((v[i]*k[i] for i in range(10)),s.zeros(5)) for v in conditions.nullspace()]
assert len(stab)==7
b0={(0,811):F(1),(1,4875):F(-1),(1,8971):F(1),(8,783):F(1)}
for K in stab:
    action=[[(b,F(s.trace(p[b]*(K*Q-Q*K)))) for b in range(14) if s.trace(p[b]*(K*Q-Q*K))] for Q in dual]
    result=BO.root_action(b0,action);expected={key:F(5*sum(K[i,i] for i in W))*v for key,v in b0.items()};expected={key:v for key,v in expected.items() if v}
    assert result==expected
# Determinants controlling the complete five transverse derivatives.
R=s.Matrix(2,2,s.symbols('r0:4'));D=s.Matrix(2,2,s.symbols('d0:4'));z0,z1,z2=s.symbols('z0 z1 z2');Z=s.Matrix([[z0,z1],[z1,z2]])
symcoords=lambda M:s.Matrix([M[0,0],M[0,1],M[1,1]])
syl=symcoords(R.T*Z+Z*R).jacobian([z0,z1,z2]);cong=symcoords(D*Z*D.T).jacobian([z0,z1,z2])
assert s.factor(syl.det()-4*R.det()*s.trace(R))==0
assert s.factor(cong.det()-D.det()**3)==0
q,tr=s.symbols('q tr',positive=True);a=q**-2
J0=s.cancel((q**3/(q*q*tr))*(a*a*q/(q*q)))
assert J0==q**-4/tr
# Hamiltonian positive-graph identities, without diagonalizing R.
aa,ab,ac,ss,st,su=s.symbols('aa ab ac ss st su')
A=s.Matrix([[aa,ab],[ab,ac]]);S=s.Matrix([[ss,st],[st,su]]);C=s.Matrix(2,2,s.symbols('c0:4'))
R=C-A*S;Fmat=S*C+C.T*S-S*A*S
H=C.row_join(-A).col_join(Fmat.row_join(-C.T));J=s.diag(1,1,-1,-1);M=J*H;U=s.eye(2).col_join(S)
assert (M*U-J*U*R).applyfunc(s.expand)==s.zeros(4,2)
yc=A.row_join(C).col_join(C.T.row_join(Fmat))
assert ((-S).row_join(s.eye(2))*yc-R.T*s.eye(2).row_join(S)).applyfunc(s.expand)==s.zeros(2,4)
T=s.eye(4);T[2:,:2]=S;Ti=s.eye(4);Ti[2:,:2]=-S
assert (Ti*H*T-R.row_join(-A).col_join(s.zeros(2).row_join(-R.T))).applyfunc(s.expand)==s.zeros(4)
# Compact metric Gram calibration, needed to obtain det Delta.
E=s.Matrix([[1,1,0,0],[s.I,-s.I,0,0],[0,0,1,1],[0,0,s.I,-s.I]])[:,[1,3,0,2]]
assert E.conjugate().T*E==2*s.eye(4)
yy=s.symbols('y0:10',real=True);Y=s.zeros(4);j=0
for i in range(4):
    for k0 in range(i,4):Y[i,k0]=Y[k0,i]=yy[j];j+=1
yc=E.inv()*Y*E.inv().T;Mc=yc[:2,2:].row_join(-yc[:2,:2]).col_join((-yc[2:,2:]).row_join(yc[2:,:2]))
assert (Mc-J*E.conjugate().T*Y*E*J/4).applyfunc(s.expand)==s.zeros(4)
assert abs(E.det())**2==16
# Same-center scaling in a common real Cholesky frame.
rho,delta=s.symbols('rho delta',positive=True);L=s.Matrix(4,4,s.symbols('l0:16'));w=s.Matrix(s.symbols('w0:4'))
F0=L.row_join(s.zeros(4,1)).col_join(w.T.row_join(s.Matrix([[delta]])))
Fr=F0*s.diag(rho,rho,rho,rho,rho**-4);e5=s.Matrix([0,0,0,0,1])
assert (Fr*Fr.T-rho*rho*F0*F0.T-(rho**-8-rho*rho)*delta*delta*e5*e5.T).applyfunc(s.expand)==s.zeros(5)
assert s.cancel((rho*rho*q)**-4/(rho*rho*tr)/(q**-4/tr))==rho**-10
# Sherman-Morrison coefficient remains bounded for all other centers.
x,alpha,beta,gamma=s.symbols('x alpha beta gamma',positive=True)
assert s.simplify(x/(1+x*alpha)-1/alpha+1/(alpha*(1+x*alpha)))==0
result=dict(stabilizer_generators=7,highest_determinant_weight=5,exact_checks=['All W-stabilizer generators act on b0 by5 trace_W','Symmetric Sylvester determinant4detR trR','Symmetric congruence determinant(detD)^3','Full transverse volume J0=q^-4/trR','Hamiltonian positive-graph block triangularization','Exact center-coordinate identity v^T=b4^T A4^-1 E4[I;S]','Positive metric Gram calibration','Common real Cholesky path and exact rho^-10 scaling','Sherman-Morrison bounded coefficient identity'],seconds=time.time()-tick,scope='Finite exact algebra used by DA_h6_mixed_growth_proof.txt. The proof separately supplies positivity, norm inequalities, compact multiplier bounds and finite-group elimination. No claim about arbitrary bounded fields or infinite CU measures.')
receipt=HERE/'DA_h6_mixed_growth.json'
if args.write_receipt:receipt.write_text(json.dumps(result,indent=2))
else:
    known=json.loads(receipt.read_text());assert {k:v for k,v in known.items() if k!='seconds'}=={k:v for k,v in result.items() if k!='seconds'}
print(json.dumps(result,indent=2))
