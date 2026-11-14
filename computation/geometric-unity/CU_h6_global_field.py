"""Exact block identities for CU's global H6 lift."""
import sympy as s
# Compact real metric calibration: the Riccati M is an actual positive Gram.
E=s.Matrix([[1,1,0,0,0],[s.I,-s.I,0,0,0],[0,0,1,1,0],[0,0,s.I,-s.I,0],[0,0,0,0,1]])
order=[1,3,4,0,2];Ei=E.inv()[order,:]
G=(Ei*Ei.T).inv();assert G.det()==16
yy=s.symbols('y0:15',real=True);Y=s.zeros(5);j=0
for i in range(5):
    for k in range(i,5):Y[i,k]=Y[k,i]=yy[j];j+=1
yreal=Ei*Y*Ei.T;Ar=yreal[:2,:2];Cr=yreal[:2,3:];Fr=yreal[3:,3:]
assert (Fr-Ar.conjugate()).applyfunc(s.expand)==s.zeros(2) and (Cr.T-Cr.conjugate()).applyfunc(s.expand)==s.zeros(2)
rows=Ei[:2,:].col_join(-Ei[3:,:]);M=Cr.row_join(-Ar).col_join((-Fr).row_join(Cr.T))
assert (rows*Y*rows.conjugate().T-M).applyfunc(s.expand)==s.zeros(4)
assert rows.rank()==4
def sym2(prefix):
    a,b,c=s.symbols(prefix+'0 '+prefix+'1 '+prefix+'2');return s.Matrix([[a,b],[b,c]])
A=sym2('A');S=sym2('S');C=s.Matrix(2,2,s.symbols('C0:4'));b=s.Matrix(s.symbols('b0:2'));v=s.Matrix(s.symbols('v0:2'));d=s.symbols('d')
R=C-A*S;F=S*C+C.T*S-S*A*S;e=S*b+R.T*v
Omega=s.zeros(4);Omega[:2,2:]=s.eye(2);Omega[2:,:2]=-s.eye(2)
H=C.row_join(-A).col_join(F.row_join(-C.T))
assert (H.T*Omega+Omega*H).applyfunc(s.expand)==s.zeros(4)
Z=sym2('Z');Ric=lambda U:F-U*C-C.T*U+U*A*U
epsilon=s.symbols('epsilon')
variation=Ric(S+epsilon*Z).diff(epsilon).subs(epsilon,0)
assert (variation+Z*R+R.T*Z).applyfunc(s.expand)==s.zeros(2)
y=A.row_join(b).row_join(C).col_join(b.T.row_join(s.Matrix([[d]])).row_join(e.T)).col_join(C.T.row_join(e).row_join(F))
n=s.eye(5);n[2,:2]=v.T;n[3:,:2]=S;ni=s.eye(5);ni[2,:2]=-v.T;ni[3:,:2]=-S
assert ni*n==s.eye(5)
yp=(ni*y*ni.T).applyfunc(s.expand)
assert yp[3:,3:]==s.zeros(2) and yp[2,3:]==s.zeros(1,2)
assert (yp[:2,3:]-R).applyfunc(s.expand)==s.zeros(2)
# Verify the parabolic factorization over Q(R,Ap,bp,q)/(q^2-4detR).
R=s.Matrix(2,2,s.symbols('r0:4'));Ap=sym2('x');bp=s.Matrix(s.symbols('u0:2'));q=s.symbols('q',nonzero=True)
D=s.diag(q,1);L=2*R*D.inv().T;a=q**-2;ell=bp/a;N=(Ap-ell*ell.T)*L.inv().T
h=s.zeros(5);h[:2,:2]=L;h[:2,2]=ell;h[:2,3:]=N;h[2,2]=a;h[3:,3:]=D
Gi=s.zeros(5);Gi[:2,3:]=s.eye(2)/2;Gi[3:,:2]=s.eye(2)/2;Gi[2,2]=1
target=Ap.row_join(bp).row_join(R).col_join(bp.T.row_join(s.Matrix([[q**-4]])).row_join(s.zeros(1,2))).col_join(R.T.row_join(s.zeros(2,1)).row_join(s.zeros(2)))
assert (h*Gi*h.T-target).applyfunc(s.cancel)==s.zeros(5)
assert s.factor(target.det())==s.factor(R.det()**2/q**4)
relation=q*q-4*R.det()
def zero_mod(expr):
    num=s.fraction(s.cancel(expr))[0]
    return s.rem(num,relation,q)==0
assert zero_mod(h.det()-1) and zero_mod(L.det()-q) and D.det()==q
assert zero_mod(target.det()-s.Rational(1,16))
base=[]
for i in range(5):
    for j in range(i,5):
        if (i,j) in [(0,1),(2,3),(4,4)]:continue
        z=s.zeros(5);z[i,j]=z[j,i]=1;base.append(z)
for i,j in [(0,1),(2,3)]:
    z=s.zeros(5);z[i,j]=z[j,i]=1;z[4,4]=-4;base.append(z)
derivatives=[]
for z in base:
    z=z[order,order];dS=z[3:,3:];dv=2*z[2,3:]
    derivatives.append([dS[0,0],dS[0,1],dv[0]/2,dS[1,1],dv[1]/2])
assert s.Matrix(derivatives).T==s.eye(14)[[0,1,3,8,9],:]
print("The symbolic positive Gram, Riccati derivative, transverse frame and block factorization checks passed.")
