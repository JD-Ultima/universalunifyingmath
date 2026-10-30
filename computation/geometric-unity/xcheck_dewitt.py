"""Independent cross-check (coordinator) of the DeWitt fiber / chimeric signatures. Run: python3 xcheck_dewitt.py"""
import itertools, numpy as np
def sym_basis(n):
    B=[]
    for i in range(n):
        for j in range(i,n):
            m=np.zeros((n,n)); m[i,j]=m[j,i]=1; B.append(m)
    return B
def sig(M):
    w=np.linalg.eigvalsh(M); return (int((w>1e-9).sum()), int((w<-1e-9).sum()), int((abs(w)<=1e-9).sum()))
for name,g in [("(1,3)",np.diag([1,-1,-1,-1.])),("(3,1)",np.diag([-1,1,1,1.])),("(4,0)",np.eye(4))]:
    gi=np.linalg.inv(g); B=sym_basis(4)
    for lam in [0,0.2,0.25,0.3,0.5,1]:
        G=np.array([[np.trace(gi@h@gi@k)-lam*np.trace(gi@h)*np.trace(gi@k) for k in B] for h in B])
        fs=sig(G); bs=sig(g)
        print(name,"lam",lam,"fiber",fs,"fiber+base",(fs[0]+bs[0],fs[1]+bs[1]),"(-fiber)+base",(fs[1]+bs[0],fs[0]+bs[1]))
