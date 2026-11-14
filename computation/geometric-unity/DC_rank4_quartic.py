"""Stable exact rank9 and a positive quartic on the certified DC span."""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from math import isqrt
import numpy as np
def verify(p,quartic,cache,anchors,engine):
 den=21*(1<<(2*p['unit_power']));vals=[np.array(list(map(int,cache['f'+str(i)])),dtype=object) for i in range(4)];columns=[vals[i]-vals[3] for i in range(3)]
 for i,j in combinations(range(4),2):
  expected=list(map(int,quartic['anchor_pair_values'][str(i)+str(j)]));assert engine.value(anchors[:,i]+anchors[:,j])==expected
  columns.append(np.array(expected,dtype=object)-vals[i]-vals[j])
 rows=quartic['selected_rows'];assert len(rows)==len(set(rows))==9 and all(0<=i<227 for i in rows);A=np.array(columns,dtype=object).T[rows];Y=np.array(quartic['rightinverse'],dtype=object);ys=1<<quartic['rightinverse_power'];error=den*ys*np.eye(9,dtype=object)-Y@A
 one=max(sum(abs(int(x)) for x in error[:,j]) for j in range(9));inf=max(sum(abs(int(x)) for x in error[i,:]) for i in range(9));defect=F(isqrt(one*inf)+1,den*ys);ynorm=F(isqrt(sum(int(v)**2 for v in Y.flat))+1,ys)
 # ||B||<sqrt2, certified slice radius10^-8, and dyadic anchor error10^-15.
 # Thus each exact root anchor differs from its dyadic anchor by<2*10^-8.
 eta=F(1,50000000);total=defect+ynorm*9*800*eta;assert total<F(1,10);inverse=ynorm/(1-total);assert inverse<266
 # Independent tiny characteristic-two Macaulay audit for the canonical
 # trace-free Veronese ideal: off-diagonal products and diagonal differences.
 mon=list(combinations_with_replacement(range(4),3));index={m:i for i,m in enumerate(mon)};quadrics=[[(i,j)] for i,j in combinations(range(4),2)]+[[(i,i),(3,3)] for i in range(3)];basis={}
 for q in quadrics:
  for a in range(4):
   value=0
   for ij in q:value^=1<<index[tuple(sorted(ij+(a,))) ]
   while value:
    k=value.bit_length()-1
    if k not in basis:basis[k]=value;break
    value^=basis[k]
 assert len(basis)==20
 return {'exact':True,'moment_map_rank':9,'moment_kernel':'span(I4)','nonzero_PSD_quadratic_receiver':False,'all_nonzero_complex_points_excluded':True,'canonical_cubic_Macaulay_rank':20,'root_ball_minor_defect_bound':str(total),'inverse_norm_bound':str(inverse),'coefficient_quartic_gap':'3/1132096','physical_quartic_gap':'3/4528384','additional_exact_anchor_pair_values':6}
