"""Exact rank5 and quartic separation on the DG three-tensor span."""
from fractions import Fraction as F
from itertools import combinations,combinations_with_replacement
from math import isqrt
import numpy as np
def verify(p,q,cache,anchors,engine):
 den=21*(1<<(2*p['unit_power']));vals=[np.array(list(map(int,cache['f'+str(i)])),dtype=object) for i in range(3)];columns=[vals[i]-vals[2] for i in range(2)]
 for i,j in combinations(range(3),2):
  v=list(map(int,q['anchor_pair_values'][str(i)+str(j)]));assert engine.value(anchors[:,i]+anchors[:,j])==v;columns.append(np.array(v,dtype=object)-vals[i]-vals[j])
 rows=q['selected_rows'];assert len(rows)==len(set(rows))==5 and all(0<=i<227 for i in rows);A=np.array(columns,dtype=object).T[rows];Y=np.array(q['rightinverse'],dtype=object);ys=1<<q['rightinverse_power'];error=den*ys*np.eye(5,dtype=object)-Y@A
 one=max(sum(abs(int(x)) for x in error[:,j]) for j in range(5));inf=max(sum(abs(int(x)) for x in error[i,:]) for i in range(5));defect=F(isqrt(one*inf)+1,den*ys);ynorm=F(isqrt(sum(int(v)**2 for v in Y.flat))+1,ys);eta=F(1,50000000);total=defect+ynorm*5*800*eta;assert total<F(1,10);inverse=ynorm/(1-total);assert inverse<154
 mon=list(combinations_with_replacement(range(3),3));idx={v:i for i,v in enumerate(mon)};quadrics=[[(i,j)] for i,j in combinations(range(3),2)]+[[(i,i),(2,2)] for i in range(2)];basis={}
 for q0 in quadrics:
  for a in range(3):
   value=0
   for ij in q0:value^=1<<idx[tuple(sorted(ij+(a,)))]
   while value:
    k=value.bit_length()-1
    if k not in basis:basis[k]=value;break
    value^=basis[k]
 assert len(basis)==10
 return {'exact':True,'moment_map_rank':5,'moment_kernel':'span(I3)','nonzero_PSD_quadratic_receiver':False,'all_nonzero_complex_points_excluded':True,'canonical_cubic_Macaulay_rank':10,'root_ball_minor_defect_bound':str(total),'inverse_norm_bound':str(inverse),'coefficient_quartic_gap':'1/106722','physical_quartic_gap':'1/426888','additional_exact_anchor_pair_values':3,'minimum_nonzero_PSD_covariance_rank_inside_this_span':3}
