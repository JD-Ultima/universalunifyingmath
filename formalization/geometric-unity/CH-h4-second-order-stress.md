# CH. A second-order stress obstruction for stationary lifts of the H4 partial value

At every value in the CB rational box, a linearly closed H4 field that satisfies grade-one scalarity and constant norm through second order must have squared covariant gradient greater than \(0.48716697\). In particular, the analytic linearly closed lifts whose first jet vanishes cannot satisfy both quadratic equations through second order. The estimate uses the physical positive metric of BX without further scaling.

## The double-divergence identity

Work at a point in a normal orthonormal frame. Write \(C_b\) for the original four-form at external index \(b\), and assume \(\delta^\nabla C=0\) and the two hook constraints. Put
\[
N=\sum_b\|C_b\|^2,\qquad E_{ij}=\langle C_i,C_j\rangle,\qquad
I_{ij}=\sum_b\langle\iota_iC_b,\iota_jC_b\rangle,\qquad S=2E-I+\tfrac17Ng.
\]
By the CB grade-one identity, grade-one scalarity is exactly \(S=0\). Differentiating the alternating identity \(C_i=\sum_b e^b\wedge\iota_i C_b\) and contracting with the covariant derivative gives \(\sum_i\nabla_i C_i=0\) from closure. The wedge and interior maps used in this manipulation are parallel tensors.

Write \(J_{ji}=\nabla_j C_i\) and use induced tensor curvature \(R_{ji}=[\nabla_j,\nabla_i]\) with bracket correction. Define
\[
\begin{split}
X&=\sum_{i,j}\langle J_{ji},J_{ij}\rangle,\quad
Y=\sum_{i,j,b}\langle\iota_iJ_{jb},\iota_jJ_{ib}\rangle,\\
E_R&=\sum_{i,j}\langle C_i,(R_{ji}C)_j\rangle,\quad
I_R=\sum_{i,j,b}\langle\iota_iC_b,\iota_j(R_{ji}C)_b\rangle.
\end{split}
\]
Expand \(\sum_{ij}(\nabla_j\nabla_i E)_{ij}\) by the product rule. The two product terms with \(\sum_i\nabla_i C_i\) vanish. The crossed first derivatives give \(X\). In the remaining second derivative, \(\sum_{ij}\langle C_i,(\nabla_j\nabla_i C)_j\rangle\), commute \(\nabla_j\) past \(\nabla_i\). The resulting covariant derivative of the external divergence is zero, leaving \(E_R\). Apply the same calculation to \(I\) with the parallel interior maps. Use \(\delta^\nabla C=0\) on its external and form indices. The resulting terms are \(Y+I_R\). Thus
\[
\operatorname{div}\operatorname{div}S-\tfrac17\Delta N=H(C)+2X-Y,\qquad H(C)=2E_R-I_R.
\]
Here \(\operatorname{div}\operatorname{div}S=\sum_{ij}(\nabla_j\nabla_iS)_{ij}\). If scalarity and constant norm hold through second order, the left side is zero.

The exchange of the derivative slot and the external slot is an orthogonal involution, so \(X\ge-\|J\|^2\). The exterior anticommutation rule gives
\[
\begin{split}
\sum_b\Big\|\sum_j e^j\wedge J_{jb}\Big\|^2
  &=\sum_{i,j,b}\langle J_{ib},\iota_i(e^j\wedge J_{jb})\rangle\\
  &=\|J\|^2-Y.
\end{split}
\]
It follows that \(2X-Y\ge-3\|J\|^2\), with equality when the two exchanged slots are skew and \(\sum_j e^j\wedge J_{jb}=0\) for every \(b\). Thus \(\|\nabla C\|^2\ge H(C)/3\) at any point satisfying the second-order conditions.

## Exact curvature bounds and scope

The verifier computes the polarization of \(H\) directly from all 91 ambient curvature matrices
\[
R_{ij}(v)=-\tfrac14[[\Phi_i,\Phi_j],v].
\]
It applies the negative transpose in each of the five covariant slots. In the physical Gram-orthogonal coframe, the two ordered curvature indices contribute inverse metric factors \(1/(g_i g_j)\). The remaining blade and external contractions contribute their respective inverse Gram factors. A separate ordered-index formula on a two-entry rational hook checks both signs, yielding \(E_R=-11/512\), \(I_R=-33/512\), and \(H=11/512\).

On the CB slice \(C(x)=c_0+\sum_j x_jc_j\), the verifier recomputes every monomial coefficient of \(H\). The center and the total variation bound for \(\|x\|_\infty\le 10^{-8}\) are
\[
H(c_0)=
\frac{2510839668219223416601135414352509}{1717986918400000000000000000000000},
\qquad
\big|H(C(x))-H(c_0)\big|\le
\frac{222230151717929555821}{2147483648000000000000000000}.
\]
Therefore \(H(C(x))>1.46150093\) on the entire box, which gives the stated gradient bound. The rigorous interval argument uses both pure and mixed monomials and is independent of any numerical field fit.

By CA there exists an analytic linearly closed extension through the CB value. Let \(s\) denote the base-point geodesic symmetry, realized by inversion of the positive symmetric matrix. It is an isometry, so the pullback commutes with the hook operations and \(\delta^\nabla\). The derivative at the fixed point is minus the identity. Thus \((C-s^*C)/2\) has the original base value and vanishing first derivative, since the covariant slot counts are five and six, respectively. The first derivatives of the pointwise quadratic equations vanish for this lift. The above estimate forces failure of their second-order scalarity/constant-norm equations. For a hypothetical solution of a homogeneous quadratic equation \(Q(C)=q_0\) with symmetry-invariant target, set \(C_0=(C-s^*C)/2\) and \(W=(C+s^*C)/2\). Averaging the two quadratic equations gives \(Q(C_0)=q_0-Q(W)\). Here \(W(o)=0\) and \(\nabla W(o)=\nabla C(o)\), so its value and first derivative make no change to the equation of \(C_0\) at the base point, while its second derivative subtracts the quadratic polarization on the original first jet. The result leaves open solutions with a nonzero first derivative through the CB value.

## Curvature-source signature on the complete BX kernel

The same curvature form is \(SO(5)\)-invariant. Its exact inertia on all 10,388 real BX directions is 8,822 positive, 1,566 negative, and zero null directions. Hence its positive sign on the CB box is a property of the chosen partial value.

For the calculation, use the 57 exact independent highest vectors from CA across the 22 irreducible types of BX. Within a type, the covariant compact transform of each highest has weight \(w\) for a fixed torus. All weights \(w\) here are nonzero. Complex conjugation gives weight \(-w\). The compact real part therefore retains the whole multiplicity pairing: the bilinear invariant form pairs \(w\) with \(-w\) and its equal-weight pairings are zero. On these real parts the physical Gram and \(H\) give the multiplicity Gram and source form, respectively, up to the same common positive factor for each type. This also follows from compact Schur orthogonality on the real irreducible factors. The comparison verifies the exact Gram and form matrices in every block and computes each inertia by rational congruence. Each sign count is repeated by the dimension of the corresponding irreducible module. For example, the ratio \(H/N\) is \(3/4\) on type \((2,0)\), \(-81/140\) on type \((4,0)\), and the two multiplicity eigenvalues on type \((1,0)\) are \(-5/8\) and \(17/40\).

## Reproduction

Run the [verifier](../../computation/geometric-unity/CH_h4_second_order_stress.py)
with the repository's pinned dependencies.

```text
python -X utf8 computation/geometric-unity/CH_h4_second_order_stress.py
python -X utf8 computation/geometric-unity/CH_h4_second_order_stress.py --check-ca
```

The default verifies the CB box, the sign calibration, all 22 multiplicity matrices, and the exact signature. The optional check also reruns the CA highest weight and overlap certificate. Both input archives are pinned by their uncompressed SHA-256 hashes.

[Review97](../../reports/reviews/REVIEW97-GU-CH-2026-11-14.md) records
the independent calculations. External peer review remains outstanding.
