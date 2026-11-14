# CB. A real H4 partial zero with scalarity, grade nine, and the invariant grade-five component

On the positive \(SL(5,\mathbb R)/SO(5)\) model of [BX](BX-h4-curvature-and-closed-fields.md), there is a real unit four-form hook \(C\) in the complete pointwise dual-curvature kernel with
\[
q_1(C)=-\frac67g,\qquad q_9(C)=0,\qquad L_5(C)=0.
\]
Here \(q\) uses the nested Clifford quadratic of BC on the original four-form, with no Clifford-volume factor. \(L_5\) is the orthogonal projection onto the single \(SO(5)\)-invariant line, evaluated against a fixed nonzero vector on that line. The grade-five coefficient at external slot \(13\) and form blade \(6172\) is strictly negative. CA gives a global real analytic linearly closed lift of this value. The displayed quadratic conditions concern the value at the base point and the separating coefficient rules out a full solution through it on the H4 branch.

## Clifford reductions

Use a positive orthonormal frame with Hermitian Clifford generators \(\gamma_a\), and regard \(C_a\) as its real degree-four Clifford symbol. The two hook traces imply both \(\sum_b\gamma_b C_b=0\) and \(\sum_b C_b\gamma_b=0\). With \(S=\sum_b C_b^2\), the nested quadratic reduces to
\[
q_a=-\frac23\left(\{\gamma_a,S\}
 -\sum_b(\gamma_b C_a C_b+C_b\gamma_a C_b+C_b C_a\gamma_b)\right).
\]
Indeed expand the nested brackets in BC and extend the sum to ordered pairs. The three ordered terms have the form \(A_{abc}C_bC_c+C_bC_cA_{abc}+C_bA_{abc}C_c\), multiplied by \(-2/3\), where \(A_{abc}=\gamma_{[a}\gamma_b\gamma_{c]}\). Insert \(\gamma_a\gamma_b\gamma_c-\delta_{ab}\gamma_c+\delta_{ac}\gamma_b-\delta_{bc}\gamma_a\) for \(A_{abc}\). In the first term put the leading triple in order \(a,c,b\), keeping the anticommutator correction. The leading triple then vanishes by the hook sum and the correction contributes \(2\gamma_a S\). For the second and third terms put the leading triple in orders \(c,b,a\) and \(b,a,c\), respectively, and apply the hook sums. The delta terms and this correction give the stated reduction.

Write \(N=\sum_c\|C_c\|^2\), \(E_{aj}=\langle C_a,C_j\rangle\), and \(I_{aj}=\sum_c\langle\iota_a C_c,\iota_j C_c\rangle\). Use normalized spinor trace against \(\gamma_j\) to extract degree one. The anticommutator contributes \(2\delta_{aj}N\). In the two mixed external terms, cyclically move the contracted \(\gamma_b\) next to \(\gamma_j\) and apply \(\{\gamma_b,\gamma_j\}=2\delta_{bj}\) and the hook sums. Each term gives \(2E_{aj}\). In the middle term,
\(\gamma_j C_b=e^j\wedge C_b+\iota_j C_b\).
The reversal signs of degrees five and three are \(+1\) and \(-1\), respectively. Its trace contribution is \(\delta_{aj}N-2I_{aj}\). Hence
\[
q_1=-\frac23Ng+\frac83E-\frac43I,\qquad \operatorname{tr}_gq_1=-12N.
\]
In particular scalarity is the tensor equation \(2E-I+Ng/7=0\). These trace and wedge calculations hold for arbitrary real four-form hooks.

The degree-nine terms are the top exterior parts. The two mixed external terms vanish by the alternating trace. The anticommutator gives \(2e^a\wedge W\) and the middle term gives \(e^a\wedge W\), where \(W=\sum_b C_b\wedge C_b\) in the orthonormal frame. Thus \(q_{9,a}=-\frac23 e^a\wedge W\). All wedge maps \(e^a\wedge \cdot\) jointly detect \(W\) in degree eight. In the covariant BC frame used by the data, \(W=\sum_b C_b\wedge C_b/g_b\), and \(q_{9,a}=-\frac23g_a e^a\wedge W\). The verifier checks the degree-one formula on a two-entry hook and the grade-nine formula on a hook with a nonzero exterior square, using independent nested Clifford products.

## Detectors and invariant equations

The affine data give \(C(x)=c_0+\sum_{j=1}^{34}x_jc_j\), with all 35 tensors rational, in the physical Gram-orthogonal BC coframe. All 35 satisfy the two hook traces and \(R^{T^*}\wedge\star C=0\). Under the coordinate permutations of the defining five-dimensional representation they transform by permutation sign. Odd permutations act on \(V=\operatorname{Sym}^2_0(\mathbb R^5)\) by the same conjugation as their negations in \(SO(5)\). The quadratic tensors are therefore invariant under this \(S_5\) action.

At permutation \(p\), the tangent character is \((f(p)^2+f(p^2))/2-1\), with \(f\) the number of fixed coordinates. Character averaging over all 120 elements yields dimensions 5 in \(\operatorname{Sym}^2 V\) and 28 in \(\Lambda^8 V\). The latter computation uses Newton's identities on the traces of powers of \(p\). Remove the metric to obtain 4 trace-free symmetric directions. The four stored trace-free matrices are invariant under the four adjacent transpositions. Reading \((0,0),(0,1),(0,7),(0,10)\) gives determinant 1. For degree eight, the verifier forms 28 exact, unnormalized \(S_5\) averages of documented rational coordinate blades. Reading the 28 chosen coordinates on these forms gives determinant \(-633318697598976\). These two invertible readings detect all of scalarity and grade nine on the affine slice.

To count the compact invariant outputs, take the \(B_2\) character of \(V=\operatorname{Sym}^2_0(\mathbb R^5)\) and multiply the characters of \(V\otimes\Lambda^r V\) by \(\prod_{\alpha>0}(1-z^{-\alpha})\) for \(\alpha=(1,-1),(1,1),(1,0),(0,1)\). The constant coefficient counts the trivial summands. This gives one each in degrees 1, 5, and 9, and zero in \(\Lambda^8 V\). Therefore the invariant degree-nine pairing already vanishes identically for four-form hooks by the wedge formula. It is also checked directly against the BN receiver.

The degree-five invariant covariant tensor is the primitive integer tensor from the antisymmetrized trace
\[
\rho_{a;i_1\ldots i_5}
=\operatorname{tr}\left(Q_a\sum_{\pi\in S_5}\operatorname{sgn}(\pi)
 Q_{i_{\pi(1)}}\cdots Q_{i_{\pi(5)}}\right),\qquad Q_a=5\Phi_a.
\]
The verifier reconstructs its 1143 nonzero entries, divides their common gcd \(50000\), and checks all ten infinitesimal invariant equations. The physical moment is \(L_5=\sum_{a,m}q_{5,a;m}\rho_{a;m}/(g_a\prod_{j\in m}g_j)\) using the primitive \(\rho\). It provides the last row. The remaining scalar condition is \(N-1=0\).

## Rational certificate and separating coefficient

The witness supplies 34 symmetric rational matrices representing the equation vector at \((1,x)\). The row order is four scalarity rows, 28 exterior-square rows, the norm row, and the invariant grade-five row. Off-diagonal matrix entries are half the corresponding monomial coefficients. There are no numerical optimizer inputs.

With the given invertible rational preconditioner \(M\), form \(T(x)=x-MF(x)\). On \(\|x\|_\infty\le 10^{-8}\), the exact contraction bound is below \(8.627\cdot10^{-6}\) and the image has sup norm below \(8.664\cdot10^{-11}\). Thus \(T\) has a fixed point within the box and \(F\) vanishes there.

For the full grade-five receiver at \((a,m)=(13,6172)\), all 630 affine monomial coefficients are reconstructed from the unprojected Clifford expansion. A separate nested multiplication computes its center value exactly. The coefficient throughout the closed box is below \(-10.36696\) in the covariant frame, and its total absolute variation from the center is below \(1.763\cdot 10^{-6}\). In particular complete grade-five vanishing fails for the exact value. A constant rescaling to any desired nonzero physical branch magnitude preserves this separation.

## Reproduction

Run the [verifier](../../computation/geometric-unity/CB_h4_partial_value.py)
from the repository with the pinned dependencies. The verifier uses
exact integer and rational arithmetic.

    python -X utf8 computation/geometric-unity/CB_h4_partial_value.py
    python -X utf8 computation/geometric-unity/CB_h4_partial_value.py --full

Both modes independently verify all rational affine geometries, the detector dimensions and nonzero minors, every scalarity, wedge, and norm matrix entry, the exact box bounds, and all 630 separating grade-five coefficients. The default reconstructs 10 selected \(L_5\) matrix entries including the center. The full mode independently computes all 630 before applying the same contraction check. Run the full mode to check the complete existence certificate.

The proof supplies one exact admissible value and a linear closed
analytic extension by [CA](CA-h4-full-analytic-lifts.md). Geodesic
symmetry $s_o$ is a global isometry with differential $-I$ at the
origin. For a chosen lift $C$, the lift $(C-s_o^*C)/2$ has the same
value and zero covariant first derivative, since $C$ has five
covariant slots and its derivative has six. The operator $\delta_\nabla$
commutes with isometries. Thus the first derivatives of the partial
quadratic equations can also be made to vanish at the point. Their
vanishing on an open set remains open.

[Review94](../../reports/reviews/REVIEW94-GU-CB-2026-11-14.md)
records the separate calculations. External peer review remains
outstanding.
