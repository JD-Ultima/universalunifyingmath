# DG. Three-value H4 covariance and an H6 rank-six example

Three linearly independent real algebraic tensors in the full
110-dimensional permutation-sign H4 curvature space have all their
homogeneous Clifford outputs cancel in a covariance. Their exact
three-dimensional span excludes every nonzero complex common zero and
has a strictly positive real quartic residual bound. Separately, an
explicit rank-six H6 hook attains the [CX](CX-hook-sparse-support.md)
external-rank bound with complete scalarity and zero top grade.

## DG1. The full permutation-sign output

Use \(q(U)=Q(U)+(6/7)\|U\|^2\gamma\) for the homogeneous H4
residual, with bilinear extension of the squared norm over the
complex numbers. The source certificate supplies 110 independent
integral tensors in the positive BC coframe and a rational basis
change. Exact arithmetic checks both hooks, coefficient curvature
and the sign on the four adjacent permutation generators.
The complete [CA](CA-h4-full-analytic-lifts.md) character inventory gives
sign multiplicity 110.

The product of two sign tensors has invariant quadratic output.
The grade-one output is symmetric and traceless by
[CD](CD-hook-clifford-identities.md). The fixed symmetric trace-free
output has dimension four. The fixed space in \(V\otimes\Lambda^5V\)
has dimension 223, and removing the exterior trace subtracts 28.
The fixed top output has dimension 28. Thus the complete output
has dimension at most \(4+223-28+28=227\), as in
[CW](CW-h4-fixed-positive-moment.md).

Choose 227 actual unit-frame output coordinates. The invertible
228-row Jacobian below has these 227 outputs followed by scalar
normalization. Its first 227 rows are linearly independent and each
is a derivative of the quadratic output. This proves that the
227 selected coordinates detect vanishing of the complete residual
on this entire permutation-sign domain.

## DG2. A real algebraic contraction certificate

Let \(B\colon\mathbb R^{110}\to W_{\mathrm{sign}}\) be the exact
basis. For \(X\in\mathbb R^{110\times3}\), use the 227 selected
coordinates of \(\sum_r q(BX_r)\) and
\(\|X\|_F^2-1\) as an equation \(f(X)=0\).
The dyadic center \(X_0\) has denominator \(2^{48}\).
Vary 228 selected entries and fix the remaining 102.
The center satisfies \(X_0^T X_0>I_3/5\).

Integer square-root enclosures bound the difference between each
exact basis column or exact center anchor and its unit-frame dyadic
proxy by \(\delta=10^{-15}\). Exact 14-bit-limb Gram products give
\(\|B\|<\sqrt2\). The individual basis columns have norm below
\(1.1\), and the three anchors have norm below one. The 228-column
dyadic derivative uses polarization of the ambient five-term
Clifford reduction, with exact normalization row.
The homogeneous polynomial and polarization bounds from
[DC](DC-h4-four-value-moment.md) give a Lipschitz constant 802 for
the derivative of \(f\). A derivative column changes by at most
\(840\delta\) when restoring the exact tensor coordinates.
The center residual changes by at most \(1600\delta\).

The certificate provides a dyadic inverse \(Y\). Let \(J_D,v_D\)
be the dyadic derivative and residual in the 228 affine variables.
Exact norm estimates give
\[
 \|Y\|<1800,\quad
 e=\|I-YJ_D\|+16\cdot840\delta\|Y\|<1/1000,\quad
 \|Y\|(\|v_D\|+1600\delta)<\tfrac12 10^{-8},\quad
 e+802\|Y\|10^{-8}<1/50.
\]
The fixed-point map \(t\mapsto t-Y f(X_0+St)\), where \(S\)
inserts the 228 entries, is a contraction into the closed radius
\(10^{-8}\) ball. The inverse and the derivative stay invertible.
The resulting isolated root is real algebraic because removing the
fixed coordinate scalings gives rational polynomial equations on the
rational affine slice. Its factor has least singular value greater
than \(2/5\). The three independent tensors \(F_r=BX_r\) therefore
give a positive semidefinite covariance of rank three with
\(\sum_r q(F_r)=0\).

## DG3. No nonzero point on the same support

On this span let \(L\colon\operatorname{Sym}^2\mathbb R^3\to
\text{output}\) be the moment map. The covariance gives \(L(I_3)=0\).
Use the traceless basis \(E_{00}-E_{22},E_{11}-E_{22}\) and
\(E_{ij}+E_{ji}\) for \(i<j\). Three exact anchor-pair evaluations
give the off-diagonal columns. At the root the anchors are within
\(\eta=2\cdot10^{-8}\) of their dyadic proxies, so each column error
is at most \(800\eta\). A verified five-row minor has inverse norm
less than 154 throughout the ball. Thus \(L\) has rank five and
kernel precisely the identity line over the reals or complex numbers.

If \(q(\sum_r x_rF_r)=0\), then \(xx^T=cI_3\).
Rank comparison gives \(c=0\) and hence \(x=0\) over the complex
numbers. The five canonical quadratic generators on this span are
the three products \(x_ix_j\) for \(i<j\) and the two differences
of diagonal squares. Their ideal contains every cubic.

For real \(x\), the traceless part
\(M=xx^T-\|x\|^2I_3/3\) has squared Frobenius norm
\((2/3)\|x\|^4\). Its coefficient vector in the stated traceless
basis satisfies \(\|M\|_F^2\le3\|t\|^2\). The minor inverse bound
therefore gives \(\|q(\sum_r x_rF_r)\|^2\ge
\|x\|^4/106722\). Since \(\|B\|<\sqrt2\) and \(\|X\|_F=1\),
every real tensor \(U\) in this span satisfies
\[
 \|q(U)\|^2\ge\|U\|^4/426888 .
\]
Every positive semidefinite quadratic receiver on this span has
trace zero by \(L(I_3)=0\), so it vanishes. All nonzero positive
semidefinite covariances in the moment kernel on this span have
rank three. The full ambient sign space may admit other ranks.

## DG4. Sharp external rank for the H6 partial equations

Split the positive 14-space as \(R\oplus T\) with orthonormal
axes \(e_0,\ldots,e_5\) in \(R\) and \(u_1,\ldots,u_8\) in \(T\).
Write \(u_{ijk}=u_i\wedge u_j\wedge u_k\) and use
\[
 \phi=u_{123}
   +\tfrac12(u_{147}+u_{246}+u_{257}+u_{345}-u_{156}-u_{367})
   +\tfrac{\sqrt3}{2}(u_{458}+u_{678}).
\]
Direct calculation gives \(\|\phi\|^2=4\) and inner contraction
Gram \((3/2)I_8\). Thus \(\eta=\sqrt2\,*_T\phi\) has squared
norm eight and form moment \(5I_8\). Choose six unit four-blades
\(\omega_b\) as the three complementary pairs generated by
\(u_{1234},u_{1256},u_{1357}\). Their total form moment is \(3I_8\).
Set \(d=(1,1,1,-1,-1,-1)\) and let \(C\) have the two identical
diagonal blocks
\[
 \begin{pmatrix}0&1&-1\\-1&0&1\\1&-1&0\end{pmatrix}.
\]
The rank-six hook is zero on external axes in \(T\), with
\[
 B_a=d_a e_a\wedge\eta+
       \sqrt3\sum_b C_{ab}e_a\wedge e_b\wedge\omega_b,
 \qquad 0\le a\le5.
\]
Every term in row \(a\) contains \(e_a\), so its alternating hook
and exterior square vanish. The contraction hook vanishes by
\(\sum_a d_a=0\) and the column sums of \(C\).
The two orthogonal \(R\)-degrees have external Grams \(8I_6,6I_6\)
and energies 48 and 36. The \(R\)-block of the form moment is
\(20I_6\), the \(T\)-block is \(48I_8\), and the mixed blocks
vanish by \(d^T C=0\). Hence \(N=84\), \(E=14I_6\oplus0_8\),
and \(I+2E=(4N/7)I_{14}\). The top formula of
[CD](CD-hook-clifford-identities.md) gives zero grade 13.
The original nested Clifford calculation verifies these equations
in each radical channel and exhibits nonzero grades five and nine.
This example concerns the algebraic hook before imposing curvature.

    python -X utf8 computation/geometric-unity/DG_h4_rank3_moment.py --full
    python -X utf8 computation/geometric-unity/DG_h6_rank6_sharpness.py

The first command rebuilds all 338 integer evaluations.
Its default mode checks the exact geometry and contraction bounds,
six representative values and all three additional quartic pairs.
