# DN. Real H4 exclusion on 960 high-Cartan directions

In the positive H4 curvature kernel, fix the rotation of the first
two defining axes. Let \(L\) be the sum of its complex weight spaces
of weights at least four, and \(W=\{z+\bar z:z\in L\}\).
The real dimension of \(W\) is 960. There is a rational real receiver
\(Y\) such that the complete homogeneous H4 residual
\(q_*(B)=Q(B)+(6/7)\|B\|^2\gamma\) satisfies
\[
 \langle Y,q_*(B)\rangle>\|B\|^2/20,qquad
 \|q_*(B)\|^2>\|B\|^4/3
 \quad\text{for all }0\ne B\in W.
\]
Compact equivariance gives the same exclusion on every compact
translate of \(W\). [DJ](DJ-hook-weight-and-commutator-bounds.md)
shows that the complex subspace \(L\), of dimension 480, consists
of zero-norm homogeneous solutions. The real claim concerns
\(L\oplus\bar L\) with compact conjugation.

## DN1. Complete source and receiver spaces

The compact centralizer of the defining rotation is
\(SO(2)\times SO(3)\). On the positive-weight part the exact
highest multiplicities are as follows.

| First weight | Second spins in order | Highest multiplicities |
| ---: | --- | --- |
| 4 | 0, 1, 2, 3, 4, 5 | 4, 21, 17, 17, 6, 2 |
| 5 | 0, 1, 2, 3, 4 | 5, 5, 9, 4, 2 |
| 6 | 1, 2, 3 | 3, 1, 1 |
| 7 | 0 | 1 |

The 98 highest vectors have compact rational real and imaginary
parts. The verifier builds the infinitesimal action directly from
the physical symmetric five-by-five matrices. It checks both hooks,
dual curvature, both Cartan actions and the \(SO(3)\) raising
equation on each vector. Exact
[BX](BX-h4-curvature-and-closed-fields.md) character branching
supplies every highest multiplicity and sums to 480.

The receiver has grades one, five and nine in the BC covariant
coframe. Its grade-one trace and all four infinitesimal centralizer
actions vanish exactly. It is also fixed by the diagonal matrix
\(s=\operatorname{diag}(1,-1,1,-1,1)\) in \(SO(5)\), which acts
as conjugation on the chosen split source tensors. Their real and
imaginary parts have opposite \(s\) parity. Thus the receiver has
no mixed real-imaginary entry on the highest lines, and rotation
invariance equates the two diagonal blocks.

The spin representations of \(SO(3)\) have real form. Schur's
lemma and orthogonality reduce the invariant quadratic inequality
to the 15 highest multiplicity matrices \(A_{w,l}\) and their
physical Grams \(G_{w,l}\). Every \(A_{w,l}-G_{w,l}/20\) has
strictly positive leading principal minors by exact integer
elimination. This also proves independence of the source basis.
The trace-free receiver has exact squared norm
\[
 \|Y\|^2=
 \frac{31603334182253213234930488845}
      {4436777100798802905238461218816}<\frac3{400}.
\]
Cauchy–Schwarz gives the full residual bound stated above.

## DN2. Exact Clifford reconstruction

Use the five-term reduction of [CD](CD-hook-clifford-identities.md)
on actual Clifford grade-four hooks. The receiver trace is zero,
so the homogeneous norm correction does not contribute.
The modular routine evaluates the adjoint of the polarized
quadratic and forms all 15 complete symmetric Gram matrices.
Separate sparse or full-spin formulas check independent entries.

Let \(d_Y,d_B\) be common raw denominators of the receiver and
source, and \(g_i\) the positive BC metric entries. A common
denominator of all matrix entries is
\[
 D=3d_Yd_B^2\prod_i\operatorname{numerator}(g_i)^3 .
\]
There are three tensor slots in a summand, two source slots and one
receiver slot. The three tensor slots can contain the same axis
at most six times. The mask relation and the paired external indices
give an even exponent for every raw metric scale.

In the terms with \(\gamma_a\) the two source external indices
are equal and the inserted generator equals the receiver external
index \(a\). The exterior blade parity gives one more occurrence
exactly on \(a\). In the terms with \(\gamma_b\) the source
external indices are \(a,b\), the receiver external index is
\(a\), and the exterior blade parity gives one more occurrence
exactly on \(b\). This proves cancellation of all half-integral
scalings. The remaining metric denominators divide the factor
in \(D\).

With normalized Clifford norm, \(\|XY\|\le\sqrt{128}\|X\|\|Y\|\).
Estimating each of the five terms gives
\(\|Q(B)\|<150\|B\|^2\). Polarization therefore bounds each
matrix entry by \(150\|Y\|\|B_i\|\|B_j\|\).
The ten distinct reconstruction primes have a 200-bit product,
more than twice the certified 190-bit absolute numerator bound.
Each prime is one modulo 120, so the required BC scale factors
have square roots in its field. Changing such roots changes only
orthonormal frame signs. The modular products use binary64 for
exact integers only, with every possible partial sum bounded
by \(2\cdot896(p-1)^2<2^{53}\) before reduction. Chinese remainder
reconstruction identifies every exact signed numerator.

    python -X utf8 computation/geometric-unity/DN_high_cartan.py
    python -X utf8 computation/geometric-unity/DN_high_cartan.py --full

Both modes check all source geometry and positivity. Full mode
regenerates all ten modular Gram images. This exclusion applies to
each stated high-Cartan support and leaves general mixed weights open.
