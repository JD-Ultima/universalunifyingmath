# DJ. Cartan support bounds and Clifford commutator energy

Consider the compact real pure H4 and H6 hooks on the positive
14-dimensional tangent space. Write \(N=\|B\|^2\). The grade-one
equation restricts the fraction of \(N\) carried by high absolute
Cartan weights. Independently, the two hook identities give a strict
lower bound on the quartic commutator energy. The full equations
strengthen that lower bound.

## DJ1. High Cartan weights

Restrict the defining real five-dimensional representation to rotation
in one coordinate plane. In a unitary weight basis, the traceless
symmetric-square tangent representation has weights
\[
 (2,1,1,1,0,0,0,0,0,0,-1,-1,-1,-2).
\]
The complex receivers in degrees \(1,5,9,13\) are subrepresentations
of \(V\otimes\Lambda^rV\), with absolute weight bounds \(4,7,7,4\),
respectively. The homogeneous residual, including its bilinear norm
term, is equivariant. On the sum of weights at least four, all its
polarized products have weight at least eight and therefore vanish.
The bilinear norm itself vanishes by weight pairing.

The exact character inventories [BX](BX-h4-curvature-and-closed-fields.md)
and [BO](BO-h6-full-curvature-kernel.md) give dimensions 480 for
this complex linear space in H4 and 573 in H6 after removing the
type \((2,0)\) by [BU](BU-h6-second-jet-obstruction.md).
Compact conjugation takes these subspaces to the negative weights.
There is no nonzero real vector in either positive weight space.
The homogeneous quadratic ideal therefore has a nontrivial complex
projective zero set on each full kernel. To study nonzero-norm complex
solutions by ideal saturation one may saturate by the bilinear norm.

For the real scalarity bounds, let \(P_j\) be the real orthogonal
projection onto the tangent slots of absolute weight \(j\), with
\((\dim P_0,\dim P_1,\dim P_2)=(6,6,2)\).
For a unitary tensor monomial, let \(n_j\) count its form slots of
absolute weight \(j\), and let \(s\) be its absolute external weight.
Set \(u_j=n_j-2\mathbf{1}_{s=j}\) for H4 and
\(u_j=n_j+2\mathbf{1}_{s=j}\) for H6. Taking the trace of scalarity
against \(P_j\) gives
\[
 \sum_{\alpha} u_j(\alpha)|B_{\alpha}|^2
 =\begin{cases}
    (N/7)\dim P_j,&\text{H4},\\
    (4N/7)\dim P_j,&\text{H6}.
  \end{cases}
\]
Here the exterior-slot moment and external Gram enter with signs
\(-2\) and \(+2\), respectively, by [CB](CB-h4-certified-partial-value.md) and
[BS](BS-h6-certified-partial-value.md).

Let \(w\) be the total signed weight of a monomial. Pointwise
incidence inequalities on all tensor slots give the following bounds.
The two H4 inequalities hold even where some \(u_j\) are negative.

| Branch | Threshold \(t\) | Majorant of \(\mathbf{1}_{|w|\ge t}\) | Upper fraction of \(N\) |
| --- | ---: | --- | ---: |
| H4 | 6 | \((2u_0+4u_1+3u_2)/7\) | \(6/7\) |
| H4 | 7 | \((u_0+2u_1+u_2)/5\) | \(4/7\) |
| H6 | 6 | \((2u_1+3u_2)/13\) | \(72/91\) |
| H6 | 7 | \(u_2/3\) | \(8/21\) |

The verifier enumerates all \(14\binom{14}{4}\) and
\(14\binom{14}{6}\) incidences and checks each majorant. In particular
real scalarity in either branch requires some mass below absolute
weight six.

## DJ2. A strict Clifford energy bound

Use the actual Hermitian Clifford inputs \(C_a\), of degree
\(k=4\) in H4 and degree \(k=8\) in H6. Let
\(S=\sum_a C_a^2\), and use the normalized Hilbert–Schmidt inner
product \(\langle U,V\rangle=\operatorname{tr}(U^\dagger V)/128\).
The Clifford identification is isometric. The hooks give
\(\sum_a\gamma_a C_a=\sum_a C_a\gamma_a=0\).
Put \(K_{ab}=[C_a,C_b]\) and \(\Gamma_{ab}=\gamma_a\gamma_b\)
for \(a\ne b\). Cyclicity and the Clifford relation give
\[
 \operatorname{tr}\left(\sum_{a,b}\gamma_a C_bC_a\gamma_b\right)/128
 =2N-\operatorname{tr}\left(\sum_{a,b}\gamma_b C_bC_a\gamma_a\right)/128
 =2N.
\]
The corresponding expression with \(C_aC_b\) in place of
\(C_bC_a\) vanishes by the same hook sums. Subtracting and grouping
unordered pairs yields \(\sum_{a<b}\langle K_{ab},\Gamma_{ab}\rangle=-N\).

Write \(D_{ab}=K_{ab}+(N/91)\Gamma_{ab}\). Orthogonality to the
array \(\Gamma\) gives
\[
 \sum_{a<b}\|K_{ab}\|^2=\frac{N^2}{91}+\sum_{a<b}\|D_{ab}\|^2.
\]
The Jacobi identity for the three Clifford components gives
\[
 \frac{N}{91}J(C)_{abc}=
 [C_a,D_{bc}]+[C_b,D_{ca}]+[C_c,D_{ab}],
 \qquad J(C)_{abc}=[C_a,\Gamma_{bc}]+\text{cyclic}.
\]
The map \(J\) is equivariant and the traceless hook is irreducible
for the orthogonal group. Thus its squared norm is a scalar multiple
of \(N\). Calibrate with \(C_0=e^0\wedge T\),
\(C_1=-e^1\wedge T\), all other inputs zero, with \(T\) a unit
\((k-1)\)-blade avoiding axes zero and one. Triples containing both
external axes cancel by Clifford multiplication. Those containing
exactly one contribute squared norm four exactly when one of the
other axes lies in \(T\) and the other lies outside it. Therefore
\[
 \sum_{a<b<c}\|J(C)_{abc}\|^2=\lambda_k N,
 \qquad \lambda_k=4(k-1)(13-k).
\]

For every spinor matrix \(D\), Hermiticity gives
\[
 \sum_a\|[C_a,D]\|^2
 \le 2\sum_a(\|C_aD\|^2+\|DC_a\|^2)
 \le 4\|S\|_{\rm op}\|D\|^2 .
\]
Applying Cauchy–Schwarz to the three cyclic terms and summing bounds
the squared right-hand side of the Jacobi identity by
\(12\|S\|_{\rm op}\sum_{a<b}\|D_{ab}\|^2\). Thus for \(N>0\),
\[
 \sum_{a,b}\|[C_a,C_b]\|^2
 \ge \left(\frac{2}{91}
      +\frac{\lambda_k}{6\cdot91^2}\frac{N}{\|S\|_{\rm op}}\right)N^2 .
\]
Every \(S\) here is positive semidefinite with trace \(64N\) on
each chirality space, since no degree-14 part occurs.
This gives \(\|S\|_{\rm op}\le64N\) without the field equations.
[CD](CD-hook-clifford-identities.md) gives smaller caps
\(2N\) for H4 and \((78/43)N\) for H6 on the full nonlinear equation.

| Branch | Bound in units of \(N^2\) on a hook | Bound on a full solution |
| --- | ---: | ---: |
| H4 | \(5833/264992\) | \(191/8281\) |
| H6 | \(2501/113568\) | \(3257/138411\) |

## DJ3. Stronger low-rank energy bound

Let \(E_{ac}=\langle C_a,C_c\rangle\) be the external Gram.
The grade-two commutator coefficients form an array
\(M_{ab,cd}=\langle K_{ab},\Gamma_{cd}\rangle\) with antisymmetry
in each index pair. Define its Ricci contraction by
\(\operatorname{Ric}(M)_{ac}=\sum_b M_{ab,cb}\).
One has exactly
\[
 \operatorname{Ric}(M)=-2E .
\]
Indeed the trace formula is
\(-\sum_b\operatorname{tr}([C_a,C_b]\gamma_c\gamma_b)/128\).
The summand with the product \(C_bC_a\) vanishes after summation
by cyclicity and the left hook sum. In the other summand use
\(\sum_b C_b\gamma_c\gamma_b=2C_c\), by the right hook sum.
Terms with equal antisymmetric indices vanish by the trace of a commutator.

For a symmetric matrix \(H\) define the array
\[
 L(H)_{ab,cd}
 =\delta_{ac}H_{bd}-\delta_{ad}H_{bc}
  -\delta_{bc}H_{ad}+\delta_{bd}H_{ac}.
\]
On unordered pairs, its inner product with \(M\) is
\(\langle H,\operatorname{Ric}(M)\rangle\).
Also \(\operatorname{Ric}(L(H))=(\operatorname{tr}H)g+12H\).
Thus the scalar and traceless parts of the orthogonal projection
onto the image of \(L\) contribute squared norms
\(N^2/91\) and \(\|E-Ng/14\|_F^2/3\), respectively.
Discarding the remaining components and summing over ordered
commutator pairs gives
\[
 \sum_{a,b}\|[C_a,C_b]\|^2\ge
 \frac{2N^2}{91}+\frac23\|E-Ng/14\|_F^2 .
\]
For external rank \(r>0\), this is at least
\(\left[2/91+(2/3)(1/r-1/14)\right]N^2\).

    python -X utf8 computation/geometric-unity/DJ_hook_bounds.py

The verifier also checks the two hook sums, trace calibration,
Ricci contraction, and Jacobi norm on four independent sparse
rational samples.
