# CJ. A universal H6 gradient lower bound

On the real pure H6 branch of the normalized flat positive observer,
every full solution satisfies the strict pointwise estimate
$$
 \|\nabla B\|^2>\frac{1731}{25000}\,\|B\|^2 .
$$
The estimate combines the H6 stress divergence with the necessary
invariant quadratic equations. Without those invariant equations,
closure, scalarity and constant norm still give the lower bound
$19\|B\|^2/300$. At the [BS](BS-h6-certified-partial-value.md) unit
partial value the corresponding necessary squared gradient is
greater than $0.55699088$.

## CJ1. The actual eight-form and the stress

With the orientation and Clifford conventions of
[BA](BA-hook-closure.md), the actual degree-eight input
$C=B\Gamma$ is the form-slot Hodge dual up to a common sign.
It is a real hook. The closure $d_\nabla B=0$ is equivalent
to $\delta_\nabla C=0$, and the star is a parallel isometry.
In a positive orthonormal frame write
$$
 N=\|C\|^2=\|B\|^2,\quad
 E_{ij}=\langle C_i,C_j\rangle,\quad
 I_{ij}=\sum_b\langle\iota_i C_b,\iota_j C_b\rangle .
$$
For actual Clifford degree eight, the grade-one stress is
$$
 q_1(C)=-\tfrac23Ng+\tfrac83E-\tfrac43I .
$$
Indeed, the trace pairing with $\gamma_j$ in the reduced
Clifford formula works as in [CB](CB-h4-certified-partial-value.md).
The wedge part of $\gamma_jC_b$ has grade nine and reversal
sign $+1$. Its contraction part has grade seven and reversal
sign $-1$. This gives the same relative coefficients as for
degree four. Equivalently,
$I(C)=Ng-I(B)$ by the Hodge identity, and substitution into
the BS grade-one formula gives the displayed expression.
As $\operatorname{tr}I=8N$, scalarity is
$$
 S=2E-I+\tfrac37Ng=0 .
$$
A full field has constant $N>0$ by its positive grade-one baseline.

## CJ2. The double divergence

The [CH](CH-h4-second-order-stress.md) double-divergence calculation
uses only the alternating hook and the form-slot codifferential.
It applies in degree eight as well. For $J_{ji}=\nabla_jC_i$ set
$$
 X=\sum_{i,j}\langle J_{ji},J_{ij}\rangle,\qquad
 Y=\sum_{i,j,b}\langle\iota_iJ_{jb},\iota_jJ_{ib}\rangle .
$$
With $R_{ji}=[\nabla_j,\nabla_i]$ acting on all nine slots, put
$$
 H(C)=2\sum_{i,j}\langle C_i,(R_{ji}C)_j\rangle
       -\sum_{i,j,b}\langle\iota_iC_b,\iota_j(R_{ji}C)_b\rangle .
$$
Then
$$
 \operatorname{div}\operatorname{div} S-\tfrac37\Delta N
 = H(C)+2X-Y .
$$
The same orthogonal exchange and exterior norm identities as
in CH give $2X-Y\ge-3\|J\|^2$. Thus a closed field preserving
scalarity and constant norm through second order satisfies
$3\|\nabla B\|^2\ge H(C)$.

## CJ3. Exact positivity and the full gradient estimate

Use the 15 curvature types and 18 exact highest vectors of
[BO](BO-h6-full-curvature-kernel.md) in the common compact real
directions of [BQ](BQ-h6-invariant-moments.md). The curvature form
$H(*B)$ is $SO(5)$-invariant. Its restriction to each multiplicity
space is defined by its exact Gram and polarized source matrices.
The [verifier](../../computation/geometric-unity/CJ_h6_stress_gradient.py)
applies the full ambient curvature matrices to the actual
eight-forms. It uses the rational form-slot star that omits
the scalar volume, then multiplies the quadratic by the
metric determinant to recover the actual physical form.
The exact Gram identity checks this normalization.

On every multiplicity block, the following two rational forms
have strictly positive leading principal minors:
$$
 H-\frac{19}{100}N,\qquad
 H+\frac{2303}{25000}L_5-\frac{8647}{50000}(384L_9)
      -\frac{5193}{25000}N .
$$
Therefore the first form is positive on the full 6,834-dimensional
curvature kernel. On $L_5=L_9=0$, the second gives
$H>(5193/25000)N$. The double divergence gives the announced
gradient estimates, since $5193/(3\cdot25000)=1731/25000$.

The verifier also reconstructs all 45 upper-triangular coefficients
of $H$ on the BS affine slice. Its center is approximately
1.670972681, with exact total variation below $3.763\cdot10^{-8}$.
The certified box gives $H>1.67097264$, so the stronger stated
gradient estimate applies to any lift preserving scalarity and
constant norm through second order there. This estimate concerns
the normalized positive observer. Dilation and homothety transport
the identities with the metric scaling of the squared gradient.

    python -X utf8 computation/geometric-unity/CJ_h6_stress_gradient.py

[Review99](../../reports/reviews/REVIEW99-GU-CJ-2026-11-14.md)
records the independent matrix assembly. The inequalities require
nontrivial variation and leave general pointwise and field
existence open. External peer review remains outstanding.
