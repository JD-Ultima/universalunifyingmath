# BO. The full H6 curvature kernel on the flat positive observer

On the flat positive Z1 observer at $\lambda=1/5$, the complete
coefficient-curvature kernel in the real $H_6$ hook has dimension 6,834.
It contains 18 irreducible copies across 15 types. This describes the
possible values allowed by the linear pointwise curvature condition
arising from [BA](BA-hook-closure.md). General closed field extensions
and the nonlinear equations require further conditions.

The [exact verifier](../../computation/geometric-unity/BO_h6_full_kernel.py)
reconstructs the full stack and checks the adjacent
[extremal vector and real witness data](../../computation/geometric-unity/BO_h6_full_kernel_witnesses.json).
[Review79](../../reports/reviews/REVIEW79-GU-BO-2026-11-14.md) records the
independent checks.

## BO1. The stack and its 185 weight blocks

Use [BG](BG-h6-fixed-values.md)'s normalized trace-metric realization
$V=\operatorname{Sym}^2_0(\mathbb R^5)$ and extend scalars to $\mathbb C$.
On the 42,042-dimensional $V^*\otimes\Lambda^6 V^*$, stack the following
maps with the indicated targets:
$$
\mathrm{Alt}\colon\Lambda^7 V^*,\qquad
\mathrm{Ctr}\colon\Lambda^5 V^*,\qquad
\mathcal R\colon V^*\otimes\Lambda^8 V^*.
$$
The last map uses the coefficient covector action $-R_{ij}^{T}$ with
$R_{ij}P=-[[P_i,P_j],P]/4$. The contraction uses the full inverse Gram
matrix. These maps are equivariant under conjugation. The first two
jointly have rank $\binom{14}{7}+\binom{14}{5}=5434$ by the hook
decomposition in BA.

To work over rational coordinates, use the five-dimensional form $K$
with nonzero entries $K_{01}=K_{10}=K_{23}=K_{32}=2$ and $K_{44}=1$.
Take symmetric matrices $C=E_{ii}$ or $E_{ij}+E_{ji}$ with $i<j$
in lexicographic pair order, omitting $(0,1),(2,3),(4,4)$, then append
$E_{01}+E_{10}-4E_{44}$ and $E_{23}+E_{32}-4E_{44}$.
The tangent basis is $P=CK$ for these 14 matrices. If $U$ has columns
$e_0+ie_1,e_0-ie_1,e_2+ie_3,e_2-ie_3,e_4$, then $U^TU=K$.
Conjugation by $U$ maps $CK$ to $UCU^T$ and preserves the trace
metric and commutator. It therefore identifies the rational model
after complexification with BG's complexified tangent model.

The diagonal Cartan has vector weights $(1,0),(-1,0),(0,1),(0,-1),(0,0)$
on the five-dimensional space. In $V$ the weights add on the symmetric
pairs, with the two appended tangent matrices of weight zero. All
tensor slots are covariant, so their Cartan weights are the negatives
of the additive labels used to index the blocks. The stack preserves
these labels. There are 185 blocks with input columns, and the largest
has 1,306 columns. The verifier constructs each modular matrix directly
from the rational metric and commutators. The ranks over $\mathbb F_{65537}$
sum to 35,208. The rational denominators are invertible at this prime.
Nonzero modular minors give a characteristic-zero rank at least
35,208, hence kernel dimension at most $42042-35208=6834$.

## BO2. Characteristic-zero witnesses

In the rational model, take the root matrices
$F_1=E_{20}-E_{13}$ and $F_2=E_{34}-2E_{42}$. Their actual Cartan
weights are $(-1,1)$ and $(0,-1)$. With the signs of both Cartan axes
reversed they are the two simple positive roots. Their commutators give
the remaining two positive root spaces. Those four root matrices,
their opposites and the two Cartan matrices are independent and span
the ten-dimensional orthogonal Lie algebra of $K$.

The verifier checks each of the 18 exact rational vectors in the
certificate. Every vector has the stated weight for the reversed
Cartan, lies in the kernel of all three maps, and is annihilated
by both simple positive root actions on all seven slots. For each
weight it also checks independence of the supplied vectors. The kernel
is an invariant finite-dimensional representation of the complex
semisimple Lie algebra. Complete reducibility and the highest-weight
classification imply that the independent highest vectors supply
independent copies of the listed types. The Weyl dimension formula in
the reversed orthogonal weight convention is
$$
d(l,m)=\frac{(l-m+1)(l+m+2)(2l+3)(2m+1)}6,
\qquad l\ge m\ge0.
$$

| Weight | Copies | Dimension each |
|:---|---:|---:|
| $(7,3)$ | 1 | 1190 |
| $(7,2)$ | 1 | 935 |
| $(7,1)$ | 1 | 595 |
| $(6,3)$ | 1 | 770 |
| $(6,2)$ | 1 | 625 |
| $(6,1)$ | 2 | 405 |
| $(5,5)$ | 1 | 286 |
| $(5,2)$ | 2 | 390 |
| $(5,1)$ | 1 | 260 |
| $(5,0)$ | 1 | 91 |
| $(4,1)$ | 2 | 154 |
| $(3,2)$ | 1 | 105 |
| $(3,0)$ | 1 | 30 |
| $(2,1)$ | 1 | 35 |
| $(2,0)$ | 1 | 14 |

Their dimensions sum to 6,834, matching the upper bound. Complexifying
a real kernel preserves its dimension. The curvature rank on the
36,608-dimensional hook is consequently 29,774. Positive constant
observer scales multiply curvature by a common nonzero factor under
the natural identification, preserving the kernel and the decomposition.

## BO3. The grade-nine sign outside the contained module

The spherical types $(2,0),(3,0),(5,0)$ are BG's contained module.
[BN](BN-h6-invariant-grade-nine.md) excludes the entire module from
the full equation by a positive invariant projection $L_9$. The
remaining kernel contains 6,699 dimensions. The adjacent certificate
also gives a real 96-entry tensor $b$ in BG's physical coordinates. The
verifier checks both hooks and the full physical curvature map, and
computes
$$
\|b\|^2=1/2048,\qquad L_9(b)/\|b\|^2=-1/384.
$$
Thus the BN quadratic has both signs on the full real kernel. The
construction supplies a point value satisfying the linear curvature
condition. Full pointwise solutions and smooth field extensions remain
open.
