# BN. An invariant grade-nine obstruction on the 135-dimensional H6 module

The entire real 135-dimensional curvature submodule found in
[BG](BG-h6-fixed-values.md) is excluded from the full equation on the
flat positive Z1 pure H6Γ branch. The result covers arbitrary mixtures
of the three harmonic summands without imposing symmetry on the value
or its field extension.

The [exact helper](../../computation/geometric-unity/BN_h6_invariant_grade9.py)
reconstructs the receiver and all three calibrations from the public
tangent basis and seeds. [Review78](../../reports/reviews/REVIEW78-GU-BN-2026-11-14.md)
records the independent checks.

## BN1. An invariant receiver

Retain [BM](BM-mixed-grade-nine-obstruction.md)'s positive signature
$(14,0)$, compact real form, fixed parallel real $W_A$, constant
$\kappa\ne0$, and
$$
T_o=i(3\kappa/104)\Phi^1,\qquad T_e=iB\Gamma.
$$
The flat observer has $\lambda=1/5$. Work first at $c_H=2,c_V=1$, with
BG's trace-metric identification of the tangent space with
$V=\operatorname{Sym}^2_0(\mathbb R^5)$. All exterior norms and contractions
use this positive metric.

Write $Q_a=5\Phi_a$ for BG's integer tangent matrices, and set
$F=\sum_a Q_a e^a$. Here $F$ is five times the tautological inclusion of
$V$ into matrices, with conjugation action. On matrix-valued forms use
matrix multiplication in the displayed order and the ordinary exterior
product. Define the covariant receiver by
$$
\Omega_9(e_a,\cdot)=10^{-10}\operatorname{tr}(Q_a F^9).
$$
Coefficients are taken in increasing-index exterior blades. In each
increasing-index slot, $F^9$ contains the signed sum over all permutations
with no factorial divisor. Equivariance of $F$, invariance of the trace
under conjugation, and the covariant first argument make $\Omega_9$
$SO(5)$-invariant. It has 1,143 nonzero coordinates, and the raw
coefficients have positive common divisor $10^{10}$. Its norm squared is
$2835/2097152$.

Let $Q_9(iB\Gamma)/i$ denote degree nine of the declared dualized
quadratic. Pair it with $\Omega_9$:
$$
L_9(B)=\left\langle\Omega_9,Q_9(iB\Gamma)/i\right\rangle.
$$
For diagonal Gram $g$, the $(a,I)$ summand has inverse weight
$1/(g_a\prod_{j\in I}g_j)$. As in BM, the helper computes $B$ times the
coordinate Clifford volume on the right and multiplies the resulting
quadratic coefficients by $\det g$ to restore the physical volume
normalization.

## BN2. Positivity on the entire contained module

Use the orthogonal decomposition $B=B_2+B_3+B_5$ of BG's module into
irreducible spherical degrees 2, 3, and 5, with dimensions 14, 30, and
91. The Clifford quadratic is equivariant, so polarization of $L_9$
defines a symmetric invariant bilinear form. Its cross terms between
distinct irreducibles vanish. On an individual irreducible, associate
to this bilinear form its self-adjoint operator with respect to the
positive metric. Its real eigenspaces are invariant, so it is a scalar.
It therefore suffices to calibrate one nonzero tensor in each summand.

The degree-two seed is BG's third fixed-kernel tensor. For the other
two, let $b$ be BC's explicit closed generator and $C$ the positive
Casimir of [BF](BF-h6-curvature-module.md). The seeds are
$(40b-Cb)/22$ and $(Cb-18b)/22$. The exact calibrations are

| Degree | Seed norm squared | $L_9$ on the seed | $L_9$ divided by norm squared |
|---:|---:|---:|---:|
| 2 | $18225/16$ | $26325/2048$ | $13/1152$ |
| 3 | $313875/88$ | $25875/1408$ | $23/4464$ |
| 5 | $455625/45056$ | $354375/5767168$ | $7/1152$ |

Hence, for every real value in the whole module,
$$
\boxed{35712 L_9(B)=403\|B_2\|^2+184\|B_3\|^2+217\|B_5\|^2.}
$$

## BN3. Field consequence

The baseline and curvature contributions have grade one. The entire
variable derivative term of the pure hook has grade seven by
[BA](BA-hook-closure.md). Thus a full solution requires
$Q_9(iB\Gamma)=0$ at each point. If a value belongs to the contained
135-dimensional module, the positivity identity forces $B=0$ there.
That value fails the grade-one baseline
$$\ell=\frac{15}{2c_V}+\frac{15\kappa^2}{52}>0.$$

The positive constant observer scales use the same natural
identifications as [BG](BG-h6-fixed-values.md) and [BJ](BJ-h6-orbit-averages.md).
In the homothetic frame the physical quadratic stress gains a common
positive factor, preserving its vanishing condition. The baseline stays
positive for every $c_V>0$.

The obstruction applies even when a candidate field leaves the
135-dimensional module elsewhere in its domain. [BO](BO-h6-full-curvature-kernel.md)
classifies the full linear curvature kernel and verifies that $L_9$
has both signs there. General nonlinear feasibility and corresponding
closed fields remain open on the stated branch. External peer review
remains outstanding.
