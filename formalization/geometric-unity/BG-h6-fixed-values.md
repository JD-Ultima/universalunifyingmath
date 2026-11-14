# BG. The fixed-value H6 curvature kernel and its nonlinear obstruction

On the flat positive Z1 model at $\lambda=1/5$, the complete $SO(4)$-fixed
$H_6$ coefficient-curvature kernel has dimension three. Its $SO(5)$ orbit
span has dimension 135. No value fixed by this $SO(4)$, or a conjugate
stabilizer, satisfies the pure $H_6\Gamma$ full field equation. This excludes
such a value at even one point of a solution; amplitudes need not be constant.
It does not exclude arbitrary values in the 135-dimensional span.

Evidence: [curvature helper](../../computation/geometric-unity/BG_h6_curvature.py),
[complete tensor certificate](../../computation/geometric-unity/BG_h6_curvature_kernel.json),
[moment helper](../../computation/geometric-unity/BG_h6_moments.py) and their
[curvature](../../computation/geometric-unity/BG_h6_curvature.out) and
[moment](../../computation/geometric-unity/BG_h6_moments.out) outputs.
Internal [Review70](../../reports/reviews/REVIEW70-GU-BG-BH-2026-11-05.md) and [Review71](../../reports/reviews/REVIEW71-GU-BG-BH-2026-11-05.md) accept BG1–BG3 and disclose shared contributions and libraries. External peer review remains outstanding.

## Setting

Use [BF's](BF-h6-curvature-module.md) flat positive four-base, reconstructed
Z1 at $\lambda=1/5$, positive constant scales, $K=SO(5)$ and $L=SO(4)$ fixing
the last unit vector. At $c_H=2,c_V=1$, the tangent is
$V=\operatorname{Sym}^2_0(\mathbb R^5)$ with BC's Gram-orthogonal basis.
All slots of the real hook
$B\in H_6\subset V^*\otimes\Lambda^6V^*$ are covariant.

For the field exclusion, retain [BA2](BA-hook-closure.md): the true dualized
equation, fixed parallel real $W_A$, signature $(14,0)$, compact
$\mathfrak u(128)$, real smooth fields, constant $\kappa\ne0$ and
$$T_o=i(3\kappa/104)\Phi^1,\qquad T_e=iB\Gamma.$$
The vector block is regular and all other even blocks are nonsingular.
The geometric kernel statement itself needs no operator or $\kappa$.

## BG1. The complete fixed-value curvature map

The map is coefficient curvature, rather than curvature on every tensor slot:
$$\mathcal R B=F^{T^*Y}\wedge B,\qquad
(\mathcal RB)_a=-\sum_{i<j,d}R_{ij}{}^d{}_a\,e^i\wedge e^j\wedge B_d.$$
The symmetric-space isometry $\Phi$ gives the full ambient curvature
$$R_{ij}\Phi_a=-\tfrac14[[\Phi_i,\Phi_j],\Phi_a].$$
The helper calibrates $Q_a=5\Phi_a$ against the Gram form and all six
stabilizer generators, then checks metric skewness, the first Bianchi
identity and $\operatorname{Ric}=-5g/4$. Horizontal, mixed and vertical
channels are all retained.

[BC's complete invariant construction](BC-affine-hooks.md) gives 30 raw
degree-six coefficient tensors. The two hook constraints have joint rank
seven, leaving $\dim(H_6)^L=23$. Hodge duality sends the raw invariant basis
isomorphically onto the invariant degree-eight target, also of dimension
30. A nonzero $30\times30$ coordinate minor certifies the output coordinates.
No closure condition is imposed on the domain.

The resulting $30\times23$ matrix has rank 20. The recorded
$20\times20$ minor is
$$-3971266945655899580006400000\ne0.$$
Three independent reconstructed nullvectors $K_0,K_1,K_2$ have zero
**full sparse** curvature output, not just zero selected coordinates.
The lower rank bound and these nullvectors therefore prove
$$\boxed{\operatorname{rank}(\mathcal R|_{(H_6)^L})=20,
\qquad \dim(\ker\mathcal R)^L=3.}$$
The adjacent JSON records all 738, 1,358 and 878 coefficients of the three
tensors. The helper reconstructs and compares this certificate without
rewriting it.

In this basis BF's positive Casimir is
$$[\mathcal C]=\begin{pmatrix}352/9&76/9&0\\20/9&170/9&0\\0&0&10\end{pmatrix},
\qquad \det(zI-\mathcal C)=(z-10)(z-18)(z-40).$$
Every component of these identities is checked. The two known BF
projections are $v_{18},v_{40}$, and set $v_{10}=K_2$. Each is real,
nonzero, $L$-fixed and in the curvature kernel. BF's sphere
matrix-coefficient argument applies to each: eigenvalues 10, 18 and 40
give harmonic degrees two, three and five, of dimensions 14, 30 and 91.
Thus
$$\boxed{\dim\operatorname{span}_\mathbb R K\cdot(\ker\mathcal R)^L=135,
\qquad \dim\ker\mathcal R\ge135.}$$
Equivalently, the full map has rank at most $36608-135=36473$.
This is not its exact full rank: further components without $L$-fixed
vectors have not been classified. The new $v_{10}$ has no demonstrated
covariantly closed field lift.

## BG2. A positive quadratic obstruction on all three fixed values

Write $B=A v_{18}+B'v_{40}+Z v_{10}$. The actual compact field is
$iB\Gamma$, of Clifford degree eight. For its quadratic term divided by
$i$, let $q_H,q_S,q_I$ be the metric-relative grade-one eigenvalues on
the multiplicity-one $L$-modules $H_4$, $\operatorname{Sym}^2_0(\mathbb R^4)$
and $I_1$. Complete quadratic polarization, including every mixed term,
gives
$$\begin{aligned}
q_I-q_H
&=\frac{1125}{45056}\big(54272A^2+960AB'+505(B')^2\big)
  +\frac{28575}{32}Z^2\\
&=\frac{1125}{45056}\left[
54272\left(A+\frac{15B'}{1696}\right)^2+
\frac{26540}{53}(B')^2\right]+\frac{28575}{32}Z^2.
\end{aligned}$$
Every square coefficient is positive. Therefore scalar $Q_1$ forces
$A=B'=Z=0$ over the reals. The direct norm and full trace also agree:
$$N=\frac{313875}{88}A^2+\frac{455625}{45056}(B')^2+
\frac{18225}{16}Z^2,\qquad 4q_H+9q_S+q_I=-\frac{52}{3}N.$$
The helper checks the sum of squares, positive principal minors, exact
phase/volume conversion and 12 targeted comparisons with the separate
nested Clifford implementation. The two engines share the declared cubic
formula; they are not independent mathematical derivations of it.

## BG3. Consequence for smooth full fields

Any full solution on this pure branch obeys BA's universal derivative
identity and hence $d_\nabla B=0$ and $\mathcal RB=0$. BA differentiates
every tensor slot and both hook identities: the derivative has only grade
seven. There is no grade-one derivative even when the amplitudes vary.
The grade-one equation is
$$\ell g+Q_1=0,\qquad \ell=\frac{15}{2c_V}+\frac{15}{52}\kappa^2>0.$$

At any point where the value is $L$-fixed, BG1 places it in the complete
three-space. BG2 forces that value to vanish; then $Q_1=0$ contradicts
$\ell>0$. Thus a full solution cannot take an $L$-fixed value at any point.
This applies to values fixed by any conjugate $kLk^{-1}$ as well, since
$\mathcal R$, $Q_1$ and the scalar baseline are $K$-equivariant.
It requires no field reconstruction or constancy of the amplitudes.

BC's affine-conjugating dilation reduces arbitrary positive $c_H,c_V$ to
a positive homothety of the normalized model. Under the induced isometry
of orthonormal tangent representations, curvature is multiplied by a
nonzero scalar, preserving its kernel; the algebraic Clifford quadratic
and its scalarity condition are natural. Hence the kernel dimension,
positive obstruction and field exclusion hold for all positive scales.
The positive baseline retains the displayed $c_V$ dependence.

The complete ambient curvature kernel, other parameters, curved bases
and coincident branches are outside BG's scope. [BH](BH-h6-degree2.md)
treats the degree-two summand, and [BN](BN-h6-invariant-grade-nine.md)
excludes the whole contained module from the full equation on this
branch.
