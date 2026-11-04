# BD. The affine-invariant pure H4 branch is empty

There is no affine-invariant full bosonic field on the pure $H_4$ branch $t=\kappa/24,u=0$ over the flat positive Z1 observer at $\lambda=1/5$, for any positive observer scales. [BC1](BC-affine-hooks.md#bc1-complete-linear-reduction-within-the-invariant-ansatz) leaves six derivative-compatible amplitudes. Eighteen necessary higher-grade quadratic equations now force all six to vanish; the grade-one baseline then gives a contradiction.

This excludes the stated homogeneous ansatz, not arbitrary smooth $H_4$ fields. Evidence: [helper](../../computation/geometric-unity/BD_h4_affine.py), [output](../../computation/geometric-unity/BD_h4_affine.out) and [record](../../metrics/geometric-unity/BD-affine-h4.yaml). Internal [Review64](../../reports/reviews/REVIEW64-GU-BD-2026-11-04.md) and [Review65](../../reports/reviews/REVIEW65-GU-BD-2026-11-04.md) accept BD1–BD3 and disclose their contributors' roles. External peer review remains outstanding.

## Setting

Retain BC's true dualized equation, fixed parallel real $W_A$, compact real $(14,0)/\mathfrak u(128)$, real smooth fields and constant real $\kappa\ne0$. The base is flat and positive, $\lambda=1/5$, and $c_H,c_V>0$ are constant. The even field is the restriction of a tensor invariant under
$$G=GL^+(4,\mathbb R)\ltimes\mathbb R^4,\qquad Y=G/SO(4).$$
The spin/Clifford identification is included in this assumption. Invariance does not mean parallelism.

Fix $T_o=i\kappa\Phi^1/24$. R/V leave precisely $T_e=iC$, with $C$ a real $H_4$ tensor and no other even component. Its actual Clifford grade is four. The vector block is regular. [BA](BA-hook-closure.md) makes the full derivative equation equivalent to closure of the dual $H_{10}$ coefficient, or to $\delta_\nabla C=0$. Metric Hodge duality is used for this derivative step only; the quadratic calculation below acts directly on $C_4$.

## BD1. Six amplitudes and eighteen necessary equations

Normalize $c_H=2,c_V=1$ using BC's affine model. The helper calls BC's public construction as a library, without running its main. It reconstructs all 20 invariant vector-valued four-forms, imposes both weighted hook conditions to obtain 13 hooks, and computes the full all-slot Nomizu derivative. Its rank-seven map leaves six independent real tensors.

Let $K_0,\ldots,K_5$ be the deterministic raw closure combinations printed by the helper. Replace the final direction by $K_5-K_4$ and normalize each resulting tensor primitively. Denote this basis by $P_j$ and write $C=\sum_{j=0}^5x_jP_j$. Every basis tensor is checked against the six infinitesimal isotropy generators, both hook traces and full dual closure. No derivative slot is omitted.

The full odd equation has derivative grade three, baseline grade one, and quadratic grades one, five and nine. Thus
$$[Q(iC)]_5=[Q(iC)]_9=0$$
are necessary equations with no derivative cancellation. The helper uses the complete differentiated cubic [BC2.2](BC-affine-hooks.md). In the non-unit covector basis, Clifford squares are inverse Gram entries and the external output contributes its Gram factor. The common compact-real phase is removed by dividing by $i$.

Eighteen genuine output coordinates are fixed in the source and printed in the output. Each entire quadratic coefficient is obtained by expanding the three gradient terms into ordered Clifford products. For a triple mask $A$, input masks $I,J$ and target mask $T$, the equality $A\mathbin{\mathrm{xor}}I\mathbin{\mathrm{xor}}J=T$ fixes the matching second blade. Exact summation retains every contribution. Three individual-input comparisons and one mixed-amplitude polarization comparison with BC's complete cubic evaluator check all eighteen coordinates: 72 coefficient equalities, including nonzero grade-five and grade-nine controls. The mixed check independently tests the assembly of cross terms.

Set $y_{ij}=x_ix_j$ for $i\le j$. For a symmetric quadratic matrix, mixed coefficients are $2Q_{ij}$, not $Q_{ij}$. The selected equations are $My=0$, where $M$ is $18\times21$ in the order
$$00,01,02,03,04,05,11,12,13,14,15,22,23,24,25,33,34,35,44,45,55.$$
Exact rational elimination gives rank 18 and kernel dimension three, with free coordinates $y_{44},y_{45},y_{55}$. The output records all rows, the nullspace and a nonzero $18\times18$ minor. These are eighteen necessary equations; no claim about the rank or completeness of all higher outputs is required.

## BD2. The actual rank-one monomials force zero

Rank 18 alone does not exclude nonzero amplitudes. Use the additional identities satisfied by actual products $y_{ij}=x_ix_j$.

The computed kernel has $y_{3j}=0$ for all $j\ne3$.

* If $x_3\ne0$, every other amplitude is zero. All three free monomials then vanish, so the whole monomial vector is zero, contradicting $y_{33}=x_3^2\ne0$. Hence $x_3=0$.
* If $x_4=0$, then $y_{44}=y_{45}=0$. The $y_{33}=0$ relation becomes
$$-\frac{150231665688711}{1931741174969578496}x_5^2=0.$$
Thus $x_5=0$, all free monomials vanish and every amplitude is zero.
* Otherwise set $z=x_5/x_4$ and divide every monomial relation by $x_4^2$. The relation for $y_{33}=0$ gives the primitive polynomial
$$F(z)=150231665688711z^2+127368305168964600z-1375492211609179736800.$$
Write $R_{ij}(z)$ for the normalized relation $y_{ij}/x_4^2$ read from the three nullspace columns. The actual rank-one identity $y_{00}y_{11}-y_{01}^2=0$ gives
$$G(z)=\operatorname{primitive}\big(R_{00}(z)R_{11}(z)-R_{01}(z)^2\big).$$
Its explicit quartic coefficients are in the output. Exact arithmetic gives $\gcd_{\mathbb Q[z]}(F,G)=1$. The helper prints integer polynomials $U,V$ and directly verifies $UF+VG=D$ for a nonzero integer $D$. This excludes the remaining chart.

Therefore the selected higher-grade equations have only $x=0$. No positive-semidefinite relaxation or sampled root search substitutes for the rank-one identities.

## BD3. Full-field contradiction and positive scales

At $C=0$, the remaining grade-one baseline is
$$
L=\frac{15}{2c_V}+\kappa t+312t^2
=\frac{15}{2c_V}+\frac7{12}\kappa^2>0.
$$
It cannot vanish. Both signs of $\kappa$ are included.

Arbitrary positive scales are recovered by
$$
g(c_H,c_V)=c_VF_s^*g(2,1),\qquad F_s(x,y)=(sx,y),\quad s=\sqrt{c_H/(2c_V)}.
$$
The dilation conjugates the affine action and preserves the invariant ansatz. Constant homothety preserves the connection and the zero hook/closure conditions. For a fixed five-covariant $C_4$, each quadratic term contains eleven Clifford covectors. A grade-$r$ output requires $(11-r)/2$ contractions; the external Gram factor contributes one power of the scale. Consequently, under $g\mapsto hg$,
$$Q_r/i\mapsto h^{(r-9)/2}Q_r/i$$
in fixed covariant coordinates. Thus grade five scales by $h^{-2}$ and grade nine by 1. Their simultaneous zero set is unchanged, so BD2 still forces zero. The actual baseline $L$ remains strictly positive. This proves the result for every $c_H,c_V>0$.

## Reproduction and boundary

Run

    python -X utf8 -I computation/geometric-unity/BD_h4_affine.py

The author run completes in 1.829 seconds with pinned Python 3.12/SymPy 1.14.0. It reconstructs the invariant/closure basis, performs the four complete-gradient comparisons above, generates all 21 coefficients of each selected quadratic and verifies every rank, chart and integer Bezout identity. Assertions require elimination success. No private input, saved basis receipt or earlier helper main is used. The calculation reuses BC's tensor and Clifford library; these convention checks are not a second independent algebra engine.

The result closes BC's six-dimensional nonlinear gap only within the affine-invariant pure $H_4$ ansatz. General smooth $H_4$ fields, curved bases, other $\lambda$, other hook copies and coincident kernels remain outside it. It supplies no field construction, literature-priority claim or physical conclusion.
