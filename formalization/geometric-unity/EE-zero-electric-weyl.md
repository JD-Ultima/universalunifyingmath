# EE. Vanishing electric Weyl curvature

At a nonzero rank-one Ricci point, the obstruction now covers every Weyl tensor whose electric part relative to the Ricci eigenline vanishes. The mixed part may be arbitrary. Three finite fibre configurations have an exact compatibility cover and full coefficient rank at those same configurations. Arbitrary higher metric jets are allowed.

This complements [CC](CC-diagonal-weyl.md), where the mixed part vanishes and the electric part is arbitrary. [DD](DD-mixed-weyl.md) retains the shorter one-direction mixed certificate. Neither condition covers a general Weyl tensor with both parts nonzero.

Code: `computation/geometric-unity/EE_zero_electric_weyl.py`, `EE_rank_blocks.py` and `EE_rank_certificates.json`; output: `EE_zero_electric_weyl.out`; record: `metrics/geometric-unity/EE-zero-electric-weyl.yaml`. The calculation uses U's reviewed assembly and the earlier referee curvature engine. [Independent internal Review29](../../reports/reviews/REVIEW29-GU-EE-2026-11-02.md) accepts EE1–EE6 and reproduces the portable computation. This is internal review, not external peer review.

## Intrinsic condition and normal form

Let $p$ be a unit vector spanning the Ricci eigenline at the base point, with metric-dual covector $p^\flat$. Assume
$$\mathrm{Ric}=\rho\,p^\flat\otimes p^\flat,\qquad \rho\ne0,
\qquad E(u,v):=W(p,u,p,v)=0\quad(u,v\perp p).\tag{EE1}$$
The sign of $p$ does not affect this condition. The electric tensor here is a Weyl contraction; no vanishing of the base Ricci contraction is assumed.

**Lemma EE1 (pure-mixed normal form).** Condition $E=0$ is equivalent to the existence of a Ricci-adapted orthonormal frame in which the only independent Weyl entries are
$$W_{0123}=\alpha,\qquad W_{0231}=\beta,\qquad
W_{0312}=-(\alpha+\beta),\tag{EE2}$$
with the usual curvature symmetries.

*Proof.* Choose $e_0=p$ and orient the perpendicular three-space. Define
$$B_{ij}=\tfrac12\epsilon_{jkl}W_{0ikl},\qquad i,j,k,l=1,2,3.$$
The zero Weyl Ricci contraction gives
$$\mathrm{Ric}^{W}_{01}=B_{23}-B_{32},\quad
\mathrm{Ric}^{W}_{02}=B_{31}-B_{13},\quad
\mathrm{Ric}^{W}_{03}=B_{12}-B_{21}.$$
Thus $B$ is symmetric. Bianchi gives $\operatorname{tr}B=0$. Since $E=0$, the spatial three-dimensional curvature tensor has Ricci tensor $-E=0$ and hence is zero. Its curvature-to-Ricci contraction is an isomorphism; the code independently checks its six-dimensional matrix determinant is $-2$. The remaining mixed tensor is recovered exactly from $B$ by $W_{0ikl}=\epsilon_{jkl}B_{ij}$.

A real symmetric tracefree matrix is diagonalizable by an $SO(3)$ change of perpendicular frame; an eigenvector sign can be adjusted to preserve orientation. Its eigenvalues are $(\alpha,\beta,-\alpha-\beta)$, giving (EE2). Conversely the tensor in (EE2) has zero Ricci contraction, satisfies Bianchi and has no electric or spatial component. The magnetic extraction, symmetry equations and trace equation are verified symbolically for a general nine-entry matrix, rather than only diagonal examples. ∎

Let $H_A,H_B$ be the pure-mixed Weyl tensors in (EE2) at $(\alpha,\beta)=(1,0),(0,1)$. With $R_0=\tfrac12(e^0\otimes e^0-\delta/6)\owedge\delta$, the actual smooth metric has normal-coordinate expansion
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\rho R_0+\alpha H_A+\beta H_B)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\tag{EE3}$$
where the symmetric cubic and quartic terms are arbitrary and the smooth remainder has zero derivatives through order four at $0$. Positivity holds on some neighbourhood, without a uniform radius over these coefficients.

Use Z1 at $\lambda=1/5$, $c_H,c_V>0$. The full field-equation statement has S.0's compact $(14,0)$ real form, gauge algebra $\mathfrak u(128)$, invariant odd field, slaved vector even field, one fixed reconstructed Shiab member with coefficients $(x,y)\ne(0,0)$, $\kappa\ne0$ and $\det M\ne0$. U's naturality and parameter maps normalize
$$\rho=\sigma\in\{1,-1\},\quad c_H=2,\quad c_V=3,\qquad
a=\alpha/|\rho|,\quad b=\beta/|\rho|.\tag{EE4}$$
The normalized cubic and quartic terms $Q_3/|\rho|^{3/2}$ and $Q_4/\rho^2$ remain arbitrary. DD is the subcase $b=-a$.

## Degree bounds and compatibility

First set the higher terms to zero and calculate on the auxiliary quadratic metric.

**Lemma EE2 (entry degrees).** At a fixed positive fibre matrix, the relaxed coefficient entries have total degree at most two in $(a,b)$, and the RHS entries total degree at most three, with parameter-independent denominators.

*Proof.* The quadratic base coefficient is affine in $(a,b)$. BB1's jet argument applies to any such curvature: the base connection's linear coefficient has degree at most one, its quadratic coefficient is zero, and its cubic coefficient has degree at most two. The observer metric and inverse at the point are fixed; its first derivatives have degree at most one and its second and third derivatives degree at most two. Curvature has degree at most two and its covariant derivative at most three. These are the tensors entering the relaxed system. The degree bounds are established before interpolation. ∎

The triangular nodes $(i,j)$, $i,j\ge0$, $i+j\le3$, are unisolvent for total degree three. The first six with $i+j\le2$ are unisolvent for degree two. Their absolute Vandermonde determinants in the code's order are $576$ and $4$. These exact interpolation spaces, rather than independent univariate fits, certify the polynomial identities below.

**Proposition EE3 (three-fibre compatibility cover).** At least one of
$$y_1=\operatorname{diag}(2,3,4,5),\quad
y_2=\operatorname{diag}(2,3,5,7),\quad
y_3=\operatorname{diag}(2,4,5,9)$$
has incompatible selected rows for every pair $(a,b)$ and either sign.

*Certificate.* At each of these stated fibres, rows $(168,337,492)$ are $(a,b,d)=(0,13,0),(1,13,1),(2,13,2)$ and have exactly the form
$$L_y(n)+q_i(a,b)\xi_{13}=r_i(a,b),\qquad i=0,1,2.\tag{EE5}$$
The common metric combination $L_y$ is independent of curvature. Every unwanted coefficient is zero at the six degree-two nodes, hence identically zero there as a function of curvature. The four additional nodes check the reconstructed coefficients. The ten-node RHS interpolation has zero cubic coefficients as an exact calculation within EE2's cubic bound.

Let $P_{\sigma,j}$ be the primitive numerator of the determinant of rows $(q_i,1,r_i)$ at $y_j$. The code prints every coefficient and RHS polynomial and all six quartic numerators. Exact Gröbner computation gives, for each sign,
$$\langle P_{\sigma,1},P_{\sigma,2},P_{\sigma,3}\rangle
=\langle1\rangle\quad\text{in }\mathbb Q[a,b].\tag{EE6}$$
Thus these compatibility polynomials have no common complex zero and in particular no common real zero. No parameter search, root sampling or fourth configuration is needed. ∎

## Same-witness rank and real completion sign

Take the 722-row cubic-invariant subsystem $M$ of BB3/DD3. At a diagonal fibre its four independent kernel columns are
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\qquad i=0,1,2,3.\tag{EE7}$$
The helper checks $MK=0$ at all six coefficient nodes, proving the identity by EE2.

**Proposition EE4 (rank at every compatibility fibre).** For both signs, every real $(a,b)$ and each of the three numerical fibres, the invariant subsystem has rank 38 and kernel $\operatorname{span}K$. The full coefficient matrix has rank 42.

*Constant block certificates.* The quadratic pure-mixed background is fixed by the even-sign subgroup $s_0s_1s_2s_3=1$. An odd coordinate sign flip reverses its Weyl part, so the full independent sign group cannot be used. As in DD2, the even-sign characters split the unknowns into a 12-column block, three complementary off-diagonal 6-column blocks and four horizontal 3-column blocks. Six exact evaluations certify the polynomial block support.

For the 12-column and 6-column blocks, constant rational left inverses satisfy $TM(a,b)=I$. For each horizontal block, two constant row combinations have rank two; (EE7) is its remaining kernel. All row indices and rational coefficients are in the compact JSON certificate. The helper verifies these identities against the freshly computed matrices, giving
$$\operatorname{rank}M=12+3\cdot6+4\cdot2=38$$
for every parameter. These are polynomial identities, rather than ranks observed at sample values.

*Real completion positivity.* Append rows $(i,O,O)$ for $i=0,1,2,3$ and $O=5,6,7$, representing $y_{01},y_{02},y_{03}$. Their restrictions to $K$ are negative diagonal matrices with positive rational prefactors times the following factors:

| fibre | $J_5$ factor | $J_6$ factor | $J_7$ factor |
|---|---|---|---|
| $y_1$ | $10+27a^2$ | $5+24b^2$ | $2+15(a+b)^2$ |
| $y_2$ | $35+54a^2$ | $7+30b^2$ | $5+42(a+b)^2$ |
| $y_3$ | $5+8a^2$ | $2+5b^2$ | $10+81(a+b)^2$ |

Every factor is strictly positive for real parameters. The JSON encodes each actual diagonal polynomial as a positive constant plus a positive semidefinite rational quadratic form, with positive diagonal scaling coefficients. The helper checks the exact polynomial identity, the positive constant, and the two-dimensional principal-minor inequalities proving semidefiniteness. Therefore each completion is invertible on $K$, and any one gives rank 42 at that same fibre.

This step is a real positivity argument. The complex compatibility unit ideals in (EE6) do not assert complex nonvanishing of the completion factors. ∎

## Higher jets and the local equation

**Proposition EE5 (arbitrary higher jets).** EE3–EE4 hold for all cubic, quartic and smooth higher terms in (EE3).

*Proof.* DD3 explicitly audits the curvature-independent version of BB3. For any fixed quadratic curvature, the cubic curvature variation at a diagonal fibre has exactly three horizontal slots and one off-diagonal vertical slot; its Ricci variation has horizontal/off-diagonal-vertical support. The diagonal-vertical mixed variation is a base gradient, eliminating that vertical slot by commutation of derivatives. This fixes $M$ and the selected three coefficient rows.

The simultaneous inversion $(x,y)\mapsto(-x,y)$ fixes every quadratic metric, reverses $Q_3$ and leaves the quartic jet fixed. Affine cubic dependence and the two-horizontal parity of the selected RHS rows eliminate their cubic variation. W2's quartic gradient calculation eliminates the remaining principal third-observer-jet variation. The smooth zero-four-jet remainder affects neither relevant observer jet. Thus the compatibility determinants remain fixed.

For every completion row $(i,O,O)$, only the coefficient of $\xi_O$ can vary with the cubic tensor. It annihilates $K$, since $\xi_V=0$ on (EE7). Quartics change no coefficient. Hence the strictly invertible restrictions in EE4 persist for every higher jet. The separate-coordinate even-sign symmetry is used only for the auxiliary quadratic rank calculation; arbitrary higher jets need not have it. ∎

**Theorem EE6.** At a point satisfying (EE1), the reconstructed Z1 relaxed system at $\lambda=1/5$, $c_H,c_V>0$, is consistent on no nonempty fibre-open patch. Under the stated S.0 assumptions, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* EE1 puts the actual curvature in the pure-mixed normal form. EE3 supplies an incompatible numerical fibre for its actual normalized real parameters. EE4–EE5 give coefficient rank 42 at that same configuration for the actual higher metric jets. Thus the augmented rank is 43. S4's rational-minor argument excludes every fibre-open consistency patch, and U's reviewed parameter maps cover every original nonzero $\rho$ and positive observer scales.

For the full equation, fix one actual nonzero augmented 43-minor at its positive coordinate fibre matrix. It remains nonzero over a base neighbourhood. Shrink that neighbourhood to preserve nonzero curvature, using nonzero Ricci at the initial point. S5 excludes Einstein pieces, while Q5 on the dense variable patch with $H\ne0$ would give relaxed consistency, contradicting the fibrewise obstruction as in U5. The continued object is the fixed full-system minor. The zero-electric hypothesis is needed only at the initial point, not throughout that neighbourhood. ∎

## Reproduction and remaining boundary

Run `python -X utf8 -I computation/geometric-unity/EE_zero_electric_weyl.py` from the repository root with the pinned dependencies. The portable default freshly computes all 60 exact node systems once in memory and passes the first six coefficient matrices per case to the rank helper. The JSON contains compact rational row combinations and real positivity certificates, with no raw node matrices or off-repository cache dependency.

The recorded fresh run passes 32 main checks and 72 rank-helper checks in 114.0 s; Review29 independently reproduces them in 113.7 s. The compact certificate file is 49,791 bytes.

This covers all zero-electric Weyl tensors relative to a nonzero Ricci eigenline, including DD's simpler subcase. CC covers the separate zero-mixed condition. [FF](FF-joint-weyl.md) separately treats the joint commuting/common-kernel family. Other cases with both parts nonzero remain open at nonzero Ricci; [GG](GG-ricci-zero-weyl.md) separately treats the Ricci-zero case with singular chiral blocks. The field statement remains within the reconstructed sector and the explicit S.0 assumptions.
