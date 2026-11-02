# FF. A joint electric and mixed Weyl plane

This result allows both electric and mixed Weyl curvature to be nonzero. It treats a two-parameter simultaneous frame orbit: the electric and magnetic Weyl matrices commute and have a common nonzero kernel vector. That common-kernel condition is essential; arbitrary commuting pairs and arbitrary Weyl tensors are not included.

Three finite fibre configurations give an exact compatibility cover. Constant rational rank certificates and an invertible completion hold at those same configurations, and the curvature-independent support argument transfers the obstruction to arbitrary higher metric jets.

Code: `computation/geometric-unity/FF_joint_weyl.py`, `FF_rank_blocks.py` and `FF_rank_certificates.json`; output: `FF_joint_weyl.out`; record: `metrics/geometric-unity/FF-joint-weyl.yaml`. The calculation uses U's reviewed assembly and the earlier referee curvature engine. [Independent internal Review32](../../reports/reviews/REVIEW32-GU-FF-2026-11-02.md) accepts FF1–FF6 and reproduces the portable calculation. This is internal review, not external peer review.

## The precise curvature orbit

At a Riemannian four-metric point, let $p$ be a unit vector spanning the nonzero Ricci eigenline and $p^\flat$ its metric-dual covector. Assume
$$\mathrm{Ric}=\rho\,p^\flat\otimes p^\flat,\qquad \rho\ne0.$$
Choose an oriented orthonormal frame of $p^\perp$ and define
$$E_{ij}=W(p,e_i,p,e_j),\qquad
B_{ij}=\tfrac12\epsilon_{jkl}W(p,e_i,e_k,e_l),\qquad i,j,k,l=1,2,3.$$
Both matrices are real symmetric and tracefree by the Weyl identities, as derived in CC1 and EE1. Reversing the chosen orientation or the sign of $p$ changes the sign of $B$ without changing the conditions below.

**Lemma FF1 (intrinsic normal form).** The conditions
$$[E,B]=0,\qquad \ker E\cap\ker B\ne\{0\}\tag{FF1}$$
are equivalent to an adapted orthonormal frame in which
$$W=\alpha E_0+\beta H,\qquad \alpha,\beta\in\mathbb R,\tag{FF2}$$
where $E_0$ has sectional entries $+1$ at $01,23$, $-1$ at $02,13$ and zero at $03,12$, with all cross-pair entries zero; and $H$ has $H_{0123}=1$, $H_{0231}=-1$, $H_{0312}=0$, with curvature symmetries.

*Proof.* Normalize a common kernel vector as $e_3$. Symmetry makes its orthogonal two-plane invariant under both matrices. Their restrictions are commuting real symmetric tracefree two-by-two matrices, hence simultaneously diagonalizable. An eigenvector sign adjustment gives an oriented perpendicular frame without changing the diagonal entries. Thus
$$E=\operatorname{diag}(\alpha,-\alpha,0),\qquad
B=\operatorname{diag}(\beta,-\beta,0).$$
This includes either matrix being zero. The electric and magnetic matrices determine the full Weyl tensor: the spatial three-dimensional tensor has Ricci $-E$ and scalar zero, hence equals $(-E)\owedge\delta$, and the mixed tensor is recovered from $B$ by $W_{0ikl}=\epsilon_{jkl}B_{ij}$. These give exactly (FF2). Conversely, extracting the two matrices from $E_0,H$ gives the displayed diagonal profiles, proving (FF1). The code checks their precise signs and common kernel. ∎

Let $R_0=\tfrac12(e^0\otimes e^0-\delta/6)\owedge\delta$. Normal coordinates give the smooth expansion
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\rho R_0+\alpha E_0+\beta H)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\tag{FF3}$$
where the cubic and quartic symmetric matrices are arbitrary and all derivatives of the smooth remainder through order four vanish at $0$. Positivity holds on a sufficiently small neighbourhood, with no uniform radius over these coefficients.

Use the reconstructed Z1 observer metric at $\lambda=1/5$, $c_H,c_V>0$. The full field statement assumes S.0's compact $(14,0)$ real form with gauge algebra $\mathfrak u(128)$, invariant odd field, slaved vector even field, one fixed reconstructed Shiab member with coefficients $(x,y)\ne(0,0)$, $\kappa\ne0$ and $\det M\ne0$. U's naturality, lifted dilation and observer homothety normalize
$$\rho=\sigma\in\{1,-1\},\quad c_H=2,\quad c_V=3,\qquad
a=\alpha/|\rho|,\quad b=\beta/|\rho|.\tag{FF4}$$
The normalized cubic and quartic terms $Q_3/|\rho|^{3/2}$ and $Q_4/\rho^2$ remain arbitrary. No relation between $a$ and $b$ is imposed.

## Degree bounds and compatibility

First set the higher terms to zero and work on the auxiliary quadratic metric.

**Lemma FF2 (entry degrees).** At a fixed positive fibre matrix, the relaxed coefficient entries have total degree at most two in $(a,b)$ and the RHS entries total degree at most three. Their denominators are parameter independent.

*Proof.* The quadratic base coefficient is affine in $(a,b)$. BB1's fixed-curvature jet argument gives a base connection linear term of degree at most one, zero quadratic term, and cubic term of degree at most two. The observer metric and inverse at the point are fixed; its first derivatives have degree at most one and its second and third derivatives degree at most two. Curvature has degree at most two and its covariant derivative at most three. These are the tensors entering the relaxed system. The bounds precede interpolation. ∎

The ten triangular nodes $(i,j)$, $i,j\ge0$, $i+j\le3$, are unisolvent for total degree three. The first six with $i+j\le2$ are unisolvent for degree two; their absolute Vandermonde determinants in the code's order are $576$ and $4$.

**Proposition FF3 (three-fibre compatibility cover).** For each sign and every real $(a,b)$, one of
$$y_1=\operatorname{diag}(2,3,4,5),\quad
y_2=\operatorname{diag}(2,3,5,7),\quad
y_3=\operatorname{diag}(2,4,5,9)$$
has incompatible selected rows.

*Certificate.* At each of these stated fibres, rows $(168,337,492)$, or $(a,b,d)=(0,13,0),(1,13,1),(2,13,2)$, reduce exactly to
$$L_y(n)+q_i(a,b)\xi_{13}=r_i(a,b),\qquad i=0,1,2.\tag{FF5}$$
Here $L_y$ is the common metric combination, independent of curvature. Every unwanted coefficient vanishes at the six degree-two nodes, hence identically within FF2's bound. The four remaining nodes check the reconstructed coefficient rows. The ten-node RHS interpolation has zero cubic coefficients as an exact calculation, rather than an assumed smaller degree bound.

Let $P_{\sigma,j}$ be the primitive numerator of the determinant of rows $(q_i,1,r_i)$ at $y_j$. All coefficients, RHS polynomials and six quartic numerators are printed in the output. For each sign, exact Gröbner computation gives
$$\langle P_{\sigma,1},P_{\sigma,2},P_{\sigma,3}\rangle
=\langle1\rangle\quad\text{in }\mathbb Q[a,b].\tag{FF6}$$
Thus the three configurations have no common compatibility exception, even over the complex numbers. Only real curvature parameters are used in the metric theorem. ∎

## Constant rank at those same fibres

Use the 722-row cubic-invariant coefficient subsystem $M$ from BB3/DD3. Its four independent kernel columns at a diagonal fibre are
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\qquad i=0,1,2,3,\tag{FF7}$$
in the $(v,\xi,n)$ ordering. Six exact node evaluations prove $MK=0$ identically by FF2.

**Proposition FF4 (same-witness rank).** For both signs and every pair $(a,b)$, $M$ has rank 38 and kernel $\operatorname{span}K$ at each of the three numerical fibres. The full coefficient matrix has rank 42 there.

*Constant block certificates.* Even coordinate sign flips preserve both fixed tensors. Odd flips reverse the mixed tensor but preserve the electric tensor, so only the even-sign subgroup is used for the quadratic rank calculation. Its characters split the 42 unknowns into a 12-column trivial block, three complementary off-diagonal 6-column blocks and four horizontal 3-column blocks. The helper checks the polynomial support at the six coefficient nodes.

For the 12-column and 6-column blocks, constant rational matrices satisfy $TM(a,b)=I$. The horizontal blocks have constant rank-two row outputs and their respective kernel columns in (FF7). The compact JSON gives all row indices and rational coefficients; the helper checks every constant identity against the freshly computed node matrices. Consequently
$$\operatorname{rank}M=12+3\cdot6+4\cdot2=38.$$
These are degree-bounded polynomial identities, not ranks inferred from samples or inherited from the diagonal $a=b$ pencil.

Append $J_7=(91,259,413,553)$, the rows $(i,7,7)$ for $i=0,1,2,3$, with coordinate $7=y_{03}$. Direct reconstruction gives the following constant restrictions, independent of both parameters and the Ricci sign:

| fibre | $J_7K$ |
|---|---|
| $y_1$ | $\operatorname{diag}(-9/125,-27/250,-18/125,-9/50)$ |
| $y_2$ | $\operatorname{diag}(-5/49,-15/98,-25/98,-5/14)$ |
| $y_3$ | $\operatorname{diag}(-49/405,-98/405,-49/162,-49/90)$ |

Every restriction is invertible. The helper verifies these identities at all six nodes, establishing rank 42 at each compatibility configuration without parameter exceptions or a determinant search. ∎

## Higher jets and the local equation

**Proposition FF5 (arbitrary higher jets).** FF3–FF4 remain valid for all higher terms in (FF3).

*Proof.* DD3 audits BB3's support and parity arguments for any fixed quadratic curvature, without requiring its electric or mixed part to vanish. Cubics vary only the mixed $HV$ second observer jet. At a diagonal vertical slot its degree-two variation is a base gradient, eliminating that slot by commuting derivatives. The remaining curvature variation has exactly three horizontal slots and one off-diagonal vertical slot; Ricci variation has horizontal/off-diagonal-vertical support. This fixes $M$ and the selected three coefficient rows.

The simultaneous inversion $(x,y)\mapsto(-x,y)$ fixes every quadratic base metric, reverses the cubic jet and preserves the quartic jet. Affine cubic dependence and the two-horizontal parity of the selected RHS rows eliminate their cubic variation. W2's quartic gradient calculation eliminates the remaining principal third-observer-jet variation. A smooth zero-four-jet remainder affects neither relevant observer jet.

In a completion row $(i,7,7)$, the only variable cubic coefficient multiplies $\xi_7$, which vanishes on $K$. Quartics change no coefficient. Thus the constant invertible completion survives all higher jets. The separate-coordinate even-sign symmetry is used only for the auxiliary quadratic rank calculation and is not imposed on the actual higher terms. ∎

**Theorem FF6.** Suppose the smooth Riemannian base has nonzero rank-one Ricci tensor at a point and its electric and magnetic Weyl matrices there satisfy (FF1). For reconstructed Z1 at $\lambda=1/5$, $c_H,c_V>0$, the relaxed system is consistent on no nonempty fibre-open patch over that point. Under the explicit S.0 assumptions, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* FF1 puts the actual curvature in the declared two-parameter frame orbit. FF3 selects an incompatible configuration for the actual normalized parameters. FF4–FF5 give coefficient rank 42 at that same fibre for the actual higher metric jets, hence augmented rank 43. S4's rational-minor argument excludes every fibre-open consistency patch. U's reviewed parameter maps transfer both ranks to the original parameters.

For the full equation, fix one actual nonzero augmented 43-minor at its positive coordinate fibre matrix. It stays nonzero over a base neighbourhood, which can also be chosen to have nonzero curvature because Ricci is nonzero at the initial point. S5 excludes Einstein pieces, while Q5 on the dense variable patch with $H\ne0$ would supply relaxed consistency. This contradicts the fibrewise minor obstruction as in U5. The continued object is the fixed full-system minor; (FF1) is required only at the initial point. ∎

## Reproduction and boundary

Run `python -X utf8 -I computation/geometric-unity/FF_joint_weyl.py` from the repository root with the pinned dependencies. The default freshly assembles 60 exact node systems once in memory and passes each case's first six coefficient matrices to the rank helper. The compact JSON stores rational row combinations and constant completion data, with no raw node matrices or off-repository cache dependency.

The recorded fresh run passes 31 main checks and 60 rank-helper checks in 107.1 s; Review32 independently reproduces them in 109.6 s. The certificate file is 46,360 bytes.

This includes the one-parameter $a=b$ joint pencil and allows arbitrary independent amplitudes in the same two fixed directions. It does not cover commuting electric and magnetic matrices without a common kernel, noncommuting pairs, arbitrary joint Weyl curvature, or Ricci-flat nonflat points. CC and EE separately cover their respective pure-part conditions without FF's common-kernel restriction.

[GG](GG-ricci-zero-weyl.md) separately treats a nonzero-Weyl Ricci-zero boundary with singular chiral blocks, using new rank certificates rather than a limit of this theorem.
