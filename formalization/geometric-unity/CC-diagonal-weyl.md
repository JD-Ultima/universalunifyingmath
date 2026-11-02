# CC. Rank-one Ricci curvature with vanishing mixed Weyl components

The obstruction extends to every diagonal Weyl tensor in a frame adapted to the nonzero rank-one Ricci tensor. Arbitrary higher metric jets are allowed. Intrinsically, the additional curvature hypothesis is that the Weyl components with one Ricci-line argument and three perpendicular arguments vanish. The full Weyl tensor need not vanish.

The proof combines four finite compatibility configurations with exact rank certificates at those same configurations. It does not assume that a sparse three-row identity verified at one fibre remains sparse at other fibres. [BB](BB-weyl-direction.md) retains the simpler one-direction certificate and supplies the higher-jet support argument used here.

Code: `computation/geometric-unity/CC_diagonal_weyl.py`, `CC_rank_blocks.py` and `CC_rank_certificates.json`; output: `CC_diagonal_weyl.out`; metric record: `metrics/geometric-unity/CC-diagonal-weyl.yaml`. The calculation uses U's reviewed relaxed-system assembly and the earlier referee curvature engine. [Independent internal Review24](../../reports/reviews/REVIEW24-GU-CC-2026-11-02.md) accepts CC1–CC6 and reproduces the fresh computation. This is internal review, not external peer review.

## Curvature condition and normalization

Use the reconstructed Z1 observer metric at $\lambda=1/5$, with $c_H,c_V>0$. The full twisted-cone conclusion uses S.0's compact $(14,0)$ real form with gauge algebra $\mathfrak u(128)$, invariant odd field, slaved vector even field, one fixed reconstructed Shiab member with Einstein coefficients $(x,y)\ne(0,0)$, $\kappa\ne0$ and $\det M\ne0$. Let $p$ span the nonzero Ricci eigenline at the base point, and write $p^\flat$ for its metric-dual covector. The intrinsic hypotheses are
$$\mathrm{Ric}=\rho\,p^\flat\otimes p^\flat,\qquad |p|=1,\quad \rho\ne0,
\qquad W(p,u,v,w)=0\quad(u,v,w\perp p).\tag{CC1}$$

**Lemma CC1 (normal form).** Condition (CC1) is equivalent to the existence of an orthonormal frame with $p=e_0$ and Weyl sectional entries
$$W_{0101}=W_{2323}=\alpha,\qquad
W_{0202}=W_{1313}=\beta,\qquad
W_{0303}=W_{1212}=-(\alpha+\beta),\tag{CC2}$$
all cross-pair components zero and the remaining entries determined by curvature symmetries.

*Proof.* The matrix $E_{ij}=W_{0i0j}$, $i,j=1,2,3$, is symmetric. Its trace is zero by the Weyl tensor's zero Ricci contraction, independently of the nonzero Ricci tensor of the base metric. The spatial curvature tensor has Ricci tensor $-E$ and scalar curvature zero. In three dimensions the Ricci contraction determines the curvature tensor, here as $(-E)\owedge\delta$. An $SO(3)$ change of perpendicular frame diagonalizes $E$ to $(\alpha,\beta,-\alpha-\beta)$. Its complementary spatial sectional entries are respectively $\alpha,\beta,-\alpha-\beta$, giving (CC2). The mixed components are zero by (CC1). Conversely, the displayed sectional tensor has no mixed components. The choice of sign of $p$ does not affect the vanishing condition. ∎

Put $R_0=\tfrac12(e_0\otimes e_0-\delta/6)\owedge\delta$. Let $\mathcal W_A,\mathcal W_B$ be the two Weyl tensors obtained from (CC2) at $(\alpha,\beta)=(1,0),(0,1)$. Both have zero Ricci contraction and satisfy Bianchi. In normal coordinates, every metric under consideration has the expansion
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\rho R_0+\alpha \mathcal W_A+\beta \mathcal W_B)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\tag{CC3}$$
where $Q_3,Q_4$ are arbitrary symmetric homogeneous cubic and quartic matrices and the smooth remainder has zero derivatives through order four at $0$. The metric is positive definite on a sufficiently small neighbourhood, without a uniform radius over these coefficients.

U's naturality, lifted base dilation and observer homothety normalize
$$\rho=\sigma\in\{1,-1\},\quad c_H=2,\quad c_V=3,\qquad
a=\alpha/|\rho|,\quad b=\beta/|\rho|.\tag{CC4}$$
The normalized higher terms are $Q_3/|\rho|^{3/2}$ and $Q_4/\rho^2$, still arbitrary. The curvature parameters $a,b$ range independently over all real numbers.

## Exact compatibility at three fibres

First set the higher terms to zero. All reflection arguments in the rank calculation concern this auxiliary quadratic background. Arbitrary higher terms need not possess these symmetries.

**Lemma CC2 (degree bound).** At a fixed fibre matrix, every entry of the auxiliary relaxed coefficient matrix has total degree at most two in $(a,b)$, and every RHS entry has total degree at most three. Denominators are independent of these parameters.

*Proof.* The quadratic base coefficient is affine in $(a,b)$. As in BB1, the base connection has a linear term of total degree at most one, zero quadratic term and a cubic term of total degree at most two. The observer metric and inverse at the point are parameter independent; its first derivatives have degree at most one and its second and third derivatives degree at most two. Curvature has degree at most two, and its covariant derivative at most three. These are the tensors entering the relaxed system. This is a total-degree argument before interpolation. ∎

The ten nodes $(i,j)$ with $i,j\ge0$, $i+j\le3$, determine all degree-three polynomials. The first six, with $i+j\le2$, determine degree two. In the code's graded monomial order, the absolute Vandermonde determinants are $576$ and $4$. Every coefficient identity below is therefore certified by six exact evaluations, and every RHS identity by ten. The four remaining nodes additionally check the reconstructed coefficients in the three certificate rows.

Use the positive fibre matrices
$$y_1=\operatorname{diag}(2,3,4,5),\quad
y_2=\operatorname{diag}(2,3,5,7),\quad
y_3=\operatorname{diag}(2,4,5,9).$$
Rows $I=(168,337,492)$ are $(a,b,d)=(0,13,0),(1,13,1),(2,13,2)$ in the ordering of U. Vertical coordinate $13$ is $y_{33}$. At each of these three specific fibres they have the exact form
$$L_{y_j}(n)+q_i(a,b)\,\xi_{13}=r_i(a,b),\qquad i=0,1,2,\tag{CC5}$$
with one common curvature-independent metric combination $L_{y_j}(n)$. Every other coefficient is zero. These are polynomial identities in curvature at the stated fibres. The code prints all nine $q_i$ polynomials and their exact coefficients.

The RHS formulas are compact:
$$r_0=s_{j,0}(3a+3b-\sigma)^2,\qquad
r_1=s_{j,1}(6b-\sigma)^2,\qquad
r_2=s_{j,2}(6a-\sigma)^2,$$
where

| fibre | $s_{j,0}$ | $s_{j,1}$ | $s_{j,2}$ |
|---|---:|---:|---:|
| $y_1$ | $1/500$ | $1/4500$ | $1/12000$ |
| $y_2$ | $5/4116$ | $1/6174$ | $1/20580$ |
| $y_3$ | $7/8748$ | $5/69984$ | $1/21870$ |

Let $P_{\sigma,j}$ be the primitive rational numerator of the determinant of the three rows $(q_i,1,r_i)$. The full quartics are printed in the output, making the following computer-assisted polynomial certificate reproducible.

**Proposition CC3 (four finite compatibility configurations).** For either sign,
$$\langle P_{\sigma,1},P_{\sigma,2},P_{\sigma,3}\rangle
=\langle(6a-\sigma)^2,(6a-\sigma)(6b-\sigma),(6b-\sigma)^2\rangle
\quad\text{in }\mathbb Q[a,b].\tag{CC6}$$
Thus their sole common complex, and hence real, zero is $a=b=\sigma/6$. Adding the configuration obtained by swapping perpendicular axes $2,3$ removes this exception.

*Certificate and symmetry.* Exact Gröbner reduction over $\mathbb Q[a,b]$ gives (CC6). The signed swap that also flips $e_0$ is in $SO(4)$ and acts on the tensors as
$$R_0\mapsto R_0,\quad \mathcal W_A\mapsto \mathcal W_A-\mathcal W_B,\quad \mathcal W_B\mapsto-\mathcal W_B.$$
It therefore sends $(a,b)$ to $(a,-a-b)$ and the exceptional point to $(\sigma/6,-\sigma/3)$. The rotated first compatibility numerator at that point is $-204984$ for $\sigma=1$ and $209736$ for $\sigma=-1$. Equivalently the exact polynomial calculation verifies
$$\langle P_{\sigma,1},P_{\sigma,2},P_{\sigma,3},
P_{\sigma,1}(a,-a-b)\rangle=\langle1\rangle.\tag{CC7}$$
The fourth configuration in the original frame is $y_4=\operatorname{diag}(2,3,5,4)$ with rows $(140,309,605)$ and diagonal vertical slot $y_{22}$. Naturality transfers the first configuration's equations and ranks to it. No general-diagonal-fibre sparsity identity is required. ∎

## Rank at those same fibres

Pointwise incompatibility alone is insufficient for the fibre-open theorem if coefficient rank is unproved. The following certificates establish rank uniformly at each of the three numerical fibres, separately from (CC6).

Let $M$ be the 722 coefficient rows fixed by cubic variations in BB3: $VVV$, $H\!D\!H$, $H\!V\!D$, $HH\!D$ and $DD\!H$, of counts $450,64,160,24,24$. Write $D=\{4,8,11,13\}$. For $y=\operatorname{diag}(d_0,d_1,d_2,d_3)$ its four independent kernel columns are
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\qquad i=0,1,2,3.\tag{CC8}$$
The code verifies $MK=0$ at the six unisolvent nodes, hence identically by CC2.

**Proposition CC4 (uniform rank).** At each $y_1,y_2,y_3$, for both signs and every real $(a,b)$, $M$ has rank 38 and kernel exactly $\operatorname{span}K$. The full coefficient matrix has rank 42.

*Character decomposition.* On the quadratic background, the four independent coordinate sign flips preserve the curvature tensor and each diagonal fibre. A horizontal coordinate has character $s_i$, a diagonal vertical coordinate character $1$, and an off-diagonal vertical coordinate $y_{ij}$ character $s_is_j$. The full sign group, rather than only its orientation-preserving subgroup, is used here. It splits the 42 unknowns into a 12-column trivial block, six 3-column off-diagonal blocks and four 3-column horizontal blocks. The exact six-node support check proves this decomposition for the polynomial matrices. It is not asserted as a symmetry of arbitrary higher jets.

*Trivial block.* A constant rational matrix $T$ satisfies $TM_0(a,b)=I_{12}$. There are 64 nonzero rows in this block. The JSON certificate records all 218 nonzero entries of $T$ for each fibre and sign, with exact row and column indices. The helper verifies the polynomial identity against the freshly reconstructed entries. One construction takes constant combinations that cancel all nonconstant coefficients, then inverts twelve independent constant output rows. Consequently this block has rank 12 everywhere.

*Horizontal blocks.* In each such block, two columns admit a constant rational left inverse. The corresponding column of (CC8) is already a kernel vector of the full three-column block. Thus every horizontal block has rank exactly two.

*Off-diagonal blocks.* Two constant row combinations give a rank-two constant output matrix with fixed one-dimensional kernel $k_O$. Two further rows, restricted to this kernel, are nonzero rational multiples of
$$u^2=(3x+\sigma)^2,\qquad v^2=(6x-\sigma)^2,$$
with $x=a$ or $b$ as recorded for that block. Their row indices, prefactors and exact combination coefficients are in the certificate. They have no common zero; a polynomial identity independent of any real-root test is
$$\frac{(20\sigma-48x)u^2+(12x+7\sigma)v^2}{27\sigma^3}=1.\tag{CC9}$$
The helper checks the scaled identity using the actual residual polynomials. Each off-diagonal block therefore has rank three everywhere. Summing the disjoint blocks gives $12+6\cdot3+4\cdot2=38$; (CC8) is the whole kernel.

*Four-row completion.* Append rows $(i,O,O)$, $i=0,1,2,3$, with $O=5,6,7$ respectively representing $y_{01},y_{02},y_{03}$. Their restrictions to $K$ are nonzero diagonal rational matrices times
$$-(3a+\sigma)^2,\qquad -(3b+\sigma)^2,\qquad
-(\sigma-3a-3b)^2,$$
respectively; all prefactors are recorded and checked at each fibre. The three linear factors sum to $3\sigma$, so cannot vanish simultaneously. At least one completion has rank four on the kernel. This proves rank 42 at each of the three compatibility fibres, for every parameter. ∎

## Higher jets and the field equation

**Proposition CC5 (higher-jet transfer).** The compatibility and full-rank conclusions hold for all higher terms in (CC3).

*Proof.* BB3 proves the required support, reflection-parity and diagonal-vertical gradient identities for any fixed quadratic curvature, including the present two-parameter plane. Its statement covers all 200 cubic and 350 quartic coefficients, with a smooth zero-four-jet remainder. The three coefficient rows and their RHS in (CC5), and the 722 coefficient rows of $M$, are unchanged. This transfers the four-configuration incompatibility and the rank-38 kernel to every higher jet.

In a completion row $(i,O,O)$, the only cubic-dependent coefficient is that of $\xi_O$. Every column in $K$ has $\xi_V=0$, so this variation vanishes on $K$. Quartics change no coefficient. Thus the completion restriction remains invertible for at least one $O$ at the same fibre. The same naturality argument transfers this conclusion to the fourth configuration. The auxiliary reflection symmetry is used to calculate the unchanged subsystem, not imposed on the arbitrary higher terms. ∎

**Theorem CC6 (fibre-open and local obstruction).** For every smooth Riemannian four-metric satisfying (CC1) at a point, the reconstructed Z1 relaxed system at $\lambda=1/5$, $c_H,c_V>0$, is consistent on no nonempty fibre-open patch over that point. In the slaved-vector sector of S.0, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* CC3 supplies one of four incompatible finite configurations for each actual pair $(a,b)$. CC4–CC5 give coefficient rank 42 at that same configuration for the actual higher metric jets. The augmented rank is therefore 43. S4's rational-minor argument excludes every fibre-open consistency patch. U's reviewed parameter maps transfer the ranks to every original $\rho,\alpha,\beta,c_H,c_V$.

For the field equation, choose one actual nonzero augmented 43-minor at that positive coordinate fibre matrix. It stays nonzero over a base neighbourhood. Nonzero Ricci at the initial point lets this neighbourhood also have nonzero base curvature. U5's argument applies: S5 excludes Einstein pieces, and the Q5 reduction on its dense variable patch with $H\ne0$ would supply relaxed consistency, contradicting the fibrewise minor obstruction. It is the fixed full-system minor that is continued, not a sparse three-row relation at arbitrary fibres. Normal coordinates and CC1 then give the stated intrinsic scope. ∎

## Reproduction and boundary

From the repository root, run `python -X utf8 -I computation/geometric-unity/CC_diagonal_weyl.py` with the pinned dependencies. The default path freshly computes the 60 exact node systems once in memory. Their first six coefficient matrices per case are passed to the rank helper, which checks the compact JSON row-combination certificates. The driver checks the Weyl identities, unisolvence, sparse rows, full RHS interpolation, both exact compatibility ideals and the signed frame action. There are no off-repository caches, parameter searches or sampled determinant-degree assumptions.

The recorded fresh run passes 32 main checks and 96 rank-helper checks in 98.3 s; Review24 independently reproduces them in 98.2 s. These are 60 curvature-engine assemblies followed by exact rational/polynomial certificate verification; no raw node matrices are stored in the repository.

[EE](EE-zero-electric-weyl.md) separately excludes every zero-electric Weyl tensor relative to the nonzero Ricci line, including DD's one-direction subcase. [FF](FF-joint-weyl.md) additionally treats commuting electric/magnetic matrices with a common nonzero kernel. Other cases with both parts nonzero, and Ricci-flat nonflat points, remain open. The field-equation conclusion concerns the reconstructed geometric sector under the stated S.0 assumptions.
