# DD. One mixed Weyl direction

This result treats one Weyl direction with zero electric part and nonzero mixed components relative to the Ricci eigenline. Its nonzero members lie outside CC's vanishing-mixed-Weyl condition. Two finite fibre configurations give incompatible relaxed equations at every amplitude, and constant rational rank certificates hold at those same fibres. Arbitrary higher metric jets are included through a separate support argument.

Code: `computation/geometric-unity/DD_mixed_weyl.py`, `DD_rank_blocks.py` and `DD_rank_certificates.json`; output: `DD_mixed_weyl.out`; record: `metrics/geometric-unity/DD-mixed-weyl.yaml`. The computation uses U's reviewed assembly and the earlier referee curvature engine. [Independent internal Review26](../../reports/reviews/REVIEW26-GU-DD-2026-11-02.md) accepts DD1–DD4 and reproduces the portable computation. This is internal review, not external peer review.

## Declared curvature and scope

Let $\theta^{ij}=e^i\wedge e^j$ in an orthonormal base frame and set
$$H=\theta^{01}\otimes\theta^{23}+\theta^{23}\otimes\theta^{01}
+\theta^{02}\otimes\theta^{13}+\theta^{13}\otimes\theta^{02}.\tag{DD1}$$
Thus $H_{0123}=1$, $H_{0231}=-1$ and $H_{0312}=0$. Pair antisymmetry and interchange hold by construction. Bianchi reduces to $1-1+0=0$. Every nonzero component has four distinct indices, so the Ricci contraction, all sectional entries and the electric entries $H_{0i0j}$ vanish. Hence $H$ is an algebraic Weyl tensor. Its mixed part is nonzero.

Write $R_0=\tfrac12(e^0\otimes e^0-\delta/6)\owedge\delta$. The declared smooth local metrics have expansion
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\rho R_0+\omega H)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\qquad \rho\ne0,\ \omega\in\mathbb R,\tag{DD2}$$
where $Q_3,Q_4$ are arbitrary symmetric homogeneous cubic and quartic matrices and the smooth symmetric remainder has zero derivatives through order four at $0$. They are positive definite on a sufficiently small neighbourhood, without a uniform radius over their coefficients. At the origin their Ricci tensor is $\rho e^0\otimes e^0$ and Weyl tensor is $\omega H$.

Use Z1 at $\lambda=1/5$, $c_H,c_V>0$. The full field-equation consequence assumes S.0's compact $(14,0)$ real form with gauge algebra $\mathfrak u(128)$, invariant odd field, slaved vector even field, one fixed reconstructed Shiab member with coefficients $(x,y)\ne(0,0)$, $\kappa\ne0$ and $\det M\ne0$. U's naturality, lifted dilation and homothety reduce the parameters to
$$\rho=\sigma\in\{1,-1\},\quad c_H=2,\quad c_V=3,\qquad
w=\omega/|\rho|.\tag{DD3}$$
The normalized higher terms $Q_3/|\rho|^{3/2}$ and $Q_4/\rho^2$ remain arbitrary.

The intrinsic scope is exactly the simultaneous frame orbit in which $\mathrm{Ric}=\rho e^0\otimes e^0$ and $W=\omega H$. No claim about all Weyl tensors with zero electric part is made. For $\omega\ne0$, the mixed Weyl tensor relative to this Ricci line is nonzero, an invariant failure of CC's hypothesis under perpendicular frame changes.

## Degree bounds and finite compatibility

First set $Q_3=Q_4=R_5=0$ to obtain the auxiliary quadratic metric. BB1's hand jet argument applies to the affine curvature $\sigma R_0+wH$: the observer metric and inverse at the origin are independent of $w$; its first derivatives have degree at most one and its second and third derivatives degree at most two. Consequently every coefficient entry has degree at most two and every RHS entry degree at most three, with $w$-independent denominators. These bounds precede interpolation.

**Proposition DD1 (two finite compatibility witnesses).** For every real $w$ and either sign, one of
$$y_1=\operatorname{diag}(2,3,4,5),\qquad
y_2=\operatorname{diag}(2,3,5,7)$$
has three incompatible rows $(168,337,492)$, or $(a,b,d)=(0,13,0),(1,13,1),(2,13,2)$.

*Certificate.* At each of these two stated fibres the rows have only the combinations
$$L_y(n)+q_i(w)\xi_{13}=r_i(w),\qquad i=0,1,2,\tag{DD4}$$
where $L_y$ is the common metric combination independent of curvature. All unwanted coefficients vanish at $w=-1,0,1$, so vanish identically within the degree bound. Four nodes $-1,0,1,2$ determine the RHS; the fourth also checks the entire reconstructed coefficient matrix.

The determinant of the three rows $(q_i,1,r_i)$ is $P_{\sigma,j}/d_j$, where $d_1=138240000$, $d_2=777924000$, and

| sign | $P_{\sigma,1}$ | $P_{\sigma,2}$ |
|---|---|---|
| $+$ | $419830w^2-3933$ | $39690w^4+2320806w^2-16867$ |
| $-$ | $-444170w^2+5667$ | $39690w^4-2441994w^2+27233$ |

The code prints all coefficient and RHS polynomials and exact Bézout coefficients $A_\sigma,B_\sigma\in\mathbb Q[w]$ satisfying
$$A_\sigma P_{\sigma,1}+B_\sigma P_{\sigma,2}=1.\tag{DD5}$$
Thus no amplitude is exceptional for both configurations. Roots of either individual compatibility polynomial are exceptions only of that certificate. ∎

## Constant coefficient-rank certificates

The three-row contradiction needs full coefficient rank at the same witness before a fibre-open argument is available. The rank calculation below is independent of the compatibility determinants.

Take BB3's 722-row cubic-invariant subsystem $M$, with classes $VVV$, $H\!D\!H$, $H\!V\!D$, $HH\!D$ and $DD\!H$. Its four independent sparse kernel columns at a diagonal fibre are
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\qquad i=0,1,2,3,\tag{DD6}$$
in the usual $(v,\xi,n)$ ordering. Exact evaluation at the three unisolvent nodes proves $MK=0$ identically in $w$.

**Proposition DD2 (rank at both compatibility fibres).** For both signs, every $w$, and either $y_1$ or $y_2$, $M$ has rank 38 and kernel exactly $\operatorname{span}K$. The full coefficient matrix has rank 42.

*Proof and certificates.* A separate coordinate sign flip reverses the mixed tensor. Only the even-sign subgroup $s_0s_1s_2s_3=1$ fixes this quadratic background at nonzero amplitude. Its characters merge each complementary pair of off-diagonal fibre coordinates. The unknowns split as
$$12\text{ diagonal-vertical columns}
\quad+\quad3\cdot6\text{ complementary off-diagonal columns}
\quad+\quad4\cdot3\text{ horizontal columns}.\tag{DD7}$$
The helper checks this support at all three coefficient nodes; the degree bound makes it a polynomial identity. Full independent sign-flip symmetry is not assumed.

For the 12-column block and each 6-column block, a constant rational matrix $T$ satisfies $TM(w)=I$ on that block. The compact certificate records the row combinations and indices, and the helper checks them against the freshly computed matrices at the three nodes. They are identities in $w$, giving ranks 12 and 6 everywhere. Each horizontal block has a constant rank-two row output and exactly its kernel column in (DD6). Therefore
$$\operatorname{rank}M=12+3\cdot6+4\cdot2=38.$$

Append the four rows
$$J_7=(91,259,413,553),\qquad(a,b,d)=(i,7,7),\quad i=0,1,2,3,$$
where coordinate $7$ is $y_{03}$. Their restriction to $K$ is constant in $w$ and independent of the Ricci sign:

| fibre | $J_7K$ |
|---|---|
| $y_1$ | $\operatorname{diag}(-9/125,-27/250,-18/125,-9/50)$ |
| $y_2$ | $\operatorname{diag}(-5/49,-15/98,-25/98,-5/14)$ |

Both are invertible. The degree-two identities are checked at all three nodes. Thus appending $J_7$ gives rank 42 at each compatibility witness for every amplitude. No polynomial determinant search or rank extrapolation from evaluated matrices is used. ∎

## Why arbitrary higher jets are allowed

**Proposition DD3 (higher-jet transfer).** DD1–DD2 remain valid for all higher terms in (DD2).

*Proof.* This step uses BB3's support, parity and gradient proof, whose independence of the fixed quadratic curvature must be checked here. With any fixed quadratic curvature, cubics change only the mixed $HV$ second metric derivatives. For a diagonal vertical slot $A=y_{aa}$, their leading degree-two variation is
$$\delta h_{jA}=\partial_jF_a,\qquad
F_a=-\frac{c_V}{d_a}(Q_{3;aa}-\lambda\operatorname{tr}Q_3).$$
Commuting third derivatives eliminates curvature variation with three horizontal slots and a diagonal vertical slot. The remaining cubic curvature variation has exactly three horizontal slots and one off-diagonal vertical slot; Ricci variation has only horizontal/off-diagonal-vertical slots. Lower observer jets and the inverse metric are fixed. None of these statements uses zero mixed curvature. This fixes $M$ and the three coefficient rows of (DD4).

The simultaneous base reflection $(x,y)\mapsto(-x,y)$ fixes every quadratic base metric, including the present mixed one, and reverses $Q_3$. The selected RHS rows have two horizontal indices and even parity. Their dependence on $Q_3$ is affine, so their cubic variation is zero. For quartics, the analogous leading degree-three mixed variation is again a gradient. W2's commuting-fourth-derivative calculation cancels its principal covariant-Ricci variation in these rows. The zero-four-jet remainder affects neither observer 2-jet nor 3-jet. Thus the RHS is unchanged as well.

In a completion row $(i,7,7)$, the cubic support excludes Riemann variation. The only variable Ricci coefficient multiplies $\xi_7$, which vanishes on $K$. Quartics change no coefficient. Consequently the constant completion restriction remains invertible. The even-sign symmetry used to calculate the quadratic rank blocks is not imposed on arbitrary higher jets. ∎

## Local consequence and remaining scope

**Theorem DD4.** Let a smooth Riemannian four-metric have an orthonormal frame at a point in which $\mathrm{Ric}=\rho e^0\otimes e^0$, $\rho\ne0$, and $W=\omega H$. For the reconstructed Z1 observer metric at $\lambda=1/5$, $c_H,c_V>0$, its relaxed system is consistent on no nonempty fibre-open patch over that point. Under the stated S.0 assumptions, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* DD1 gives an incompatible configuration for the actual amplitude. DD2–DD3 give coefficient rank 42 at that same configuration for the actual higher metric jets. Hence the augmented rank is 43. S4's rational-minor argument excludes every fibre-open consistency patch. U's reviewed parameter maps cover all original nonzero $\rho$, all $\omega$ and positive observer scales. Normal coordinates put every smooth metric with the declared pointwise frame condition into (DD2).

For the field equation, fix one actual nonzero augmented 43-minor at that positive coordinate fibre matrix. It stays nonzero over a base neighbourhood. Shrink the neighbourhood to preserve nonzero curvature, using nonzero Ricci at the initial point. S5 excludes Einstein pieces, and Q5 on its dense variable patch with $H\ne0$ would give relaxed consistency. This contradicts the fibrewise obstruction, as in U5. The continued object is the fixed full-system minor, not a sparse three-row relation at arbitrary fibres. ∎

The portable default command is `python -X utf8 -I computation/geometric-unity/DD_mixed_weyl.py`. It freshly assembles 16 exact node systems once in memory and checks the compact rational certificates against those matrices. There are no off-repository caches or raw stored node matrices. The result treats one declared mixed-Weyl pencil and its simultaneous Ricci-frame orbit; general mixed Weyl tensors and Ricci-flat nonflat points remain outside its scope.

The recorded fresh run passes 26 main checks and 40 rank-helper checks in 27.6 s; Review26 independently reproduces them in 26.3 s. Its separate mixed cubic/quartic control changes hundreds of other entries while preserving the protected subsystem, selected rows and completion. The JSON certificate is 22,252 bytes and contains rational row combinations and the constant completion data.
