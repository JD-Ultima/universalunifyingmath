# GG. Ricci-zero points with singular chiral Weyl blocks

[HH](HH-ricci-zero.md) extends the obstruction to all nonzero-curvature Ricci-zero points. This note preserves the simpler singular-block proof used there.

The obstruction reaches a nonflat part of the Ricci-zero boundary. At a Riemannian four-metric point with $\mathrm{Ric}=0$, suppose $W\ne0$ and both chiral Weyl blocks are singular. An exact positive compatibility polynomial and new normalized rank certificates exclude fibre-open relaxed solutions. The complete flat quadratic control retains its lower rank; zero-curvature points remain outside the theorem. Arbitrary higher metric jets are allowed at the covered nonzero-curvature points.

Unlike FF, there is no distinguished Ricci eigenline. The intrinsic normal form uses independent rotations of the two chiral three-spaces. The coefficient certificates below are calculated anew at Ricci zero; no nonzero-Ricci rank result is extended by taking a limit.

Code: `computation/geometric-unity/GG_ricci_zero_weyl.py`, `GG_rank_blocks.py`, `GG_rank_certificates.json` and `GG_chiral_action.py`; output: `GG_ricci_zero_weyl.out`; record: `metrics/geometric-unity/GG-ricci-zero-weyl.yaml`. The observer calculation uses U's reviewed assembly and the earlier referee curvature engine. The small chiral-action checker has no observer-engine dependency. [Independent internal Review35](../../reports/reviews/REVIEW35-GU-GG-2026-11-02.md) accepts GG1–GG6 and reproduces the portable computation. This is internal review, not external peer review.

## Intrinsic normal form

Fix an orientation of the positive four-dimensional tangent space, with
$$\theta_i=e^0\wedge e^i,\qquad
(\tau_1,\tau_2,\tau_3)=(e^{23},e^{31},e^{12}),\qquad
\eta_i^\pm=(\theta_i\pm\tau_i)/\sqrt2.$$
These are orthonormal Hodge chiral bases. Use the curvature-operator convention $\langle\mathcal R(e^{ij}),e^{kl}\rangle=R_{ijkl}$.

**Lemma GG1 (Ricci-zero chiral orbit).** In a positive real four-space,
$$\mathrm{Ric}=0,\quad\det W^+=\det W^-=0$$
is equivalent to the existence of an oriented orthonormal frame with
$$R=\alpha E_0+\beta H,\qquad \alpha,\beta\in\mathbb R,\tag{GG1}$$
where $E_0,H$ are FF's fixed electric and mixed tensors. Nonzero curvature is equivalent to $(\alpha,\beta)\ne(0,0)$.

*Proof.* In the paired basis $(\theta,\tau)$, a Ricci-zero curvature operator has the form
$$\mathcal R=\begin{pmatrix}E&B\\B&E\end{pmatrix},\qquad
E=E^\top,\ B=B^\top,\quad\operatorname{tr}E=\operatorname{tr}B=0.$$
Indeed $\mathrm{Ric}_{00}=\operatorname{tr}E$, the spatial Ricci matrix of general paired blocks $(E,B;B^\top,D)$ is $E+(\operatorname{tr}D)I-D$, the mixed Ricci entries give the antisymmetric part of $B$, and Bianchi gives $\operatorname{tr}B=0$. Thus $D=E$. The chiral blocks are $W^\pm=E\pm B$.

The action $SO(4)\to SO(3)_+\times SO(3)_-$ is onto. Explicitly, let $J_{ab}e_a=e_b$, $J_{ab}e_b=-e_a$. For cyclic $(i,j,k)$, the generators $(J_{0i}+J_{jk})/2$ and $(-J_{0i}+J_{jk})/2$ act on the chiral bases as $(L_i,0)$ and $(0,L_i)$, where $L_i e_j=\epsilon_{ijk}e_k$. The differential is onto; the portable checker verifies all six actions and the original paired-generator determinant $-8$. The group image is an open subgroup of the connected target, hence the whole target.

Each real symmetric tracefree singular three-by-three block has eigenvalues $(\lambda,-\lambda,0)$. Its eigenbasis can be made oriented by changing an eigenvector sign. Independent rotations of the two chiral spaces therefore make the blocks $\lambda_+D,\lambda_-D$, with $D=\operatorname{diag}(1,-1,0)$. But
$$E_0:(D,D),\qquad H:(D,-D),$$
so $\alpha=(\lambda_++\lambda_-)/2$, $\beta=(\lambda_+-\lambda_-)/2$ gives (GG1). The converse follows from these blocks. Zero blocks are permitted. Changing the four-space orientation interchanges the two blocks, preserving the condition that both determinants vanish. This is a pointwise statement; no Ricci-flat neighbourhood or preferred line is assumed. ∎

Normal coordinates give the actual smooth metric
$$g_{ij}(x)=\delta_{ij}-\tfrac13(\alpha E_0+\beta H)_{iajb}x^ax^b
+Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\tag{GG2}$$
where the symmetric cubic and quartic terms are arbitrary and all derivatives of the smooth remainder through order four vanish at $0$. The metric is positive definite near $0$, without a uniform radius over these coefficients.

Use Z1 at $\lambda=1/5$, $c_H,c_V>0$. The full equation has S.0's compact $(14,0)$ real form, gauge algebra $\mathfrak u(128)$, invariant odd field, slaved vector even field, one fixed reconstructed Shiab member with coefficients $(x,y)\ne(0,0)$, $\kappa\ne0$ and $\det M\ne0$.

## Covariance without a Ricci hypothesis

**Lemma GG2 (normalization and the complete witness).** Positive base dilation and observer rescaling normalize any nonzero pair in (GG1), while preserving coefficient rank and augmented inconsistency together at one positive fibre point. This does not require nonzero Ricci.

*Proof.* For $F_s(x)=sx$, let $g'=s^{-2}F_s^*g$. Its quadratic curvature is $R'=s^2R$, its cubic and quartic coefficients are $s^3Q_3,s^4Q_4$, and its first metric derivative at the origin is zero. Constant base homothety does not change the Levi-Civita connection. The natural lift
$$\Phi_s(x,y)=(sx,s^{-2}y)$$
therefore has the same Z1 pullback identity as U, independently of Ricci. Fibre dilation $\phi_q(x,y)=(x,qy)$ changes $c_H$ to $qc_H$ and leaves $c_V$ fixed. Set
$$\gamma=3/c_V>0,\qquad q=2c_V/(3c_H)>0.$$
Then
$$ (\Phi_s\circ\phi_q)^*(\gamma h_g(c_H,c_V))=h_{g'}(2,3).\tag{GG3}$$
Under $h\mapsto\gamma h$, the relaxed unknowns transform as $(v,\xi,n)\mapsto(v/\gamma,\xi,n/\gamma^2)$, as follows directly from the connection being fixed and the raised Ricci tensor scaling by $\gamma^{-1}$. This invertible transformation, together with pullback, preserves both coefficient and augmented ranks.

For $\alpha\ne0$, choose $s=|\alpha|^{-1/2}$. The normalized pair is $(\operatorname{sign}\alpha,t)$ with $t=\beta/|\alpha|$. For $\alpha=0$, $\beta\ne0$, choose $s=|\beta|^{-1/2}$ and normalize the pure-mixed pair to $(0,\operatorname{sign}\beta)$. All transformed higher jets remain arbitrary.

A normalized witness at $y_*$ is transported to the original positive fibre matrix $q s^{-2}y_*$. The entire coefficient/rank and incompatibility witness is transferred together. No assertion that the original fixed numeric matrix $y_*$ works at every amplitude is needed. At $(\alpha,\beta)=(0,0)$ this normalization is unavailable. ∎

## Compatibility and the flat origin

**Proposition GG3 (positive compatibility polynomial).** For the auxiliary quadratic metric in (GG2), at
$$y_* =\operatorname{diag}(2,3,5,7),$$
the three selected rows are incompatible for every real nonzero $(\alpha,\beta)$.

*Degree proof and certificate.* The coefficient total degree is at most two and the RHS total degree at most three before interpolation. BB1's proof uses only an affine normal quadratic curvature and a parameter-independent observer metric/inverse at the point, so applies unchanged when Ricci is zero. Curvature has degree at most two and its covariant derivative at most three. The normalized $t$ families consequently retain coefficient degree at most two.

Six triangular nodes determine all coefficient entries and ten determine the RHS within these bounds. Rows $(168,337,492)$ have only
$$L_{y_*}(n)+q_i(\alpha,\beta)\xi_{13}=r_i(\alpha,\beta),\qquad i=0,1,2,$$
where $L_{y_*}$ is the common metric combination independent of curvature. The driver prints the exact coefficient and RHS formulas and checks them at five additional normalized nodes. Their compatibility determinant is
$$\det[(q_i,1,r_i)]
=\frac{529\alpha^4+1052\alpha^2\beta^2+735\beta^4}{14406000}.\tag{GG4}$$
The numerator is a positive weighted sum of $(\alpha^2)^2$, $(\alpha\beta)^2$ and $(\beta^2)^2$, hence is strictly positive for every real nonzero pair. This is an exact polynomial/SOS certificate, not a sampled sign test.

At the origin it vanishes. Fresh evaluation of the complete flat quadratic metric gives coefficient/augmented ranks $28/28$ at $y_*$. No origin exclusion, continuity across its rank defect, or unit-ideal claim including that point is made. The pure-mixed line is retained and checked separately below. ∎

## New normalized coefficient ranks

**Proposition GG4 (rank at the normalized compatibility witness).** At $y_*$, the cases $\alpha=\pm1$, $\beta=t$ have coefficient rank 42 for every real $t$. The two pure-mixed normalized cases $(0,\pm1)$ have coefficient/augmented ranks $42/43$.

*Polynomial rank certificates.* These matrices are computed anew with $\mathrm{Ric}=0$. The 722-row cubic-invariant subsystem has the four independent kernel columns
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\qquad i=0,1,2,3.$$
For each sign of $\alpha$, exact checks at $t=-1,0,1$ prove the degree-two identities $MK=0$ and the even-sign character support. The unknowns split into a 12-column block, three complementary off-diagonal 6-column blocks and four horizontal 3-column blocks.

New constant rational left inverses give ranks 12 and 6 on the first blocks. The horizontal blocks have constant rank-two row outputs and their respective kernel columns. The compact JSON records the coefficients and row indices. Identities checked at the three unisolvent nodes give invariant rank $12+3\cdot6+4\cdot2=38$ for every $t$, rather than a rank inference from those values.

The old $J_7K$ completion is zero at Ricci zero. Instead $J_5=(61,229,383,523)$, or $(i,5,5)$, has
$$J_5K=-(35+6t^2)\operatorname{diag}(1/175,3/350,1/70,1/50).\tag{GG5}$$
Its real scalar factor is strictly positive, so it completes the rank to 42 for every real $t$. Together with (GG4), this gives augmented rank 43 at that same normalized fibre.

For $(0,\pm1)$, separate exact rational elimination gives invariant/full/augmented ranks $38/42/43$, with invertible
$$J_5K=\operatorname{diag}(-6/175,-9/175,-3/35,-3/25).$$
These are fixed-case certificates; their nonzero pure-line amplitudes are covered by GG2's positive covariance, not by interpolation from cases. ∎

## Higher jets and the local equation

**Proposition GG5 (higher-jet transfer).** GG3–GG4 apply to every higher jet in (GG2) after the normalization of GG2.

*Proof.* DD3's support proof requires only a fixed normal quadratic curvature, not nonzero Ricci. Cubics change only the mixed $HV$ second observer jet. The diagonal-vertical gradient cancellation removes that vertical slot, leaving only three-horizontal/one-off-diagonal-vertical curvature variation and horizontal/off-diagonal-vertical Ricci variation. Thus the 722-row subsystem and the selected coefficient rows are unchanged.

Simultaneous base inversion fixes every quadratic metric, reverses the cubic jet and gives even parity to the selected two-horizontal RHS rows. Their affine cubic variation is therefore zero. The quartic gradient cancellation and the zero-four-jet remainder argument are likewise independent of Ricci. In a completion row $(i,5,5)$ only the cubic coefficient multiplying $\xi_5$ varies; it annihilates $K$. Quartics change no coefficient. Thus the rank-42 and incompatible normalized witness survives all higher jets, without imposing baseline even-sign symmetry on them. GG2 transfers the entire actual witness back together. ∎

**Theorem GG6.** At a Riemannian four-metric point satisfying
$$\mathrm{Ric}=0,\qquad W\ne0,\qquad \det W^+=\det W^-=0,$$
the reconstructed Z1 relaxed system at $\lambda=1/5$, $c_H,c_V>0$, is consistent on no nonempty fibre-open patch. Under the explicit S.0 assumptions, the twisted cone has no $\Upsilon=0$ solution on an open set meeting that fibre.

*Proof.* GG1 supplies the fixed-plane frame. GG2 normalizes the actual nonzero curvature pair. GG3–GG5 supply a rank-42 incompatible witness at one and the same normalized positive fibre, and GG2 transports it to one original positive fibre. Its augmented rank is 43, so S4 excludes every fibre-open consistency patch.

For the full equation, fix one actual nonzero augmented 43-minor and its coordinate fibre matrix. It remains nonzero over a base neighbourhood. Here $R=W\ne0$ at the initial point, so shrink that neighbourhood to keep Riemann curvature nonzero. Nonzero Ricci is not needed. A constant-field open piece would make Z1 Einstein; S5(b) for the definite four-base would force flatness on its open projection, a contradiction. Q5's metric-independent dense variable patch with $H\ne0$ then gives the forbidden relaxed-system patch, exactly as in U5. The pointwise Ricci-zero and chiral hypotheses need not persist on the neighbourhood. The fixed full-system minor is the object continued there. ∎

## Reproduction and boundary

Run `python -X utf8 -I computation/geometric-unity/GG_ricci_zero_weyl.py` from the repository root with the pinned dependencies. The default assembles the required union of 15 exact systems at $y_*$ once in memory: ten compatibility nodes, normalized rank nodes for both signs, and the second pure-mixed sign. It then verifies compact new rank data and the separate portable chiral-action controls, with no scratch paths, stored raw matrices or cache dependency.

The recorded fresh run passes 10 main checks, 28 rank checks and 6 chiral checks in 20.9 s; Review35 independently reproduces them in 20.9 s. Its separate nonunit control at $\alpha=-2$, $\beta=3$, $c_H=5$, $c_V=7$ verifies ranks $42/43$ at the transported fibre matrix $(28/15)y_*$. The rank certificate is 10,906 bytes.

Ricci-zero Weyl tensors with a nonsingular chiral block are outside this note’s theorem; [HH](HH-ricci-zero.md) covers them. Zero-curvature points and indefinite real forms remain outside the present results. This is a pointwise curvature condition in the reconstructed sector, not a classification of all Ricci-flat metrics or a separate weighted base-metric equation.
