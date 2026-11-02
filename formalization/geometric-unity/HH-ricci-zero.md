# HH. Nonflat Ricci-zero points

At every point of a positive Riemannian four-metric with zero Ricci tensor and nonzero Riemann curvature, the reconstructed Z1 relaxed system has no nonempty fibre-open consistency patch. The proof combines GG's singular chiral-block case with a finite spectral-frame identity and new normalized coefficient-rank certificates. All higher metric jets are allowed. The pointwise Ricci-zero condition need not persist on a neighbourhood.

The geometric relaxed-system statement, the twisted-cone consequence and the all-member corollary are distinguished below. The latter uses S16(b) as well as the new twisted-cone exclusion. These are results in the reconstructed invariant/slaved sector, not statements about arbitrary even fields or non-invariant odd fields.

Code: `computation/geometric-unity/HH_ricci_zero.py`, `HH_rank_blocks.py`, `HH_rank_certificates.json`, `HH_orbit_identity.py` and `HH_orbit_certificate.json`; output: `HH_ricci_zero.out`; record: `metrics/geometric-unity/HH-ricci-zero.yaml`. [Independent internal Review36](../../reports/reviews/REVIEW36-GU-HH-2026-11-02.md) accepts HH1–HH6 and freshly reproduces the portable computation. This is internal AI-assisted mathematical review. External peer review and Lean formalization are not claimed.

## Chiral frames and the division into cases

Fix an orientation at the base point and use GG1's orthonormal bases
$$\theta_i=e^0\wedge e^i,\qquad
(\tau_1,\tau_2,\tau_3)=(e^{23},e^{31},e^{12}),\qquad
\eta_i^\pm=(\theta_i\pm\tau_i)/\sqrt2.$$
The curvature-operator convention is $\langle\mathcal R(e^{ij}),e^{kl}\rangle=R_{ijkl}$. For $\mathrm{Ric}=0$, the chiral blocks $W^\pm$ are independent real symmetric tracefree three-by-three matrices. GG1 proves that $SO(4)\to SO(3)_+\times SO(3)_-$ is onto; its six explicit infinitesimal actions and the connected-group argument are unchanged here.

**Lemma HH1 (frame choice).** Every Ricci-zero curvature tensor has an oriented orthonormal frame in which
$$W^+=\operatorname{diag}(p,q,-p-q),\qquad
W^-=\operatorname{diag}(r,s,-r-s).\tag{HH1}$$
In the paired $(\theta,\tau)$ basis its operator is
$$\mathcal R=\begin{pmatrix}E&B\\B&E\end{pmatrix},\quad
E=\tfrac12\operatorname{diag}(p+r,q+s,-p-q-r-s),\quad
B=\tfrac12\operatorname{diag}(p-r,q-s,-p-q+r+s).\tag{HH2}$$
All 36 independent permutations of the two eigenvalue triples are realized by oriented four-space frame changes. If at least one block is nonsingular, the orientation can first be chosen so that $W^+$ is nonsingular; then the first eigenvalue $p$ is nonzero in every one of these 36 frames.

*Proof.* Apply the spectral theorem independently to the two real symmetric blocks. Every orthogonal spectral permutation can be made determinant one by changing an eigenvector sign, without changing its diagonal eigenvalues. The surjective chiral action realizes the resulting two rotations simultaneously. Reversing four-space orientation interchanges the two chiral blocks.

This last reversal is compatible with the oriented 14-dimensional observer geometry. A base linear change $A$ has determinant $(\det A)^{-5}$ on $\operatorname{Sym}^2(T^*X)$; its natural total-space lift has determinant $(\det A)^{-4}>0$. In particular an orthogonal reflection preserves the total-space orientation used by the reconstructed Hodge operation. No claim that it preserves the four-space chirality labels is needed. ∎

If both chiral blocks are singular and curvature is nonzero, GG1–GG6 already apply, including the pure-mixed and one-zero-block cases. The remaining argument concerns a nonsingular $W^+$ and retains all four amplitudes in (HH1).

The actual smooth metric in these normal coordinates is
$$g_{ij}(x)=\delta_{ij}-\tfrac13R_{iajb}x^ax^b
 +Q_{3;ij}(x)+Q_{4;ij}(x)+R_{5;ij}(x),\tag{HH3}$$
with arbitrary symmetric cubic and quartic terms and a smooth remainder whose derivatives through order four vanish at zero. Positivity holds on a sufficiently small neighbourhood, with no uniform radius over the jets. Throughout the geometric calculation, Z1 has $\lambda=1/5$ and $c_H,c_V>0$.

## A finite exact compatibility identity

**Proposition HH2 (finite frame cover).** For every nonzero real quadruple in (HH1), at least one of the 36 spectral frames gives incompatible selected relaxed-system rows at the numeric positive fibre
$$y_* = \operatorname{diag}(2,3,5,7),\qquad (c_H,c_V)=(2,3).$$

*Degree bound before interpolation.* BB1/GG3's argument requires only affine normal quadratic curvature and a parameter-independent observer metric/inverse at the origin. The coefficient entries have total degree at most two and the RHS entries at most three in $(p,q,r,s)$. Fifteen simplex nodes are unisolvent for the former space and 35 for the latter. Their evaluation determinants have absolute values $16$ and $84934656$. The code constructs these spaces explicitly; interpolation is not used to guess the degree bounds.

At rows $(168,337,492)$, or $(0,13,0),(1,13,1),(2,13,2)$, the equations have the form
$$a_i(p,q,r,s)\xi_{13}+N_{y_*}(n)=b_i(p,q,r,s),\quad i=0,1,2,$$
where the common $N_{y_*}$ is independent of curvature. Define the actual determinant, with its stated sign and scale,
$$P(p,q,r,s)=\det[(a_i,1,b_i)]_{i=0}^2.\tag{HH4}$$
The freshly reconstructed polynomial has degree four and 35 homogeneous terms; the output prints it in full. Its sign varies, so no pointwise positivity of $P$ is assumed.

Set
$$X=p^2+q^2+(p+q)^2,\qquad Y=r^2+s^2+(r+s)^2.$$
For each $\sigma\in S_3\times S_3$, let $\sigma$ act by permuting the two triples in (HH1). The compact certificate gives the exact polynomial identity
$$\sum_{\sigma\in S_3\times S_3}\sigma(P)\,\sigma(M)=(X+Y)^3,\tag{HH5}$$
where the single quadratic seed is
$$M=-\frac{78400}{16071864927}\left[
18385622537(p^2+r^2)+38447522213(pq+rs)
+8874071033(q^2+s^2)-1663777752pr\right].$$
Direct rational multiplication verifies all 84 degree-six coefficient slots, including the 40 zero-target slots. The helper receives the freshly computed $P$ and stores only the seed and frame permutations; it contains no saved source determinant. A negative control rejects an incorrectly doubled $P$.

For real nonzero curvature, $X+Y>0$. Hence not all $\sigma(P)$ vanish. HH1 realizes the selected frame geometrically, giving the claimed incompatible selected rows there. This is a finite exact cover, with no exceptional norm-ratio cones and no inference from evaluated signs. It supplies compatibility only; coefficient rank is established next. ∎

## New normalized rank certificates

**Proposition HH3 (rank for a nonzero chosen first eigenvalue).** At $y_*$ and $(c_H,c_V)=(2,3)$, fix $p=\epsilon\in\{1,-1\}$ and allow arbitrary real $q,r,s$. The coefficient matrix has rank 42.

*Proof.* These matrices are computed anew at Ricci zero. For each sign, ten unisolvent total-degree-two nodes in $(q,r,s)$, with evaluation determinant of absolute value $8$, verify all the following polynomial identities. This use of the hand degree bound is distinct from computing ranks at ten samples.

The 722-row cubic-invariant subsystem has four independent kernel columns
$$K_i=e_{28+i}+\frac{24d_i}{5}e_{14+i},\quad (d_0,d_1,d_2,d_3)=(2,3,5,7).\tag{HH6}$$
The even-sign base reflections preserve every diagonal electric/magnetic quadratic curvature in (HH2) and this diagonal fibre. Their characters split the columns into a 12-column trivial block, three complementary off-diagonal 6-column blocks and four horizontal 3-column blocks. Full independent coordinate reflections are not asserted. The code verifies this support and $AK=0$ at all ten nodes.

Constant rational row combinations give a left inverse for the 12-column block and two of the 6-column blocks. Each horizontal block has a constant rank-two row output and its one-dimensional part of (HH6). The remaining 6-column block, with columns $(5,12,19,26,33,40)$, has a constant rank-four output $C$ with kernel basis
$$Q=\begin{pmatrix}
0&0\\0&0\\12/5&0\\0&72/175\\1&0\\0&1
\end{pmatrix}.$$
All combinations, row indices and rational coefficients are in the compact rank JSON. To complete this block, write
$$Q_1=31q^2+22qs+31s^2,\qquad
Q_2=29(q+\epsilon)^2+2(q+\epsilon)(r+s)+29(r+s)^2.$$
For the block matrix $A_3$, four rows of $A_3Q$ are exactly

| Original row | Row of $A_3Q$ |
|---|---|
| 776 | $-27(Q_1,0)/14000$ |
| 791 | $-3(Q_2,0)/784$ |
| 958 | $81(0,Q_1)/245000$ |
| 1043 | $9(0,Q_2)/13720$ |

After dividing by these nonzero constants and adding the corresponding rows, one obtains $h_\epsilon I_2$, where
$$\begin{aligned}
h_\epsilon=Q_1+Q_2
={}&\frac{35}{3}+48(q+7\epsilon/12)^2+20s^2+28(r+s)^2\\
&+11(q+s)^2+(q+r+s+\epsilon)^2>0.
\end{aligned}\tag{HH7}$$
Thus a vector killed by the full block is first in $\ker C=\operatorname{im}Q$ and then has zero two coordinates by (HH7). This block also has full rank six. No large parameter determinant or Groebner calculation is required for the published certificate.

The invariant subsystem therefore has exact rank $12+3\cdot6+4\cdot2=38$ and kernel precisely (HH6). For $J_5=(61,229,383,523)$, the rows $(i,5,5)$, direct degree-two identities give
$$J_5K=-(41r^2+58\epsilon r+41)
\operatorname{diag}(1/700,3/1400,1/280,1/200).\tag{HH8}$$
Its scalar is
$$41(r+29\epsilon/41)^2+840/41>0.$$
These four rows complete the coefficient rank to 42. The new data and identities are valid for every real $(q,r,s)$ with $p=\epsilon$; no nonsingularity of the other eigenvalues is required by the rank certificate itself. ∎

## One selected witness through normalization and higher jets

**Proposition HH4 (same-witness transfer).** If at least one chiral block is nonsingular, HH2 and HH3 give a rank-42 incompatible witness at one and the same positive fibre of the actual metric, for every $c_H,c_V>0$ and every higher jet in (HH3).

*Proof.* Choose the nonsingular block to be $W^+$, using HH1 if necessary. HH2 chooses one spectral frame with $P\ne0$. Its first eigenvalue $p$ is nonzero. Only now normalize by GG2's Ricci-independent dilation: for $F_a(x)=ax$, put $g'=a^{-2}F_a^*g$ with $a=|p|^{-1/2}$. Then
$$R'=R/|p|,\qquad
(p',q',r',s')=(\operatorname{sign}p,q/|p|,r/|p|,s/|p|).$$
The determinant is homogeneous of degree four, so $P'=P/|p|^4\ne0$. At this same normalized numeric fibre $y_*$, HH3 gives coefficient rank 42 and HH2 gives inconsistency. Its augmented rank is therefore 43.

For clarity, GG2's observer normalization uses
$$\gamma=3/c_V>0,\qquad q_0=2c_V/(3c_H)>0,\qquad
\Phi_a(x,y)=(ax,a^{-2}y),\quad\phi_{q_0}(x,y)=(x,q_0y),$$
and
$$ (\Phi_a\circ\phi_{q_0})^*(\gamma h_g(c_H,c_V))=h_{g'}(2,3).$$
The relaxed unknowns under $h\mapsto\gamma h$ transform invertibly as $(v,\xi,n)\mapsto(v/\gamma,\xi,n/\gamma^2)$. Pullback and this scaling preserve coefficient and augmented ranks. The whole normalized witness is transported to the original positive matrix $q_0a^{-2}y_*$, followed by the chosen frame change. Rank and incompatibility are not taken from different fibres. No assertion that the original numeric $y_*$ works at all amplitudes is used.

The arbitrary higher jets are handled before transporting the witness. DD3/GG5's support proof is independent of the chosen fixed quadratic curvature and of nonzero Ricci. Cubics leave the 722-row subsystem and the three selected coefficient rows fixed. Total base inversion, which fixes every quadratic metric, and the affine cubic dependence make the selected RHS variation zero. Quartic diagonal-vertical gradient cancellation leaves those RHS rows fixed; terms of order at least five change no required observer third jet. For $J_5$, the only cubic coefficient variation multiplies $\xi_5$, so it annihilates the horizontal kernel (HH6); quartics change no coefficient. Hence the actual higher-jet system still has rank 42 and the same incompatible rows at the normalized witness. These higher jets need not have the baseline even-sign symmetry used to compute the rank blocks. ∎

## Fibre-open and full-equation consequences

**Theorem HH5 (geometric relaxed system).** Let a smooth positive Riemannian four-metric have
$$\mathrm{Ric}_g(x_0)=0,\qquad R_g(x_0)\ne0.$$
At $\lambda=1/5$ and $c_H,c_V>0$, the reconstructed Z1 relaxed system is consistent on no nonempty open subset of the positive fibre over $x_0$.

*Proof.* If both chiral blocks are singular, use GG6. Otherwise HH4 supplies an actual augmented rank-43 witness. A nonzero augmented minor at that witness is a nonzero rational function on the fibre, by S4. Consistency on any fibre-open patch would make every augmented 43-minor zero there and hence identically zero, a contradiction. All cases of nonzero Ricci-zero curvature are covered. ∎

**Corollary HH6 (full equation in the stated sector).** Impose S.0: the compact $(14,0)$ real form, gauge algebra $\mathfrak u(128)$, invariant odd field, slaved vector even field with $\det M\ne0$, one fixed reconstructed Shiab member with coefficients $(x,y)\ne(0,0)$, and $\kappa\ne0$. Use positive Z1 and $c_H,c_V>0$.

(a) For the twisted-cone family (C1), $\Upsilon=0$ has no solution on an open domain meeting the fibre over a point satisfying HH5.

(b) The same exclusion holds for every other nonzero reconstructed member under these assumptions, by S16(b). This part is a corollary of that earlier classification, not an extension of Q5 to additional members.

*Proof of (a).* Fix the actual nonzero augmented 43-minor and its coordinate fibre point from HH5. Its nonvanishing persists over a sufficiently small base neighbourhood. Since $R_g(x_0)\ne0$, shrink further so Riemann curvature remains nonzero there. A constant-field open piece would make Z1 Einstein and S5(b) would force a flat base on its open projection, a contradiction. Q5's twisted-cone reduction gives the dense variable patch with $H\ne0$ and then the forbidden relaxed-system patch, as in U5. The Ricci-zero or spectral conditions need not persist on this neighbourhood. For (b), S16(b) forces $\lambda=1/5$ for a nonzero member and forces a flat base on the projection for every non-C1 member. Such a projection cannot contain $x_0$. No arbitrary even-field or non-invariant odd-field conclusion is made. ∎

## Reproduction and boundary

Run `python -X utf8 -I computation/geometric-unity/HH_ricci_zero.py` from the repository root with the pinned dependencies. The portable default assembles the union of 45 exact node systems once in memory: 35 compatibility nodes, including the ten $p=+1$ rank nodes, and ten additional $p=-1$ rank nodes. It then checks the compact rational row data, the finite polynomial identity using the freshly computed determinant, and a scale-rejection control. No raw matrices, caches or off-repository files are needed. GG's both-singular theorem remains an explicit reviewed dependency, rather than being numerically re-established by these nodes.

The recorded fresh author run passes 9 main checks, 28 rank checks and 6 orbit checks in 76.4 seconds. Review36 independently reproduces these checks in 77.6 seconds and reruns the direct GG dependency in 21.1 seconds. The rational rank data are 94,410 bytes and the orbit data 6,395 bytes.

Review36 also checks a nonunit witness at chiral parameters $(-2,3,5,-7)$, $(c_H,c_V)=(5,7)$ and transported fibre $(28/15)y_*$, obtaining ranks $42/43$. Its normalized non-grid control applies a deterministic rule across the cubic/quartic spaces, giving 180 of 200 cubic and 311 of 350 quartic coefficients nonzero: 1,188 other coefficient entries and 636 RHS entries change, while the protected rows, $J_5K$ and ranks $42/43$ are preserved. These controls use the accepted referee engine. The review does not supply a second metric/curvature engine, independently reconstruct every imported geometric calculation, or establish physical evidence or novelty.

The complete flat quadratic origin remains a separate control with coefficient/augmented ranks $28/28$. The theorem does not exclude zero-Riemann-curvature points, assert that arbitrary higher jets at those points retain rank 28, or cover indefinite real forms. It is a pointwise curvature obstruction in the stated reconstructed sector, not a classification of all Ricci-flat metrics or of unrestricted fields.
