# AR. Geometry of the remaining Lambda5 branch

The opposite single $W_A$ Lambda5 conic has no zero-$u$ points in a smooth solution domain. A solution requires $\lambda=1/5$, base Ricci rank at most one, and avoidance of the local-parallel set. A Ricci–Hessian change of scalar transfers the earlier geometric relaxed-system obstructions to this branch. The Euler and conservation identities also give a necessary one-variable profile equation. These restrictions do not establish existence or exclude every remaining field.

The [helper](../../computation/geometric-unity/AR_opposite_lambda5_geometry.py) and [output](../../computation/geometric-unity/AR_opposite_lambda5_geometry.out) record 34 new exact controls. [Metadata](../../metrics/geometric-unity/AR-opposite-lambda5-geometry.yaml) records dependencies. [Independent internal Review48](../../reports/reviews/REVIEW48-GU-AR-2026-11-03.md) accepts AR1–AR5 and independently reproduces the helper. External peer review remains outstanding.

## Setting

Retain [AQ](AQ-opposite-lambda5-conservation.md)'s positive-Z1, compact-real, true-equation and single-kernel assumptions, now on an entire connected open opposite-conic domain, initially allowing $w=0$. The base is a smooth positive four-metric; $c_H,c_V>0$ and $\kappa\ne0$ are constants, $\lambda<1/4$, $\det M\ne0$, and coincident non-vector factors are excluded. The odd field is $T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma)$ and the even field is the complete slaved vector plus the declared five-form kernel. All chart weights and derivatives remain included.

Use
$$w^2=f(a)=\frac{(1+\frac{704}{3}a)(1+160a)}{1280},\qquad
\delta=1+224a-3328f(a),\qquad C=\frac{24}{\delta}da=d\chi.$$
Write $g_Y$ for the observer metric, $g$ for the base metric, $c=|C|^2$, $h=\operatorname{div}C$, and
$$c_0=4c_V(1-4\lambda)>0,\quad T=\kappa^2U+Vc,\quad
L=2h+T/6+Wc/12.$$
The functions $U,V,W$ are AQ's explicit polynomials. Its metric equation and trace reversal are
$$G-24(\nabla C-hg_Y)+Tg_Y+WC\otimes C=0,\qquad
\operatorname{Ric}=24\nabla C+Lg_Y-WC\otimes C.\tag{AR0}$$
Here and below unlabelled curvature tensors are those of $g_Y$. [N1](N-z1-neighbourhood.md) gives the Euler field $E$ with $|E|^2=c_0$, $\nabla E=P_H/2$, and $\operatorname{Ric}(E,\cdot)=-E^\flat/c_0$.

## AR1. No endpoint or locally constant patch

Every full solution in this setting has $w\ne0$ everywhere. Also $\{da\ne0\}$ is open dense.

*Proof.* [AN](AN-opposite-endpoint-obstruction.md) excludes an open subset with $w=0$: its connected components would have one constant odd endpoint. Thus $w\ne0$ is dense. The tensor on the left of (AR0) is smooth on the actual domain, including $w=0$, since $\delta\ne0$. AQ proves its vanishing where $w\ne0$, so continuity proves (AR0) everywhere. This extends a smooth metric identity on the actual domain, not a divided form-norm formula.

At $w=0$, differentiating $w^2=f(a)$ twice gives $da=0$ and $\nabla^2a=2\,dw\otimes dw/f'(a)$. Consequently
$$\nabla C=\frac{48}{\delta f'}dw\otimes dw.$$
The exact endpoint data are

| $a$ | $\delta$ | $f'$ | $\delta f'$ | $U$ |
|---|---|---|---|---|
| $-1/160$ | $-2/5$ | $-7/120$ | $7/300$ | $47/3200$ |
| $-3/704$ | $1/22$ | $7/120$ | $7/2640$ | $877/185856$ |

Thus $\nabla C$ is positive semidefinite, $h\ge0$ and $C=0$. Trace reversal in (AR0) gives
$$\operatorname{Ric}=24\nabla C+(2h+\kappa^2U/6)g_Y>0,$$
contradicting the negative Euler Ricci direction. This also covers $dw=0$ contacts.

On the two real nonzero-$w$ arcs, $a<-1/160$ or $a>-3/704$, the factored $U=(208a+1)(21840a+113)/480$ is strictly positive: both its roots lie in the intervening forbidden interval. A constant-$a$ open patch would have $C=h=0$ and $\operatorname{Ric}=\kappa^2U g_Y/6>0$, again impossible. The critical set therefore has empty interior. ∎

## AR2. An exact Ricci–Hessian bridge

On a small $a$-interval choose $\mu(a)>0$ and $\psi(a)$ by
$$\mu_a=\frac{576k(a)}{\delta(a)}\mu,\qquad
\psi_a=\frac{576\mu}{\delta},\qquad
k(a)=-\frac{W}{576}=\frac{16(676a-53)}{2565}.$$
Then $\psi_a\ne0$, so $\mu$ is a smooth function of $\psi$, and
$$d\psi=24\mu C,\qquad
\nabla^2\psi=\mu(\psi)\operatorname{Ric}+\nu g_Y,\qquad \nu=-\mu L.\tag{AR2}$$
Indeed $d\mu=k\,d\psi$, so the additional Hessian term is $576k\mu C\otimes C=-\mu WC\otimes C$. This argument does not require $da\ne0$ or that $\nu$ depend only on $\psi$.

The exact scalar identities are
$$\mu_\psi'=k,\qquad \mu\mu_{\psi\psi}''=\frac{169\delta}{23085}=Q(k),\qquad
a=\frac{2565k+848}{10816},$$
$$Q(k)=-\frac{(57k+20)(141075k+48668)}{200070}.$$
Commuting (AR2) gives precisely [S.0](S-z1-off-section.md)'s geometric relaxed system with
$$v=\operatorname{grad}\psi/\mu=24C^\sharp,\qquad
\xi=kv^\flat,\qquad n=\operatorname{grad}\nu/\mu.$$
No original slaved-vector field corollary is imported; this is a direct identification of the compatibility equations.

## AR3. The metric parameter and base Ricci rank

A solution requires $\lambda=1/5$ and $\operatorname{rank}\operatorname{Ric}_g\le1$ throughout its base projection.

**Auxiliary continuation.** Fix a fibre meeting an actual regular patch. Every augmented 43-minor of the relaxed 42-column system vanishes there and hence identically as a rational function of its positive fibre matrix $y$. Thus its full-rank points are consistent. Choose a radial full-rank point where its unique $v$ is nonzero, and one 42-row minor nonzero there. This minor's nonzero locus meets the actual fibre-open patch. On a smaller product of a base neighbourhood and an actual fibre patch, its Cramer fields agree with the actual $(v,\xi,n)$ and $|v|^2\ne0$.

Recover the scalar by $k=\xi(v)/|v|^2$ and the affine inverse in AR2. Finite base derivatives of these Cramer expressions remain rational in $y$, with smooth base-jet coefficients. They agree with the actual derivatives on that product. Clearing denominators at each nearby base point therefore continues both $da=\delta v^\flat/576$ and the full metric equation (AR0), including its base derivatives, to a neighbourhood of the selected radial point. This product argument supplements S4's fibre-derivative formulation.

The reconstructed finite scalar cannot hit a root $a_0$ of $\delta$ there. Along every small positive-fibre path it solves the locally Lipschitz ODE $a'=\delta(a)v^\flat(\dot\gamma)/576$. Uniqueness would force the constant root solution on a fibre neighbourhood, hence make the rational scalar identically $a_0$ on that fibre. This contradicts $\delta\ne0$ on its original patch. Thus the reconstruction gives genuine local Ricci–Hessian germs near a radial interval. Neither $\omega$ nor real conic inequalities are extended to that interval.

**Excluding $\lambda\ne1/5$.** Put $b_Z=(1-5\lambda)/(c_V(1-4\lambda))$ and $B=4b_Zc_0=16(1-5\lambda)$. If $b_Z\ne0$, S13 gives full rank at radial points except at most one; consistency forces an Einstein base, and its $v_E=2b_Z-e^{-\theta}\sigma_g/2$ is nonzero except at most one further point, where $\sigma_g=\mathrm{scal}_g/c_H$. Apply the reconstruction on an interval avoiding these exceptions. Replacing $g$ by $e^{\theta_0}g$ places that interval in [Q7](Q-z1-remaining-members.md)'s notation without changing the connection.

Q7 with $\sigma_g=0$ forces $k=1/B$ on a nontrivial $\psi$-interval, contradicting $Q(k)\ne0$. If $\sigma_g\ne0$, its persistence equations give
$$X=\rho^2=1-\frac1{Bk},\qquad
\rho=1+\frac{4k^2X}{Q(k)},\qquad
\rho=e^{-\theta}\sigma_g/(4b_Z).$$
Here $k$ varies on an interval. The expression for $\rho$ is rational in $k$, but its square cannot equal $(Bk-1)/(Bk)$: the latter has a simple pole at zero, whereas rational squares have even pole orders. The equivalent cleared polynomial has nonzero constant coefficient $947429689600B$. Hence $b_Z=0$.

**Bounding base Ricci.** At $\lambda=1/5$, a projected rank-at-least-two Ricci point would give an open base rank condition meeting the dense projected regular set. S13 then supplies full-rank radial reconstruction. Consistency forces nonzero Einstein Ricci and $v_E=-e^{-\theta}\sigma_g/2\ne0$. Q6(c) forces $k=0$ along the radial interval, while the continued differential identity gives $E(a)=\delta c_0v_E/576\ne0$. The affine inverse for $a$ makes this impossible. Continuity of the two-by-two Ricci minors gives the rank bound on the whole base projection.

The geometric Q6/Q7/S13 input retains its previously recorded computer-assisted N6(f) dependency. The new scalar controls do not rebuild that input or supply a new metric engine.

## AR4. Transferred obstructions and the parallel set

The bridge supplies relaxed consistency at every actual point, including spatial critical points. Consequently its domain cannot meet a fibre excluded by the geometric results [CC](CC-diagonal-weyl.md), [EE](EE-zero-electric-weyl.md), [FF](FF-joint-weyl.md) or [HH](HH-ricci-zero.md): CC/EE/FF concern their specified nonzero rank-one Ricci/Weyl families, while HH concerns Ricci-zero, nonflat points. Their exact hypotheses and higher-jet scopes are unchanged. The transfer uses their no-fibre-open-consistency statements, not their original field corollaries.

The domain also avoids [S2](S-z1-off-section.md)'s actual local-parallel set $\Pi$. If it met $\Pi$, S3/F2.4 would identify a neighbourhood with the section of a positive parallel base metric. Q6 forces that metric Einstein on an open base neighbourhood; AR3's rank bound makes its Ricci tensor zero. HH excludes every nonzero-curvature point there, so the base is flat. The geometric flat-base lemma Q12 at $\lambda=1/5$ then forces $d\psi=0$ on the observer neighbourhood, contradicting AR1. No unrestricted field equivalence is asserted.

## AR5. Euler variation and a necessary scalar profile

The set $\{E(a)\ne0\}$ is open dense. Moreover
$$\operatorname{Ric}(C^\sharp)=-C^\sharp/c_0\tag{AR5a}$$
everywhere, with the identity trivial at $C=0$.

*Proof.* Set $p=C(E)$ and $C_H=C\circ P_H$. Differentiating $p$ and using (AR0) gives
$$dp=-\frac{1/c_0+L}{24}E^\flat+\frac{Wp}{24}C+\frac12C_H.$$
If $p$ vanished on an open patch, this identity would give $C_H=0$ and $L=-1/c_0$. The needed all-fibre horizontal trace formula is
$$\operatorname{tr}_H\operatorname{Ric}
=\frac1{c_H}\operatorname{tr}(y^{-1}\operatorname{Ric}_g)
-\frac4{c_0}-2|\mathcal A|^2,\qquad
\mathcal A(X,Y)=\tfrac12[X^H,Y^H]^V.$$
To derive it, split the vertical horizontal-second-fundamental tensor into its symmetric part $b_H$ and antisymmetric part $\mathcal A$. Koszul gives the ordered horizontal curvature trace
$c_H^{-1}\operatorname{tr}(y^{-1}\operatorname{Ric}_g)-|m|^2+|b_H|^2-3|\mathcal A|^2$; the ordered symmetric/antisymmetric cross sum cancels. F's orthogonal-splitting formula gives the mixed trace $\operatorname{div}m-|b_H|^2+|m|^2+|\mathcal A|^2$. N1 gives $m=-2E/c_0$ and $\operatorname{div}m=-4/c_0$. Norms sum over ordered horizontal pairs; no Riemannian-submersion assumption on Z1 is used.

Since $C_H=0$, $\operatorname{tr}_H\nabla C=-C(m)=2p/c_0=0$. Thus (AR0) implies the metric-only equality $\operatorname{tr}_H\operatorname{Ric}=-4/c_0$. S4 continues it rationally along the positive fibre. At $y=g$, $\mathcal A=0$; its first fibre variation also makes $d|\mathcal A|^2=0$, so the equality forces $\operatorname{Ric}_g=0$. Positivity then gives $\mathcal A=0$ for every $y>0$. Every base curvature endomorphism is therefore skew for every positive $y$, which forces it to vanish. The base is flat on the projection. N2's horizontal Ricci block and (AR0) now imply $(\nabla C)_{HH}=0$. N1's symmetric horizontal Koszul map from the vertical covector $C$ is injective, so $C=0$ on that open patch, contradicting AR1. Hence $p\ne0$ is dense.

On regular patches AQ gives $\operatorname{Ric}(C)=\zeta(a)C$. Self-adjointness with the Euler identity yields $(\zeta+1/c_0)p=0$. Density and smoothness prove (AR5a), including critical points. ∎

On a regular patch, write $c=c(a)$ as in AQ. Eliminating $h$ between (AR5a) and AQ's Bianchi identity gives the implicit linear equation
$$\boxed{D c'+\left[V'+W'-\frac{2W}{\delta}D\right]c
+\kappa^2\left[U'-\frac{2W}{\delta}U\right]
-\frac{12(W-48)}{c_0\delta}=0,}\tag{AR5b}$$
$$D=V-11W/2=-\frac{64(343408a+6383)}{285}.$$
Its zero $a=-491/26416$ lies on the admissible left arc; it is retained, not divided away. Here $c>0$ on a regular patch, and the actual form/Gram constraints remain required. The profile equation supplies neither arbitrary initial-data admissibility nor a solution of the remaining form and metric equations.

## Evidence and boundary

Run `python -X utf8 -I computation/geometric-unity/AR_opposite_lambda5_geometry.py`. The recorded 34 controls check scalar transformations, endpoints, Euler transport, ordered frame contractions and the profile identity. Universal covariant, rational-continuation and density arguments are hand proofs; the source geometric certificates retain their own evidence. No old suite, rank certificate or metric engine was rerun for this helper.

The remaining opposite-conic system is still open. [AT](AT-opposite-profile-crossing.md) subsequently excludes critical points at the profile equation's exceptional level; [AU](AU-opposite-ricci-spectrum.md) gives spectral and finite critical-value restrictions. These conclusions apply to an entire declared single-kernel domain; they do not classify a field touching or leaving that stratum, coincident kernels, $\det M=0$, other operators or indefinite real forms. No new physical prediction, global solution or literature-first claim is asserted.
