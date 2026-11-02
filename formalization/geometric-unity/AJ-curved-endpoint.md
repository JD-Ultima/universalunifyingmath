# AJ. A singular endpoint forces a flat base

This note classifies one of the four constant odd endpoints left by AF/AG without assuming a flat base. On positive Z1 at $\lambda=1/5$, the first Lambda3 branch at $t/\kappa=-1/400,u=0$ admits only the known zero-even solution over a flat base, at $c_V\kappa^2=150000/11$. No curvature-trace sign or analyticity assumption is imposed.

The proof first establishes a metric rigidity criterion, then applies the endpoint's grade-one equation. It combines the reconstructed geometry with a spectral elimination; it is not a theorem about arbitrary fields in the source proposal.

## Setting and inputs

Let $g$ be a smooth Riemannian four-metric and $Y=\operatorname{Met}_{+}(X)$. Use reconstructed Z1 with $\lambda=1/5$ and constant $c_H,c_V>0$. Write
$$G=\widehat{\operatorname{Ric}^Y}-\tfrac12\operatorname{scal}_Y I,$$
where the hat raises an index with the positive observer metric. Thus $G$ is a self-adjoint endomorphism of the fourteen-dimensional observer tangent space.

The imported geometric results are:

- [F4](F-zorro-reduction.md), proved in [K](K-z1-ricci-general-n.md), gives the Ricci blocks along a metric section. Replacing $g$ by $qg$ preserves its connection and hence Z1, so the same formula applies at $y=qg(x)$ for every $q>0$.
- [N1(d)](N-z1-neighbourhood.md) gives the universal identity $\widehat{\operatorname{Ric}^Y}E=-E/c$, with nonzero Euler field $E=y\partial_y$ and $c=4c_V/5$.
- [S4(a)](S-z1-off-section.md) makes the metric and curvature components rational in the fibre coordinates, with denominators powers of $\det y$.
- [S5(b)](S-z1-off-section.md) proves that Z1 Einstein on an open observer set over a definite base forces that base to be flat on its projection. Its proof uses N10's curvature pairing and its exact coefficient computation.

## AJ1. A metric rigidity criterion

Let $U\subset Y$ be nonempty and open and $h$ a real constant. Set $B=G-hI$ and $\tau=\operatorname{tr}B$. If
$$3300B^2-380\tau B+9\tau^2I=0\quad\text{on }U,\tag{1}$$
then
$$h=\frac{15}{2c_V},\qquad B=0\text{ on }U,$$
and $g$ is flat on $\pi(U)$.

*Proof.* Fix $x\in\pi(U)$. Each coordinate component of (1) is rational in $y$. It vanishes on the nonempty fibre-open set $U\cap Y_x$. Clearing powers of $\det y$ gives a polynomial vanishing on an open subset of the positive cone, so (1) holds on its entire positive fibre. Only this metric identity is continued.

Choose $g$-normal coordinates at $x$ and their induced horizontal/vertical frame along the radial ray. At $y=qg(x)$, F4 gives the endomorphism blocks
$$G=\frac{15}{2c_V}I+\frac1{qc_H}
\operatorname{diag}\!\left(\operatorname{Ric}_g^\sharp-\frac{R_g}{2}I_4,
-\frac{R_g}{2}I_{10}\right).\tag{2}$$
Taking $q\to\infty$ in (1) yields its scalar value at $B=\delta I$, where $\delta=15/(2c_V)-h$. In dimension fourteen the left side is $-256\delta^2 I$. Hence $\delta=0$.

Now $h=6/c$. Since $\operatorname{tr}G=-6\operatorname{scal}_Y$, N1(d) implies
$$\tau=-6\operatorname{scal}_Y-14h,\qquad
BE=\left(-\frac1c-\frac{\operatorname{scal}_Y}{2}-h\right)E
=\frac{\tau}{12}E.$$
Applying (1) to $E$ gives
$$\left(\frac{3300}{144}-\frac{380}{12}+9\right)\tau^2E
=\frac14\tau^2E=0.$$
Thus $\tau=0$, since $E\ne0$. Equation (1) reduces to $B^2=0$. Self-adjointness in the positive observer metric gives $B=0$. Consequently $G=hI$ and Z1 is Einstein on $U$. S5(b) supplies flatness on $\pi(U)$. ∎

No solution field, projector, inequality or choice of frame is extended to the radial ray in this argument. Smooth base metrics suffice; the rational dependence concerns the fibre variable alone.

## AJ2. Classification at the first branch's $-1/400$ endpoint

Keep the true dualized equation, compact real form $(14,0)$ with gauge algebra $\mathfrak u(128)$, fixed parallel real $W_A$, real fields and constant $\kappa\ne0$. Let a connected open observer domain remain in the single $\Lambda^3:\mathsf c,\mathsf w\Gamma$ kernel, away from every other non-vector singular factor, with regular vector slaving $\det M\ne0$. At the endpoint
$$T_o=-\frac{i\kappa}{400}\Phi^1,\qquad u=0,$$
a full solution exists exactly when $g$ is flat on the projected domain and
$$c_V\kappa^2=\frac{150000}{11}.$$
Its even field is necessarily $T_e=0$.

*Proof.* The odd field is constant, so the slaved vector part vanishes identically. Locally the remaining kernel field is $k=s\mathsf w(\omega)\Gamma$ with real $s\ne0$. [AA6–AA7](AA-conic-traces.md) give the full-equation hook condition and make the real three-form $\omega$ zero or simple. AA4 includes the full covariant derivative, including weight derivatives; its output has only grades nine and thirteen, so it cannot alter the grade-one equation.

Put $z=s^2|\omega|^2$. [AI2](AI-flat-lambda3.md) establishes the metric-independent quadratic tensor
$$[Q(k)]_1/i=180zI-\frac{320}{3}zP,$$
where $P$ is the orthogonal three-plane projector at nonzero $\omega$, and $zP=0$ at its zeros. The curved baseline contributes $G+(\kappa t+312t^2)I$. At this endpoint $\kappa t+312t^2=-11\kappa^2/20000$. Therefore, with $h=11\kappa^2/20000$,
$$B=G-hI=-180zI+\frac{320}{3}zP.$$
Its eigenvalues are $-220z/3$ on the three-plane and $-180z$ on the perpendicular eleven-plane. Its trace is $\tau=-2200z$, and its two roots are $\tau/30$ and $9\tau/110$. Multiplying their polynomial by 3300 gives (1), including at $z=0$.

AJ1 gives $h=15/(2c_V)$, $B=0$ and a flat base. The trace then forces $z=0$; the positive form norm gives $k=0$. The scale relation follows. Conversely, AI3 checks the complete equation for this constant odd field and $T_e=0$ over flat Z1 at that scale, including regularity of the other blocks. Both signs of $\kappa$ are allowed. ∎

## Verification, prior work and remaining scope

[AJ_curved_endpoint.py](../../computation/geometric-unity/AJ_curved_endpoint.py) and its [recorded output](../../computation/geometric-unity/AJ_curved_endpoint.out) provide eleven exact symbolic controls. They check the spectral elimination, radial limit, an independent section scalar-curvature identity, Euler coefficients and scale. A nonzero spectral control shows why the Euler step is necessary; an indefinite nilpotent control records the role of positivity. These controls verify algebraic combinations of imported geometric lemmas, not a second metric implementation. The universal proof is above. [Internal Review37](../../reports/reviews/REVIEW37-GU-AJ-2026-11-02.md) records the independent assessment.

Run from the repository root:

    python -X utf8 -I computation/geometric-unity/AJ_curved_endpoint.py

This closes one endpoint within the repository's reconstruction. A bounded literature comparison on November 2, 2026 considered the source draft, [Nguyen–Polya's response](https://files.timothynguyen.org/geometric_unity.pdf) and [Atik's audit](https://kleis.io/docs/papers/geometric_unity_audit.pdf), together with searches for GU singular-field and observer-metric classifications. No matching endpoint theorem was identified in that comparison. This is not an exhaustive novelty search or external peer review.

[AK](AK-lambda3-classification.md) subsequently excludes the other three endpoints without a curvature-trace sign and classifies both branches across all positive Z1 parameters. AJ retains the flatness proof for the surviving endpoint; AI3 supplies its full-equation converse. Lambda5 fields, coincident kernels, other operators, indefinite real forms, arbitrary odd fields and a separate weighted base-metric variation remain outside AJ2. No new physical evidence follows.
