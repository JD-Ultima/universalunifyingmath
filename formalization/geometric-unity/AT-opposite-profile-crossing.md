# AT. The exceptional opposite-conic profile level

At the zero of AR's leading profile coefficient, boundedness selects one analytic scalar germ with a forced norm and slope. An actual field cannot have a critical point at this level. Consequently any actual crossing is regular and obeys an exact scale bound. The positive scalar germs below that bound do not establish existence of a metric, five-form or full solution.

Code: [AT helper](../../computation/geometric-unity/AT_opposite_profile_crossing.py) and [output](../../computation/geometric-unity/AT_opposite_profile_crossing.out). Record: [AT metadata](../../metrics/geometric-unity/AT-opposite-profile-crossing.yaml). [Independent internal Review50](../../reports/reviews/REVIEW50-GU-AT-2026-11-03.md) accepts AT1–AT3 and independently reproduces the helper. External peer review remains outstanding.

## Setting and scalar equation

For statements about actual fields, retain all [AR](AR-opposite-lambda5-geometry.md) hypotheses: a smooth positive four-base and positive Z1 with constant $c_H,c_V>0$; compact observer signature $(14,0)$ and $\mathfrak u(128)$; the true bosonic equation and fixed parallel real $W_A$; real fields and constant $\kappa\ne0$; $\det M\ne0$; and an entire connected open domain in the single opposite $\Lambda^5:\mathsf w,\mathsf c\Gamma$ conic, away from coincident non-vector factors. AR proves $w\ne0$, $\lambda=1/5$ and its other necessary geometric restrictions. Put
$$c_0=4c_V(1-4\lambda)>0,\quad C=d\chi=\frac{24}{\delta(a)}da,
\quad c=|C|^2,\quad t=\kappa^2c_0>0,\quad y=c_0c.$$
All derivatives and norms use the actual observer metric. Use AR's coefficients
$$\begin{gathered}
\delta=-\tfrac8{15}(208a+1)(880a+3),\qquad
U=\tfrac{(208a+1)(21840a+113)}{480},\\
V=-\tfrac{64(402896a+1719)}{285},\qquad
W=\tfrac{1024(53-676a)}{285},\qquad D=V-\tfrac{11}{2}W,
\end{gathered}$$
and set
$$P=V'+W'-2WD/\delta,\qquad Q=U'-2WU/\delta,\qquad
R=-12(W-48)/\delta.$$
Primes here denote derivatives in $a$. On a regular patch, AR5b is
$$Dy'+Py+tQ+R=0.\tag{AT0}$$
The exceptional level and its coefficients are
$$a_*=-\frac{491}{26416},\quad D=d_1(a-a_*),\quad
d_1=-\frac{21978112}{285},\quad\delta(a_*)=-\frac{4939648}{241935}\ne0.$$
This is on the real left arc. The scalar statements below require only (AT0) on a sufficiently small interval about $a_*$ avoiding the zeros of $\delta$.

## AT1. Unique bounded scalar germ

For every fixed real $t$, (AT0) has exactly one real analytic germ at $a_*$. Every bounded solution on either punctured side is its restriction and extends uniquely to that germ.

*Proof.* Write $x=a-a_*$ and
$$p(x)=P(a_*+x)/d_1,\qquad q(x;t)=-[tQ(a_*+x)+R(a_*+x)]/d_1.$$
These functions are analytic, and the equation becomes
$$xy'+p(x)y=q(x;t),\qquad p_0=p(0)=153/127>0.$$
Define
$$v(x)=\exp\!\left(\int_0^x\frac{p(s)-p_0}{s}\,ds\right),\qquad
z(x)=\int_0^1u^{p_0-1}v(ux)q(ux;t)\,du.$$
The first integrand has a removable singularity; $v$ is analytic and nonzero near zero, with $v(0)=1$. The second integral converges because $p_0>0$. If $vq=\sum b_nx^n$, uniform convergence on a smaller interval gives
$$z(x)=\sum_{n\ge0}\frac{b_n}{n+p_0}x^n.$$
Thus $y=z/v$ is analytic. Integration of the derivative of $u^{p_0}v(ux)q(ux;t)$ gives $xz'+p_0z=vq$, proving the equation. This establishes convergence as well as the formal coefficient recurrence.

The difference of two solutions is $A_+x^{-153/127}/v(x)$ for $x>0$ or $A_-|x|^{-153/127}/v(x)$ for $x<0$. A nonzero such term is unbounded. Boundedness therefore forces the difference to vanish on either side. ∎

The forced value and slope are
$$y_0(t)=\frac{95(3241880613-6380808980t)}{259536139978752},\tag{AT1a}$$
$$y_1(t)=\frac{103523628088901644t+54385536931208595}{788642454954332160}.\tag{AT1b}$$
Indeed $y_0=q(0;t)/p_0$ and $y_1=[q'(0;t)-p'(0)y_0]/(1+p_0)$. With
$$\tau_*:=\frac{3241880613}{6380808980},$$
the germ is positive on a small two-sided interval when $0<t<\tau_*$. For $t>\tau_*$ its value is negative. At equality its value is zero and its slope is positive, so only its right side is positive near zero. Hence no positive bounded two-sided crossing exists for $t\ge\tau_*$. This is a scalar existence and uniqueness statement.

## AT2. No actual critical point at the exceptional level

An actual AR solution has no point with $a=a_*$ and $C=0$.

*Proof.* Suppose such a point exists and put $Z=C^\sharp$, $A=(\nabla Z)$ there. The smooth identity $\mathrm{Ric}(Z)=-Z/c_0$ from AR5a differentiates at $Z=0$ to
$$ (\mathrm{Ric}+I/c_0)A=0.$$
AR0 at this point gives $\mathrm{Ric}=24A+LI$, where
$L=2\operatorname{tr}A+\kappa^2U(a_*)/6$. Therefore
$$24A^2+(L+1/c_0)A=0.\tag{AT2a}$$
Self-adjointness for the positive observer metric implies that $A$ has only the eigenvalues zero and
$\eta=-(L+1/c_0)/24$. It cannot vanish: otherwise $\mathrm{Ric}=\kappa^2U(a_*)I/6>0$, contrary to N1's negative Euler Ricci direction. A double zero root would also force $A=0$. Thus $m=\operatorname{rank}A\in\{1,\ldots,14\}$ and
$$\eta=-\frac{\kappa^2U(a_*)/6+1/c_0}{24+2m}<0.\tag{AT2b}$$
N1 gives a nonzero Euler field $E$ with $|E|^2=c_0$, $\nabla_EE=0$ and $\mathrm{Ric}(E)=-E/c_0$. Consequently $A(E)=\eta E$.

Follow its actual short integral curve $\gamma(s)$ through the point, using parallel transport to compare vectors. Smoothness gives
$$Z(\gamma(s))=\eta sE+O(s^2),\quad
y(\gamma(s))=\eta^2c_0^2s^2+O(s^3),$$
$$a(\gamma(s))-a_*=\frac{\delta(a_*)\eta c_0}{48}s^2+O(s^3).\tag{AT2c}$$
The last coefficient is positive. For small $s>0$, the curve consists of regular points and $a$ is strictly increasing. AQ's local profiles therefore give, after inverting $a(s)$, a solution of (AT0) for $x>0$, bounded with limit zero. AT1 applies from this one side; no prior profile regularity at the critical point is assumed. It forces $t=\tau_*$, and the Taylor ratio in (AT2c) gives
$$c_0\eta=\frac{\delta(a_*)y_1(\tau_*)}{48}
=-\frac{53607532260267}{929045787488000}.$$
But (AT2b) would then require
$$m=-\frac{109962574640284}{53607532260267}<0,$$
which contradicts its positive rank. Only the scalar ODE germ is analytically extended; observer fields remain on their actual curve. ∎

## AT3. Necessary condition for an actual crossing

Every actual point with $a=a_*$ is regular and satisfies
$$0<\kappa^2c_0<\tau_*.$$
Its local profile is the unique positive germ of AT1, with the norm and slope (AT1a–b).

*Proof.* AT2 excludes a critical point. On a regular patch, AQ gives a smooth positive profile $c(a)$ on a two-sided interval about its attained value. AT1 supplies the conclusion. ∎

The result does not prove that any permitted germ is the Gram tensor of a real five-form, meets AS's full compatibility or is realized by a Z1 metric. It does not exclude critical points at other values of $a$ or classify the remaining conic.

## Evidence

Run `python -X utf8 -I computation/geometric-unity/AT_opposite_profile_crossing.py`. The recorded run passes 13 exact controls in 0.091 seconds with Python 3.12 and SymPy 1.14.0; Review50 independently reproduces them in 0.092 seconds. The standalone helper checks the profile normalization, forced coefficients, integrating-factor algebra, scale signs, critical Hessian trace and Euler-curve ratio. It imports no previous helper or metric engine. The integral convergence, bounded uniqueness, positive spectral argument and passage along an actual Euler curve are hand proofs.
