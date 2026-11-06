# BJ. Closed H6 orbit averages and their nonlinear obstruction

Finite real signed-measure averages of BC's closed $H_6$ seed cannot solve the full bosonic equation on any nonempty open domain of the flat positive Z1 model at $\lambda=1/5$. This includes every finite real linear combination of its $SL(5,\mathbb R)$ translates and a constructive 121-dimensional family of non-affine closed fields. The measure must be fixed globally and have finite total variation; arbitrary closed hook fields remain outside this result.

The proof uses the [global model AM](AM-global-flat-models.md), [closure identity BA](BA-hook-closure.md), [closed seed and exact nonlinear moments BC](BC-affine-hooks.md), and [spherical value modules BF](BF-h6-curvature-module.md). The new [scalar helper](../../computation/geometric-unity/BJ_sphere_boundary.py) and [output](../../computation/geometric-unity/BJ_sphere_boundary.out) check only the optional sphere constants in BJ4. The group-action, measure, analytic and nonlinear arguments below are hand proofs. Internal [Review72](../../reports/reviews/REVIEW72-GU-BI-BJ-2026-11-06.md) and [Review73](../../reports/reviews/REVIEW73-GU-BI-BJ-2026-11-06.md) accept BJ1–BJ4 with contribution-aware disclosures. External peer review remains outstanding.

## Setting

Retain BA's true dualized equation, fixed parallel real $W_A$ with V's normalization, observer signature $(14,0)$, compact $\mathfrak u(128)$, real smooth fields and constant real $\kappa\ne0$. The base is flat and positive, $\lambda=1/5$, and $c_H,c_V>0$ are constant. Restrict to the pure branch
$$
T_o=i\frac{3\kappa}{104}\Phi^1,\qquad T_e=iB\Gamma,\qquad B\in H_6.
$$
Here $B$ is a seven-covariant tensor, the Clifford generators are Hermitian and $\Gamma^2=-1$. By R/V the even equation has exactly this hook kernel and a regular vector block. BA makes $d_\nabla B=0$ necessary, and its derivative term has grade seven, disjoint from grade one.

First normalize $c_H=2,c_V=1$ and use the full simply connected model
$$
Y=S/K,\qquad S=SL(5,\mathbb R),\quad K=SO(5),\quad o=K.
$$
The natural tensor action is $g\cdot B=(g^{-1})^*B$. These connected isometries preserve orientation, $\Gamma$ and the natural $W_A$ Clifford operations; a spin-lift sign acts trivially on the bosonic adjoint tensors. Let $b$ be BC's global affine-invariant real closed $H_6$ field and $B_*=b(o)$. Its affine group is the connected positive hyperplane parabolic
$$
P=\left\{(\det A)^{-1/5}
\begin{pmatrix}A&a\\0&1\end{pmatrix}:\det A>0\right\}\subset S.
$$
AM's matrix realization verifies its affine action. It is transitive on $Y$, so $S=PK=KP$. Its compact intersection $L=P\cap K=SO(4)$ fixes the last unit vector $e_4$; thus $K/L=S^4$, with an oriented normal, rather than an unoriented axis. Since $P$ fixes the entire field $b$,
$$S\cdot b=K\cdot b$$
as sets of global tensor fields.

## BJ1. Exact global profiles for finite-measure averages

For $n=ke_4\in S^4$, put
$$b_n=k\cdot b,\qquad B(n)=b_n(o)=kB_*.$$
These are well-defined independently of the representative $k$: $L$ fixes the entire seed field. Let $\mu$ be any fixed finite real signed Borel measure on $S^4$, with $|\mu|(S^4)<\infty$, and define
$$
B_\mu(y)=\int_{S^4}b_n(y)\,d\mu(n).\tag{BJ1.1}
$$
The smooth orbit map descends to the compact parameter space $K/L$. Every local tensor derivative is uniformly bounded there on a compact coordinate neighborhood. Integration against finite total variation therefore gives a smooth real hook field and commutes with covariant differentiation. Consequently
$$d_\nabla B_\mu=0.$$
No unique measure representation or injectivity is required. Atomic and singular measures are included; atomic measures include all finite real linear combinations of $S$-translates of $b$.

At $y=so$, evaluate in the $s$-isometric frame by applying $s^{-1}$. Factor $s^{-1}k=k'p$ with $p\in P$. The oriented hyperplane normal of this coset is
$$
n'=\frac{(s^{-1}k)^{-T}e_4}{\|(s^{-1}k)^{-T}e_4\|}
=\frac{s^Tn}{\|s^Tn\|}.
$$
Because $p$ fixes the global field, $[(s^{-1})\cdot b_n](o)=k'B_*=B(n')$. Hence the exact profile is
$$
\boxed{[(s^{-1})\cdot B_\mu](o)
=\int_{S^4}B\!\left(\frac{s^Tn}{\|s^Tn\|}\right)d\mu(n).}\tag{BJ1.2}
$$
There is no parabolic character and no Jacobian: the field is $P$-fixed and the original integration variable $n$ is unchanged. The transpose agrees with AM's inverse-congruence action $Q\mapsto s^{-T}Qs^{-1}$.

BF's componentwise spherical argument gives vector-valued homogeneous harmonic polynomials $P_3,P_5$ such that
$$
B(n)=P_3(n)+P_5(n),\qquad B(-n)=-B(n),\qquad \|B(n)\|=\|B_*\|.
$$
This is a statement about the value map, requiring no global field-Casimir polynomial for $b$. The integrand in (BJ1.2) is therefore
$$
\frac{P_3(s^Tn)}{(n^Tss^Tn)^{3/2}}
+\frac{P_5(s^Tn)}{(n^Tss^Tn)^{5/2}}.\tag{BJ1.3}
$$

## BJ2. A collapse limit excludes the full equation

There is an axis $u\in S^4$ whose equator has zero total-variation measure:
$$|\mu|(\{n:n\cdot u=0\})=0.$$
Indeed normalized sphere measure $\sigma$ and Tonelli give
$$
\int_{S^4}|\mu|(u^\perp)\,d\sigma(u)
=\int_{S^4}\sigma(n^\perp)\,d|\mu|(n)=0.
$$
Each great equator has zero sphere measure, so almost every axis is suitable. Total variation is essential here; signed cancellation is not being used.

Choose $k\in K$ with $ke_4=u$ and set
$$
s_{a,u}=k\operatorname{diag}(a,a,a,a,a^{-4})k^T,\qquad a>0.
$$
Every finite $a$ defines an actual positive point of the same metric manifold and an isometric transport to $o$. For $n\cdot u\ne0$,
$$\frac{s_{a,u}^Tn}{\|s_{a,u}^Tn\|}\longrightarrow\operatorname{sign}(n\cdot u)u\qquad(a\downarrow0).$$
The integrand is bounded by $\|B_*\|$. Dominated convergence with respect to $|\mu|$, oddness of $B$, and (BJ1.2) yield
$$
[(s_{a,u}^{-1})\cdot B_\mu](o)\longrightarrow\beta_\mu(u)B(u),
\qquad\beta_\mu(u)=\int_{S^4}\operatorname{sign}(n\cdot u)\,d\mu(n)\in\mathbb R.
\tag{BJ2.1}
$$
The limit is in one finite tensor fiber after isometric transport. No field or metric is extended to the degenerate endpoint, and no derivative-jet limit is needed.

The necessary grade-one equation on the normalized model is
$$
\ell g_Y+Q_1(iB_\mu\Gamma)/i=0,\qquad
\ell=\frac{15}{2}+\frac{15}{52}\kappa^2>0.
$$
If it held globally, its transported limit would be
$$\ell g_o+\beta_\mu(u)^2Q_1(iB(u)\Gamma)/i=0.\tag{BJ2.2}$$
The pointwise quadratic is continuous, homogeneous of degree two and $K$-equivariant. BC's exact metric-relative seed moments satisfy
$$q_I(B_*)-q_H(B_*)=\frac{5700375}{4096}\ne0.$$
Since $B(u)=kB_*$, compare the $k$-transforms of these metric-unit directions in (BJ2.2). Their difference forces $\beta_\mu(u)=0$. The remaining equation then forces $\ell=0$, a contradiction. Thus every field (BJ1.1) fails the full equation somewhere on $Y$. A zero hemisphere coefficient or even a zero field is included: it fails the positive baseline. No injectivity of the hemisphere transform or nonzero coefficient assertion is needed.

## BJ3. Any open restriction and all positive scales

Every fixed-measure field (BJ1.1) is real analytic. To prove this without assuming analyticity of arbitrary solutions, choose a local positive upper Cholesky section on AM's determinant-one positive matrix model:
$$Q=R(Q)^TR(Q),\qquad s(Q)=R(Q)^{-1},\qquad\det R=1.$$
Positive-pivot recursion makes this section analytic. Near a fixed positive $Q_0$, $n^Ts(Q)s(Q)^Tn$ has a uniform positive lower bound for real $n\in S^4$. Shrink a complex coordinate neighborhood so its real part stays uniformly positive for every such $n$. Equation (BJ1.3) then has a single holomorphic branch of each denominator power and uniform local holomorphic bounds. Its numerators and tensor-frame conversion are analytic too. Finite-total-variation integration preserves these extensions; no density or absolute continuity of $\mu$ is required.

The grade-one residual is consequently an already-defined analytic tensor on connected $Y$. If the restriction of this particular $B_\mu$ solved the full equation on any nonempty open subset, the analytic identity theorem would make its grade-one residual vanish everywhere, contradicting BJ2. The domain need not reach the collapse curve. No arbitrary smooth field or its domain is continued.

For general positive scales let $h=c_V$ and $D_s(x,y)=(sx,y)$, $s=\sqrt{c_H/(2c_V)}$. BC gives
$$g(c_H,c_V)=hD_s^*g(2,1).$$
Transfer the fields and group actions by covariant pullback. Constant homothety preserves the connection, hook and closure kernels. For a fixed seven-covariant tensor its metric-relative quadratic moments scale by $h^{-7}$, so the limiting difference becomes
$$h^{-7}\beta_\mu(u)^2\frac{5700375}{4096}.$$
The actual baseline remains strictly positive,
$$\ell_h=\frac{15}{2h}+\frac{15}{52}\kappa^2>0.$$
The same contradiction and analyticity argument apply. Thus for all positive $c_H,c_V$ and either sign of nonzero $\kappa$, no nonempty open restriction of a finite real signed-measure orbit average solves the declared pure-branch full equation.

## BJ4. A constructive 121-dimensional closed family

This supplement uses classical compact averaging and matrix-coefficient orthogonality. For their induced-representation and Peter–Weyl setting, see [Woit's notes](https://www.math.columbia.edu/~woit/LieGroups-2012/peterweyl.pdf), pp. 1, 4–7. No new general harmonic-analysis method is claimed.

BF gives $B_*=v_{18}+v_{40}$ in two orthogonal irreducible real modules $W_{18},W_{40}$ of dimensions $d_{18}=30,d_{40}=91$. Put $N_\eta=\|v_\eta\|^2>0$ and $v_\eta(n)=kv_\eta$ when $n=ke_4$. For constant real $w=w_{18}+w_{40}$ define
$$
F_w=\sum_{\eta=18,40}\frac{d_\eta}{N_\eta}
\int_{S^4}\langle w_\eta,v_\eta(n)\rangle b_n\,d\sigma(n).\tag{BJ4.1}
$$
These are smooth-density cases of (BJ1.1). Haar invariance gives $k\cdot F_w=F_{kw}$. At $o$, the same-module averaging operator is positive self-adjoint and equivariant, hence scalar by real irreducibility; its trace is $N_\eta$. It is therefore $(N_\eta/d_\eta)I$. The cross-module operator is zero because it intertwines distinct Casimir eigenvalues. Consequently
$$F_w(o)=w.$$
Thus (BJ4.1) is an injective $K$-equivariant right inverse to evaluation, constructing a genuine global 121-dimensional family of real closed hooks. It does not assert that every closed extension of a value in this module has these profiles.

For the axis $e_4$, normalized zonals and coordinate density are
$$
Z_3(t)=\frac{7t^3-3t}{4},\qquad Z_5(t)=\frac{33t^5-30t^3+5t}{8},\qquad
d\sigma_t=\frac34(1-t^2)\,dt.
$$
Their signed moments are $-1/16$ and $3/128$. Hence the coefficient in (BJ2.1) for (BJ4.1) is
$$
\alpha(w)=-\frac{15}{8N_{18}}\langle w_{18},v_{18}\rangle
+\frac{273}{128N_{40}}\langle w_{40},v_{40}\rangle.
$$
In particular $\alpha(v_{18})=-15/8$, $\alpha(v_{40})=273/128$ and $\alpha(B_*)=33/128$. Along the axis-$e_4$ curve the affine seed transports constantly to $B_*$, whereas $F_{B_*}$ has this different limiting multiplier despite the same value at $o$. These optional constants distinguish the actual profiles; they are not the finite-measure no-go proof.

Every nonzero $F_w$ is non-affine. Indeed if $P$ fixed $F_w$, its transported value along that curve in $P$ would be constantly $w$. The limit forces $w=\alpha(w)B_*=cB_*$. But $\alpha(cB_*)=33c/128$, so $c=33c/128$ and $w=0$, hence $F_w=0$ by the injective construction. This proves the non-affine assertion without a classification of global $S$-modules.

The scalar helper checks zonal normalization, orthogonality, the two sphere eigenvalues, signed moments and these multipliers. It imports no previous helper and verifies no global profile, dominated-convergence argument or nonlinear equation. Those are the hand arguments BJ1–BJ3.

## Scope

The measure is fixed, real and of finite total variation on the full model sphere. Point-dependent measures, positive-order distributions such as derivatives of Dirac masses, infinite-variation expressions and unspecified completions of the orbit span are not covered. There is no claim about arbitrary smooth closed $H_6$ fields, arbitrary 135-module mixtures, all fields with values in a pointwise curvature-kernel module, other extensions sharing one value, curved bases, other metric parameters, coincident kernels, quotient descent or other operators/real forms. No physical conclusion or literature priority is asserted.
