# T. The residual required by an extra non-vector even field

For the current singular-field results and open questions, see the [research synthesis](../../reports/SYNTHESIS.md#geometric-unity).

Date: 2026-11-01, overnight continuation. This note extends R10 to nonconstant invariant odd fields. It derives a necessary scalar condition, excludes two $W_A$ singular lines under a curvature sign condition, and rules out nonconstant motion along its six positive hook lines.

Code: `computation/geometric-unity/T_singular_residual.py` → `.out`. It uses the earlier independent Clifford implementation `xcheck_round14.py`, with no imports from J, R, or the original field-equation evaluator. Result record: `metrics/geometric-unity/T-singular-residual.yaml`. Claims: GU-054 (context); gaps GU-G14, GU-G15. Internal cross-review: [T1–T5](../../reports/reviews/REVIEW13-GU-T-2026-11-02.md). No physical prediction is made.

## Setting

Fix the metric and Levi-Civita spin connection on $U\subset Y$, and one member of H7's parallel four-parameter Shiab family. Use the true, dualized equation of D5 and J.0 in signature (14,0), with compact gauge real form $\mathfrak u(128)$. All fields below are real in this gauge real form. Assume $\kappa\ne0$ and $\det M(t,u)\ne0$ on $U$.

The odd field is invariant but need not be constant:
$$T_o=i(t\Phi^1+u\Phi^1\Gamma),$$
where $t,u$ denote I.0's real primed coefficients. Let $v:=T_e^{\mathcal V}=\Sigma(dt,du)^\top$ be the uniquely slaved field in the vector isotypic component $\mathcal V$ (J1). Write $\mathcal W:=\mathcal V^\perp$, the sum of the other even isotypic components of R1. Both bundles are $\nabla^0$-parallel.

Define the **baseline residual**
$$E_v:=Y_{\rm odd}(T_o,v).$$
The subscript records the chosen even field; it does not mean a projection of the residual. In particular, $E_v$ may have all the odd components listed in J2. No baseline field equation is assumed unless stated.

Use the positive pointwise pairing
$$\langle A,B\rangle=-\sum_a\operatorname{tr}(A_aB_a),\qquad\operatorname{tr}(1)=1.$$
Here $\operatorname{tr}$ is the normalized spinor matrix trace; $\operatorname{tr}_{g_Y}$ below denotes a metric tensor trace. A different positive normalization rescales every pairing and norm together. Let $\mathcal B$ denote I2's symmetric polarization of $Q$, so $Q(A)=\mathcal B(A,A)$. The cubic potential makes $\langle A,\mathcal B(B,C)\rangle$ symmetric in all three fields.

## T1. Two orthogonality identities

For any $k\in\Gamma(\mathcal W)$,
$$\langle T_o,{*}^{-1}\circledast(d_0k)\rangle=0,
\qquad \langle T_o,2\mathcal B(v,k)\rangle=0.\tag{T1}$$

**Proof.** Put $\mathcal D_\xi A={*}^{-1}\circledast(\xi\wedge A)$. D2 and R4 give
$$\langle T_o,\mathcal D_\xi w\rangle=-\langle w,\mathcal D_\xi T_o\rangle=0
\quad(w\in\mathcal W),$$
because $\mathcal D_\xi T_o\in\mathcal V$. In a normal frame, each derivative $\nabla^0_a k$ still lies in $\mathcal W$, so the first identity follows term by term. This calculation does not differentiate $T_o$; $t,u$ may vary.

For the second identity, cubic symmetry gives
$$\langle T_o,2\mathcal B(v,k)\rangle
=\langle k,2\mathcal B(T_o,v)\rangle
=\langle k,(L_{T_o}-\kappa)v\rangle=0.$$
R2 says $L_{T_o}$ preserves $\mathcal V$, which is orthogonal to $\mathcal W$. ∎

## T2. Residual and kernel norm

Every field solving the even equation has, by R4,
$$T_e=v+k,\qquad k\in\mathcal W\cap\ker L_{T_o}.$$
For such a field the following identity holds, whether or not its odd equation is zero:
$$\left\langle T_o,Y_{\rm odd}(T_o,v+k)\right\rangle
=\langle T_o,E_v\rangle-\frac\kappa2\|k\|^2.\tag{T2.1}$$
Consequently every full solution satisfies
$$\boxed{\frac{\langle T_o,E_v\rangle}{\kappa}=\frac12\|k\|^2\ge0.}\tag{T2.2}$$

**Proof.** Pair $L_{T_o}k=\kappa k+2\mathcal B(k,T_o)=0$ with $k$. Cubic symmetry gives
$$\langle T_o,Q(k)\rangle=-\frac\kappa2\|k\|^2.$$
The odd equation at $v+k$ differs from the baseline by
$$Y_{\rm odd}(T_o,v+k)-E_v
={*}^{-1}\circledast(d_0k)+2\mathcal B(v,k)+Q(k).$$
Apply T1 and the preceding pairing. ∎

Equation (T2.2) is necessary, not sufficient. A scalar residual of the allowed sign supplies no solution of the other odd components. It only fixes the pointwise size that an extra field would need.

## T3. Rigidity for a fixed invariant odd field

If $\langle T_o,E_v\rangle=0$ throughout $U$, every full solution with this same $T_o$ has $k=0$. Thus it can only be the baseline, which must then satisfy $E_v=0$.

In particular, if the slaved baseline $(T_o,v)$ already solves the full equation, it admits no other real even field with the same invariant $T_o$. Unlike R10, $T_o$ and $v$ may be nonconstant. This is uniqueness of the even field with the odd field fixed; it does not classify possible odd fields.

**Proof.** T2 and positivity give $k=0$ at every point. Since this holds on $U$, its derivatives also vanish, leaving the baseline odd equation. ∎

No integration, boundary condition, analyticity, or compactness of the manifold is used in T1–T3. The positive pairing of the compact gauge real form is essential.

## T4. Two nonconstant $W_A$ singular strata are excluded

Take the real member $W_A$, and suppose $t$ is constant. Its source has no $du$ column (J1(d)); hence $v=0$ even if $u$ varies. I6 in the real Euclidean variables gives
$$Y_{\rm odd}(T_o,0)=i\left(G^Y_{ab}e^a\otimes\gamma^b+F_1\Phi^1+F_{13}\Phi^1\Gamma\right),$$
$$F_1=\kappa t+312t^2+104u^2,\qquad F_{13}=u(\kappa+208t).$$
The two invariant copies have pairing 14 and are mutually orthogonal. Therefore
$$\langle T_o,E_0\rangle
=t\operatorname{tr}_{g_Y}G^Y+14(\kappa+312t)(t^2+u^2).\tag{T4.1}$$

If $\operatorname{tr}_{g_Y}G^Y\ge0$ and $a:=t/\kappa<-1/312$, then
$$\frac{\langle T_o,E_0\rangle}{\kappa}
=a\operatorname{tr}_{g_Y}G^Y+14(1+312a)(t^2+u^2)<0.$$
This contradicts T2. Hence there is no full solution with this constant $t$, arbitrary smooth $u$, and arbitrary real even field, on an open set where the stated metric and $\det M\ne0$ conditions hold.

In particular this excludes both $W_A$ non-vector $\Lambda^7$ singular lines
$$t=-\kappa/176,\qquad t=-\kappa/208$$
under the same curvature-trace condition. R2 places a kernel in the $c$ or $w$ copy of $\Lambda^7$ along these lines for every $u$. The argument excludes additions in the entire non-vector kernel, including any extra components at intersections of singular factors.

Here $\det M\ne0$ is automatic: for $d=176$ or $208$ its factors are
$$\left(\kappa^2(1-224/d)-3328u^2\right)
\left(\kappa^2(1-544/d)-3328u^2\right)>0.$$
Both factors are strictly negative for real $u$ and $\kappa\ne0$.

The curvature-trace condition is an additional geometric hypothesis, not a consequence of the operator calculation. For the flat-base Z1 model at $\lambda=1/5$, $c_V>0$, F7 gives $G^Y=15g_Y/(2c_V)$ and $\operatorname{tr}_{g_Y}G^Y=105/c_V>0$, so the exclusion applies. This example concerns varying odd fields on a singular stratum, beyond R10's already solved constant backgrounds.

## T5. The positive $W_A$ hook lines admit no nonconstant $u$

Let $U$ be connected and take $W_A$ with constant $t/\kappa$ equal to one of
$$\frac3{136},\quad\frac1{40},\quad\frac3{104},\quad\frac3{88},\quad\frac1{24},\quad\frac3{56}.$$
These are R6's six positive hook singular lines. Every full solution under the standing hypotheses has **constant $u$** on $U$. This statement needs no sign condition on $\operatorname{tr}_{g_Y}G^Y$.

**Proof.** Again $v=0$. Put $z=i\Phi^1\Gamma$. The hook blocks in R2 are independent of $u$ for $W_A$, so $\mathcal B(k,z)=0$ for every hook field $k$. Cubic symmetry then gives
$$\langle z,Q(k)\rangle=0$$
whenever $k$ lies in the sum of the hook isotypes. Also $\mathcal D_\xi z=0$ for every $\xi$ by J1's zero $du$ column. The skew symbol therefore gives $\langle z,{*}^{-1}\circledast(d_0k)\rangle=0$. On a pure-hook kernel the odd Γ-partner trace consequently requires
$$0=\langle z,E_0\rangle=14u(\kappa+208t).\tag{T5.1}$$
The coefficient $\kappa+208t$ is nonzero on all six lines, so $u=0$ at every such point.

For fixed $t,\kappa$, the $\Lambda^3$ and $\Lambda^5$ determinants of $W_A$ are genuine quartics in $u$: their leading coefficients are $126877696/9$ and $1638400$ (R6(b)). Their real roots form a finite set $F$. The $\Lambda^7$ blocks are invertible because $t/\kappa>0$. Thus, at a point with $u\notin F$, the non-vector kernel is entirely hook, and (T5.1) forces $u=0$. Every full solution has
$$u(U)\subset F\cup\{0\}.$$
A continuous function on a connected set with finite image is constant. ∎

If the constant value is outside $F$, it must be zero. This does **not** exclude full solutions with constant $u=0$ and a hook field, or constant nonzero values at coincident kernel factors. Their odd equations remain open unless another result applies. T5 excludes nonconstant motion along these lines; it does not solve their constant-field cases.

## Checks and remaining questions

The exact script checks the compact real-form phases, cubic symmetry, both identities in T1, the kernel pairing, and the scalar coefficient in T4. It reproduces R10's zero-residual control at $t=-1/400$, $u=0$, $c_V=150000/11$. Four cases on the two singular lines check both signs of $\kappa$. The diagonal $G^Y$ used for those coefficient checks is algebraic test data, not an asserted observer metric. For T5 it checks the zero $\mathcal B(k,z)$ and Γ-current on all six hook copies, the Γ-residual coefficient, the six distinct hook ratios, and the quartic determinants' nonzero leading coefficients.

An independent internal cross-review checked T1–T5, including the nonconstant-field derivative argument, the cubic factor of two, T4's symbolic scalar expansion, and all six T5 ratios. This remains an internal mathematical result under the stated reconstruction.

The result does not cover signature (7,7), $\det M=0$, non-invariant odd fields, or a variable/nonparallel Shiab operator. Nonconstant motion on general singular curves remains open. The positive hook lines illustrate a limit of T2's sign test: their scalar residual is positive when $\operatorname{tr}_{g_Y}G^Y\ge0$, so that test alone does not exclude them. T5 forces $u$ to be constant there; the remaining odd tensor equations for those constant values still need to be solved or obstructed.
