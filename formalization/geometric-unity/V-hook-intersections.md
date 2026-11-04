# V. Intersections on the positive hook lines and a grade-five obstruction

For the later exclusion of the entire first pure hook over flat-base positive Z1 at $\lambda=1/5$, see [AX](AX-hook-curvature.md). The other pure cases and all coincident factors retain their separate scope.

This note resolves the finite intersection set used in T5. It also separates the remaining odd equations on the first hook branch and gives an even-kernel field that passes both invariant trace equations but fails the full odd equation.

Code: `computation/geometric-unity/V_hook_intersections.py` → `.out`. It uses T's compact-real-form helpers and the earlier independent Clifford engine `xcheck_round14.py`, not R or J's evaluator. Record: `metrics/geometric-unity/V-hook-intersections.yaml`. Internal review: `reports/reviews/REVIEW16-GU-V-2026-11-02.md`.

## Setting

Keep T's fixed metric and spin connection, signature (14,0), compact gauge real form $\mathfrak u(128)$, true field equation, $\kappa\ne0$ and $\det M\ne0$. Fix the real member $W_A$, whose coefficient in I/J's code normalization is $-i$. Thus
$$T_o=i(t\Phi^1+u\Phi^1\Gamma),\qquad t,u\in\mathbb R.$$
Use the positive normalized pairing $\langle A,B\rangle=-\sum_a\operatorname{tr}(A_aB_a)$, $\operatorname{tr}(1)=1$.

Assume $t/\kappa=a$ is one of the six positive hook ratios in R6. T5 forces $u$ to be constant on each connected full-solution domain. Since $dt=0$, the vector-sector slaved field is zero; every solution of the even equation is $T_e=k\in\mathcal W\cap\ker L_{T_o}$. None of the catalogue entries below is asserted to solve the odd equation.

## V1. Exact intersection catalogue

Put $r=(u/\kappa)^2$. The four non-hook singular factors give
$$\begin{aligned}
r_{3c}&=(1+\tfrac{112}3a)(1+400a)/(11264/3),&
r_{3w}&=(1+240a)(1+\tfrac{272}3a)/(11264/3),\\
r_{5c}&=(1+96a)(1+\tfrac{832}3a)/1280,&
r_{5w}&=(1+\tfrac{704}3a)(1+160a)/1280.
\end{aligned}$$
The suffixes refer to the paired copies $(c,w\Gamma)$ or $(w,c\Gamma)$ in R2.

| $a$ | hook copy | $r_{3c}$ | $r_{3w}$ | $r_{5c}$ | $r_{5w}$ |
|---|---|---|---|---|---|
| $3/136$ | $H_2\Gamma$ | $15531/3255296$ | $963/191488$ | $6413/369920$ | $1617/73984$ |
| $1/40$ | $H_4\Gamma$ | $29/5120$ | $343/56320$ | $2023/96000$ | $103/3840$ |
| $3/104$ | $H_6\Gamma$ | $13203/1903616$ | $14523/1903616$ | $441/16640$ | $7373/216320$ |
| $3/88$ | $H_6$ | $12075/1362944$ | $13635/1362944$ | $1081/30976$ | $639/14080$ |
| $1/24$ | $H_4$ | $1219/101376$ | $43/3072$ | $113/2304$ | $2231/34560$ |
| $3/56$ | $H_2$ | $1413/78848$ | $11931/551936$ | $4773/62720$ | $1273/12544$ |

**Proposition V1.** Each row has four distinct positive values, none colliding with $\det M=0$. At $u=\pm\kappa\sqrt{r_{3c}}$ or $\pm\kappa\sqrt{r_{3w}}$, the non-vector kernel is the listed hook copy plus one $\Lambda^3$ copy (dimension 364). At the two analogous $\Lambda^5$ values it is the hook copy plus one $\Lambda^5$ copy (dimension 2002). At $u=0$ only the hook copy is singular. No $\Lambda^7$ or second hook copy occurs on these lines.

Thus a full solution confined to one of these lines has at most nine possible constant odd pairs: $u=0$ and the eight signed intersection values. This is a list of necessary candidates, not an existence result.

**Proof.** R2 is homogeneous in $(\kappa,t,u)$, so dividing by nonzero $\kappa$ reduces the kernel calculation to $\kappa=1$, including negative $\kappa$. At $\theta=0$, the $\Lambda^3$ and $\Lambda^5$ matrices split into the two paired two-by-two blocks. Their determinants are $\kappa^2(1+da)(1+ea)-pq\,u^2$, with $(d,e,pq)$ given in the formulas above. At a zero each paired block has rank one, since $1+da$ and $1+ea$ are positive. Exact rational comparisons in the table separate all four values and the two vector roots $(1+224a)/3328$, $(1+544a)/3328$.

The other paired block is regular because its squared root differs. The $\Lambda^7$ eigenvalues are $\kappa(1+176a)$ and $\kappa(1+208a)$, both nonzero. The six hook ratios are distinct, so exactly one hook copy is singular. At $u=0$, all non-hook diagonal entries are nonzero. These computations also show $\det M\ne0$ at every listed candidate. T5 supplies constancy of $u$. ∎

## V2. An additional norm requirement at each intersection

At a nonzero intersection write $k=k_H+k_\Lambda$ according to V1. Let $(d,e,pq)$ label its paired block, and set $A=1+da$, $B=1+ea$. Every full solution satisfies
$$\|k_\Lambda\|^2
=14\kappa^2(1+208a)\frac{A+B}{pq}>0.\tag{V2.1}$$
Its hook norm is then constrained by
$$\|k_H\|^2
=2a\operatorname{tr}_{g_Y}G^Y
 +28(1+312a)(\kappa^2a^2+u^2)-\|k_\Lambda\|^2\ge0.\tag{V2.2}$$
This strengthens T2's nonnegative-total-norm requirement at those points. It does not settle the remaining tensor equations.

**Proof.** Write $L=\kappa D+uC$ on the paired $\Lambda$ block, where $D$ has eigenvalues $A,B$ on its two orthogonal copies. Here $C=2\mathcal B(\cdot,i\Phi^1\Gamma)$ is self-adjoint for the positive pairing and preserves the isotypic splitting. The hook part has $Ck_H=0$. The Γ-partner derivative current vanishes by $\mathcal D_\xi(i\Phi^1\Gamma)=0$ for $W_A$. Pairing the full odd equation with $i\Phi^1\Gamma$ therefore gives
$$14u\kappa(1+208a)+\tfrac12\langle k_\Lambda,Ck_\Lambda\rangle=0.$$
The kernel equation gives $u\langle k_\Lambda,Ck_\Lambda\rangle=-\kappa\langle k_\Lambda,Dk_\Lambda\rangle$, so
$$\langle k_\Lambda,Dk_\Lambda\rangle=28u^2(1+208a).$$
Pair each of its two component kernel equations with that component. Self-adjointness makes the cross terms equal, giving $A\|k_1\|^2=B\|k_2\|^2$. Hence
$$\|k_\Lambda\|^2=14u^2(1+208a)(1/A+1/B).$$
Use $u^2=\kappa^2AB/(pq)$ for (V2.1), and T2/T4 for (V2.2). ∎

## V3. The grade-five equation is not automatic

Now restrict to the pure first hook branch
$$t=3\kappa/136,\qquad u=0,\qquad k\in H_2\Gamma.$$
Its real compact fields have Clifford grade 12. The full odd equation is exactly the following three graded equations:
$$\begin{array}{ll}
\text{grade 1:}&i\big(G^Y+ (\kappa t+312t^2)g_Y\big)_{ab}e^a\otimes\gamma^b+[Q(k)]_1=0,\\
\text{grade 5:}&[Q(k)]_5=0,\\
\text{grade 11:}&{*}^{-1}\circledast(d_0k)=0.
\end{array}\tag{V3.1}$$
In particular, the derivative of $k$ is still an equation to satisfy. It has not been set to zero.

**Grade proof.** Write a grade-12 element as $B_2\Gamma$, with $B_2$ a bivector. Here $\Gamma^2=-1$; Γ commutes with even grades and anticommutes with odd grades. For bivectors $B_2,C_2$ and a trivector $A_3$,
$$[B_2\Gamma,C_2\Gamma]=-[B_2,C_2],\qquad
\{A_3,B_2\Gamma\}=[A_3,B_2]\Gamma.$$
The first expression has grade 2; the second has grade 11, since a bivector commutator preserves grade 3. Thus $k\wedge k$ has grade 2, and the derivative term has grade 11. Its covariant derivative preserves this Clifford decomposition because the spin connection preserves Γ and grade.

The first quadratic-gradient term is a grade-3 anticommutator with the grade-2 product, giving grades 1 and 5. In the adjoint-gradient term the intermediate element is $A'_3\Gamma$; its commutator with $B_2\Gamma$ is
$$[A'_3\Gamma,B_2\Gamma]=-\{A'_3,B_2\},$$
again of grades 1 and 5. Scalar phases in the real form do not change these grade statements. The curvature and invariant background terms have grade 1. ∎

For a concrete obstruction, define the hook tensor $h(p,q;r)$, with distinct indices, by
$$h_p=\gamma_q\gamma_r,\qquad h_q=\gamma_p\gamma_r,$$
and all other components zero. The contraction $\sum_a\gamma^ah_a$ is zero, so $h\in H_2$. Let
$$k_0=i\big(h(0,1;2)+h(0,3;4)\big)\Gamma.$$
It is a real anti-Hermitian field in the kernel, with $\|k_0\|^2=4$. Direct exact evaluation of the cubic gradient gives
$$[Q(k_0)]_5=4i\,w(E_{1234})\ne0,\qquad
w(E_{1234})_a=[\gamma_aE_{1234}]_5.\tag{V3.2}$$
Every nonzero real multiple of this element fails the grade-five equation, regardless of its derivative jet. No derivative term or curvature term can cancel that grade. This rejects these specified kernel elements, not arbitrary hook fields.

The distinction from the scalar tests can be realized on the flat-base Z1 model. Take $\lambda=1/5$, $\kappa=1$, $c_V=36414/30599$ and $c_H=2c_V$. F7 gives
$$\operatorname{tr}_{g_Y}G^Y=152995/1734.$$
At $t=3/136$, $u=0$, the field $k_0$ solves the even equation and passes both invariant odd trace equations. It still fails (V3.2). A local orthonormal frame defines such a smooth field; it need not be parallel, and all its derivative terms remain in grade 11. This is an explicit failure of sufficiency for the scalar tests, not a new solution.

## Evidence and remaining work

The script checks all 24 squared intersections, vector-block regularity, multiplicities, and the norm formulas with exact rationals. Four separate Clifford controls on the first hook line cover all four paired block types. It verifies (V3.2), anti-Hermitian phases, the even-kernel condition, both scalar traces in the flat-model control, and all 14 derivative symbols of the displayed element.

The catalogue reduces the six line branches to finitely many constant odd pairs. It does not solve the full odd equations for those pairs. V3.1 displays the pure first-hook tensor and derivative constraints; AX later excludes that entire branch over flat-base positive Z1 at $\lambda=1/5$. The other pure hook cases and curved-base first-hook problem remain open. The coincident hook–$\Lambda^3/\Lambda^5$ cases, nonconstant motion on other singular curves, $\det M=0$, and indefinite real forms also remain open.
