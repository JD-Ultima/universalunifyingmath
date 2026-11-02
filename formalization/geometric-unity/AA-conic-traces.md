# AA. Scalar constraints on the nonconstant exterior conics

This note reduces two invariant traces of the full equation on each of the four $W_A$ exterior conics. The extra kernel norm is algebraic, and the remaining scalar equation has a simple change of variable. On every negative branch, a nonnegative curvature trace excludes a local maximum of $t/\kappa$. The untraced odd equations remain unsolved.

Code: `computation/geometric-unity/AA_conic_traces.py` → `.out`. It uses T's compact-real-form helpers and the earlier independent Clifford engine, with no imports from the J/R evaluator. Record: `metrics/geometric-unity/AA-conic-traces.yaml`. Independently reviewed in `reports/reviews/REVIEW20-GU-AA-2026-11-02.md`; all 23 assertions were reproduced. This is internal review, not external peer review. GU-054 supplies source-claim context; the results here are AA1–AA3.

## Setting and notation

Keep the true dualized equation, real $W_A$, signature $(14,0)$, compact gauge real form $\mathfrak u(128)$, the positive normalized $-\operatorname{Tr}$ pairing, constant $\kappa\ne0$ and $\det M\ne0$. The odd field is
$$T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\qquad a=t/\kappa,\quad w=u/\kappa.$$
Work on a smooth arc of one conic, away from all other non-vector singular factors. Its even kernel is one exterior copy. The parameters are:

| exterior pair | $d$ | $e$ | $p$ | $q$ | $P=pq$ |
|---|---:|---:|---:|---:|---:|
| $\Lambda^3:c,w\Gamma$ | $112/3$ | $400$ | $352/3$ | $32$ | $11264/3$ |
| $\Lambda^3:w,c\Gamma$ | $240$ | $272/3$ | $32$ | $352/3$ | $11264/3$ |
| $\Lambda^5:c,w\Gamma$ | $96$ | $832/3$ | $48$ | $80/3$ | $1280$ |
| $\Lambda^5:w,c\Gamma$ | $704/3$ | $160$ | $80/3$ | $48$ | $1280$ |

In the compact real bases of these two copies, R2's block is
$$\kappa\begin{pmatrix}A&-pw\\-qw&B\end{pmatrix},\quad
A=1+da,\quad B=1+ea,\quad w^2=f(a):=\frac{AB}{P}.$$
Put
$$h=A+B,\qquad \delta=1+224a-3328f(a),\qquad s=|da|^2,\qquad
\Delta=\operatorname{div}\operatorname{grad},\qquad \mathcal G=\operatorname{tr}_{g_Y}G^Y.$$
On a real conic, $h\ne0$: if $A+B=0$ and $AB\ge0$, both would vanish, contrary to $d\ne e$. Slaving requires $\delta\ne0$ and the other factor $1+544a-3328f(a)\ne0$.

The even field is $T_e=v+k$, with $k$ in this single non-vector kernel. J1 gives, in the real compact phases,
$$v=P_2(C)+P_{14}(B_0),\qquad C=\frac{24}{\delta}\,da,\qquad B_0=208wC.$$
There is no $du$ source. Write $N=\|k\|^2\ge0$; this is the kernel norm only.

## AA1. The two traces and T2's residual

Define
$$F(a)=a+312a^2+104f(a),\qquad
H(a)=\frac{d+e+2de\,a}{2h}=\frac{P f'(a)}{2h}.$$
For every full solution the following two equations hold:
$$0=\mathcal G+14\kappa^2F(a)+312\operatorname{div}C+1456|C|^2+H(a)N,\tag{AA1.1}$$
$$0=w\left[14\kappa^2(1+208a)-\frac{24920064}{\delta^2}s-\frac{P}{h}N\right].\tag{AA1.2}$$
The first is the $i\Phi^1$ trace; the second is the $i\Phi^1\Gamma$ trace. In particular neither contains a derivative of $k$.

The baseline residual $E_v$ of T has the explicit contraction
$$\boxed{\frac{\langle T_o,E_v\rangle}{\kappa}
=a\mathcal G+14\kappa^2(1+312a)(a^2+f)
+\frac{7488a}{\delta}\Delta a
+\frac{7488[-112a+3328(af'-f)]}{\delta^2}s.}\tag{AA1.3}$$
T2 requires this expression to equal $N/2$. No $\Delta u$ or derivative of $k$ appears. Since $\kappa$ is constant, $s=|dt|^2/\kappa^2$ and $\Delta a=\Delta t/\kappa$ give the same formula in the original variable.

*Proof.* The independent differentiated cubic gives, for arbitrary real vector covectors,
$$\langle i\Phi^1,Q(v)\rangle=1456|C|^2,\qquad
\langle i\Phi^1\Gamma,Q(v)\rangle=-208\langle B_0,C\rangle.$$
The derivative traces are $312\operatorname{div}C$ and zero. Also $\mathcal D_\xi(i\Phi^1\Gamma)=0$ for every $\xi$. T1's skew-symbol argument makes both non-vector derivative traces zero. Cubic symmetry and preservation of the vector isotype make both traces of $2\mathcal B(v,k)$ zero. These statements concern traces; the derivative and cross tensors themselves need not vanish.

Let $N_1,N_2$ be the squared norms of the two grade components of $k$. Self-adjointness and the kernel equations give $AN_1=BN_2$. Differentiating the cubic in the two invariant odd directions gives
$$\langle i\Phi^1,Q(k)\rangle=\tfrac12(dN_1+eN_2)=H N,\qquad
\langle i\Phi^1\Gamma,Q(k)\rangle=-\frac{Pw}{h}N.$$
These formulas also hold at a conic endpoint $w=0$. They use no derivative of a chosen exterior-basis field. The invariant-background traces are $\mathcal G+14\kappa^2F$ and $14\kappa^2w(1+208a)$, proving (AA1.1–2).

Finally, $\operatorname{div}C=24(\Delta a/\delta-\delta's/\delta^2)$ and $\delta'=224-3328f'$. Pair the baseline traces with $a,w$ and use $w^2=f$. This gives (AA1.3). T2 supplies the required norm equality. ∎

## AA2. Algebraic norm and a scalar equation

On an arc with $w\ne0$, the Gamma trace fixes the kernel norm:
$$\boxed{N=\frac{h}{P}\left[14\kappa^2(1+208a)-\frac{24920064}{\delta^2}|da|^2\right].}\tag{AA2.1}$$
Substituting into (AA1.1) gives
$$\boxed{\Delta a-\frac{\delta'}{2\delta}|da|^2
=-\frac{\delta}{7488}\left[\mathcal G+\kappa^2\mathscr R(a)\right],}\tag{AA2.2}$$
where the explicit quadratic polynomial is
$$\mathscr R(a)=14[a+312a^2+104f(a)]+7f'(a)(1+208a).$$
These identities are valid at $a=0$ as well: their derivation uses the direct $\Phi^1$ trace, without dividing T2 by $a$. At $w=0$, use the undivided traces (AA1.1–2); the divided norm formula is not an extra endpoint condition unless justified by continuity from $w\ne0$.

*Proof.* Equation (AA2.1) is (AA1.2) divided by $w$. The norm contribution in (AA1.1) becomes
$$HN=7\kappa^2f'(1+208a)-\frac{12460032f'}{\delta^2}|da|^2.$$
Combining it with the derivative and vector-quadratic terms leaves gradient coefficient
$$\frac{7488(112-\delta')-12460032f'}{\delta^2}
=-\frac{3744\delta'}{\delta^2},$$
which proves (AA2.2). The script also checks that (AA1.3) minus $N/2$ is $a$ times the substituted first trace, as required by T2. ∎

On the right real branch, $a\ge\max(-1/d,-1/e)$, both $h$ and $1+208a$ are positive. Positivity of $N$ therefore yields the necessary gradient bound
$$|da|^2\le\frac{14\kappa^2(1+208a)\delta^2}{24920064}.\tag{AA2.3}$$
This bound and the scalar equation do not establish existence.

## AA3. A sign obstruction on every left branch

Suppose $\mathcal G\ge0$ and the field remains on a left branch
$$a\le\alpha:=\min(-1/d,-1/e).$$
The four endpoints are $-3/112,-3/272,-1/96,-1/160$, respectively. There is no full solution attaining an interior local maximum of $a$. In particular a full solution cannot reach its left endpoint $w=0$.

*Proof away from the endpoint.* Every left branch has
$$a<-1/312,\quad 1+208a<0,\quad f\ge0,\quad f'<0,\quad h<0,\quad \delta<0.$$
Indeed $\alpha<-1/312$, $1+224\alpha<0$ and $1+208\alpha<0$ for all four listed endpoints. Thus $F=a(1+312a)+104f>0$ and $\mathscr R>0$. Define an increasing scalar coordinate $z(a)$ by
$$z'(a)=(-\delta(a))^{-1/2}.$$
Equation (AA2.2) gives
$$\Delta z=\frac{\sqrt{-\delta}}{7488}[\mathcal G+\kappa^2\mathscr R]>0.$$
A smooth function with strictly positive Laplacian cannot have an interior local maximum. Since $z$ is increasing, the same holds for $a$.

*Endpoint.* If $w=0$ on this branch, $a=\alpha$ is a local maximum, so $da=0$ and $\Delta a\le0$. In (AA1.1), $312\operatorname{div}C=7488\Delta a/\delta\ge0$, $F>0$, $H>0$, and $N\ge0$. Its RHS is then strictly positive. This uses the undivided first trace and does not assume that (AA2.1) holds at the endpoint. ∎

For flat-base Z1 at $\lambda=1/5,c_V>0$, $\mathcal G=105/c_V>0$, so this bounded obstruction applies. It allows nonconstant motion without a local maximum on a left arc; it does not exclude every nonconstant field there. No compactness, integration or boundary condition is used.

## Evidence and remaining equations

The exact script checks both vector traces, all fourteen zero Gamma symbols, all four compact exterior pairs, their symbolic kernels, both kernel traces, cross-trace cancellation, and both currents on each pair. It then verifies the norm elimination, T2 identity and rational branch-sign premises symbolically. Its output prints $f$, $\delta$ and $\mathscr R$ for every conic. All 23 assertions passed in about two seconds.

The untraced odd equation is still
$$E_v+{*}^{-1}\circledast(d_0k)+2\mathcal B(v,k)+Q(k)=0.$$
In particular $d_0k$ is retained in its other odd components. The scalar reduction does not classify exterior tensor fields or solve their differential equation. Intersections, $\det M=0$, other real forms, non-invariant odd fields and a variable Shiab operator remain outside this result.
