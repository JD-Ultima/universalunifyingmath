# AA. Traces and derivatives on the exterior conics

This note reduces two invariant traces of the full equation on each of the four $W_A$ exterior conics. The extra kernel norm is algebraic, and the remaining scalar equation has a simple change of variable. On every left branch, a nonnegative curvature trace excludes a local maximum of $t/\kappa$. For the two $c/w\Gamma$ conics, it also gives the derivative term explicitly through weighted exterior derivatives and codifferentials. The untraced odd equations remain unsolved.

Code: `computation/geometric-unity/AA_conic_traces.py` → `.out`, with supporting derivative identities in `AA_exterior_derivative.py` → `.out`. Both use T's compact-real-form helpers and the earlier independent Clifford engine, with no imports from the J/R evaluator. Record: `metrics/geometric-unity/AA-conic-traces.yaml`. AA1–AA3 were independently reviewed in `reports/reviews/REVIEW20-GU-AA-2026-11-02.md`; all 23 assertions were reproduced. AA4–AA5 were independently reviewed in `reports/reviews/REVIEW23-GU-AA-DERIV-2026-11-02.md`, including all 14 supporting assertions. This is internal review, not external peer review. GU-054 supplies source-claim context; the results here are AA1–AA5.

## Setting and notation

All derivatives, form degrees and norms below are on the fourteen-dimensional observer domain $U\subset Y$, with metric $g_Y$.

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

## AA4. The derivative on the two $c/w\Gamma$ conics

Use $m=3,5$ for form degree; $p,q$ retain their off-diagonal block meanings in the table. Write $\mathsf c,\mathsf w$ for R1's contraction and exterior lifts, so the scalar $w=u/\kappa$ is distinct from $\mathsf w$. For a real $m$-form $\omega$, let
$$\varepsilon_3=1,\qquad\varepsilon_5=i,\qquad
k=\varepsilon_m[\mathsf c(\alpha)+\mathsf w(\beta)\Gamma],\qquad
\alpha=r_m\omega,\quad\beta=s_m\omega.$$
Here $r_m,s_m$ are real scalar kernel weights. The subscript distinguishes them from AA1's gradient scalar $s=|da|^2$. Both displayed lifts have the same compact phase $\varepsilon_m$: their grades are $(2,10)$ for $m=3$ and $(4,8)$ for $m=5$.

On the first or third row of the table, the kernel conditions are
$$Ar_m-pws_m=0,\qquad Bs_m-qwr_m=0.$$
Two smooth kernel charts are $(r_m,s_m)=(pw,A)$ and $(B,qw)$. At least one is nonzero at every conic point, including $w=0$, since $A,B$ cannot both vanish. A real $m$-form $\omega$ then parametrizes the single kernel copy locally. The other two table rows, $w/c\Gamma$, are not included in this derivative reduction.

Let $\delta_\nabla=-\sum_j\iota_{e_j}\nabla^0_{e_j}$ be the Riemannian codifferential, distinct from the scalar $\delta(a)$. Put $\mathcal D(k)={*}^{-1}W_A(d_0k)$ in the real-member conventions used above. The full covariant derivative is included. Then
$$\boxed{\begin{aligned}
\mathcal D(k)=i\varepsilon_m\big[&2m\,\mathsf w(d\alpha)-2(m-2)\mathsf c(\delta_\nabla\alpha)\\
&+\{2(12-m)\mathsf w(d\beta)-2(14-m)\mathsf c(\delta_\nabla\beta)\}\Gamma\big].
\end{aligned}}\tag{AA4.1}$$
Thus the coefficient pairs $(d,\delta_\nabla)$ are $(6,-2),(18,-22)$ for $m=3$ and $(10,-6),(14,-18)$ for $m=5$, for the $c$ and $w\Gamma$ copies respectively.

Variable weights must be differentiated:
$$d(r_m\omega)=r_m d\omega+dr_m\wedge\omega,\qquad
\delta_\nabla(r_m\omega)=r_m\delta_\nabla\omega-\iota_{\operatorname{grad}r_m}\omega,\tag{AA4.2}$$
and likewise for $s_m$. On a conic arc with $w\ne0$, a sign of $w=\pm\sqrt{f(a)}$ is fixed locally and the weights may be treated as functions of $a$, giving $dr_m=r_m'(a)da$. At a zero-$w$ endpoint, $w$ need not be a smooth function of $a$. Use a smooth conic chart and the general differentials in (AA4.2); no division by $w$ is needed for (AA4.1).

*Proof.* The lifts and $\Gamma$ are Levi-Civita parallel. In a normal frame the real-member operator is
$$\mathcal D(k)_a=i\sum_{b,c}\{\gamma_{abc},\nabla^0_bk_c\},$$
where $\gamma_{abc}$ is the alternating Clifford triple. For $\mathsf c(\omega)$ the anticommutator has grades $m-2,m+2$. For $\mathsf w(\omega)\Gamma$, move $\Gamma$ to the right: its anticommutator becomes a commutator before multiplication by $\Gamma$, with the same two pre-$\Gamma$ grades.

The decomposition
$$T^*\otimes\Lambda^m=\Lambda^{m+1}\oplus\Lambda^{m-1}\oplus\mathcal H_m$$
has one Cartan summand $\mathcal H_m$. For $m=3,5$, it does not occur in $T^*\otimes(\Lambda^{m-2}\oplus\Lambda^{m+2})$ or its right-$\Gamma$ image. The middle-grade splitting when $m+2=7$ does not introduce $\mathcal H_5$. Equivariance therefore annihilates that summand and factors the derivative through $d\omega$ and $\delta_\nabla\omega$.

The coefficients follow by counting complete basis-jet outputs. On $e^0\otimes(e^1\wedge\cdots\wedge e^m)$, only $d\omega$ is nonzero. The $c$ copy has $m$ contributing terms; the $w\Gamma$ copy has $14-m-2$. Each contributes twice the exterior lift, and all other output components vanish by Clifford parity. At $m=5$, the two equivalent grade-seven exterior copies are fixed by this full-output calibration and irreducibility of $\Lambda^6$; absence of $\mathcal H_5$ alone would not fix their coefficients. On $e^0\otimes(e^0\wedge\cdots\wedge e^{m-1})$, only the codifferential is nonzero. The counts are $m-2$ and $14-m$, with a minus sign from $\delta_\nabla=-\iota\nabla$. This gives (AA4.1); the ordinary covariant product rule gives (AA4.2). No parallelness of $\omega$ or the weights is assumed. ∎

## AA5. Injectivity of the derivative symbol

Freeze the smooth conic weights $r_m,s_m$, not both zero, and consider the principal derivative symbol in $\omega$ alone. For a real covector $\xi$, replacing $\nabla\omega$ by $\xi\otimes\omega$ in (AA4.1) gives
$$\begin{aligned}
\|\sigma_\xi\omega\|^2={}&4(13-m)[m^2r_m^2+(12-m)^2s_m^2]\,|\xi\wedge\omega|^2\\
&+4(m-1)[(m-2)^2r_m^2+(14-m)^2s_m^2]\,|\iota_{\xi^\sharp}\omega|^2.
\end{aligned}\tag{AA5.1}$$
Both coefficients are strictly positive for $m=3,5$. Since the two form norms sum to $|\xi|^2|\omega|^2$, the symbol is injective for every real $\xi\ne0$.

*Proof.* With the displayed compact phases, or equivalently the ordinary positive Clifford coefficient norm, the squared norms of the contraction and exterior lifts on real $j$-forms are $j$ and $14-j$ times the form norm. For $m=3$ all four output grades in (AA4.1) are distinct. For $m=5$ the two exterior lifts both have grade seven. On that real compact space, right multiplication by $\Gamma$ is orthogonal and squares to $-1$, hence is skew-adjoint. Its cross pairing with the same field is zero. All remaining cross terms have distinct grades. Applying these facts to (AA4.1) proves (AA5.1). Derivatives of the weights multiply $\omega$ without differentiating it and do not enter this symbol. ∎

This is a statement about the derivative operator in the declared kernel parametrization with the conic field fixed. It does not establish ellipticity of the coupled full system, solve $Q(k)$, or supply a solution of its odd tensor equations.

## Evidence and remaining equations

The exact script checks both vector traces, all fourteen zero Gamma symbols, all four compact exterior pairs, their symbolic kernels, both kernel traces, cross-trace cancellation, and both currents on each pair. It then verifies the norm elimination, T2 identity and rational branch-sign premises symbolically. Its output prints $f$, $\delta$ and $\mathscr R$ for every conic. All 23 assertions passed in about two seconds.

The supporting derivative script checks 112 complete coordinate-jet outputs, mixed covariant jets, the variable-weight product rules, the compact phases and the positive symbol norm. A nonzero Cartan jet is annihilated on both declared copies and has nonzero output on both excluded $w/c\Gamma$ copies, for each degree. All 14 assertions passed in 0.2 seconds. The hand grade/counting proof supplies the universal statement; these exact controls check its conventions and coefficients.

The untraced odd equation is still
$$E_v+{*}^{-1}\circledast(d_0k)+2\mathcal B(v,k)+Q(k)=0.$$
In particular $d_0k$ is retained in its other odd components. The scalar reduction does not classify exterior tensor fields or solve their differential equation. Intersections, $\det M=0$, other real forms, non-invariant odd fields and a variable Shiab operator remain outside this result.
