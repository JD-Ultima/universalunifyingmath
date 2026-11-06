# BI. Cubic and mixed H6 value obstructions

Every nonzero value in the 30-dimensional harmonic-cubic H6 module fails
the necessary scalar grade-one equation. A separate positive-square
identity excludes a specified mixed cubic/quintic family with discrete
symmetry. Neither statement classifies arbitrary mixed values.

Evidence: [cubic helper](../../computation/geometric-unity/BI_h3_stress.py),
[mixed helper](../../computation/geometric-unity/BI_mixed_h6.py) and their
[cubic](../../computation/geometric-unity/BI_h3_stress.out) and
[mixed](../../computation/geometric-unity/BI_mixed_h6.out) outputs.
Internal [Review72](../../reports/reviews/REVIEW72-GU-BI-BJ-2026-11-06.md) and [Review73](../../reports/reviews/REVIEW73-GU-BI-BJ-2026-11-06.md) accept BI1–BI3 and disclose shared contributions. The complementary degree-two result is [BH](BH-h6-degree2.md).

## Setting

Use [BF/BG's](BG-h6-fixed-values.md) flat positive Z1 model at
$\lambda=1/5$, initially $c_H=2,c_V=1$. Its tangent is
$V=\operatorname{Sym}^2_0(\mathbb R^5)$ with trace metric. Let
$W_{18},W_{40}\subset H_6$ be the contained real harmonic-degree-three
and degree-five modules, of dimensions 30 and 91. All hook slots are
covariant. Their polynomial intertwiners are normalized by
$$\iota_3\left(t^3-\tfrac37t|x|^2\right)=v_{18},\qquad
\iota_5\left(t^5-\tfrac{10}{11}t^3|x|^2+\tfrac5{33}t|x|^4\right)=v_{40},
\quad t=x_4.$$

For field consequences, retain [BA2](BA-hook-closure.md): the true
dualized equation, fixed parallel real $W_A$, signature $(14,0)$,
compact $\mathfrak u(128)$, real smooth fields, constant $\kappa\ne0$ and
$$T_o=i(3\kappa/104)\Phi^1,\qquad T_e=iB\Gamma.$$
The necessary grade-one equation is
$$\ell g_Y+Q_1/i=0,\qquad \ell=\frac{15}{2c_V}+\frac{15}{52}\kappa^2>0.$$
BA's universal differentiated-hook identity leaves only derivative grade
seven, including for variable amplitudes.

## BI1. The whole harmonic-cubic module

Write a real harmonic cubic as $f_T(x)=T_{abc}x_ax_bx_c$, where $T$ is
fully symmetric and trace free. Put $(T_p)_{ab}=T_{abp}$,
$A=\sum_pT_p^2$ and $N=\operatorname{tr}A=\sum T_{abc}^2$.
For $S\in V$, the symmetric grade-one quadratic has the universal form
$$\begin{aligned}
q_T(S)={}&\alpha N\operatorname{tr}(S^2)+\beta\operatorname{tr}(AS^2)\\
&+\gamma\sum_p\operatorname{tr}(T_pS)^2
+\delta\sum_p\operatorname{tr}(ST_pST_p),
\end{aligned}\tag{BI1.1}$$
with the physical compact phase $i\iota_3(T)\Gamma$.

**Completeness of this list.** Precompose with the orthogonal trace-free
projections to regard $q$ as an invariant polynomial on two symmetric
three-tensors and two symmetric two-tensors. Its coefficient has ten
standard-vector slots. Since $-I$ acts trivially on an even tensor power
in dimension five, $SO(5)$ invariance implies $O(5)$ invariance. The
classical orthogonal invariant-tensor theorem gives metric pairings;
see Goodman, [Theorem 7.4 and Corollary 7.5, printed p. 26](https://sites.math.rutgers.edu/~goodman/pub/nuslect.pdf#page=26).
The real form follows by complexification, or the equivalent compact-real
version of that theorem.

Self-pairings vanish by trace-freeness. If $k$ pairs join the two $T$
factors and $j$ pairs join the two $S$ factors, counting slots gives
$k=j+1$. The possibilities $(k,j)=(3,2),(2,1),(1,0)$ give respectively
the first term, the second term and the two last terms in (BI1.1).
These exhaust the pairings. Thus calibration determines a universal
identity, rather than fitting selected samples.

Four exact controls, with determinant $-4096/60025\ne0$, give
$$ (\alpha,\beta,\gamma,\delta)=
\left(-\frac{9124875}{968},\frac{4095875}{484},
-\frac{268625}{484},\frac{1031625}{484}\right).\tag{BI1.2}$$
The controls use the zonal cubic and $x_0x_1x_2$. The latter is obtained
by three calibrated rotations of the zonal seed. All 28 diagonal entries,
a nonzero symmetric off-diagonal mixture and two separate nested-Clifford
targets check the conventions. The common norm ratio
$\|\iota_3(T)\|^2/N=2197125/352$ also agrees with the full trace identity.

**Rank obstruction.** Define
$$L_T:V\to\mathbb R^5,\qquad (L_TS)_p=\operatorname{tr}(T_pS),
\qquad K_T=L_T^*L_T.$$
This is positive semidefinite and has rank at most five. In indices,
$K_{ab,cd}=\sum_pT_{abp}T_{cdp}$, with pair symmetry and zero pair traces.
Equation (BI1.1) is a linear transform $F$ of this operator.

For completeness, the relevant operator decomposition is explicit. Set
$$P_{ab,cd}=\tfrac12(\delta_{ac}\delta_{bd}+\delta_{ad}\delta_{bc})
-\tfrac15\delta_{ab}\delta_{cd}.$$
For trace-free symmetric $H$, the operator $K_H(S,S)=\operatorname{tr}(HS^2)$ is
$$ (K_H)_{ab,cd}=\tfrac14(\delta_{ac}H_{bd}+\delta_{ad}H_{bc}
+\delta_{bc}H_{ad}+\delta_{bd}H_{ac})
-\tfrac15(\delta_{ab}H_{cd}+H_{ab}\delta_{cd}).$$
For an arbitrary pair-symmetric, pair-trace-free $K$, take its cross trace
$A_{ac}=\sum_bK_{ab,cb}$ and $N=\operatorname{tr}A$. Subtract
$K_0=(N/14)P$ and $K_2=K_H$, with $H=20(A-NI/5)/27$.
The remainder has every trace zero. Its full symmetrization is $K_4$;
the remainder $K_{22}$ has zero full symmetrization. This constructs a
direct decomposition $K=K_0+K_2+K_4+K_{22}$ without a representation
multiplicity assumption.

The trace-free projection of the reshuffle
$J(K)(S,S)=K_{ac,bd}S_{ab}S_{cd}$ acts on those components by
$3/10,1/10,1,-1/2$, respectively. Substituting $P,K_H$, full symmetry
and zero symmetrization proves these identities. The helper checks all
625 entries of nonzero controls. The three nonscalar multipliers of $F$
are therefore
$$\frac{27}{20}\beta+\gamma+\frac{\delta}{10}
=\frac{21455875}{1936},\qquad
\gamma+\delta=\frac{190750}{121},\qquad
\gamma-\frac\delta2=-\frac{142625}{88}.$$
All are nonzero. Scalar $q_T$ consequently forces $K_T$ scalar.
A nonzero scalar operator on the 14-dimensional $V$ has rank 14,
contradicting $\operatorname{rank}K_T\le5$. Thus $K_T=0$ and
$N=\operatorname{tr}K_T=0$ force $T=0$ over the real positive metric.

Hence no nonzero value in the entire $W_{18}$ has scalar symmetric
grade-one stress. Only the symmetric part was needed; no untested
antisymmetric component was assumed to vanish.

## BI2. A mixed family beyond the fixed-value obstruction

Set $r^2=\sum_{i=0}^4x_i^2$ and define the harmonic power sums
$$P_3=\sum_i\left(x_i^3-\tfrac37x_ir^2\right),\qquad
P_5=\sum_i\left(x_i^5-\tfrac{10}{11}x_i^3r^2+\tfrac5{33}x_ir^4\right).$$
Let $U_3=\iota_3(P_3)$ and $U_5=\iota_5(P_5)$. Their values are fixed by
the even coordinate permutations $A_5\subset SO(5)$. The helper constructs
them as sums of five quarter-rotated zonal seeds, checks both hooks and
the full ten-term Casimir, and verifies their zero infinitesimal
stabilizers. Since the two summands are distinct, every nonzero
$aU_3+bU_5$ also has zero infinitesimal stabilizer. This family is therefore
outside all conjugate-$SO(4)$ fixed-value cases of BG.

All fourteen metric-relative grade-one diagonal moments are computed by
exact polarization, including the mixed term. For $B=aU_3+bU_5$, the
difference between the displayed diagonal and off-diagonal tangent-axis
moments is
$$\begin{aligned}
q_D-q_O&=\frac{375}{495616}(892928a^2+163520ab+15975b^2)\\
&=\frac{375}{495616}\left[
892928\left(a+\frac{2555b}{27904}\right)^2+\frac{7402175}{872}b^2\right].
\end{aligned}$$
Every coefficient is positive. Scalar stress forces $a=b=0$.
The norm, complete diagonal trace and three independent nested-Clifford
targets also agree. These are necessary diagonal tests; no full
off-diagonal tensor or irreducible ten/four splitting is asserted.
Common $SO(5)$ rotations preserve the obstruction. Independently rotated
components and arbitrary cubic/quintic mixtures remain outside this family.

## BI3. Field consequence and limits

At any point of a full solution on the declared branch, a nonzero value
in $W_{18}$ fails BI1, and a nonzero value in BI2's family fails BI2.
Zero fails $\ell>0$. Thus neither kind of value can occur even at one
point. There is no assumption of constant amplitudes or a closed-field lift.

As in BG/BH, affine-conjugating dilation and positive homothety transfer
these results to all constant $c_H,c_V>0$ at this fixed $\lambda$.
Metric-relative quadratic moments share the positive factor $c_V^{-7}$;
scalarity and the positive baseline are preserved.

Arbitrary mixed values in the contained 135-dimensional curvature module,
the whole degree-five summand, other ambient curvature components,
curved bases, other parameters and coincident kernels remain open.
[BJ](BJ-h6-orbit-averages.md) addresses a separately specified global
closure family. These are conditional mathematical results, with no
physical or literature-priority claim. External peer review is outstanding.
