# U. An obstruction for a declared rank-one-Ricci metric family

**Later extension.** [Z](Z-cubic-rank-one.md) removes the higher-jet restrictions at a nonzero rank-one-Ricci, zero-Weyl point. This note retains the shorter certificate and proof for its declared family.

At $\lambda=1/5$, S15 leaves the twisted cone open away from the parallel set when the base Ricci tensor has rank at most one. This note excludes a family of complete local metric polynomials, including both signs and every nonzero magnitude of the rank-one Ricci eigenvalue. It does not classify arbitrary metrics with rank-one Ricci curvature.

The proof combines two exact augmented minors with coordinate covariance. Code: `computation/geometric-unity/U_rank_one_frontier.py`; certificates and checks: `U_rank_one_frontier.out`; record: `metrics/geometric-unity/U-rank-one-frontier.yaml`. The computation imports only the earlier referee engine `xcheck_round10.py`, not the S or N computation. The result is computer-assisted. Internal review: `reports/reviews/REVIEW14-GU-U-2026-11-02.md`.

## The family and the necessary linear system

Work near $0\in\mathbb R^4$ with Euclidean metric $\delta$. For a unit vector $p$ put
$$P_p=\tfrac12(p\otimes p-\delta/6),\qquad R_p=P_p\owedge\delta,$$
where $\owedge$ is the Kulkarni–Nomizu product with the convention of S and the referee curvature engine. For $\rho\in\mathbb R$ define the complete local metric
$$g_{\rho,p;ij}(x)=\delta_{ij}-\frac{\rho}{3}(R_p)_{iajb}x^ax^b.\tag{U1}$$
It is positive definite on a sufficiently small ball. At $0$, its Ricci tensor is $\rho p\otimes p$ and its Weyl tensor is zero. Formula (U1), not just these two curvature conditions, is the hypothesis: it fixes every connection jet used below. In particular its cubic and quartic metric coefficients vanish, while covariant higher curvature derivatives need not vanish.

Let $h_{\rho,p}(c_H,c_V)$ be the Zorro metric Z1 built from the Levi-Civita connection of (U1), with
$$\lambda=\tfrac15,\qquad c_H,c_V>0.$$
Thus its signature is $(14,0)$. Let $F_0$ be the positive-definite metric fibre over $0$.

At a point of this fibre, S's relaxed system is
$$R_h(X,Y)v-\big[\xi(X)\widehat{\mathrm{Ric}}_hY-\xi(Y)\widehat{\mathrm{Ric}}_hX\big]
-\big[h(n,X)Y-h(n,Y)X\big]
=(\nabla_X\widehat{\mathrm{Ric}}_h)Y-(\nabla_Y\widehat{\mathrm{Ric}}_h)X.\tag{U2}$$
Its unknowns are two vectors $v,n$ and one covector $\xi$, giving 42 scalar unknowns and 1,274 equations. A line equation $\nabla^2r=\mu(r)\mathrm{Ric}_h+\nu h$ with $\mu(r)\ne0$ implies (U2), as in S.0. The relaxation does not require $\xi$ to be proportional to $v^\flat$.

## Covariance removes the continuous parameters

**Lemma U1.** Consistency and both matrix ranks of (U2) are invariant under base coordinate changes and constant nonzero rescalings of the observer metric.

*Proof.* Coordinate changes act by invertible changes of tensor bases, on both equations and unknowns. For $h'=\alpha h$, the connection and mixed Riemann tensor are unchanged, while
$$\widehat{\mathrm{Ric}}_{h'}=\alpha^{-1}\widehat{\mathrm{Ric}}_h,
\qquad\nabla\widehat{\mathrm{Ric}}_{h'}=\alpha^{-1}\nabla\widehat{\mathrm{Ric}}_h.$$
The invertible change of unknowns
$$v'=\alpha^{-1}v,\qquad\xi'=\xi,\qquad n'=\alpha^{-2}n\tag{U3}$$
multiplies every side of (U2) by $\alpha^{-1}$. This proves rank and consistency invariance for the relaxed linear system itself, independently of whether a line-equation solution exists. ∎

**Lemma U2.** Every metric $h_{\rho,p}(c_H,c_V)$ with $\rho\ne0$ and $c_H,c_V>0$ is, up to a constant positive homothety and an explicitly lifted base/fibre coordinate change, equivalent to one of
$$h_{+,e_0}(2,3),\qquad h_{-,e_0}(2,3).$$

*Proof.* Rotations are natural for (U1); $SO(4)$ is transitive on unit vectors, so take $p=e_0$. Under $F_s(x)=sx$, $s>0$,
$$F_s^*g_{\rho,e_0}=s^2g_{\rho s^2,e_0}.$$
Constant rescaling of the base metric preserves its Levi-Civita connection. Consequently the natural bundle lift
$$\Phi_s(x,y)=(sx,s^{-2}y)$$
satisfies $\Phi_s^*h_{\rho,e_0}(c_H,c_V)=h_{\rho s^2,e_0}(c_H,c_V)$. The direction of the fibre factor is important: $y$ is a covariant tensor, transported by the inverse base derivative.

For $\phi_q(x,y)=(x,qy)$, $q>0$, the defining formula $h=c_Hy(dx,dx)+c_VG_\lambda|_y(Dy,Dy)$ gives
$$\phi_q^*h(c_H,c_V)=h(qc_H,c_V),\qquad \alpha h(c_H,c_V)=h(\alpha c_H,\alpha c_V).$$
Choose
$$s=|\rho|^{-1/2},\qquad\alpha=3/c_V,\qquad q=\frac{2c_V}{3c_H}.$$
Then
$$\Phi:=\Phi_s\circ\phi_q,\qquad
\Phi^*\big(\alpha h_{\rho,e_0}(c_H,c_V)\big)=h_{\operatorname{sgn}\rho,e_0}(2,3).\tag{U4}$$
This proves the parameter reduction. ∎

The code checks (U4) on two nonunit controls, $(\rho,c_H,c_V)=(4,5,7)$ and $(-9,3/2,5/2)$. The pulled-back metric, covariant Ricci tensor, mixed Riemann tensor and covariant derivative of Ricci agree exactly, and both relaxed rank pairs are $(42,43)$. These checks confirm conventions; the covariance proof establishes the statement for all parameters.

## Exact certificates and their consequence

**Proposition U3 (two signed certificates).** At $y_*=\operatorname{diag}(2,3,4,5)$, both canonical metrics in U2 have relaxed coefficient rank 42 and augmented rank 43.

*Certificate.* `U_rank_one_frontier.out` records the same set of 43 zero-based row indices for each sign, together with the full nonzero rational determinant of that augmented minor. Rows are indexed by $(a,b,d)$ in lexicographic order, $0\le a<b<14$, $0\le d<14$; columns are $(v^i,\xi_i,n^i,\mathrm{rhs})$. The computation assembles these entries directly from the earlier referee metric and curvature engine, selects independent rows by exact rational elimination, and evaluates the resulting square determinant over $\mathbb Q$. A nonzero 43-minor forces augmented rank 43; deleting its last column leaves coefficient rank at least 42 and at most 42. ∎

**Theorem U4 (the declared family).** For every $\rho\ne0$, every unit $p$, and every $c_H,c_V>0$, (U2) is inconsistent at a point of $F_0$ and is consistent on no nonempty open subset of $F_0$. In particular the line equation with $\mu(r)\ne0$ holds on no open subset of $Y$ whose projection contains $0$.

*Proof.* Lemmas U1–U2 transfer U3 to the entire family. If $p=Oe_0$, one witness is
$$y=\frac{2c_V|\rho|}{3c_H}\,O\operatorname{diag}(2,3,4,5)O^\top.$$
At this point an augmented 43-minor is nonzero. By S4, its entries are rational in the fibre coordinates with denominators powers of $\det y$. Consistency on a nonempty real fibre-open set would make every augmented 43-minor vanish there, hence identically on the fibre, contradicting the witness. An open line-equation domain projecting onto $0$ would supply just such a fibre-open consistent set. ∎

**Corollary U5 (stability and the twisted cone).** The full-rank inconsistent witness persists under sufficiently small perturbations of the admissible complete finite Levi-Civita connection jet at $0$. For the metrics (U1), in the hypotheses of S.0, the twisted cone admits no $\Upsilon=0$ solution on an open set meeting $F_0$. In particular, this retains S.0's invariant odd field, slaved vector even field, $\kappa\ne0$ and $\det M\ne0$ assumptions.

*Proof.* At fixed positive $y_0$, the determinant is a continuous rational function of the finite connection jet, on the open set where the metric and its denominators are nondegenerate. Nonvanishing therefore persists on a neighbourhood of the admissible jet; no quantitative or parameter-uniform radius is asserted.

For a fixed metric (U1) with $\rho\ne0$, trivialize the metric bundle over a coordinate ball. Choose the witness $y_0$ of U4 and one fixed set of 43 rows giving a nonzero augmented minor at $(0,y_0)$. Keep that coordinate matrix and those rows fixed. Smoothness permits a smaller base neighbourhood $V$ on which the same minor is nonzero at every $(x,y_0)$; $y_0$ remains positive definite. Shrink $V$ further so $R_g\ne0$ everywhere, possible since $R_g(0)\ne0$. U4's rational-minor argument now excludes a fibre-open consistent patch in **every** fibre over $V$. The witness need not lie in a candidate solution's domain.

Suppose a twisted-cone solution has a nonempty open domain over $V$. On any component of $\mathrm{int}\,Z$, where $Z=\{dt=du=0\}$, the slaved even field vanishes and the odd equation makes Z1 Einstein (Q5(b)). S5 would force a flat base on its open projection, contradicting $R_g\ne0$ on $V$. Hence $\mathrm{int}\,Z=\varnothing$, so $U^*=\{dt\ne0\}$ is dense and nonempty, and $u=0$ by Q5(a) and continuity. These parts of Q5's argument do not require a point on $s(X)$.

Choose a patch in $U^*$. The rational coefficient $H(t)$ of Q5(c) is not identically zero; since $dt\ne0$, $t$ ranges over an interval and a smaller nonempty patch has $H(t)\ne0$. Its denominators are regular under $\det M\ne0$. The local reparametrization can be chosen with $F'\ne0$, so $\mu=xF'/H\ne0$. Q5's reduced equation then implies (U2), contradicting the fibre-open obstruction over $V$. Finally, any solution domain meeting $F_0$ has a nonempty open intersection with $\pi^{-1}(V)$, to which this argument applies. ∎

The exception $\rho=0$ is the flat base. At $\lambda=1/5$, Z1 is Einstein there and (U2) has its zero solution; Q15 describes the separate algebraic existence conditions for the constant fields. No exclusion at $\rho=0$ is claimed.

## Scope and next question

The family theorem is stronger than S15's single specified jet but remains bounded. It uses the full normal quadratic metric (U1), not the assertion $\operatorname{rank}\mathrm{Ric}=1$ alone. The Weyl tensor, cubic metric coefficients and quartic metric coefficients were not free parameters. The stability corollary gives an open neighbourhood of each certified jet, without an explicit uniform bound.

The higher-jet restrictions of this proof are removed by [W](W-quartic-rank-one.md) and [Z](Z-cubic-rank-one.md). Later certificates cover the nonzero rank-one-Ricci Weyl families of [CC](CC-diagonal-weyl.md), [EE](EE-zero-electric-weyl.md) and [FF](FF-joint-weyl.md), and every nonzero-curvature Ricci-zero point in [HH](HH-ricci-zero.md). Those results retain their stated field-sector assumptions. The [synthesis](../../reports/SYNTHESIS.md#geometric-unity) lists the remaining joint-Weyl and zero-curvature cases.
