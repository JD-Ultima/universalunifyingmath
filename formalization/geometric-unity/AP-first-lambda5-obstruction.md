# AP. No first Lambda5 conic domains

The first $W_A$ Lambda5 conic admits no smooth full bosonic solution on a nonempty positive-Z1 observer domain. At a nonzero-$u$ point, its higher algebraic equations force a simple five-form, contradicting [AH2's Gram-rank requirement](AH-lambda5-gram.md#ah2-rank-and-support-at-a-nonzero-w-point). This obstruction retains the slaved vector field and all derivatives. [AH9](AH-lambda5-gram.md#ah9-the-first-pair-constant-endpoint-domains-are-empty) supplies the remaining constant endpoints.

The helper is [AP_first_lambda5_obstruction.py](../../computation/geometric-unity/AP_first_lambda5_obstruction.py), with its recorded output beside it. [Independent internal Review46](../../reports/reviews/REVIEW46-GU-AP-2026-11-03.md) accepts AP1–AP3 and independently reproduces all 40 controls. External peer review remains outstanding.

## Setting

Keep the true reconstructed equation, compact real $(14,0)$ and $\mathfrak u(128)$, fixed parallel real $W_A$, real fields, constant $\kappa\ne0$, invariant odd field and $\det M\ne0$. Work in a single first $\Lambda^5:c,w\Gamma$ kernel, away from every other non-vector singular factor:
$$
T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\qquad
k=i[r\,\mathsf c(\omega)+s\,\mathsf w(\omega)\Gamma],
$$
where $\omega$ is a smooth real five-form and $r,s$ are smooth real chart weights, not both zero. They satisfy
$$
w^2=\frac{(1+96a)(1+\frac{832}{3}a)}{1280},\qquad
(1+96a)r=48ws,\qquad
(1+\tfrac{832}{3}a)s=\tfrac{80}{3}wr.
$$
The full even field is
$$
T_e=v+k,\qquad v=P_2(C)+P_{14}(208wC),\qquad
C=\frac{24}{\delta}da,\quad \delta=1+224a-3328w^2.
$$
No constant-odd, parallel-form or curvature-trace assumption is imposed.

## AP1. Universal higher algebra

For any real five-form and nonzero real weight pair $(r,s)$, the equations $[Q(k)]_5=[Q(k)]_9=0$ hold exactly when $\omega$ is zero or simple.

*Proof.* Write $\Omega$ for its Clifford image, $P_j$ for grade projection,
$$
F_j=P_j(\Omega^2),\qquad
Y_{j,a}=P_j(\Omega\gamma_a\Omega),\qquad
E_{j,a}=P_j(\gamma_a\Omega^2).
$$
All Gamma products below are **right multiplication**, with $\Gamma^2=-1$. The full identities are
$$
\begin{aligned}
Q_{5,a}/i={}&(\tfrac{28}{3}r^2-36s^2)Y_{5,a}
+(-\tfrac83r^2+\tfrac{112}{3}s^2)E_{5,a}\\
&+rs(\tfrac{104}{3}Y_{9,a}-\tfrac{136}{3}E_{9,a})\Gamma,\\
Q_{9,a}/i={}&(\tfrac{40}{3}r^2-40s^2)Y_{9,a}
+(-24r^2-\tfrac{16}{3}s^2)E_{9,a}\\
&+rs(\tfrac{128}{3}Y_{5,a}-24E_{5,a})\Gamma.
\end{aligned}\tag{AP1a}
$$
Their exhaustive basis certificate checks every output slot on $\omega=Ae^{01234}+Be^J$ for all six intersection sizes $0,\ldots,5$, with independent symbolic real $A,B,r,s$. Any coordinate five-blade pair occupies at most ten axes. A signed permutation carries its representative to that pair; an unused-axis reflection repairs its determinant without changing either input. Thus $SO(14)$ equivariance, which preserves $\Gamma$, and polarization cover every quadratic coefficient, including diagonal terms.

Reversion gives $\Omega^2=N+F_4+F_8$, $N=|\omega|^2$, and sandwich grades $1,5,9$. Using $\sum_a\gamma_a\Omega\gamma_a=-4\Omega$ and $\sum_a\gamma_aE_{5,a}=10F_4$, $\sum_a\gamma_aE_{9,a}=6F_8$:
$$
\begin{aligned}
\sum_a\gamma_aQ_{5,a}/i
 &=(-64r^2+\tfrac{1552}{3}s^2)F_4-\tfrac{1232}{3}rsF_8\Gamma,\\
\sum_a\gamma_aQ_{9,a}/i
 &=(-\tfrac{592}{3}r^2+128s^2)F_8-\tfrac{1232}{3}rsF_4\Gamma.
\end{aligned}\tag{AP1b}
$$
These separate grades $4,6$ and $8,10$. If $rs\ne0$, both square components vanish. If one weight is zero, both ordinary trace factors remain nonzero, giving the same conclusion.

Now $E_5=E_9=0$. The equations for $(Y_5,Y_9\Gamma)$ have matrix
$$
\begin{pmatrix}
\frac{28}{3}r^2-36s^2&\frac{104}{3}rs\\
-\frac{128}{3}rs&\frac{40}{3}r^2-40s^2
\end{pmatrix},
\qquad
\det=\frac{32}{9}(35r^4+176r^2s^2+405s^4)>0.
$$
Hence the sandwich is vector-valued. For $N>0$, conjugation by $\Omega$, with $\Omega^{-1}=\Omega/N$, is an orthogonal involution of the vector space. In an orthonormal eigenbasis, a grade-five monomial commutes with a generator exactly when its support contains that axis. Thus every nonzero coefficient must have support equal to the entire $+1$ eigenspace; it has dimension five and only one five-blade can occur. This is AH9's elementary simplicity argument. Conversely a simple form has scalar square and vector sandwich, so (AP1a) vanishes. ∎

## AP2. Nonzero-$w$ points are impossible

The full odd equation retains $D(v)$, $D(k)$, $Q(v)$ and $2B(v,k)$. Nevertheless its grade-five/nine equations are precisely those of AP1.

Indeed [AA4](AA-conic-traces.md) gives
$$
D(k)=-10\mathsf w(d(r\omega))+6\mathsf c(\delta_\nabla(r\omega))
+[-14\mathsf w(d(s\omega))+18\mathsf c(\delta_\nabla(s\omega))]\Gamma,
$$
of grades $3,7,11$, including every weight derivative. The exact mixed support is also
$$
\operatorname{grades}\bigl(2B(P_2(C)+P_{14}(B_0),k)\bigr)\subset\{3,7,11\}
$$
for independent arbitrary covectors $C,B_0$. Its complete bilinear certificate has two coordinate types: the vector axis lies inside or outside the five-blade. Symbolic independent vector coefficients and kernel weights check all output slots. The same signed-permutation argument and linearity prove universality. Vector derivatives and quadratics have no grade five/nine, as checked on the complete equal/unequal coordinate-pair types. Neither does the curvature/invariant background.

At $w\ne0$, AH2 gives $N>0$ and contraction-Gram rank at least thirteen. AP1 makes $\omega$ nonzero simple, with rank five: a contradiction. No derivative has been set to zero, and no metric identity or inequality is continued off the actual domain. ∎

## AP3. The entire first conic is empty on positive Z1

Take a smooth positive four-base, reconstructed Z1, constant $c_H,c_V>0$ and $\lambda<1/4$, with the standing field assumptions. If a full solution existed, AP2 would force $w=0$ everywhere. The conic then forces $a=-1/96$ or $-3/832$, constant on each connected domain. AH9 excludes both constant endpoint domains. Therefore no nonempty full solution domain exists in this single conic.

Run from the repository root:

    python -X utf8 -I computation/geometric-unity/AP_first_lambda5_obstruction.py

Its 40 exact controls share T's accepted Clifford engine. The complete basis arguments and hand determinant/simplicity proof supply universality. This is no independent metric-engine or physical test. The opposite conic, coincident kernels, $\det M=0$, other operator members, non-invariant odd fields and indefinite real forms remain outside this theorem.
