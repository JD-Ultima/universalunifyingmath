# AS. The opposite moving Lambda5 system

The opposite conic has an exact first-order formulation. Within its ansatz, the full bosonic equation is equivalent to the grade-one metric equation, the entire Gamma equation, a determined five-form jet and an algebraic compatibility tensor. This retains the actual form Gram; [AQ's eliminated metric equation](AQ-opposite-lambda5-conservation.md) alone does not give a converse. [AW](AW-lambda5-conic-classification.md) subsequently proves that no smooth actual domain satisfies this full system under the stated positive-Z1 hypotheses.

The [helper](../../computation/geometric-unity/AS_opposite_lambda5_system.py) and its [output](../../computation/geometric-unity/AS_opposite_lambda5_system.out) use the accepted T/AL Clifford operations. [Independent internal Review49](../../reports/reviews/REVIEW49-GU-AS-2026-11-03.md) accepts AS1–AS3 and independently reproduces all 52 controls. External peer review remains outstanding.

## Setting

Retain AQ's smooth positive four-base, reconstructed Z1, constant $c_H,c_V>0$, $\lambda<1/4$, observer signature $(14,0)$, compact $\mathfrak u(128)$, true dualized equation, fixed parallel real $W_A$, real fields, constant $\kappa\ne0$ and $\det M\ne0$. Stay in one $\Lambda^5:\mathsf w,\mathsf c\Gamma$ kernel away from coincident non-vector factors, on a smooth nonzero-$w$ arc:
$$
T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\quad
k=r\mathsf w(\omega)+s\mathsf c(\omega)\Gamma,\quad T_e=v+k,
$$
$$
w^2=\frac{(1+\frac{704}{3}a)(1+160a)}{1280},\quad
(1+\tfrac{704}{3}a)r=\tfrac{80}{3}ws,\quad(1+160a)s=48wr.
$$
The real smooth weights $r,s$ are nonzero. With all derivatives retained,
$$
v=P_2(C)+P_{14}(B),\quad B=208wC,\quad
C=\frac{24}{\delta}da,\quad \delta=1+224a-3328w^2.
$$
All norms and derivatives use the observer metric $g$. Put $J=\nabla\omega$, $n=|\omega|^2$, $S_{ab}=\langle\iota_a\omega,\iota_b\omega\rangle$, $c=|C|^2$ and $h=\operatorname{div}C$.

## AS1. Complete derivative, mixed and quadratic channels

Use [AL's](AL-opposite-lambda5-system.md) full jet operator
$$
K=I-\mathsf c\operatorname{Alt}-\mathsf w\operatorname{Ctr},\qquad
K^{-1}=I-\tfrac15\mathsf c\operatorname{Alt}-\tfrac19\mathsf w\operatorname{Ctr}.
$$
Its eigenvalues are $-5,-9,1$; Cartan jets survive. Right multiplication by $\Gamma$, with $\Gamma^2=-1$, gives
$$
D(k)_5/i=-16K(rJ+dr\otimes\omega),\qquad
D(k)_9\Gamma/i=8K(sJ+ds\otimes\omega).
$$
Define, for a covector $D$,
$$
\mathcal X(D)_a=D_a\Omega,\quad
\mathcal Y(D)_a=\iota_a(D\wedge\omega),\quad
\mathcal Z(D)_a=e^a\wedge\iota_D\omega,
$$
so $K\mathcal X=\mathcal X-\mathcal Y-\mathcal Z$. The entire $M=2\mathcal B(v,k)$ has grades five/nine:
$$
\begin{aligned}
M_5/i&=r[304\mathcal X(C)-\tfrac{208}{3}\mathcal Y(C)+32\mathcal Z(C)]
-\tfrac{16}{3}sK\mathcal X(B),\\
M_9\Gamma/i&=s[-160\mathcal X(C)+\tfrac{80}{3}\mathcal Y(C)]
-\tfrac{32}{3}rK\mathcal X(B).
\end{aligned}\tag{AS1a}
$$
For $\Omega$ the Clifford image of $\omega$, set
$$
F_j=P_j(\Omega^2),\quad Y_{j,a}=P_j(\Omega\gamma_a\Omega),\quad
E_{j,a}=P_j(\gamma_a\Omega^2).
$$
The complete higher quadratic identities are
$$
\begin{aligned}
Q_5/i={}&(52r^2-\tfrac{44}{3}s^2)Y_5+(48r^2+\tfrac83s^2)E_5
+rs(56Y_9+40E_9)\Gamma,\\
Q_9/i={}&(\tfrac{152}{3}r^2-\tfrac{40}{3}s^2)Y_9
+(\tfrac{16}{3}r^2-\tfrac{56}{3}s^2)E_9
+rs(\tfrac{160}{3}Y_5+\tfrac{56}{3}E_5)\Gamma .
\end{aligned}\tag{AS1b}
$$
Only grades $1,5,9,13$ occur in $Q(k)$. Its **entire**, unsymmetrized converted grade-thirteen tensor is
$$
[Q(k)]_{13}\Gamma/i=\tfrac{160}{3}rsn\,g-\tfrac{304}{3}rsS.\tag{AS1c}
$$

*Complete basis proof.* The helper checks every output slot on $\omega=Ae^{01234}+Be^I$ for all six blade-intersection sizes, with independent symbolic $A,B,r,s$. This includes every coefficient of (AS1b) and all 196 tensor entries of (AS1c), without symmetric projection. At most ten axes are occupied, allowing an unused-axis determinant repair to $SO(14)$; equivariance and polarization cover all quadratic coefficients. The linear derivative and bilinear mixed identities have the two exhaustive coordinate types: derivative/vector axis inside or outside the five-blade. Independent vector coefficients cover arbitrary $C,B$. These basis arguments retain all weight derivatives.

## AS2. Equivalent full system and conditional uniqueness

Set $Z_5=(Q_5+M_5)/i$ and $Z_9=(Q_9+M_9)\Gamma/i$. Within the stated ansatz, the full equation is equivalent to these four tensor equations:

1. AQ1.1, with the actual $n,S$:
   $$
   0=G-24(\nabla C-hg)+224C\otimes C+(b+88c+qn)g+\ell S,
   $$
   where $b=\kappa^2(a+312a^2+104w^2)$,
   $q=112r^2/3+40s^2$ and $\ell=320r^2/3-32s^2$.
2. The entire AH2 equation:
   $$
   [-\kappa^2w(1+208a)+3328wc+\tfrac{160}{3}rsn]g
   =3328wC\otimes C+\tfrac{304}{3}rsS.
   $$
3. The determined full jet:
   $$
   J=\frac1{16r}K^{-1}(Z_5)-\frac{dr}{r}\otimes\omega.
   $$
4. Its algebraic compatibility:
   $$
   sZ_5+2rZ_9+16K((r\,ds-s\,dr)\otimes\omega)=0.\tag{AS2}
   $$

*Proof of both directions.* J1 and the R2 conic kernel solve the entire even equation. In the odd equation, J2 and AQ give the complete grade-one expression above. Since $C,B$ are closed and $B=208wC$, the vector sector has no grades three/eleven. Its converted grade thirteen is exactly
$-\kappa^2w(1+208a)g+16(\langle C,B\rangle g-C\otimes B)$; $D(v)$ contributes no Gamma term. AS1c adds the entire remaining tensor, giving AH2 without a lost antisymmetric equation. All other odd terms have grades five/nine, namely
$$
-16K(rJ+dr\otimes\omega)+Z_5=0,\qquad
8K(sJ+ds\otimes\omega)+Z_9=0.
$$
Invertibility of $K$ and $r,s\ne0$ makes these equivalent to items 3–4. Thus no Clifford grade or full compatibility component is omitted. ∎

For fixed background and fixed smooth odd field and chart weights, item 3 is a smooth linear-plus-quadratic equation in $\omega$. Existing smooth solutions with the same point value agree on a connected common domain by pathwise ODE uniqueness. This gives neither existence nor uniqueness across different backgrounds or odd fields.

Freezing those weights, the full jet map has norm squared
$$
(256r^2+64s^2)\|KJ\|^2\ge(256r^2+64s^2)\|J\|^2.
$$
Its principal symbol is injective at every nonzero real covector. No ellipticity or well-posedness of the coupled system follows.

## AS3. Traced compatibility and its limit

Conic differentiation gives
$$
r\,ds-s\,dr=\frac{7rs}{240w^2}da
=\frac{7\delta rs}{5760w^2}C.
$$
Put $\mathcal X=\mathcal X(C)$ and similarly for $\mathcal Y,\mathcal Z$, and define
$$
\vartheta=\frac{3328}{3}w\frac{s^2+4r^2}{rs}-\frac{7\delta}{360w^2}.
$$
Then (AS2) is
$$
\mathcal Q+rs[-16(\mathcal X+\mathcal Y-2\mathcal Z)-\vartheta K\mathcal X]=0,
$$
$$
\begin{aligned}
\mathcal Q={}&-\tfrac43s(41r^2+11s^2)Y_5+\tfrac83s(4r^2+s^2)E_5\\
&+\tfrac83r(38r^2+11s^2)Y_9\Gamma+\tfrac83r(4r^2+s^2)E_9\Gamma.
\end{aligned}
$$
Its contraction/exterior traces give
$$
F_4=\alpha\,\iota_C\omega,\quad F_8\Gamma=\zeta\,C\wedge\omega,\qquad
\alpha=-\frac{3r(304+9\vartheta)}{16(61r^2+16s^2)},\quad
\zeta=\frac{3s(-112+5\vartheta)}{16(64r^2+19s^2)}.
$$
The denominators are positive. These imply
$$
(\Omega-\tfrac12\alpha C-\tfrac12\zeta C\Gamma)^2
=n+\tfrac14(\alpha^2+\zeta^2)c,
$$
using $\{C,\Omega\}=2\iota_C\omega$ and
$\{C\Gamma,\Omega\}=-2(C\wedge\omega)\Gamma$.
This shifted inhomogeneous element need not be a five-blade. The full compatibility tensor, including its hook part, remains essential; its traces establish no simplicity or existence claim.

Run from the repository root:

    python -X utf8 -I computation/geometric-unity/AS_opposite_lambda5_system.py

The exact controls share T's Clifford engine, not a second metric evaluator. [AR](AR-opposite-lambda5-geometry.md) excludes zero-$w$ contacts, so the divisions here are used only where $w\ne0$. [AW](AW-lambda5-conic-classification.md) closes the actual field problem by pairing the full normal compatibility and using conservation. [AV](AV-opposite-gram-transport.md)'s question about cubic implication from zero-order data remains separate. Coincident kernels, $\det M=0$, other operators, non-invariant odd fields and indefinite real forms are outside this formulation.