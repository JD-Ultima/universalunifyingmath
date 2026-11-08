# BL. A real mixed solution of the grade-one equations

The cubic and quintic H6 modules admit a nonzero real mixed value with
scalar **full** grade-one stress. Its existence follows from an exact
rational contraction certificate, rather than a small numerical residual.
An overall real amplitude then satisfies the grade-one field equation.
This is a pointwise result; it constructs no full field solution.

Evidence: the [cross-term calibration](../../computation/geometric-unity/BL_mixed_stress.py),
[complete finite-family reduction](../../computation/geometric-unity/BL_symmetric_mixed.py)
and [real-zero certificate](../../computation/geometric-unity/BL_real_certificate.py),
with their respective [cross](../../computation/geometric-unity/BL_mixed_stress.out),
[reduction](../../computation/geometric-unity/BL_symmetric_mixed.out) and
[certificate](../../computation/geometric-unity/BL_real_certificate.out) outputs.
[Review76](../../reports/reviews/REVIEW76-GU-BL-2026-11-07.md) accepts BL1–BL3
and records contributions, shared engines and fresh reproduction.

## Setting and normalization

Use the flat positive Z1 model at $\lambda=1/5$, initially $c_H=2,c_V=1$.
Let $E=\mathbb R^5$, $V=\operatorname{Sym}^2_0E$ and
$\mathcal H_d=\operatorname{Sym}^d_0E$, with ordinary ordered-index
Euclidean tensor inner products. The dimensions of $V,\mathcal H_3,
\mathcal H_5$ are 14, 30 and 91. Use the real $SO(5)$ intertwiners of
[BI](BI-h6-cubic-and-mixed.md) and [BK](BK-quintic-gram-and-symmetry.md):
$$\iota_3(z_3)=v_{18},\qquad
z_3=t^3-\tfrac37t|x|^2,$$
$$\iota_5(z_5)=v_{40},\qquad
z_5=t^5-\tfrac{10}{11}t^3|x|^2+\tfrac5{33}t|x|^4,\quad t=x_4.$$
For real harmonic tensors $T\in\mathcal H_3$ and $U\in\mathcal H_5$, put
$$B=\iota_3(T)+\iota_5(U),\qquad
N=\|B\|^2=\frac{2197125}{352}\|T\|^2+
\frac{1366875}{32768}\|U\|^2.$$
The mixed norm term vanishes because the two Casimir modules are orthogonal.

Field statements retain [BA2](BA-hook-closure.md): the true dualized
equation, fixed parallel real $W_A$, signature $(14,0)$, compact
$\mathfrak u(128)$, real smooth fields, constant $\kappa\ne0$ and
$$T_o=i(3\kappa/104)\Phi^1,\qquad T_e=iB\Gamma.$$
The variable hook derivative has only grade seven. In grade one the
equation is
$$\ell g_Y+q_B=0,\qquad q_B=Q_1(iB\Gamma)/i,
\qquad \ell=\frac{15}{2c_V}+\frac{15\kappa^2}{52}>0.$$
The result below concerns this algebraic equation at a point.

## BL1. The complete mixed grade-one tensor

### Completing the pure-cubic formula

BI1 calibrated the symmetric pure-cubic stress. Its extension to the
**whole** grade-one tensor requires a short additional argument. Polarize
the two output matrices $S,R\in V$ before applying the classical
orthogonal invariant-tensor theorem recorded in BI. Canonical harmonic
projection gives ten standard-vector slots: $3+3+2+2$. Since $-I$ acts
trivially on an even tensor power in dimension five, the $SO(5)$ invariant
coefficient extends to an $O(5)$ invariant coefficient and is a sum of
metric pairings. Self-pairings vanish by trace-freeness.

If $k$ pairs join the two cubic tensors and $j$ join $S$ to $R$, valence
counting gives $k=j+1$. The possibilities $(3,2),(2,1),(1,0)$ give exactly
$$\|T\|^2\operatorname{tr}(SR),\quad
\operatorname{tr}(A_TSR),\quad
\sum_p\operatorname{tr}(T_pS)\operatorname{tr}(T_pR),\quad
\sum_p\operatorname{tr}(ST_pRT_p),$$
where $(T_p)_{ij}=T_{ijp}$ and $A_T=\sum_pT_p^2$. Every form is symmetric
in $S,R$. For the second, $A_T$ is symmetric and
$\operatorname{tr}(A_T[S,R])=0$; the last uses cyclic trace. Identical
input factors identify the other placements. Thus BI's calibrated formula
is the full pure-cubic tensor. BK1 already supplies the corresponding
full pure-quintic formula.

### The cross term

Use the polarization without a factor of one half:
$$q_{3,5}=q_{\iota_3(T)+\iota_5(U)}-q_{\iota_3(T)}-q_{\iota_5(U)}.$$
Define
$$C_{ab}=\sum_{ijk}T_{ijk}U_{ijkab},\qquad
Z_{ab,cd}=\sum_{ij}T_{ija}U_{ijbcd},$$
$$v_1=\operatorname{tr}(CSR),\qquad
v_2=Z_{ab,cd}S_{ab}R_{cd},\qquad
v_3=Z_{ab,cd}R_{ab}S_{cd}.$$
Here $C$ is symmetric and trace free, and $Z$ is symmetric in its last
three slots. Repeated indices run from zero through four.

The twelve-slot polarized pairing argument is complete. If $k$ pairs
join $T$ to $U$ and $j$ join $S$ to $R$, then $k=j+2$. Only
$(k,j)=(3,1),(2,0)$ are possible. The first gives $v_1$; the second
attaches the remaining cubic index to either receiver and gives $v_2,v_3$.
These three forms are independent. Use the real tensors associated with
$\operatorname{Re}(x_0+ix_1)^3$ and $\operatorname{Re}(x_0+ix_1)^5$,
zero-extended to $E$. Set $S=\operatorname{diag}(1,-1,0,0,0)$ and
$R=\operatorname{diag}(1,1,-2,0,0)$. The receiver pairs $(S,R)$, $(R,S)$
and $(E_{02}+E_{20},E_{02}+E_{20})$ give rows
$(8,0,8),(8,8,0),(4,0,0)$, whose determinant is $-256$. Invariance
therefore permits an antisymmetric output component.

The ordered calibration uses inputs $(z_3,z_5)$ and $(z_3,D_0z_5)$,
where $D_0=x_0\partial_4-x_4\partial_0$. Its receiver pairs are
$(\Phi_0,\Phi_0)$ for the first input, and $(\Phi_0,\Phi_{13})$,
$(\Phi_{13},\Phi_0)$ for the second. The exact pattern matrix is
$$\begin{pmatrix}
2/11&-4/77&-4/77\\
-32/55&-32/77&-384/385\\
-32/55&-384/385&-32/77
\end{pmatrix},\qquad \det=-6144/33275\ne0.$$
It gives the universal formula
$$\boxed{q_{3,5}(S,R)=-\frac{2131875}{2816}v_1+
\frac{2874375}{11264}(v_2+v_3).}\tag{BL1.1}$$
The equality of the last two coefficients is an actual operator
cancellation, not an imposed symmetric projection. Twelve additional
ordered entries and four entries from a separate nested-Clifford
implementation check the calibration and phase conventions. Hook traces,
full Casimir identities and both intertwiner norm ratios are also checked.
These controls calibrate the universal form; the pairing argument proves
its completeness.

For clarity about output channels, fully symmetrize $Z$:
$$W_{abcd}=\tfrac14(Z_{abcd}+Z_{bacd}+Z_{cabd}+Z_{dabc}).$$
Then $W_{aa cd}=C_{cd}/2$. Put
$$H_{abcd}=W_{abcd}-\tfrac1{18}
(\delta_{ab}C_{cd}+\delta_{ac}C_{bd}+\delta_{ad}C_{bc}
+\delta_{bc}C_{ad}+\delta_{bd}C_{ac}+\delta_{cd}C_{ab}).$$
This is fully symmetric and trace free, and
$v_2+v_3=2H(S,R)+(4/9)v_1$. Hence
$$q_{3,5}=\frac{2874375}{5632}H(S,R)
-\frac{453125}{704}\operatorname{tr}(CSR).$$
The cross term has only harmonic-degree-four and degree-two output. It
has zero operator trace and no Young-$(2,2)$ component.

Together with the two pure formulas, this supplies the full mixed tensor.
Their weighted trace is
$$\operatorname{tr}_Vq_B=-\frac{52}{3}N.$$
Consequently any scalar value satisfies
$$\boxed{q_B=-\frac{26}{21}Ng_Y.}\tag{BL1.2}$$
The trace identity alone does not imply scalarity.

## BL2. A complete six-amplitude family and an exact real zero

Let $m_\lambda$ be the sum over distinct coordinate monomials with exponent
partition $\lambda$. The complete coordinate-permutation invariant harmonic
cubic space has basis
$$H_0=4m_3-3m_{21},\qquad H_1=m_{111}.$$
Its three monomial orbits have Laplacian rank one. BK3 supplies the complete
four-dimensional quintic space with basis
$$\begin{aligned}
P_0&=4m_5+5m_{41}-10m_{32},& P_1&=m_{41}-2m_{221},\\
P_2&=2m_{311}-9m_{2111},& P_3&=m_{11111}.
\end{aligned}$$
Set
$$f_T=a_0H_0+a_1H_1,\qquad
f_U=a_2P_0+a_3P_1+a_4P_2+a_5P_3.$$
The symmetry condition is on these polynomials; no invariance of a global
field is assumed. Even coordinate permutations lie in $SO(5)$. The ambient
H6 action corresponds to the determinant-twisted natural $O(5)$ action on
these odd harmonic polynomials.

The reduction helper evaluates every ordered receiver pair in BG's actual
Gram-orthogonal basis $\Phi_a$, whose Gram is $g_a\delta_{ab}$. For each
of the 21 monomials $a_ia_j$, $i\le j$, it subtracts the scalar trace part
$$q_{ab}-g_a\delta_{ab}\,\frac1{14}\sum_cq_{cc}/g_c.$$
After division only by nonzero rational constants, the 196 equations give
four distinct integer coefficient rows of rank four. Their receiver labels
are $(0,0),(0,1),(0,7),(0,10)$. Thus the four quadrics are equivalent to
all scalarity equations in this entire family. This is a finite exact
coefficient calculation, not a rank inferred from numerical samples.

The integer rows are recorded explicitly in the reduction source and
output, in the lexicographic unordered-pair order $(i,j)$, $0\le i\le j<6$.
Their canonical JSON serialization with separators `(',', ':')` has SHA-256
`8d365144715aac0307b26b0051d5c59610bc6ead63a903f66d727d005788d6bb`.
The reduction reconstructs and verifies those rows; the real-zero helper
uses that same exact coefficient data without reading a timing-bearing
output file.

The positive norm in this family is
$$\begin{aligned}
N={}&\frac{2197125}{352}\left(140a_0^2+\tfrac53a_1^2\right)\\
&+\frac{1366875}{32768}
\left(380a_2^2+40a_2a_3+8a_3^2+33a_4^2+\tfrac1{120}a_5^2\right).
\end{aligned}$$
The $(a_2,a_3)$ block has determinant $2640>0$; all remaining displayed
coefficients are positive. Thus $N>0$ for every nonzero amplitude vector.

### Rational contraction certificate

Fix $a_4=-1/100$, $a_5=1$, and write the four quadrics as a polynomial
map $F:\mathbb R^4\to\mathbb R^4$ in $x=(a_0,a_1,a_2,a_3)$. Take
$$c=10^{-24}(4780996876601263938643,-17430528795097021097941,
-33861696489766133108980,289635016937541269176233),$$
$$r=10^{-10},\qquad \mathcal B=\{x:\|x-c\|_\infty\le r\}.$$
The center was selected with a numerical guide. Neither an approximate
residual nor a floating-point enclosure is a premise of the certificate.
Every subsequent quantity is rational.

Let $J=DF$, $R=J(c)^{-1}$ and $G(x)=x-RF(x)$. The helper checks
$\det J(c)\ne0$ and $RJ(c)=I$. With the constant Hessians
$D_{\ell jk}=\partial_j\partial_kF_\ell$, define
$$\eta=\|RF(c)\|_\infty,\qquad
L=\max_i\sum_{j,k}\left|\sum_\ell R_{i\ell}D_{\ell jk}\right|,
\qquad q=rL.$$
Because $F$ is quadratic,
$$DG(c+h)_{ij}=-\sum_k\left(\sum_\ell R_{i\ell}D_{\ell jk}\right)h_k.$$
Thus $\|DG\|_\infty\le q$ on the whole box. Exact rational inequalities give
$$0<\eta<10^{-23},\qquad 0<L<10^9,\qquad 0<q<1/10<1,
\qquad \eta+qr<r.$$
The output records the complete rational bounds and positive self-map
margin. For $x\in\mathcal B$,
$$\|G(x)-c\|_\infty\le\eta+q\|x-c\|_\infty<r.$$
Hence $G$ maps the closed complete box into itself and is a contraction.
Banach's theorem gives a unique fixed point in this box. Since $R$ is
invertible, it is exactly a real common zero of $F$. Uniqueness is only
on this fixed-coordinate slice and in this box, not in the whole family.

All four variable intervals exclude zero, and $a_5=1$. Both harmonic
components are therefore nonzero. The center itself is not asserted to be
a zero, nor are the exact zero coordinates asserted to be rational.

## BL3. Grade-one normalization and the unresolved equations

Let $B_*$ be the mixed value supplied by BL2, with $N_* =\|B_*\|^2>0$.
By (BL1.2),
$$B=\sqrt{\frac{21\ell}{26N_*}}\,B_*$$
satisfies the complete grade-one field equation. This is an exact
existence statement about a real point value under the declared
intertwiner and field conventions.

The natural affine conjugacy and positive homothety used in BK transfer
this pointwise conclusion to all constant $c_H,c_V>0$ at $\lambda=1/5$.
Both modules are seven-covariant hooks, so their squared norm and
metric-relative quadratic scale by the same positive factor $c_V^{-7}$.
The positive target $\ell$ uses the actual scales.

This theorem establishes grade one only. [BM](BM-mixed-grade-nine-obstruction.md)
separately proves that this certified box fails grade nine, so the
candidate cannot be a full pointwise solution. No full local or global
solution is constructed. [BJ4](BJ-h6-orbit-averages.md#bj4-a-constructive-121-dimensional-closed-family)
already gives a closed hook extension of every value in this cubic/quintic
module; such an extension need not satisfy the nonlinear equations.
The general pure-quintic Gram problem and
arbitrary mixed-family classification remain open. This result shows why
the separate pure-module exclusions cannot be extended to all mixed
values at grade one. It supplies no physical or literature-priority claim.
External peer review remains outstanding.
