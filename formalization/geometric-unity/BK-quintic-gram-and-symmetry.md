# BK. The quintic Gram condition and two symmetry-class exclusions

On the pure 91-dimensional quintic H6 module, the complete grade-one equation
reduces to an explicit real Gram identity and a norm normalization. The
general Gram problem remains open. Two full permutation-invariant subspaces
are excluded exactly: a four-dimensional coordinate-permutation space and
the one-dimensional standard-simplex space.

Evidence: [Gram helper](../../computation/geometric-unity/BK_quintic_gram.py),
[symmetry helper](../../computation/geometric-unity/BK_symmetric_quintics.py)
and their [Gram](../../computation/geometric-unity/BK_quintic_gram.out) and
[symmetry](../../computation/geometric-unity/BK_symmetric_quintics.out) outputs.
Internal [Review74](../../reports/reviews/REVIEW74-GU-BK-2026-11-06.md) and [Review75](../../reports/reviews/REVIEW75-GU-BK-2026-11-06.md) accept BK1–BK4 and disclose shared contributions and implementations.

## Setting

Use the flat positive Z1 model at $\lambda=1/5$, initially $c_H=2,c_V=1$.
Put $E=\mathbb R^5$, $V=\operatorname{Sym}^2_0E$ and
$W=\operatorname{Sym}^3_0E$, of dimensions 14 and 30. Tensor inner products
are ordered-index Euclidean contractions. [BF](BF-h6-curvature-module.md)
identifies the real quintic module $W_{40}\subset H_6$ with harmonic quintics
$$f_T(x)=T_{ijklm}x_ix_jx_kx_lx_m,\qquad T=T_{(ijklm)},\quad T_{aaijk}=0.$$
Normalize its intertwiner by
$$\iota_5\left(t^5-\tfrac{10}{11}t^3|x|^2+\tfrac5{33}t|x|^4\right)=v_{40},
\qquad t=x_4.$$

For field consequences retain [BA2](BA-hook-closure.md): the true dualized
equation, fixed parallel real $W_A$, signature $(14,0)$, compact
$\mathfrak u(128)$, real smooth fields, constant $\kappa\ne0$ and
$$T_o=i(3\kappa/104)\Phi^1,\qquad T_e=iB\Gamma,\quad B\in W_{40}.$$
The variable hook derivative has only grade seven. Hence grade one requires
$$\ell g_Y+Q_1/i=0,\qquad \ell=\frac{15}{2c_V}+\frac{15}{52}\kappa^2>0.$$

## BK1. The complete grade-one criterion

For each ordered triple $J=(a,b,c)$ let $(T_J)_{ij}=T_{ijabc}$, and set
$$A_T=\sum_JT_J^2,\qquad N=\operatorname{tr}A_T=\|T\|^2.$$
For $S,R\in V$, the full grade-one tensor of $i\iota_5(T)\Gamma$, divided by
$i$, is the symmetric bilinear form
$$\begin{aligned}
q_T(S,R)={}&\alpha N\operatorname{tr}(SR)+\beta\operatorname{tr}(A_TSR)\\
&+\gamma\sum_J\operatorname{tr}(T_JS)\operatorname{tr}(T_JR)
+\delta\sum_J\operatorname{tr}(ST_JRT_J),
\end{aligned}\tag{BK1.1}$$
where
$$ (\alpha,\beta,\gamma,\delta)=
\left(-\frac{2311875}{32768},\frac{2244375}{32768},
\frac{1873125}{32768},\frac{860625}{16384}\right).$$

**Why there is no missing antisymmetric part.** Polarize the two output
arguments before setting them equal. As in [BI1](BI-h6-cubic-and-mixed.md),
the orthogonal invariant-tensor theorem reduces the fourteen standard-vector
slots to metric pairings. Trace-free projection gives a canonical extension
to the ambient tensor spaces; $-I$ acts trivially on fourteen slots, so
$SO(5)$ invariance extends to $O(5)$ here.
If $k$ pairs join the two five-tensors and $j$ join the two output matrices,
valence counting gives $k=j+3$. The possibilities $(5,2),(4,1),(3,0)$ give
the norm term, the $A_T$ term and the two remaining terms in (BK1.1).
The two input $T$ factors are identical. Each resulting bilinear form is
symmetric in $S,R$: the only possibly ordered trace satisfies
$\operatorname{tr}(A_TSR)=\operatorname{tr}(A_TRS)$ for real symmetric
matrices; the last term is symmetric by cyclic trace. Thus the complete
grade-one tensor, not only its symmetric projection, has this form.

A four-row exact calibration with determinant $131072/49413375\ne0$ fixes
the coefficients. Three short Lie-word quintics supply 42 diagonal
controls, weighted hook and full Casimir checks. Both ordered entries of
a nonzero off-diagonal control agree with the formula and a separate
nested-Clifford implementation. The physical phase and determinant factor
for $iB\Gamma$ are retained. These finite controls calibrate the universal
contraction formula; its completeness is the preceding hand argument.

Define
$$L_T:V\to W,\qquad (L_TS)_{abc}=T_{ijabc}S_{ij},\qquad K_T=L_T^*L_T.$$
Its four-index tensor and trace are
$$K_{ij,kl}=\sum_{abc}T_{ijabc}T_{klabc},\qquad \operatorname{tr}K_T=N.$$
The explicit operator decomposition and reshuffle from BI1 apply unchanged.
The three nonscalar multipliers are
$$\frac{20300625}{131072},\qquad \frac{3594375}{32768},
\qquad \frac{253125}{8192},$$
all nonzero. Therefore
$$\boxed{q_T\text{ is scalar}\quad\Longleftrightarrow\quad
L_T^*L_T=\frac{N}{14}I_V.}\tag{BK1.2}$$
This is an identity of positive real tensor Grams, not a pointwise gradient
norm condition for $f_T$.

The scalar multiplier is $-5923125/8192$, and
$\|\iota_5(T)\|^2=(1366875/32768)N$. On (BK1.2),
$$q_T=-\frac{26}{21}\|B\|^2g_Y.$$
Consequently the entire grade-one equation is equivalent to (BK1.2) and
$$\boxed{\|B\|^2=\frac{21}{26}\ell.}\tag{BK1.3}$$
If a nonzero Gram solution exists, a real amplitude can meet this positive
normalization. Higher quadratic grades and covariant closure still need
their own equations.

## BK2. An equivalent constrained-isometry problem

Define $a_i:W\to V$ by $(a_iU)_{ab}=U_{iab}$. For
$f_S(x)=x^TSx$, its adjoint is harmonic creation:
$$f_{a_i^*S}(x)=x_if_S(x)-\tfrac27|x|^2(Sx)_i.$$
Trace-free projection proves this formula; its two Laplacian terms are
$4(Sx)_i$ and $14(Sx)_i$. Ordered-index contraction gives
$\sum_i a_i^*a_i=I_W$.

A nonzero solution of (BK1.2) exists exactly when there is a real isometry
$$L:V\hookrightarrow W\quad\text{such that every }a_iL\text{ is self-adjoint}.
\tag{BK2.1}$$
To check the symmetry condition, use the canonical pair-trace-free coefficient
tensor $T_{abc,de}$ of a map $V\to W$. Self-adjointness of $a_aL$ is precisely
$T_{abc,de}=T_{ade,bc}$. Its existing $S_3\times S_2$ symmetries and the new
permutation $(2\ 4)(3\ 5)$ generate $S_5$: conjugating $(1\ 2)$ gives the
cross transposition $(1\ 4)$. Hence all five indices are symmetric; the
existing traces then give complete trace-freeness. The converse is immediate.
Rescaling a nonzero Gram solution normalizes $L^*L=I_V$ and $N=14$.

The codomain has dimension 30, so the cubic rank obstruction does not apply.
Indeed the zonal quintic already has Gram rank 14, with eigenvalues
$28/1089$, $4/3267$ and $140/1089$ of multiplicities 4, 9 and 1.
Its failure is unequal eigenvalues, not insufficient rank.
Equations (BK2.1) remain a genuine real feasibility problem in the
91-dimensional harmonic-quintic subspace of $\operatorname{Hom}(V,W)$.
Dimension counting or the existence of arbitrary isometric embeddings
does not settle that intersection.

## BK3. Two complete permutation-invariant subspaces are excluded

### Coordinate permutations in five dimensions

Let $m_\lambda$ be the sum of distinct coordinate monomials of exponent
partition $\lambda$. The seven partitions of degree five give the complete
$S_5$-invariant polynomial space. Its Laplacian has rank three. The whole
four-dimensional harmonic kernel has basis
$$\begin{aligned}
F_0&=4m_5+5m_{41}-10m_{32},&F_1&=m_{41}-2m_{221},\\
F_2&=2m_{311}-9m_{2111},&F_3&=m_{11111}.
\end{aligned}$$
For $f=\sum_{j=0}^3u_jF_j$, the complete 625-slot Gram deviation has coefficient
rank four on the ten quadratic amplitude monomials. Thus all its equations
are equivalent to the following four, not just sampled receivers:
$$\begin{aligned}
0={}&46800u_0^2+5280u_0u_1+384u_1^2-2070u_2^2-u_3^2,\\
0={}&1200u_0^2+240u_0u_1+240u_0u_2+16u_1^2+72u_1u_2+54u_2^2-3u_2u_3,\\
0={}&600u_0^2+480u_0u_1-240u_0u_2+56u_1^2+216u_1u_2+264u_2^2-9u_2u_3,\\
0={}&660u_0u_2+84u_1u_2-2u_1u_3+111u_2^2.
\end{aligned}$$
Every nonzero amplitude vector has a highest nonzero coordinate $u_j$.
Set higher coordinates to zero and rescale $u_j=1$. For these four
strata, the helper reconstructs rational polynomial identities
$\sum_i h_if_i=1$, with multiplier degrees 3, 1, 0 and 0 for $j=3,2,1,0$.
It expands each identity exactly. Thus there is no nonzero complex
polynomial Gram solution in this subspace, and in particular none over
the reals. No numerical root exclusion is used.

This also covers the whole subspace fixed by even coordinate permutations
$A_5\subset SO(5)$. An $A_5$-invariant polynomial's odd-permutation part
would be alternating and divisible by the degree-ten Vandermonde product;
it vanishes in degree five.

### The standard simplex in five dimensions

Use $E=\{y\in\mathbb R^6:\sum y_i=0\}$ with projector $\Pi=I-J/6$.
Write $p_j=\sum y_i^j$ and $r^2=p_2$. Restricting the seven degree-five
monomial orbit sums to this hyperplane gives exactly
$\operatorname{span}\{p_5,p_2p_3\}$. Their Laplacians are
$\Delta p_5=(50/3)p_3$ and $\Delta(r^2p_3)=22p_3$.
Hence the entire standard-$S_6$ invariant harmonic space is the line
$$f=p_5-\tfrac{25}{33}r^2p_3.$$
Even permutations $A_6$ have the same fixed space in this degree, because
the alternating part would be divisible by a degree-fifteen Vandermonde.

The helper constructs all ordered five-tensor entries in the six-coordinate
ambient model and verifies trace-freeness and annihilation of the all-ones
direction. Ambient contractions therefore equal those on $E$; no radical
orthonormal basis is needed. Its norm is $N=10/33$, so (BK1.2) requires
$c=N/14=5/231$. For the projector
$$P_{ab,cd}=\tfrac12(\Pi_{ac}\Pi_{bd}+\Pi_{ad}\Pi_{bc})-\tfrac15\Pi_{ab}\Pi_{cd},$$
one exact coefficient already contradicts that requirement:
$$K_{00,00}-cP_{00,00}=\frac{35}{3267}-\frac{25}{2079}
=-\frac{10}{7623}\ne0.$$
All 1,296 ambient coefficients are also compared. Scaling or orthogonal
rotation cannot repair this mismatch, so the whole declared line and its
rotates are excluded.

## BK4. Field scope and the remaining question

No nonzero value in either fixed subspace of BK3 can satisfy grade one;
zero fails $\ell>0$. Such a value cannot occur even at one point of a full
solution on the declared pure branch. Common $SO(5)$ rotations preserve
the conclusion, including conjugate copies of these finite groups.

Affine-conjugating horizontal dilation and positive homothety transfer the
criteria and exclusions to all constant $c_H,c_V>0$ at this fixed $\lambda$.
The seven-covariant hook norm and metric-relative quadratic share the
positive factor $c_V^{-7}$, preserving (BK1.2)–(BK1.3) with the actual
baseline above.

The general real problem (BK1.2), equivalently (BK2.1), is unresolved.
No nonzero general quintic candidate or universal quintic exclusion follows
from these fixed-space classifications. Arbitrary mixed harmonic values,
higher quadratic grades, closed-field reconstruction, curved bases, other
parameters and coincident kernels retain their separate scope. Classical
invariant theory is used as recorded in [BI](BI-h6-cubic-and-mixed.md);
no physical or literature-priority claim is made. External peer review is
outstanding.
