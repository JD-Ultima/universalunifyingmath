# BM. The certified mixed candidate fails grade nine

The real mixed value proved to exist in [BL](BL-mixed-grade-one-existence.md)
satisfies the complete grade-one equation after amplitude normalization.
It cannot satisfy the full pointwise field equation: one grade-nine
coefficient is nonzero throughout BL's certified rational box. This
excludes that candidate, not every mixed value or every scalarity zero
in the six-amplitude family.

[BN](BN-h6-invariant-grade-nine.md) subsequently excludes the entire
135-dimensional BG module by an invariant grade-nine projection.

Evidence: the [exact reconstruction and coefficient helper](../../computation/geometric-unity/BM_mixed_grade9.py)
and its [output](../../computation/geometric-unity/BM_mixed_grade9.out).
[Review77](../../reports/reviews/REVIEW77-GU-BM-2026-11-07.md) records the
internal review, fresh reproduction and shared dependencies.

## Setting

Retain BL's flat positive Z1 model, $\lambda=1/5$, compact $(14,0)/
\mathfrak u(128)$ real form, fixed parallel real $W_A$, constant
$\kappa\ne0$ and
$$T_o=i(3\kappa/104)\Phi^1,\qquad T_e=iB\Gamma.$$
Initially $c_H=2,c_V=1$. Let $H_0,H_1,P_0,\ldots,P_3$ be BL's complete
coordinate-permutation invariant polynomial bases and define
$$B(a)=a_0\iota_3(H_0)+a_1\iota_3(H_1)
+\sum_{j=0}^3a_{j+2}\iota_5(P_j).$$
Intertwiners, tensor norms, ordered contractions and the physical phase
are exactly those of BL. No affine invariance of a field is imposed.

## BM1. Reconstructing the six actual H6 tensors

The grade-one contraction formula alone does not determine a higher-grade
coefficient. The helper therefore reconstructs all six tensors from the
public zonal seeds, rather than treating polynomial coefficients as
Clifford coefficients.

Set $c=4/5$, $s=3/5$. Successive rotations in the $(i,4)$ planes give the
rational unit directions
$$\begin{aligned}
n_0&=(0,0,0,0,1),\\
n_1&=(3/5,0,0,0,4/5),\\
n_2&=(3/5,12/25,0,0,16/25),\\
n_3&=(3/5,12/25,48/125,0,64/125),\\
n_4&=(3/5,12/25,48/125,192/625,256/625).
\end{aligned}$$
For a harmonic-degree-$d$ module, the infinitesimal rotation $D$ has
frequencies $0,\pm1,\ldots,\pm d$. On that module, $\exp(\theta D)$ is
therefore the real polynomial $p_d(D)$ determined by
$$p_d(ik)=(c+is)^k,\qquad -d\le k\le d.$$
Its coefficients are rational. Splitting even and odd powers reduces
their construction to the nonsingular interpolation systems
$$\sum_{m=0}^d e_m(-k^2)^m=\cos(k\theta),\quad 0\le k\le d,$$
$$\sum_{m=0}^{d-1}o_m(-k^2)^m=\sin(k\theta)/k,\quad 1\le k\le d.$$
The helper verifies all frequencies exactly. In these four planes,
$D=2A_i$, where $A_i$ is BC's induced covariant connection action. Full frequency
interpolation is necessary after successive rotations; an odd-frequency
zonal shortcut would not justify the later steps.

Coordinate transpositions act on $V=\operatorname{Sym}^2_0E$ by matrix
conjugation. The helper constructs their tangent matrix $X$ using the
actual Gram-orthogonal $\Phi_a$ and verifies $X^2=I$ and $X^TgX=g$.
On covectors, a source $e^a$ maps to $\sum_b(X^{-1})_{ab}e^b$; the
external slot and all six form slots receive this action.

There is a parity issue. The ambient H6 action corresponds to the
determinant-twisted natural $O(5)$ action on odd harmonic polynomials,
whereas conjugation on $V$ identifies $P$ and $-P$. A polynomial fixed by
ordinary coordinate permutations thus corresponds to the ambient sign
character. A reflection of $E$ has four negative off-diagonal directions
on $V$, so its action on $V$ has determinant $+1$ and preserves $\Gamma$.
The required Reynolds operator is the **signed** ambient projection
$$\mathcal R_{S_n}^{\rm sign}
=\frac1n\left(I-\sum_{j=0}^{n-2}\rho((j,n-1))\right)
\mathcal R_{S_{n-1}}^{\rm sign},\qquad n=2,\ldots,5.$$
This coset recurrence equals the 120-term signed group average. Under
the intertwiner, it matches ordinary polynomial averaging.

The polynomial calculation independently reconstructs the averaged zonals
at the directions above. In the BL bases, directions $(0,2)$ give cubic
coordinate determinant $1728/546875\ne0$; directions $(0,1,2,4)$ give
quintic determinant $71230958862336/7212162017822265625\ne0$.
The corresponding exact inverse matrices recover all six basis tensors
from the rotated and averaged zonal seeds. Their mapping is thus checked
from polynomials, not inferred from norm matching.

The reconstruction additionally verifies both weighted hook traces, the
full Casimir values 18/40, every adjacent-transposition sign identity and
the entire six-by-six norm Gram of BL. Actual grade-one coefficients for
one cubic basis tensor, one quintic basis tensor and their cross term
agree with the BI/BK/BL formulas. The cross value is
$-8623125/1408\ne0$, checking the relative intertwiner phase and
orientation that norm data alone would not fix.

## BM2. One exact grade-nine coefficient

Convert each six-form hook to the actual field coefficient $B\Gamma$ by
right Clifford multiplication. In BC's fourteen-dimensional basis, take
external receiver $a=0$ and form mask
$$1279=2^0+2^1+\cdots+2^7+2^{10}.$$
It has degree nine. The corresponding coefficient of $Q(iB(a)\Gamma)/i$
is
$$\boxed{\frac5{1672704}\,h(a),}\tag{BM2.1}$$
where
$$\begin{aligned}
h(a)={}&570343710720a_0a_1-51352963200a_0a_4-3736212480a_0a_5\\
&+35548651520a_1^2+49612780800a_1a_2+5095537920a_1a_3
+1183248000a_1a_4\\
&+11396602800a_2a_4+157894110a_2a_5+337252410a_3a_4\\
&-8742492a_3a_5+1884938121a_4^2.
\end{aligned}$$
The helper obtains all 21 amplitude coefficients by exact quadratic
polarization of the reconstructed tensors. It retains BC's actual
inverse-Gram Clifford products, the phase $iB\Gamma$ and the determinant
conversion factor. One separate nested-Clifford implementation agrees
at $(a_0,\ldots,a_5)=(2,-3,1,0,1,-1)$. That control checks the coefficient
implementation; it is not the real-box exclusion proof.

Fix $a_4=-1/100$, $a_5=1$, and let $c,r$ be BL's exact rational center
and radius. Write the restriction as $\widehat h(x)$ in the four variable
amplitudes. Exact arithmetic gives
$$\widehat h(c)=
-\frac{1661891521323587018793941654514756913923866165313977}
{31250000000000000000000000000000000000000000}<-10^7.$$
Because this is quadratic, the whole box satisfies
$$|\widehat h(c+v)-\widehat h(c)|\le
r\sum_j|\partial_j\widehat h(c)|+
\frac{r^2}{2}\sum_{j,k}|\partial_j\partial_k\widehat h|=\Delta,$$
with the exact bound
$$\Delta=\frac{11998004659139876640488812845599}
{7812500000000000000000000000000}<2.$$
Thus $\widehat h<0$ everywhere in BL's rational box. In particular,
the exact real scalarity zero certified there has nonzero grade-nine
quadratic stress. No approximation to that zero is used.

## BM3. Field consequence and scope

On this branch, curvature and the invariant baseline have grade one.
[BA](BA-hook-closure.md) proves that the entire variable hook derivative
has grade seven. Neither can cancel (BM2.1). Therefore no value in this
fixed-slice box can be the even field of a full solution on the declared
branch, even at one point.

The conclusion persists under nonzero real field rescaling and common
$SO(5)$ rotations: quadratic homogeneity and equivariance preserve
nonvanishing. In particular it survives BL's positive grade-one amplitude
normalization. The natural orthonormal identifications for BL's positive
constant observer scales preserve the zero/nonzero tensor condition.

BL's grade-one existence theorem remains valid. BM shows why it cannot
be promoted to a full-equation solution. Other scalarity zeros, even in
the same six-amplitude family, have not been excluded. No complete mixed
classification, closed-field construction, curved-base extension or
physical conclusion follows. External peer review remains outstanding.
