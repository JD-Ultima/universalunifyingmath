# Internal review of the exterior derivative identities

November 2, 2026. Review of AA4–AA5 in `formalization/geometric-unity/AA-conic-traces.md` and `computation/geometric-unity/AA_exterior_derivative.py` with its output. Both results are valid under their stated assumptions. Review20 continues to cover AA1–AA3 only. This is internal review, not external peer review.

## Copies, phases and connection

The identities concern the two $\mathsf c/\mathsf w\Gamma$ conics, in degrees $m=3,5$. They do not extend the derivative reduction to the other two $\mathsf w/\mathsf c\Gamma$ conics. The helper supplies nonzero Cartan-jet outputs on those excluded copies.

In compact $(14,0)$, $\Gamma^2=-1$. The two input grades are $(2,10)$ for $m=3$, both anti-Hermitian, and $(4,8)$ for $m=5$, both Hermitian before multiplying by $i$. Thus $\varepsilon_3=1$ and $\varepsilon_5=i$ are correct. The real-member derivative has the displayed additional factor $i$, and its outputs remain anti-Hermitian. These phases are required when using the positive normalized $-\operatorname{Tr}$ pairing.

The Clifford lifts and volume element are parallel for the fixed Levi–Civita spin connection. A formula proved for arbitrary jets in a normal orthonormal frame therefore holds covariantly. It does not assume that the exterior form, the conic weights or a kernel basis is parallel.

## Independent coefficient count

I checked the full derivative on arbitrary coordinate jets, rather than inferring it from the finite controls. Let only $\nabla_j\omega=e^I$ be nonzero, with $|I|=m$. Output component $a=j$ vanishes. For $a\ne j$, the contraction copy gives, after removing the common phase,

$$\sum_{c\notin\{a,j\}}
\{\gamma_a\gamma_j\gamma_c,\iota_c e^I\}.$$

For the partner copy, moving $\Gamma$ to the right replaces the anticommutator by a commutator with $e^c\wedge e^I$. The complete output counts are

| membership of $a,j$ in $I$ | contraction copy | exterior-partner copy |
|---|---:|---:|
| neither | $m$ terms | $14-m-2$ terms |
| both | $m-2$ terms | $14-m$ terms |
| exactly one | $0$ | $0$ |

Each contributing term is twice the same Clifford product $\gamma_a\gamma_j e^I$. In the first case this is the exterior output; in the second it is the double contraction. The convention $\delta_\nabla=-\sum_j\iota_j\nabla_j$ supplies the negative codifferential coefficients. This yields exactly

$$\mathcal D(\varepsilon_m\mathsf c(\alpha))
=i\varepsilon_m[2m\,\mathsf w(d\alpha)-2(m-2)\mathsf c(\delta_\nabla\alpha)],$$

$$\mathcal D(\varepsilon_m\mathsf w(\beta)\Gamma)
=i\varepsilon_m[2(12-m)\mathsf w(d\beta)-2(14-m)\mathsf c(\delta_\nabla\beta)]\Gamma.$$

Linearity covers every covariant jet. In particular the Cartan summand, where both exterior derivative and codifferential vanish, is annihilated on these two copies.

There is a representation-theoretic subtlety at $m=5$: the grade-seven target contains two equivalent exterior copies. Absence of $\mathcal H_5$ proves factorization but does not alone select the lift of $d\alpha$. The complete decomposable-jet output fixes both possible coefficients; irreducibility of $\Lambda^6$ extends that calibration. The coordinate count above also checks the entire output directly, leaving no additional middle-grade intertwiner untested.

## Weights and conic endpoints

Both kernel charts $(r_m,s_m)=(pw,A)$ and $(B,qw)$ satisfy the block equations. At least one is nonzero at every conic point, including an endpoint, because $A,B$ cannot both vanish. In a valid chart, the kernel form is represented without a zero parametrizing vector.

The covariant product rules retain both $dr_m\wedge\omega$ and $-\iota_{\operatorname{grad}r_m}\omega$, and the corresponding $s_m$ terms. On a nonzero-$w$ arc, the weights can be written smoothly as functions of $a$. At $w=0$, $w=\pm\sqrt{f(a)}$ need not be smooth in $a$; a smooth conic chart and the general differentials are required. Even when $da=0$ at an endpoint, a weight derivative involving $dw$ need not vanish. The derivative identities themselves involve no division by $w$.

## Positive symbol norm

With the displayed compact phases, the lift norms are the ordinary positive coefficient norms: contraction of a $j$-form multiplies squared norm by $j$, and the exterior lift by $14-j$. For $m=3$, the four derivative output grades are $1,5,9,13$, so the cross terms vanish by grade orthogonality.

For $m=5$, the output grades are $3,7,7,11$. The two middle terms in the frozen-weight symbol are the same field $X=\mathsf w(\xi\wedge\omega)$ and $X\Gamma$. On this real compact grade-seven space, right multiplication by $\Gamma$ is orthogonal, squares to $-1$, and is therefore skew-adjoint. Hence $\langle X,X\Gamma\rangle=0$. The cancellation uses the same field in the two terms; it is specific to the frozen-weight symbol.

The resulting two coefficients are

$$C_w=4(13-m)[m^2r_m^2+(12-m)^2s_m^2],$$
$$C_c=4(m-1)[(m-2)^2r_m^2+(14-m)^2s_m^2].$$

Both are positive for a valid nonzero chart. Together with
$|\xi\wedge\omega|^2+|\iota_{\xi^\sharp}\omega|^2=|\xi|^2|\omega|^2$, this proves injectivity for real nonzero $\xi$.

This symbol freezes the conic field and differentiates only $\omega$. Weight derivatives are lower-order terms for that operator. They cannot be discarded when studying the coupled unknowns. No ellipticity, existence or sufficiency statement for the full nonlinear system follows.

## Reproduction and remaining equation

`python -X utf8 -I computation/geometric-unity/AA_exterior_derivative.py` passed all **14 checks** in 0.2 seconds with Python 3.12 and SymPy 1.14.0. Deterministic output matched the recorded file. The checks include 112 representative coordinate-jet full outputs, mixed jets, Cartan and excluded-copy controls, weighted chain rules, phases and the symbol norm. The helper uses the accepted Clifford engine through T; it is not a second Clifford implementation. The coordinate counting and norm argument above were checked independently.

The full odd equation still contains $\mathcal D(k)$, $2\mathcal B(v,k)$ and $Q(k)$. AA4 gives the first term explicitly on the two declared copies; it does not set it to zero. AA1–AA3 remain necessary scalar conditions, and neither they nor AA4–AA5 solve the remaining tensor equations. Intersections, $\det M=0$, other real forms and variable Shiab operators remain outside these results.
