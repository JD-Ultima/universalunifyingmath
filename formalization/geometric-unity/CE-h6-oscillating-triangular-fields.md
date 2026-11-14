# CE. A small-frequency obstruction for oscillating triangular fields

Take the compact real pure H6 branch of
[BW](BW-h6-triangular-orbit-averages.md). Write $Y=G/K$ with
$G=SL(5,\mathbb R)$ and $K=SO(5)$. Let $S=AN$ be the positive
upper triangular subgroup. There is a neighborhood of zero in the
space of imaginary triangular characters such that no finite
measure superposition of the corresponding rotated closed fields
solves the full equation on a nonempty open set. The result permits
an arbitrary finite measure over the frequencies and rotations
jointly. The existence of closed fields at some nonzero frequencies
follows from the affine character family of
[BR](BR-h6-closed-field-lifts.md).

## CE1. Closed fields and a common positive receiver

Use the real positive frame spaces of BV and their
complexifications. A real frequency is $\zeta\in\mathbb R^5$ with
$\sum_j\zeta_j=0$. If $s\in S$, let $d(s)$ be the vector of the
logarithms of its positive diagonal entries. Then
$\chi_\zeta(s)=\exp(i\zeta\cdot d(s))$ is a unitary character.
Let $\alpha_\zeta$ be the invariant one-form whose value at the
identity is the differential of $s\mapsto \zeta\cdot d(s)$.
Let $W(\zeta)$ be the complex hook kernel of
$$
 D_S+i\alpha_\zeta\wedge,
$$
where $D_S$ is BR's full invariant exterior covariant derivative.
The field
$$
 B_{\zeta,b}(sK)=\chi_\zeta(s)\,s\cdot b,\qquad b\in W(\zeta),
$$
is analytic and closed by the product rule. Conversely, the displayed
kernel describes its linear closure. In particular
$W(0)=W_{\mathbb C}$, $W(-\zeta)=\overline{W(\zeta)}$, and the
even diagonal sign group preserves every $W(\zeta)$.

Let $L$ be the real polarization of the exact BV quadratic
receiver $F_y$, and write $H(v,w)=L(v,\overline w)$ for its
Hermitian extension. On $W(0)$ the real positivity of $F_y$
gives $H(v,v)>c_0\|v\|^2$ for some $c_0>0$ by separating real
and imaginary parts. There exist $\varepsilon,c>0$ such that
$$
 H(v,v)\ge c\|v\|^2
 \qquad(v\in W(\zeta),\ \|\zeta\|\le\varepsilon).
$$
Indeed, otherwise a sequence of normalized counterexamples with
$\zeta\to0$ has a convergent subsequence. The continuous
linear stack puts the limiting unit value in $W(0)$, where
it has strict positive margin. No constant-rank assumption is
needed.

## CE2. A radial Fourier limit

Let $\mu$ be a finite positive measure on pairs $(\zeta,k)$
with $\|\zeta\|\le\varepsilon$ and $k\in K$. Let
$b(\zeta,k)\in W(\zeta)$ be a measurable complex tensor with
$\int\|b\|\,d\mu<\infty$. This includes finite vector measures
by taking $\mu$ to be total variation. Consider the real part
of the integral of $k\cdot B_{\zeta,b(\zeta,k)}$. Because the
linear kernels respect conjugation, replacing the representing
vector measure by its average with its conjugate under
$(\zeta,k)\mapsto(-\zeta,k)$ makes the integral real.

At $gK$ use $g^{-1}k=Q(g,k)R(g,k)$ as in BW. The transported
integrand equals
$$
 Q(g,k)\cdot b(\zeta,k)\ \chi_\zeta(R(g,k)^{-1}).
$$
This formula is locally holomorphic in complex coordinates on
the observer, using the same QR branches as BW and the locally
holomorphic logarithms of the positive diagonal at real
arguments. Compactness of $k$ and of the bounded frequencies
gives uniform local bounds on the kernels and their derivatives.
The integral is therefore real analytic and linearly closed.

Choose $u\in K$ by Haar Fubini as in BW, using the variation
measure $\|b\|\mu$ on rotations. Write $\Delta_j(k)$ for the
leading $j$-minor of $u^Tk$ and set $\Delta_0=\Delta_5=1$.
For this choice all are nonzero almost everywhere. Choose
$\lambda_0>\cdots>\lambda_4$ with sum zero and use the ray
$g(t)=u\operatorname{diag}(e^{-\lambda_0t},\dots,e^{-\lambda_4t})u^T$.
The norm of the wedge of the first $j$ columns of $g(t)^{-1}k$
is the product of the first $j$ positive QR diagonals. The
leading Plücker coefficient gives
$$
 \log R_{jj}(g(t),k)=\lambda_jt+\ell_j(k)+o(1),
 \qquad
 \ell_j(k)=\log|\Delta_{j+1}(k)/\Delta_j(k)| .
$$
Here the individual column indices are zero-based. In the same
frames $Q(g(t),k)\to uD(k)$ for $D(k)\in M$, with $M$ the
even diagonal sign subgroup.

After transporting back by $g(t)$ and then $u$, the ray
value satisfies
$$
 v(t)=\int e^{-it\zeta\cdot\lambda}\,d\nu(\zeta)+o(1)
$$
in coefficient norm, where $\nu$ is the finite vector measure
on frequencies obtained by integrating
$D(k)\cdot b(\zeta,k)\exp(-i\zeta\cdot\ell(k))$ over the
rotations. The error estimate follows by dominated convergence
on $\|b\|\mu$ because the real QR actions are orthogonal and the
character phases have unit modulus. For each frequency atom,
$\nu(\{\zeta\})\in W(\zeta)$ since this is a closed linear
subspace invariant under $D(k)$.

## CE3. Cesàro identities for the joint measure

Choose $\lambda$ in the open strict ordering chamber so that
$(\zeta-\eta)\cdot\lambda\ne0$ for almost every pair of
distinct frequencies under the product of the finite variation
marginal with itself. This is possible by Fubini because for
each fixed $\zeta\ne\eta$ the equality is a proper hyperplane.
It also excludes collisions of distinct atoms. Denote the
Fourier integral in CE2 by $A(t)$.

For any bounded sesquilinear form $h$, Fubini expresses the
Cesàro average of $h(A(t),A(t))$ as a double measure integral
with factor
$$
 \frac1T\int_0^T e^{-it(\zeta-\eta)\cdot\lambda}\,dt.
$$
The factor has modulus at most one and tends to zero wherever
the two frequencies have different phases. The product measure
on the equal-frequency diagonal is supported on the countable
frequency atoms, since the contribution at a given $\zeta$ is
the singleton mass in the other coordinate. Dominated
convergence therefore yields
$$
 \lim_{T\to\infty}\frac1T\int_0^T h(A(t),A(t))\,dt
   =\sum_{\zeta}h(\nu(\{\zeta\}),\nu(\{\zeta\})).
$$
Both $A(t)$ and $v(t)$ are uniformly bounded, so the same identity
holds with $v(t)$. This argument works for diffuse measures even
when their support has many coincident projected frequencies.
Only the pair measure of those coincidences enters the average.

If the full equation holds on an open set, the real analytic
residual vanishes everywhere. After the transports, the
equivariant full equation gives $F_y(v(t))=0$ and positive
constant squared norm $R^2$. Since $v(t)$ is real,
$H(v(t),v(t))=F_y(v(t))$. Apply the identity first to $H$ and
then to the positive Hermitian metric. CE1 shows
$$
 0=\sum_\zeta H(\nu(\{\zeta\}),\nu(\{\zeta\}))
  \ge c\sum_\zeta\|\nu(\{\zeta\})\|^2=cR^2>0,
$$
a contradiction. Thus the entire declared small-frequency
family is excluded. The observer dilation and homothety
transfer the conclusion in the same group character
coordinates. Arbitrary imaginary frequencies and general
closed fields remain open. Internal checks are recorded in
[Review92](../../reports/reviews/REVIEW92-GU-CE-2026-11-14.md).
External peer review remains outstanding.
