# DI. Decaying and flat-periodic lifts of all H6 values

Every one of the 6,820 admissible H6 values is realized by a
global analytic linearly closed field that vanishes at infinity.
Every value also has a bounded analytic linear lift invariant
under a sufficiently widely spaced lattice of translations in
a flat, for any lattice rank from one to four. The periodic
fields constructed here decay away from their flat. A nonzero
prescribed value persists along the translation orbit.

## DI1. Gaussian triangular profiles and value interpolation

For a rational polynomial chain
\(p(z)=\sum_{j=0}^d z^jp_j\) with
\((A+zD)p(z)=0\), where \(A\) is triangular covariant
closure and \(D=e^{13}\wedge\), the profile
\(b_f(t)=\sum_j p_j f^{(j)}(t)\) solves
\(A b_f+D b_f'=0\) for any entire \(f\). Take \(f_j(t)\)
to be \(e^{-t^2}\) times the degree-\(d\) Taylor polynomial
of \(e^{t^2}t^j/j!\). Then \(f_j^{(\ell)}(0)=\delta_{j\ell}\)
for \(0\le\ell\le d\). Each \(p_j\) is therefore the base
value of one entire Gaussian profile.

The [verifier](../../computation/geometric-unity/DI_h6_axial_gaussian_lifts.py)
checks eight full chains from the axial sign blocks of
[DS](DS-h6-axial-multisign.md), including all their hooks
and closure equations. Its chosen degree-zero evaluations
have nonzero highest-weight pairing minors in all 13 types
of the 6,534-dimensional BR complement. The pairing includes
the identity and one rational compact rotation. All minors
are Gaussian rational reductions modulo 65,537, with
denominators prime to it. Nonzero minors prove the
characteristic-zero multiplicity ranks.

Select the resulting finite list of base values \(r_\alpha\)
and corresponding Gaussian fields \(E_\alpha\). Write \(U\)
for their compact orbit span, precisely the BR complement
because all \(r_\alpha\) are untwisted triangular closed
values. The operator
$$
 S_U=\sum_\alpha\int_{SO(5)}
         (kr_\alpha)\otimes(kr_\alpha)\,dk
$$
is positive definite on \(U\) by that span. The field
$$
 F_v(Y)=\sum_\alpha\int_{SO(5)}
      \langle kr_\alpha,S_U^{-1}v\rangle
      (k\cdot E_\alpha)(Y)\,dk
$$
has \(F_v(o)=v\) for \(v\in U\). Its values everywhere are
in the full admissible module by closure and BU. At other
points they need not be confined to \(U\).
Add the pure-type right inverse from
[DB](DB-h6-bounded-analytic-lifts.md) to obtain a linear
seed map on the full admissible module with exact base value.

## DI2. Global decay and analytic control

For an upper triangular representative \(s\) with
\(Y=ss^T\), its last diagonal satisfies \(s_{44}^2=Y_{44}\).
For a compact rotation \(k\) and \(a=ke_4\) in zero-based
indexing, the rotated scalar coordinate is
$$
 \chi_k(Y)=\tfrac58\log(a^TYa).
$$
The other moving-frame factor is orthogonal. Let
\(\Lambda=\lambda_{\max}(Y)\) and \(u\) be a unit eigenvector.
Outside the strip \(|a\cdot u|\le\Lambda^{-1/4}\), the
argument of the logarithm is at least \(\Lambda^{1/2}\).
The strip has normalized measure at most
\((3/2)\Lambda^{-1/4}\). A Gaussian times any fixed
polynomial decays faster than \(\Lambda^{-1/4}\) there
after substituting this lower bound on the logarithm.
Every bounded-weight compact average above therefore has
norm bounded by a constant times \(\Lambda^{-1/4}\),
uniformly on the whole observer. The trace-metric estimate
\(\log\Lambda\ge d(o,Y)/\sqrt{20}\) gives exponential
distance decay. Together with DB this proves the first claim.

The same estimates hold on a common complex radius at
every real base \(g_0K\) in the charts
\(g_0\exp(Z/2)K_{\mathbb C}\).
At \(Z=0\), factor \(g_0^{-1}k=Q_0R_0\) by real QR.
For complex small \(Z\), factor
\(\exp(-Z/2)Q_0=Q_ZR_Z\) by the holomorphic branch near
the real QR factorization. The matrix
\(Q_0^T\exp(-Z)Q_0\) stays uniformly near the identity
because \(Q_0\) is compact. Its upper triangular
Cholesky factor and its inverse are uniformly bounded and
holomorphic, as is \(Q_Z\). The resulting upper factor
is \(R_ZR_0\), so the character coordinate differs from
the real one by a uniformly bounded holomorphic shift.
Gaussian polynomial factors retain their decay in the real
coordinate under this bounded complex shift. Thus the
compact integral is analytic, with uniform holomorphic
distance decay on these charts. Closure passes under it.

## DI3. Periodization and its scope

Let \(\gamma_n\) be translations by \(Ln\),
\(n\in\mathbb Z^d\), in a fixed rank-\(d\) flat through
\(o\). Such flats exist for \(1\le d\le4\).
Sum the translates of the full seed map. Its Gaussian
part converges normally on local complex charts by DI2.
The pure \((5,5)\) part converges locally uniformly and
is analytic by the elliptic argument in
[DE](DE-h6-elliptic-periodic.md).
The sum preserves exact closure and the hooks.

Nearest flat projection is nonexpansive, so with flat
projection coordinate \(s(x)\),
\(d(x,\gamma_n o)\ge\|s(x)-Ln\|\). The corresponding
exponential lattice sum is uniformly bounded in \(s(x)\).
At the origin its nonzero translate tail tends to zero
as \(L\to\infty\), for example by
\(\|n\|_2\ge\|n\|_1/2\) and geometric product sums.
The evaluation map of the periodized seed on the full
6,820-dimensional module is therefore invertible for
large \(L\). Invert it to prescribe any base value.
Translation invariance preserves this norm along the orbit.

If \(r(x)\) is the distance from the flat, then
\(d(x,\gamma_n o)\ge(r(x)+\|s(x)-Ln\|)/2\).
Thus the constructed periodic field decays exponentially
away from its flat. A ray whose endpoint lies outside
the visual boundary of a closed flat has distance tending
to infinity from that flat. Indeed bounded distance along
a subsequence, and the corresponding flat projections,
produce a flat ray with the same endpoint by convexity.
For a regular fixed ray direction, the endpoints of that
direction type in any fixed flat form a finite Weyl orbit.
Hence a generic compactly rotated regular ray escapes every
one of finitely many specified flats.

Finite sums of the constructed periodic fields therefore
vanish along such a ray. In a mixture with a bounded analytic
[BW](BW-h6-triangular-orbit-averages.md) field and the decaying
DB fields, one may choose the same ray satisfying BW's
genericity condition. The limit of the transported mixture
is the BW limit, excluded by the BV positive receiver and
the positive constant norm of a full solution. Analyticity
extends this ansatz exclusion from any open set. Arbitrary
periodic fields and infinite mixtures over flats remain
outside that argument.

    python -X utf8 computation/geometric-unity/DI_h6_axial_gaussian_lifts.py

[Review111](../../reports/reviews/REVIEW111-GU-DS-DH-DI-2026-11-14.md)
records the independent checks. External peer review remains outstanding.
