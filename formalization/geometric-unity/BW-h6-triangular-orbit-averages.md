# BW. No finite-measure average of triangular closed H6 fields solves the full equation

The [BR](BR-h6-closed-field-lifts.md) family provides global analytic
closed fields on the flat positive observer, and finite sums of its
compact rotations realize every point value in a 6,534-dimensional
submodule. The full equations exclude every field obtained by a fixed
finite signed measure of those rotated fields, even on a nonempty open
subset. The key inputs are a diagonal flag limit and
[BV](BV-h6-triangular-positive-obstruction.md)'s exact quadratic
obstruction to the entire triangular family.

## BW1. Fields from finite vector measures

Write the observer as $Y=G/K$ with $G=SL(5,\mathbb R)$ and $K=SO(5)$,
and let $S=AN$ be the positive upper triangular determinant-one subgroup.
Let $W$ be BR's 952-dimensional value space at $o=K$, and $B_b$ the
$S$-invariant field with $B_b(o)=b$. The group $M$ of determinant-one
diagonal signs belongs to $K$, normalizes $S$, and preserves $W$.

For a finite real $W$-valued Borel measure $\mu$ on compact $K$, define
$$
B_\mu=\int_K k\cdot B_{d\mu(k)} .
$$
This is a finite-dimensional vector-measure integral: choose a basis of
$W$, integrate each continuous field against the corresponding scalar
signed measure, and sum. Its definition is independent of that basis.
Finite real sums are the atomic cases.

For $x=gK$ use the representative $g$ that is symmetric positive
definite with determinant one. Let $g^{-1}k=Q(g,k)R(g,k)$ be the QR
factorization with $Q(g,k)\in K$ and $R(g,k)\in S$. Then
$gQ(g,k)=kR(g,k)^{-1}$ is in $kS$. This and the $S$-invariance of
$B_b$ give the exact transported profile at $o$,
$$
g^{-1}\cdot B_\mu(x)=\int_K Q(g,k)\cdot d\mu(k).
$$

These fields are real analytic. Near any fixed real positive $g_0$, the
leading column Gram determinants of $g^{-1}k$ are strictly positive
uniformly over compact $K$. After shrinking a complex coordinate
neighborhood of $g_0$, each determinant stays in a right half-plane.
Use its holomorphic positive square-root branch in the QR recurrence,
with the Gram products continued bilinearly. On a smaller neighborhood
all integrands have uniform bounds. Integration against finite variation
therefore preserves holomorphicity by the Cauchy formula, giving real
analyticity. The two hook traces and the linear exterior covariant
derivative pass under the integral by locally uniform derivative bounds.

## BW2. A common limiting flag

Fix $\lambda_0>\cdots>\lambda_4$ with sum zero. For $u\in K$ let
$$
g(t)=u\operatorname{diag}(e^{-\lambda_0t},\ldots,e^{-\lambda_4t})u^T .
$$
For any fixed $k$, each leading principal minor of $u^Tk$ of orders
one through four is real analytic in $u$ and is one at $u=k$.
Its zero set has Haar measure zero. Applying Fubini to Haar measure
and the total variation $|\mu|$ gives one $u$ with all four minors
nonzero for $|\mu|$-almost every $k$.

In the $u$ frame, the wedge of the first $j$ columns of $g(t)^{-1}k$
has coefficient on the ordered blade with indices $I$ equal to
$$
e^{t\sum_{i\in I}\lambda_i}\Delta_{I,j}(u^Tk),
$$
where $\Delta_{I,j}$ uses rows $I$ and the first $j$ columns. The
strict ordering of the exponents gives a unique largest term at
$I=(0,\ldots,j-1)$. Thus the normalized wedge converges to the
leading coordinate blade times the sign of that minor. The positive
diagonal of $R(g,k)$ ensures the same oriented wedge for the first
$j$ columns of $Q(g,k)$. In succession the oriented wedges force
$Q(g(t),k)\to uD(k)$, where $D(k)$ is diagonal. Its $j$-th sign is
the ratio of the leading signs at orders $j$ and $j-1$ in
one-based indexing, with both boundary signs at orders zero and
five equal to one. Hence $D(k)\in M$.

The integrands are orthogonal actions on a positive tensor space.
Dominated convergence gives a limiting transported value in $uW$,
since $D(k)$ preserves $W$.

## BW3. The full equation

On the declared branch the pointwise quadratic equations are
equivariant and the grade-one baseline is a positive constant in
every physical frame. If the full equation held on a nonempty
open set, the analytic pointwise tensor residual would vanish
throughout the connected observer. It would also vanish for
the transported limit in $uW$. Rotating back, BV's positive
receiver forces the limiting value to be zero. That value fails
the positive baseline, giving a contradiction.

Every arbitrary $G$ translate of an individual $B_b$ is a compact
translate because $G=KS$ and $S$ fixes $B_b$. Thus finite sums of
these more general translates are included. Positive observer scales
transfer by BC's dilation and homothety. The averaging uses a fixed
finite total-variation measure with values in $W$.
[BZ](BZ-h6-triangular-limit-separation.md) extends the limit to
finite rotation derivatives and distributions with support in
a common regular QR cell. [CE](CE-h6-oscillating-triangular-fields.md)
allows sufficiently small imaginary character frequencies in a
joint finite rotation measure. General closed fields remain open.
External peer review remains outstanding.
