# BZ. Distributional triangular fields and a uniform separation bound

On the pure H6 branch of [BW](BW-h6-triangular-orbit-averages.md),
compactly supported rotation distributions have the same obstruction
when their support lies in one regular QR cell. In particular, no
finite linear combination of finite-order rotation derivatives of
triangular closed fields solves the full equation on an open set.
There is also an explicit uniform approximation gap on the entire
observer for fields in these families.

## BZ1. Distributional profiles

Keep BW's notation for $G,K,S,W$ and its analytic profile
$Q(g,k)\cdot b$, with $b\in W$. For a real $W$-valued distribution
$T$ of order $m$ on compact $K$, pair this profile against $T$.
Equivalently, write $T$ in a basis of $W$ and pair each kernel
component with the corresponding scalar distribution. The joint
analyticity from BW gives a uniform complex $g$ neighborhood on
each compact parameter chart. By reducing that neighborhood,
the Cauchy formula in the parameter variables gives uniform
bounds on all derivatives through order $m$. Pairing against
$T$ preserves holomorphicity in $g$. The resulting field is
globally real analytic. Its linear hook and exterior closure
equations follow by moving the required derivatives through
the distribution pairing.

Suppose there exists $u\in K$ such that all four leading minors
$m_j(k)$ of $u^Tk$, of sizes $1\le j\le4$, are nonzero on
$\operatorname{supp}T$. They are uniformly bounded away from
zero on a compact neighborhood of the support inside that cell.
Choose the strictly ordered diagonal ray from BW. Normalize each
exterior vector of its first $j$ columns after dividing by the
largest exponential factor. It takes the form
$$
  m_j(k)e_0\wedge\cdots\wedge e_{j-1}+E_j(t,k),
  \qquad E_j(t,k)=O(e^{-ct})
$$
for some $c>0$, uniformly on small complex parameter
neighborhoods of the support. Use the bilinear Gram norm there.
Its square is $m_j(k)^2+O(e^{-2ct})$. On each of finitely
many smaller neighborhoods, its positive real square root
extends holomorphically for all sufficiently large $t$.
The normalized wedge converges exponentially to the leading
coordinate wedge times the locally constant sign of $m_j$.
Cauchy estimates give convergence for every fixed finite
number of parameter derivatives.

The frame columns follow from the oriented wedges. Specifically,
contract the normalized wedge of the first $j$ orthonormal
columns against the ordered blade of the first $j-1$ columns.
This recovers column $j$ and uses only polynomial operations
on the coefficients. For the final column use the orientation.
Thus $Q(g(t),k)$ converges in $C^m$ near the support to
$uD(k)$, with locally constant $D(k)$ in the even diagonal
sign group $M$. A cutoff equal to one near the support permits
distribution pairing in these neighborhoods. The limit of the
transported field values belongs to $uW$ since every $D(k)$
preserves $W$. For any finite set of rotation parameters,
one common $u$ exists by avoiding finitely many minor zero
sets as in BW. This covers finite derivatives of point masses,
and hence finite sums of parameter derivatives of arbitrary
fixed finite order.

If the exact full equation held on a nonempty open set, real
analyticity would extend it to the connected observer. Its
transported limit in $uW$ would contradict BV's positive receiver
and the positive baseline.

## BZ2. An explicit separation

The exact [BV verifier](../../computation/geometric-unity/BV_h6_an_positive.py)
also checks the row-sum bound
$D^TD\prec(10001/10000)I$ for each reconstructed dyadic basis.
Together with $\|B-D\|_{\rm op}<10^{-8}$ this gives
$B^TB\prec(101/100)I$ for each exact closed basis. The polarized
receiver $F_y(v)=\langle y,q^0(v)\rangle$ therefore satisfies
$$
  F_y(w)>\frac1{26}\|w\|^2\quad (0\ne w\in W),\qquad
  |F_y(v,w)|\le1650\|v\|\|w\|.
$$
The first bound uses BV's coefficient lower bound
$399/10000$ and the orthogonal sign blocks. The second uses
the universal polarized coefficient estimate and $\|y\|<11$.

Let $z$ be a nonzero full pointwise zero and let $w$ be its
orthogonal projection onto $W$. Write $R=\|z\|$ and
$\delta=\|z-w\|$. If $w=0$ then $\delta=R$.
Otherwise $F_y(z)=0$ and the two bounds give
$$
  (R^2-\delta^2)/26 < 1650(2R\delta+\delta^2).
$$
Here the right-hand side follows by polarization and
$\|w\|\le R$. For $\delta\le R/100000$ this inequality
contradicts
$$
  1650\left(2/100000+1/100000^2\right)
  <(1-1/100000^2)/26.
$$
The same estimate holds on every rotation of $W$ by
orthogonal equivariance. For the full pointwise equation
$R^2=(21/26)(15/(2c_V)+15\kappa^2/52)>0$ by the grade-one
trace. Consequently every global full solution at those fixed
scales has uniform distance at least $R/100000$ from any
field with values converging to one rotated $W$ along a ray in
transported frames. This applies to both BW's finite-measure
family and the distributions of BZ1, and extends to their
global uniform closure by the triangle inequality.

[Review90](../../reports/reviews/REVIEW90-GU-BZ-2026-11-14.md)
records the checks. General distributions whose support meets
every candidate exceptional minor set, and approximation merely
on compact observer subsets, remain open. External peer review
remains outstanding.
