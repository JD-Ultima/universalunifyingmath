# CD. Clifford contractions and operator bounds for H4 and H6

The real pure H4 and H6 quadratic outputs obey an exterior trace
identity in every degree. A Clifford contraction also gives a
positive operator bound at any full solution. These statements use
the real hook conditions on a positive 14-dimensional frame.

## CD1. Exterior and contraction traces

Let $\gamma_a$ be Hermitian Clifford generators. Write $C_a=B_a$
for H4 and $C_a=B_a\Gamma$ for H6, with $B$ a real hook and
$\Gamma^2=-1$. Then the homogeneous Clifford degree of $C$ is
$k=4$ or $8$, each $C_a$ is Hermitian, and
$\sum_a\gamma_aC_a=\sum_a C_a\gamma_a=0$. Set
$$
 S=\sum_a C_a^2,\quad N=\operatorname{tr}(S)/128,\quad
 J(X)=\sum_a \gamma_aX\gamma_a,\quad
 V_0=\sum_{a,b}\gamma_a C_b C_a\gamma_b .
$$
In this section superscripts in parentheses denote Clifford
projections, including $V_0^{(r)}$ for the grade-$r$ part of the
operator $V_0$.

From the complete reduced quadratic of [BS](BS-h6-certified-partial-value.md)
and [CB](CB-h4-certified-partial-value.md), contraction gives
$$
 \sum_a\gamma_a q_a=-\frac23\big((2k-2)S+J(S)-V_0\big).
$$
In this calculation
$\sum_{a,b}\gamma_a\gamma_b C_a C_b=2S$
by the hook sums and the Clifford relation.
The sandwich term gives
$\sum_{a,b}\gamma_a C_b\gamma_a C_b=(14-2k)S$.
Exchanging $a,b$ after reversal shows that $V_0$ is Hermitian.
Thus the displayed contraction is Hermitian, and the reverse
product gives its same adjoint.

For an odd homogeneous symbol, left minus right Clifford action
of $\gamma_a$ is twice exterior multiplication. As the possible
grades of $q_a$ are $1,5,9,13$, their distinct exterior degrees give
$$
 \sum_a e^a\wedge q_a^{(r)}=0,\qquad r=1,5,9,13 .
$$
In particular, degree one is symmetric. Since
$J(S^{(r)})=(14-2r)S^{(r)}$, the contraction traces are
$$
 \sum_a\iota_a q_a^{(r+1)}
 =-\frac23\big((14+2k-2-2r)S^{(r)}-V_0^{(r)}\big),
 \qquad r=0,4,8,12 .
$$
Cyclicity gives $V_0^{(0)}=2N$, so the grade-one trace is
$-12N$ for H4 and $-52N/3$ for H6. The middle contractions
retain both terms. Exact controls certify that their two terms
need not be proportional.

## CD2. The orthogonal output span

For $d=4,6$ put $\Omega_{2d}=\sum_b B_b\wedge B_b$.
The top exterior terms give
$$
 q_a^{(9)}=-\frac23e^a\wedge\Omega_8\quad\text{in H4},
 \qquad
 q_a^{(13)}=2e^a\wedge\Omega_{12}\quad\text{in H6}.
$$
In either case the top contraction is $-4S^{(2d)}$, or
equivalently $V_0^{(2d)}=-2S^{(2d)}$. For a tensor of
form degree $r$ whose exterior trace vanishes, the contraction
summand is the exterior product with its contraction trace
divided by $15-r$. Subtracting it leaves the two-trace-free
orthogonal hook $H_r$.

Consequently the complete linear spans of the quadratic outputs
on the full real hook domains are, respectively,
$$
 \mathbb R\oplus\operatorname{Sym}_0^2V
 \oplus\Lambda^4V\oplus H_5\oplus\Lambda^8V
$$
and
$$
 \mathbb R\oplus\operatorname{Sym}_0^2V
 \oplus\Lambda^4V\oplus H_5\oplus\Lambda^8V
 \oplus H_9\oplus\Lambda^{12}V .
$$
Here $V=\mathbb R^{14}$ with its positive inner product.
The harmonic tensor classification gives the listed irreducible
orthogonal hooks. See [Goodman and Wallach, chapter 10](https://webhomes.maths.ed.ac.uk/~v1ranick/papers/goodwallx.pdf).
The summands are inequivalent irreducibles of $O(14)$.
Form-slot Hodge star identifies $H_9$ with the determinant
twist of $H_5$. The restriction of $H_5$ to $SO(14)$ is
irreducible, since its orthogonal harmonic Young diagram is
off the middle degree. An intertwiner between the two full
orthogonal extensions would have to be scalar on this
restriction, and a reflection then forces that scalar to
be zero. The quadratic H6 output remains
$O(14)$-equivariant because the determinant factor from the
Clifford volume occurs twice. The verifier exhibits a nonzero
projection onto every summand. Complete reducibility and
inequivalence then prove the spans. Their dimensions are 28,133
and 52,248, respectively. Removing the scalar trace reduces
each dimension by one. These dimensions concern all hooks.

## CD3. An operator sum of squares

Let $A_{ab}=\iota_a C_b-\iota_b C_a$ and
$P=\sum_{a,b}A_{ab}^\dagger A_{ab}$. Each $A_{ab}$ is
anti-Hermitian. Using the hook sums and
$\iota_a C_b=(\gamma_a C_b-C_b\gamma_a)/2$ gives
$$
 2P=-(14-4k+2)S+J(S)-V_0,\qquad
 \sum_a\gamma_a q_a=-\frac23\big((14-2k)S+2P\big).
$$
For the first identity, expand $P=-\sum_{ab}A_{ab}^2$.
The two squared interior sums use $J(S)$ and
$\sum_a\gamma_aC_b\gamma_a=(14-2k)C_b$. The crossed
interior sum is $-(V_0+2S)/4$ by the hook sums and
the Clifford relation.

On antisymmetric spinor-valued arrays let
$(TA)_a=\sum_b\gamma_b A_{ab}$ with the inner product summing
over all ordered pairs. Its adjoint is
$(T^*c)_{ab}=(\gamma_b c_a-\gamma_a c_b)/2$.
On $\gamma$-traceless $c$, direct contraction gives $TT^*c=6c$.
The hook sums give $T(A)=-(k-1)C$. Subtracting the
$T^*$-component and completing the square yields
$$
 P-\frac{(k-1)^2}{6}S
 =\sum_{a,b}U_{ab}^\dagger U_{ab}\succeq0,\qquad
 U_{ab}=A_{ab}+\frac{k-1}{12}(\gamma_b C_a-\gamma_a C_b).
$$
The spinor argument applies to every spinor and therefore gives
the operator identity.

At a full solution the higher grades vanish and grade one is
$-(6/7)N\gamma_a$ for H4 or $-(26/21)N\gamma_a$ for H6.
Substituting gives
$$
 \begin{array}{ll}
 \mathrm{H4}:&P+3S=9N\,\mathbf1,\qquad S\preceq2N\,\mathbf1,\\
 \mathrm{H6}:&P-S=13N\,\mathbf1,\qquad S\preceq(78/43)N\,\mathbf1.
 \end{array}
$$
The upper-cap eigenspace is the common kernel of all $U_{ab}$.
The Hermitian real coefficient expansion of $S$ has degrees
$0,4,8,12$. In particular the chirality trace difference is
zero, so each 64-dimensional chirality space has trace $64N$.
For $N>0$ the two caps give ranks at least 64 and 72 in
the full spin space.

For $p>1$ the convex scalar maximum at fixed trace and cap uses
as many cap eigenvalues as possible and at most one intermediate
eigenvalue in each chirality space. It gives
$$
 \frac{\operatorname{tr}(S^p)}{128N^p}\le
 \begin{cases}
 2^{p-1}&\text{in H4},\\
 \big(35(78/43)^p+(22/43)^p\big)/64&\text{in H6}.
 \end{cases}
$$
For equality the H4 blocks require 32 cap eigenvalues and
32 zeros each. The H6 blocks require 35 cap eigenvalues,
one $22N/43$ eigenvalue and 28 zeros each.

## Verification

    python -X utf8 computation/geometric-unity/CD_clifford_hook_identities.py

The helper compares full quadratic outputs of four rational
hook samples against the original nested Clifford formula.
It checks every exterior and contraction identity, both top
formulas and the full operator sum of squares.
[Review95](../../reports/reviews/REVIEW95-GU-CD-2026-11-14.md)
records the checks. General pointwise and field solutions
remain open. External peer review remains outstanding.
