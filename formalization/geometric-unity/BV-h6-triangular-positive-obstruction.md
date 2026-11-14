# BV. A positive receiver excludes the full triangular closed H6 family

On the fixed flat positive branch, the entire 952-dimensional space
of closed $H_6$ fields invariant under the simply transitive positive
upper triangular subgroup is excluded from the full field equation.
The proof gives an exact linear combination of necessary quadratic
coefficients that is positive on every nonzero real value in this
space.

The [verifier](../../computation/geometric-unity/BV_h6_an_positive.py)
reconstructs the rational field columns and the exact dyadic
certificate. [Review86](../../reports/reviews/REVIEW86-GU-BV-2026-11-14.md)
records the independent checks.

Retain [BR](BR-h6-closed-field-lifts.md)'s frame and sign characters
and [BS](BS-h6-certified-partial-value.md)'s complete physical
quadratic formula. Let
$$
q^0(B)=q(B)+(26/21)\|B\|^2\gamma
$$
be the homogeneous necessary quadratic in the orthonormal coframe.
Here $\gamma$ denotes the canonical degree-one tensor. The added
term removes the grade-one trace. On the full equation, all its
degrees vanish.

## BV1. A shared invariant receiver

The adjacent data give a dyadic receiver $y$ in the positive unit
coframe, with denominator $2^{36}$ and integer numerators. Its
grade-one trace is zero. Every coefficient has equal sign character
in its external and exterior slots, so it is invariant under the
even sign group. Consequently, for input values of distinct
characters the symmetric polarization of $\langle y,q^0(B)\rangle$
vanishes. The same receiver can therefore be verified on the sixteen
individual character blocks.

The receiver was located with a numerical search on the exact
closed spaces. The proof uses its fixed dyadic coefficients and
the exact inequalities below.

## BV2. Exact evaluation on rounded matrices

For each sign block, let $Z$ be a full rational closed coordinate
map from BR's stack. The certificate selects rational linear
combinations $B=ZC$ and checks membership in the exact stack.
With positive square-root intervals for the Gram factors, it
constructs a dyadic unit-frame matrix $D$ with denominator $2^{36}$.
The errors and norms satisfy
$$
\|B-D\|_{\rm F}<10^{-8},\qquad \|D\|_{\rm op}<2,\qquad \|y\|<11.
$$
Here $B$ in the norm denotes its exact conversion to the unit frame.
The first inequality uses integer square roots and the second uses
exact Gram row sums. The receiver norm bound is a sum of squares.

Let $H(D)$ be the symmetric matrix of the quadratic
$x\mapsto\langle y,q^0(Dx)\rangle$ using the polynomial
extension of BS1 on the rounded tuple. Its integer multiple
$J=3\cdot 2^{108}H(D)$ has no input-dependent denominator,
because the grade-one correction pairs to zero and the remaining
Clifford quadratic has denominator three in the unit frame.
To evaluate it, use the four even-sign Clifford lifts. Their
binary commutation form is nondegenerate of rank four. Selecting
a symplectic Pauli coordinate basis puts their matrix algebra
in a four-dimensional factor of the spin space. In a fixed
input character, each external component has one fixed Pauli
factor on those two qubits. Products multiply the scalar phases
there and use 32-dimensional Clifford matrices on the remaining
factor. This reorganization preserves the exact trace pairing.
The verifier checks the phase conventions and the complete targeted
quadratic contraction.

For each block, compute $J$ modulo six distinct primes, with
the required Pauli phases realized in the prime field. The
pairwise products of reduced integer matrices remain exactly
representable in binary64. The inner dimensions and prime sizes
give a strict bound below $2^{53}$ for every dot product.
Chinese remaindering and an explicit bound on $|J_{ij}|$ determine
the unique signed integer matrix. Fraction-free elimination
checks all leading minors of $25J-3\cdot2^{108}I$. They are
strictly positive for every block. Thus
$$
H(D)\succ \frac1{25}I .
$$

## BV3. Transfer to the exact closed spaces

Write $\|\cdot\|_{\rm c}$ for the coefficient norm, the
Hilbert–Schmidt norm divided by $\sqrt{128}$. For a product,
$\|XY\|_{\rm c}\le\sqrt{128}\|X\|_{\rm c}\|Y\|_{\rm c}$.
In the full BS1 quadratic, the anticommutator, the two external
mixed terms, and the middle sandwich give
$$
\|q^0(v)\|_{\rm c}
\le\left(\frac{10}{3}\sqrt{14\cdot128}
+\frac{26}{21}\sqrt{14}\right)\|v\|^2
<150\|v\|^2.
$$
The estimate holds for the extension to arbitrary grade-eight
tuples, with the same norms. For either mixed term, Cauchy–Schwarz
in the contracted index gives a component bound
$\sqrt{14\cdot128}\|v_a\|\|v\|$. The other three contributions
use the triangle inequality and $\sum_a\|v_a\|^2=\|v\|^2$.
Polarization and the parallelogram identity give a corresponding
bilinear bound $150\|v\|\|w\|$.

For unit $x$, the difference between the exact closed square and
the rounded square has receiver pairing of absolute value at most
$$
150\cdot 11\cdot 10^{-8}(4+10^{-8}) < 10^{-4}.
$$
Therefore $H(B)\succ(399/10000)I$ on each selected coordinate
system. Positivity implies their independence. Each system has
the complete dimension of its closed block by BR, so the exact
closed maps span it. Combining the blocks proves positivity on
the entire real 952-dimensional space.

Every full field must satisfy the homogeneous quadratic, giving
zero pairing with $y$. Within the declared triangular family,
this forces the zero field value, which fails the positive
grade-one baseline. The result includes arbitrary mixtures of
the sixteen triangular sign characters. Finite sums of fields
rotated by different compact group elements have a larger value
space. [BW](BW-h6-triangular-orbit-averages.md) excludes those
fields from the full equation via a limiting value in a rotation
of the invariant space. General full fields and physical
identification remain open.
