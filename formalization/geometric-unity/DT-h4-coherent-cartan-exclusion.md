# DT. Real exclusion on an H4 two-weight slice

In the complexification of the 10,388-dimensional H4 curvature
kernel, let \(L\) be the subspace with two Cartan weights \((3,3)\)
annihilated by the \(SO(3)\) rising operator on the last three
defining axes. Use the \(SO(3)\) convention of
[DN](DN-h4-high-cartan-real-exclusion.md). The complex dimension
of \(L\) is 28. On the 56-dimensional real space
\(W=\{z+\bar z:z\in L\}\), one rational receiver gives
\[
 \langle Y,q_*(B)\rangle>\tfrac{19}{20}\|B\|^2,qquad
 \|q_*(B)\|^2>\tfrac15\|B\|^4
 \quad(0\ne B\in W).
\]
Here \(q_*\) is the complete homogeneous H4 Clifford residual.
Equivariance extends the exclusion to every compact translate
of \(W\).

The two-weight complex slice itself consists of homogeneous
quadratic zeros. On \(L\), every polarized value would have
weight \((6,6)\). The maximum sum of the two weights of any
H4 receiver is ten, as follows by counting the largest
one-plane weights of the equal-speed double rotation.
The norm-correction term has the same weight restriction.
Compact conjugation adds a negative-weight part in \(W\),
whose mixed products contribute to the zero-weight receivers.

## DT1. The complete slice and its receiver matrix

The witness supplies compact rational real and imaginary parts of
28 complex highest vectors. Each pair passes both hook equations,
dual curvature, both Cartan eigenrelations and the raising equation.
The [BX](BX-h4-curvature-and-closed-fields.md) character inventory
gives 28 for the difference between the full \((3,3)\) and \((3,4)\)
weight multiplicities. This is the required highest multiplicity.

The receiver is fixed by both compact rotations and by
\(s=\operatorname{diag}(1,-1,1,-1,1)\). Its grade-one trace is
zero. In the source, \(s\) is even on real parts and odd on
imaginary parts, so its invariant quadratic form has no cross
terms between them. Circle invariance equates the diagonal
blocks. The exact 28-by-28 real Gram of the receiver is \(A\),
and the physical norm Gram in the same basis is \(G\).
Exact integer elimination certifies \(20A-19G>0\). The receiver
squared norm is
\[
 \frac{665957917357558952343007011503833913}
      {163557351043847070298710634370433024}<\tfrac92 .
\]
Cauchy–Schwarz gives a squared residual lower coefficient
\(361/1800>1/5\).

## DT2. Reproduction

The verifier uses DN's rational denominator and integer arithmetic
construction. Ten independent primes give a 200-bit modulus
strictly above twice the 193-bit signed numerator bound.
The full-spin Clifford evaluator checks the anchor and a mixed pair.

    python -X utf8 computation/geometric-unity/DT_coherent_slice.py
    python -X utf8 computation/geometric-unity/DT_coherent_slice.py --full

Both modes verify every geometric and positivity claim.
The default checks one complete modular image. Full mode
reconstructs all ten. The real domain uses the highest-weight
slice of the indicated \(SO(3)\) type. The full type and its
sum with DN's larger high-Cartan domain remain open.
