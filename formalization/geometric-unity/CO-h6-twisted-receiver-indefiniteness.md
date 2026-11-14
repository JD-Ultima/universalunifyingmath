# CO. Exact finite-frequency indefiniteness of the BV receiver on a twisted H6 kernel

At the imaginary character with differential \(4i\,e^{13}\), the
[BV](BV-h6-triangular-positive-obstruction.md) receiver has both
Hermitian signs on the trivial diagonal-sign block of the exact
twisted hook and closure kernel from [CE](CE-h6-oscillating-triangular-fields.md).
The certificate gives actual linearly closed field amplitudes.

## Frequency convention and witnesses

Use the real additive character from CL,
\[
\chi(s)=-\tfrac14\sum_{j=0}^{3}\log s_{jj}+\log s_{44}.
\]
The unitary character is \(e^{4i\chi(s)}\), corresponding to a CE log-diagonal frequency \(\zeta=(-1,-1,-1,-1,4)\) in zero-based axis order. Its derivative at the identity is \(4i\,e^{13}\). The real covector \(4e^{13}\) has squared physical norm \(5\).

Write \(W(4)\) for the complex hooked tensor space with
\[
(D_S+4i\,e^{13}\wedge)b=0.
\]
Then \(e^{4i\chi(s)}s\cdot b\) is a globally analytic complex linearly closed field on the positive \(SL(5,\mathbb R)/SO(5)\) model. The real part is also closed. Both the alternating and contraction hook equations remain real-linear on its components.

The archive has one negative witness \(b=u+iv\), with 2638 integer covariant coefficients each for \(u\) and \(v\), and one positive real witness with 762 integer coefficients. The verifier independently checks the complete hook equations and coefficient-curvature wedge for every component. It verifies \(D_Su=4e^{13}\wedge v\) and \(D_Sv=-4e^{13}\wedge u\) in the exact BR coframe. The positive witness lies in the CL carrier intersection, so both \(D_S b\) and \(e^{13}\wedge b\) vanish on it. The exact character of each witness is the trivial diagonal-sign block.

## Physical receiver and exact signs

Let \(F_y\) be the real homogeneous BV receiver, and let \(H_y\) be its Hermitian extension. For real \(u,v\),
\[
H_y(u+iv,u+iv)=F_y(u)+F_y(v).
\]
We use exactly the BV receiver vector and its same rational dyadic coefficients in the physical orthonormal frame. The verifier pins the BV data by SHA-256.

For each component it first multiplies the covariant degree-six tensor by the covariant Clifford volume; the resulting degree-eight polynomial is converted to the unit-volume normalization by the factor \(\det g\) in the BN receiver evaluator. Then it adds the scalar adjustment \((26/21)\|b\|^2\gamma\) as in BV. For a receiver label \((a,m)\), its orthonormal coefficient is its covariant coefficient divided by \(\sqrt{g_a\prod_{j\in m}g_j}\). This conversion yields finitely many square-root classes with rational weights. The verifier computes every weighted coefficient exactly, checks selected Clifford entries by independent nested multiplication, and encloses the residual radicals by rational integer square-root bounds.

For the negative witness the normalized Hermitian value is strictly between \(-0.314427\) and \(-0.314426\). The second witness has strictly positive Hermitian value. Positivity also follows for that witness from BV because it already belongs to \(W(0)\).

This shows the BV receiver cannot extend the CE positivity proof
over the full frequency range. A different receiver would be
needed at this frequency. General nonlinear fields remain open.
[CL](CL-h6-oscillating-carriers.md) gives an all-frequency exclusion
on the smaller common carrier space.

## Reproduction

Run the [verifier](../../computation/geometric-unity/CO_h6_twisted_receiver.py)
with the repository's pinned dependencies.

```text
python -X utf8 computation/geometric-unity/CO_h6_twisted_receiver.py
```

[Review102](../../reports/reviews/REVIEW102-GU-CO-2026-11-14.md)
records the exact checks. External peer review remains outstanding.
