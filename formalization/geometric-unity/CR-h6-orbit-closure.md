# CR. Compact-orbit first-jet obstruction for the H6 partial value

On the positive H6 observer, the first derivative of a field taking values in one compact \(SO(5)\) orbit must lie in the tangent bundle to that orbit within the tensor bundle. At the certified BS value, covariant closure annihilates every such first derivative. The positive second-order curvature source of CJ then excludes a local closed field whose values stay in the orbit of the BS point and satisfy its scalarity and constant norm.

The closure condition on compact-orbit first derivatives is injective on a nonempty Zariski-open subset of the complete pointwise H6 curvature kernel. It is not an injectivity statement for all points of that kernel, and it does not say that every solution of the full quadratic equations lies in one compact orbit.

## The orbit-closure operator

Let \(V=T_oY\) with its physical positive metric, and let \(K_6\) be the BO pointwise real hook and coefficient-curvature kernel. At \(b\in K_6\) define
\[
A_b\colon V^*\otimes\mathfrak{so}(5)\longrightarrow V^*\otimes\Lambda^7V^*,
\qquad
A_b\!\left(\sum_i\theta^i\otimes J_i\right)_a
=\sum_i\theta^i\wedge(J_i\cdot b)_a.
\]
The derivative index is separate from the external coefficient index \(a\). The action \(J_i\cdot b\) acts in all seven covariant slots of \(b\). In this display any covariant basis \(\theta^i\) may be used, with the input still measured by the induced metric on \(V^*\) and the Euclidean unit-plane metric on \(\mathfrak{so}(5)\). The verifier calibrates the BC generator convention against the ten physical plane rotations of the defining five-dimensional representation. The six internal rotations need no change, and the four remaining BC rotations are multiplied by two. Order and sign choices of the orthonormal generating planes make no change to the Gram positivity bound.

The induced connection on the real tensor bundle preserves the compact reduction \(SO(5)\) and the parallel kernel \(K_6\). A local section whose values remain in a fixed \(SO(5)\) orbit in this reduction has covariant derivatives in the orbit tangent. This follows by taking a compact gauge local to the point and differentiating the orbit-valued section, including the connection contribution in the Lie algebra. Thus write \(\nabla_i B=J_i\cdot b\) at that point. Exterior covariant closure gives \(A_b(\sum_i\theta^i\otimes J_i)=0\). If \(A_b\) is injective then \(\nabla B=0\) there. If the action has a nontrivial stabilizer, the appropriate weaker condition is injectivity on \(V^*\otimes T_b(SO(5)b)\), since one may quotient the stabilizer in the domain.

The dimension of the domain used here is \(14\cdot10=140\). In fixed rational coframes every entry is linear in \(b\). The verifier computes one exact 140-by-140 coefficient minor at the BS center and confirms its determinant is nonzero modulo 65537. This proves that the injective set is a nonempty Zariski-open subset of \(K_6\). The rank calculation is equally valid over a linear subspace containing that center. It does not give a codimension greater than one for the rank-drop locus.

## Uniform estimate on the BS box

Use the 140 columns of \(A_b\) in the physical BC covariant coframe. Its output coefficient at index \((a,m)\) has squared metric weight \(1/(g_a\prod_{j\in m}g_j)\). The Gram on the input coefficients is the diagonal matrix \(D\) that repeats \(1/g_i\) on the ten Lie generators at each derivative index \(i\). At the BS rational center the verifier computes the entire Gram \(A_b^*A_b\) in coefficient coordinates and checks 140 positive rational pivots of \(A_b^*A_b-\tfrac12D\). In particular the center singular bound exceeds \(\sqrt{1/2}\) relative to the stated domain metric.

For any variation \(\Delta b\), each unit-plane generator has action norm at most 2 on \(V\), since \(V=\operatorname{Sym}_0^2(\mathbb R^5)\) and its action is induced from two copies of the defining representation. On the seven covariant slots of \(b\) this gives action norm at most 14. Cauchy–Schwarz over the ten Lie generators and the norm \(\sqrt7\) of exterior alternation \(V^*\otimes\Lambda^6 V^*\to\Lambda^7V^*\) yield
\[
\|A_{\Delta b}\|_{\mathrm{op}}\le 14\sqrt{70}\,\|\Delta b\|<120\|\Delta b\|.
\]
The box has sup-coordinate radius \(10^{-8}\). Using the inverse Gram, the verifier evaluates the eight affine direction squared norms and finds rational upper bounds for their square roots. Therefore the whole box perturbation has operator norm at most \(2445963/250000000000\), by the triangle inequality. This gives the uniform singular bound \(\sigma_{\min}(A_b)>3/5\) throughout the BS box, including its certified partial point. The bound is with respect to Lie generators, rather than a choice of normalization of the orbit tangent basis. The conclusion of injectivity is independent of this normalization.

## Consequence for a field confined to one orbit

The map \(A_b\) is \(SO(5)\)-equivariant, so injectivity extends to the entire orbit of any tested value. A covariantly closed field with values in one such orbit has \(\nabla B=0\) wherever it is defined. The orbit has constant norm. If it also has scalar grade one, the double-divergence identity of CJ requires the curvature source to vanish at a point with zero covariant first derivative. Instead CJ gives a strictly positive curvature source for every nonzero real H6 point value in \(K_6\). This proves the exclusion for the BS orbit, whose point value has scalar grade one and norm one, and for any other nonzero grade-one scalar orbit on which the orbit closure operator is injective.

In particular a full solution through a regular algebraic zero of the complete pointwise equations would be locally excluded if the kernel of the differentiated quadratic equations were exactly the orbit tangent and the orbit closure operator were injective there. The latter holds generically on \(K_6\); that generic statement does not specify the geometry of the unknown nonlinear zero set. Neither such a zero nor that kernel condition is part of this result.

## A real curvature-compatible counterexample to universal injectivity

Take the first rational CL carrier \(b\), with 222 nonzero covariant coefficients. It lies in \(K_6\) and satisfies genuine AN closure. Every form blade of \(b\) involves \(e^{13}\). The physical compact generator \(J_{01}\) fixes \(e^{13}\), so all form blades of \(J_{01}\cdot b\) also involve \(e^{13}\). This orbit derivative is nonzero: the covariant coefficient at \((a,m)=(4,8433)\) is exactly \(-48\). Hence \(e^{13}\otimes(J_{01}\cdot b)\) is a nonzero tangent-orbit first derivative with zero symbol. This remains a counterexample after quotienting the stabilizer. It is an algebraic first-jet statement and need not correspond to an actual field confined to the orbit of that carrier.

## Verification

Run from the root of the repository. Python, NumPy, and SymPy suffice.

```text
python -B computation/geometric-unity/CR_h6_orbit_closure.py
python -B computation/geometric-unity/CR_h6_orbit_closure.py --full
```

The default checks the complete generator calibration, selected nonzero exact minor, BS direction lengths and perturbation estimate, and the exact CL counterexample. The full option also verifies the shifted Gram at the BS center, completing the box bound.

[Review106](../../reports/reviews/REVIEW106-GU-CR-CS-CV-2026-11-14.md) records the independent checks. External peer review remains outstanding.
