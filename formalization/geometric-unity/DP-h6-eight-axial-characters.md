# DP. A common H6 axial obstruction on eight sign characters

For the positive observer, consider the [DS](DS-h6-axial-multisign.md)
axial triangular fields \(B(sK)=s\cdot b(\chi(s))\) whose values
belong to the sum of sign masks \(0,3,5,6,9,15,17,27\).
One exact common positive receiver gives
\[
 F_y(b)>\|b\|^2/1250
\]
for every nonzero local closed value in that sum. It excludes full
solutions on any open axial interval. Its positive Gram holds on
the complete 1,165-dimensional local value space.

## DP1. Five additional complete blocks

The new characters are \(3,5,6,9,17\). Exact character averaging
in [BO](BO-h6-full-curvature-kernel.md)'s inventory gives a
423-dimensional curvature kernel in each. The triangular-closure
pencil \(A+zD\) is the one in DS. Reconstruct and check all
polynomial coefficients in its kernel, with chain counts
\((3,27,19,6)\) in degrees zero to three. Their coefficient span
has dimension \(138\), and its image under \(D\) has dimension
\(83\). The exact quotient joint rank has lower bound
\(570=2(423-138)\), verified modulo 65,537 using a saturated
modular curvature kernel as in DS. This forces the two quotient
columns to be individually injective with disjoint images. Any
local closed \(C^1\) profile therefore has value in the coefficient
span.

| Sign mask | Curvature dimension | Local value dimension | Chain image dimension | Quotient joint rank lower |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 389 | 107 | 64 | 564 |
| 3, 5, 6, 9, 17, each | 423 | 138 | 83 | 570 |
| 15 | 443 | 187 | 105 | 523 |
| 27 | 443 | 181 | 110 | 524 |

The three DS blocks retain their complete local mode classification.
In particular the eleven additional real exponential modes of
character 15 are included in the new Gram.

## DP2. One receiver on the combined spaces

Use the same orthonormalization, rational rounding metric bounds
and exact Gram reconstruction as DS with the separately archived
DP receiver and eight exact bases. The smallest of the eight
certified rational margins is greater than \(1/1250\). The receiver
is invariant under the diagonal sign group, so its polarization
vanishes between distinct characters. This proves the inequality
on their entire direct sum without restricting the relative sizes
or phases of components. Coefficients on the positive observer use
the full homogeneous Clifford residual.

The DS escaping-ray argument also applies to finite-variation
mixtures of compact rotations of bounded profiles in this enlarged
domain, under the same measurability and uniform integrability
hypotheses. The result covers these eight sign characters.

    python -X utf8 computation/geometric-unity/DS_h6_axial_multisign.py --eight-signs --character 6 --geometry-only
    python -X utf8 computation/geometric-unity/DS_h6_axial_multisign.py --eight-signs --full

The selected command checks the complete linear claim in character
6. The full command verifies every coefficient and every modular
Gram image across all eight characters.
