# DM. The complete H6 first-derivative tableau

Let \(V\) be the 14-dimensional tangent representation, and let
\(K\) be the 6,820-dimensional H6 value space after the universally
forbidden \((2,0)\) type is removed by
[BU](BU-h6-second-jet-obstruction.md). The exterior symbol
\(d\colon V\otimes K\to V\otimes\Lambda^7V\) has rank 17,251.
Its kernel, the first-derivative tableau \(g_1\), has dimension
78,229. Contraction maps \(g_1\) onto a 15,953-dimensional
subspace of \(V\otimes\Lambda^5V\). The joint closed and
coclosed tableau therefore has dimension 62,276.

Every smooth closed hook has covariant derivative in \(g_1\).
The tableau is the complete space satisfying these linear
first-order conditions. Actual first derivatives can also face
higher-order compatibility conditions.

## DM1. The complete decomposition

Use the exact [BO](BO-h6-full-curvature-kernel.md) fiber inventory.
The 95,480-dimensional tensor product \(V\otimes K\) has
38 highest output types. The verifier constructs their exact
multiplicities using the B2 Racah–Speiser rule and checks the total
dimension of every input tensor product.

The portable archive supplies all highest tensors in the
raw BO coordinates \((i,a,M)\), with derivative index \(i\),
coefficient index \(a\) and six-form blade mask \(M\).
For every tensor the verifier checks its total weight, both raising
equations, both hooks and coefficient-curvature equation on every
derivative slice, and the fiber Casimir. The Casimir values of
distinct fiber types in the inventory are distinct. Independence
of the raw highest columns verifies independence across multiple
copies. The exact character multiplicities give completeness.

## DM2. Exterior and contraction ranks

In each highest block apply exterior multiplication and contraction
directly to the raw tensors, with the BO wedge signs and inverse
metric. Reduce the exterior image matrix over the rationals.
Independently check the rank of a displayed rational kernel basis
and its zero exterior product. Restrict contraction to that entire
kernel and compute its rational image rank.

Equivariance and compact complete reducibility extend the highest
space computations to the full representations. Weighting each
rank by the Weyl dimension gives the stated totals.
The 38 block ranks and kernel multiplicities are in the
[verifier receipt](../../computation/geometric-unity/DM_h6_first_tableau.json).

The calculation recovers DK's three forced-zero derivative
channels from the \((5,5)\) fiber. Further forced-zero channels
include both \((5,2)\) fibers into output \((5,3)\), the
\((5,1)\) fiber into \((4,0)\) and \((3,1)\), the \((3,2)\) fiber
into \((3,3)\), and the \((2,1)\) fiber into \((1,1)\).
Contraction on \(g_1\) has substantial nonzero image, so the
algebraic closure condition alone leaves coclosure freedom.

    python -X utf8 computation/geometric-unity/DM_h6_first_tableau.py
    python -X utf8 computation/geometric-unity/DM_h6_first_tableau.py --full

The first command checks five representative blocks.
The second checks all 38 blocks, the full rational membership,
the exact kernel basis and the complete rank totals.
