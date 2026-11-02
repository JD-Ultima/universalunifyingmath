# Internal review of AA7: the three-form hook zero locus

Scope: AA7 in `formalization/geometric-unity/AA-conic-traces.md`, its metric record and `AA_three_form_simple.py/.out`. This is an independent internal proof and implementation review, not external peer review or a formal proof-assistant verification. Reviews 20, 23 and 25 retain their earlier scopes.

## Assessment

AA7 is valid as stated for a real three-form on a positive-definite fourteen-dimensional space. Its hook equation holds exactly for zero or decomposable three-forms. The application to AA6 remains a necessary pointwise condition on its single Lambda3 c/wGamma kernel branch. No sufficiency for the other algebraic or differential field equations follows.

The proof was audited independently before packaging and checked again against the final note. The orthogonal complement in the last step explicitly excludes X as well as the two contraction-plane vectors; this ensures that its residual forms contain no x component. No mathematical narrowing was required.

## Critical proof steps

**Annihilator and positivity.** If a nonzero Y annihilates omega, contraction anticommutation gives `i_Y F=0`. The hook equation gives `Y-flat wedge F=0`; contracting it yields `|Y|^2 F=0`. Positivity is used here with the original metric. A null annihilator in an indefinite metric would not justify that conclusion.

**Exact wedge lemma.** In a Darboux basis on a symplectic support of dimension `2 ell`, the adjoint/wedge commutator is `(ell-j)I` on degree j. Counting each pair's absent, double and single occupancy gives coefficients +1, -1 and 0. Its norm identity therefore makes wedge multiplication injective for `j<ell`. The auxiliary positive Darboux norm proves an algebraic injectivity statement; it is never substituted for the original metric in F or in the hook equation.

**All ranks under F nonzero.** Nonzero F implies that at least one of its summands, a contraction square, is nonzero. This uses no positivity assertion about a sum of four-forms. For that unit X, the contraction has rank `2 ell>=4`, and the part of the hook equation without x is `theta wedge psi=0`.

- For `ell>=4`, injectivity kills every component of psi. The kernel of a two-form on the thirteen-dimensional perpendicular space is nonzero, giving an annihilator.
- For `ell=3`, components with a kernel index vanish. The form is supported in dimension at most seven, again leaving an annihilator.
- For `ell=2`, the only possible kernel-dependent part maps the nine-dimensional kernel into primitive two-forms on a symplectic four-space. The latter dimension is five: the six-dimensional space of two-forms maps nontrivially onto the one-dimensional fourth exterior power. The map therefore has a nonzero kernel, giving an annihilator.

These cases are exhaustive. Each contradicts F nonzero by the first step, so F vanishes and the hook equation reduces to Z=0.

**Simple-form step.** Contracting each equation `Z_X=0` by X makes every contraction square zero. A nonzero contraction consequently has rank two. In the final orthogonal splitting, the three-kernel-index residual vanishes first. The square of contraction by the first plane vector has the unique x-containing term displayed in the note, forcing sigma to vanish; the second forces tau to vanish. The remaining form is `u wedge v wedge (lambda*x+nu)`, with nonzero lambda. Conversely, a decomposable three-form has Z=F=0. This proves both directions, rather than only an obstruction on selected examples.

## Exact reproduction and limits

The fresh command

```text
python -X utf8 -I computation/geometric-unity/AA_three_form_simple.py
```

passed all 24 reported checks in under a second with Python 3.12. They include all 625 relevant basis commutators, the primitive dimension, projector and trace conventions, simple and non-simple controls, an exact nonzero seven-form identity and the isolated term in the last proof step. The helper uses standard-library rational exterior algebra and imports no Clifford or curvature engine. These controls support signs and conventions; the universal conclusion comes from the rank/support proof above.

The explicit seven-dimensional comparison form has zero hook with F nonzero, while its embedding in dimension fourteen has nonzero hook. This supports the stated dimension boundary; it is not a classification in other dimensions. Complex forms, indefinite metrics, smooth global choices of simple factors, integrability and PDE sufficiency are not promoted by AA7. Other kernel branches and coincident kernels retain their separate assumptions and open questions.
