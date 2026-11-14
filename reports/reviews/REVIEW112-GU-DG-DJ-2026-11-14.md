# Review112. H4 sign covariance and hook bounds

- The independent full DG run regenerated all 338 integer Clifford
  evaluations in 459 seconds. The contributor run recomputed the
  same evaluations in 625 seconds. Each run checked the 110 source
  columns, both hooks, full curvature, all four permutation sign
  generators, all 227 original Clifford receiver controls, and the
  contraction and quartic inequalities.
- A separate full 128-spinor matrix implementation compared all 227
  coordinates on anchor `f0` and anchor-plus-basis probe `p281`
  modulo 65537. It used the ambient reduced expression on dyadic
  proxies without a parity split, limb expansion or Walsh transform.
  Both probes agreed with the archived integers.
- The DG rank-six hook passed a direct nested Clifford calculation in
  each of four independent radical channels. All grade-one and
  grade-thirteen coefficients have their stated values. Grades five
  and nine retain nonzero coefficients.
- DJ checked all 14,014 H4 and 42,042 H6 incidence slots, exact
  Weyl coefficients, and the receiver weight bounds. Four independent
  rational hook controls checked both hook sums, scalar commutator
  trace, the full grade-two Ricci contraction, and the Jacobi norm.

The covariance counts three independent tensor values in its support.
Its minimum-rank statement is confined to that support. The H6 partial
rank example imposes the hook, scalarity and top equations. It leaves
the coefficient-curvature test and both intermediate grades outside
its scope. External peer review remains outstanding.
