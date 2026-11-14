# DB. Bounded analytic H6 lifts from compact averages

Every value in the 286-dimensional real (5,5) module has a
globally bounded real analytic closed-hook lift which vanishes at
infinity. Consequently every value in K6820 has a globally bounded
real analytic closed-hook lift, using [BR](BR-h6-closed-field-lifts.md) for the other 6,534 directions.
The converse value inclusion follows from BU, already for local C2
fields. This theorem is entirely linear.

Let B be the normalized real CU field, b=B(o), N0=||b||, and K=SO5
with normalized Haar measure. Thus N0^2=3/131072. Write B_k=k.B and
a_k=k e5 for its center. DA proves the pointwise majorant
 ||B_k(Y)|| <= N0 g_Y(a_k),
 g_Y(a)=(a^T Y^(-1)a)^(-5/4), detY=1.
Here det(Y restricted to a^perp)=a^T Y^(-1)a is the cofactor identity.
## Exact endpoint integral and uniform decay
For any real positive definite5x5 matrix M, polar coordinates in the
Gaussian integral give
 integral_S4 (a^T M a)^(-5/2) d sigma(a) = det(M)^(-1/2),
where sigma is normalized surface measure. Indeed integration in the
radial variable gives Gamma(5/2)/2 times the sphere integral, and the
same Gaussian equals pi^(5/2)/sqrt(detM); normalize at M=I.
Taking M=Y^(-1) gives ||g_Y||_(L2(S4))=1. Consequently
 integral_K ||B_k(Y)||^2 dk <= N0^2 for every Y.

These kernels tend to zero in L1 uniformly as Y escapes to infinity.
For a quantitative estimate let Lambda=lambda_max(Y^(-1))>=1 and
let u be a unit eigenvector of Lambda. On |a.u|>delta,
g_Y(a)<=Lambda^(-5/4)delta^(-5/2). The normalized S4 coordinate
density is (3/4)(1-t^2), so sigma{|a.u|<=delta}<=3delta/2.
Cauchy-Schwarz and the exact L2 identity therefore give
 ||g_Y||_1 <= sqrt(3delta/2)+Lambda^(-5/4)delta^(-5/2).
Choosing delta=Lambda^(-5/12) yields
 ||g_Y||_1 <= (1+sqrt(3/2))Lambda^(-5/24).
Bounded Lambda defines a compact subset of positive determinant-one
matrices: all eigenvalues of Y lie between Lambda^(-1) and Lambda^4.
Thus the decay is uniform at infinity, not just along a selected ray.
## L2 densities and the center-marginal extension
For h in L2(K), define H_h(Y)=integral_K h(k) B_k(Y) dk (real h),
or take the real part for complex h. Cauchy-Schwarz gives
 ||H_h||_infinity <= N0 ||h||_2.
For bounded h the L1 decay estimate shows H_h vanishes uniformly at
infinity. General L2 h can be approximated by bounded functions in L2;
the uniform norm bound passes this C0 property to the limit.

More generally let mu be a finite signed or complex measure on K.
Assume the pushforward of its total variation to the center sphere
has density w in L2(S4). Disintegration is not needed: directly
 ||integral B_k dmu|| <= N0 integral_S4 w(a)g_Y(a)d sigma(a)
 <= N0||w||_2.
Bounded truncations of w and the same uniform L2 bound prove C0.
The measure may be singular in the fibers over each center. The
condition concerns the TOTAL-VARIATION marginal, so signed fiber
cancellation is not silently used to justify an estimate.
For an L2 Haar density h, Jensen's inequality for conditional
expectation gives ||w||_2<=||h||_2, recovering the first statement.

Every such integral is globally real analytic, closed, hooked, and
pure-(5,5). To justify analytic regularity, the integrand is jointly
real analytic on K times Y. Near any fixed Y0, compactness of K gives
a common complex neighborhood and a uniform holomorphic bound in
local tensor coordinates, by finitely many analytic charts. The
finite total variation of the measure permits integration and
differentiation there. Thus the integral is locally holomorphic and
its real restriction is analytic. Linear differential and algebraic
equations pass under the integral. This is stronger than a merely
smooth dominated-convergence argument.
## Explicit lift through every (5,5) value
For v in the real irreducible module V55, define
 F_v(Y)=(286/N0^2) integral_K <k b,v> B_k(Y) dk.
The operator integral (k b) tensor (k b) dk commutes with SO5 and
is symmetric. Irreducibility makes it scalar, and its trace N0^2
gives scalar N0^2/286. Hence F_v(o)=v exactly.
Also integral_K <k b,v>^2 dk=N0^2||v||^2/286, so
 ||F_v||_infinity <= sqrt(286)||v||.
The density is bounded and analytic on K, therefore F_v is C0 with
the explicit bound
 ||F_v(Y)|| <= 286(1+sqrt(3/2))||v|| Lambda(Y)^(-5/24).
No optimality of either bound is claimed.
## Nonlinear scope of mixed families
Let G be a BW bounded analytic compact-measure AN field. Let H be
any of the CU averages above. Also allow finitely many compact
rotated CU fields multiplied by arbitrary entire functions of their
own transverse coordinates. Suppose the resulting global analytic
field satisfies the full nonlinear equations on a nonempty open set.
Real analytic continuation extends those polynomial tensor identities
to all connected Y; in particular the prescribed positive constant
norm gives global boundedness. DA, with bounded background G+H,
forces every finite same-center CU group to vanish identically.
The remaining field is G+H, with H in C0.

Along the generic radial ray used in BW, the transported G has a
limit in the single rotated AN value family to which the BV positive
receiver applies. Since H vanishes in norm, the transported sum has
the same limit. The nonlinear equations and positive prescribed norm
pass to this limit. BV excludes this nonzero limiting AN value,
giving a contradiction. Thus this entire mixed global ansatz is
excluded even as a solution on an open set.

This does not exclude arbitrary bounded pure-(5,5) fields: DB actually
constructs many. It also does not represent arbitrary analytic fields,
cover all singular center measures, or rule out all nonlinear fields.

## Reproduction

    python -X utf8 computation/geometric-unity/DB_h6_bounded_average.py

[Review107](../../reports/reviews/REVIEW107-GU-DA-DB-2026-11-14.md)
records the independent checks. External peer review remains outstanding.
