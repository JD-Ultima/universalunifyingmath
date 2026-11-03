# AW. The two single Lambda5 conics on positive Z1

The opposite single $W_A$ Lambda5 conic has no smooth full bosonic solution on a nonempty positive-Z1 observer domain. A normal component of the full AS compatibility equation, paired with the actual form, forces a rational norm profile incompatible with AR's Euler/conservation profile. Together with [AP](AP-first-lambda5-obstruction.md), this excludes both single Lambda5 conics under the stated assumptions. It does not classify zero-order form data or coincident kernels.

Code: [AW helper](../../computation/geometric-unity/AW_lambda5_conic_classification.py) and [output](../../computation/geometric-unity/AW_lambda5_conic_classification.out). Record: [AW metadata](../../metrics/geometric-unity/AW-lambda5-conic-classification.yaml). [Review53](../../reports/reviews/REVIEW53-GU-AW-2026-11-03.md) accepts the full proof and [Review54](../../reports/reviews/REVIEW54-GU-AW-2026-11-03.md) records the phase and scope audit, disclosing its contributor's role in the convention controls. Both reproduce all 34 controls. These are internal AI-assisted reviews; external peer review remains outstanding.

## Setting

Retain [AS/AQ](AS-opposite-lambda5-system.md)'s true reconstructed dualized equation, smooth positive four-base and positive Z1 with constant $c_H,c_V>0$ and $\lambda<1/4$, compact real $(14,0)$ and $\mathfrak u(128)$, fixed parallel real $W_A$, real fields, constant $\kappa\ne0$ and $\det M\ne0$. On an entire connected open domain, stay in the single opposite $\Lambda^5:\mathsf w,\mathsf c\Gamma$ kernel away from every other non-vector singular factor:
$$T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\quad
w^2=f(a)=\frac{(1+704a/3)(1+160a)}{1280},$$
$$k=r\mathsf w(\omega)+s\mathsf c(\omega)\Gamma,\quad
(1+704a/3)r=\tfrac{80}{3}ws,\quad(1+160a)s=48wr.$$
The full even field includes the slaved vector part
$$T_e=P_2(C)+P_{14}(208wC)+k,\qquad
C=\frac{24}{\delta}da,\quad\delta=1+224a-3328f(a).$$
All derivatives of the odd fields, form and chart weights remain included. [AR1](AR-opposite-lambda5-geometry.md) proves $w\ne0$ everywhere and that $da\ne0$ is open dense. It is therefore enough to exclude a regular patch. On it set
$$\eta=r\omega,\quad\tau=s/r,\quad x=\tau^2>0,\quad c=|C|^2>0.$$
This canonical change retains $\nabla\eta=r\nabla\omega+dr\otimes\omega$. The real chart gives
$$b(a):=\frac w\tau=\frac{1+160a}{48},\qquad
x(a)=\frac{9(1+704a/3)}{5(1+160a)}.$$

## AW1. A normal moment identity

Choose an oriented orthonormal frame with $C=\rho e^0$, $\rho=\sqrt c>0$, and write
$$\eta=e^0\wedge\phi+\psi,\quad
\mathsf U=\phi_4,\quad\mathsf V=\psi_5J,\quad
J=\Gamma_{13},\quad\Gamma=\Gamma_{14}=e^0J.$$
Here $J$ is central in the tangent Clifford algebra, $J^2=1$, and it anticommutes with $e^0$. Put $\Omega=e^0\mathsf U+\psi$. Use scalar Clifford coefficient $\operatorname{tr}=\operatorname{sc}$ in the ambient fourteen-dimensional algebra, equivalently its normalized cyclic matrix trace. Grade orthogonality gives
$$\operatorname{tr}\mathsf U=\operatorname{tr}\mathsf V=
\operatorname{tr}(\mathsf U\mathsf V)=0,\qquad
m=\operatorname{tr}\mathsf U^2=|\phi|^2,\quad
\ell=\operatorname{tr}\mathsf V^2=|\psi|^2.$$
No single-chirality odd-dimensional trace is substituted for this ambient trace.

Use AS3's canonical coefficients
$$\begin{gathered}
\theta=\frac{3328}{3}b(x+4)-\frac{7\delta}{360f},\qquad
\alpha=-\frac{3(304+9\theta)}{16(61+16x)},\quad
\zeta=\frac{3\tau(-112+5\theta)}{16(64+19x)},\\
A=-\tfrac43\tau(41+11x),\quad B=\tfrac83(38+11x),\quad
E=\tfrac83\tau(4+x),\quad F=\tfrac83(4+x),\\
p=E\alpha+16\tau,\qquad q=32\tau-F\zeta.
\end{gathered}$$
AS3 gives $F_4=\alpha\rho\mathsf U$ and $F_8=\zeta\rho\mathsf V$: the latter uses $e^0\psi\Gamma=-\mathsf V$. Thus
$$\Omega^2=n+\alpha\rho\mathsf U+\zeta\rho\mathsf V,\qquad n=|\eta|^2.$$
Let $M_{ij}=\operatorname{tr}(\mathsf U^i\mathsf V^j)$. Since
$\Omega^2=\mathsf U^2+\psi^2+e^0[\mathsf U,\psi]$ and $\psi^2=\mathsf V^2$, pairing this square with $\mathsf U,\mathsf V$ gives
$$M_{30}+M_{12}=\alpha\rho m,\qquad
M_{21}+M_{03}=\zeta\rho\ell.\tag{AW1a}$$
The normal term has zero scalar trace in these pairings. No commutation or simultaneous diagonalization is needed.

The full normal AS compatibility, with all Gamma products on the right, is
$$A Y_{5,0}+B(Y_{9,0}\Gamma)+\rho p\,e^0\mathsf U-\rho q\,\psi=0.\tag{AW1b}$$
Indeed $X_0=\rho(e^0\mathsf U+\psi)$, $Y_0=\rho\psi$, $Z_0=\rho e^0\mathsf U$ and $K(X)_0=0$. The square traces give $E_{5,0}=\alpha\rho e^0\mathsf U$ and $E_{9,0}\Gamma=\zeta\rho\psi$.

Direct multiplication retains the essential anticommutator:
$$\Omega e^0\Omega=e^0(\mathsf U^2-\psi^2)+\{\mathsf U,\psi\}.$$
Reversion fixes a five-form, so its real form pairing is $\operatorname{sc}(XY)$. Ambient cyclicity and grade projection give the four pairings

| Pattern | Paired with $e^0\mathsf U$ | Paired with $\psi$ |
|---|---|---|
| $Y_{5,0}$ | $M_{30}-M_{12}$ | $2M_{12}$ |
| $Y_{9,0}\Gamma$ | $-2M_{21}$ | $M_{21}-M_{03}$ |

For the negative factor two, $\Gamma e^0\mathsf U=-J\mathsf U$; cyclicity combines both anticommutator terms. Terms with one remaining normal $e^0$ have zero scalar trace. Also $\langle e^0\mathsf U,\psi\rangle=0$. Pairing (AW1b), then eliminating all third moments with (AW1a), proves
$$\boxed{(A\alpha+p)m=(B\zeta+q)\ell.}\tag{AW1}$$
This is a universal normal trace identity. It uses neither a bound on joint spectral values nor an example asserted to realize AS.

## AW2. The forced rational norm

The actual canonical Gram equation and its trace give
$$\beta=\frac{624b}{19},\quad
n=\kappa^2\nu_0+\nu_1c,\quad
\nu_0=\frac7{120}b(1+208a),\quad\nu_1=-\frac{2704}{15}b,$$
$$m=\frac{5n-13\beta c}{14},\qquad
\ell=\frac{9n+13\beta c}{14}.$$
Define the rational functions
$$\begin{gathered}
u_0=5\nu_0/14,\quad u_1=(5\nu_1-13\beta)/14,\qquad
v_0=9\nu_0/14,\quad v_1=(9\nu_1+13\beta)/14,\\
r_1=(A\alpha+p)/\tau=-4(11+3x)\alpha+16,\qquad
r_2=(B\zeta+q)/\tau=\tfrac{16}{3}(17+5x)(\zeta/\tau)+32,\\
L_*(a)=r_1u_1-r_2v_1,\qquad K_*(a)=r_1u_0-r_2v_0.
\end{gathered}$$
Equation (AW1) is exactly
$$L_*c+K_*\kappa^2=0.\tag{AW2a}$$
The coefficient $L_*$ is not identically zero: its exact rational value at zero is $-4672630976/37698489$. Its possible zero levels are finite after clearing denominators. On a regular patch $a$ is a local coordinate, so a smaller nonempty patch avoids those levels and all rational denominators. There
$$c=\kappa^2\mathcal F(a),\qquad
\mathcal F=-K_*/L_*.\tag{AW2}$$
No exceptional level is silently divided away, and no assumption about realizability of a prescribed Gram is used.

## AW3. No opposite-conic domain

Write the scalar coefficients from AQ/AR as
$$\begin{gathered}
U=\frac{(208a+1)(21840a+113)}{480},\quad
V=-\frac{64(402896a+1719)}{285},\quad
W=\frac{1024(53-676a)}{285},\quad D=V-11W/2,\\
P=V'+W'-2WD/\delta,\quad Q=U'-2WU/\delta,\quad
R_0=-12(W-48)/\delta.
\end{gathered}$$
On the actual regular patch, AR5b is
$$Dy'+Py+tQ+R_0=0,\quad
y=c_0c,\quad t=\kappa^2c_0>0,\quad c_0=4c_V(1-4\lambda)>0.$$
Substitute (AW2). Since $y=t\mathcal F(a)$, this requires
$$t\mathcal H(a)+R_0(a)=0,\qquad
\mathcal H=D\mathcal F'+P\mathcal F+Q.\tag{AW3a}$$
These rational functions would agree on an open interval of scalar values. Clearing denominators makes this a polynomial identity. The following small exact certificates, derived from the displayed rational definitions, refute it:

| Rational quantity | Exact value |
|---|---|
| $\mathcal F(a)/a$ at infinity | $-35575543/31723328$ |
| $\mathcal H(a)/a$ at infinity | $3138842736144/15324925\ne0$ |
| $R_0(a)/a$ at infinity | $0$ |
| $\mathcal H(0)$ | $35018717715292549351/51860784957622875>0$ |
| $R_0(0)$ | $20296/19>0$ |

The leading rational coefficient would force $t=0$. Independently, the rational identity evaluated at zero would give a sum of two positive terms for $t>0$. Both are tests of a rational identity; neither continues an actual field, its domain or positivity to zero, infinity or across excluded levels.

Thus no regular opposite-conic solution patch exists. AR1 supplies a nonempty regular patch in every proposed nonempty entire conic domain and excludes moving zero-$w$ contacts. Consequently the entire opposite single-conic domain is empty. This argument uses AR1 and AR5, not AR3's base-rank classification or the earlier geometric matrix-rank certificates.

## AW4. Both single Lambda5 conics

Under the standing positive-Z1, compact-real and true-equation hypotheses, a full solution cannot remain on an entire connected open domain in either single $W_A$ Lambda5 conic away from coincident factors. The opposite pair follows from AW3; the first $c/w\Gamma$ pair is AP3. All endpoints are included in this domain statement through AR1 and AP3's endpoint argument.

This excludes actual smooth domains in the declared sectors. It does not prove that the remaining cubic tensor in AS/AV follows algebraically from the zero-order equations, nor decide real zero-order admissible feasibility. Coincident kernels, $\det M=0$, other operators, non-invariant odd fields and indefinite real forms remain outside the theorem. No physical or literature-first claim is made.

## AW5. Pointwise exclusion of an isolated Lambda5 addition

Let $T$ be a smooth full solution on an open observer domain under the standing background, true-equation, fixed-$W_A$, invariant-odd and $\det M\ne0$ hypotheses. The odd parameters need not stay globally on a conic. At a point away from all coincident non-vector singular factors, including the other Lambda5 factor, its even field has no nonzero extra Lambda5 component.

*Proof.* If its smooth Lambda5 isotypic projection were nonzero at that point, it would remain nonzero on a smaller neighbourhood. [R2 and R4](R-even-operator.md) give no Lambda5 source in the even equation, so its homogeneous copy block annihilates that projection. The other factors remain invertible after shrinking. Thus one fixed Lambda5 determinant factor vanishes throughout the neighbourhood and every other non-vector component is zero. Its two-by-two copy block cannot be the zero matrix when $\kappa\ne0$, so its rank-one kernel has local smooth charts, including at the endpoints. These put the full field in exactly one of the two single-conic ansatzes above. On a connected smaller patch, AP3 or AW3 is a contradiction. This does not exclude an odd-field conic contact or crossing whose Lambda5 component is zero. ∎

## Evidence

Run `python -X utf8 -I computation/geometric-unity/AW_lambda5_conic_classification.py` from the repository root. The recorded author run passes 34 exact controls in 0.544 seconds with Python 3.12 and SymPy 1.14.0; Reviews53–54 reproduce them in 0.538 and 0.531 seconds. Twelve checks calibrate the full normal Clifford pairings on three declared noncommuting inputs; these are convention controls, not asserted AS data. Twenty-two scalar checks verify the moment elimination, actual Gram coefficients, forced rational norm, direct AQ-to-AR profile reduction and both contradiction certificates. Cyclicity, grade orthogonality, the actual Gram reduction and the open-interval rational-identity argument supply the universal proof. The helper imports no earlier check or metric engine. Earlier AS/AQ/AR inputs retain their recorded evidence; an independent reproduction of every imported calculation is not claimed.
