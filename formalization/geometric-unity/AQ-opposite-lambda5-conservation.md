# AQ. Conservation on the opposite moving Lambda5 conic

For the single opposite $W_A$ conic, the grade-one and symmetric Gamma equations reduce to a metric equation for the odd scalar. Contracted Bianchi forces local isoparametricity where its gradient is nonzero. These identities, the [AR](AR-opposite-lambda5-geometry.md) Euler argument and [AS](AS-opposite-lambda5-system.md)'s full form equations now support [AW](AW-lambda5-conic-classification.md)'s completed actual-domain exclusion.

The [helper](../../computation/geometric-unity/AQ_opposite_lambda5_conservation.py) and [recorded output](../../computation/geometric-unity/AQ_opposite_lambda5_conservation.out) check the coefficients. Record: [AQ metadata](../../metrics/geometric-unity/AQ-opposite-lambda5-conservation.yaml). [Independent internal Review47](../../reports/reviews/REVIEW47-GU-AQ-2026-11-03.md) accepts the proof and scope; external peer review remains outstanding.

## Setting

Keep a smooth positive four-base and reconstructed positive Z1 with constant $c_H,c_V>0$ and $\lambda<1/4$; observer signature $(14,0)$ and compact $\mathfrak u(128)$; the true dualized bosonic equation; fixed parallel real $W_A$; real fields; constant $\kappa\ne0$; and $\det M\ne0$. Stay in the single $\Lambda^5:\mathsf w,\mathsf c\Gamma$ kernel, away from coincident non-vector factors. Write
$$T_o=i\kappa(a\Phi^1+w\Phi^1\Gamma),\qquad
w^2=f(a)=\frac{(1+\frac{704}{3}a)(1+160a)}{1280},\qquad w\ne0.$$
The full even field is
$$T_e=v+k,\quad v=P_2(C)+P_{14}(B_0),\quad
C=\frac{24}{\delta}da=d\chi,\quad B_0=208wC,\quad
\delta=1+224a-3328f(a),$$
$$k=r\mathsf w(\omega)+s\mathsf c(\omega)\Gamma,\qquad
(1+\tfrac{704}{3}a)r=\tfrac{80}{3}ws,\quad (1+160a)s=48wr.$$
Here $\omega$ is a smooth real five-form, $\Gamma^2=-1$, and both smooth real chart weights are nonzero. Arbitrary common chart scales can be absorbed into $\omega$, with all derivatives retained. Since $w\ne0$, $dw=f'(a)da/(2w)$; therefore $C$ and $B_0$ are closed. The local potential $\chi$ exists because $24/\delta$ depends only on $a$.

All tensors and derivatives below use the observer metric $g$. Put
$$n=|\omega|^2,\quad S_{ab}=\langle\iota_a\omega,\iota_b\omega\rangle,
\quad c=|C|^2,\quad h=\operatorname{div}C=\Delta\chi,$$
and let $G$ and $\mathrm{Ric}$ be its actual Einstein and Ricci tensors. Metric identification is understood in $\mathrm{Ric}(C)$ and $S(d\ell)$. Thus $\operatorname{tr}S=5n$ and $\Delta=\operatorname{div}\operatorname{grad}$.

## AQ1. Full grade-one stress and metric reduction

Every full solution satisfies
$$0=G-24(\nabla C-hg)+224C\otimes C
+[b+88c+qn]g+\ell S,\tag{AQ1.1}$$
where
$$b=\kappa^2F(a),\quad F=a+312a^2+104f,
\qquad q=\tfrac{112}{3}r^2+40s^2,\quad
\ell=\tfrac{320}{3}r^2-32s^2.$$

*Proof.* The invariant baseline contributes $G+bg$. J2 gives the complete derivative contribution $-24(\nabla C-hg)$ and vector quadratic $224C\otimes C+88cg$. In particular $D(v)$ is retained. The grade-one parts of $D(k)$ and $2\mathcal B(v,k)$ vanish: their input $V\otimes\Lambda^5$ has none of the scalar, alternating two-tensor or symmetric tracefree two-tensor summands of $V\otimes V$. This includes every form and weight derivative.

The grade-one part of $Q(k)/i$ is $qn g+\ell S$. The only metric contractions of two identical five-forms into a two-tensor are $ng$ and $S$, both symmetric; twelve tensor slots are below fourteen, excluding an epsilon contribution. A complete unit-simple calibration with symbolic $r,s$ fixes both coefficients, including the absence of a cross term. This argument concerns grade one, not the higher quadratic channels. ∎

AH1's symmetric right-Gamma equation, also retaining derivatives and mixed terms in their proper channels, gives
$$rsn=\frac{w}{240}[14\kappa^2(1+208a)-43264c],\qquad
S=\alpha g-\beta C\otimes C,\quad
\beta=\frac{624w}{19rs},\quad\alpha=\frac{5n+\beta c}{14}.\tag{AQ1.2}$$
Substitution into (AQ1.1) yields
$$\boxed{G-24(\nabla C-hg)+[\kappa^2U(a)+V(a)c]g
+W(a)C\otimes C=0,}\tag{AQ1.3}$$
with
$$U=\frac{(208a+1)(21840a+113)}{480},\qquad
V=-\frac{64(402896a+1719)}{285},\qquad
W=\frac{1024(53-676a)}{285}.$$
These coefficients are chart-independent: the kernel equations give $wr/s=(1+160a)/48$ and $ws/r=3(1+704a/3)/80$.

## AQ2. Contracted Bianchi compatibility

Before eliminating the form, the complete divergence of (AQ1.1) is
$$0=q\,dn+\ell\operatorname{div}S+n\,dq+S(d\ell)+db+200\,dc
+224hC-24\mathrm{Ric}(C).\tag{AQ2.1}$$
In particular no variable-weight derivative is suppressed. For arbitrary $J_a=\nabla_a\omega$,
$$ (\operatorname{div}S)_b
=\langle\operatorname{Ctr}J,\iota_b\omega\rangle
+\tfrac12(dn)_b-\langle\iota_b\operatorname{Alt}J,\omega\rangle,$$
where $\operatorname{Alt}J=d\omega$ and $\operatorname{Ctr}J=-\delta_\nabla\omega$. This is the covariant product rule, valid for the full jet, including its Cartan part. No AN constant-endpoint cancellation is imported.

Put $T=\kappa^2U+Vc$. Tracing (AQ1.3) gives $\operatorname{tr}G=-312h-14T-Wc$. Since $\mathrm{Ric}=G-(\operatorname{tr}G)g/12$ and $\nabla_C C=dc/2$,
$$\mathrm{Ric}(C)=12\,dc+[2h+T/6-11Wc/12]C.\tag{AQ2.2}$$
Using $\operatorname{div}(\nabla C-hg)=\mathrm{Ric}(C)$, contracted Bianchi becomes exactly
$$\boxed{\begin{aligned}
0={}&(V+W/2-288)\,dc\\
&+\left[(W-48)h+\kappa^2(\delta U'/24-4U)
+c\{\delta(V'+W')/24+22W-4V\}\right]C.
\end{aligned}}\tag{AQ2.3}$$
Primes denote derivatives in $a$. This identity keeps the Hessian terms through $h$ and $dc$; it does not force $n$ constant.

## AQ3. Local isoparametricity and Ricci alignment

On every sufficiently small solution patch with $da\ne0$, there are smooth one-variable functions $\rho,\sigma$ such that
$$|da|^2=\rho(a),\qquad\Delta a=\sigma(a).\tag{AQ3}$$
Also $\mathrm{Ric}(\operatorname{grad}a)=\zeta(a)\operatorname{grad}a$ for a local smooth function $\zeta$.

*Proof.* Wedge (AQ2.3) with $da$. Up to sign, its coefficient of $dc$ is
$$A_c=288-W/2-V=\frac{32(816608a+5155)}{285}.$$
Its one zero, $a_*=-5155/816608$, lies on the real left branch. It cannot be divided away. On a regular patch $a$ is a coordinate; the complement of this one level is dense, so $dc\wedge da=0$ extends to it by continuity. Hence $c=c(a)$ locally.

Where $W\ne48$, (AQ2.3) now makes $h$ a function only of $a$. The coefficient
$$W-48=\frac{16(2537-43264a)}{285}$$
has just one zero. Its complement is again dense on the regular patch. Consequently $dh\wedge da=0$ extends across that level, so $h=h(a)$ locally there too. No division at either exceptional level is needed for the conclusion.

With $H(a)=24/\delta$, the identities $c=H^2|da|^2$ and $h=H\Delta a+H'|da|^2$ prove (AQ3). Equation (AQ2.2) then proves the Ricci alignment. ∎

Here “locally isoparametric” means precisely the two scalar identities (AQ3), following the definition in [Qi-Ming Wang, *Isoparametric Functions on Riemannian Manifolds. I*, Math. Ann. 277 (1987), 639–646, first page](https://gdz.sub.uni-goettingen.de/fulltext/PPN235181684_0277/00000645?_format=xml). No global foliation, constant principal curvatures or space-form classification is inferred. AW uses these necessary identities in its exclusion of the entire declared single-conic sector; the local derivation here is retained.

## Reproduction and limits

Run `python -X utf8 -I computation/geometric-unity/AQ_opposite_lambda5_conservation.py`. The recorded run passes 23 exact controls in 1.724 seconds with Python 3.12 and SymPy 1.14.0: thirteen stress/derivative controls and ten scalar elimination/Bianchi controls, restricted to the opposite conic. The invariant-contraction and covariant arguments supply universality; the density argument is a hand proof. The helper shares T's accepted Clifford engine and AH's form operations, and does not reconstruct a metric or prove PDE existence.

The result concerns only the stated nonzero-$w$ single opposite conic. It does not cover coincident kernels, $\det M=0$, other operators, indefinite real forms or arbitrary non-invariant odd fields. No field, inequality or spectral shape is continued along a fibre.
