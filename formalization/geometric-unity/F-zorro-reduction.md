# F. The Zorro metric on $Y=\mathrm{Met}(X)$ and the 14→4 reduction

Claims: GU-006, GU-021, GU-024, GU-054. Gaps: GU-G02, GU-G05, GU-G06, GU-G14. Code, all exact unless stated:

- `computation/geometric-unity/F_jets.py` (library);
- `F_fiber_curvature.py` → `.out` (FC1–FC4);
- `F_zorro_reduction.py` → `.out` (R0–R5);
- `F_zorro_einstein.py` → `.out` (E1–E4).

Metric records: `metrics/geometric-unity/F-fiber-geometry.yaml`, `F-zorro-reduction.yaml`. Conventions follow A (signature $(a,b)$ = ($a$ positive, $b$ negative)). Notation: $n=\dim X$, $N=n(n+1)/2$, $d=n+N=n(n+3)/2$ (= 14 for $n=4$). Curvature: $R_{abcd}=g(R(\partial_c,\partial_d)\partial_b,\partial_a)$, $\mathrm{Ric}_{bd}=g^{ac}R_{abcd}$; the unit sphere has $\mathrm{Ric}=(n-1)g$ (control R0).

## F.0 Source assertions (layer a)

All quotations are from [W21]; page = PDF page = printed page.

| # | Location | Text (verbatim, abridged with …) |
|---|---|---|
| a1 | Def. 3.1, (3.1), p. 16 | "local Riemannian embeddings … into a Riemannian manifold $Y^d$ … inducing a pull back metric $g_X=\iota^*(g_Y)$ and defining a normal bundle $N_\iota^{d-n}\subset T(Y)$ with metric" |
| a2 | p. 17 | "Note that in this case, Y doesn't usually have a pre-existing metric. Instead the choice of a section ι of the metric bundle induces a metric on Y in the portion of Y that lies above U ⊂ X." |
| a3 | (3.7), p. 18 | "$H^*=\pi^*(T^*(X))$ … which at the point $g\in Y$ carries an induced metric $g$" |
| a4 | (3.8)–(3.9), p. 19 | Frobenius form $\langle A,B\rangle_y=\mathrm{Tr}_y(A\cdot B)$ on $V$; "a freedom to assign either a (0,1) or (1,0) signature to the trace component" |
| a5 | p. 20 | "the fibers at a point g ∈ Y impart metrics to … the vertical $V_g$ and horizontal sub-bundles $H^*_g$, taking those sub-spaces to be orthogonal" |
| a6 | (3.15), p. 21 | "When a metric is chosen below on X compatibly with V and H, the Chimeric bundle and Tangent are given a natural isomorphism via the Levi-Civita connection" |
| a7 | §3.4, (3.16)–(3.17), p. 21 | "what is necessary to split the cyclic long exact sequence … is in fact a connection on the space Y viewed as a bundle over X … a metric choice below on X leads ultimately to the choice of a connection above over Y"; diagram $\gimel\to\aleph_\gimel$, $g_\aleph\to A_g$; "each observation of Y via a choice of metric ג on X actually induces a metric and connection on Y" |
| a8 | p. 22 | "The subset of metrics that we are considering … are equivalent to the space of connections that can arrise on X as Levi-Civita connections"; "the issue of what upstairs gravity waves are possible doesn't arise" |
| a9 | §3.5, p. 22 | "we must combine these choices of fiber metrics in Y with metrics pulled up from the base space of signature (1, 3)" |
| a10 | (3.39), (3.42), pp. 25–26 | $\iota^*(TY)=TU\oplus N_\iota$; "(Upstairs) Metric Via the Zorro diagram"; $V=V^*=N_\iota$, $H=H^*=TX$ |
| a11 | (4.5), p. 28 | "$D\gimel:TX\to TY$ splits the tangent bundle on Y along ג(X) ⊂ Y into the image of Dג and its 10 dimensional complement $V^{10}$ … $\mathrm{Im}(D\gimel)=(H^4_{\gimel_x})^*$" |
| a12 | §5.1.1, p. 31 | "a derived field $\aleph=\aleph_\gimel$ representing the Levi-Civita connection across all bundles on which it is induced from the metric ג" |
| a13 | §6, p. 36 | "we begin with a metric ג, which determines a connection Aleph $\aleph_\gimel$ on the bundle TX. At that point, the connection ℵ determines a metric $g_\aleph$ on TY given the metric data already present on the horizontal and vertical tangent sub-bundles HY and VY. This, in turn, determines a Levi-Civita connection $\nabla^{g_\aleph}=\nabla^0$" |
| a14 | (9.1), (9.4), pp. 43–44 | $I^B_1:\mathcal G\times\mathrm{MET}(X^{1,3})\to\mathbb R$; the pairing in (9.4) is labelled "Zorro Metric" |
| a15 | (9.9)–(9.10), p. 45 | $\Upsilon_\omega=S_\omega-T_\omega=0$ annotated "$S_\omega\leftrightarrow R_{\mu\nu}-\frac s2g_{\mu\nu}$", "$T_\omega\leftrightarrow\Lambda g_{\mu\nu}$" |
| a16 | §11, p. 49 | "To interpret differential equations governing fields native to Y back on X means pulling back not only bundles but so-called Jets or Sprays of sections on Y." |
| a17 | §12.3, p. 56 | "our subset of 'metrics' under consideration is really equivalent to the space of connections that split the long repeating exact sequence … the Zorro construction … uses a connection to derive the highly restricted metrics on Y" |

**Not specified by the draft:**

- the horizontal metric at a point $y\neq\gimel(x)$: a3/a13 give the point $y$ itself, a9 gives the base metric;
- $\lambda$ (a4 "assigns" a sign);
- the scale $c_V/c_H$;
- whether $\Upsilon=0$ is imposed on all of $Y$ or only along $\gimel(X)$;
- how a 13-form equation on $Y^{14}$ is read as the 2-tensor equation (9.10) on $X^4$.

No derivation of (9.10) is given.

## F.1 Reconstruction (layer b)

**R-Z0.** $X^n$ carries a metric $g=\gimel$ of signature $(p,q)$, and $Y=\mathrm{Met}_{(p,q)}(X)$ has coordinates $(x^\mu,y_{ij})$ with $i\le j$. The section is $s=g$, and $\aleph$ is the Levi-Civita connection of $g$.

**R-Z1 (horizontal distribution).** $\aleph$ induces a linear connection on $\mathrm{Sym}^2T^*X$. Its connection map gives
$$Dy=dy-(\Gamma_\mu^{T}y+y\Gamma_\mu)\,dx^\mu,\qquad (\Gamma_\mu)^k{}_i=\Gamma^k_{\mu i},$$
and $H^\aleph=\ker D$. By Lemma F0 there is no metric-free alternative.

**R-Z2 (metric).** Two readings:

- **Z1 (primary):**
  $$g_\aleph=c_H\,y_{\mu\nu}\,dx^\mu dx^\nu+c_V\,G_\lambda|_y(Dy,Dy),\qquad G_\lambda|_y(h,k)=\mathrm{tr}(y^{-1}hy^{-1}k)-\lambda\,\mathrm{tr}(y^{-1}h)\mathrm{tr}(y^{-1}k).$$
  Equivalently $g_\aleph=\Theta_\aleph^*g_C$. Here $\Theta_\aleph=(D,\ \flat_y\circ\pi_*):TY\to C=V\oplus H^*$, and $g_C=c_VG_\lambda\oplus c_Hy^{-1}$ is the chimeric metric R-A3 of A. So Z1 is literally a3 + a5 + a6 + a13: the chimeric metric transported by the Levi-Civita splitting. It depends on $g$ only through $\aleph$, which is the property claimed in a8/a17.
- **Z2 (alternative):** the horizontal block is $c_H\,g_{\mu\nu}(x)$, the base metric pulled back, following a9. This is the candidate listed under GU-G06 in `gaps/geometric-unity.md`. It is the standard Kaluza–Klein/Sasaki-type choice. It violates a8/a17 because $g$ enters beyond $\aleph$.

Z1 and Z2 agree on $s(X)$ to zeroth order. They differ in their first vertical derivatives, so their curvatures on $s(X)$ differ.

**Parameters.** $\lambda\neq1/n$ and $c_H,c_V\neq0$.

- The draft's $(7,7)$ requires $\lambda>\tfrac14$ and $\operatorname{sign}(c_Vc_H)=+1$ for base $(1,3)$ (Theorem A4).
- The literal reading of a4 ("flip the sign of the trace part, keep its size") is $G(h,h)=|h_0|^2-\tfrac1n(\mathrm{tr}\,y^{-1}h)^2$, i.e. $\lambda_D=2/n$ ($=\tfrac12$ for $n=4$).
- The standard DeWitt value $\lambda=1$ is also tested.
- $c_H=1$ makes $s^*g_\aleph=g$ (a1).

**R-Z3 (how 4D gravity is read off).** Three candidate readings:

- **(P) Projection:** $s^*(G^Y)(X,Y):=G(g_\aleph)(s_*X,s_*Y)$. This is what one gets by Hodge-dualising the 13-form equation to a $T^*Y$-valued 1-form and pulling it back (a16). By Theorem C5, under R-S $\circledast F_{\nabla^0}=-{*}G(g_\aleph)$.
- **(A) Reduced action:** $S[g]=\int_X s^*(s_Y)\,dvol_{s^*g_\aleph}$.
- **(E) Field equation on $Y$ along $s(X)$:** $G(g_\aleph)+\Lambda_Yg_\aleph=0$ on $s(X)$.

The full Einstein–Hilbert action $\int_Y s_Y\,dvol$ is not defined. The fiber has infinite volume (Prop. F1: it contains a line). In Z2, $s_Y$ is invariant under fiber scaling $y\mapsto e^ty$, and in Z1 that scaling only rescales $c_H$. So the integrand does not decay along the fibers. (Made precise 2026-10-30, round 11, `O-gimel-variation.md`: with the fibre-invariant measure every nonzero fibre integrand of the invariant sector diverges (O5); over an indefinite base no natural weight function converges, and the convergent natural weights are measures on the scaling ray through $s(X)$ (O6, O7).)

## F.2 Results (layer c)

Imported lemmas:

- **L-F1** O'Neill's formulas for (pseudo-)Riemannian submersions: B. O'Neill, *Michigan Math. J.* 13 (1966) 459; Besse, *Einstein Manifolds* (1987) 9.28–9.37; the semi-Riemannian case in O'Neill, *Semi-Riemannian Geometry* (1983). Algebraically identical for nondegenerate fibers.
- **L-F2** The Gauss equation for nondegenerate submanifolds: O'Neill (1983) 4.5.
- **L-F3** Walczak's mixed-scalar-curvature formula for complementary orthogonal distributions: P. Walczak, *Colloq. Math.* 58 (1990) 243; V. Rovenski, *Foliations on Riemannian Manifolds and Submanifolds* (1998):
  $$\mathrm{S}_{\rm mix}=\mathrm{div}(H_1+H_2)-|h_1|^2-|h_2|^2+|H_1|^2+|H_2|^2+|T_1|^2+|T_2|^2.$$
  It is a local frame identity, and we use it in indefinite signature. Its output is re-verified exactly in R2/R3.
- **L-F4** Curvature of a symmetric space: $R(X,Y)Z=-[[X,Y],Z]$ on $\mathfrak p$ (Helgason, *Differential Geometry, Lie Groups and Symmetric Spaces* (1978) IV.4). Fixed-point sets of isometries are totally geodesic (Kobayashi, *Transformation Groups in Differential Geometry* (1972) II.5).
- **L-F5** For a connection on a vector bundle $E$, $[X^H,Y^H]=[X,Y]^H-(R(X,Y)e)^{\rm vert}$ (Kobayashi–Nomizu I, II.5/III.1).
- **L-F6** $\delta\int(R-2\Lambda)\,dvol_g=-\int\langle G+\Lambda g,\delta g\rangle\,dvol_g$ (Besse 4.17).

For the DeWitt geometry of $\mathrm{Met}(\mathbb R^n)$ in the positive-definite case, cf. B. DeWitt, *Phys. Rev.* 160 (1967) 1113; D. Freed and D. Groisser, *Michigan Math. J.* 36 (1989) 323; O. Gil-Medrano and P. Michor, *Quart. J. Math.* 42 (1991) 183. Our checks are independent and cover every signature and every $\lambda$.

**Lemma F0 (no natural connection on $\mathrm{Met}(X)$).** No $\mathrm{Diff}(X)$-natural horizontal distribution on $\mathrm{Met}(X)\to X$ exists.

*Proof.* A horizontal $n$-plane at $y\in Y_x$ is a 1-jet $j^1_xh$ of local metrics with $h(x)=y$. Diffeomorphisms $\varphi$ with $\varphi(x)=x$ and $j^1_x\varphi=\mathrm{id}$ act on such 1-jets by $\partial_kh_{ij}\mapsto\partial_kh_{ij}+y_{lj}a^l_{ik}+y_{il}a^l_{jk}$, where $a^l_{ik}=\partial_i\partial_k\varphi^l(x)$ is arbitrary and symmetric in $(i,k)$. The map $a\mapsto(\cdot)$ is injective: put $b_{jik}=y_{jl}a^l_{ik}$; then $b_{jik}+b_{ijk}=0$ with $b$ symmetric in its last two indices forces $b=0$ (the Christoffel argument). Since both sides have dimension $n\cdot N$, the map is bijective. So the action on 1-jets over $y$ is simply transitive and has no fixed point. ∎

This is why the draft needs $\aleph$ (a7).

**Proposition F1 (the fiber; FC1–FC4).** Let $F=(\mathrm{Met}_{(p,q)}(\mathbb R^n),c_VG_\lambda)$ with $\lambda\ne1/n$.

1. $F$ is a pseudo-Riemannian symmetric space. The geodesic symmetry at $y_0$ is $S(z)=y_0z^{-1}y_0$ (FC3).
2. $F=\mathbb R\times SL(n,\mathbb R)/SO(p,q)$ isometrically, with $c_VG_\lambda=c_V\,n(1-n\lambda)\,ds^2+c_V\,\mathrm{tr}((u^{-1}du)^2)$, where $y=e^su$ and $|\det u|=1$ (FC4). The $\mathbb R$-factor changes sign at $\lambda=1/n$, which is Theorem A1.
3. The curvature is independent of $\lambda$: $c_VG_\lambda(R(h_c,h_d)h_b,h_a)=-\tfrac{c_V}4\mathrm{tr}([[H_c,H_d],H_b]H_a)$ with $H=y^{-1}h$ (FC2, all $N^4$ components for $n=2,3$).
4. $\mathrm{Ric}_F=-\tfrac n4\,G_{1/n}$ (the traceless projection), independent of $\lambda$ and $c_V$. Also $s_F=-n(n-1)(n+2)/(8c_V)$, which is $-9/c_V$ for $n=4$. FC1 checks this for every signature with $p\ge q$, $n\le4$, plus non-diagonal $y_0$. Since $G_\lambda$ is even in $y$, $(p,q)$ and $(q,p)$ are isometric.

For $(p,q)=(n,0)$ this is the noncompact-type space $SL(n)/SO(n)$, with sectional curvature $\le0$.

*Proof.* Item 2 follows from $y^{-1}dy=ds\,I+u^{-1}du$ with $\mathrm{tr}(u^{-1}du)=0$. Then 1, 3 and 4 follow from L-F4 on $\mathfrak{sl}(n)/\mathfrak{so}(p,q)$ for the flat-times-symmetric product. The constants are fixed by the exact computation. ∎

**Proposition F2 (the section; R1 in R2/R3).** For both readings, every $g$, every $\lambda$ and every $c$:

1. $s_*T_xX=H^\aleph_{g(x)}$, because the vertical part of $ds$ is $\nabla g=0$.
2. $s^*g_\aleph=c_H\,g$. So a1 holds with $c_H=1$.
3. $H\perp V$, hence $N_s=V$ isometrically with metric $c_VG_\lambda$. This is a10/a11 and GU-024.
4. In Z1, $g_\aleph$ depends on $g$ only through $\aleph$, and every $\aleph$-parallel metric is a section with horizontal image. (a8/a17 hold in Z1 and fail in Z2.)
5. In Z1, $\nabla^0:=\Theta_\aleph\nabla^{g_\aleph}\Theta_\aleph^{-1}$ is $g_C$-metric and torsion-free relative to the soldering form $\Theta_\aleph$, by construction. This closes the check requested in GU-G06. In Z2, $\Theta_\aleph$ is an isometry onto $(C,g_C)$ only along $s(X)$.

**Theorem F3 (Z2 along the section).**

- $\pi:(Y,g_\aleph)\to(X,c_Hg)$ is a pseudo-Riemannian submersion with totally geodesic fibers ($T\equiv0$).
- Its $A$-tensor $A_XY=-\tfrac12(R(X,Y)\cdot y)^{\rm vert}$ vanishes on $s(X)$.
- $\Phi(y)=g\,y^{-1}g$ is an involutive isometry with fixed set $s(X)$.

Consequently $s(X)$ is totally geodesic, every curvature component with an odd number of vertical slots vanishes there, and on $s(X)$:
$$\mathrm{Ric}^Y_{HH}=\mathrm{Ric}_g,\quad \mathrm{Ric}^Y_{HV}=0,\quad \mathrm{Ric}^Y_{VV}=-\tfrac n4G_{1/n},\quad s_Y=\frac{R(g)}{c_H}-\frac{n(n-1)(n+2)}{8c_V}.$$

*Proof.*

- $T\equiv0$: horizontal flows act on fibers by $y\mapsto P^*y$ with $P\in GL$, which are $G_\lambda$-isometries, so $T=0$ by Koszul.
- $A$: L-F5 gives the formula. $R(X,Y)\in\mathfrak o(g)$ annihilates $g$, so $A=0$ on $s(X)$.
- $\Phi$: parallel transport preserves $g$, so $\Phi$ commutes with it and preserves $H$ and $\pi^*g$. On each fiber $\Phi$ is the geodesic symmetry at $g_x$ (Prop. F1).
- The curvature values follow from L-F1 with $T=0$ and $A|_{s(X)}=0$ (Besse 9.36–9.37): $r(U,V)=\hat r(U,V)$, $r(X,Y)=\check r(X,Y)$ and $s=\check s+\hat s$. $\mathrm{Ric}_{HV}=0$ by $\Phi$-parity. ∎

Verification:

- R2: $n=2$, general metric, with the full $5$-dimensional Riemann tensor, Gauss equation, fiber block and mixed sum.
- R3: $n=3,4,5$; signatures $(3,0),(2,1),(1,3),(3,1),(4,0),(2,2),(5,0)$; non-normal coordinates; $\nabla\mathrm{Ric}\neq0$.
- R4: off the section, $s_Y=R/c_H+s_F-|A|^2$ with $|A|^2\neq0$, an exact check of L-F1 and of the $A$ formula.

**Theorem F4 (Z1 along the section, the draft-literal reading).** On $s(X)$, with $\mu:=\lambda/(1-n\lambda)$:

(a) The fibers are totally geodesic. The second fundamental form of $s(X)$ is
$$\mathrm{II}(X,Y)=-\tfrac{c_H}{2c_V}\big[X^\flat\odot Y^\flat+\mu\,g(X,Y)\,g\big],$$
with mean curvature $-g/(2c_V(1-n\lambda))$. So $s(X)$ is not totally geodesic, and $\nabla^0$ along $s(X)$ does not reduce to $\mathrm{Spin}(1,3)\times\mathrm{Spin}(6,4)$: its off-diagonal block is $\pm\mathrm{II}$, of size $c_H/c_V$.

(b) $s_Y=R(g)/c_H+\sigma_1/c_V$ with
$$\sigma_1=-\tfrac{n(n-1)(n+2)}8-\tfrac n4\Big[\tfrac{1+\lambda}{1-n\lambda}+\tfrac{n+1}2\Big];\qquad n=4:\ \sigma_1=-\tfrac{23}2-\tfrac{1+\lambda}{1-4\lambda}.$$

(c) The Ricci blocks are
$$\mathrm{Ric}^Y_{HH}=\mathrm{Ric}_g+\frac{c_H}{4c_V(n\lambda-1)}\,g,\qquad \mathrm{Ric}^Y_{HV}=0,\qquad \mathrm{Ric}^Y_{VV}=-\tfrac n4G_{1/n}-\tfrac14G_0=-\tfrac{n+1}4\,G_{1/(n+1)}.$$

*Proof.*

- (a) By Koszul, $2g_\aleph(\nabla_XY,U)=g_\aleph([X,Y],U)-U\,g_\aleph(X,Y)$ for basic $X,Y$. On $s(X)$, $[X^H,Y^H]^V=0$. Also $U(c_Hy(X,Y))=c_HU(X,Y)$, and $G_\lambda^{-1}=G_0^{-1}+\mu\,y\otimes y$ (Sherman–Morrison).
- (b) With a block-diagonal metric, $s=s_{HH}+2s_{HV}+s_{VV}$.
  - $s_{VV}=s_F$ (fibers totally geodesic).
  - $s_{HH}=R/c_H-|H|^2+|\mathrm{II}|^2$ (L-F2).
  - $s_{HV}=\mathrm{div}\,H_H-|\mathrm{II}|^2+|H_H|^2$ (L-F3, with $h_V=H_V=T_V=0$ and $T_H=A=0$ on $s(X)$).
  - The horizontal mean-curvature field is $H_H=-E/(2c_V(1-n\lambda))$ everywhere on $Y$, where $E=y\,\partial_y$ is the Euler field. Then $\mathrm{div}\,E=n/2$, because the volume density is homogeneous of degree $n/2-N$ in $y$.
  - Finally $|H|^2=n/(4c_V(1-n\lambda))$ and $|\mathrm{II}|^2=(N+n\mu)/(4c_V)$.
- (c), sketch (completed to a proof for every $n$ in `K-z1-ricci-general-n.md`, Theorem K, 2026-10-30):
  - Write $\mathrm{Z1}=\mathrm{Z2}+\delta$, $\delta=c_H(y-g)(\pi_*,\pi_*)$, which vanishes on $s(X)$.
  - On $s(X)$ the only non-zero parts of $\nabla\delta$ and $\nabla^2\delta$ are $(\nabla_U\delta)=c_HU$ and $\nabla^2_{U,V}\delta=\tfrac{c_H}2(Ug^{-1}V+Vg^{-1}U)$. These use the fiber Levi-Civita connection, which is $\lambda$-independent by Prop. F1(2).
  - So the corrections to Theorem F3 are algebraic in $g$. They can be computed in the flat model $c_Hy\,dx^2+c_VG_\lambda(dy,dy)$, where L-F1 applied to the submersion $(x,y)\mapsto y$ (flat fibers, $A=0$) gives $\mathrm{Ric}_{\mu\nu}=-\tfrac12\Delta h+\tfrac12(\partial h\,h^{-1}\partial h)-\tfrac14\partial\log|\det h|\cdot\partial h=\frac{c_H}{4c_V(n\lambda-1)}y$ and $\mathrm{Ric}_{AB}^{\rm mixed}=-\tfrac12\mathrm{tr}(h^{-1}\nabla^2_{AB}h)+\tfrac14\mathrm{tr}(h^{-1}\partial_Ah\,h^{-1}\partial_Bh)=-\tfrac14G_0$.
  - $\mathrm{Ric}_{HV}=0$ by parity under $\Phi$ ($\nabla\delta$ is odd). ∎

Status: (a) and (b) are proven. (c) is proven for $n=2$ (R2: symbolic, for every metric) and $n=1$ (E1: the hyperbolic plane, $K=-1/(4c_V(1-\lambda))$). For $n=3,4,5$ it is verified exactly on random metrics in every tested signature (R3). (Updated 2026-10-30: (c) is now proven for every $n$, every signature and every $\lambda\neq1/n$ in K (Theorem K), which completes the sketch above; K4 adds independent exact checks for $n=2,\dots,6$.)

**Corollary F5 (4D equations: $\lambda$- and route-dependence; R5).** Both (P) and (A) give the 4D vacuum Einstein equation with a cosmological term, $G(g)+\Lambda g=0$:

- (P): $s^*G^Y=G(g)+\Lambda_Pg$.
- (A): $S[g]=c_H^{n/2}\int(R/c_H+\sigma/c_V)\,dvol_g$, whose Euler–Lagrange equation is $G(g)+\Lambda_Ag=0$ with $\Lambda_A=-c_H\sigma/(2c_V)$ (L-F6).

| | $\Lambda_A$ | $\Lambda_P$ | $n=4$ |
|---|---|---|---|
| Z2 | $\frac{c_Hn(n-1)(n+2)}{16c_V}$ | $=\Lambda_A$ | $\frac{9c_H}{2c_V}$, for every $\lambda$ |
| Z1 | $\frac{c_Hn}{16c_V}\big[n^2+2n-1+\frac{2(1+\lambda)}{1-n\lambda}\big]$ | $\Lambda_A+\frac{c_H}{4c_V(n\lambda-1)}$ | $\Lambda_A=\frac{5c_H}{4c_V}\frac{5-18\lambda}{1-4\lambda}$, $\Lambda_P=\frac{3c_H}{2c_V}\frac{4-15\lambda}{1-4\lambda}$ |

Values for Z1 at $n=4$, in units of $c_H/c_V$:

| $\lambda$ | $\Lambda_A$ | $\Lambda_P$ |
|---|---|---|
| $0$ (Frobenius) | $25/4$ | $6$ |
| $\tfrac12$ (draft-literal) | $5$ | $21/4$ |
| $1$ (DeWitt) | $65/12$ | $11/2$ |

The zeros are at $\lambda=5/18$ for (A) and $\lambda=4/15$ for (P), both on the $(7,7)$ side but different from each other. Under the $(7,7)$ sign rule the physical constant is $\Lambda_{\rm phys}=-\Lambda$ in $(+,-,-,-)$ with $c_H/c_V>0$; the mirror convention gives the same.

- Z2: $\Lambda_{\rm phys}=-9|c_H/c_V|/2<0$, the anti-de Sitter sign, for every $\lambda$.
- Z1: de Sitter sign only for $\tfrac14<\lambda<4/15$ (P) or $<5/18$ (A).

In 2D, $G\equiv0$, so (P)/(A) have no solutions unless $\Lambda=0$.

*Remark (2026-10-30, round 11).* Reading (A) is what the ℷ-variation of (9.4) gives at first order in $T$ on the invariant sector, for the weight localised at $s(X)$ (`O-gimel-variation.md`, Cor. O11). For a ray weight μ, Z2 gives $\Lambda_A$ for every μ and Z1 gives $(m_2/m_1)\Lambda_A$. So where Z1 separates (A) from (P), the ℷ-variation selects (A).

**Theorem F6 (Einstein condition along the section; E2).** Let $n\ge2$. On $s(X)$, $\mathrm{Ric}(g_\aleph)=\kappa\,g_\aleph$ holds iff the reading is Z1, $\lambda=1/(n+1)$ and $\mathrm{Ric}(g)=0$. Then $\kappa=-(n+1)/(4c_V)$. Equivalently, (E) holds on $s(X)$ with $\Lambda_Y=-(d-2)(n+1)/(8c_V)$, which is $-15/(2c_V)$ for $n=4$.

At $\lambda=1/5$ ($n=4$) the fiber is $(7,3)$, so $g_\aleph$ has signature $(8,6)$, $(10,4)$ or a mirror (Theorem A4), never $(7,7)$. For $\lambda>\tfrac14$, and in Z2 for every $\lambda$, $g_\aleph$ is nowhere Einstein on $s(X)$. In particular $G(g_\aleph)\neq0$ on $s(X)$ for every $g$, $\lambda$ and $c$, because $\mathrm{Ric}^Y_{VV}\neq0$.

*Proof.* Use Theorems F3/F4, together with the fact that $G_\alpha\propto G_\lambda$ iff $\alpha=\lambda$ (compare the trace and traceless parts). ∎

E2 checks this on a Ricci-flat-at-the-point base ($\eta-\tfrac13Wxx$, with $W$ an exact random Weyl tensor). The condition holds for Z1 at $\lambda=1/5$. It fails at $\lambda=\tfrac12$ and at $\lambda=1$, fails for Z2, and fails for a non-Ricci-flat base. **Off the section** it fails for the same Weyl base (E2(v): 116 of 196 components non-zero). So for this base (E) holds only along $s(X)$, not on a neighbourhood in $Y$. For flat ℵ at $\lambda=1/(n+1)$ it holds on all of $Y$ (F7), and for $n=4$ over a definite base it holds on a neighbourhood of a point $p\in s(X)$ only if $\lambda=\tfrac15$ and ℵ is flat near $\pi(p)$ (`N-z1-neighbourhood.md`, N10). (Corrected 2026-10-30, round 10; the text had "So (E) holds only along $s(X)$" without the qualification; see `reports/CORRECTIONS.md`.)

**Proposition F7 (flat model = $SL(n+1)/SO$; E4).** For flat $\aleph$ and $\lambda=1/(n+1)$, Z1 with $c_H=2c_V$ is locally isometric to $SL(n+1,\mathbb R)/SO(p',q')$ with metric $c_V\,\mathrm{tr}((M^{-1}dM)^2)$. The map is $M=L^T\mathrm{diag}(u,1/\det u)L$, $L=\begin{psmallmatrix}1&x\\0&1\end{psmallmatrix}$, $y=\det(u)\,u$; other values of $c_H$ follow by rescaling $x$. For base $(1,3)$ this is $SL(5)/SO(3,2)$, of dimension 14 and signature $(8,6)$, and its Einstein constant $-(n+1)/(4c_V)$ (Prop. F1 at $n+1$) is $\kappa$. Theorem F6 says that a curved $\aleph$ preserves the Einstein condition along $s(X)$ exactly when $g$ is Ricci-flat.

## F.3 Conditional and conjectural (layer d)

**Proposition F8 (the draft's equation; conditional on R-S and Theorem C5).**

1. At the Levi-Civita point ($\varepsilon=1$, $\varpi=0$, so $T_\omega=0$), $\Upsilon=-{*}(G^Y_{ab}e^a\otimes\gamma^b)$. Its Hodge-dual pull-back is $s^*({*}\Upsilon)=-(G(g)+\Lambda_Pg)_{\mu\nu}dx^\mu\otimes\gamma^\nu$, with no $\gamma^A$ components because $\mathrm{Ric}_{HV}=0$. This is the precise content of (9.10)'s "$S_\omega\leftrightarrow R_{\mu\nu}-\frac s2g_{\mu\nu}$". The "$\Lambda g$" appears from the fiber geometry with $T=0$. However, $\Upsilon\ne0$ there for every $\gimel$ (Theorem F6).
2. For $\varpi=t\,\Phi^1$ (the soldering form, $T_\omega=t\Phi^1$), $d_{\nabla^0}\Phi^1=0$ (zero torsion). Then $F_{A_\omega}=F_{\nabla^0}+t^2\,\Phi^1\wedge\Phi^1$, the last term being the spin curvature of a constant-curvature tensor. By C5,
   $$\Upsilon=-{*}\big((G^Y+\Lambda_Y(t)g_\aleph)_{ab}e^a\otimes\gamma^b\big),\qquad \Lambda_Y(t)=\alpha t^2+\beta\kappa_1t,$$
   with $\alpha,\beta$ fixed by R-S normalisations that we have not computed.
3. So, along $s(X)$, $\Upsilon=0$ is (E). By Theorem F6 it is equivalent to 4D vacuum GR with $\Lambda_4=0$ exactly when Z1 holds, $\lambda=1/(n+1)$ and $\Lambda_Y(t)=-15/(2c_V)$ has a real root $t$. Otherwise it has no solution of this form.

Status: steps 1 and 3 follow from F3–F6 together with C5. Step 2 is a derivation sketch whose constants and sign conditions are open. The ε-dependence and $\varpi\notin\mathbb R\Phi^1$ are not explored. This bears on the Shiab normalisation studied separately (files `D-*`), which we did not touch.

**Conjecture F9.** For curved $\aleph$, no $\lambda$ makes $g_\aleph$ Einstein on an open set of $Y$ unless $\aleph$ is flat. Evidence: E2(v) and the non-vanishing of $A$ off $s(X)$. Not proven.

*Remark (2026-10-30, round 8).* In Z2 with $n=4$, $g_\aleph$ is Einstein on no open subset of $Y$, for every $\lambda\neq\tfrac14$ and every base metric, flat or curved (`L-nonradial-z2.md`, L2(a) and the remark after it). So F9 holds for Z2 at $n=4$; its open content is Z1, where the flat model F7 is Einstein.

*Remark (2026-10-30, round 10).* For Z1 at $n=4$: over a flat base, $\mathrm{Ric}=\beta_Sg_\aleph+b(g_\aleph|_H+c\,ds^2)$, so $g_\aleph$ is Einstein on an open set iff $\lambda=\tfrac15$ (`N-z1-neighbourhood.md`, N2). Over a definite base, Einstein on a neighbourhood of a point of $s(X)$ forces $\lambda=\tfrac15$ and $R_g=0$ there (N10: the order-$\varepsilon^2$ fibre pairing $\alpha|W|^2$, $\alpha=45/(8c_H^2c_V)$). Open: indefinite bases, and open sets away from $s(X)$.

## F.4 Grading (layers c–e)

| Item | Status |
|---|---|
| Lemma F0, Prop. F1, Prop. F2, Theorem F3 (Z2), Theorem F4(a)(b), Cor. F5 for Z2 and for $\Lambda_A$ (Z1), Prop. F7 | **Proven** (with L-F1…L-F6), each also checked exactly |
| Theorem F4(c) (Z1 Ricci blocks); hence $\Lambda_P$ (Z1) and the Z1 half of Theorem F6 | **Proven** for every $n$ (K, Theorem K and Cor. K4; 2026-10-30). Exact checks: R2, E1, R3, K4 ($n=2,\dots,6$) and, for F6, a Weyl base at $n=4$ (E2). The Z2 half of F6 (never Einstein) is proven. |
| Prop. F8 (link to (9.7)–(9.10)) | **Conditional** on R-S/C5; constants $\alpha,\beta$ open |
| Conjecture F9 | Conjectured; for Z2 at $n=4$ proven in a stronger form (no Einstein open set for any base; `L-nonradial-z2.md` L2); for Z1 at $n=4$ proven near $s(X)$ over definite bases (`N-z1-neighbourhood.md` N10); over flat bases its exception occurs exactly at $\lambda=\tfrac15$ (N2) |
| Which horizontal metric (Z1 vs Z2), $\lambda$, $c_V/c_H$, readout (P)/(A)/(E), domain of $\Upsilon=0$ | **Undefined in the draft**. The text favours Z1 (a3, a6, a8, a13, a17). |

Claims:

- **GU-006:** C. Well-posed under R-Z; a1 holds with $c_H=1$; resolves GU-G02.
- **GU-021:** C. The Zorro metric is Z1 $=\Theta_\aleph^*g_C$; resolves GU-G06.
- **GU-024:** V.
- **GU-054:** C. In the (P) and (A) senses, (9.10) holds with a fiber-induced $\Lambda$ that is $\lambda$-independent in Z2 and $\lambda$- and route-dependent in Z1. As a $Y$-equation it is exact 4D vacuum GR only for Z1 with $\lambda=\tfrac15$, which is incompatible with the draft's $(7,7)$. GU-G14 is closed under R-Z, and conditional on R-S for the link to $\Upsilon$.

**Testable (layer e).** No small parameter appears. In Z2 at $n=4$, $|\Lambda|=9/(2\ell_F^2)$ with $\ell_F^2=|c_V/c_H|$. Matching the observed $\Lambda=1.1056\times10^{-52}\,\mathrm m^{-2}$ needs $\ell_F\approx2.0\times10^{26}$ m (21 Gly; E3). A Planck-scale fiber overshoots by a factor $1.6\times10^{122}$, and the sign is anti-de Sitter. In Z1 only a tuning of $\lambda$ near $5/18$ (A) or $4/15$ (P) can make $\Lambda$ small.

**Coordinator reproduction (2026-10-30).** All three F scripts were re-run. Outputs are identical to the committed `.out` files apart from timing annotations. `F_zorro_reduction.py` needs about 27 min on 4 shared cores.
