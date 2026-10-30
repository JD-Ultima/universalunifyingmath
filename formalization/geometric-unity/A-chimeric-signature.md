# A. Observerse fiber metrics, the chimeric signature, and the spin structure on $C$

Claims: GU-010…015, GU-022, GU-023 (`claims/geometric-unity.md`). Code: `computation/geometric-unity/fiber_signature.py` → `fiber_signature.out`. Metric record: `metrics/geometric-unity/A-chimeric-signature.yaml`. Imported lemmas L-A1…L-A4 are listed in `sources/geometric-unity.md`.

## A.0 Source assertions (layer a)

- [W21 §3.1, (3.8)–(3.9)] The vertical tangent space of $Y=\mathrm{Met}(X^n)$ at $g$ carries the Frobenius form $\langle A,B\rangle_g=\mathrm{tr}_g(A\cdot B)$. Its traceless part has signature $\bigl(in-i^2,\tfrac{n^2+2i^2-2in+n-2}{2}\bigr)$, which is $(3,6)$ for $(i,n)=(1,4)$. The trace part may be "assigned" either $(0,1)$ or $(1,0)$. With the second choice the fiber signature is $(4,6)$. Fn. 2 calls this "one of the few non-forced choices".
- [W21 §3.5, p. 22] Combined with the base metric of signature $(1,3)$, the possibilities are $(4,10)$, $(8,6)$, $(5,9)$ or $(7,7)$; $(9,5)$ is mentioned as an alternative. "We will assume that the metric on $Y$ is split with signature $(7,7)$."
- [W21 (3.10)] $C=V\oplus H^*$, where $H^*=\pi^*T^*X$ carries the metric induced by $g$ at the point $g\in Y$.
- [W21 (3.32)–(3.33)] $P_H=P_{\widetilde{\mathrm{Fr}}(C^{7,7})}\times_{\rho_D}H$ over $Y$, whose fiber over $X$ is written $\widetilde{GL}(4,\mathbb R)/\mathrm{Spin}(1,3)$.

## A.1 Reconstruction (layer b)

**R-A1.** $X$ is a smooth $n$-manifold. $Y=\mathrm{Met}_{(p,q)}(X)\subset \mathrm{Sym}^2T^*X$ is the open subbundle of nondegenerate forms of signature $(p,q)$. Its vertical bundle is $V=\pi^*\mathrm{Sym}^2T^*X$, because the fiber is an open subset of a vector space, and $H^*=\pi^*T^*X$.

**R-A2. Natural fiber metrics.** A *natural* fiber metric is a field of symmetric bilinear forms on $V$ that is invariant under the fiberwise action of $GL(n,\mathbb R)$, $g\mapsto A^{T}gA$, $h\mapsto A^{T}hA$. Here it is the DeWitt family
$$G_\lambda(h,k)\big|_g=\mathrm{tr}(g^{-1}hg^{-1}k)-\lambda\,\mathrm{tr}(g^{-1}h)\,\mathrm{tr}(g^{-1}k),\qquad \lambda\in\mathbb R,$$
up to an overall factor $c_V\neq0$. Lemma A2 shows nothing else is natural. $\lambda=0$ is the source's Frobenius form. Flipping the sign of the trace part, the source's free choice, corresponds to $\lambda>1/n$ (Theorem A1).

**R-A3. Chimeric metric.** $g_C=c_V\,G_\lambda\oplus c_H\,g^{-1}$ on $C=V\oplus H^*$, with $c_V,c_H\in\mathbb R^\times$. Lemma A3 shows that naturality forces $V\perp H^*$.

**Signature convention.** $(a,b)$ means $a$ positive and $b$ negative directions. "Base $(1,3)$" means $g$ has one positive eigenvalue. The source lists the Frobenius signature in the order (negative, positive) in (3.9) and §3.5, but in (4.6) it writes $\mathrm{Spin}(6,3)\times\mathrm{Spin}(0,1)$, which is (positive, negative). Check A6 confirms that its formula matches ours with the order reversed.

## A.2 Results (layer c)

**Theorem A1 (fiber signature).** Let $g$ have signature $(p,q)$, $n=p+q\ge1$, and $P=\tfrac{p(p+1)}2+\tfrac{q(q+1)}2$. Then
$$\mathrm{sig}\,G_\lambda=\begin{cases}(P,\;pq) & \lambda<1/n,\\ (P-1,\;pq+1)& \lambda>1/n,\end{cases}$$
and for $\lambda=1/n$ the form is degenerate with kernel $\mathbb R g$. In particular, for $n=4$: base $(1,3)$ or $(3,1)$ gives $(7,3)$ for $\lambda<\tfrac14$ and $(6,4)$ for $\lambda>\tfrac14$; base $(4,0)$ gives $(10,0)$ and $(9,1)$.

*Proof.* $G_\lambda$ is $GL(n)$-invariant (verified symbolically, A3) and even in $g$. By L-A1 we may take $g=\eta=\mathrm{diag}(\epsilon_1,\dots,\epsilon_n)$. In the basis $E_{ii}$, $E_{ij}+E_{ji}$ ($i<j$), the off-diagonal vectors are mutually orthogonal and orthogonal to the $E_{ii}$, with $G_\lambda(E_{ij}+E_{ji},E_{ij}+E_{ji})=2\epsilon_i\epsilon_j$. That gives $\binom p2+\binom q2$ positive and $pq$ negative directions. On $\mathrm{span}\{E_{ii}\}$ the Gram matrix is $I-\lambda\epsilon\epsilon^T$, with eigenvalue $1$ on $\epsilon^\perp$ (multiplicity $n-1$) and $1-\lambda n$ on $\epsilon=\eta$. Adding up gives the claim. ∎ Exact checks: the symbolic spectra are $\{1-4\lambda,1^{(3)},\pm2\}$ (A1). The inertia was computed by two independent exact methods for 7 values of $\lambda$, including $\lambda=\tfrac14$, on diagonal and random non-diagonal $g$ (A2). The closed form was checked for every $(p,q)$ with $2\le n\le 6$ (A5).

**Lemma A2 (naturality ⇒ DeWitt).** For $n\ge3$, every $GL(n)$-invariant symmetric bilinear form on $V$ equals $a\,\mathrm{tr}(g^{-1}hg^{-1}k)+b\,\mathrm{tr}(g^{-1}h)\mathrm{tr}(g^{-1}k)$.
*Proof.* $GL(n)$ acts transitively on the fiber, so an invariant form is determined by an $O(p,q)$-invariant form on $T_g=\mathrm{Sym}^2=\mathbb R g\oplus\mathrm{Sym}^2_0$. $\mathrm{Sym}^2_0$ is absolutely irreducible and nontrivial (L-A2), so $\mathrm{Hom}_{O(p,q)}(\mathbb Rg,\mathrm{Sym}^2_0)=0$ and its invariant forms are multiples of the trace form. ∎

**Lemma A3 (orthogonality forced).** Every $GL(n)$-invariant bilinear pairing $V\times H^*\to\mathbb R$ vanishes, because $-I\in O(p,q)$ acts by $+1$ on $\mathrm{Sym}^2$ and by $-1$ on $T^*_x$. ∎

**Theorem A4 (when the chimeric signature is $(7,7)$; $n=4$).** Take base signature $(1,3)$, $(3,1)$ or $(4,0)$, and $g_C=c_VG_\lambda\oplus c_Hg^{-1}$ with $\lambda\ne\tfrac14$. Then $\mathrm{sig}(g_C)=(7,7)$ iff all three of the following hold:
1. $g$ is Lorentzian;
2. $\lambda>\tfrac14$;
3. $c_Hg^{-1}$ has exactly one eigenvalue with the sign of $c_V$. Equivalently, $\mathrm{sign}(c_Vc_H)=+1$ for base $(1,3)$ and $-1$ for base $(3,1)$.

For $\lambda<\tfrac14$ the possible signatures are $(8,6)$, $(10,4)$, $(6,8)$, $(4,10)$. For $\lambda>\tfrac14$ they are $(7,7)$, $(9,5)$, $(5,9)$. For a Riemannian base they are $(14,0)$, $(10,4)$, $(13,1)$, $(9,5)$ and their mirrors, never $(7,7)$.
*Proof.* Add the signatures from Theorem A1 to those of $\pm g^{-1}$. All 24 cases are listed in A4, with exact $14\times14$ cross-checks. ∎

So the source's list in §3.5, $\{(4,10),(8,6),(5,9),(7,7)\}$ plus $(9,5)$, is correct once its trace "assignment" is realised as $\lambda>\tfrac14$. Its "(4,6)" is our $(6,4)$ written in the reversed order. The $(7,7)$ choice is therefore not forced by the geometry: it needs the DeWitt parameter on the supermetric-type side $\lambda>1/n$ and a relative sign between the vertical and horizontal blocks. The standard DeWitt value $\lambda=1$ qualifies.

**Corollary A5 (classification of neutral chimeric signatures, all $n$).** With $g_C$ as in R-A3, $\lambda\neq1/n$ and $n\ge2$, the signature is neutral $(N,N)$ exactly for base $(p,q)\in\{(1,3),(3,1),(2,2)\}$, where $N=7$, and $(p,q)\in\{(2,3),(3,2)\}$, where $N=10$. All cases have $\lambda>1/n$.
*Proof.* Write $d=p-q$ and use $P-pq=\tfrac12(d^2+p+q)$. For $\lambda>1/n$ with horizontal block $(p,q)$, neutrality is $d^2+3p-q-4=0$, i.e. $p=\tfrac12(4-d-d^2)$. Requiring $p,q\ge0$ leaves $d\in\{-2,-1,0,1\}$, giving $(1,3),(2,3),(2,2),(1,0)$. For $\lambda<1/n$ it is $d^2+3p-q=0$, giving only $n\le1$. The swapped horizontal block is the mirror image. A5 confirms this by exhaustive search up to $n=40$. ∎

**Proposition A6 (topology of the fiber; the space $Y$).** The fiber $GL(4,\mathbb R)/O(1,3)$ of $\mathrm{Met}_{(1,3)}(X)$ is homotopy equivalent to $\mathbb{RP}^3$. The fiber written in (3.33) is $\widetilde{GL}^+(4,\mathbb R)/\mathrm{Spin}_0(1,3)=GL^+(4)/SO_0(1,3)$, which is homotopy equivalent to $S^3$ and is the connected double cover: the space of time-oriented Lorentzian forms.
*Proof.* Fix a Euclidean $\delta$. The self-adjoint operator $\delta^{-1}g$ has exactly one eigenvalue of the minority sign (L-A1). The map $g\mapsto$ (its eigenline) $\in\mathbb{RP}^3$ is continuous, and its fibers $\cong\mathbb R_{>0}\times\mathrm{PosDef}(3)$ are convex. Using the eigen*ray* instead gives the $S^3$ version. The preimage of $SO_0(1,3)$ in $\widetilde{GL}^+(4)$ is connected because $\pi_1(SO(3))\to\pi_1(GL^+(4))$ is an isomorphism. ∎
The source therefore uses two different spaces under the name $Y$: $\mathrm{Met}(X)$ in Def. 3.1 and its double cover in (3.33). See gap GU-G03.

**Proposition A7 (spin structure on $C$).** For $n=4$ and either choice of $Y$: $w_1(C)=0$ and $w_2(C)=\pi^*w_2(TX)$, and $\pi^*$ is injective on $H^2(\,\cdot\,;\mathbb Z_2)$. Hence $C\to Y$ admits a spin structure, and so the bundle $P_{\widetilde{\mathrm{Fr}}(C)}$ of (3.32) exists, **iff $X$ is spin**. Spin structures on $C$ are a torsor over $H^1(Y;\mathbb Z_2)$. This is $H^1(X;\mathbb Z_2)\oplus\mathbb Z_2$ for $Y=\mathrm{Met}(X)$ and $H^1(X;\mathbb Z_2)$ for the time-oriented $Y$. In the first case even a $Y$ over an $X$ with a unique spin structure has two spin structures on $C$.
*Proof.* $V=\pi^*\mathrm{Sym}^2T^*X$ and $H^*=\pi^*T^*X$. By the splitting principle (L-A4) and the computation in A7: $w(\mathrm{Sym}^2)$ has $w_1=w_1(X)$ and $w_2=w_1(X)^2$, and $w(C)$ has $w_1=0$ and $w_2=w_2(X)$. For injectivity, Prop. A6 gives a fiber homotopy equivalence $Y\simeq P(TX)$ or $S(TX)$. Leray–Hirsch for $P(TX)$ and the Gysin sequence for $S(TX)$ (Milnor–Stasheff §12, Hatcher 4.D) make $\pi^*$ injective in degrees $\le2$. Both also give the stated $H^1$. ∎

## A.3 Status (layers d, e)

- GU-011, GU-013: proven as stated, modulo the order convention and with the trace choice realised by $\lambda>1/n$.
- GU-012: the free choice is exactly the choice of a side of $\lambda=1/n$. It is natural, since both sides are $GL$-invariant.
- GU-014: the canonical part is $C$ with $g_C$, given the data $(\lambda,c_V,c_H)$. The isomorphism $C\cong TY$ needs a horizontal distribution, so it is not canonical.
- GU-015: the fiber written in (3.33) is the double cover of the fiber of $\mathrm{Met}(X)$ (Prop. A6).
- GU-022, GU-023: $P_H$ exists iff $X$ is spin. The topological data include $(\lambda,c_V,c_H)$ and the choice of $Y$ (Prop. A7).
- Conjectural or physical content: none. The selection of $(7,7)$ over $(9,5)$ is stated as a choice ("more balanced"). The fact that a neutral chimeric signature essentially singles out $n=4$ (with $n=5$ the only other case, and there $\dim V=15$ is odd) is recorded as pattern GU-P01.
- Failure conditions: Theorem A1 would fail if some inertia computation disagreed with the closed form; none of the 7×3×3 + 2·25 exact checks do. Theorem A4 depends on R-A2. A non-natural fiber metric, for example one that depends on a background structure, is outside its scope.

## Independent cross-check (coordinator)

The fiber and chimeric signature tables were recomputed independently, before this note was read, with a separate script. The script builds the numerical Gram matrix of $\mathrm{tr}(g^{-1}hg^{-1}k)-\lambda\,\mathrm{tr}(g^{-1}h)\mathrm{tr}(g^{-1}k)$ on a basis of $\mathrm{Sym}^2\mathbb R^4$ for $g$ of signature $(1,3)$, $(3,1)$, $(4,0)$ and $\lambda\in\{0,0.2,0.25,0.3,0.5,1\}$. Every fiber signature agrees: degenerate at $\lambda=\tfrac14$. Every block-sign combination agrees as well: $(7,7)$ only for a Lorentzian base with $\lambda>\tfrac14$ and the stated relative sign, and $(13,1)$ for a Riemannian base with $\lambda>\tfrac14$. The script is `computation/geometric-unity/xcheck_dewitt.py`.
