# Lynchpin, "wave conjugations", "where n bubbles meet" (TH-020 – TH-025)

Code: `computation/howard/lynchpin_geometry.py` → `lynchpin_geometry.out`. Lean: `pentagon_faces_ge_12`, `no_six_pentagon_sphere` in `lean/HowardArith.lean`.

## 0. What the sources give (layer a)

| Datum | Source | Fidelity |
|---|---|---|
| "a Lynchpin structure may include six pentangular areas, and one or more of the pentangular areas may include a propulsion device"; structures combine into a "compound Lynchpin propulsion structure" | US 11,117,065 B2 abstract (TH-S6) | **Search-engine summary of Justia text**; full patent not obtained |
| Two Lynchpins joined "along three mated flanges"; "tetrahedrally supported Lynchpin structure 700 (FIG. 7)"; Fig. 3B embodiment "may include four pyramidal inner spaces" | same | same (paraphrase) |
| WO 2021/138203 claim: "pentangular structure vehicle" of two "six-sided 3-D structures" | TH-S6 | summary |
| "one on the top and then three in groups of three around that are spaced at 120 degrees apart" (Speaker C, describing the configuration) | JRE #2152 75:07 | transcript |
| "If I take where the twelve bubbles meet, that's where I made the linchpin from" | JRE #2152 36:58 | transcript |
| "This is what happens when four bubbles meet. This is the negative space where they can't touch each other. This is hydrogen … it has four contractive poles" | 31:28–31:54 | transcript |
| 8 bubbles: "eight contractive poles, but it only has six magnetic poles"; 6 bubbles: "30 poles"; 12 bubbles: "four unaccounted electrical poles"; 24 bubbles: "all of the electrical potential has been accounted for" | 32:35–34:29 | transcript |
| Design patent D1,002,750 "Circular internal lynchpin structure" | TH-S7 | title only (design patents claim drawings) |

**Well-posedness verdict on the sources alone:** **not well-posed.** No dimensions, angles or incidence data for the pentangular areas were obtained. We therefore test explicit reconstructions.

## 1. Reconstruction L0: a closed polyhedron with six pentagonal faces — impossible (c)

**Theorem G1.** Let a polyhedral 2-sphere (any cellular map on $S^2$) have all faces pentagons and all vertex degrees $\ge3$. Then $F\ge12$, with equality only for the dodecahedral map.
*Proof.* $2E=5F$, $2E\ge3V$ and $V-E+F=2$ give $F\ge12$. ∎ *Lean:* `pentagon_faces_ge_12`, `no_six_pentagon_sphere` (omega).
So "six pentangular areas" cannot be the full face set of a closed convex solid. The six areas must be walls, fins or partial faces, as in L-T.

## 2. Reconstruction L-T: six walls on the tetrahedral cone (b)

**Definition.** Let $V_1,\dots,V_4$ be the vertices of a regular tetrahedron centred at $O$, with $V=(\pm1,\pm1,\pm1)$ and an even number of minus signs. The *Lynchpin cone* is $\mathcal T=\bigcup_{i<j}W_{ij}$, where $W_{ij}=\{sV_i+tV_j:s,t\ge0\}$. This is the cone over the tetrahedron's 1-skeleton.
A *Lynchpin* is $\mathcal T\cap K$ for a $T_d$-invariant convex body $K\ni O$.

Why this reading:
- it has six planar walls ("six pentangular areas");
- its complement has four chambers ("four pyramidal inner spaces");
- it is "tetrahedrally supported";
- one ray is "on the top" and three surround it at 120°;
- it is exactly the junction "where four bubbles meet" in a foam (see below).

**Proposition T1** (exact unless noted):
1. The rays meet pairwise at $\theta_T=\arccos(-\tfrac13)=109.471220634°$. Seen from the top ray, the other three are $120°$ apart in azimuth.
2. Each wall plane $\operatorname{span}(V_i,V_j)$ is a mirror plane of $T_d$, and these six are all of them (each reflection permutes $\{V_k\}$).
3. Along every ray, three walls meet with pairwise dihedral angles $120°$.
4. $\mathbb R^3\setminus\mathcal T$ is the disjoint union of the four open trihedral cones over the faces.
   *Proof:* $\sum V_k=0$, so the $V_k$ positively span $\mathbb R^3$. The face-cones therefore cover $\mathbb R^3$, their interiors are disjoint, and their boundaries are the $W_{ij}$.
   *Check:* 200 000 random directions each lie in exactly one open face-cone.
5. (Imported: Plateau 1873; Taylor 1976.) Items 1–3 are Plateau's laws. $\mathcal T$ is Taylor's tetrahedral minimal cone, one of the only three tangent-cone types of soap-film-like $(M,0,\delta)$-minimal sets in $\mathbb R^3$. The other two are the plane and the $Y$-cone.

   So under L-T the Lynchpin is the four-bubble vertex of a foam. **This is a known object, not a new one.**

**Proposition T2 (pentagons).**
1. A regular pentagon's interior angle is $108°$. A wall's apex angle is $\theta_T$, a mismatch of $1.47122063449°$. *Six regular pentagons cannot meet at $O$ along the four rays.*
2. With apex angle $\theta_T$, a mirror-symmetric equilateral pentagon exists. Its angles are $109.4712°$, $106.8132°$ (×2) and $108.4512°$ (×2).
3. **Explicit well-posed instance:** take $K=\{x:V_k\cdot x\le3\ \forall k\}\cap[-2,2]^3$. Every wall becomes the planar pentagon
   $O,(1,1,1),(2,\tfrac12,\tfrac12),(2,-\tfrac12,-\tfrac12),(1,-1,-1)$,
   with edges $\sqrt3,\tfrac{\sqrt6}2,\sqrt2,\tfrac{\sqrt6}2,\sqrt3$ and angles $109.47°,90°,125.26°,125.26°,90°$ (sum $540°$).
   In general $K=\{V_k\cdot x\le h\}\cap[-g,g]^3$ yields pentagons for $h/3<g<h$, a one-parameter family.

**Proposition T3 ("circular" version).** Take $K=B(0,r)$. Each wall is a circular sector of angle $\theta_T$. The total wall area is $3\theta_T r^2=5.73189970875\,r^2$, and the density relative to a flat disc is $3\theta_T/\pi=1.82452034391$ (the $Y$-cone's is $3/2$).

**The idea inside the claim.** Two of Howard's statements have a precise content here:
- "the platonic solids were averages and approximations" (16:57);
- the dodecahedron becomes "naturally curved" when 20 cells are combined (36:50).

A cell with 12 pentagonal neighbours (the "twelve bubbles") would need flat faces with $108°$ corners and dihedral angle $\arccos(-1/\sqrt5)=116.565°$. Plateau requires $109.47°$ and $120°$, so flat-faced dodecahedra cannot satisfy Plateau's laws, and real foam faces curve. This frustration is classical (Kelvin 1887; Weaire & Phelan 1994). It is real, but it is not new.

## 3. Reconstruction L-B: "negative space where n bubbles meet" (b, c)

Place $n$ equal tangent unit balls at the vertices of a polyhedron with edge length 2. The "negative space" is the convex hull minus the balls.

- $n=4$ (tetrahedron):
  - volume of the hull: $2\sqrt2/3=0.942809$;
  - interstice: $\tfrac{2\sqrt2}{3}+\tfrac{4\pi}3-4\arccos\tfrac13=0.207762$;
  - covered fraction: $\sqrt2\,\Omega=0.779636$, where $\Omega=3\arccos\tfrac13-\pi=\arccos\tfrac{23}{27}$ (Rogers' simplex bound, 1958).
- $n=8$ (cube): the interstice is $8-4\pi/3=3.811210$, and the central gap has radius $\sqrt3-1=0.732051$.

**"Pole" counts.** The stated counts have no definition (TH-G10). They partly match (vertices, faces):
- tetrahedron $(4,4)$ ↔ "four contractive poles";
- cube $(8,6)$ ↔ "eight contractive … six magnetic".

No candidate cluster reproduces the counts for $n=6,12,24$, including octahedron, icosahedron, cuboctahedron, truncated solids and snub cube (table in `.out`). For example, the octahedron gives $(6,8)$ against the stated "30 poles". **Verdict:** not reproducible from any vertex-transitive cluster we tried. A further reconstruction would need an independently specified counting rule before these counts could be used as a test (P-TH-4).

## 4. Status of "wave conjugations", "supersymmetric systems", "predicts every star"

These terms are not defined in any source we obtained, and no construction rule is given. They are not formalizable (TH-G08, TH-G09). The prediction "when we put this in the solar apex … it predicts every star on the galactic plane" (149:53) would be testable if the construction and a star catalogue were specified. Neither is (TH-G11).

## References (imported)

- Plateau, J. (1873). *Statique expérimentale et théorique des liquides soumis aux seules forces moléculaires.* Gauthier-Villars.
- Taylor, J. E. (1976). The structure of singularities in soap-bubble-like and soap-film-like minimal surfaces. *Ann. of Math.* 103, 489–539.
- Rogers, C. A. (1958). The packing of equal spheres. *Proc. London Math. Soc.* (3) 8, 609–620.
- Thomson, W. (Lord Kelvin) (1887). On the division of space with minimum partitional area. *Phil. Mag.* 24, 503–514.
- Weaire, D. & Phelan, R. (1994). A counter-example to Kelvin's conjecture on minimal surfaces. *Phil. Mag. Lett.* 69, 107–110.
