# Geometric Unity (E. Weinstein): sources

Access date for everything below: 2026-10-30. Network policy in this environment: `geometricunity.org`, `geometricunity.nyc3.digitaloceanspaces.com`, `web.archive.org`, `scribd.com`, `fliphtml5.com`, `semanticscholar.org`, `theportal.wiki`, `timothynguyen.org` and `arxiv.org` could not be reached (DNS failure or proxy 403). `github.com` was reachable.

## 1. Primary sources

### [W21] *Geometric Unity: Author's Working Draft, v 1.0*, E. Weinstein, dated Thursday April 1st, 2021 (69 pp.). **Obtained via a third-party mirror. Not verified against the original host.**

| Field | Value |
|---|---|
| Original URL (unreachable) | `https://geometricunity.nyc3.digitaloceanspaces.com/Geometric_Unity-Draft-April-1st-2021.pdf` |
| Copy used | `github.com/CarbonEdge/moonstar-physics`, commit `9f39e97b0c6d1a9c61a325f4b63da0d6c6495d89` (2026-10-11), path `papers/pdfs/geometric-unity.pdf` (identical copy at `docs/papers/geometric-unity.pdf`) |
| sha256 | `3f28d742234a9841fc8e51ff172053200aa3eddf3ece38154a3328b9ebd186d4` |
| Size, pages | 2,087,649 bytes; 69 pages; US Letter |
| PDF metadata | Creator TeX; Producer pdfTeX-1.40.21; CreationDate = ModDate = Thu Apr 1 16:26:28 2021 UTC |
| Consistency checks | Title page, date and table of contents match the bibliographic data quoted by search engines and Wikipedia ("Author's Working Draft, v 1.0", April 1, 2021). The internal timestamps match the release date. |
| Status | Treated as the primary text. All Claim-Register locations (section, equation, page) refer to this file. If a hash of the original becomes available and differs, every GU claim location must be re-checked; log that in `reports/CORRECTIONS.md`. |
| Rights | The front-matter footnote says the work "may not be built upon, renamed, or profited from without express permission of the author". This repository does not redistribute the PDF. We quote short passages for analysis and criticism and transcribe tables only to compare them with our computations. |

How it was read: `pdftotext -layout` (poppler) for text. Pages 19, 22, 41, 43, 44 and 51–53 were also rendered with `pdftoppm` at 80–220 dpi to check equations and tables that the text layer garbles. The Shiab equation (9.3) on p. 43 and the bars on SU(3) labels in the p. 53 table were read from the renders.

Self-declared status (quoted for context, not evidence): the author describes the text as "a draft of work in progress", "stitched together from entirely heterogeneous sources", with "multiple sign conventions" in Sec. 6 (fn. 7), a Shiab operator whose original choice "cannot yet now locate" (fn. 10), and a deformation-complex diagram that "may contain some inconsistancies" (p. 49).

### [W13] Oxford lectures, May 2013. **Not obtained.**
Search-result snippets list a Martin Wood Lecture on 23 May 2013 and a seminar "Geometric Unity" on 31 May 2013 (maths.ox.ac.uk event pages). Those pages were not fetched, so the dates are unverified. No video or transcript was found. The theportal.wiki entry is reported, again only in a snippet, to have no edited transcript, and the host is blocked here anyway. Nothing in this repository relies on the content of these lectures.

### [W20] "A Portal Special Presentation – Geometric Unity: A First Look", YouTube, 1 April 2020 (the draft's ref. [10]). **Not obtained.** We had neither the video nor a transcript (theportal.wiki is blocked).

### [W92] E. Weinstein, *Extension of self-dual Yang–Mills equations across the 8th dimension*, PhD thesis, Harvard, 1992 (the draft's ref. [8]). **Not obtained**, and not used.

## 2. Secondary sources: leads only, not evidence

| Ref | Item | Access | Use |
|---|---|---|---|
| [NP21] | T. Nguyen, T. Polya, "A Response to Geometric Unity" (2021) | **Not obtained** (timothynguyen.org blocked; no mirror found) | Search snippets describe it only as a critique in which the undefined Shiab operator is central. We do not attribute any specific argument to it. The objections we test are the obvious ones raised by the draft itself: well-definedness and gauge covariance of the Shiab operator, the complexification Spin(7,7)→U(64,64), and consistency of the field equations. Phase-2 results on these are in `formalization/geometric-unity/C-gauge-shiab.md` and `B-spinor-branching.md`. |
| — | `CarbonEdge/moonstar-physics` `reviews/geometric-unity/*` | Present in the mirror repository | An AI-generated review. Not read for content and not used. |

## 3. Imported lemmas (standard results used in `formalization/geometric-unity/`)

| ID | Statement | Reference |
|---|---|---|
| L-A1 | Sylvester's law of inertia. | any linear algebra text, e.g. Lang, *Algebra*, XV §4 |
| L-A2 | For $n\ge 3$ the traceless symmetric square $\mathrm{Sym}^2_0(\mathbb C^n)$ is an irreducible $SO(n,\mathbb C)$-module (highest weight $2\omega_1$). An identity-component real form $SO_0(p,q)$ is Zariski-dense in $SO(n,\mathbb C)$, so a real $SO_0(p,q)$-invariant bilinear form on $\mathrm{Sym}^2_0$ is a multiple of the trace form (Schur). | Fulton–Harris, *Representation Theory*, §19.5; Borel, *Linear Algebraic Groups*, 18.3 |
| L-A3 | A real symmetric matrix has only real eigenvalues. For a real-rooted polynomial, Descartes' rule of signs counts the positive roots exactly. | standard |
| L-A4 | Splitting principle; $w$ of tensor products of line bundles; Whitney sum formula. | Milnor–Stasheff, *Characteristic Classes*, §4, §7 |
| L-B1 | $\mathrm{Cl}(p,q)\cong\mathbb R(2^{n/2})$ when $p-q\equiv 0,2 \pmod 8$, with $n=p+q$ even. Real Majorana–Weyl spinors exist iff $p-q\equiv 0\pmod 8$. Weyl spinors are complex conjugate to each other iff $p-q\equiv \pm 2 \pmod 8$, and quaternionic iff $p-q\equiv 4\pmod 8$. | Lawson–Michelsohn, *Spin Geometry*, I.4; Atiyah–Bott–Shapiro, *Topology* 3 (1964) |
| L-B2 | Weights of the vector and half-spin modules of $D_m$; characters determine finite-dimensional representations of a connected reductive group; branching to a regular subalgebra is computed by restricting weights. | Fulton–Harris §20, §25 |
| L-B3 | A real irreducible representation whose complexification is $W\oplus\overline W$ with $W\not\cong\overline W$ irreducible is of complex type. Its commutant is $\mathbb C$, and it is the realification of $W$. | Bröcker–tom Dieck, *Representations of Compact Lie Groups*, II.6 (the algebraic statement does not need compactness) |
| L-B4 | The standard chain $\mathrm{Spin}(10)\supset\mathrm{Spin}(6)\times\mathrm{Spin}(4)=SU(4)\times SU(2)_L\times SU(2)_R\supset SU(3)\times SU(2)_L\times U(1)_Y$ with $Y=\tfrac12(B-L)+T_{3R}$. The $16$ contains one SM generation plus $\nu^c$. | Pati–Salam, PRD 10 (1974) 275; Slansky, Phys. Rep. 79 (1981) 1, Tables |
| L-C1 | The gauge transformation of a connection is $A\cdot h=h^{-1}Ah+h^{-1}dh$. Covariant derivative of a section of $\mathrm{Ad}\,P$: $d_{A}h=dh+[A,h]$. Curvature covariance $F_{A\cdot h}=h^{-1}F_Ah$. | Kobayashi–Nomizu I, II.5 |
| L-C2 | Integration by parts for a first-order operator $L=K\circ\nabla$ whose coefficient tensor $K$ is $\nabla$-parallel. Formal self-adjointness is equivalent to a pointwise symbol identity. | Lawson–Michelsohn II.5 (formal adjoints) |
| L-C3 | Algebraic curvature tensors satisfy $R=W+P\owedge g$ (Kulkarni–Nomizu product with the Schouten tensor), where $W$ is totally traceless. | Besse, *Einstein Manifolds*, 1.G (1.110–1.116) |
