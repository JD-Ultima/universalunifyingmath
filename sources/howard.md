# Sources — Terrence Howard, "Terryology" / "Grand Unified Theory"

Accessed 2026-10-30 from a sandbox with restricted egress. Most news, patent and video hosts were refused by the proxy or did not resolve: rollingstone.com, wikipedia.org, patents.google.com, USPTO ppubs, justia.com, grantland.com, archive.org, youtube.com. Only WebSearch result summaries, GitHub and raw.githubusercontent.com were usable.

The source register distinguishes interviews, transcript material, self-published writing, social-media statements, and patent records. Acquisition status is listed per item; several primary texts were unavailable. Publication type alone does not decide the mathematical claims. Howard mentions "three papers that's been peer reviewed" (JRE #2152 136:09), but these were not located (TH-G12).

Timestamps "mm:ss" follow the transcript file, so 124:37 means 2:04:37. Downloaded files are kept outside the repo, under the session scratchpad (`howard-dl/`), and are not committed.

| ID | Item | Date | Access | sha256 / pin | Use |
|---|---|---|---|---|---|
| TH-S1 | *Rolling Stone*, "Terrence Howard's Dangerous Mind" (profile; URL `rollingstone.com/tv-movies/tv-movie-features/terrence-howards-dangerous-mind-37057/`) | Sept 2015 | **Original not obtained** (proxy 403). Howard's quotes were obtained verbatim as reprinted by CBS News, "Terrence Howard challenges the times tables" (2015-09-15). That article is archived sentence-split in the WMT16 corpus file `newsdev2016.en`, lines 205–230. | `newsdev2016.en` sha256 `c59454c1…86f3e0` (rsennrich/wmt16-scripts @ `9695851c`) | H1, SQ, EFF; "Terryology"; building-block patent pointer; "Twenty years from now…" |
| TH-S2 | *The Joe Rogan Experience* #2152, "Terrence Howard" (YouTube `g197xdRZsW0`) | 2024-05-18 | **Third-party machine transcript with speaker labels** from GitHub `zreese/jre-excerpts`, `content/jre-2152--terrence-howard--2024-05-18.md`, 2141 lines covering 00:02–188:44. Not checked against the audio, so speaker attribution and digits may be wrong. | sha256 `3cb1de0a…a2c7c6` (repo @ `e4338df1`) | Main verbatim source for TH-001…TH-025 |
| TH-S3 | *JRE* #2171, "Eric Weinstein & Terrence Howard" | 2024-07-01 | **Not obtained.** Known only from search summaries. Reported: Howard called 1×1=2 "just a metaphor for challenging the status quo". | — | Status caveat on TH-001 |
| TH-S4 | Oxford Union address | Dec 2017 (reported) | **No transcript or recording obtained.** Howard's own account is in TH-S2 124:01–125:45: he showed the √2 loop on calculators. | — | TH-005 |
| TH-S5 | Book *One Times One Equals Two*, "initial public draft", at `tcotlc.com` (site title "Terrence Howard : One Times One Equals Two"). Secondary reports describe "a free 162-page online book". Scribd preview `OTOET-PREVIEW-044-NOV7-2019`. | ≤2019 | **Not obtained** (host unreachable). Page refs from TH-S2: p. 64 (Platonic solids), p. 134 (flower "juxtapose"). | — | **The key missing primary text** (TH-G01) |
| TH-S6 | US 11,117,065 B2 "Systems and methods for lynchpin structure applications", inventor T. Dashon Howard. Filed 2020-01-03, granted 2021-09-14. Pub. US 2021/0205727 A1; WO 2021/138203. | 2021 | **Full text not obtained.** Abstract and figure captions are known only through WebSearch summaries of Justia/sumobrain pages (paraphrase, not verbatim). | — | TH-020 (six pentangular areas; four pyramidal inner spaces; tetrahedral support) |
| TH-S7 | US D1,002,750 "Circular internal lynchpin structure" (design) | 2023-10-24 | Title only (via search) | — | TH-020 (L-T.7) |
| TH-S8 | US 2015/0079872 A1; US 9,259,660 B2; US 9,427,676 B2, "Systems and methods for enhanced building block applications" (T. Dashon Howard). Claim 1 of '660 is reported to recite a building block with a piezoelectric element and an LED. A related "All-Shape" family (US 9,192,875; 9,168,465) is reported to have tetrahedron-edge flanges with magnets. | 2015–2016 | Summaries only | — | TH-021 (blocks may be "platonic solids") |
| TH-S9 | R. J. Lipton & K. W. Regan, "Does 1 Times 1 Equal 2?", *Gödel's Lost Letter* | 2018-05-22 | Search summary only | — | **Commentary (lead, not evidence).** Prior art for the $2AB$ renormalisation and the $2\times2$ tension. A commenter there states that Howard keeps $2\times2=4$ (unverified). |
| TH-S10 | Tweet `x.com/terrencehoward/status/925754491881877507`: "This is the proof to the World of Science and Mathematics that 1x1=2. #mathematics #worldofscience @rollingstone #proof" | 2017-11-01 | Text from the search-result title. Attachment not obtained. | — | Dating of the "proof" |
| TH-S11 | KTLA red-carpet interview, Emmys 2019 ("open up the flower of life and find the real wave conjugations"), quoted by *The Root* | 2019 | Search summary | — | TH-023 context |
| TH-S12 | Academia.edu, "Tetrahedron with Six Pentagons Geometry and the Lynchpin Concept of Terrence Howard" (third party) | — | Title and summary only | — | Lead for reconstruction L-T (pentagon vertex at the tetrahedron centroid; 109.47° vs 108°) |
| TH-S13 | Third-party notes on TH-S2, GitHub `evilusean/Pysean/.../JRE-TerrenceHoward-Notes.md` | — | Obtained | sha256 `e1814efa…de474` | Cross-check only (frequency list) |

Critical responses mentioned in search results (N. deGrasse Tyson 2024 video; S. Hossenfelder 2024; blogs) were **not used as evidence** and contributed no premises.

## Imported theorems and data (cited in `formalization/howard/*.md`)

| Result | Reference | Used in |
|---|---|---|
| Cauchy equation: regular additive maps $\mathbb R\to\mathbb R$ are linear | Aczél 1966, Ch. 2 | A3 on $\mathbb R$ |
| Non-measurable additive maps exist | Hamel 1905 | A3 remark |
| Continuous, strictly increasing, associative operations on an interval are conjugate to + | Aczél 1949/1966 §6.2; Craigen & Páles 1989 | A7 |
| Near-ring definition | Pilz 1983 | A2 corollary |
| Tropical and log semirings | Maclagan & Sturmfels 2015; Litvinov 2007 | §5 |
| Plateau's laws; Taylor's classification of minimal cones | Plateau 1873; Taylor 1976 | T1.5 |
| Rogers' simplex bound $\sqrt2\,\Omega$ | Rogers 1958 | L-B |
| Kelvin problem; Weaire–Phelan | Kelvin 1887; Weaire & Phelan 1994 | T2 remark |
| $c=299\,792\,458$ m/s exact; 1 mi = 1609.344 m | BIPM SI Brochure 2019; 1959 agreement | N3 |
| Pythagorean tuning, cents | Benson 2006 | N4 |
| Element data | mendeleev 1.3.0 (PyPI), bundled data | N4 |
