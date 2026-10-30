"""Curated data for the round-3 mass-cell audit (masscells.py).

Every entry carries a verbatim anchor: (statement number, string). masscells.py checks each string
against the statement unit (box + following prose, bridge and derivation, segmented as in
scorecard_extract.py, whitespace-normalised). FILE_ANCHORS are checked against a tex file instead
(narrative text that belongs to no statement unit). Nothing here is computed; it is transcription.
"""

# ---------------------------------------------------------------- trace anchors (task 1)
TRACE = [
    # the electron cell 67
    ("152.1", "27 + 40 + 12 = 79, \\qquad \\operatorname{rank}(e) = 79 - 12 = 67"),
    ("152.1", "so the unit of the mass readout stands at 201 at the hub reading and 67 at the dock reading, the two offices parted by the factor nine at every body"),
    ("185.12", "The electron at its two-dock web has the rank $67$ (its measure $79$ less the two struck (2,3) docks) and the dressing $3$ from its one in-leg, so $S_1=67\\cdot 3=201$."),
    ("185.12", "whose act readout is the measure of the link, the product of its two end dimensions"),
    ("186.1", "S_1=(79-12)\\cdot 1 = 67"),
    ("176.2", "the 79 closed records of the one link fall into 48 mixed records, 28 at two qutrit slots and 3 at two qubit slots, and no union of these sets stands at 27, 40 or 12, the weights of the three sets of the roster electron, so the electron is forced at no partition of the closed records"),
    # the muon and tau cells
    ("152.1", "27 + 2 \\cdot 40 + 6 + 4 = 117, \\qquad 27 + 3 \\cdot 40 + 0 + 8 = 155, \\qquad \\operatorname{rank}(\\mu) = 107, \\quad \\operatorname{rank}(\\tau) = 147"),
    ("152.1", "the struck docks at measures six and zero and the struck two-qubit leaf links at four and eight"),
    ("185.12", "S_1 = 201,\\ 963,\\ 3969,"),
    ("186.1", "$(117-10)\\cdot 9=963$"),
    ("186.1", "$(155-8)\\cdot 27=3969$"),
    ("186.1", "reaches its two tied hubs at three steps with both in-legs unpassed"),
    ("20.1", "reaches two hubs at the third step at once through two in-legs, a tie, both passed"),
    ("20.1", "p_{e}=(79-12)\\cdot 1=67,\\quad p_{\\mu}=(79-10)\\cdot 1=69,\\quad p_{\\tau}=(24-8)\\cdot 27=432,\\quad p_{p}=(319-18)\\cdot 9=2709,\\quad p_{n}=(319-18)\\cdot 27=8127"),
    ("20.1", "r_{e}=79,\\qquad r_{\\mu}=79,\\qquad r_{\\tau}=24,\\qquad r_{p}=319,\\qquad r_{n}=319"),
    ("195.1", "The muon stands at its web of fifteen links, the measure 117 and the dock set 10. The tau stands at its web of twenty links, the measure 155 and the dock set 8."),
    ("195.1", "hubs H; H2; in-legs h-a; H2-B; dressing 9; ports H-a3; w2-C; w1-b3"),
    ("195.1", "hubs H; H2; H3; in-legs h-a; H2-B; H3-C; dressing 27; ports H-a3; w1-b3; w2-c3"),
    ("189.1", "the third lepton stands at the web of the library and at no sink of the census at either hub-count reading"),
    ("189.1", "Open is the tau as a sink of the saturation: the tau of the library is a web the census reaches at no hub-count-two and no hub-count-three sink."),
    ("189.1", "The electron and the muon of the fixed library are each matched piece by piece by a sink of the saturation census"),
    ("189.1", "the first powers spread over ninety-two values from 95 to 3321"),
    ("154.1", "their dressings are 3 at 199 sinks, 9 at 88 and 27 at 13, and their first powers run from 369 to 4617 at 63 distinct values"),
    # 319 and the proton cell 2871
    ("176.2", "the support rank of the one qutrit-qubit link is the Born count over both ends, 319 as the sum over the 172 record letters"),
    ("175.3", "101 = 3 \\cdot 9 + 3 \\cdot 6 + 3 \\cdot 12 + 2 \\cdot 8 + 4"),
    ("175.3", "the center $J$ with a (3,3) leg to each quark, an up dock at each quark, the hub $H$ with a (3,4) leg from each quark and the pair, ten sites and twelve links"),
    ("195.1", "The proton stands at the sink of the bound-respecting reading with nine links, the measure 81 and the dock set 18, its entered link capped at the cap of its genre."),
    ("195.1", "dark links 6; hubs s7; dressing 9"),
    ("195.1", "in-legs s3-s7; s4-s7; s5-s7; dressing 27; ports s0-s3; s1-s4; s2-s5"),
    ("20.1", "The proton, entered along the dock at its first site, reaches the core and the hub at the second step through one in-leg passed and leaves two in-legs not passed."),
    ("186.1", "S_1=(319-0)\\cdot 9+(0-0)\\cdot 1=2871"),
    ("186.1", "every in-leg the read leaves unpassed dresses the entry segment, the dock set is struck at an uncapped entry alone, and a capped link beyond the first cap reads nothing"),
    ("186.1", "The test applied to each read is the equality of this sum with the derived cell of the body"),
    ("186.1", "meets every cell at 8 of its 36 readings, every meeting reading at the strike and the beyond records adopted here"),
    ("186.1", "Striking the dock set at a capped entry moves the proton to $2709$"),
    ("186.1", "moves the proton to $8613$"),
    ("186.1", "Dressing the arm segment, the passed in-leg at the stop dressed at every read, moves the deuteron by the factor three and the electron to $201$"),
    ("20.1", "the candidate cell two thousand eight hundred seventy-one at the proton"),
    ("20.1", "Over the grid of sixty readings, four of the register, three of the records and five of the path, no one reading meets the four derived cells at once"),
    ("20.1", "The proton is met at six readings, every one at the capped register of its entered link with nothing removed"),
    ("20.1", "so the passed records at the leptons and the passed leg at the proton are objects of two rules"),
    ("21.1", "the read stopping at the capped dock link dresses the geodesic leg that the candidate cell of the cap census needs undressed"),
    ("32.1", "but the chain is a reading at two steps, the typing of which web the measuring read enters and the bound at which its count stops, each still owing a derivation"),
    ("175.3", "no sink of the site-bounded closure, 44341 sinks at wave four, and no sink of the clean twin reads the proton figure of the comparison table at any reading tried, so the candidate is refuted at the bounded closures"),
    ("175.4", "asks whether the readout principle applied to such a capped web reproduces the two proton cells of the construction"),
    ("179.1", "S_{1} = 319\\cdot 3^{u}"),
    ("179.1", "Striking the typing to uncapped moves the ground office from 2871 to 747"),
    # the mass law
    ("185.9", "The mass is the order-two chain sum at the common support, the first power squared at the up quark and that square multiplied by one plus the order-two body term at the twelve bodies with a dark link."),
    ("195.1", "$\\mathrm{office}=\\mathrm{pair}+\\mathrm{closure}$"),
    ("195.1", "S_1=67,\\ 963,\\ 3969,\\ 2871,\\ 2871"),
    ("195.1", "at the proton and at the neutron through the deuteron the closure part is zero"),
    ("176.2", "the mass at the first plus second order $T_{1}(1 + G_{2}^{2})$"),
    ("176.2", "the two bodies entering at a (2,3) port carrying one $G_{2}$ each"),
    ("200.3", "so the exponent is one and the word reads at the order two and no lower"),
    ("176.1", "The proton in the electron unit stands outside the comparison by orders at every reading the construction has produced"),
    # faces and the dressing law (task 3)
    ("155.1", "The two readings of the mass law, the dock reading and the hub reading, are one law up to the unit absorbing one leg at every established result, the ratio form carries no reading"),
    ("155.1", "the mixed face the hub numerator over the dock unit"),
    ("155.1", "\\text{dock: } \\frac{M \\cdot 3^{k}}{67}, \\qquad \\text{hub: } \\frac{M \\cdot 3^{k+1}}{201}"),
    ("155.1", "The in-leg count is one function of the web alone, invariant under every reading"),
    ("155.1", "The strike rule that a dock is struck exactly when it is the probe attach site reproduces the struck set at no tabulated lepton result."),
    ("198.3", "The dressing of every body and of every joined content is a function of its hub in-legs alone, the product over the distinct hub in-legs of the far-end dimension, three at a qutrit in-leg and four at a ququart in-leg counted once"),
    ("185.6", "the product over the in-legs of the far-end dimension equals the dressing"),
    ("157.1", "with the electron at its dock-struck first power sixty-seven in place of the hub first power two hundred and one"),
    ("157.1", "At the electron the dock-spelled sum is the hub-spelled sum over the factor nine"),
    ("176.2", "\\Bigl(\\frac{963}{201}\\Bigr)^{2} = \\frac{103041}{4489}, \\qquad \\Bigl(\\frac{3969}{201}\\Bigr)^{2} = \\frac{1750329}{4489}"),
    ("21.1", "\\text{first power at the stop clause} & 67 & 621 & 432 & 8613"),
    ("181.5", "8613 = 319 \\cdot 27, \\qquad 25839 = 319 \\cdot 81, \\qquad 2871 = 319 \\cdot 9"),
    ("181.5", "top sink & the first reached top & 2871 & 8242641"),
]

FILE_ANCHORS = [
    ("parts/09-cell-3-3-strong.tex", "\\section{The proton-to-electron ratio: a candidate}"),
    ("parts/09-cell-3-3-strong.tex", "Theorem~176.1 reads the mass ratio at the second-order read, a candidate stated as a candidate"),
]

# ---------------------------------------------------------------- printed second-order objects
# (exact; a + b sqrt2 + c sqrt3 + d sqrt6 as rational strings) with their anchors
CLOSURES = {   # 195.1 closure parts (dock spelling of the electron), added to the pair part S^2
    "e":   (("13708610839/477757440000", "661/36450", "5977541/23887872000", "0"),
            ("195.1", "\\mathrm{closure}_{e}&=\\tfrac{13\\,708\\,610\\,839}{477\\,757\\,440\\,000}+\\tfrac{661}{36\\,450}\\sqrt{2}+\\tfrac{5\\,977\\,541}{23\\,887\\,872\\,000}\\sqrt{3}")),
    "mu":  (("27214907303/5898240000", "661/225", "6301957/294912000", "0"),
            ("195.1", "\\mathrm{closure}_{\\mu}&=\\tfrac{27\\,214\\,907\\,303}{5\\,898\\,240\\,000}+\\tfrac{661}{225}\\sqrt{2}+\\tfrac{6\\,301\\,957}{294\\,912\\,000}\\sqrt{3}")),
    "tau": (("40721203767/655360000", "1983/50", "6626373/32768000", "0"),
            ("195.1", "\\mathrm{closure}_{\\tau}&=\\tfrac{40\\,721\\,203\\,767}{655\\,360\\,000}+\\tfrac{1983}{50}\\sqrt{2}+\\tfrac{6\\,626\\,373}{32\\,768\\,000}\\sqrt{3}")),
    "p_massrow": (("732073/122880000", "0", "5069/1536000", "0"),
            ("O.1", "the closure zero at the act read and $\\tfrac{732073}{122880000} + \\tfrac{5069}{1536000}\\sqrt{3}$ at the mass-row weight")),
}
G2 = (("79/20736", "1/432", "7/3456", "1/6912"),
      ("176.2", "G_{2} = \\frac{79}{20736} + \\frac{1}{432}\\sqrt{2} + \\frac{7}{3456}\\sqrt{3} + \\frac{1}{6912}\\sqrt{6}"))

# ---------------------------------------------------------------- printed cells (task 2, the printed space)
# (value, label, anchor). Value is the first power as printed; '214.6' is irrational and given as a+c sqrt3.
PROTON_CELLS = [
    (2871, "186.1/195.1 derived cell; 20.1, 175.4 'candidate cell'", ("186.1", "S_1=(319-0)\\cdot 9+(0-0)\\cdot 1=2871")),
    (8613, "21.1 stop clause (hub numerator); 175.1 'candidate cell 8613'; 186.1 merge", ("175.1", "the candidate cell 8613 equals the first power at no sink of the bounded closure at any reading")),
    (2709, "20.1 adopted reading (strike 18); 186.1 strike at a capped entry", ("20.1", "p_{p}=(319-18)\\cdot 9=2709")),
    (2727, "175.1 seed, hub spelling", ("175.1", "the first power 2727 at the hub spelling and 101 at the dock spelling")),
    (101, "175.1 seed, dock spelling", ("175.1", "the first power 2727 at the hub spelling and 101 at the dock spelling")),
    (2187, "175.3 first catalogued sink, hub spelling", ("175.3", "reads 2187 at the hub spelling of the measure law")),
    (3108, "175.3 Born count over the open end sets of the junction legs", ("175.3", "the Born count over their open end sets is 3108")),
    (25839, "175.4 counted reading", ("175.4", "25839 at the counted reading and 8613 at the excluded reading")),
    (729, "180.1 share-face register 81 at dressing 9", ("180.1", "its first power 729 against 2871")),
    (747, "179.1 typing struck to uncapped", ("179.1", "Striking the typing to uncapped moves the ground office from 2871 to 747")),
    (1701, "195.1 mass-register row (63 x 27)", ("195.1", "the read BP: the class touching the ququart; the first power 1701")),
    (255, "185.12 proton star (rank 85, dressing 3)", ("185.12", "at the proton star (the rank $85$, the dressing $3$, the first power $255$)")),
    (85, "195.1 hadron-offices row, proton star", ("195.1", "the read d-ps of the hadron offices row: the first power 85")),
    (201, "185.12 site product of the proton", ("185.12", "and at the site product of the proton, $201$")),
    ("214.6", "214.6 S_p at the hub unit 201", ("214.6", "S_p=\\tfrac{460\\,867}{3072}+\\tfrac{25}{1536}\\sqrt{3}")),
]
ELECTRON_UNITS = [
    (67, "dock reading (rank, dressing 1)", ("152.1", "201 at the hub reading and 67 at the dock reading")),
    (201, "hub reading (rank x 3)", ("152.1", "201 at the hub reading and 67 at the dock reading")),
    (79, "measure, no strike, dressing 1 (175.3 control)", ("175.3", "The electron control parts at the readings striking no dock, 237 and 79")),
    (237, "measure, no strike, dressing 3 (175.3 control)", ("175.3", "The electron control parts at the readings striking no dock, 237 and 79")),
]
MUON_CELLS = [
    (963, "185.12/186.1/195.1 (107 x 9)", ("186.1", "$(117-10)\\cdot 9=963$")),
    (69, "20.1 adopted reading", ("20.1", "p_{\\mu}=(79-10)\\cdot 1=69")),
    (621, "21.1 stop clause", ("21.1", "\\text{first power at the stop clause} & 67 & 621 & 432 & 8613")),
]
TAU_CELLS = [
    (3969, "185.12/186.1/195.1 (147 x 27)", ("186.1", "$(155-8)\\cdot 27=3969$")),
    (432, "20.1 adopted reading; 21.1", ("20.1", "p_{\\tau}=(24-8)\\cdot 27=432")),
    (4185, "154.1 the tau sink of the census (star kind)", ("154.1", "the nearest of the star kind the tau sink itself at 4185")),
    (4212, "154.1 generation-step web", ("154.1", "9336401 over 13739282 at the first power 4212")),
    (3996, "154.1 nearest triangle-plus-one sink", ("154.1", "the nearest candidate of the triangle-plus-one kind is the sink at the first power 3996")),
]
# printed numbers deliberately left out of the printed space, with reasons (shown in the output for audit)
EXCLUDED = [
    ("116.1", "first powers 67, 71, 75, 79 / 963 ... 1071 / 3969 ... 4293 / 2871: bodies with serially inserted (2,2) links at the metric, not the bodies' own webs"),
    ("175.3", "10773, 133: the neutron reading of the seed (another body)"),
    ("186.1, O.1", "8127, 2873 and the neutron-through-deuteron 2871: neutron readings (round 2, SC-39/40)"),
    ("W.1", "1647, 1971: deuteron at the joined-web office formula (another body)"),
    ("176.2", "'the electron at its dock-struck rank 67 against 350': 350 is not identified as an electron first power (ambiguous)"),
    ("154.1", "the 300 three-hub sinks of the tau census: only the range 369 to 4617 and 63 distinct values are printed, not the values"),
    ("181.5", "the arm and heavy-sink cells 319, 8613, 25839, 2871 belong to quark bodies; used as a cross-check, not as proton readings"),
]

# ---------------------------------------------------------------- the electron web figure of 152.1 is parsed from the tex.
# The muon, tau and proton webs are built from these printed descriptions (masscells.py, build_webs):
WEB_NOTES = {
    "mu":  "core triangle A,B,C ((3,3) x3); hubs H (in-leg from A), H2 (in-leg from B), each with three (2,4) qubit legs and one "
           "(2,2) qubit pair (152.1); dock w2-C (2,3) and leaf w1-b3 (2,2), both struck (152.1 'struck docks at measures six "
           "and zero and the struck two-qubit leaf links at four and eight'; 195.1 row 'in-legs h-a; H2-B ... ports ... w2-C; w1-b3')",
    "tau": "core triangle A,B,C; hubs H, H2, H3 with in-legs from A, B, C; three (2,4) legs and one (2,2) pair each; leaves w1-b3, "
           "w2-c3 (2,2), both struck (152.1; 195.1 row 'in-legs h-a; H2-B; H3-C ... ports H-a3; w1-b3; w2-c3'); entered at w1 "
           "along the pendant pair (20.1, 186.1)",
    "p_seed": "175.3: centre J with a (3,3) leg to each quark Q1..Q3, a (2,3) dock n_i at each quark, hub H with a (3,4) leg from "
              "each quark, and the pair: two qubits f1, f2 on (2,4) legs at H joined by a (2,2) link (101 = 3*9 + 3*6 + 3*12 + 2*8 + 4)",
    "p_sink": "195.1 sink of record: quarks s3,s4,s5 with docks s0-s3, s1-s4, s2-s5 (2,3), in-legs s3-s7, s4-s7, s5-s7 (3,4) to "
              "the hub s7 and three (3,3) legs to the junction s6 (8 sites, 9 links, measure 81 = 18 + 36 + 27)",
}
