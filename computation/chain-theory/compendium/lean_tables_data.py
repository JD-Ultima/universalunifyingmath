"""Curated data for lean_tables.py (round 7: the 37 table-check rows and the 11 opaque-only rows).

Every quoted string is checked verbatim by lean_tables.py: each clause against the box body in compendium/tex
(whitespace-normalised, as in round 6), each check anchor against the `checks` statement of the named data module.

The rule is round 6's (09-lean-audit-at-scale.md §3 and §8): each box is split into clauses; a clause has a
counterpart when a tied declaration states its finite content under the library's own definitions of its words;
MATCH = every clause has a counterpart, PARTIAL = some clause does, NONE = no clause does. As applied here and in round 6,
a clause counts when its counterpart states some part of its finite content; PART_ONLY (at the end) lists the credited
clauses that are stated only in part, and lean_tables.py [10] prints the strict reading without them (REVIEW4-CT).

Counterpart kinds used here:
  ENTERED   the check evaluates the relation the clause asserts over entered per-instance values of the quantities the
            clause is about (round 6's ENTERED kind: LemmaB_zero_transfer evaluates the entered step triples). Counts
            as a counterpart.
  RECORDED  the check confirms that the entered list contains a string that records the clause's conclusion: a verdict
            (YES, EQUAL, MONOTONE, NONE, ...), a count written as a string ("5 of 5", "13"), or a per-instance verdict
            flag counted by a filter. The library defines none of the words of the key; the string is the claim.
            Not a counterpart under the rule; reported as the label-inclusive upper bound.
A RECORDED match needs a lexical link: the key or value carries a number or a content word of the clause.
Whole-box verdict cells (`cell.C-XXXX = YES`, a letter for the whole box) are not matched to clauses.
"""

E, R = "ENTERED", "RECORDED"
KV = lambda k, v: f'r.key == "{k}" && r.value == "{v}"'

# row -> list of (verbatim clause, [(declaration (last name component), check anchor, kind)], note)
ROWS = {
    "3.1": [
        ("every near basis direction whose transport is nonnull fixes exactly one far direction, the traceless part of its projector, and frees the two transverse to it, at all eleven nonnull directions", [], ""),
        ("while the one null direction frees the far end", [], ""),
        ("The washed read along $Z$ fixes the one far letter $Z$ and frees $X$ and $Y$", [],
         "entered in D525257_HTOR (OPEN-READ-OF-RECORD, cap read Z: fixed set Z, 1 fixed, 2 unfixed), which is tied to 57.1, not to 3.1"),
        ("the three open reads $X$, $XZ$ and $XZ^{2}$ fix none and free three, their far directions a trine", [],
         "HTOR rows X, XZ, XZ2: fixed NONE, 0 fixed, 3 unfixed; check_4 counts 3 with far_washed = NO; tied to 57.1"),
        ("the fixing is the same at every move of the identification group on the support",
         [("row_535360_NDFE", KV("eqfi.record.moves.covariance", "YES"), R)], "label 'record.moves.covariance'"),
        ("the far marginal is invariant on the state at every read", [], ""),
        ("On this reading the premise that the far end is fixed and freed direction by direction is a tautology of the definition of transport to the far end, drafted as Lean text with its Lean proof still to be written", [],
         "the box itself says the Lean proof is not written; the table enters eqfi.lean.remet = YES"),
    ],
    "11.1": [
        ("the on-shell condition of the literature, residue one at the physical pole and the vertex at zero momentum transfer, is the unit lemma of the compendium, met as an exact identity by the construction of the readout unit and of the read map", [], ""),
        ("the readout unit of the one link is an exact element at every derived cell, the trace form with the centre valued one, the identity reading one at the register dimension, the order step one tick at the pair state, the bare two-point function and the bare coupling both one",
         [("row_535412_OLUS", KV("cell.C-UNIT", "YES"), R)], "cell named UNIT"),
        ("the read map product of the one-link record web is the identity on its labels with one pole at the unit, its residue the derived pole residue",
         [("row_535412_OLUS", KV("cell.C-POLE", "YES"), R)], "cell named POLE"),
        ("the residue is normalised at the readout unit at both registers with the vertex at zero displacement",
         [("row_535412_OLUS", KV("cell.C-NORM", "YES"), R)], "cell named NORM"),
        ("the dressed residue being the coupling the scheme defines and never the normalisation it fixes", [], ""),
    ],
    "13.1": [
        ("The record phase group is finite cyclic at every link type and every end",
         [("row_535680_RPHG", KV("eqpq.letter", "FINITE"), R)], "letter FINITE; 'cyclic' not recorded"),
        ("the phases of the admitted matchings, the entries and eigenvalues of the admitted cap bases and the stabiliser of the unit are roots of unity of the ring of the twenty-fourth root at every end", [], ""),
        ("generating cyclic groups of order four at the qubit end, twelve at the qutrit end and eight at the ququart end, the (2,3) link at order twelve", [], "no order is entered"),
        ("every derived holonomy, the Berry result, the classical-limit sequences and the Chern result, lies inside its group", [], ""),
        ("the per-vector phases act by conjugation on every holonomy, so that the characteristic polynomial is free of every phase variable and the phase entry of the definition census is typed a gauge", [], ""),
        ("The reading on which an antilinear holonomy is read as a linear matrix is printed beside it.", [], ""),
    ],
    "28.1": [
        ("This chapter establishes that measurement is the wash (decoherence) at one act, and reads the probe at the register.", [], ""),
        ("At every one of the two hundred twenty-eight one-act reads of the network list the readout at the open end, computed as exact arithmetic in the cyclotomic field, equals the reduced state of the pair at the far end at every direction", [],
         "the 14 entered WAMR rows sum to 228 acts (own check), which no check of WAMR states; no readout or reduced state is entered"),
        ("with the Born square the unit at every read", [], ""),
        ("and is unmoved under every change confined to a washed end",
         [("row_526019_WAMR", "rows.all (fun r => r.wash_moved_acts == 0) = true", R)],
         "an entered per-body count of acts moved by the wash, 0 at the 14 bodies: the conclusion as entered"),
        ("so measurement is the wash at one act and the washed ends are the apparatus", [], ""),
    ],
    "28.2": [
        ("The star-pair and one-link apparatus webs are built at every admitted reading for every measurement of the item.",
         [("row_527964_MWSP", KV("eqsp.built", "5 of 5"), R), ("row_527964_MWSP", KV("eqsp.readings.built", "10 of 10"), R),
          ("row_527964_MWSP", KV("eqsp.letter", "BUILT"), R)], "labels 'built'"),
    ],
    "28.5": [
        ("The classical limit is the zero-density register", [], ""),
        ("at every admitted web the sign-record density counts the pairs of reads commuting only up to a sign against the register size", [], ""),
        ("the webs at zero density are listed at every reading, the gluon and the Higgs bodies among them", [], ""),
        ("the semiclassical expansion is the expansion in that density", [], ""),
        ("the first power of every sink re-met as the sum of its terms grouped by the sign-pair count", [], ""),
    ],
    "57.1": [
        ("At the first body at which the twice-spin read from the open qubit ends alone parted from the name, the docked up quark at the reading that keeps the nonzero labels", [],
         "PQCD enters equs.record.reading = DOCK, but no check reads it"),
        ("a web of three sites and two links at the genre (2,3) with one open ket end and one open bra end of dimension two and no open qutrit or ququart end", [], ""),
        ("the label content re-built at a second route from its open ends alone, the weight route with no pairwise Clebsch-Gordan rule, reads zero and two at four states with the multiplicities one and one", [], ""),
        ("equal to the rule at that body at every one of its three readings, the low at zero against zero, the nonzero at two against two and the maximal at two against two", [],
         "PQCD's checked counts are '8 of 8' and '5 of 8', not three readings"),
        ("the arrow spaces of its two links are printed at the dimension six each, the ket end tensor the bra end, with both weight faces as the catalogue carries them", [], ""),
        ("the body has no ququart end whose frame could part", [], ""),
        ("The construction of the open-end label therefore stands re-met at its first differing body by a second route", [], ""),
        ("the parting the previous open result named is a parting from the name alone, the body up at the twice-spin one by name, which differs from the second route at every reading and is an input to no cell", [], ""),
        ("the deciding table is re-derived in Lean", [], "both tied tables are entered literal lists; none is derived"),
    ],
    "65.1": [
        ("This chapter establishes that the gauge group is a corollary of one construction.", [], ""),
        ("The stabiliser of the Hermitian unit inside the isometries of the split, spanned by the adjoint matrices of the derived directions, applied at the three end kinds of the listed webs, yields the three factors as exact objects", [],
         "GSTB enters cells C-CONS, C-ENDS, C-GRPS, C-MTCH = YES; no factor, dimension or rank"),
        ("at the hub-leg qubit copies the register of dimension two gives dimension three and rank one with the lift centre of order two, the weak factor", [], ""),
        ("at the qutrit ends open at the hub leg the register of dimension three gives dimension eight and rank two with the lift centre of order three, the colour factor", [], ""),
        ("at the dock qubit ends the $u(1)$ of the unit has the slice representative equal to the centre generator of the group result", [], ""),
        ("while the stabiliser at the dock register is the rotation group, spatial and compared with no factor", [], ""),
        ("the kernel of the covering is re-met from the three register centres", [], ""),
        ("The reading beside, the stabiliser inside the isometries whole, parts from its factor first at the hub-leg copy.", [], ""),
    ],
    "65.2": [
        ("At the ququart hub ends, the hub pair read as one qudit of dimension four by the Axiom, the same construction yields the Lie algebra of dimension fifteen and rank three with the lift $SU(4)$ and the centre cyclic of order four",
         [("row_530426_UB4H", KV("cell.C-ALGB", "YES"), R), ("row_530426_UB4H", KV("cell.C-HUBE", "YES"), R)],
         "cells named ALGB and HUBE; dimension 15, rank 3 and order 4 are not entered"),
        ("the adopted carrier phase datum, the fourth root, is a generator of that centre, so the datum is derived at the hub centre",
         [("row_530426_UB4H", KV("cell.C-DATM", "YES"), R)], "cell named DATM"),
        ("the hub-pair weights at the two qubit copies with the $u(1)$ weight at the dock slice equal the doublet of the Higgs entry by entry, so the Higgs representation is an output of the hub-pair weights",
         [("row_530426_UB4H", KV("cell.C-HIGG", "YES"), R)], "cell named HIGG"),
        ("the readings beside parting at the centre or at a weight", [], ""),
    ],
    "145.1": [
        ("The excluded count between two sinks of the saturation placed along a chain of the construction is a definite function of the separation at every separation the chain admits, with one sign",
         [("row_531942_CASM", KV("cas.SINK.ARC.COUNT.QUTRIT.s0-s1.N10.definite", "YES"), R),
          ("row_531942_CASM", KV("eqca.first.undefined", "NONE"), R)],
         "per-separation labels 'definite' and 'one.sign' at N = 2, 3, 4, 10, 11, 12 of one chain"),
        ("so the Casimir effect is a composed cell of the second-law and chain results", [], ""),
        ("at the standing readings, the chain, the web between and the count face, the excluded records are defined at all forty-four chains", [],
         "the detail entries cover one chain (QUTRIT s0-s1) at six separations and stop at the 60th entry; 44 is not entered"),
        ("the nonzero differences carry one sign at every one",
         [("row_531942_CASM", KV("eqca.first.sign.change", "NONE"), R),
          ("row_531942_CASM", KV("cas.SINK.ARC.COUNT.QUTRIT.s0-s1.N10.one.sign", "YES"), R)], "labels 'one.sign', 'first.sign.change'"),
    ],
    "181.5": [
        ("every arm of the two catalogued stars and every heavy sink is an office by the body-alone rule", [], ""),
        ("the six arms at the first power 319 and the office 101761 each", [], ""),
        ("the up and the down equal with their ratio one", [], "QLHL's checked up and down counts are 9 and 124 (another object)"),
        ("the heavy sink offices the cap squared at the powers of the dressing", [], ""),
        ("every ratio among them an exact rational cell", [], ""),
        ("each star re-met as the sum of its arms by the cross-term theorem", [], ""),
    ],
    "186.1": [
        ("The measuring read of a body is cut into segments at every capped link and at every hub it meets", [], ""),
        ("every in-leg the read leaves unpassed dresses the entry segment, the dock set is struck at an uncapped entry alone, and a capped link beyond the first cap reads nothing", [], ""),
        ("Under that reading one formula, the sum over the segments of the records less the strike multiplied by the dressing, returns the derived cells at the five bodies and at the neutron through the deuteron", [], ""),
        ("restates the one rule of the strike-rule result at every read with no body-specific term", [], ""),
    ],
    "200.6": [
        ("Every Higgs and carrier office of the compendium is an exact structure by the body-alone rule at the roster words", [], ""),
        ("the closed-genre objects and the free genres are read at their own sinks", [], ""),
        ("every office is an exact cell at every reading the results admit", [], ""),
        ("every ratio among them is an exact rational or ring cell, at all 52 ratios", [], ""),
        ("the Higgs faces close 2 of 2, the carriers 6 of 6", [], "HGHL's counts are 28 of 28, 4 of 4, 48 of 48, 0 of 4, 8 of 28"),
        ("no sink the rule does not close and no ratio with no second read is found", [], ""),
    ],
    "203.5": [
        ("This theorem, closing the open result of the previous releases at the literal closure reading", [], ""),
        ("is stated once, as Theorem~170.1 of the chapter Beta decay and the boson widths, and this chapter cites it there without restating it", [],
         "a pointer to 170.1, an opaque-only row; ENPC records census counts 3, 8, 35"),
    ],
    "217.1": [
        ("The object of this chapter is the boundary count (the entanglement entropy across a cut) of an interface at its counting face.", [], ""),
        ("The theorem states that this count is non-decreasing under every one-act read whose open end lies at the boundary, and under every one-act read beside, at each of the thirteen bodies and the neutral body",
         [("row_525335_SLBC", "rows.all (fun r => r.W_target_site - r.W_source_site == r.W_diff_site) = true", E),
          ("row_525335_SLBC", "rows.all (fun r => decide (0 <= r.W_diff_site)) = true", E)],
         "W_target_site >= W_source_site at every entered read; the 14 bodies are entered, not checked"),
        ("at all 228 one-act reads of the fourteen bodies (two per link at 114 links, 78 with an open far end, 150 washed)",
         [("row_525335_SLBC", "rows.length = 228", E)],
         "228 entered reads; 114 links, 78 OPEN, 150 WASHED and 14 bodies follow from the entered rows (own check) but are not checked"),
        ("the count rises or stays at the site bound, at the free-direction bound and at the directions read",
         [("row_525335_SLBC", "rows.all (fun r => decide (0 <= r.W_diff_site)) = true", E),
          ("row_525335_SLBC", "rows.all (fun r => decide (0 <= r.R_diff_site)) = true", R)],
         "site bound: W checked with its consistency (W is constant in all 228 rows, own check); R_diff_site >= 0 is an entered difference whose relation to the counts is not checked; the free-direction fields are not checked"),
        ("no decreasing act is found",
         [("row_525335_SLBC", "rows.all (fun r => decide (0 <= r.W_diff_site)) = true", E)], ""),
        ("the running count along the order coordinate is non-decreasing", [], "run_site is entered, not checked"),
        ("The result is a corollary of the record clause (a read is a further record and no read removes one) and of the definition of an open end", [], ""),
        ("so the second law is a theorem at the counting face with the arrow of the order-coordinate theorem as its direction, the words entropy and arrow entering by name", [], ""),
    ],
    "217.2": [
        ("The parent theorem fixes the count: at every one of the 60 deciding cells of the horizon composite and the universe web the count at the counting face equals the site bound of the valence-law question at 60 of 60", [],
         "SLBC's rows are the 228 reads of 14 bodies; no horizon composite or universe web is entered"),
        ("and the free-direction bound at 0 of 60", [], ""),
        ("so the first premise of the Einstein-shaped equation is a theorem of the count at the composite", [], ""),
    ],
    "219.1": [
        ("Sixteen thermal names of the literature", [("row_527268_ENSB", KV("eqen.entries", "16"), R)], "count label"),
        ("are read as readings taken at the caps of the joined networks with a bath body", [], ""),
        ("the three ensembles, the partition sums, the bath readout and the zeroth and third laws are exact objects at the joined networks", [], ""),
        ("the three distributions follow the sign of the own-site rotation with every parting named, and no Gibbs ratio is asserted", [], ""),
        ("Thirteen names are derived and three are glossed, none open",
         [("row_527268_ENSB", KV("eqen.derived", "13"), R), ("row_527268_ENSB", KV("eqen.open", "0"), R)],
         "eqen.gloss = 3 is entered but not checked"),
        ("a second reading at zero freedom, reading every derived name against the textbook object it stands for, typed five as proof sketches, where the textbook definition is instantiated at the established object, and eleven as namings", [], ""),
        ("both typings printed and the deciding letter the same at both", [], ""),
    ],
    "220.2": [
        ("Restricted to the forward moves (182940 forward edges of 185496 at the large truncated closure, the 2556 inverse edges excluded)", [], ""),
        ("W is monotone along every forward edge",
         [("row_530247_TOP2", KV("eqtf.mono", "MONOTONE"), R), ("row_530247_TOP2", KV("CGB2.eqtf.mono", "MONOTONE"), R)], ""),
        ("strictly increasing at the forward pairing edges",
         [("row_530247_TOP2", KV("eqtf.strict", "STRICTLY-INCREASING"), R)], ""),
        ("and falling at the inverse edges alone", [], ""),
        ("the class changes along the forward edges at 9377 edges at 5 of 12 kinds",
         [("row_530247_TOP2", KV("CGB2.eqtf.transition.kinds.present", "5 of 12"), R)], "9377 is in PHTR (tied to 221.1), not here"),
        ("the first from the charged class to the dark class at the first wave by the attach move at the genre (2,3)", [], ""),
        ("and no web changes its determined class along a forward path",
         [("row_530247_TOP2", KV("eqtf.clss", "FIXED"), R), ("row_530247_TOP2", KV("eqtf.first.change", "NONE"), R)], ""),
        ("so the order parameter is the pair of the class and W at every sink, a theorem at the forward saturation", [], ""),
    ],
    "220.3": [
        ("At every forward edge of the twelve readable saturations (322516 webs, 138129 sinks, 570349 forward edges)", [],
         "DST2's keys name 12 saturations (own count); no count is entered"),
        ("the pair the move joins is, by the Axiom, a maximally entangled pair of its common support", [], ""),
        ("so its fidelity to the flat pair is one exactly at every link", [], ""),
        ("monotone along every forward path", [("row_532023_DST2", KV("eqdj.CGB2.mono.cell", "YES"), R)], "'mono.cell = YES' at each of the 12"),
        ("every sink at full fidelity", [("row_532023_DST2", KV("eqdj.CGB2.sink.cell", "YES"), R)], "'sink.cell = YES' at each of the 12"),
        ("the distillation protocols of the literature, each a rise of the fidelity from below full, are located as no path of the saturation", [], ""),
    ],
    "221.1": [
        ("The order parameter of the preceding chapter, the pair of the class of a sink and the register count W, is re-met along every forward edge of the large closure (101973 webs, 44341 sinks, 182940 forward edges) at the per-web and order tables of the forward saturation", [], ""),
        ("W strictly increasing along the forward edges", [("row_531680_PHTR", KV("ordp.strict.CGB2", "STRICTLY-INCREASING"), R)], ""),
        ("and tabulated beside at the proton seeded saturation (218340 webs), the bound-respecting seeded saturation (93481 webs, W unread at its table) and the library seeds (64 seeds at 1203 sink entries)", [], ""),
        ("at the large closure the class changes at 9377 forward edges (9270 at the fourth wave, 107 below it)",
         [("row_531680_PHTR", KV("eqpt.changes", "9377"), R)], "9270 and 107 are not entered"),
        ("the critical edge running from a charged parent at W two to a dark child at W four at the first wave by the attach move at the genre (2,3), so the critical W is four with the charged and the dark classes on either side",
         [("row_531680_PHTR", KV("eqpt.crit.classes", "CHARGED>DARK"), R)], "the classes only; W two and four are not entered"),
        ("at the proton seed the critical W is 93312 (from the class beyond two to the neutral class), at the library seeds 5184 (charged to dark), and at the bound-respecting seed the critical measure is 21, each closure its own", [], ""),
        ("fourteen bodies rebuilt from the caps stand eleven at a class boundary and none below the critical W of the large truncated closure", [], ""),
        ("the ratio two of the critical edge equals the twirl at the support at the dimension two alone", [], ""),
        ("The identification with the phase transition of the literature is a correspondence by name, and every count is a count of a truncated closure.", [], ""),
    ],
    "223.2": [
        ("with W taken as the entropy, the entropy production is strictly increasing at every one of the 182940 forward edges of the large truncated closure",
         [("row_530584_ENT2", KV("entw.strict", "STRICTLY-INCREASING"), R), ("row_530584_ENT2", KV("entw.sample.count", "182940"), R)], ""),
        ("the second law holds at the strict reading along every forward path and the name entropy is keyed",
         [("row_530584_ENT2", KV("entw.path", "STRICTLY-INCREASING"), R)], ""),
        ("the boundary count of the second law, read beside at the same edges, is typed two objects and not the entropy",
         [("row_530584_ENT2", KV("entw.slaw.typing", "TWO-OBJECTS"), R)], ""),
        ("the fluctuation relation has 179 defined pairs of 2989 at the sinks graph, the two candidate functions meeting the ratio at 0 and at 6 of them", [], "179 and 2989 are not entered"),
        ("so the fluctuation name stays silent, a count of a truncated closure with no ensemble",
         [("row_530584_ENT2", KV("fluc.function.name", "NONE"), R)], ""),
        ("the asymmetric Onsager pairs are exactly the pairs at a non-zero production difference at every one of the 43",
         [("row_530584_ENT2", KV("onsg.setequal", "EQUAL"), R)],
         "the same table enters and checks onsg.asym.pairs = \"24\"; the box prints 43"),
        ("so the Onsager name is keyed", [], ""),
    ],
    "224.1": [
        ("A chain of two, three or four copies of the electron body joined at the dock by the pairing move, the join graph a path or a cycle, is one network", [], ""),
        ("whose direct record sum is an exact element at every one of the 5 arrays", [("row_527373_LATW", KV("eqlw.arrays.built", "5 of 5"), R)], ""),
        ("every record graded five ways", [], "eqlw.gradings.closed = \"5 of 5\" is entered but not checked"),
        ("the readout is periodic under the cycle shift at 2 of 2 cycles at four tests",
         [("row_527373_LATW", KV("eqlw.cycles.periodic", "2 of 2"), R), ("row_527373_LATW", KV("eqlw.letter", "PERIODIC"), R)], "'four tests' not entered"),
        ("so the lattice and the Bloch theorem are corollaries by name", [], ""),
    ],
    "224.2": [
        ("The periodic network is admitted by the grammar at both cells of the item, the library star with its dock", [], ""),
        ("rings of N docked stars at three dock kinds, every join link a (2,2) pairing link between two docks, the translation an automorphism", [], ""),
        ("54 networks built, 17 periodic, 29 admitted", [], ""),
        ("the first link the dock pairing refuses is the like-kind link of the chain of three", [],
         "LWE3 checks STAR.first.refused = NONE and QUTRIT.first.refused = NONE (keys not defined; not matched)"),
        ("the 17 periodic networks are absent from both truncated closures", [], ""),
    ],
    "224.3": [
        ("The periodic sink exists at a named N with its cell an output", [], ""),
        ("at the library star the sink of the ring of N docked stars under the pairing move at the docks is periodic at the even N two, four and six and refused at the odd N three and five (the charge carried at one open dock)", [], ""),
        ("its cell the lattice cell with its free dock paired at the half turn (sites 11, measure 91)", [], ""),
        ("at the two-dock qutrit the ring is its own sink at every N from two to six with the lattice limit cell exactly (sites 3, measure 16)",
         [("row_530477_LWE3", KV("eqlw.QUTRIT.N3.of.record", "s0-s1.RING"), R), ("row_530477_LWE3", KV("eqlw.QUTRIT.N4.of.record", "s0-s1.RING"), R)],
         "the ring recorded as the sink at N = 3 and 4 only"),
        ("40 saturations terminated at 82 sinks, the parity law at every one", [], ""),
    ],
    "225.1": [
        ("At the periodic networks of the preceding chapter, the rings of N docked stars at the three dock kinds and of N docked qutrits, re-built by the same join at N two to six", [], ""),
        ("every link is a register of the step generator tiled by the translation", [], ""),
        ("the cycle-order spectrum of the ring is the spectrum of the cell repeated at every translation eigenvalue",
         [("row_530533_ANDB", KV("web.STAR.n1-n2.RING.2.tiled.spectrum.eq", "YES"), R),
          ("row_530533_ANDB", KV("web.STAR.n1-n2.RING.2.spectrum.translation.equal", "YES"), R)], "per-ring labels"),
        ("20 networks built, 20 periodic, 20 tiled, the spectrum equal to the tiling and translation-equal at every register at every one",
         [("row_530533_ANDB", KV("web.QUTRIT.s0-s1.RING.3.tiled.spectrum.eq", "YES"), R)],
         "per-ring labels at 18 rings; the list stops at the 60th entry inside QUTRIT RING.4; RING.5 and RING.6 of the qutrit are absent; 20 is not entered"),
        ("with no parting N or register at either cell",
         [("row_530533_ANDB", KV("eqbn.STAR.first.parting", "NONE"), R), ("row_530533_ANDB", KV("eqbn.QUTRIT.first.parting", "NONE"), R)], ""),
        ("The library star carries 11 registers and the two-dock qutrit 2, the join link the (2,2) register at every ring", [], ""),
        ("the band at the star is the set of orders one to five at every N with the multiplicity per cell constant across N, at the qutrit the set of orders one, two and four likewise", [], ""),
        ("the limit per cell equals the limit table of the lattice result", [], ""),
        ("and the dispersion reading is flat at 10 of 10 tables, with no value and no material constant", [], ""),
    ],
    "226.1": [
        ("At every ring of the lattice (the 12 rings of the lattice tables and 16 more built at N five to eight)", [],
         "the entered per-ring keys are 4 cells x N 2..8 = 28 (own count), consistent with 12 + 16"),
        ("the N join links are N (2,2) pairing links between dock-labelled pendant qubits of consecutive copies, disjoint and carried to the next by the translation", [], ""),
        ("so the chain of N pairs is admitted at all 28 rings with no parting N or register",
         [("row_530610_PCHS", KV("rings.face.A", "28 of 28"), R), ("row_530610_PCHS", KV("rings.face.B.first", "NONE"), R)], ""),
        ("every pair register is 4 with the pair ket a fixed point of the step", [], ""),
        ("the chain step is the product of the N pair steps at order two with the translation commuting with it", [], ""),
        ("the ring holonomy at the cycle is trivial at the letters Z and X everywhere and at the letter XZ trivial at every even N and nontrivial at every odd N", [], ""),
        ("the photon chain precedent is re-met at the same pairing move", [], ""),
    ],
    "230.1": [
        ("Every ready name of the computation vocabulary is a corollary of the text or of a theorem, or a gloss named beside an object the construction has, and every gated one names its gate",
         [("row_527242_PRCM", KV("gate.all", "MET"), R)], "label 'gate.all'"),
        ("The bit is the qubit pair ket of the Axiom; the gate is the act at one link; the channel is the wash at one act; the measurement is the wash at a cap; the classical register is the washed record; teleportation is the read passed along the pair together with its record; the Schmidt rank is the common support and the purification the Born face; the cost of a computation is the record count of a read", [], ""),
        ("the circuit is the closure of the act from a root, with the saturated networks its halting states and the stop condition of the nuclear roster its halting rule", [], ""),
        ("and the Turing machine stands open at the universality question", [], ""),
        ("Of the thirty-six names (fifteen ready, twenty-one gated) nine are derived, fourteen glossed and thirteen open at their named gates",
         [("row_527242_PRCM", KV("typing.total", "36"), R), ("row_527242_PRCM", KV("typing.rows", "36 of 36"), R)],
         "9, 14, 13, '15 of 15', '21 of 21' are entered but not checked"),
        ("the cost of a read is printed at 18 of 18 networks with the sum equal to the orders and at twelve closures (six stopped, six at the web cap)", [], ""),
        ("The words bit, gate, circuit, channel, teleportation, Turing machine and complexity class enter by name.", [], ""),
    ],
    "234.1": [
        ("At the sixteen networks of the item and at the sixty-four library networks the records of every link are counted by the closed count law at the slot polarity by order zero to three and by kind", [],
         "the counts are entered per network and polarity; no count law is in Lean"),
        ("the orders and the kinds summing to the records at every link and network",
         [("row_525516_NCRE", "rows.all (fun r => r.order0 + r.order1 + r.order2 + r.order3 == r.records) = true", E),
          ("row_525516_NCRE", "rows.all (fun r => r.same_site + r.one_link + r.cross_link == r.records) = true", E)],
         "at every entered network-polarity row (246); no per-link row is entered"),
        ("(64 of 64 and 18 of 18 defined)",
         [("row_525516_NCRE", 'rows.filter (fun r => r.group == "KIT" && r.polarity == "SLOT" && r.defined == "YES")).length = 64', R)],
         "a count of per-row 'defined = YES' flags; 18 follows from the entered rows (13 BODY + 3 NEUTRAL + 2 STAR, own check) but is not checked"),
        ("the record sums under the lowering at the depths zero, one and two grow with the ratio of successive sums non-increasing", [], ""),
        ("the finiteness beyond the second depth open", [], ""),
        ("and the growth with the order terminates at the third order at every depth", [], "order3 = 0 in all 246 entered rows (own check), not checked"),
        ("so the cost of a read at a network at an order is the record count at that order", [], ""),
    ],
    "234.3": [
        ("At every one of the 5707 applications every class of the records at the links touching no move site lands whole in the class of the same name at the result", [], ""),
        ("the kind of every one of the 270584 untouched links unmoved, with no parting",
         [("row_528475_DURA", KV("eqau.first.parting", "NONE"), R)], "label 'first.parting'"),
        ("so the adiabatic theorem is derived at the untouched records by a sketch of five lines from the locality of the kind rule, the move action on the grade and the placement classes", [], ""),
        ("the slow-change clause fenced from the sketch", [], ""),
    ],
    "235.1": [
        ("The third clause of the read line is a corollary of the Axiom word maximally entangled with the first two clauses", [], ""),
        ("the far set is an orthogonal basis of the far support at every genre and orientation", [], ""),
        ("and the orientation sign is needed at no read", [], ""),
    ],
    "235.2": [
        ("At the 341 closed paths of the deciding set of the library the transport composed around the path at the identity matching is an exact unitary on the register of the starting end", [], ""),
        ("defined in the ring of the twenty-fourth root of unity", [], ""),
        ("nontrivial at the 333 odd paths and trivial at the 8 even ones", [], ""),
        ("the Berry phase by name is the eigenphases at the identity matching, a gloss", [], ""),
    ],
    "235.3": [
        ("At every one of the 1692 registers of the step-generator result the step generator equals the ring quotient multiplied by the cycle-index operator",
         [("row_529015_YOGC", KV("dict.number.operator", "1692 of 1692"), R)], ""),
        ("a number operator with an equally spaced spectrum whose spacing is the quotient of that register",
         [("row_529015_YOGC", KV("dict.spacing.equal", "1692 of 1692"), R)], ""),
        ("and whose level count is the cycle order", [], "dict.levels.greatest = 5 is checked; no cycle order is entered"),
    ],
    "235.4": [
        ("At every carrier carrying a qutrit-qubit register the carrier phase datum equals the primitive fourth roots of unity in the holonomy spectrum of the transport-generator four-cycle of that register",
         [("row_530067_HRCP", KV("eqcr.WORD.PRIMITIVE-ROOT.carriers.equal", "3 of 3"), R), ("row_530067_HRCP", KV("eqcr.reading.record", "PRIMITIVE-ROOT"), R)],
         "the same table checks COMPOSED and PHASE-CONTENT readings at 0 of 2 and 0 of 3, not printed in the box"),
        ("the four-cycle present at the W, the Z and the photon words and absent at the gluon", [],
         "the table's 'NO-CYCLE' entries are at r-HIG (the Higgs word), not at the gluon"),
    ],
    "D.3": [
        ("The hub step is one move on the pair state", [], ""),
        ("the unit of the hub, attached as one body at the site of the step by its leg and its odd leg with the far ends open, re-makes every derived hub step",
         [("row_532499_HPK2", KV("eqhp.all.equal", "16 of 16"), R)], ""),
        ("the sixteen deciding webs standing at their grammar, dressing and direct sum",
         [("row_532499_HPK2", KV("eqhp.webs", "16"), R), ("row_532499_HPK2", KV("eqhp.grammar.equal", "16 of 16"), R),
          ("row_532499_HPK2", KV("eqhp.dressing.equal", "16 of 16"), R), ("row_532499_HPK2", KV("eqhp.directsum.equal", "16 of 16"), R),
          ("row_526261_RLRO", 'rows.filter (fun r => r.equal == "YES")).length = 16', R)],
         "RLRO: 16 entered webs (13 BODY, 3 BESIDE), each flagged equal = YES"),
    ],
    "D.4": [
        ("The pair state of the one qutrit-qubit link is closed at both orientations", [], ""),
        ("opened by the first act at the qutrit end with the far qubit end open and at the qubit end with the far qutrit end closed", [], ""),
        ("every listed body with its conjugate at every end and every arrow reversed is a web of the same grammar and the same direct sum, the law of the root", [],
         "RLRO's per-web flags do not name the conjugate; ATSR records readings R1-R3 at sinks (14), not the root"),
    ],
    "H.3": [
        ("Every seeded saturation terminates under the valence bound and the site bound", [],
         "ERMT records 14 demonstration runs, 7 'STAYS-BOUNDED' and 7 'GROWS-TO-THE-DEMONSTRATION-CAP' (max sites 40 to 42); LSRC records loop counts 0"),
        ("since the moves map the webs of at most $S$ sites into themselves and such webs are finitely many", [], "a finiteness argument; no Lean statement of it"),
    ],
    "I.7.1": [
        ("Decoherence is the wash at a sub-web", [], ""),
        ("at every tabulated sub-web carrying a shared record, 1753 of the 1755 sub-webs in 85 equivalence classes, the wash of every shared record at the complement, the act of the no-signalling result summed over its outcomes with the complement copy discarded and read as a channel on the state of the sub-web, is a twirl of the shared end", [], ""),
        ("a completely positive, unital, trace-preserving and idempotent map at every combination of the site dimension and the complement support", [], ""),
        ("the off-diagonal cells of the reduced state in the standing cap basis reach zero at the first record at every equivalence class, exact at a finite count and never asymptotic", [], ""),
        ("the record basis of the act is the pointer basis the wash leaves fixed at every dimension", [], ""),
    ],
    # ------------------------------------------------------------------ the 11 opaque-only rows
    "62.1": [
        ("the order parameter (the order parameter W), the product over the links of a web of the dimension of the common support, is preserved by the group moves, the step and the holonomy, at every gate edge by construction, a gate changing no dimension", [], ""),
        ("and raised by every transfer (the transfer), the attach, the cut, the hub-forming, the hub steps and the joins, at every read record edge with none non-rising, the increment per kind the ratio census", [], ""),
        ("the reverse-marked moves lower it at every read edge", [], ""),
        ("so every kind is partitioned into preserving, raising and lowering by the data alone, with no mixed kind", [], ""),
        ("every count is a count over the enumerated part carrying the caveat of its closure, and the theorem rests on the weakest of its parents, an open result among them", [], "the box says it rests on an open result"),
    ],
    "64.1": [
        ("Every readout of the one qutrit-qubit link is invariant under every matching of the identification group as a corollary of that clause", [], ""),
        ("the group is re-met as exact objects, the fourth-root monomial groups at the orders four, thirty-two, three hundred eighty-four and six thousand one hundred forty-four", [], ""),
        ("the whole group at the Lie dimensions four and nine", [], ""),
        ("the identification algebra at the six link types with the one readable type carrying the acting centre", [], "at (2,3) this is restated by AGUQ, tied to section 9, not to 64.1"),
        ("the one hundred seventy-two records are rebuilt and every readout is invariant on the state under every matching, every Lie generator, the Cayley unitaries, the chain unitaries and the Fourier elements", [], ""),
        ("The group is named as the gauge group of the Standard Model by correspondence.", [], ""),
    ],
    "64.2": [
        ("the face group of the qutrit-qubit tether, the image of the local unitaries $U(2)$ x $U(3)$ on the link, is the connected determinant slice $S(U(2)$ x $U(3))$ with Lie algebra of dimension twelve, derived algebra of dimension eleven and a one-dimensional centre spanned by the generator of the identification-algebra result", [], ""),
        ("isomorphic to $(SU(2)$ x $SU(3)$ x $U(1))/\\mathbb{Z}_6$ by an explicit covering whose kernel is the cyclic group of order six", [], ""),
        ("the uniqueness cell at the dimension twelve and the rank four picks the qutrit-qubit type at both", [], ""),
        ("the sixteen-row law of the dart holds at every row", [], ""),
    ],
    "73.1": [
        ("This chapter establishes that the isospin sign is the side of the acting site (the acting site) at every hub reached through its qutrit", [], ""),
        ("and that the slot rule is not the sink orientation read at the cut end", [], ""),
        ("every hub the carry enters through its qutrit leg carries the far side at its odd-leg qubit copy, equal to the letter of the network list there, at all five hubs of the ladder so entered and at every hub of the three lepton-seeded sinks", [], ""),
        ("the tau's second hub is entered through its qubit leg, its odd-leg copy reading the near side as a corollary of its leaf attachment, the one hub outside the domain", [], ""),
        ("The open result that the sign is one function of the side at every hub was entered at that hub and dissolves at the re-wording.", [], ""),
    ],
    "170.1": [
        ("Beta decay is read as the release at the plus-two arm at the literal closure reading of the decay graph", [], ""),
        ("at that reading the catalogued neutron web is a non-sink carrying its dock edge at the plus-two arm", [], ""),
        ("the hub transfer at the dock of the arm, the hub unit of the released body attached at the dock by its in-leg together with the kept-leg hub at the same dock, adds with its end the cross-core object of the released body by an explicit isomorphism at dimensions and links", [], ""),
        ("read through the hub at the order two of its (2,4) in-leg with the dressing of that link at both dressing readings, the readout 2916 at the hub spelling and the cell 1296 over 17689 at both spellings", [], ""),
        ("so the beta rate cell is an output of the rate law at a non-sink", [], ""),
        ("at the stop condition of the pairing move the neutron carried no forward edge, and the object the dock move added, a dock of one site at one link of measure 6, parted from the released body and from the cross-core object", [], ""),
        ("at the literal closure reading the forward edge is present and the hub transfer adds the cross-core object, while the dock move of the previous reading remains the first admitted edge at the dock at which the release parts, its object a dock, printed beside", [], ""),
        ("the deciding table is re-derived in Lean", [], "the two tied declarations are open-result records; no table is tied"),
    ],
    "229.1": [
        ("The four bodies joined to their measuring bodies at the attachment site, with the joined read opened at the dock or nowhere under the transport of reference, give the unit and the register as outputs of the readout", [], ""),
        ("every one of those readings is equal at two routes: the unit and the register at the three unit readings of the eleven joined webs and of the six body rows at the two ladders, fifty-one readings in all", [], ""),
        ("read equal at the two routes at 51 of 51, 33 of 33 at the web level and 18 of 18 at the body level", [], ""),
        ("At 2 of 4 bodies the outputs re-meet the established faces, the electron and the muon at the dock at the register 16 against the dock-struck readings 67 and 107 at 16", [], ""),
        ("and at 2 of 4 they part, the tau first, its joined read opening at neither the hub nor the dock", [], ""),
        ("the hub unit is reached by no joined read, so the hub-or-dock question of the previous releases resolves to the dock wherever the joined read opens at all", [], ""),
        ("The first parting of the previous open result, the tau at the library ladder, re-reads at the second route as the same output, neither at no register against the dock-struck reading 147 at the register 16, a parting of the inherited face and not of the construction, which is thereby closed", [], ""),
        ("at the derived ladder the tau carries the two (2,3) docks of the electron and re-meets the dock at 16, the ladder reading 3 of 4", [], ""),
        ("the deciding table is re-derived in Lean", [], "the one tied declaration is an open-result record; no table is tied"),
    ],
    "232.2": [
        ("The candidate read from the readings table of that result is the Hadamard reached at the invariant cut at order two", [], ""),
        ("the Hadamard tensor the identity is a member of the finite gate group at the pair register (the order 2304 unchanged)", [], ""),
        ("so no move beyond the twenty-eight is admitted, the gate set is closed under the candidate and is the limit of the move set, and the open result is a residue of the finite group", [], ""),
        ("the infinite order belongs to the invariant-cut rotations alone, beyond the register face at every seeded sink carrying a qubit-qubit link (46412 of 46563 sinks), candidates the saturation cannot admit under the termination theorem", [], ""),
    ],
    "233.1": [
        ("The distilling line of the text is a corollary of the Axiom and of the qudit line", [], ""),
        ("the output of an act exists only as a configuration of maximally entangled pairs", [], ""),
        ("the act on two qudits is a projection onto a pair on the common support (the rank identity at 69 of 69 operators)", [], ""),
        ("and the pair lives on the common support alone", [], ""),
        ("the four steps are derived with no posed step", [], ""),
    ],
    "233.2": [
        ("The first clause of the read line is a washed clause, a corollary at the level of the known", [], ""),
        ("the basis-matched restriction of the pair projections is the family of maximally entangled projections, characterised by the certainty the known clause names (a read at every basis if and only if the reduced state is flat, the Gram identity at 120 of 120)", [], ""),
        ("the basis of that clause is a label beyond the projection, washed at the matching group", [], ""),
        ("and the read cap, the dimension plus one, is the count of mutually unbiased bases", [], ""),
    ],
    "233.4": [
        ("The re-typing of the 1653 subdivisions of the census as paths of four one-act moves is forced by the text at the pair-ket unit of the first reading, every subdivision decomposing into one read and three joins", [], ""),
        ("under that unit the one-act clause is a theorem at every instance", [], ""),
        ("at the unit of the projection with its common-support byproduct the subdivision is one inverse act at 1448 instances and the path alone at 205, both readings printed and neither chosen", [], "1448 + 205 = 1653"),
    ],
    "I.16.2": [
        ("The seventy-five standard-model entries of the textbook corpus stand at the grades listed in the table that follows.", [], ""),
    ],
}

# The 16 rows of round 6: labels in their table-check declarations (not graded in round 6), matched to round 6's clauses
# (lean_resolve_data.ROWS, by clause index). Every one is RECORDED; none of these modules has a check of kind REL.
LABELS16 = {
    "2.6": {0: [("row_525542_TTCF", 'rows.filter (fun r => r.eq == "YES")).length = 32', R)]},
    "2.9": {1: [("row_532049_DIM1", 'r.cell == "C-EQD1" && r.value_E1 == "YES" && r.value_E2 == "YES" && r.letter == "EXCLUDED"', R)]},
    "9.3": {2: [("row_528882_LTGO", 'r.cell == "C-EQOT" && r.cell_value == "YES" && r.letter == "EXACT"', R)]},
    "18.2": {3: [("row_531081_GPOE", 'r.cell == "C-EQGL" && r.value_E1 == "YES" && r.value_E2 == "YES" && r.letter == "MEASURE"', R)]},
    "74.1": {0: [("row_527909_NVCI", KV("cell.C-KERN", "YES"), R)],
             2: [("row_527909_NVCI", KV("cell.C-CHRG", "YES"), R)]},
    "74.4": {0: [("row_526972_NRGE", KV("eqer.junctions.equal", "577 of 577"), R)]},
    "78.3": {2: [("row_527777_AJBP", KV("eqab.order2.identity", "YES"), R)],
             3: [("row_527777_AJBP", KV("eqab.order3.zero", "YES"), R)]},
}

# Data modules whose content is that of a box, by the box's own label or by content (judgment where stated)
OWN_BOX = {
    "D535360_NDFE": ("3.1", "label", "box label 'box:I-1-D535360-ndfe-t1:1'"),
    "D535412_OLUS": ("11.1", "label", "box label 'box:I-2-D535412-olus-t1:1'"),
    "D535680_RPHG": ("13.1", "label", "box label 'box:I-2-D535680-rphg-t1:1'"),
    "D536716_SSFP": ("186.1", "label", "box label 'box:VIII-1-D536716-ssfp-t1:1'"),
    "D536583_PSYO": ("I.7.1", "label", "section label 'mod:II-10-D536583-psyo-t1'"),
    "D525257_HTOR": ("3.1", "content", "the OPEN-READ-OF-RECORD rows enter 3.1's sentence on the reads Z, X, XZ, XZ2"),
    "D532499_HPK2": ("D.3", "content", "16 webs equal at grammar, dressing and direct sum; the code letters H, P, K occur in order in D.3's title"),
}
OWN_BOX_ANCHORS = {
    "D525257_HTOR": "the three open reads $X$, $XZ$ and $XZ^{2}$ fix none and free three",
    "D532499_HPK2": "the sixteen deciding webs standing at their grammar, dressing and direct sum",
}

# Strict reading (REVIEW4-CT, R4-CT-8 and R4-CT-9). The verdicts above (and round 6's in lean_resolve_data.py) count a
# clause when its counterpart states some part of its finite content. These credited clauses are stated only in part;
# under the strict reading (a clause counts only when all its finite content is stated) they do not count.
# (row, clause number from 1) -> the part that is stated, and what is not
PART_ONLY = {
    ("9.4", 3): "LemmaB.zeroTransfer over the 172 entered step triples: the third triple (0,0,0) is the conclusion "
                "'zero transfer' as entered; second = minus first is a relation over entered triples, read as "
                "'zero register displacement' by position only (the library's displacement of the object of record, "
                "LemmaB.massless_faces, is (-1,1,0) and is not used)",
    ("196.2", 2): "12 open sites at the kit tau are computed; the 17 sites are only in derived_bodies_counts, outside "
                  "`checks`",
    ("217.1", 3): "only rows.length = 228 is checked; 14 bodies, 114 links, 78 open and 150 washed are not",
    ("217.1", 4): "only the site bound is checked (W); the free-direction fields are not checked, and the directions "
                  "read (R_diff_site >= 0) are RECORDED",
    ("234.1", 2): "the sums are checked at the 246 network-and-polarity rows (82 networks x 3 polarities); no per-link "
                  "row is entered",
}
