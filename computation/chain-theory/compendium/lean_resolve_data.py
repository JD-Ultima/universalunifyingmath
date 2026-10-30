"""Curated data for lean_resolve.py (round 6: the restated AlphaDerivation theorems and the kernel modules).

Every quoted string here is checked verbatim by lean_resolve.py against the source: Lean anchors against the
declaration text in compendium/lean, box anchors against the box body in compendium/tex (whitespace-normalised).

Content kinds of a conjunct (what its proof establishes):
  COMPUTED      the kernel evaluates functions defined in the library on inputs entered in the library
                (matrices, letter lists, web and record lists, step lists); the proposition is about the values
  GENERAL       proved for all values of its variables by a tactic proof
  ENTERED       a property that an entered table has as entered (the check reads the table's own fields)
  DEFINITIONAL  holds by unfolding a definition that is written to yield it (a case split on the object)
  +LITERAL      the conjunct also contains equalities between entered constants and literals (true by rfl)
"""

# conjunct (AlphaDerivation.Statements.<name>) -> (kind, note, [anchors in the Proofs-level declaration text])
CONJ = {
    "OrbitReading_orbit_reading": (
        "COMPUTED+LITERAL",
        "typing of the 172 entered webs against the 172 entered records: 79 odd-orbit, 79 even-orbit, 11 fixed, "
        "3 re-entry webs; re-entry webs are commutators +(x,y) -(y,x) of two letters at slot 1 (the dimension-2 slot); "
        "the odd orbits, their swaps and the fixed webs cover the 169 slot-pair labels; literal part 2*79+11+3 = 172",
        ["(websOf .reentry).all isCommutatorWeb = true",
         "t1.2.1.slot == 1 && t1.2.2.slot == 1 && t2.2.1 == t1.2.2 && t2.2.2 == t1.2.1",
         "2 * Expected.orbitRecords + Expected.fixedRecords + Expected.reentryCount = Expected.recordCount"]),
    "OrbitReading_no_fixedness_asked": (
        "COMPUTED",
        "matrix action of the entered generators: web R001 is not fixed by the third qutrit generator; the web at "
        "index 168 (R169) is fixed by every qutrit and qubit generator",
        ["fixedUnder (LemmaA.gens3.getD 2 []) (identity 2) (webs.getD 0",
         "(webs.getD (Expected.recordCount - Expected.reentryCount - 1) ⟨\"\", .fixed, []⟩).id = Expected.recordFixedUnderGroupId"]),
    "IdAlg_closure": (
        "COMPUTED", "the 12 entered 5x5 Gaussian-integer matrices (su(2)+su(3)+u(1) block basis) close under the bracket",
        ["closureHolds = true ∧ basis.length = 12"]),
    "IdAlg_centre_one_dimensional": (
        "COMPUTED", "the centre of that 12-dimensional algebra (genre (2,3)) is one-dimensional: TU in the kernel of the "
        "adjoint constraint and an 11x11 minor nonzero",
        ["ZM.det (ZM.minor adjointConstraint centreMinorRows (List.range 11)) ≠ 0"]),
    "IdAlg_tu_central": (
        "COMPUTED", "TU = diag(3i,3i,-2i,-2i,-2i) commutes with every basis element",
        ["basis.all (fun b => isZero (comm n TU b)) = true"]),
    "IdAlg_forms": (
        "COMPUTED", "Killing form: TU in its kernel, rank 11; trace form nondegenerate, entry (11,11) = -30",
        ["ZM.det traceForm ≠ 0 ∧ ZM.entry traceForm 11 11 = -30"]),
    "IdAlg_centre_acts": (
        "COMPUTED", "rho(TU) = i Id_6 on C^2 (x) C^3; at (2,2) the entered element TU22 = diag(i,i,-i,-i) acts as 0 "
        "(TU22 is entered, not computed as the (2,2) centre)",
        ["rho TU = smul G.iu (identity 6)", "isZero (rho22 TU22) = true",
         "def TU22 : GM := mk 4 (fun i j => if i == j then (if i < 2 then ⟨0, 1⟩ else ⟨0, -1⟩) else G.zero)"]),
    "LemmaB_readable_and_acting": (
        "COMPUTED", "readable genres from alphabet {2,3,4} and read cap 6 (both definitions): [(2,2),(2,3)]; "
        "'abelian' (pairScalar = (d2-d1)/gcd != 0) selects [(2,3)]",
        ["def alphabet : List Nat := [2, 3, 4]", "def readCap : Nat := 6",
         "def abelian (g : Nat × Nat) : Bool := pairScalar g != 0"]),
    "LemmaB_centre_generators": (
        "COMPUTED+LITERAL", "centralFunctional (2,3) = (3,-2) -> [3,3,-2,-2,-2]; (2,2) -> (0,0) by the case g.1 == g.2; "
        "literal part: [1,1,-1,-1] = Expected.centreGenerator22 (an entered list equal to a literal list)",
        ["(List.replicate 2 (1 : Int) ++ List.replicate 2 (-1 : Int)) = Expected.centreGenerator22",
         "if g.1 == g.2 then (0, 0) else"]),
    "LemmaB_spin_content": (
        "COMPUTED+LITERAL", "angular-momentum labels from entered frame contents (2 -> [1], 3 -> [0,0,0], 4 -> [0,2]): "
        "cut register at (2,3) is 8 copies of label 2 (spin one), dim 24; at (2,2) labels [0,2,4], dim 9; "
        "literal parts 2*2 = 1+3 and 1*(1+1) = 2",
        ["| 2 => [1] | 3 => [0, 0, 0] | 4 => [0, 2] | _ => []", "2 * 2 = 1 + 3 ∧ 1 * (1 + 1) = 2"]),
    "LemmaB_massless_faces": (
        "COMPUTED", "no closed non-backtracking walk on the open one-link graph (2,3): rest energy 0, displacement (-1,1,0); "
        "the closed (2,3) form and the open (2,2) form are not 'massless'",
        ["closedReads objectOfRecord = [] ∧ restEnergy objectOfRecord = 0"]),
    "LemmaB_unique_object": (
        "COMPUTED", "of the 4 readable objects ((2,2),(2,3) x open, closed) only (2,3).open satisfies "
        "abelian && conserved && massless && spinOne",
        ["readableObjects.filter allFour = [objectOfRecord] ∧ readableObjects.length = 4",
         "def allFour (o : Obj) : Bool := abelian o.genre && conserved o.genre && massless o && spinOne o.genre"]),
    "LemmaB_objectOfRecord_electromagneticLike": (
        "COMPUTED", "the structure ElectromagneticLike bundles five Booleans defined in LemmaB (abelian, conserved, massless, "
        "spinOne, zeroTransfer) at (2,3).open; 'conserved' is 'generatorBlocksScalar && centralFunctional != (0,0)', and "
        "generatorBlocksScalar is true at every pair (own check), so 'conserved' and 'abelian' both reduce to d1 != d2",
        ["def conserved (g : Nat × Nat) : Bool := generatorBlocksScalar g && centralFunctional g != (0, 0)",
         "(List.replicate g.1 a ++ List.replicate g.2 b).length == g.1 + g.2"]),
    "Theorem5_correspondence_at_structure": (
        "DEFINITIONAL+COMPUTED", "five opaque hypotheses, none used (underscore binders; the proof term does not mention "
        "them); onShellCoupling o is defined as (registerSum, register) when o is the object of record and (0,1) otherwise, "
        "so its value at (2,3).open is the closed form by definition; registerSum = (79,48,42,3) is computed from the "
        "172 entered records",
        ["if o == LemmaB.objectOfRecord then (registerSum, register) else (Q4.zero, 1)",
         "theorem correspondence_at_structure (_cap : ReadCapWoronowiczStormer)",
         "opaque ReadCapWoronowiczStormer : Prop"]),
    "Census_weights_from_construction": (
        "COMPUTED+LITERAL", "exchange pairs (x,y) from the entered records' parity and kinds by the end rule (same: (-p,1), "
        "mixed: (1,-p), closure: from ClosureSign); squared weight = 3^[x=-1] * 2^[y=-1]; histogram (1,79),(2,48),"
        "(3,42),(6,3); sector counts also from the formulas in a = 8, b = 3. The record-wise identity "
        "((2-x)(3-y))/2 = 3^[x=-1] 2^[y=-1] holds for every (x,y) in {+-1}^2 (own check), so it does not depend on the "
        "records; literal part Expected.endMap = \"X:centre:3;Y:far:2\"",
        ["records.all (fun r => squaredWeightProduct r == squaredWeightCount r) = true",
         "def squaredWeightProduct (r : Rec) : Int := let (x, y) := exchanges r; ((2 - x) * (3 - y)) / 2",
         "| .same => (-(p r), 1)", "Expected.endMap = \"X:centre:3;Y:far:2\""]),
    "IdAlg_tensor_commutant_scalars": (
        "COMPUTED", "the commutant of rho(algebra) on C^6 is the scalars: identity in the kernel of 432 linear constraints "
        "on 36 unknowns, a 35x35 minor nonzero",
        ["ZM.det (ZM.minor tensorConstraint tensorMinorRows tensorMinorCols) ≠ 0"]),
    "Theorem1_register": (
        "COMPUTED", "products of 2 and 3: whole = 36, cut = 24, dimR = 24*36*36 = 31104, one-link register 6",
        ["dimR 2 3 = 24 * 36 * 36 ∧ dimR 2 3 = Expected.closedFormDen"]),
    "Theorem1_readout_is_trace_form": (
        "COMPUTED+LITERAL", "the conjunction IdAlg.forms and IdAlg.tensor_commutant_scalars, plus two entered constants "
        "equal to literals (Expected.tensorCommutantDim23 = 1, Expected.killingRank23 = 11, proved `by decide`); "
        "the statement does not relate the two constants to the computed ranks",
        ["⟨IdAlg.forms, IdAlg.tensor_commutant_scalars, by decide, by decide⟩",
         "Expected.tensorCommutantDim23 = 1 ∧ Expected.killingRank23 = 11"]),
    "LemmaB_zero_transfer": (
        "ENTERED", "the 172 entered step triples: third component (0,0,0) in every entry as entered, second component the "
        "negative of the first; plus steps.length = records.length",
        ["steps.all fun s => s.2.2.2 == (0, 0, 0) &&"]),
    "Theorem5_sub_add_sub": ("GENERAL", "(x - y) + (y - z) = x - z over Int", ["(x - y) + (y - z) = x - z"]),
    "Theorem5_telescoping": (
        "GENERAL", "vertex W c [(a,c'),(c',b),(c',c')] = vertex W c [(a,b)] for every functional W: a sum of "
        "differences telescopes", ["def subdivided (a c b : Nat) : List (Nat × Nat) := [(a, c), (c, b), (c, c)]"]),
    "Theorem5_one_link_iff": ("GENERAL", "count + (W a - W b) = count iff W a = W b",
                              ["vertex W count [(a, b)] = count ↔ W a = W b"]),
    "Theorem5_line_of_record": (
        "COMPUTED+LITERAL", "evaluations of the unit and wall functionals at the 9 alphabet pairs and 6 subdivisions "
        "(wall vertex at (3,2) = -6 - 6 = -12); literal parts Expected.wardRecordIdentity = \"172 of 172\", "
        "Expected.wardLinkPaths = \"48 of 48\" (entered strings equal to literals)",
        ["Expected.wardRecordIdentity = \"172 of 172\" ∧", "Expected.wardLinkPaths = \"48 of 48\""]),
    "Theorem5_ward_identity": (
        "GENERAL", "the conjunction of telescoping, one_link_iff, unit_washed and one evaluation",
        ["⟨telescoping, one_link_iff, unit_washed, line_of_record.1⟩"]),
    "Theorem5_unit_washed": (
        "GENERAL", "for the constant functional every link difference is 0, so the vertex equals the count",
        ["def unitFunctional : Functional := fun _ => 1"]),
    "Theorem5_readings_washed": (
        "COMPUTED", "41 readings enumerated (32 path, 8 single, 1 open); orderings 1/8/4/24; the unit functional holds at "
        "every reading and read class; the wall functional fails at one link (-12). That the unit functional's identity "
        "holds on every path follows from unit_washed (own check identityAt_unit)",
        ["readings.all (fun rd => readClasses.all fun ends => recordHolds unitFunctional rd ends) = true"]),
}

# the restated chapter declarations: name -> (box the content and the code point to, evidence)
OWN_BOX = {
    "ChainTheory.Foundations.Ch01.row_524038_ORKC": ("2.2", "C", "code = initials of 'Orbit Reading ... Knowledge Clause'; "
        "the box's 'invariant under every matching' is no_fixedness_asked"),
    "ChainTheory.StandardModel.Ch10.row_518160_AGUQ": ("9.2", "T", "code = 'Abelian Genre UniQueness'"),
    "ChainTheory.StandardModel.Ch10.row_520119_GWDC": ("18.3", "C", "code = initials of 'Gain ... Washed Direction Count'; "
        "18.3: 'The gain squared of a record is the product, over the signed exchanges a record survives, of the dimension "
        "of the end each exchange acts at'"),
    "ChainTheory.Foundations.Ch03.row_521104_OUTR": ("-", "-", "no box of section 56 states a tensor commutant, a register "
        "or a trace form"),
    "ChainTheory.Foundations.Ch05.row_518830_MDFN": ("74.2", "C", "74.2: 'the momentum transfer composes to the zero signed "
        "step along the direction register at every record'"),
    "ChainTheory.Foundations.Ch05.row_521443_WDCL": ("78.1", "C", "code = 'Ward iDentity ... Conservation ... Link'; 78.1 "
        "states one_link_iff and telescoping in words"),
    "ChainTheory.Foundations.Ch05.row_524064_WSWC": ("78.2", "C", "78.2: 'the identity holds at every one of the forty-one "
        "readings of the records the definitions admit'"),
    "ChainTheory.StandardModel.Ch12.row_527992_ESRL": ("196.2", "T", "196.2: 'Every site of the tau is typed under the "
        "definition of open'"),
}

# box anchors (verbatim) for the own-box attribution above
OWN_BOX_ANCHORS = {
    "2.2": "invariant under every matching",
    "18.3": "The gain squared of a record is the product, over the signed exchanges a record survives, of the dimension of the end each exchange acts at",
    "74.2": "the momentum transfer composes to the zero signed step along the direction register at every record",
    "78.1": "the vertex equals the count exactly when the per-dimension sums are constant",
    "78.2": "the identity holds at every one of the forty-one readings of the records the definitions admit",
    "196.2": "Every site of the tau is typed under the definition of open",
}

# per row: list of (verbatim clause of the box, [counterpart conjuncts or kernel checks], note)
K = "ChainTheory.Kernel.D527992_ESRL."
ROWS = {
    "2.1": [
        ("is a corollary of the Axiom at both clauses", [], "not a finite statement; nothing in Lean"),
        ("fix jointly with certainty whichever way the basis directions are matched across its pairs",
         ["OrbitReading_no_fixedness_asked"], "one entered web (R169) fixed by every generator, R001 not"),
        ("two reads commuting only up to a sign are a further record with their order as the sign, by the holonomy sign at dimension two",
         ["OrbitReading_orbit_reading"], "3 re-entry webs are commutators at the dimension-2 slot"),
    ],
    "2.3": [
        ("The one-act clause is a theorem of the derived act", [], ""),
        ("every entry of the census of moves then carries at most one act", [], ""),
        ("a subdivision being a path of four one-pair acts", [], ""),
    ],
    "2.6": [
        ("all thirty-two derived cells with a phase input read the same at both of its faces", [], ""),
        ("the re-spelling of the tenth from the ninth is an instrument", [], ""),
    ],
    "2.9": [
        ("A pair with an end of dimension one is a product state", [], ""),
        ("excludes the one-dimensional end", [], ""),
    ],
    "2.10": [
        ("No second copy of a pair exists", [], ""),
        ("a read acts once on the one copy of each pair the qudit meets", [], ""),
        ("has no statistics over fresh copies of the pair to draw on", [], ""),
    ],
    "9.2": [
        ("Exactly one readable type carries an acting centre", ["LemmaB_readable_and_acting"],
         "readable [(2,2),(2,3)] under cap 6; 'acting' as abelian = pairScalar != 0"),
        ("the centre of the identification algebra has dimension one at (2,2)", [],
         "not computed: Expected.centreDims = [1, 1] is defined and used by no theorem"),
        ("and at (2,3)", ["IdAlg_centre_one_dimensional", "IdAlg_tu_central"], "computed"),
        ("it acts on the pair space as a scalar only at (2,3)", ["IdAlg_centre_acts"], "i Id_6 at (2,3); 0 at (2,2) (entered TU22)"),
        ("whose open read alone is abelian, conserved, massless and of spin one",
         ["LemmaB_unique_object", "LemmaB_objectOfRecord_electromagneticLike"], "by the Boolean definitions of the four words"),
    ],
    "9.3": [
        ("is a function of the link type alone", [], "onShellCoupling is defined on objects; no proposition states it"),
        ("re-meeting the derived fine-structure cell at (2,3)", ["Theorem5_correspondence_at_structure"],
         "by the case-split definition of onShellCoupling"),
        ("printing the other types as exact faces", [], "onShellCoupling is (0,1) at every other object"),
    ],
    "9.4": [
        ("The open (2,3) link carries no mass face", ["LemmaB_massless_faces"], "computed on the one-link graph"),
        ("its two-read record at the qubit end carries spin one", ["LemmaB_spin_content"], "from the entered frame contents"),
        ("zero register displacement equals zero momentum transfer at the one-link star",
         ["LemmaB_objectOfRecord_electromagneticLike"], "its zeroTransfer field: a property of the entered step table"),
    ],
    "9.5": [
        ("The passed-hub count is monotone against the hub dressing", [], ""),
        ("at every pair of bodies and every face", [], ""),
    ],
    "9.6": [
        ("The one-link tables are defined at all six types by the machinery of the support-rank result", [], ""),
        ("the (2,3) table agreeing with that result at every cell", [], "no one-link table is restated (OneLinkWord_* restated nowhere)"),
    ],
    "18.2": [
        ("The multiplicity of every record is the product of the dimensions of the ends it passes, one at no end, two at the qubit end, three at the qutrit end and six at both",
         ["Census_weights_from_construction"], "the squared weight is defined as that product; the checked identity holds for every exchange pair"),
        ("at all one hundred seventy-two records of the qutrit-qubit link", ["Census_weights_from_construction"],
         "records.all over the 172 entered records; histogram sums to 172"),
        ("its four end-set sums agreeing with the support-rank result", [],
         "stated in GainLemma.end_set_sums, which no chapter declaration restates"),
        ("a link-product measure equal to the Born count at every end set and every link type", [],
         "GainLemma.born_count_is_a_measure (not restated); only (2,3) anywhere"),
    ],
    "56.1": [
        ("the forty-eight matchings of the five coordinates that the construction admits", [], ""),
        ("0 of 322714 records at 0 of 3072 pairs of a web and a matching over the 64 catalogued webs", [], "322714 occurs nowhere in the Lean tree"),
        ("carries 319332 of those records", [], "319332 occurs nowhere in the Lean tree"),
        ("384 maps against 48", [], ""),
        ("the deciding table re-derived in Lean at the forty-eight matchings", [], "no such declaration found"),
    ],
    "74.1": [
        ("The left kernel of the adopted move matrix over the grade vector has rank zero at both faces", [], ""),
        ("the lepton number as the dark-triangle count is broken at a freed-at-cut column", [], ""),
        ("the dart charge at the pure open-leg column", [], ""),
    ],
    "74.4": [
        ("Baryon number read as the unsigned junction count is conserved at every interaction and moved by the non-interactions alone", [], ""),
        ("the structure counts at the moves are broken at the pinned release", [], ""),
        ("the triangle count is broken with the join count not additive", [], ""),
        ("the lepton number as the dark-triangle count under the hub step is unequal at the freed-at-cut attachments alone", [], ""),
    ],
    "78.3": [
        ("every record of the direct sum over records is classified into the body, the cross records across the join and the apparatus self-terms", [], ""),
        ("the dock-struck first power of the electron, sixty-seven, with the sign plus", [], ""),
        ("at order two an exact element of the field with the sign plus and a nonzero cross term", [], ""),
        ("the order-three term zero", [], ""),
    ],
    "196.2": [
        ("Every site of the tau is typed under the definition of open, the far end of the entered link open exactly when the near end is the larger",
         [K + "check_1"], "every site of the 5 built bodies typed by openAt (69 rows, computed)"),
        ("twelve of the seventeen sites are open attach sites", [K + "check_6"],
         "12 open sites at the kit tau (computed); 17 sites only in derived_bodies_counts, outside `checks`"),
        ("the joined face of the tau is defined at one of them", [], ""),
        ("the body alone at the first order and the second-order object an exact element of the field", [], ""),
    ],
}
