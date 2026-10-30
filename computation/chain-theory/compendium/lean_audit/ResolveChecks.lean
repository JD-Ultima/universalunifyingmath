import AlphaDerivation

/-!
Round 6 checks (ours, not part of the source). Compiled with Lean 4.35.0 against Build B of the compendium
(`lean_build.py`, LEAN_PATH = Build B's `out`). Each theorem is stated about the source's own definitions.
-/

namespace Round6
open AlphaDerivation

/-- `generatorBlocksScalar`, one of the two clauses of `LemmaB.conserved`, is true at every pair:
the length of `replicate g.1 a ++ replicate g.2 b` is always `g.1 + g.2`. -/
theorem generatorBlocksScalar_always (g : Nat × Nat) : LemmaB.generatorBlocksScalar g = true := by
  unfold LemmaB.generatorBlocksScalar
  split
  simp

/-- On every pair of dimensions 1..12, `abelian`, `conserved` and `d1 != d2` are the same Boolean. -/
theorem abelian_conserved_distinct :
    ((List.range 12).map (· + 1)).all (fun a => ((List.range 12).map (· + 1)).all (fun b =>
      LemmaB.abelian (a, b) == LemmaB.conserved (a, b) && LemmaB.abelian (a, b) == (a != b))) = true := by
  decide +kernel

/-- The record-wise identity checked by `Census.weights_from_construction` holds for every exchange pair in
{+1,-1}^2, whatever the record: ((2 - x)(3 - y))/2 = 3^[x = -1] * 2^[y = -1]. -/
theorem weight_identity (x y : Int) (hx : x = 1 ∨ x = -1) (hy : y = 1 ∨ y = -1) :
    ((2 - x) * (3 - y)) / 2 = (if x == -1 then Census.dimX else 1) * (if y == -1 then Census.dimY else 1) := by
  rcases hx with rfl | rfl <;> rcases hy with rfl | rfl <;> decide

/-- `onShellCoupling` is the closed form at the object of record by its definition, and (0, 1) at every other
object. -/
theorem onShellCoupling_elsewhere (o : LemmaB.Obj) (h : o ≠ LemmaB.objectOfRecord) :
    Theorem4.onShellCoupling o = (Q4.zero, 1) := by
  unfold Theorem4.onShellCoupling Theorem4.twoPointRatio
  simp [h]

/-- The conclusion of `Theorem5.correspondence_at_structure`, without any of its five opaque hypotheses. -/
theorem correspondence_without_hypotheses :
    LemmaB.ElectromagneticLike LemmaB.objectOfRecord ∧
    Theorem4.onShellCoupling LemmaB.objectOfRecord = Theorem4.closedFormOfRecord :=
  ⟨LemmaB.objectOfRecord_electromagneticLike, Theorem4.value_at_structure LemmaB.objectOfRecord_electromagneticLike⟩

/-- The unit functional satisfies the identity on every path whatsoever (from `Theorem5.unit_washed`), so the
41-reading enumeration in `readings_washed` cannot fail at the identity step for it. -/
theorem identityAt_unit (path : List (List Theorem5.Link)) :
    Theorem5.identityAt Theorem5.unitFunctional path = true := by
  simp [Theorem5.identityAt, Theorem5.unit_washed]

/-- Literal parts of the restated conjuncts: each holds by `rfl` (an entered constant equals a literal). -/
theorem literal_parts :
    Expected.wardRecordIdentity = "172 of 172" ∧ Expected.wardLinkPaths = "48 of 48" ∧
    Expected.endMap = "X:centre:3;Y:far:2" ∧ Expected.tensorCommutantDim23 = 1 ∧ Expected.killingRank23 = 11 ∧
    (List.replicate 2 (1 : Int) ++ List.replicate 2 (-1 : Int)) = Expected.centreGenerator22 :=
  ⟨rfl, rfl, rfl, rfl, rfl, rfl⟩

/-- The step table behind `LemmaB.zero_transfer` has the zero transfer (0, 0, 0) entered in every one of its 172 rows. -/
theorem steps_entered_zero : steps.length = 172 ∧ steps.all (fun s => s.2.2.2 == (0, 0, 0)) = true := by
  decide +kernel

end Round6

#print axioms Round6.generatorBlocksScalar_always
#print axioms Round6.abelian_conserved_distinct
#print axioms Round6.weight_identity
#print axioms Round6.onShellCoupling_elsewhere
#print axioms Round6.correspondence_without_hypotheses
#print axioms Round6.identityAt_unit
#print axioms Round6.literal_parts
#print axioms Round6.steps_entered_zero
