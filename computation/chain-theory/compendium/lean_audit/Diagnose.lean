import AxiomDeclaration.Structure
import AxiomDeclarationTwo.Structure

/-!
Why `AxiomDeclaration.Theorems` and `AxiomDeclarationTwo.Theorems` fail under Lean 4.35.0 (ours; compiled against
Build B's `out`, which holds both `Structure` modules unchanged). The five statements are those of
AxiomDeclaration/Theorems.lean:42, :52 and AxiomDeclarationTwo/Theorems.lean:42, :52, :56, verbatim. Under 4.35.0
instance resolution does not unify `Fin w.qudits` with `Fin n` for a plain `def` world `w`, so `decide` finds no
`Decidable` instance; `rfl` (componentwise for the conjunctions) proves all five. The instance exists for the type
as stated:
-/
open AxiomDeclaration

#synth DecidableEq (List (Fin worldThree.qudits))

theorem a42 : worldThree.sharedQudits = [(1 : Fin 3)] := by
  fail_if_success decide
  rfl

theorem a52 : worldThree.passedTo ⟨(1 : Fin 3), (0 : Fin 3), 0⟩ = [(0 : Fin 3), (2 : Fin 3)] := by
  fail_if_success decide
  rfl

namespace AxiomDeclarationTwo
open AxiomDeclaration (Pair World pairAdmitted worldThree worldOne)

theorem b42 : worldTwoThree.sharedEnd ⟨(0 : Fin 2), (1 : Fin 2)⟩ = some (1 : Fin 2) := by
  fail_if_success decide
  rfl

theorem b52 : worldTwoHub.genres = [(2, 4)] ∧ worldTwoHub.sharedEnd ⟨(0 : Fin 2), (1 : Fin 2)⟩ = some (0 : Fin 2) ∧ worldTwoHub.support ⟨(0 : Fin 2), (1 : Fin 2)⟩ = 2 := by
  fail_if_success decide
  exact ⟨rfl, rfl, rfl⟩

theorem b56 : worldTwoSix.shared (4 : Fin 6) = true ∧ worldTwoSix.sharedByDimension (4 : Fin 6) = true ∧ worldTwoSix.sharedEnd ⟨(0 : Fin 6), (2 : Fin 6)⟩ = some (2 : Fin 6) := by
  fail_if_success decide
  exact ⟨rfl, rfl, rfl⟩

end AxiomDeclarationTwo

#print axioms a42
#print axioms a52
#print axioms AxiomDeclarationTwo.b42
#print axioms AxiomDeclarationTwo.b52
#print axioms AxiomDeclarationTwo.b56
#print axioms AxiomDeclaration.World.sharedQudits
