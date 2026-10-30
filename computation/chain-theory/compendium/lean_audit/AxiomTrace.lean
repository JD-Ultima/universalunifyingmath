import Lean
import ChainTheory

/-!
Where `Classical.choice` enters under Lean 4.35.0 (ours): a shortest chain of constants, through types and values,
from five declarations of the source to `Classical.choice`. At 4.32.1 the source reports no `Classical.choice` at all.
-/
open Lean Elab Command

/-- a shortest chain of constants from `root` to `target` through types and values -/
def chainTo (env : Environment) (target root : Name) : List Name := Id.run do
  let mut parent : NameMap Name := {}
  let mut seen : NameSet := ({} : NameSet).insert root
  let mut queue := #[root]
  let mut i := 0
  while i < queue.size do
    let n := queue[i]!
    i := i + 1
    if n == target then
      let mut path := [n]
      let mut cur := n
      while cur != root do
        cur := (parent.get? cur).getD root
        path := cur :: path
      return path
    match env.find? n with
    | none => pure ()
    | some ci =>
      let cs := ci.type.getUsedConstants ++ ((ci.value? (allowOpaque := true)).map Expr.getUsedConstants |>.getD #[])
      for c in cs do
        if !seen.contains c then
          seen := seen.insert c
          parent := parent.insert c n
          queue := queue.push c
  return []

#eval show CommandElabM Unit from do
  let env ← getEnv
  for root in [``AlphaDerivation.LemmaB.abelian, ``AlphaDerivation.LemmaB.zeroTransfer, ``AlphaDerivation.Theorem5.identityAt,
               ``AlphaDerivation.Statements.Census_weights_from_construction, ``ChainTheory.Kernel.D527992_ESRL.checks] do
    logInfo m!"{root}: {chainTo env ``Classical.choice root}"
