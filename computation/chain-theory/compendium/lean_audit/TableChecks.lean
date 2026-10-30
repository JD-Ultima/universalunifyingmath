import ChainTheory.Data.D526019_WAMR
import ChainTheory.Data.D525335_SLBC
import ChainTheory.Data.D525516_NCRE
import ChainTheory.Data.D525257_HTOR
import ChainTheory.Data.D526261_RLRO

/-!
Round 7 checks (ours, not part of the source). Compiled with Lean 4.35.0 against Build B of the compendium
(`lean_build.py`, LEAN_PATH = Build B's `out`). Each theorem is a statement about the source's entered `rows` lists
that the source's own `checks` do not make: printed numbers that follow from the entered rows, and properties of the
entered values. `lean_tables.py` computes the same quantities from the text of the rows.
-/

set_option maxRecDepth 100000

namespace Round7
open ChainTheory.Data

/-- 28.1 and 217.1: the 14 entered rows of D526019_WAMR (tied to section 28) sum to 228 acts, 114 links,
78 open and 150 washed ends. Its `checks` state none of these sums. -/
theorem wamr_sums :
    (D526019_WAMR.rows.map (·.acts)).foldl (· + ·) 0 = 228 ∧
    (D526019_WAMR.rows.map (·.links)).foldl (· + ·) 0 = 114 ∧
    (D526019_WAMR.rows.map (·.c_open)).foldl (· + ·) 0 = 78 ∧
    (D526019_WAMR.rows.map (·.washed)).foldl (· + ·) 0 = 150 := by decide

/-- 217.1: in the 228 entered reads of D525335_SLBC, 78 have an open end and 150 are washed, and each link is entered
at the two orientations AB and BA (114 each). Not stated by its `checks`. -/
theorem slbc_counts :
    (D525335_SLBC.rows.filter (fun r => r.c_open == "OPEN")).length = 78 ∧
    (D525335_SLBC.rows.filter (fun r => r.c_open == "WASHED")).length = 150 ∧
    (D525335_SLBC.rows.filter (fun r => r.orientation == "AB")).length = 114 ∧
    (D525335_SLBC.rows.filter (fun r => r.orientation == "BA")).length = 114 := by decide

/-- 217.1: the W counts that the source's `checks` relate (W_target_site - W_source_site = W_diff_site >= 0) do not
change at any entered read; the free-direction W difference is 0 as well. -/
theorem slbc_W_constant :
    D525335_SLBC.rows.all (fun r => r.W_diff_site == 0 && r.W_diff_free == 0 &&
      r.W_target_site == r.W_source_site) = true := by decide

/-- 217.1: the R counts rise at 54 entered reads; R_target_site - W_source_site = R_diff_site at every row, a relation
the source's `checks` do not state (they check only 0 <= R_diff_site). -/
theorem slbc_R :
    D525335_SLBC.rows.all (fun r => r.R_target_site - r.W_source_site == r.R_diff_site) = true ∧
    (D525335_SLBC.rows.filter (fun r => decide (0 < r.R_diff_site))).length = 54 := by decide

/-- 234.1: the order-three count is 0 in every entered row of D525516_NCRE, and 18 SLOT rows lie outside the library
group KIT (13 BODY, 3 NEUTRAL, 2 STAR). Neither is stated by its `checks`. -/
theorem ncre_rows :
    D525516_NCRE.rows.all (fun r => r.order3 == 0) = true ∧
    (D525516_NCRE.rows.filter (fun r => r.polarity == "SLOT" && r.group != "KIT")).length = 18 ∧
    (D525516_NCRE.rows.filter (fun r => r.polarity == "SLOT" && r.group == "BODY")).length = 13 ∧
    (D525516_NCRE.rows.filter (fun r => r.polarity == "SLOT" && r.group == "NEUTRAL")).length = 3 ∧
    (D525516_NCRE.rows.filter (fun r => r.polarity == "SLOT" && r.group == "STAR")).length = 2 := by decide

/-- 3.1 content in a table tied to 57.1: in D525257_HTOR the open read Z fixes one far direction and leaves two, and
the three open reads X, XZ, XZ2 fix none and leave three. -/
theorem htor_open_reads :
    (D525257_HTOR.rows.filter (fun r => r.reading == "OPEN-READ-OF-RECORD" && r.cap_read == "Z" &&
      r.fixed_set == "Z" && r.fixed_count == 1 && r.unfixed_count == 2)).length = 1 ∧
    (D525257_HTOR.rows.filter (fun r => r.reading == "OPEN-READ-OF-RECORD" && r.fixed_set == "NONE" &&
      r.fixed_count == 0 && r.unfixed_count == 3)).length = 3 := by decide

/-- D.3: the 16 entered webs of D526261_RLRO are 13 BODY and 3 BESIDE (the source checks both counts). -/
theorem rlro_webs :
    D526261_RLRO.rows.length = 16 ∧
    (D526261_RLRO.rows.filter (fun r => r.standing == "BODY" && r.equal == "YES")).length = 13 := by decide

end Round7

#print axioms Round7.wamr_sums
#print axioms Round7.slbc_counts
#print axioms Round7.slbc_W_constant
#print axioms Round7.slbc_R
#print axioms Round7.ncre_rows
#print axioms Round7.htor_open_reads
#print axioms Round7.rlro_webs
