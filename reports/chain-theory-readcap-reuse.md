# Reusing the operator-theoretic read cap

November 3, 2026. We freshly compiled the reusable positive-map theorem in
[SproutSeeds/math-lab](https://github.com/SproutSeeds/math-lab/tree/48d14da643ec5f6e40ea0e2ed9f2a5eb9c8d5fda)
and connected it to the author's unchanged `ChainAlpha.ReadCap` binder. This
reuses classical mathematics and existing formal proofs; it supplies no new
selection argument or physical confirmation of Chain Theory.

`ChainWoronowicz.woronowicz` proves that every positive complex-linear map
from the two-by-two to the three-by-three complex matrix algebra is decomposable.
`ChainWoronowicz.readcap` then proves, for both end dimensions at least two,
that every positive map in either direction is decomposable exactly when the
dimension product is at most six. Positivity quantifies over all positive
semidefinite inputs, including the reality of their output quadratic forms.
Decomposability means an exact finite Kraus sum of completely positive maps and
completely positive maps after transpose. This is a universal operator theorem,
not a finite matching-group test or a numerical SDP experiment.

The [adapter source](../proofs/lean/chain-theory/readcap-reuse/ReadCapReuse.lean)
contains five named theorems. It passes the existing cap directly into the
author's abstract binder, proves its dimension-level equivalence with
`readableByReads`, and applies the original `selection_under_read_cap` theorem
under its explicit readability and unequal-dimension assumptions. Two controls
retain an important distinction: the weaker lowercase `ChainAlpha.readable`
admits `(3,3)`, whereas operator-theoretic readability does not. It must not be
silently substituted for the operator predicate. The equivalence with
`readableByReads` is restricted to dimensions at least two.

The exact source pins are math-lab
`48d14da643ec5f6e40ea0e2ed9f2a5eb9c8d5fda`, Mathlib
`1b3c41969602424bd2f04b3d4f731734e5b00cde`, and the author's
`a28db8998bce5bd5ef11195c05f13c86cd9dac5c`. Lean is the official Linux
4.35.0-rc4 release, compiler commit
`c29b6dda4f7c20e3eeaa717c4e565663c5cfa364`. Its published archive SHA256,
the eight transitive dependency revisions, all 54 local module hashes, the
unchanged author-file hash and the public adapter's exact copy hash are recorded
in [verification.json](../proofs/lean/chain-theory/readcap-reuse/verification.json).
The local closure comprises nine MathLab modules and 45 vendored PSI modules;
the latter's upstream provenance is recorded separately in that receipt.

The initial run freshly compiled those **54 local modules in 65 seconds**.
Mathlib used its official master cache: all **9,010** required artifacts were
present, and Lake's no-build check reported **9,024 jobs up to date**. We did
not rebuild all Mathlib or run math-lab's entire 512-theorem suite. The local
build completed with 94 inherited style warnings and no proof error. The
unchanged author file and initial adapter probe also compiled successfully.
The exact public adapter source was subsequently copied, hash-checked and
freshly compiled once in **3 seconds**, using those verified dependency outputs.
No upstream proof source is copied into this repository.

All four requested Woronowicz declarations and all five public adapters report
only `propext`, `Classical.choice` and `Quot.sound`. There is no custom
Woronowicz/PPT axiom, `sorryAx`, or undisclosed cap premise. The
[selected actual transcript](../proofs/lean/chain-theory/readcap-reuse/ReadCapReuse.out)
contains the complete public-adapter compiler output; full setup and upstream
build logs remain private. Its header distinguishes those records. Reproduction
requires the pinned external MathLab project, its dependencies and an unchanged
author `ChainAlpha.lean` copy; this adapter is not a standalone core-Lean file.
The receipt records the commands and zero exit statuses.

For reproduction, prepare an LF checkout of the pinned MathLab project and its manifest-pinned dependencies with the stated Lean release. Place the unchanged author `ChainAlpha.lean` and this adapter at its root. After acquiring the official Mathlib cache, run the receipt's commands in that project:

```bash
lake --no-cache --no-build build +Mathlib:olean
lake --no-cache --fail-fast build +MathLab.ChainTheory.Woronowicz:olean
lake --no-cache env lean -o .lake/build/lib/lean/ChainAlpha.olean ChainAlpha.lean
lake --no-cache env lean -o .lake/build/lib/lean/ReadCapReuse.olean ReadCapReuse.lean
```

The repository preserves the checked adapter/transcript as LF; external copies must retain the recorded bytes. These commands check the selected target, not every MathLab module.

The trust boundary includes the official compiler, Lean's kernel, imported
Mathlib cache artifacts and the stated axioms. This verifies the mathematical
cap and its formal interface. It does not derive the model's single-copy
certainty or allowed read operations from that cap, justify the record/weight
clauses, or identify the computed candidate with the photon or a measured
coupling. Those semantic and physical questions remain separate.
