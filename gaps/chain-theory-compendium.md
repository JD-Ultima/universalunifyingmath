# Gaps: Chain Theory compendium, release 0.1

Source commit `a28db899`. Each gap names what is missing, where, its effect, and what would close it.
Claim IDs refer to `claims/chain-theory-compendium.md`.

## G1. Identification maps from model quantities to physical quantities are not defined

| Missing map | Affected claims | Effect | Would close it |
|---|---|---|---|
| record-sum order n -> power of alpha/pi | CT-C015 (a_e), CT-C021 | no a_e coefficient can be read off; under the literal map T3 = 0 contradicts C3 (9e4 sigma) | a definition of the expansion parameter and the printed order-1 coefficient in it |
| "the electron unit" -> m_e c^2 (or anything else) | CT-C016 (hydrogen), CT-C026 (masses) | with m_e c^2, E1 is ~1e-9.2 of the Rydberg energy; 176.1 concedes m_p is "outside the comparison by orders" | a stated unit and its derivation |
| count metric / rung -> energy scale | CT-C020 (running) | the printed geometric running (1/36 per rung) cannot be compared with log running | the scale map |
| tick -> seconds | CT-C022 (decays) | rates are counts per tick; the source states width ratios part "by many orders" | a tick unit by derivation |
| office ratio -> mixing angle | CT-C018 | under the source's stated correspondence cos^2 theta_W = 1/2, 1067 sigma off | either a different correspondence or a different number |

## G2. Quantities used but not printed or not re-derived

- 141.1: the "closure" part of T2 ("an exact element of the field of the two roots"), the link list of
  the two-dock body (13 multisets of eight measures give the printed 67), and the eight self-terms
  (only one appears, in the Lean data table).
- 143.1: the two body registers 2^30 3^18 and 2^31 3^17 ("re-derived on the page not"); the "cut gain"
  that sets the fine structure; the Lamb interval's value.
- 65.1: the end-kind assignment of the three factors is stated, not derived.

## G3. Grades inconsistent with the stated chain rule

The front matter states that a composed claim carries the lowest status in its chain. 66 of 100
THEOREMs list an OPEN premise directly (robust to designator resolution); 84 list an OPEN or CONDITIONAL
premise. Examples among the physics theorems: 64.1, 64.2, 65.1 (gauge), 74.1 (conservation), 83.1
(three generations), 145.1 (Casimir, resting on CONDITIONAL 143.1), 217.1, 220.3, 221.1, 225.1. Also,
premise lists label CONDITIONAL statements "Theorem" (e.g. "Theorem 124.1", "Theorem 141.1").
Would close it: per-statement logical premises (not whole chapters), and re-grading under the rule.

## G4. Circular premise chains

Under the index's own designator mapping, 124.1 (alpha) and 16.1, 58.1 list each other (via the
chapter-form premises of 124.1); one strongly connected component of 22 statements. Would close it:
statement-level premises for 124.1 and 200.1-200.4.

## G5. Cross-reference ambiguity

82 of 162 chapter designators (e.g. "5-2") have two readings under two numberings (`designators.tsv`);
Part VII prose uses one ("Chapters 5-2 and 5-3 are the dipole chapters") while the index maps 5-2 to
the gauge-group chapter. 3 premise references (238.1, 239.1, 240.1) point to appendices with no boxed
statement. Only 43.4% of premise edges are echoed in the premise's `used_by`.

## G6. Lean coverage does not reach the claims

For 141.1, 143.1, 64.1, 64.2, 65.1 the named Lean declarations are opaque Props, open-result records,
or `decide` proofs about hand-entered `String` tables. No arithmetic of these chapters (T1 = 309,
E1 = 67/(2^28 3^15), E_n = E1/n^2, the Z6 kernel) is in Lean. Across the 64 Lean-tied theorems, 11
cite only opaque Props/records; 137 of 150 data modules hold only `String` fields. The build record and
axiom report were not re-verified (not built). Would close it: Lean statements whose types mention the
mathematical objects (matrices, groups, rationals) and whose proofs compute them.

## G7. Selection inputs that are not derived from the axiom

- Dimension alphabet {2,3,4}: inherited from "the sixty-four listed webs" (9.1 derivation), i.e. from an
  input network list (H.2, "The webs of the network list as seeds", is OPEN).
- Read cap d1 d2 <= 6 (9.1, CONDITIONAL): uniqueness of the (3,2) genre fails at cap 8 or 9.
- Adopted readings: W1 in 141.1 (W2, W3, PWONCE part); whole vs cut faces in 143.1 (factor 2, both
  printed, none selected); "the adopted reading" in 83.1 (readings beside give count 0, differing counts, or no bound).

## G8. Abstract vs boxes

- "parity, time reversal, and CPT" as conditional theorems: 61.1 is graded THEOREM; 60.1 is OPEN.
- "as a theorem, conservation laws": 74.1 proves the conserved-functional kernel is trivial; the
  conservation of energy-momentum and charge-like functionals is CONDITIONAL (74.2, 74.3).
- "the Chern number" as a theorem: the instrument chapter's table is titled "The Chern number, awaited";
  the flux quantization premise is OPEN (I.51.4).
- "the electron anomalous moment through third order", "the hydrogen spectrum", "the running of the
  couplings": none of the three chapters outputs a number comparable with the measured quantity
  (CT-C015, CT-C016, CT-C020).
- README "Tables 56" vs 54 in the index.

## G9. Physical content absent from the gauge-group theorem

The theorems give the local-unitary group of one link. Gauge fields, spacetime locality of the
symmetry, matter representations (66.x OPEN) and hypercharges are not derived. The Z6 global form is
consistent with all observed fields but is not experimentally established (CT-C013).
