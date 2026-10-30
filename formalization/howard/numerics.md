# Howard numerical statements (TH-005, TH-006, TH-008 – TH-011)

All values are printed by `computation/howard/numerics.py`; the output is in `numerics.out`. Exact statements use sympy. Floating-point statements give their precision. Every quote comes from the JRE #2152 transcript (TH-S2), at the timestamps shown. That transcript is third-party and machine-made (see `sources/howard.md`).

## N1. The $\sqrt2$ "loop" (TH-005) — JRE #2152 124:37–125:38

(a) Howard: "the square root of two cubed is the same value as the square root of two times two … divided by two, it gets back to 1.414 … That's a loop, that's x cubed being equal to two x, which is equal to x plus x. That's unnatural equation." He gives the digits $\sqrt2=1.414213562373095$ and $(\sqrt2)^3=2.82842712174619$.

(c) **Proposition N1.**
1. $(\sqrt2)^3=2\sqrt2$ and $(\sqrt2)^3/2=\sqrt2$, exactly. Every arithmetic step he states is true in standard arithmetic.
2. For $x\neq0$, $x^3=2x\iff x^2=2$. So the "loop" is the defining equation of $\sqrt2$ multiplied by $x$. It is not special to 2: $\sqrt c$ is a fixed point of $x\mapsto x^3/c$ for every $c>0$.
3. The fixed point is repelling, since $f'(\sqrt2)=3$.
   - **float64:** iterating $x\mapsto x^3/2$ from $\mathrm{fl}(\sqrt2)$ (error $9.67\times10^{-17}$) leaves $|x-\sqrt2|\le10^{-3}$ at step 28. The first-order prediction is $\log(10^{-3}/\varepsilon_0)/\log3=27.28$.
   - **mpmath, dps 30 / 50:** the iteration escapes at steps 59 / 101 (predicted 59.01 / 100.4).
   - **mpmath, dps 10–120:** the iteration escapes at 64 precisions. At the other 47, including dps 100, the rounded map fixes the rounded $\sqrt2$ exactly.

   So the calculator "loop" is real for the three iterations he performs. It is a rounding artefact for long runs, and with typed input 1.414 it drifts at once (1.413573, 1.412293, …).
4. The quoted $(\sqrt2)^3$ digits differ from the true value $2.8284271247461902909$ in the 9th decimal (…712**1**74… vs …712**4**74…). This is either a transcription error or a misreading. We have not checked it against the audio.

**Verdict.** The arithmetic holds exactly. The claim that the equation is "unnatural" is evaluative, not mathematical. Under the relabelling of arithmetic Thm A4, the loop is invariant (`arithmetic.md` §4), so it is not caused by $1\cdot1=1$.
(e) "every quantum mechanics … equation has that square root of two in there" is not well-posed: no class of equations is specified. → TH-G05.

## N2. "The square root of a number added to itself would equal more than that number squared … that's what happens with most of the numbers" (TH-006) — 39:50

(b) There are two readings:
- R1: $2\sqrt n>n$, which holds iff $n\in(0,4)$.
- R2: $x+x>x^2$, which holds iff $x\in(0,2)$.

(c) Each solution set is a bounded interval. Among the integers $1..10^6$, R1 holds for exactly 3 ($n=1,2,3$) and R2 for exactly 1 ($x=1$). So under natural density on $\mathbb Z_{>0}$ or $\mathbb R_{>0}$, **"most numbers" is false under both readings.** It is true at $n=2$.

## N3. Speed-of-light figures (TH-008) — 33:25, 126:59

| Quoted | Exact reference | Difference |
|---|---|---|
| "299,752,400 something" (m/s implied) | $c=299\,792\,458$ m/s (exact, SI 2019) | 40 058 m/s (relative $1.336\times10^{-4}$) |
| "186,624 miles a second" | $c=18737028625/100584=186\,282.397051$ mi/s (international mile 1609.344 m, exact) | 341.60 mi/s (relative $1.834\times10^{-3}$) |
| "300 million km/second squared" | $299\,792.458$ km/s | wrong units and magnitude |

**Pattern.** $186\,624=432^2=2^8 3^6$ exactly (P-TH-1). The quoted figure is a 3-smooth square near $c$. It is not a measured value.

## N4. Element "tones" (TH-009) — 09:06, 10:22–10:28

(a) "Same key of e. 40.5." Then 81 Hz (carbon), 162 Hz (silicon), 324 Hz (cobalt), "doubled in each octave". Rhodium is named at 12:25, with no figure given.

(c) **Exact identification.** The numbers are $40.5\cdot2^k=3^4\cdot2^{k-1}$. In particular $324=256\cdot\tfrac{81}{64}$, the Pythagorean major third above scientific-pitch C = 256 Hz. So the sequence is the Pythagorean E in the $C=256$ convention.
- It is $-29.81$ cents from 12-TET $E_1$ at A4 = 440 Hz (41.2034 Hz).
- $81/64$ is $+7.82$ cents above the equal-tempered third.

The figures come from a tuning convention. They are not measurements.

(e) **Test of "same tone, doubled each octave"** against element data (mendeleev 1.3.0 bundled data). Successive ratios for H→C→Si→Co→Rh:

| Quantity | Successive ratios |
|---|---|
| $Z$ | 6.0, 2.333, 1.929, 1.667 |
| Atomic weight | 11.92, 2.338, 2.098, 1.746 |
| First ionization energy | 0.828, 0.724, 0.967, 0.946 |

For Li→Na→K: $Z$ 3.667, 1.727; atomic weight 3.313, 1.701.

None of these quantities doubles. The claim remains testable only if "tone" is defined as a specific spectral observable, for example a named line (NIST ASD). It was not, and line data were not reachable here. → TH-G06.

## N5. "These 64 circles overlapping" (Flower of Life) (TH-010) — 29:01

(c) Consider circles centred on a triangular lattice, with radius equal to the spacing, truncated by a centred disc. The achievable counts begin $1,7,13,19,31,37,43,55,61,73,\dots$. The centred hexagonal numbers are $1,7,19,37,61,91,\dots$. **64 is in neither list.**
$64=4^3=2^6$. A 64-circle figure would need a non-centred or non-planar arrangement. → TH-G07.

## References

- BIPM, *The International System of Units*, 9th ed. (2019): $c$ exact.
- International Yard and Pound Agreement (1959): 1 mi = 1609.344 m exactly.
- Benson, D. J. (2006). *Music: A Mathematical Offering*, CUP (Pythagorean tuning; cents).
- mendeleev 1.3.0 (PyPI), element data bundled with the package.
