# What we built and tested

Imagine someone gives us two sets of building instructions and says they might help explain how the universe works. We are checking the instructions, one step at a time.

Our goal is to find clear rules that other people can check. We do not yet have a demonstrated “source code” for the universe.

## Our math workshop

We built a workshop with five kinds of tools:

- **A source shelf:** records of where each idea came from and which version we read.
- **Careful instructions:** definitions, assumptions and written mathematical arguments.
- **Math calculators:** Python programs that recount objects, calculate formulas, try alternatives and look for contradictions.
- **A proof checker:** Lean files that let a computer check certain mathematical statements, step by step.
- **An inspection notebook:** reviews, corrections and a list of questions we have not answered.

The programs and their recorded answers are saved so another person can repeat the work.

## Two sets of instructions

### Chain Theory

One construction makes 172 mathematical pieces and uses them to calculate a particular number. We counted the pieces again and got the same number.

Then we changed some of the rules to see which choices matter. The answer depends on choices that still need an independent reason. Getting the number again establishes that we can follow the recipe; showing why nature uses it would require more evidence.

We also rebuilt an existing computer-checked proof of the mathematical read cap and connected it to the author's dimension-selection theorem. This checks those precise definitions; explaining how the proposed physical reads implement them still needs work.

We checked the larger compendium too: whether its steps depend on unfinished steps, and whether its computer proofs prove the claims written beside them. Some proofs establish useful mathematical facts. Many check narrower statements, such as whether an entered statement appears in a supplied table.

### Geometric Unity

This proposal uses complicated geometry. We wrote down explicit versions of its objects and equations, checked how the parts fit, and repaired a sign problem. We also found that the printed equation does not follow from the proposed action under the stated assumptions.

For two specified field families, we found every allowed solution in the positive Z1 version of the model. One allows only a fixed field over flat four-dimensional starting geometry, at one exact parameter relation. The other allows no solution. The larger observer geometry is curved in this example.

The harder Lambda5 family has two specified branches. We have now ruled out both in the declared positive Z1 model, including fields that change across the geometry. For the second branch, its field equations force one formula for the form’s size, while conservation requires an incompatible formula. Both requirements cannot hold together. This closes those cases; other field families and different assumptions still need investigation.

For another candidate, called the first pure hook, we checked all 896 possible directions in a simpler comparison model. An exact parameter calculation and a rescaling argument extend that result to every smooth positive starting geometry, for every parameter value that keeps the specified Z1 metric positive. The equations force this field to be zero. The remaining equation then requires a positive value for a curvature quantity that the geometry makes negative. That rules out this candidate under those assumptions. Other hook families and their intersections still need investigation.

For fields with a specified affine symmetry, we ruled out three more hook families over flat starting geometry. In one, an exact calculation leaves only two possible choices, and each makes a required equation positive instead of zero. This holds for every parameter that keeps the specified Z1 metric positive. In the other two families at one specified parameter, the derivative test leaves one and six directions; their nonlinear equations then give contradictions. Fields with less symmetry remain open.

We also identified at least 121 independent tensor directions at a point that the curvature test allows. They still need to satisfy the nonlinear equations and fit together from point to point before they can form a working field.

We also extended a known family of working fields across complete global models, including versions with a periodic four-dimensional starting space. Their geodesics—the geometry's straightest paths—have no artificial stopping point. For fixed model parameters, an exact formula gives zero, one or two constant field configurations in the regular family. These are mathematical examples within our reconstruction; distinct physical states have not been established.

Separate results exclude several families of curved starting geometry, including points with zero Ricci curvature but nonzero total curvature. Those are different ways of measuring curvature. Every result has a named set of assumptions; the [synthesis](SYNTHESIS.md#geometric-unity) gives the precise map.

## Existing work we can use

Other people provide useful programs and computer-checked mathematics that overlap ours. We have reproduced some of that work, including the read-cap proof above. The [prior-work survey](RELATED-WORK.md) gives credit and distinguishes what we reran from what we only inspected. None of the work we assessed proves either entire physical theory.

## What “the tests passed” means

The [verification guide](../computation/README.md#completed-validation) shows which Chain Theory and Geometric Unity checks were run, their results, and their limits. The [README](../README.md#reproduce-selected-checks) explains how to run the current selected suite.

Some tests confirm an answer. Others confirm that a proposed rule fails, or that a known error is detected. Passing means the program produced the expected, checked result.

This suite covers selected work. The recorded Lean checks and manuscript builds have their own scope; the selected suite does not rebuild them. Internal reviews include AI-assisted work, and external peer review remains outstanding.
