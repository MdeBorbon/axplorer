# AGENTS.md

## Scope

These instructions apply to the repository when working on the `mutations` branch, whose purpose is to explore Hirzebruch line arrangements using Axplorer and geometry-aware mutation search.

Do not merge this work to another branch unless the user explicitly asks. Make conservative commits to `mutations` as experiments or infrastructure improvements are validated.

## Mathematical objective

For a complex projective line arrangement `A` with `N` lines, let `v(L)` be the number of distinct intersection points (vertices) lying on a line `L`.

The Hirzebruch property is

    N = 3k,
    v(L) = k + 1 for every L in A.

Equivalently define

    delta_L = 3 v(L) - N - 3,
    D(A) = sum_L delta_L^2.

Then `D(A) = 0` exactly when the arrangement has the Hirzebruch property.

The research goal is to discover genuinely new complex Hirzebruch arrangements, or to extract structural information from persistent near-solutions. The computational goal is not merely to optimise a benchmark score.

## Known reflection arrangements are not discoveries

Known rank-three complex reflection arrangements may satisfy the Hirzebruch property. Rediscovering them is useful as validation but is not a new mathematical result.

In particular, configurations projectively equivalent to

- `G(m,m,3)`;
- `G(m,1,3)`;

must be recognised and excluded as discovery targets.

The current projective-equivalence classifier for these imprimitive families is implemented in:

    src/envs/projective_equivalence.py

An `UNCLASSIFIED` zero is NOT automatically a counterexample to Hirzebruch's conjecture. The classifier does not yet cover every exceptional rank-three Shephard--Todd reflection arrangement. If an unclassified zero appears, stop large-scale experimentation and analyse it carefully before making any claim.

## Current experiment architecture

Important files:

- `MUTATIONS.md` -- overview of the mutation-search idea and initial validation commands.
- `EXPERIMENTS.md` -- chronological experimental log. Keep this updated.
- `src/envs/mutations.py` -- base mutation environment and reflection seeds.
- `src/envs/mutations_tabu.py` -- exact-set tabu search for known reflection endpoints.
- `src/envs/mutations_explore.py` -- exploratory objective and one-line swap moves.
- `src/envs/projective_equivalence.py` -- projective-equivalence classification.
- `tools/mutation_probe.py` -- initial validation probe.
- `tools/mutation_tabu_probe.py` -- first tabu probe.
- `tools/mutation_explore_probe.py` -- current exploratory probe.

The base candidate pool for the first experiments is generated from a reflection seed by taking the seed lines together with all projective lines through pairs of seed vertices.

Current mutation moves include:

- one-line-for-one-line swaps at fixed `N`;
- add three lines;
- delete three lines;
- three-out/three-in replacements.

Changes in `N` should respect the condition that a Hirzebruch arrangement has `N` divisible by 3.

## Established validation results

The following facts have already been checked experimentally and should not be rediscovered repeatedly as if new:

1. Starting from `G(4,4,3)`, the original mutation environment reliably finds the three coordinate lines and reaches `G(4,1,3)`.
2. Starting from `G(4,1,3)`, the reverse three-line deletion recovers `G(4,4,3)`.
3. A tabu run that only excluded the exact coordinate realizations became trapped near the known reflection transition, showing that exact-set tabu is too weak.
4. The exploratory environment found a 15-line defect-zero configuration with candidate indices

       [3, 4, 5, 7, 9, 11, 13, 19, 24, 33, 39, 50, 53, 58, 61]

   and multiplicity profile

       {2: 12, 3: 16, 6: 3}.

   This arrangement was checked and is projectively equivalent to `G(4,1,3)`. It is therefore validation, not a discovery.
5. Projective-equivalence classification for the imprimitive families has now been incorporated into the exploratory search so that exact reflection zeroes are penalised up to projective equivalence, not merely by exact candidate indices.

## Experimental workflow

Before changing code, read `EXPERIMENTS.md` and the relevant environment files to understand what has already been tried.

For each experiment:

1. State the mathematical or algorithmic question being tested.
2. Use the smallest laptop-appropriate run that can answer it.
3. Record the exact command and parameters in `EXPERIMENTS.md`.
4. Record at least:
   - line-count range;
   - best defect `D`;
   - number/fraction of bad lines;
   - multiplicity profile of the best candidate;
   - exact zeroes found;
   - classification of zeroes;
   - important qualitative observations.
5. Only scale a run after a smaller experiment shows that additional computation is likely to be informative.
6. Commit validated code changes and the corresponding experiment-log update together when practical.

Do not run large searches blindly. If a search repeatedly falls into the same near-solution pattern, inspect the geometry of that obstruction and improve the mutation moves or objective before simply increasing the number of trials.

## Exact zero protocol

Whenever `D = 0` is found:

1. Save the complete candidate-index set and concrete line equations/coordinates.
2. Compute the multiplicity profile.
3. Compare projectively against all currently implemented known reflection families.
4. Compare the full intersection lattice/incidence graph against known candidates when practical.
5. If it is a known reflection copy, record it as validation and continue searching.
6. If it is `UNCLASSIFIED`, stop and investigate before further large runs.

For an unclassified zero, inspect at least:

- full incidence lattice;
- multiplicity vector `(t_2,t_3,...)`;
- per-line multiplicity signatures;
- automorphism group or obvious symmetries when practical;
- projective equivalence to known exceptional rank-three reflection arrangements;
- numerical/projective stability of the incidences;
- whether the zero depends on approximate floating-point coincidences.

Do not describe an unclassified zero as a counterexample, new arrangement, or theorem until those checks have been carried out.

## Near-zero protocol

Near-solutions can be mathematically informative. Prefer recording the per-line error vector or at least the number of bad lines, not only the scalar defect.

For `N = 3k`, the target on every line is `v(L) = k+1`. A candidate with only a few bad lines may reveal a rigid geometric obstruction and deserves inspection even if its total `D` is not the smallest across different values of `N`.

The exploratory objective therefore prioritises the fraction of defective lines, with defect magnitude secondary. Preserve this size-neutral principle unless an experiment gives a specific mathematical reason to change it.

## Search directions

After validating the current projective-orbit exclusion, reasonable directions include:

- more trials and moderate step increases for `m=4`;
- other imprimitive reflection seeds `G(m,m,3)` and `G(m,1,3)`;
- wider size radii;
- improved geometry-aware one-line or multi-line mutations;
- dynamically enlarging the candidate pool from promising near-solutions;
- specially aligned unions of reflection arrangements followed by deletion/mutation search;
- extending the known-arrangement classifier to exceptional rank-three Shephard--Todd reflection arrangements;
- eventually handing a productive local-search distribution to the Axplorer transformer rather than relying only on hand-written annealing.

Do not attempt all of these at once. Use the experimental log to choose the next test from the evidence of the previous one.

## Axplorer usage

The eventual purpose of the local mutation search is to provide Axplorer with a useful population and local sculptor. Do not spend indefinitely optimising the hand-written annealer if it already produces diverse, mathematically meaningful near-solutions.

Once the local search is demonstrably exploring beyond known reflection basins, it is appropriate to run the Axplorer training loop and let the model learn global patterns from the selected population.

Keep first transformer runs small and CPU-friendly on the user's laptop. Increase population size, epochs, or model size only when the smaller run provides evidence that the learned generator is adding value over the local search alone.

## Numerical care

The project currently uses complex floating-point representatives and approximate projective keys. Exact or near-exact incidences must be treated cautiously.

For any potentially important zero:

- recompute incidences independently;
- check sensitivity to tolerances;
- when possible recover exact algebraic coordinates or verify the claimed incidences symbolically;
- never rely solely on rounded projective keys for a mathematical claim.

## Communication and stopping conditions

Continue routine computational work autonomously, but stop and report to the user when any of the following occurs:

- an `UNCLASSIFIED` exact zero appears;
- a persistent near-zero pattern suggests a clear mathematical obstruction;
- the current search architecture appears systematically unable to leave a known basin;
- a substantial new search direction would require changing the mathematical model rather than just implementation details;
- a computation would become materially expensive on the laptop;
- there is ambiguity about whether a candidate is genuinely new.

When reporting, give the mathematical facts first: `N`, defect, number of bad lines, multiplicity profile, classification status, and candidate indices/coordinates. Avoid overclaiming.

## Code quality

Keep the experimental code readable and auditable. Prefer small functions with explicit mathematical meaning. Add concise comments where an implementation choice encodes a mathematical fact.

Do not silently change the definition of the Hirzebruch defect, projective equivalence, candidate pool, or discovery criterion. Document any such change in `EXPERIMENTS.md` and explain why it is mathematically justified.
