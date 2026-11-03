# Hirzebruch mutation experiments

Chronological record of computational experiments on the `mutations` branch.

Keep this file concise but complete enough that another research session can continue without reconstructing the history from terminal scrollback.

## 2026-11-03 — Initial three-line mutation validation

### Question

Can the mutation environment recover the known reflection transition

    G(4,4,3) -> G(4,1,3)

by adding the three coordinate lines?

### Command

```bash
python tools/mutation_probe.py --reflection_m 4 --seed_family gmm --trials 20 --steps 200
```

### Result

Seed:

- family: `G(4,4,3)`;
- lines: 12;
- candidate pool: 63;
- allowed sizes: 9..15;
- defect: 0;
- multiplicities: `{3: 16, 4: 3}`.

Observed outcomes:

- defect histogram: `{0: 19, 45: 1}`;
- distinct non-seed zeroes: 1.

The unique non-seed zero was the known `G(4,1,3)` extension:

- lines: 15;
- novelty: 3;
- multiplicities: `{2: 12, 3: 16, 6: 3}`;
- selected candidate indices:

      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 18, 19, 42]

### Interpretation

The projective geometry, defect function, candidate pool, and three-line mutation mechanism are working. This is a validation result, not a discovery.

---

## 2026-11-03 — Exact-set tabu experiment

### Question

If the exact coordinate realizations of `G(4,4,3)` and `G(4,1,3)` are made repulsive, does pure local mutation search find another zero or explore larger line counts?

### Command

```bash
python tools/mutation_tabu_probe.py
```

Default parameters:

- seed: `G(4,4,3)`;
- allowed sizes: 12..18;
- trials: 50;
- steps: 1000.

### Result

- defect histogram: `{27: 42, 45: 8}`;
- line-count histogram: `{12: 50}`;
- distinct non-tabu zeroes: 0.

Best candidate:

- defect: 27;
- lines: 12;
- novelty: 2;
- multiplicities: `{2: 7, 3: 13, 5: 2}`;
- selected candidate indices:

      [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 18]

### Structural analysis

For `N=12`, the target is `v(L)=5`. Since

    D = 9 * sum_L (v(L)-5)^2 = 27,

exactly three lines have one extra vertex and the remaining nine lines are correct.

The candidate is a partial version of the known reflection extension: candidate `18` is a coordinate line while one seed line has been removed. The three bad lines are the remaining members of the associated pencil.

### Interpretation

The tabu penalty excluded only exact candidate-index sets. The search still hugged the known `G(4,4,3) <-> G(4,1,3)` basin and never explored `N=15` or `N=18` effectively.

This motivated two changes:

1. fixed-size one-line swaps;
2. a size-neutral exploratory objective based primarily on the fraction of defective lines rather than raw `D`.

---

## 2026-11-03 — Size-neutral exploratory search

### Question

Can a size-neutral objective plus one-line swaps leave the original reflection basin and find other exact or near-exact Hirzebruch arrangements?

### Command

```bash
python tools/mutation_explore_probe.py
```

Default parameters:

- seed: `G(4,4,3)`;
- candidate pool: 63;
- allowed sizes: 12..18;
- trials: 50;
- steps: 1500.

### Result

Defect histogram:

    {0: 1, 27: 4, 45: 1, 189: 1, 315: 1, 405: 1,
     468: 2, 486: 1, 576: 37, 792: 1}

Line-count histogram:

    {12: 48, 15: 2}

Bad-line histogram:

    {(0, 15): 1, (3, 12): 4, (4, 12): 37,
     (5, 12): 7, (6, 15): 1}

One non-tabu zero was found:

- lines: 15;
- novelty: 15;
- defect: 0;
- multiplicities: `{2: 12, 3: 16, 6: 3}`;
- selected candidate indices:

      [3, 4, 5, 7, 9, 11, 13, 19, 24, 33, 39, 50, 53, 58, 61]

Best `N=12` candidate:

- bad lines: 3/12;
- mean-square vertex error: 0.25;
- defect: 27;
- novelty: 2;
- multiplicities: `{2: 7, 3: 13, 5: 2}`;
- selected indices:

      [0, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 19]

### Classification of the 15-line zero

The zero was checked by comparing the full incidence structure and then solving for a projective transformation in the dual projective plane.

Conclusion:

    projectively equivalent to G(4,1,3)

It is therefore a disguised projective copy of the known reflection arrangement, not a new Hirzebruch arrangement.

### Interpretation

This is stronger validation than the initial experiment. The search can move far from the original coordinate realization (`novelty=15`) and reconstruct the same reflection arrangement elsewhere in the 63-line candidate pool.

However, exact-set tabu is insufficient. Known reflection arrangements must be excluded up to projective equivalence.

---

## 2026-11-03 — Projective-orbit exclusion added

### Code change

Added projective-equivalence machinery in

    src/envs/projective_equivalence.py

for the imprimitive rank-three reflection families

    G(m,m,3), G(m,1,3).

The exploratory environment now checks exact zeroes for projective equivalence to these families and penalises known reflection copies inside the local-search energy.

The expensive equivalence test is only invoked when `D=0`.

### Important limitation

The automatic classifier does NOT yet cover all exceptional rank-three Shephard--Todd reflection arrangements.

An exact zero labelled `UNCLASSIFIED` must therefore be investigated before any claim of novelty.

### Next experiment

Run a broader exploratory search after pulling the projective-orbit exclusion:

```bash
python tools/mutation_explore_probe.py --trials 100 --steps 2500
```

Primary questions:

1. Does the search now continue past projective copies of `G(4,1,3)`?
2. Are any exact zeroes found that are unclassified by the current imprimitive-family classifier?
3. If no new zero is found, what are the best near-zero configurations for `N=12`, `15`, and `18`?
4. Do near-zeroes repeatedly exhibit the same local obstruction, suggesting a new geometry-aware mutation move?

If the local search produces diverse near-zeroes but no genuinely new zero, the next stage should likely be a small Axplorer transformer run rather than indefinite tuning of the hand-written annealer.

---

## 2026-11-03 — Orbit validation, audited probes, and a two-pencil trap

Questions: does orbit exclusion penalise the disguised reflection copy; does
widening the size range escape the lower boundary; what causes repeated D=576?
All work stayed on `mutations`. No objective, pool, or mutation rule changed.

Commands (the first ran before audit instrumentation):

```bash
python tools/mutation_explore_probe.py --trials 6 --steps 300
python -u tools/mutation_explore_probe.py --trials 8 --steps 700 --output experiments/mutations-orbit/m4-r6.jsonl
python -u tools/mutation_explore_probe.py --trials 6 --steps 700 --size_radius 9 --min_lines 15 --output experiments/mutations-orbit/m4-r9-min15.jsonl
python tools/mutation_near_analysis.py experiments/mutations-orbit/m4-r6.jsonl --defect 576 > experiments/mutations-orbit/two-pencil-analysis.json
python -m unittest tools.test_mutation_audit
```

Common search parameters: `m=4`, `gmm`, pool 63, initial mutations 4,
temperature 0.15, candidate top-k 20, greedy probability 0.25.
The first two runs used radius 6, minimum 12 (allowed N=12,15,18).
Runs use entropy-seeded NumPy generators; commands reproduce the setup, not
the stochastic trajectories. Saved coordinates preserve the actual candidates.

### Orbit exclusion and recording

The previously reported disguised 15-line zero is not exact-set tabu, but is
classified `G(4,1,3)` and receives exactly +10 energy (total 9.9999765625).
Its indices, coordinates and profile are in `orbit-validation.json` in the
artifact directory above. Regression tests verify classification, penalty,
deduplicated auditing, coordinates, and save-before-stop behavior with an
injected uncertain classifier result (this is not a real unclassified zero).
Both tests passed.

The probe previously only checked final endpoints. It now saves every distinct
zero evaluated through the energy function, including rejected proposals and
exact-set tabu zeroes, and stops immediately on `UNCLASSIFIED`. JSONL files also
save parameters and every trial endpoint, with equations, per-line errors,
defect, bad-line fraction and multiplicities. An omitted `--output` creates a
timestamped file; existing files are never overwritten. Production training
behavior is unchanged; the observer belongs to the probe subclass.

### Search results

- Initial 6x300 run: all endpoints N=12; defects
  `{81:1,162:1,270:1,504:1,576:1,999:1}`. No final zero.
  Best by bad-line fraction then MSE: D=576, bad=4/12,
  profile `{2:32,4:1,8:1}`, indices
  `[2,3,9,10,15,21,22,34,36,47,51,52]`.
  Intermediate zeroes were not audited in this initial run.
- Audited 8x700 run: all endpoints N=12; defects
  `{27:1,45:1,486:1,504:1,576:3,792:1}`.
  Best D=27, bad=3/12, profile `{2:7,3:13,5:2}`, indices
  `[0,2,3,4,5,6,7,8,9,10,11,19]`; errors
  `[0,0,0,1,0,0,1,0,0,1,0,0]`. This repeats the known reflection perturbation.
  One distinct evaluated zero: `G(4,4,3)`, indices `0..11`,
  profile `{3:16,4:3}`, bad=0/12. No final or unclassified zero.
- Wider run, allowed N=15,18,21: all six endpoints N=15; defects
  `{126:2,135:1,243:1,360:1,477:1}`.
  Best D=126, bad=5/15, profile `{2:17,3:11,4:5,5:1,6:1}`, indices
  `[0,1,2,4,5,6,7,8,9,10,11,18,20,31,42]`; errors
  `[0,0,1,0,0,2,0,1,0,0,0,0,2,2,0]`.
  One distinct evaluated zero: `G(4,1,3)`, indices
  `[0,1,2,3,4,5,6,7,8,9,10,11,18,19,42]`, profile
  `{2:12,3:16,6:3}`, bad=0/15. No final or unclassified zero.
  Larger allowed sizes did not produce larger retained endpoints; this does
  not establish that the trajectories never visited them.

### Structural finding and stopping decision

Three of eight audited N=12 endpoints had D=576, bad=4/12 and profile
`{2:32,4:1,8:1}`. They consist of an eight-line pencil and a four-line pencil,
with 32 distinct cross-intersections. For the analysed representative:

- indices `[2,6,23,24,29,39,44,49,53,55,56,59]`;
- eight-line pencil `[2,6,23,29,39,44,55,56]` at `[1:-1:1]`;
- four-line pencil `[24,49,53,59]` at `[1:-i:i]`;
- errors `[0,0,0,4,0,0,0,4,4,0,0,4]`;
- eight lines have five vertices; the other four have nine vertices.

Small Gaussian-rational reconstruction of its saved coefficients independently
confirms the complete incidence lattice and errors using exact SymPy arithmetic.
Exact equations and incidence sets are saved in `two-pencil-analysis.json`.
This verifies the reconstructed concrete arrangement without rounded keys;
it is not a claim that arbitrary floating-point reconstruction is certified.
All 612 one-line swaps in the fixed pool were checked: none reduces the current
energy; equal-energy neighbors exist. No zero occurred in that neighborhood.

The mechanism generalises: two pencils of sizes 2k and k with distinct generic
cross-intersections have N=3k, 2k good lines, k lines with vertex error k,
bad fraction 1/3 and D=9k^3. Thus a stable bad fraction can conceal a large
structured defect. This explains a repeated attraction of the current objective,
not an obstruction to Hirzebruch arrangements in general.

Stopped further searches to report this persistent structural trap, as required.
No other seeds or transformer runs were attempted. Next proposed experiment:
test bounded coordinated swaps that break the pencil structure, comparing escape
and retained geometry with the existing moves while preserving the size-neutral
objective. Merely increasing trial counts is not supported by these results.

---

## 2026-11-03 — Exact four-for-four enumeration: no completion

Question: can the saved two-pencil arrangement reach N=12, D=0 by four
replacements? No stochastic search was run. The original N=12 arrangement has
D=576, four bad lines, and profile `{2:32,4:1,8:1}`. Its indices and equations
are unchanged from the preceding entry.

```bash
python -u tools/mutation_four_escape.py --output experiments/mutations-four-escape-verified
python -m unittest tools.test_mutation_four_escape
```

The script rebuilds the m=4 candidate pool over Q(i) and checks exact set equality
with all 63 reconstructed, indexed candidates. Floating-point coordinates assign
the established indices only; all incidence and enumeration decisions use exact
Gaussian rationals. It also checks the input's twelve exact equations.
The output directory must be new. Two regression tests passed, covering the
complete enumeration and exact certification of a known projective equivalence.

### Obstructions and completeness

Write A=`[2,6,23,29,39,44,55,56]` for the eight-line pencil at P=`[1:-1:1]`,
and B=`[24,49,53,59]` for the four-line pencil at Q=`[1:-i:i]`.

1. With at most three deletions, a surviving B-line still meets at least five
   surviving A-lines plus Q, so already has at least six vertices. Adding lines
   cannot remove vertices. Exact enumeration confirms this for all
   `1+12+66+220=299` subsets with zero through three deletions.
2. Among four deletions, any mixed choice leaves a B-line with at least six
   vertices. Exactly 424 of the 495 deletion sets are excluded this way.
3. Deleting all of B leaves eight concurrent lines. Any new line not through P
   has eight distinct intersections with these lines, exceeding the target five.
   If every new line passes through P, all twelve lines have only one vertex.
   This proves impossibility even for additions outside the candidate pool.
   Computationally, of the 51 pool lines outside the original arrangement,
   50 have eight intersections with A; only candidate 7 passes through P.
4. The remaining 70 deletion sets retain four A-lines and all of B. Each of
   these eight lines has five vertices: its pencil centre and four grid points.
   Every added line must preserve those vertex sets. A line through neither
   centre must be a collinear grid transversal, one point per row and column.
   A distinct line through either centre must pass through both, since otherwise
   it creates new vertices on the opposite pencil. Thus the only exception to
   the grid-transversal rule is PQ, candidate 7, with equation `y+z=0`.

### Exhaustive grid results

For all 70 survivor sets, the script enumerates admissible pool lines, verifies
their row-column permutations, and independently constructs all transversals
from all 24 permutations (1,680 exact collinearity checks in total).

| Surviving A-lines | Only transversal | Equation | Columns in B order |
|---|---|---|---|
| `[23,29,39,56]` | 51 | `x+(1-i)y/2+(-1+i)z/2=0` | `[3,0,1,2]` |
| `[23,44,55,56]` | 19 | `x=0` | `[2,1,0,3]` |

The other 68 survivor sets admit no transversal. The unrestricted geometric
enumeration gives exactly the same result: there are no additional transversals
outside the pool. Hence there are **zero transversal quartets**, and also zero
quartets formed from PQ plus three transversals. The quartet-test loops are
exhaustive but empty, rather than a sampled failure to find a zero.

No exact zeroes, unclassified zeroes, or new near-zeroes were produced. Counts
of concurrent 3-net, two-pair, and other zero patterns are all zero. There is no
best completed candidate or multiplicity profile to report for this enumeration.
The script saves exact coordinates and full lattices before classifying any zero,
certifies a returned projective correspondence symbolically, and stops on an
unclassified zero; these branches found no candidates in this run.

`experiments/mutations-four-escape-verified/enumeration.json` contains all exact
pool equations, deletion certificates, all 70 grids and admissible lists,
unrestricted transversal equations, and the empty zero list.

Conclusion: **no four-for-four escape exists for this concrete arrangement,
even allowing arbitrary new complex projective lines**. The lower bound on a
fixed-N escape is at least five replacements; existence at five is not established.
Next proposed investigation: derive the surviving-line constraints for five
deletions before considering another finite enumeration. No further search run.

---

## 2026-11-03 — Exact five-for-five obstruction

Question: can five replacements of the saved N=12 two-pencil arrangement yield
D=0? Source remains D=576, bad=4/12, profile `{2:32,4:1,8:1}`. No stochastic
search, pool expansion, or production objective change was made.

```bash
python tools/mutation_five_escape.py --output experiments/mutations-five-escape.json
python -m unittest tools.test_mutation_five_escape
```

The certificate reads the exact source equations and previous four-replacement
certificate, records their SHA-256 hashes, and uses Q(i) arithmetic exclusively.
The script's initial transversal-count calculation was first checked interactively
using the same 24 injections per core; the command above is its saved, reproducible
implementation. This is an unrestricted projective calculation, not a pool search.

### Geometric reduction

Of the 792 five-deletion subsets, 448 already leave a line with more than five
vertices: types (2,3) and (3,2). In type (2,3) only one B-line survives, so Q is
no longer a vertex of the core: that line has **six**, not seven, vertices.
This corrects the formula `9-a` in this boundary case without changing exclusion.
The eight subsets of type (1,4) leave seven concurrent A-lines. An added line
off P has at least seven vertices, whereas all additions through P give a pencil.
Thus the only remaining cores are 280 of type (4,1) and 56 of type (5,0).

Unify both cases as three **saturated** lines S through C and four **unsaturated**
lines U through D. Their vertex counts are `[5,5,5,4,4,4,4]`. C and D exchange
the roles of P and Q between the two cases. Each added line is of one of these
three types:

1. **Off C:** it must meet each S-line at a grid point. Those three points lie
   on distinct U-lines; otherwise the added line equals a surviving U-line.
   Enumerate the `4*3*2=24` injections from S to U. Two points determine the
   line and exact collinearity tests the third. All such lines lie over Q(i),
   even though arbitrary complex additions are allowed.
2. **CD:** the centre connector creates no new vertex on any survivor.
3. **Through C but different from CD:** it creates exactly one new point on
   every U-line. That point cannot be an old grid point (which would force it
   to equal an S-line) or D. Therefore at most **one** such line can be added:
   two distinct C-lines meet each U-line at two distinct new points.

If there are t off-C transversals, there can consequently be at most `t+2`
added lines. Five additions require t>=3. Even the apparent continuous freedom
in type 3 would become finite in a completion: an off-C transversal misses one
U-line and creates a new point X there, forcing the type-3 line to be CX if used.
Here the stronger counting obstruction already settles every core.

### Complete counts and certificate

| Deletion type | Cores with t=0 | t=1 | t=2 | t>=3 |
|---|---:|---:|---:|---:|
| (4,1) | 232 | 46 | 2 | 0 |
| (5,0) | 46 | 10 | 0 | 0 |

All 336 cores and 8,064 injection collinearity checks are recorded. Each core's
certificate contains saturated/opposite indices, exact transversal equations,
their three grid incidences, omitted U-line, new point there, and upper bound
on added lines. In every case that upper bound is at most four. Consequently
the set of admissible five-line completions is empty: zero quintets need full
defect evaluation, zero exact zeroes occur, and no classification is applicable.
No new near-zero is produced or ranked. The script fails explicitly if the
counting bound is ever inconclusive rather than reporting a false exclusion.

Regression checks cover the 792 deletion subsets and 336 cores, independently
construct all off-centre transversals from pairs of old saturated-line vertices,
and check the two-new-centre-lines obstruction on an exact synthetic grid.
All three tests passed.

**Conclusion:** combined with the previous certificate, the replacement distance
from this concrete arrangement to any 12-line Hirzebruch arrangement is at least
**six**, allowing arbitrary new complex lines. Existence at six is not established.
The exact five-for-five problem is settled; objective assessment follows separately.

### Multiplicity-excess assessment after the exact result

```bash
python tools/mutation_multiplicity_assess.py --output experiments/mutations-multiplicity-assessment.json
```

This deterministic diagnostic evaluates all 612 fixed-N one-line swaps using
exact pool intersections and rational score arithmetic. It changes no search
code and runs no annealing. Its hardcoded score constants match the current
exploratory energy (0.02 times MSE, with the existing novelty tie-breaker).

For a non-pencil Hirzebruch arrangement, any vertex p of multiplicity m has a
line avoiding p. That line meets the m incident lines at distinct points, proving
`m <= k+1`. More quantitatively, when `m>k+1`, all N-m lines avoiding p are bad,
and `D >= 9(N-m)(m-k-1)^2`. Thus multiplicity excess detects a real geometric
obstruction, though it is not sufficient for the Hirzebruch property and its
vanishing says nothing about novelty or reflection equivalence.

Let `M=sum_p max(0,m_p-(k+1))^2`. The trap has M=9 and relative excess
`M/(k+1)^2=9/25`. The diagnostic compares
`E_lambda = existing_energy + lambda*M/(k+1)^2` at N=12. Normalization avoids
using an unscaled quadratic term when comparing different line counts; it does
not itself enforce primary bad-fraction ordering.

- No swap lowers M without increasing the bad-line fraction.
- At lambda=0, 1/50, 1/20, 1/10, 1/4, or 1, no swap improves E_lambda.
- The first strict improvement requires
  `lambda > 286397/192000 = 1.491651041666...`.
  A witness removes candidate 2 and adds 19: N=12, D=387, bad=8/12,
  profile `{2:27,3:4,4:1,7:1}`, M=4. Its indices, error vector, and exact energy
  are saved in the assessment artifact; its equations are in the saved pool.
- At lambda=2, eight swaps improve the weighted energy; all have eight bad lines,
  twice the starting count. No exact zero occurs in this neighborhood.

Assessment: M is a justified structural diagnostic and a useful secondary
tie-breaker, but a modest secondary penalty alone does not untrap this state.
A sufficiently strong additive penalty would sacrifice the intended bad-line
fraction priority. A strict secondary/lexicographic term cannot prefer these
eight-bad-line moves over the four-bad-line starting state.

The geometric information is also relevant to proposals: the current guided
single-swap move preferentially removes bad lines, which here belong to the
four-line pencil. The excessive multiplicity is supported by the eight *good*
lines instead. Prefer testing geometry-guided coordinated moves targeting that
pencil, with multiplicity excess recorded as a diagnostic, before enabling an
energy change. A direct zero requires at least six replacements, but this does
not preclude a sequence of smaller moves through intermediate arrangements.
Production energy and moves remain unchanged. Next proposed work is a bounded
deterministic study of such moves, not larger annealing runs or an assumed
six-replacement completion.
