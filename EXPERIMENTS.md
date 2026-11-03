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
