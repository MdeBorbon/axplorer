# Hirzebruch mutation experiments

This branch explores Hirzebruch arrangements by mutating known rank-three complex reflection arrangements rather than fixing the number of lines in advance.

## Objective

For an arrangement A of N projective lines, let v(L) be the number of vertices on a line L. Define

    delta_L = 3 v(L) - N - 3,
    D(A) = sum_L delta_L^2.

For a nontrivial arrangement, D(A)=0 exactly when N=3k and every line has k+1 vertices, i.e. A has the Hirzebruch property.

The Axplorer score uses D as the primary objective. Among arrangements with equal D it rewards distance from the starting reflection arrangement.

## Reflection seeds

The first implementation uses the imprimitive rank-three reflection arrangements

- `gmm`: G(m,m,3), with 3m lines x=zeta^a y, y=zeta^a z, z=zeta^a x;
- `gm1`: G(m,1,3), obtained by adding the coordinate lines x=0, y=0, z=0.

This gives a built-in validation test: G(m,m,3) -> G(m,1,3) is a three-line extension, and the reverse operation is a three-line deletion. Both endpoints have Hirzebruch defect zero.

## Candidate lines and mutation moves

Starting from the seed, the environment forms a finite candidate pool consisting of

1. all seed lines;
2. every projective line through a pair of vertices of the seed arrangement.

The selected arrangement is allowed to change size. Local moves are:

- add three lines;
- delete three lines;
- replace three lines by three others.

The default `size_radius=3` searches the seed line count and the adjacent sizes differing by three. Increase this later to explore longer mutation paths.

A geometric prior ranks additions by how many seed-line intersections are reused at existing seed vertices. Deletions are ranked by the number of double vertices on a seed line. For m>=4 these criteria identify the three coordinate lines in the known G(m,m,3) <-> G(m,1,3) transition.

The candidate-pair intersection table is precomputed, so the inner local-search loop evaluates combinatorial defects using integer arrays rather than repeated complex cross-products.

## Laptop sanity checks

After switching to this branch, first check the known three-line extension:

```bash
python tools/mutation_probe.py --reflection_m 4 --seed_family gmm --trials 20 --steps 200
```

A successful run should report a non-seed defect-zero configuration labelled `known-neighbor`, with 15 lines.

Then check the reverse deletion:

```bash
python tools/mutation_probe.py --reflection_m 4 --seed_family gm1 --trials 20 --steps 200
```

This should recover the 12-line G(4,4,3) arrangement as a defect-zero `known-neighbor`.

## First Axplorer smoke test

Once the probe works, use a deliberately small CPU run:

```bash
python train.py --env_name hirzebruch_mutations --exp_name mutations-g4 --reflection_m 4 --seed_family gmm --size_radius 3 --local_search_steps 200 --gensize 500 --max_epochs 5 --max_steps 200 --num_samples_from_model 500 --pop_size 500 --ntest 50 --batch_size 32 --n_layer 4 --n_embd 256 --save_best true --cpu true
```

Do not scale this up yet. The first purpose is to verify that Axplorer learns the known G(4,4,3) -> G(4,1,3) move and that generation/local search remain stable on the laptop.

## Next research step

After the validation run, the known seed and known reflection neighbor should be put on a tabu list so that score-zero rediscoveries of them do not dominate the population. Then increase `size_radius`, introduce further reflection seeds, and finally add specially aligned unions of reflection arrangements to the candidate pool.
