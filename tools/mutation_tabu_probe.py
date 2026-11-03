#!/usr/bin/env python3
"""Search away from the known G(m,m,3) and G(m,1,3) endpoints."""

import argparse
from collections import Counter
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.envs.mutations_tabu import HirzebruchTabuMutationDataPoint


def configure(args):
    HirzebruchTabuMutationDataPoint.configure(
        args.reflection_m,
        args.seed_family,
        args.size_radius,
        args.min_lines,
        args.steps,
        args.initial_mutations,
        args.temperature,
        args.candidate_top_k,
        args.greedy_move_probability,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--reflection_m", type=int, default=4)
    parser.add_argument("--seed_family", choices=["gmm", "gm1"], default="gmm")
    parser.add_argument("--size_radius", type=int, default=6)
    parser.add_argument("--min_lines", type=int, default=12)
    parser.add_argument("--steps", type=int, default=1000)
    parser.add_argument("--trials", type=int, default=50)
    parser.add_argument("--initial_mutations", type=int, default=3)
    parser.add_argument("--temperature", type=float, default=30.0)
    parser.add_argument("--candidate_top_k", type=int, default=20)
    parser.add_argument("--greedy_move_probability", type=float, default=0.20)
    args = parser.parse_args()

    configure(args)
    cls = HirzebruchTabuMutationDataPoint

    print(
        f"seed={args.seed_family}(m={args.reflection_m}), "
        f"seed_lines={cls.SEED_COUNT}, pool={cls.POOL_SIZE}, "
        f"allowed_sizes={cls.SIZE_MIN}..{cls.SIZE_MAX}"
    )
    print("tabu:", sorted(cls.TABU_LABELS.values()))

    defect_histogram = Counter()
    size_histogram = Counter()
    zeroes = {}
    best = None

    for _ in range(args.trials):
        d = cls(cls.POOL_SIZE, init=True)
        d.local_search(improve_with_local_search=True)

        defect = d.defect()
        selected = d.selected_indices()
        n_lines = len(selected)
        defect_histogram[defect] += 1
        size_histogram[n_lines] += 1

        if d.is_tabu():
            continue

        record = (
            defect,
            n_lines,
            d.novelty(),
            dict(sorted(d.multiplicity_profile().items())),
            selected,
        )

        if best is None or (defect, -d.novelty()) < (best[0], -best[2]):
            best = record

        if defect == 0:
            zeroes[d.features] = record

    print("defect histogram:", dict(sorted(defect_histogram.items())))
    print("line-count histogram:", dict(sorted(size_histogram.items())))
    print(f"distinct non-tabu zeroes: {len(zeroes)}")

    for number, record in enumerate(zeroes.values(), start=1):
        defect, n_lines, novelty, profile, selected = record
        print(
            f"[{number}] NEW ZERO: lines={n_lines}, novelty={novelty}, "
            f"multiplicities={profile}"
        )
        print("    selected candidate indices:", selected)

    if best is not None and not zeroes:
        defect, n_lines, novelty, profile, selected = best
        print(
            f"best non-tabu candidate: defect={defect}, lines={n_lines}, "
            f"novelty={novelty}, multiplicities={profile}"
        )
        print("    selected candidate indices:", selected)


if __name__ == "__main__":
    main()
