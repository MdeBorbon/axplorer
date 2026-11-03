#!/usr/bin/env python3
"""Probe the size-neutral exploratory Hirzebruch mutation search."""

import argparse
from collections import Counter, defaultdict
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from src.envs.mutations_explore import HirzebruchExploreDataPoint


def configure(args):
    HirzebruchExploreDataPoint.configure(
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
    parser.add_argument("--steps", type=int, default=1500)
    parser.add_argument("--trials", type=int, default=50)
    parser.add_argument("--initial_mutations", type=int, default=4)
    parser.add_argument("--temperature", type=float, default=0.15)
    parser.add_argument("--candidate_top_k", type=int, default=20)
    parser.add_argument("--greedy_move_probability", type=float, default=0.25)
    args = parser.parse_args()

    configure(args)
    cls = HirzebruchExploreDataPoint

    print(
        f"seed={args.seed_family}(m={args.reflection_m}), seed_lines={cls.SEED_COUNT}, "
        f"pool={cls.POOL_SIZE}, allowed_sizes={cls.SIZE_MIN}..{cls.SIZE_MAX}"
    )
    print("tabu:", sorted(cls.TABU_LABELS.values()))

    defect_histogram = Counter()
    size_histogram = Counter()
    bad_fraction_histogram = Counter()
    zeroes = {}
    best_by_size = {}

    for _ in range(args.trials):
        d = cls(cls.POOL_SIZE, init=True)
        d.local_search(improve_with_local_search=True)

        defect = d.defect()
        selected = d.selected_indices()
        n_lines = len(selected)
        bad_fraction, mse, errors = d.line_error_stats()
        n_bad = sum(err != 0 for err in errors)

        defect_histogram[defect] += 1
        size_histogram[n_lines] += 1
        bad_fraction_histogram[(n_bad, n_lines)] += 1

        if d.is_tabu():
            continue

        record = (
            bad_fraction,
            mse,
            defect,
            n_lines,
            d.novelty(),
            n_bad,
            dict(sorted(d.multiplicity_profile().items())),
            selected,
        )

        old = best_by_size.get(n_lines)
        if old is None or record[:3] < old[:3]:
            best_by_size[n_lines] = record

        if defect == 0:
            zeroes[d.features] = record

    print("defect histogram:", dict(sorted(defect_histogram.items())))
    print("line-count histogram:", dict(sorted(size_histogram.items())))
    print("(bad lines, total lines) histogram:", dict(sorted(bad_fraction_histogram.items())))
    print(f"distinct non-tabu zeroes: {len(zeroes)}")

    for number, record in enumerate(zeroes.values(), start=1):
        bad_fraction, mse, defect, n_lines, novelty, n_bad, profile, selected = record
        print(
            f"[{number}] NEW ZERO: lines={n_lines}, novelty={novelty}, "
            f"multiplicities={profile}"
        )
        print("    selected candidate indices:", selected)

    print("best non-tabu candidate by line count:")
    for n_lines in sorted(best_by_size):
        bad_fraction, mse, defect, _, novelty, n_bad, profile, selected = best_by_size[n_lines]
        print(
            f"  N={n_lines}: bad={n_bad}/{n_lines}, mse={mse:.4f}, "
            f"D={defect}, novelty={novelty}, multiplicities={profile}"
        )
        print("      selected candidate indices:", selected)


if __name__ == "__main__":
    main()
