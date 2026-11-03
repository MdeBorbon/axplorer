#!/usr/bin/env python3
"""Fast sanity check for the Hirzebruch mutation environment."""

import argparse
from collections import Counter

from src.envs.mutations import (
    HirzebruchMutationDataPoint,
    _projective_key,
    hirzebruch_defect,
    multiplicity_profile,
    reflection_arrangement,
)


def line_key_set(lines):
    return {_projective_key(line) for line in lines}


def configure(args):
    HirzebruchMutationDataPoint.configure(
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
    parser.add_argument("--size_radius", type=int, default=3)
    parser.add_argument("--min_lines", type=int, default=6)
    parser.add_argument("--steps", type=int, default=200)
    parser.add_argument("--trials", type=int, default=20)
    parser.add_argument("--initial_mutations", type=int, default=2)
    parser.add_argument("--temperature", type=float, default=20.0)
    parser.add_argument("--candidate_top_k", type=int, default=12)
    parser.add_argument("--greedy_move_probability", type=float, default=0.35)
    args = parser.parse_args()

    configure(args)
    cls = HirzebruchMutationDataPoint

    seed = reflection_arrangement(args.reflection_m, args.seed_family)
    neighbor_family = "gm1" if args.seed_family == "gmm" else "gmm"
    known_neighbor = reflection_arrangement(args.reflection_m, neighbor_family)
    known_neighbor_keys = line_key_set(known_neighbor)

    print(
        f"seed={args.seed_family}(m={args.reflection_m}), "
        f"seed_lines={len(seed)}, pool={cls.POOL_SIZE}, "
        f"allowed_sizes={cls.SIZE_MIN}..{cls.SIZE_MAX}"
    )
    print(
        f"seed defect={hirzebruch_defect(seed, args.min_lines)}, "
        f"multiplicities={dict(sorted(multiplicity_profile(seed).items()))}"
    )
    print(
        f"known neighbor={neighbor_family}(m={args.reflection_m}), "
        f"lines={len(known_neighbor)}, "
        f"defect={hirzebruch_defect(known_neighbor, args.min_lines)}"
    )

    hits = {}
    defect_histogram = Counter()

    for _ in range(args.trials):
        d = cls(cls.POOL_SIZE, init=True)
        d.local_search(improve_with_local_search=True)
        defect = d.defect()
        defect_histogram[defect] += 1

        if defect != 0 or d.novelty() == 0:
            continue

        selected = d.selected_indices()
        selected_lines = [cls.CANDIDATE_LINES[i] for i in selected]
        keys = line_key_set(selected_lines)
        label = "known-neighbor" if keys == known_neighbor_keys else "other-zero"
        hits[d.features] = (
            label,
            len(selected),
            d.novelty(),
            dict(sorted(d.multiplicity_profile().items())),
            selected,
        )

    print(f"defect histogram: {dict(sorted(defect_histogram.items()))}")
    print(f"distinct non-seed zeroes: {len(hits)}")

    for number, hit in enumerate(hits.values(), start=1):
        label, n_lines, novelty, profile, selected = hit
        print(
            f"[{number}] {label}: lines={n_lines}, novelty={novelty}, "
            f"multiplicities={profile}"
        )
        print("    selected candidate indices:", selected)


if __name__ == "__main__":
    main()
