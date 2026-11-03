from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import argparse
from itertools import combinations
from math import gcd

import numpy as np

from src.envs.hirzebruch import count_c4

N = 39
DEGREE = 6
UNITS = tuple(u for u in range(N) if gcd(u, N) == 1)


def is_sidon(S):
    """Return True iff all ordered nonzero differences a-b mod 39 are distinct."""
    seen = set()
    for a in S:
        for b in S:
            if a == b:
                continue
            d = (a - b) % N
            if d in seen:
                return False
            seen.add(d)
    return True


def affine_canonical(S):
    """Canonical representative under x -> u*x + a, u a unit mod 39."""
    images = []
    for u in UNITS:
        for a in range(N):
            images.append(tuple(sorted((u * x + a) % N for x in S)))
    return min(images)


def build_matrix(S):
    data = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        for s in S:
            data[i, (i + s) % N] = 1
    return data


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=Path("hirzebruch39_cyclic_classes"))
    args = parser.parse_args()

    sidon_sets = []
    for tail in combinations(range(1, N), DEGREE - 1):
        S = (0,) + tail
        if is_sidon(S):
            sidon_sets.append(S)

    classes = {}
    for S in sidon_sets:
        canonical = affine_canonical(S)
        classes.setdefault(canonical, 0)
        classes[canonical] += 1

    reps = sorted(classes)

    print(f"Sidon 6-subsets containing 0: {len(sidon_sets)}")
    print(f"Affine classes: {len(reps)}")

    args.output.mkdir(parents=True, exist_ok=True)

    for idx, S in enumerate(reps, start=1):
        data = build_matrix(S)
        assert np.all(data.sum(axis=0) == DEGREE)
        assert np.all(data.sum(axis=1) == DEGREE)
        assert count_c4(data) == 0

        stem = args.output / f"class_{idx:02d}"
        np.save(stem.with_suffix(".npy"), data)
        np.savetxt(stem.with_suffix(".txt"), data, fmt="%d")
        with open(stem.with_suffix(".sidon.txt"), "w", encoding="utf-8") as handle:
            handle.write(" ".join(map(str, S)) + "\n")

        print(f"{idx:02d}: S={S} normalized-count={classes[S]}")

    if len(sidon_sets) != 2640:
        raise RuntimeError(f"expected 2640 normalized Sidon sets, found {len(sidon_sets)}")
    if len(reps) != 21:
        raise RuntimeError(f"expected 21 affine classes, found {len(reps)}")


if __name__ == "__main__":
    main()
