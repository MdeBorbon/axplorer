from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import argparse
import numpy as np

from src.envs.hirzebruch import count_c4

N = 39
S = (0, 1, 3, 7, 12, 22)


def cyclic_start():
    data = np.zeros((N, N), dtype=np.uint8)
    for i in range(N):
        for s in S:
            data[i, (i + s) % N] = 1
    return data


def try_c4free_switch(data, rng):
    r1, r2 = map(int, rng.choice(N, size=2, replace=False))
    c1 = int(rng.choice(np.flatnonzero(data[r1])))
    c2 = int(rng.choice(np.flatnonzero(data[r2])))
    if c1 == c2 or data[r1, c2] or data[r2, c1]:
        return False

    data[r1, c1] = 0
    data[r2, c2] = 0
    data[r1, c2] = 1
    data[r2, c1] = 1

    if count_c4(data) == 0:
        return True

    data[r1, c2] = 0
    data[r2, c1] = 0
    data[r1, c1] = 1
    data[r2, c2] = 1
    return False


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--count", type=int, default=20)
    parser.add_argument("--accepted-between", type=int, default=200)
    parser.add_argument("--max-attempts", type=int, default=2_000_000)
    parser.add_argument("--seed", type=int, default=0)
    parser.add_argument("--output", type=Path, default=Path("hirzebruch39_samples"))
    args = parser.parse_args()

    rng = np.random.default_rng(args.seed)
    data = cyclic_start()
    args.output.mkdir(parents=True, exist_ok=True)

    assert np.all(data.sum(axis=0) == 6)
    assert np.all(data.sum(axis=1) == 6)
    assert count_c4(data) == 0

    accepted = 0
    attempts = 0
    saved = 0

    while saved < args.count and attempts < args.max_attempts:
        attempts += 1
        if not try_c4free_switch(data, rng):
            continue
        accepted += 1
        if accepted % args.accepted_between != 0:
            continue

        saved += 1
        path = args.output / f"candidate_{saved:04d}.npy"
        np.save(path, data)
        np.savetxt(path.with_suffix(".txt"), data, fmt="%d")
        print(f"saved {path} after {accepted} accepted switches ({attempts} attempts)")

    print(f"accepted switches: {accepted}")
    print(f"attempts: {attempts}")
    print(f"saved candidates: {saved}")


if __name__ == "__main__":
    main()
