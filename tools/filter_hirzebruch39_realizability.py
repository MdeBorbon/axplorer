from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import argparse
from functools import reduce
from itertools import combinations
from math import gcd

import numpy as np

from src.envs.hirzebruch import count_c4


def primitive(v):
    if v == (0, 0, 0):
        return v
    g = reduce(gcd, (abs(x) for x in v if x != 0))
    v = tuple(x // g for x in v)
    for x in v:
        if x != 0:
            if x < 0:
                v = tuple(-y for y in v)
            break
    return v


def cross(u, v):
    return primitive(
        (
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        )
    )


def dot(u, v):
    return sum(a * b for a, b in zip(u, v))


def load_matrix(path):
    path = Path(path)
    if path.suffix == ".npy":
        data = np.load(path)
    else:
        data = np.loadtxt(path, dtype=np.uint8)
    data = np.asarray(data, dtype=np.uint8)
    if data.shape != (39, 39):
        raise ValueError(f"expected a 39 x 39 matrix, got {data.shape}")
    return data


def incidence_sets(data):
    line_points = {i: set(map(int, np.flatnonzero(data[i]))) for i in range(39)}
    point_lines = {j: set(map(int, np.flatnonzero(data[:, j]))) for j in range(39)}
    return line_points, point_lines


def pair_line(a, b, point_lines):
    common = point_lines[a] & point_lines[b]
    if len(common) > 1:
        raise ValueError("matrix is not C4-free")
    return next(iter(common)) if common else None


def triple_forced_noncollinear(triple, point_lines):
    pair_lines = []
    for a, b in combinations(triple, 2):
        line = pair_line(a, b, point_lines)
        if line is not None:
            pair_lines.append(line)
    return len(set(pair_lines)) >= 2


def seed_candidates(point_lines):
    candidates = []
    for seed in combinations(range(39), 4):
        if not all(triple_forced_noncollinear(triple, point_lines) for triple in combinations(seed, 3)):
            continue
        joins = sum(pair_line(a, b, point_lines) is not None for a, b in combinations(seed, 2))
        candidates.append((joins, seed))
    candidates.sort(reverse=True)
    return [seed for _, seed in candidates]


def same_projective(u, v):
    return primitive(u) == primitive(v)


def propagate(data, line_points, point_lines, seed):
    standard = ((1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 1))
    points = {j: p for j, p in zip(seed, standard)}
    lines = {}
    steps = []

    def contradiction(message):
        return {
            "status": "contradiction",
            "seed": seed,
            "message": message,
            "points": len(points),
            "lines": len(lines),
            "steps": len(steps),
        }

    while True:
        changed = False

        for i in range(39):
            if i in lines:
                continue
            known = sorted(line_points[i] & points.keys())
            if len(known) < 2:
                continue
            a, b = known[:2]
            line = cross(points[a], points[b])
            if line == (0, 0, 0):
                return contradiction(f"distinct prescribed points P_{a}, P_{b} coincide")
            for k, old in lines.items():
                if k != i and same_projective(line, old):
                    return contradiction(f"distinct prescribed lines L_{i}, L_{k} coincide")
            lines[i] = line
            steps.append(("line", i, a, b, line))
            changed = True

            for j, point in points.items():
                incident = dot(line, point) == 0
                required = bool(data[i, j])
                if incident != required:
                    relation = "incidence" if required else "non-incidence"
                    return contradiction(
                        f"{relation} fails for L_{i} and P_{j}: L_{i}={line}, P_{j}={point}, value={dot(line, point)}"
                    )

        for j in range(39):
            if j in points:
                continue
            known = sorted(point_lines[j] & lines.keys())
            if len(known) < 2:
                continue
            a, b = known[:2]
            point = cross(lines[a], lines[b])
            if point == (0, 0, 0):
                return contradiction(f"distinct prescribed lines L_{a}, L_{b} coincide")
            for k, old in points.items():
                if k != j and same_projective(point, old):
                    return contradiction(f"distinct prescribed points P_{j}, P_{k} coincide")
            points[j] = point
            steps.append(("point", j, a, b, point))
            changed = True

            for i, line in lines.items():
                incident = dot(line, point) == 0
                required = bool(data[i, j])
                if incident != required:
                    relation = "incidence" if required else "non-incidence"
                    return contradiction(
                        f"{relation} fails for L_{i} and P_{j}: L_{i}={line}, P_{j}={point}, value={dot(line, point)}"
                    )

        if not changed:
            break

    if len(points) == 39 and len(lines) == 39:
        return {
            "status": "realized",
            "seed": seed,
            "message": "all points and lines constructed with all incidences and non-incidences verified",
            "points": 39,
            "lines": 39,
            "steps": len(steps),
        }

    return {
        "status": "stalled",
        "seed": seed,
        "message": "join/intersection propagation did not determine the full configuration",
        "points": len(points),
        "lines": len(lines),
        "steps": len(steps),
    }


def validate_abstract(data):
    if not np.all(data.sum(axis=0) == 6):
        raise ValueError("not every point has degree 6")
    if not np.all(data.sum(axis=1) == 6):
        raise ValueError("not every line has degree 6")
    if count_c4(data) != 0:
        raise ValueError("incidence graph is not C4-free")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("path", help="candidate .npy or text incidence matrix")
    parser.add_argument("--max-seeds", type=int, default=500)
    args = parser.parse_args()

    data = load_matrix(args.path)
    validate_abstract(data)
    line_points, point_lines = incidence_sets(data)
    seeds = seed_candidates(point_lines)

    print(f"candidate: {args.path}")
    print(f"projectively normalizable seed quadruples found: {len(seeds)}")
    if not seeds:
        print("UNRESOLVED: no seed quadruple was certified in general position by the prescribed incidences")
        return

    best_stalled = None
    for idx, seed in enumerate(seeds[: args.max_seeds], start=1):
        result = propagate(data, line_points, point_lines, seed)
        if result["status"] == "contradiction":
            print("NONREALIZABLE (exact projective certificate)")
            print(f"seed: {result['seed']}")
            print(f"constructed before contradiction: {result['points']} points, {result['lines']} lines")
            print(f"construction steps: {result['steps']}")
            print(result["message"])
            return
        if result["status"] == "realized":
            print("REALIZED over Q (hence over C)")
            print(f"seed: {result['seed']}")
            print(result["message"])
            return
        if best_stalled is None or result["points"] + result["lines"] > best_stalled["points"] + best_stalled["lines"]:
            best_stalled = result

    print("UNRESOLVED by exact join/intersection propagation")
    print(f"seeds tried: {min(len(seeds), args.max_seeds)}")
    if best_stalled is not None:
        print(
            f"best closure used seed {best_stalled['seed']} and reached "
            f"{best_stalled['points']} points, {best_stalled['lines']} lines"
        )
    print("This candidate should be passed to a symbolic determinant solver.")


if __name__ == "__main__":
    main()
