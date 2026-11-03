from math import gcd
from functools import reduce


N = 39
S = {0, 1, 3, 7, 12, 22}
SEED = (0, 1, 2, 5)


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


def build_incidence():
    line_points = {i: {(i + s) % N for s in S} for i in range(N)}
    point_lines = {j: {i for i in range(N) if j in line_points[i]} for j in range(N)}
    return line_points, point_lines


def main():
    line_points, point_lines = build_incidence()

    print("Seed points:", SEED)
    print("Pairwise prescribed joining lines:")
    for a in range(len(SEED)):
        for b in range(a + 1, len(SEED)):
            p, q = SEED[a], SEED[b]
            common = sorted(point_lines[p] & point_lines[q])
            print(f"  P_{p}, P_{q}: L_{common[0]}")
            if len(common) != 1:
                raise RuntimeError("seed pair does not determine a unique prescribed line")

    for i, pts in line_points.items():
        if len(pts & set(SEED)) >= 3:
            raise RuntimeError(f"three seed points lie on prescribed line L_{i}")

    print("No three seed points lie on one prescribed line.")
    print("Therefore in every projective realization the four seed points are in general position.")
    print("Normalize:")
    print("  P_0 = [1:0:0]")
    print("  P_1 = [0:1:0]")
    print("  P_2 = [0:0:1]")
    print("  P_5 = [1:1:1]")

    points = {
        0: (1, 0, 0),
        1: (0, 1, 0),
        2: (0, 0, 1),
        5: (1, 1, 1),
    }
    lines = {}
    construction = []

    while True:
        changed = False

        for i in range(N):
            if i in lines:
                continue
            known_points = sorted(line_points[i] & points.keys())
            if len(known_points) < 2:
                continue

            a, b = known_points[:2]
            line = cross(points[a], points[b])
            if line == (0, 0, 0):
                raise RuntimeError(f"coincident constructed points P_{a}, P_{b}")
            lines[i] = line
            construction.append(("line", i, a, b, line))
            changed = True

            for j in sorted(line_points[i] & points.keys()):
                value = dot(lines[i], points[j])
                if value != 0:
                    print()
                    print("CONTRADICTION")
                    print(f"Constructed L_{i} from P_{a} and P_{b}: L_{i} = {line}")
                    print(f"Prescribed incidence requires P_{j} in L_{i}.")
                    print(f"But P_{j} = {points[j]} and L_{i}(P_{j}) = {value} != 0.")
                    print(f"Construction steps used: {len(construction)}")
                    print("Hence this cyclic 39_6 incidence structure has no realization over C.")
                    return

        for j in range(N):
            if j in points:
                continue
            known_lines = sorted(point_lines[j] & lines.keys())
            if len(known_lines) < 2:
                continue

            a, b = known_lines[:2]
            point = cross(lines[a], lines[b])
            if point == (0, 0, 0):
                raise RuntimeError(f"coincident constructed lines L_{a}, L_{b}")
            points[j] = point
            construction.append(("point", j, a, b, point))
            changed = True

            for i in sorted(point_lines[j] & lines.keys()):
                value = dot(lines[i], points[j])
                if value != 0:
                    print()
                    print("CONTRADICTION")
                    print(f"Constructed P_{j} from L_{a} and L_{b}: P_{j} = {point}")
                    print(f"Prescribed incidence requires P_{j} in L_{i}.")
                    print(f"But L_{i} = {lines[i]} and L_{i}(P_{j}) = {value} != 0.")
                    print(f"Construction steps used: {len(construction)}")
                    print("Hence this cyclic 39_6 incidence structure has no realization over C.")
                    return

        if not changed:
            break

    print("No contradiction found.")
    print(f"Constructed {len(points)} points and {len(lines)} lines.")


if __name__ == "__main__":
    main()
