from itertools import combinations, product

import numpy as np

from src.envs.mutations import _projective_key, arrangement_incidence


def line_signatures(lines):
    """Incidence signature of each line: sorted multiplicities of its vertices."""
    _, vertices = arrangement_incidence(lines)
    signatures = [[] for _ in lines]
    for rec in vertices.values():
        mult = len(rec["lines"])
        for i in rec["lines"]:
            signatures[i].append(mult)
    return [tuple(sorted(sig)) for sig in signatures]


def _is_projective_frame(vectors, tol=1e-8):
    """Four points of the dual plane with no three collinear."""
    s0, s1, s2, s3 = vectors
    S = np.column_stack((s0, s1, s2))
    if abs(np.linalg.det(S)) < tol:
        return False
    a = np.linalg.solve(S, s3)
    return bool(np.all(np.abs(a) > tol))


def _choose_source_frame(lines, signatures, tol=1e-8):
    """Choose a projective frame whose signature classes are as restrictive as possible."""
    counts = {sig: signatures.count(sig) for sig in set(signatures)}
    best = None
    best_cost = None
    for idx in combinations(range(len(lines)), 4):
        vectors = [lines[i] for i in idx]
        if not _is_projective_frame(vectors, tol=tol):
            continue
        cost = 1
        for i in idx:
            cost *= counts[signatures[i]]
        if best is None or cost < best_cost:
            best = idx
            best_cost = cost
    return best


def find_projective_equivalence(source_lines, target_lines, tol=1e-7):
    """
    Find B in PGL(3,C) sending the source line set to the target line set.

    Line equations are treated as points in the dual projective plane.  The
    returned 3x3 matrix acts on line coefficient vectors.  ``None`` means that
    no projective equivalence was found.
    """
    if len(source_lines) != len(target_lines):
        return None

    n = len(source_lines)
    source_signatures = line_signatures(source_lines)
    target_signatures = line_signatures(target_lines)
    if sorted(source_signatures) != sorted(target_signatures):
        return None

    frame = _choose_source_frame(source_lines, source_signatures, tol=tol)
    if frame is None:
        return None

    s0, s1, s2, s3 = [source_lines[i] for i in frame]
    S = np.column_stack((s0, s1, s2))
    a = np.linalg.solve(S, s3)

    target_by_signature = {}
    for j, sig in enumerate(target_signatures):
        target_by_signature.setdefault(sig, []).append(j)

    candidate_lists = [target_by_signature[source_signatures[i]] for i in frame]
    target_keys = {_projective_key(line) for line in target_lines}

    for target_frame in product(*candidate_lists):
        if len(set(target_frame)) != 4:
            continue

        t0, t1, t2, t3 = [target_lines[j] for j in target_frame]
        T = np.column_stack((t0, t1, t2))
        if abs(np.linalg.det(T)) < tol:
            continue
        b = np.linalg.solve(T, t3)
        if np.any(np.abs(b) < tol):
            continue

        B = T @ np.diag(b / a) @ np.linalg.inv(S)
        if abs(np.linalg.det(B)) < tol:
            continue

        transformed = {_projective_key(B @ line) for line in source_lines}
        if transformed == target_keys:
            return B

    return None


def classify_imprimitive_reflection_copy(lines, max_m=50):
    """
    Classify a line arrangement against the rank-three imprimitive reflection
    families G(m,m,3) and G(m,1,3), up to projective equivalence.
    """
    from src.envs.mutations import reflection_arrangement

    n = len(lines)
    candidates = []

    if n % 3 == 0:
        m = n // 3
        if 2 <= m <= max_m:
            candidates.append((f"G({m},{m},3)", reflection_arrangement(m, "gmm")))

        m = n // 3 - 1
        if 2 <= m <= max_m:
            candidates.append((f"G({m},1,3)", reflection_arrangement(m, "gm1")))

    for label, known in candidates:
        B = find_projective_equivalence(known, lines)
        if B is not None:
            return label, B

    return None, None
