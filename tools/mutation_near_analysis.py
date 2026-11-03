#!/usr/bin/env python3
"""Inspect a saved near-solution using exact Gaussian rational reconstruction.

This is intended for the m=4 pool; it rejects coordinates that cannot be
reconstructed with small rational real and imaginary parts.
"""
import argparse
from collections import Counter
from fractions import Fraction
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.envs.mutations_explore import HirzebruchExploreDataPoint as C


def rational(x):
    f = Fraction(x).limit_denominator(64)
    if abs(float(f) - x) > 1e-10:
        raise ValueError("Coordinate is not a small rational")
    return sp.Rational(f.numerator, f.denominator)


def analyse(path, defect):
    records = [json.loads(s) for s in path.read_text().splitlines()]
    p = records[0]
    r = next(r for r in records if r.get("D") == defect)
    lines = [sp.Matrix([rational(a) + sp.I*rational(b) for a, b in line])
             for line in r["lines"]]
    vertices = {}
    for i in range(len(lines)):
        for j in range(i):
            v = [sp.simplify(x) for x in lines[i].cross(lines[j])]
            pivot = next(x for x in v if x != 0)
            key = tuple(sp.simplify(x/pivot) for x in v)
            vertices.setdefault(key, set()).update((i, j))
    counts = [sum(i in s for s in vertices.values()) for i in range(len(lines))]
    errors = [v - len(lines)//3 - 1 for v in counts]
    assert errors == r["errors"]
    assert 9*sum(e*e for e in errors) == defect

    C.configure(p['reflection_m'], p['seed_family'], p['size_radius'],
                p['min_lines'], 0, p['initial_mutations'], p['temperature'],
                p['candidate_top_k'], p['greedy_move_probability'])
    d = C(C.POOL_SIZE)
    d.data[r['indices']] = 1
    original = d.data.copy()
    energy = d._local_energy()
    improving = 0
    best = None
    for remove in r['indices']:
        for add in np.flatnonzero(original == 0):
            d.data = original.copy()
            d.data[remove] = 0
            d.data[add] = 1
            e = d._local_energy()
            improving += e < energy - 1e-12
            if best is None or e < best['energy']:
                best = dict(remove=remove, add=int(add), energy=e,
                            D=d.defect(), errors=d.line_error_stats()[2].tolist())
            if d.defect() == 0:
                # Do not silently continue analysis past an unexpected zero.
                raise RuntimeError('Zero in swap neighborhood; inspect immediately')
    return dict(indices=r['indices'], N=r['N'], D=defect, errors=errors,
                exact_lines=[[str(x) for x in line] for line in lines],
                exact_multiplicities=dict(Counter(map(len, vertices.values()))),
                incidence=[dict(point=[str(x) for x in v],
                                lines=sorted(r['indices'][i] for i in s))
                           for v, s in vertices.items()],
                single_swaps=len(r['indices'])*(C.POOL_SIZE-len(r['indices'])),
                improving_single_swaps=int(improving), energy=energy, best_swap=best)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('audit', type=Path)
    parser.add_argument('--defect', type=int, required=True)
    args = parser.parse_args()
    print(json.dumps(analyse(args.audit, args.defect), indent=2))
