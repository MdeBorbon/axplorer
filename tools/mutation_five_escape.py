#!/usr/bin/env python3
"""Certify the five-replacement obstruction over Q(i), without a pool restriction.

Three saturated pencil lines force an off-centre addition to be a transversal
of three grid points on distinct opposite-pencil lines. Through their centre,
at most one new line besides the centre connector can be added: two create two
new vertices on each unsaturated survivor. Thus at least three transversals are
necessary for five additions. The saved geometry has at most two in every core.
"""
import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.mutation_four_escape import K, cross, exact, strings

SOURCE = Path('experiments/mutations-orbit/two-pencil-analysis.json')
PREVIOUS = Path('experiments/mutations-four-escape-verified/enumeration.json')


def incident(line, point):
    return sum(a*b for a,b in zip(line, point)) == K.zero


def line_vertices(lines):
    return [{cross(line, other) for other in lines if other != line} for line in lines]


def transversals(saturated, opposite):
    """All off-centre lines preserving the three saturated lines' vertices.

    There are 4*3*2 injections into the opposite pencil. The first two points
    determine a line, so no continuous or non-Q(i) solutions can be missed.
    """
    result = {}
    for permutation in permutations(range(4), 3):
        points = [cross(s, opposite[j]) for s,j in zip(saturated, permutation)]
        line = cross(points[0], points[1])
        if incident(line, points[2]):
            assert line not in result
            result[line] = permutation
    return result


def certify(source=SOURCE, previous=PREVIOUS):
    saved = json.loads(source.read_text())
    prior = json.loads(previous.read_text())
    assert prior['complete'] and not prior['zeroes'] and len(prior['cases']) == 70
    lines = {i:exact(v) for i,v in zip(saved['indices'],saved['exact_lines'])}
    A, B = prior['large_pencil'], prior['small_pencil']
    assert set(A+B) == set(lines)
    assert all(lines[i] == exact(prior['exact_pool'][i]) for i in lines)
    P, Q = cross(lines[A[0]], lines[A[1]]), cross(lines[B[0]], lines[B[1]])
    assert all(incident(lines[i], P) for i in A)
    assert all(incident(lines[i], Q) for i in B)
    assert len({cross(lines[i],lines[j]) for i in A for j in B}) == 32
    connector = cross(P,Q)
    assert all(not incident(lines[i], Q) for i in A)
    assert all(not incident(lines[i], P) for i in B)

    deletion_histogram = Counter()
    for deleted in combinations(sorted(lines), 5):
        remaining = sorted(set(lines)-set(deleted))
        a = len(set(deleted).intersection(A))
        b = 5-a
        counts = [len(v) for v in line_vertices([lines[i] for i in remaining])]
        if (a,b) in [(2,3),(3,2)]:
            assert max(counts)>5
            deletion_histogram['over_target'] += 1
        elif (a,b)==(1,4):
            assert counts == [1]*7
            deletion_histogram['seven_line_pencil'] += 1
        else:
            assert (a,b) in [(4,1),(5,0)]
            deletion_histogram['survivor_core'] += 1

    records = []
    for kind in ['4,1','5,0']:
        cores = (((list(b),list(a)) for a in combinations(A,4) for b in combinations(B,3))
                 if kind == '4,1' else ((list(a),B) for a in combinations(A,3)))
        for S, U in cores:
            saturated, opposite = [lines[i] for i in S], [lines[i] for i in U]
            counts = [len(v) for v in line_vertices(saturated+opposite)]
            assert counts == [5]*3+[4]*4
            C = cross(saturated[0],saturated[1])
            candidates = transversals(saturated, opposite)
            rows = []
            for line, permutation in candidates.items():
                assert not incident(line,C) and line not in lines.values()
                omitted = next(j for j in range(4) if j not in permutation)
                extra = cross(line, opposite[omitted])
                old_vertices = line_vertices(saturated+opposite)[3+omitted]
                assert extra not in old_vertices
                rows.append(dict(line=strings(line), opposite_indices=[U[j] for j in permutation],
                                 omitted_line=U[omitted], new_vertex=strings(extra)))
            upper_bound = len(candidates)+2  # connector + at most one other C-line
            if upper_bound >= 5:
                raise RuntimeError('Counting certificate inconclusive; finite completion analysis needed')
            records.append(dict(type=kind, saturated=S, opposite=U, core_counts=counts,
                                transversals=rows, max_additions=upper_bound,
                                five_subsets_possible=0))
    assert len(records) == 336
    hist = {kind:dict(Counter(len(r['transversals']) for r in records if r['type']==kind))
            for kind in ['4,1','5,0']}
    return dict(complete=True, N=12, replacement_distance_lower_bound=6,
                scope='arbitrary complex projective added lines',
                source_sha256=sha256(source.read_bytes()).hexdigest(),
                previous_sha256=sha256(previous.read_bytes()).hexdigest(),
                original_indices=sorted(lines), exact_lines={i:strings(v) for i,v in lines.items()},
                centre_connector=strings(connector), deletion_histogram=dict(deletion_histogram),
                transversal_histograms=hist, injection_checks=336*24,
                cases=records, five_line_completions_tested=0, zeroes=[])


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = certify()
    with args.output.open('x') as out:
        json.dump(result,out,indent=2)
        out.write('\n')
    print('Deletion cases:',result['deletion_histogram'])
    print('Transversal counts:',result['transversal_histograms'])
    print('No five-line completions; replacement distance >= 6.')
