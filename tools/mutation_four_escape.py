#!/usr/bin/env python3
"""Exact, finite 4-for-4 enumeration for the saved two-pencil arrangement."""
import argparse
from collections import Counter
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp
from sympy import QQ_I as K

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from src.envs.mutations import build_candidate_pool, reflection_arrangement, _projective_key
from src.envs.projective_equivalence import classify_imprimitive_reflection_copy
from tools.mutation_near_analysis import rational


def norm(v):
    pivot = next((x for x in v if x), None)
    return tuple(x / pivot for x in v) if pivot else None


def cross(a, b):
    return norm((a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]))


def exact(v):
    return norm(tuple(K.from_sympy(sp.sympify(x)) for x in v))


def strings(v):
    return [str(K.to_sympy(x)) for x in v]


def pool_exact():
    # Rebuild the geometric pool over Q(i), independently of rounded keys.
    seed = []
    for a in range(4):
        w = sp.I**a
        seed.extend(exact(v) for v in [(1,-w,0), (0,1,-w), (-w,0,1)])
    vertices = list(dict.fromkeys(cross(a,b) for a,b in combinations(seed,2)))
    generated = set(seed) | {cross(a,b) for a,b in combinations(vertices,2)}
    generated.discard(None)
    # Floating-point pool supplies only the established candidate-index labels.
    floating = build_candidate_pool(reflection_arrangement(4, 'gmm'))
    pool = [exact([rational(z.real)+sp.I*rational(z.imag) for z in v]) for v in floating]
    assert len(pool) == len(set(pool)) == 63 and set(pool) == generated
    return pool


def lattice(indices, table):
    vertices = {}
    for i,j in combinations(indices,2):
        vertices.setdefault(table[i,j], set()).update((i,j))
    counts = [sum(i in s for s in vertices.values()) for i in indices]
    return counts, vertices


def exact_equivalence(pool, indices, label, numeric_matrix):
    """Use the numerical classifier's correspondence, then certify it exactly."""
    m = int(label.split('(')[1].split(',')[0])
    family = 'gmm' if label.split(',')[1] == str(m) else 'gm1'
    source = [] if family == 'gmm' else [sp.eye(3)[:,i] for i in range(3)]
    for a in range(m):
        w = sp.expand_complex(sp.exp(2*sp.pi*sp.I*a/m))
        source.extend(sp.Matrix(v) for v in [(1,-w,0),(0,1,-w),(-w,0,1)])
    target = [sp.Matrix([K.to_sympy(z) for z in pool[i]]) for i in indices]
    keys = {_projective_key(np.array(v,dtype=complex).ravel()): j for j,v in enumerate(target)}
    mapping = [keys[_projective_key(numeric_matrix @ np.array(v,dtype=complex).ravel())] for v in source]
    for frame in combinations(range(len(source)),4):
        S = sp.Matrix.hstack(*(source[i] for i in frame[:3]))
        if S.det() == 0:
            continue
        a = S.inv()*source[frame[3]]
        if any(sp.simplify(x)==0 for x in a):
            continue
        T = sp.Matrix.hstack(*(target[mapping[i]] for i in frame[:3]))
        b = T.inv()*target[mapping[frame[3]]]
        B = sp.simplify(T*sp.diag(*(b[i]/a[i] for i in range(3)))*S.inv())
        assert sp.simplify(B.det()) != 0
        for i,v in enumerate(source):
            assert all(sp.simplify(x)==0 for x in (B*v).cross(target[mapping[i]]))
        return dict(matrix=[[str(x) for x in B.row(i)] for i in range(3)],
                    source_to_candidate=[indices[j] for j in mapping])
    raise AssertionError('No projective frame')


def enumerate_escape(source, output):
    saved = json.loads(source.read_text())
    pool = pool_exact()
    for i,v in zip(saved['indices'],saved['exact_lines']):
        assert pool[i] == exact(v)
    table = {(i,j): cross(pool[i],pool[j]) for i in range(63) for j in range(63) if i!=j}
    original = set(saved['indices'])
    A = next(v['lines'] for v in saved['incidence'] if len(v['lines'])==8)
    B = next(v['lines'] for v in saved['incidence'] if len(v['lines'])==4)
    P,Q = table[A[0],A[1]],table[B[0],B[1]]
    connector = pool.index(cross(P,Q))
    outside = sorted(set(range(63))-original)
    deletion_audit = Counter()
    for r in range(5):
        for deleted in combinations(sorted(original),r):
            remain = sorted(original-set(deleted))
            counts,_ = lattice(remain,table)
            if max(counts)>5:
                deletion_audit[f'{r}:already_over_target'] += 1
            else:
                deletion_audit[f'{r}:not_excluded_by_survivors'] += 1
                assert r==4 and (set(deleted)<=set(A) or set(deleted)==set(B))
    off_centre = {i: len({table[i,j] for j in A}) for i in outside}
    assert set(off_centre.values()) <= {1,8}
    # Any off-centre addition has >=8 vertices. All-centre additions have v=1.
    result = dict(complete=False, source=str(source), pool_size=63,
                  exact_pool=[strings(v) for v in pool], large_pencil=A, small_pencil=B,
                  connector=connector, deletion_audit=dict(deletion_audit),
                  small_pencil_deletion_intersection_counts=off_centre,
                  cases=[], zeroes=[], quartet_defects={}, connector_defects={})
    output.mkdir(parents=True,exist_ok=False)
    summary = output/'enumeration.json'
    def save():
        summary.write_text(json.dumps(result,indent=2)+'\n')
    hist, connector_hist = Counter(), Counter()
    for survivors in combinations(A,4):
        core = sorted(list(survivors)+B)
        counts,vertices = lattice(core,table)
        assert counts == [5]*8
        grid = {table[i,j] for i in survivors for j in B}
        assert len(grid)==16 and P not in grid and Q not in grid
        admissible = [i for i in outside if all(table[i,j] in grid for j in core)]
        # Independent exhaustive construction from the 24 row-column bijections.
        geometric_transversals = set()
        for perm in permutations(B):
            points = [table[i,j] for i,j in zip(survivors,perm)]
            line = cross(points[0],points[1])
            if all(sum(x*y for x,y in zip(line,p)) == K.zero for p in points):
                geometric_transversals.add(line)
        assert {pool[i] for i in admissible} == geometric_transversals.intersection(pool)
        # Independent row/column-permutation check of each transversal.
        transversals = []
        for i in admissible:
            permutation = [next(c for c,j in enumerate(B) if table[i,r]==table[r,j]) for r in survivors]
            assert sorted(permutation)==list(range(4))
            transversals.append(dict(index=i,columns=permutation))
        preserving = [i for i in outside if all(table[i,j] in vertices for j in core)]
        assert set(preserving)==set(admissible)|{connector}
        case = dict(survivors=list(survivors),transversals=transversals,
                    geometric_transversals=[strings(v) for v in sorted(geometric_transversals,key=str)],
                    grid=[[strings(table[i,j]) for j in B] for i in survivors],
                    quartets=0, connector_quartets=0)
        result['cases'].append(case)
        quartets = [(q,False) for q in combinations(admissible,4)]
        quartets += [(tuple(sorted((connector,)+q)),True) for q in combinations(admissible,3)]
        for added,has_connector in quartets:
            case['connector_quartets' if has_connector else 'quartets'] += 1
            indices = sorted(core+list(added))
            counts,incidence = lattice(indices,table)
            D = sum((3*v-15)**2 for v in counts)
            (connector_hist if has_connector else hist)[D] += 1
            if D:
                continue
            pairpoints = {table[i,j] for i,j in combinations(added,2)}
            pairings = []
            for j in added[1:]:
                left = (added[0],j)
                right = tuple(i for i in added if i not in left)
                if all(table[a,b] in grid for a in left for b in right):
                    pairings.append([list(left),list(right)])
            pattern = 'concurrent_3_net' if len(pairpoints)==1 else ('two_pairs' if pairings else 'other')
            if pattern=='concurrent_3_net':
                assert next(iter(pairpoints)) not in grid
                assert all(sum(table[i,r]==g for i in added for r in survivors)==1 for g in grid)
            zero = dict(indices=indices,added=list(added),survivors=list(survivors),
                        N=12,D=0,bad_lines=0,counts=counts,pattern=pattern,pairings=pairings,
                        exact_lines=[strings(pool[i]) for i in indices],
                        multiplicities=dict(Counter(map(len,incidence.values()))),
                        incidence=[dict(point=strings(p),lines=sorted(ls)) for p,ls in incidence.items()],
                        classification='PENDING')
            result['zeroes'].append(zero)
            save()  # Save full exact evidence before invoking the classifier.
            label,matrix = classify_imprimitive_reflection_copy([np.array([complex(K.to_sympy(x)) for x in pool[i]]) for i in indices])
            zero['classification'] = label or 'UNCLASSIFIED'
            if label:
                zero['exact_equivalence'] = exact_equivalence(pool,indices,label,matrix)
            save()
            print('ZERO',label or 'UNCLASSIFIED',indices,pattern,flush=True)
            if label is None:
                print('STOP: unclassified zero saved; enumeration incomplete',flush=True)
                return result
    result.update(complete=True,quartet_defects=dict(hist),connector_defects=dict(connector_hist),
                  transversal_count_histogram=dict(Counter(len(c['transversals']) for c in result['cases'])),
                  zero_patterns=dict(Counter(z['pattern'] for z in result['zeroes'])),
                  zero_classifications=dict(Counter(z['classification'] for z in result['zeroes'])))
    save()
    print('complete',len(result['cases']),'cases; transversal counts',dict(Counter(len(c['transversals']) for c in result['cases'])),
          'quartet defects',dict(hist),'connector defects',dict(connector_hist),flush=True)
    return result


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--source',type=Path,default=Path('experiments/mutations-orbit/two-pencil-analysis.json'))
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    enumerate_escape(args.source,args.output)
