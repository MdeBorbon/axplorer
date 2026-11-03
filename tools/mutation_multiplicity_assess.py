#!/usr/bin/env python3
"""Deterministic scoring diagnostic; does not alter the production objective."""
import argparse
from collections import Counter
from fractions import Fraction as F
from itertools import combinations
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from tools.mutation_four_escape import cross, exact, lattice


def assess():
    prior = json.loads(Path('experiments/mutations-four-escape-verified/enumeration.json').read_text())
    source = json.loads(Path('experiments/mutations-orbit/two-pencil-analysis.json').read_text())
    pool = [exact(v) for v in prior['exact_pool']]
    table = {(i,j):cross(pool[i],pool[j]) for i,j in combinations(range(63),2)}
    def stats(indices):
        indices = sorted(indices)
        counts,vertices = lattice(indices,table)
        profile = Counter(map(len,vertices.values()))
        errors = [v-5 for v in counts]
        bad = sum(e!=0 for e in errors)
        defect = 9*sum(e*e for e in errors)
        excess = sum(number*max(0,m-5)**2 for m,number in profile.items())
        novelty = len(set(indices)^set(range(12)))
        energy = F(bad,12)+F(defect,9*12*50)-F(novelty,10000*64)
        return dict(indices=indices, D=defect, bad_lines=bad, errors=errors,
                    multiplicities=dict(profile), excess=excess, energy=str(energy))
    original = set(source['indices'])
    start = stats(original)
    candidates = []
    for remove in sorted(original):
        for add in sorted(set(range(63))-original):
            r = stats(original-{remove}|{add})
            r.update(remove=remove,add=add)
            if r['excess'] < start['excess']:
                r['lambda_threshold'] = str((F(r['energy'])-F(start['energy'])) /
                                           F(start['excess']-r['excess'],25))
            candidates.append(r)
    thresholds = [r for r in candidates if 'lambda_threshold' in r]
    first = min(thresholds,key=lambda r:F(r['lambda_threshold']))
    sweep = {}
    for weight in [F(0),F(1,50),F(1,20),F(1,10),F(1,4),F(1),F(2)]:
        old = F(start['energy'])+weight*F(start['excess'],25)
        better = [r for r in candidates if F(r['energy'])+weight*F(r['excess'],25)<old]
        sweep[str(weight)] = dict(improving_swaps=len(better),
                                 bad_line_counts=dict(Counter(r['bad_lines'] for r in better)))
    return dict(N=12,neighborhood_size=len(candidates),baseline=start,
                diagnostic='existing_energy + lambda * multiplicity_excess / 25',
                first_improvement_threshold=first,
                lower_excess_without_worse_bad_fraction=sum(r['bad_lines']<=start['bad_lines'] for r in thresholds),
                weight_sweep=sweep,
                neighbors=[{k:r[k] for k in ('remove','add','D','bad_lines','excess','energy')}
                           for r in candidates])


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    result=assess()
    with args.output.open('x') as f:
        json.dump(result,f,indent=2)
        f.write('\n')
    print(json.dumps({k:v for k,v in result.items() if k!='neighbors'},indent=2))
