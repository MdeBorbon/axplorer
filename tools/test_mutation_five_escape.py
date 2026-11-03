"""Independent incidence enumeration and continuous-pencil boundary checks."""
from collections import Counter
from itertools import product
import unittest

from tools.mutation_five_escape import certify, cross, exact, incident, line_vertices


class FiveEscapeTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.result = certify()

    def test_complete_counts_and_bound(self):
        r = self.result
        self.assertTrue(r['complete'])
        self.assertEqual(r['deletion_histogram'],
                         {'seven_line_pencil':8,'over_target':448,'survivor_core':336})
        self.assertEqual(r['transversal_histograms'],
                         {'4,1':{0:232,1:46,2:2},'5,0':{0:46,1:10}})
        self.assertEqual(Counter(c['type'] for c in r['cases']), {'4,1':280,'5,0':56})
        self.assertEqual(r['injection_checks'],8064)
        self.assertTrue(all(c['max_additions']<5 for c in r['cases']))
        self.assertEqual(r['replacement_distance_lower_bound'],6)
        self.assertEqual(r['zeroes'],[])

    def test_independent_two_saturated_line_construction(self):
        # No row-column injections: join every pair of distinct old vertices
        # on the first two saturated lines, then check all saturated incidences.
        lines = {i:exact(v) for i,v in self.result['exact_lines'].items()}
        for c in self.result['cases']:
            S = [lines[i] for i in c['saturated']]
            U = [lines[i] for i in c['opposite']]
            C = cross(S[0],S[1])
            old = line_vertices(S+U)[:3]
            found = set()
            for p,q in product(old[0],old[1]):
                line = cross(p,q)
                if line is None or line in U or incident(line,C):
                    continue
                if all(cross(line,s) in vertices for s,vertices in zip(S,old)):
                    found.add(line)
            self.assertEqual(found,{exact(t['line']) for t in c['transversals']})

    def test_two_extra_centre_lines_exceed_unsaturated_budget(self):
        S = [exact((0,1,-j)) for j in range(3)]
        U = [exact((1,0,-j)) for j in range(4)]
        extras = [exact((0,1,-10)),exact((0,1,-11))]
        self.assertEqual([len(v) for v in line_vertices(S+U)][3:],[4]*4)
        self.assertEqual([len(v) for v in line_vertices(S+U+extras)][3:7],[6]*4)


if __name__=='__main__':
    unittest.main()
