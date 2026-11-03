"""Independent checks of the finite enumeration and exact equivalence witness."""
import contextlib
from collections import Counter
import io
import json
from pathlib import Path
import tempfile
import unittest

import numpy as np
from sympy import QQ_I

from tools.mutation_four_escape import enumerate_escape, exact_equivalence, pool_exact
from src.envs.projective_equivalence import classify_imprimitive_reflection_copy


class FourEscapeTests(unittest.TestCase):
    def test_exhaustive_certificate(self):
        with tempfile.TemporaryDirectory() as tmp, contextlib.redirect_stdout(io.StringIO()):
            output = Path(tmp)/'result'
            r = enumerate_escape(Path('experiments/mutations-orbit/two-pencil-analysis.json'), output)
            self.assertEqual(json.loads(json.dumps(r)),
                             json.loads((output/'enumeration.json').read_text()))
        self.assertTrue(r['complete'])
        self.assertEqual(len(r['cases']), 70)
        self.assertEqual(r['transversal_count_histogram'], {0:68, 1:2})
        self.assertEqual(dict(Counter(len(c['geometric_transversals']) for c in r['cases'])),
                         {0:68, 1:2})
        self.assertEqual(r['zeroes'], [])
        self.assertEqual([(c['survivors'],c['transversals'][0]['index'])
                          for c in r['cases'] if c['transversals']],
                         [([23,29,39,56],51), ([23,44,55,56],19)])
        self.assertEqual(r['deletion_audit']['4:already_over_target'], 424)

    def test_known_zero_has_exact_projective_witness(self):
        pool = pool_exact()
        indices = list(range(12))
        lines = [np.array([complex(QQ_I.to_sympy(x)) for x in pool[i]]) for i in indices]
        label, B = classify_imprimitive_reflection_copy(lines)
        self.assertEqual(label, 'G(4,4,3)')
        witness = exact_equivalence(pool, indices, label, B)
        self.assertEqual(sorted(witness['source_to_candidate']), indices)


if __name__ == '__main__':
    unittest.main()
