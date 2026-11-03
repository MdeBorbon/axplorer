"""Regression checks for evaluated-zero auditing and immediate stopping."""
import unittest
from unittest.mock import patch

from tools.mutation_explore_probe import AuditedExploreDataPoint as C, UnclassifiedZero


class MutationAuditTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        C.configure(4, 'gmm', 6, 12, 0, 4, .15, 20, .25)

    def setUp(self):
        self.records = []
        C.journal = staticmethod(self.records.append)
        C.seen_zeroes = set()
        self.d = C(C.POOL_SIZE)
        self.d.data[[3, 4, 5, 7, 9, 11, 13, 19, 24, 33, 39, 50, 53, 58, 61]] = 1

    def test_disguised_orbit_penalty_and_audit(self):
        self.assertFalse(self.d.is_tabu())
        self.assertEqual(self.d.defect(), 0)
        energy = self.d._local_energy()
        self.assertAlmostEqual(energy - self.d._exploration_energy_without_tabu(), 10)
        self.d._local_energy()
        self.assertEqual(len(self.records), 1)
        r = self.records[0]
        self.assertEqual(r['classification'], 'G(4,1,3)')
        self.assertEqual(r['multiplicities'], {2: 12, 3: 16, 6: 3})
        self.assertEqual(len(r['lines']), 15)
        self.assertEqual(r['errors'], [0]*15)

    def test_uncertain_classification_saves_before_stopping(self):
        with patch.object(self.d, 'reflection_orbit_label', return_value=None):
            with self.assertRaises(UnclassifiedZero):
                self.d._local_energy()
        self.assertEqual(self.records[-1]['classification'], 'UNCLASSIFIED')
        self.assertEqual(self.records[-1]['indices'], self.d.selected_indices())
        self.assertEqual(len(self.records[-1]['lines']), 15)


if __name__ == '__main__':
    unittest.main()
