import json
from pathlib import Path
import unittest

from early_pole_native_d4_exact_partition import DATE,run


class D4ExactFiveBands(unittest.TestCase):
    def test_full_budget_and_unique_panel_partition(self):
        r=Path(__file__).resolve().parent/'results'
        Q=run(r)
        self.assertEqual(Q,json.loads((r/f'early-pole-native-d4-exact-partition-{DATE}.json').read_bytes()))
        self.assertEqual(sum(Q['counts']),111)
        self.assertEqual([z['sign'] for z in Q['bands']],[1,1,1,-1,-1])
        for m in Q['by_mesh']:
            self.assertLess(m['held_net_partition_closure'],5e-33)
            self.assertAlmostEqual(m['central_share_of_positive_total'],.83917413925,places=10)
            self.assertLess(m['opposing_upper_fraction_of_positive_total'],6.2e-5)
