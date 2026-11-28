import json
from pathlib import Path
import unittest

from early_pole_native_cold_exact_paired_bands import DATE,run


class ExactMatchedPartition(unittest.TestCase):
    def test_complete_archives_replay_and_resolve_common_held_corridor(self):
        p=Path(__file__).resolve().parent/'results'
        X=run(p)
        self.assertEqual(X,json.loads((p/f'early-pole-native-cold-exact-paired-bands-{DATE}.json').read_bytes()))
        self.assertEqual(X['limiting_charts'],[4,4])
        self.assertEqual(X['limiting_spans'],[[5e-6,8e-6],[8e-6,1e-5]])
        self.assertGreater(X['common_open_corridor'][1],2.64e-9)
        self.assertTrue(all(r['D2_minus_D4_neutral_lower']>6e-10 for r in X['rows']))
