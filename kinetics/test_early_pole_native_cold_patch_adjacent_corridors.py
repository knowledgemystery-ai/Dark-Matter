import json
from pathlib import Path
import unittest

from early_pole_native_cold_patch_adjacent_corridors import DATE,run


class PatchAdjacentCorridors(unittest.TestCase):
    def test_old_negative_sign_is_now_resolved_on_adjacent_exact_span(self):
        root = Path(__file__).resolve().parent/'results'
        Z = run(root)
        self.assertEqual(Z,json.loads((root/f'early-pole-native-cold-patch-adjacent-corridors-{DATE}.json').read_bytes()))
        v = Z['rows'][4]
        self.assertEqual(v['q'],[8e-6,1e-5])
        self.assertTrue(v['D4_both_edges_exact_after_patch'])
        self.assertEqual(v['held_signs_D2_D4'],[-1,-1])
        n=Z['original_diagnostic_occupation']
        self.assertTrue(v['D4_prior_neutral_hull'][0]<n<v['D4_prior_neutral_hull'][1])
        self.assertLess(v['held_nets_D2_and_D4'][1][1],-5e-22)

    def test_documented_corridor_has_two_distinct_limiting_bands(self):
        root = Path(__file__).resolve().parent/'results'
        Z = json.loads((root/f'early-pole-native-cold-patch-adjacent-corridors-{DATE}.json').read_bytes())
        self.assertEqual(Z['limiting_band_q'],[[5e-6,8e-6],[3e-6,5e-6]])
        self.assertEqual(Z['limiting_levels'],[4,2])
        self.assertGreater(min(Z['adjoining_D4_neutral_gap_lower_above_q3']),5e-11)
        n = Z['original_diagnostic_occupation']
        self.assertTrue(Z['common_strict_occupation_corridor'][0]<n<Z['common_strict_occupation_corridor'][1])
