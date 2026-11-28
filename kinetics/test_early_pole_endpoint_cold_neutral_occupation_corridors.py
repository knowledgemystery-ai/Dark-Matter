import unittest
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_neutral_occupation_corridors import (
    corridor, overlap, root_bounds, run, sign)


class HeldOccupationCorridors(unittest.TestCase):
    def test_partial_panel_positivity_and_ratio(self):
        e=[[0,1],[1,2],[2,3]]
        f=np.array([[6,1,5],[9,6,3],[8,3,5]])
        self.assertEqual(root_bounds(e,f,1,2),[2,2])
        self.assertEqual(root_bounds(e,f,.5,2.5),[6/13,10/3])

    def test_possible_is_not_sufficient_and_intervals_are_open(self):
        c=corridor([2,4],[1,3])
        self.assertEqual(c['d2_loss_d4_gain_possible_within_enclosures'],[1,4])
        self.assertIsNone(c['d2_loss_d4_gain_sufficient'])
        self.assertEqual(sign(3,[1,3]),0)
        self.assertEqual(sign(1,[1,3]),0)
        self.assertIsNone(overlap([1,2],[2,3]))

    def test_exact_admitted_meshes_disagree_in_full_and_center_but_never_upper(self):
        root=Path(__file__).resolve().parent/'results'
        a=run(root)
        self.assertIsNotNone(a['simultaneous_full_and_central_sufficient'])
        self.assertIsNone(a['simultaneous_full_and_upper_possible'])
        self.assertTrue(a['both_upper_negative_throughout_full_opposite_sign_corridor'])
        self.assertGreater(a['gap_from_full_corridor_ceiling_to_first_upper_neutral'],1e-11)
