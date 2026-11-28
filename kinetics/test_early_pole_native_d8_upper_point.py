import json
import unittest
from pathlib import Path

import numpy as np

from early_pole_native_d8_upper_point import DATE, select


class D8UpperPoint(unittest.TestCase):
    def test_disjoint_median_enclosures_are_not_used(self):
        with self.assertRaises(ArithmeticError):
            select(dict(fine_mesh=dict(charts=[dict(raw_median_possible=[1,2]),dict(raw_median_possible=[3,4])])))

    def test_archived_upper_point_has_complete_independent_mesh_controls(self):
        p=Path(__file__).resolve().parent/'results'
        load=lambda k:json.loads((p/f'early-pole-native-{k}-{DATE}.json').read_bytes())
        q,w=select(load('cold-weight-target'))
        r=load('d8-upper-point')
        self.assertEqual((q,w),(r['q'],r['intersection_of_raw_median_enclosures']))
        self.assertEqual([t['native_level'] for t in r['results']],[2,4,8])
        for t in r['results']:
            self.assertEqual(len(t['meshes']),2)
            for z in t['meshes']:
                np.testing.assert_allclose(z['stock'],z['real_pauli_mirror'],rtol=5e-15,atol=0.)
        self.assertTrue(r['results'][2]['inverse_geometric_cosine_fraction']<.83)
        self.assertLess(r['d8_over_d4'][0],.3)
