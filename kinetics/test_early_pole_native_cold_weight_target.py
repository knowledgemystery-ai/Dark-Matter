import unittest
from pathlib import Path

import numpy as np

from early_pole_native_cold_weight_target import mass_interval, median_interval, run


class ColdWeightTarget(unittest.TestCase):
    def test_support_bounds_do_not_assume_intrapanel_shape(self):
        B=np.array([[0.,1.],[1.,2.],[2.,3.]])
        f=np.array([.2,.6,.2])
        np.testing.assert_allclose(mass_interval(B,f,.5,2.5),[.6,1.])
        self.assertEqual(median_interval(B,f),[1.,2.])

    def test_exact_chart_ledgers_separate_raw_and_weighted_targets(self):
        p=Path(__file__).resolve().parent/'results'
        a=run(p,p/'early-pole-native-cold-bose-reweight-2026-11-28.json')['fine_mesh']
        for i in a['charts']:
            self.assertLess(i['central_raw_fraction'][1],.013)
            self.assertGreater(i['central_weighted_fraction'][0],.8)
            self.assertGreater(i['upper_raw_fraction'][0],.85)
            self.assertGreater(i['raw_median_possible'][0],i['weighted_median_possible'][1])
