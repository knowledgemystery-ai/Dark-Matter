import unittest
from pathlib import Path
import numpy as np
from early_pole_native_cold_shape_attribution import moment_interval,run

class ShapeAttribution(unittest.TestCase):
    def test_signed_measure_reverses_negative_endpoints(self):
        np.testing.assert_allclose(moment_interval(np.array([[1.,2.],[2.,3.]]),np.array([1.,-1.])),[-2.,0.])
    def test_exact_shared_shape_decomposition(self):
        p=Path(__file__).resolve().parent/'results'
        out=run(p,p/'early-pole-native-cold-bose-reweight-2026-11-28.json')['by_mesh'][1]['forward_number']
        self.assertTrue(out['all_bin_components_same_sign'])
        self.assertAlmostEqual(out['native_TV'],out['kernel_component_TV']+out['occupation_component_TV'],places=13)
        self.assertGreater(out['occupation_mean_shift_interval'][0],0)
        self.assertLess(out['native_mean_shift_interval'][0],0)
