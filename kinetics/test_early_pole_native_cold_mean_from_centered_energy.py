import unittest
from pathlib import Path
import numpy as np
from early_pole_native_cold_mean_from_centered_energy import conditional_q,run

class EnergyMean(unittest.TestCase):
    def test_chord_and_jensen_bound_nonnegative_endpoint_mixture(self):
        q=np.array([[1.,2.]])
        u=np.sqrt(q*q+4)-2
        lo,hi=conditional_q(q,np.array([1.]),u.mean(axis=1),2.)[0]
        self.assertAlmostEqual(lo,1.5,places=11)
        self.assertGreater(hi,1.5)
        with self.assertRaises(ArithmeticError):
            conditional_q(q,np.array([1.]),np.array([2.]),2.)

    def test_independent_centered_moments_resolve_both_signs(self):
        p=Path(__file__).resolve().parent/'results'
        r=run(p,p/'early-pole-native-cold-shape-attribution-2026-11-28.json')['by_mesh']
        for x in r:
            for z in x.values():
                for k in ('raw_native_mean_difference','weighted_native_mean_difference',
                          'occupation_only_mean_difference','kernel_only_mean_difference'):
                    self.assertGreater(z[k][0],0)
