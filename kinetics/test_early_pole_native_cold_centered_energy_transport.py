import unittest
from pathlib import Path
import numpy as np
from early_pole_native_cold_centered_energy_transport import integral, run

class MomentTransport(unittest.TestCase):
    def test_monotone_enclosures_are_integrated_in_intervals(self):
        q=np.array([0.,1.,2.])
        a=np.array([[0.,.8,1.],[0.,.8,1.]])
        b=np.array([[0.,.2,1.],[0.,.2,1.]])
        np.testing.assert_allclose(integral(q,a,b),[0.,1.6])

    def test_exact_centered_forward_bounds_improve_support_only(self):
        p=Path(__file__).resolve().parent/'results'
        r=run(p,p/'early-pole-native-cold-mean-from-centered-energy-2026-11-28.json')
        for x in r['by_mesh']:
            a=x['forward']['weighted']
            lo,hi=a['previous_support_only_distance_bounds']
            l,h=a['distance_with_independent_mean_lower']
            self.assertGreater(l,15*lo)
            self.assertLess(h,.8*hi)
