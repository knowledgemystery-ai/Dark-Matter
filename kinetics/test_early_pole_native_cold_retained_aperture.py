import unittest
from pathlib import Path

from early_pole_native_cold_retained_aperture import p, run, sep


class NativeColdRetainedAperture(unittest.TestCase):
    def test_stable_separation_agrees_with_direct_momentum_difference(self):
        M,m,ep,t=.05,.003,26.,-.0014
        ed=ep-t
        self.assertAlmostEqual(sep(t,M,m,ep,ed),p(ed,m)-p(ep,M),places=13)

    def test_only_d4_inverse_has_a_clipped_cosine_interval(self):
        p=Path(__file__).resolve().parent/'results'
        out=run(p,p/'early-pole-native-d2d4-cold-cross-q41-2026-11-28.json')['rows']
        self.assertTrue(all(r['first_clipped_q'] is None for r in out[:-1]))
        self.assertTrue(1.2e-5<out[-1]['first_clipped_q']<1.3e-5)
        self.assertTrue(.92<out[-1]['samples'][-1]['geometric_fraction_of_cosine_interval']<.93)
