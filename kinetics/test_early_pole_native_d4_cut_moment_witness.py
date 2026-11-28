import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from early_pole_native_d4_cut_moment_witness import DATE, bump_means, center, mixture, run


class D4MomentWitness(unittest.TestCase):
    def test_bumps_have_means_in_the_right_order(self):
        q = [.01, .02, .04]
        m = .25
        h = bump_means(q, m, 16)
        self.assertTrue(np.all(np.diff(h)>0.))
        self.assertLess(h[1],center(q[1],m))
        self.assertGreater(h[2],center(q[1],m))
        p = mixture(h,h[1],(0,3))
        self.assertAlmostEqual(np.dot(p,h),h[1])

    def test_archived_moments_permit_both_signs(self):
        r = Path(__file__).parent/'results'
        with TemporaryDirectory() as d:
            Z = run(r,Path(d))
        self.assertEqual(Z,json.loads((r/f'early-pole-native-d4-cut-moment-witness-{DATE}.json').read_bytes()))
        self.assertTrue(all(v['witnesses'][0]['net']<-3e-22<v['witnesses'][1]['net'] for v in Z['rows']))
        shifts = [r['nearby_aperture_sensitivity'] for r in Z['rows']]
        self.assertLess(abs(shifts[0]['lower_cut_at_zero']-shifts[1]['lower_cut_at_zero']),2e-17)
        for s in shifts:
            (q0,lo),(q,zero),(q1,hi) = s['lower_cut_and_conservative_net_lower']
            self.assertTrue(q0 < q < q1 and lo < 0. < hi)
            self.assertAlmostEqual(zero,0.,places=33)
