import json
import unittest
from pathlib import Path
import numpy as np
from fd_portal_entropy_cubic_short_tangent import DATE

class CubicShortTangent(unittest.TestCase):
    def test_completed_estimate_is_stable_below_one_on_all_four_windows(self):
        p=Path(__file__).resolve().parent/'results'/f'fd-portal-entropy-cubic-short-tangent-{DATE}.json'
        r=json.loads(p.read_bytes())
        c=np.asarray(r['initial_rebased_fractional_columns'])
        np.testing.assert_allclose(c.T@c,np.eye(2),rtol=0.,atol=1e-14)
        self.assertEqual(len(r['windows']),4)
        self.assertTrue(all(z['ordered_gain_minus_one'][-1]<0 and
                            z['two_setting_defect']<abs(z['ordered_gain_minus_one'][-1])/100 for z in r['windows']))
