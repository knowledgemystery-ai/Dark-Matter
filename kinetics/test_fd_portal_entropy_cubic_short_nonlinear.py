import json
import unittest
from pathlib import Path

import numpy as np

from fd_portal_entropy_cubic_short_nonlinear import DATE

class CubicNonlinearWitness(unittest.TestCase):
    def test_archived_controls_are_negative_for_each_setting_and_window(self):
        f=Path(__file__).resolve().parent/'results'/f'fd-portal-entropy-cubic-short-nonlinear-{DATE}.json'
        r=json.loads(f.read_bytes())
        t=np.asarray(r['tangent_and_two_direct_gain_minus_one'])
        self.assertEqual(t.shape,(3,4,2))
        self.assertLess(np.max(t[:,:,-1]),0.)
        self.assertLess(r['max_nonlinear_vs_tangent_defect'],2e-8)
        for Q in map(np.asarray,r['local_rebased_planes']):
            np.testing.assert_allclose(Q.T@Q,np.eye(2),atol=2e-14)
