import json
import unittest
from pathlib import Path
import numpy as np
from fd_portal_entropy_short_nonlinear_witness import DATE, endpoint_columns

class NonlinearWitness(unittest.TestCase):
    def test_symmetric_trajectories_recover_fractional_columns(self):
        Q=np.eye(3)[:,:2]
        ctrl=np.array([[1.,2.,3.],[1.5,2.5,3.5]])
        r=endpoint_columns(lambda v:ctrl*v,Q,1e-5,ctrl)
        np.testing.assert_allclose(r,np.broadcast_to(Q,r.shape),atol=1e-10)
    def test_completed_two_settings_agree_including_knot_crossing(self):
        p=Path(__file__).resolve().parent/'results'/f'fd-portal-entropy-short-nonlinear-witness-{DATE}.json'
        r=json.loads(p.read_bytes())
        x=np.asarray(r['ordered_and_finite_gain_minus_one'])[:,:,-1]
        np.testing.assert_equal(np.sign(x),np.tile([1,1,-1,1],(3,1)))
        self.assertLess(r['maximum_tangent_comparison_defect'],3e-11)
