import unittest
import numpy as np
from fd_portal_entropy_short_tangent_windows import gain_difference

class ShortTangent(unittest.TestCase):
    def test_resolves_small_gain_defects(self):
        v=np.array([-2e-6,3e-6])
        np.testing.assert_allclose(gain_difference(np.eye(2),np.diag((1+v)**2)),v,atol=3e-16)
    def test_fixed_change_of_plane_coordinates_does_not_change_gains(self):
        A=np.diag([1.,3.]); B=np.diag([.9,3.2]); Q=np.array([[2.,1.],[1.,3.]])
        np.testing.assert_allclose(gain_difference(Q.T@A@Q,Q.T@B@Q),gain_difference(A,B),atol=2e-15)
