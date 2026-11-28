import unittest

import numpy as np

from fd_portal_entropy_plane_rotation_threshold import cone_distance, plane_angle


class PlaneRotation(unittest.TestCase):
    def test_known_tilted_conical_boundary(self):
        t=np.deg2rad(10)
        k=np.diag([-1.,-1.,1.])
        n=np.array([np.sin(t),0.,np.cos(t)])
        l,h=cone_distance(k,n,16384)['degrees_interval']
        self.assertLessEqual(l,35.)
        self.assertGreaterEqual(h,35.)
        self.assertLess(h-l,.02)
        self.assertAlmostEqual(plane_angle(np.eye(3)[:,:2],n),10.,places=12)

    def test_reference_plane_already_growing_is_rejected(self):
        with self.assertRaises(ArithmeticError):
            cone_distance(np.diag([-1.,-1.,1.]),np.eye(3)[0],128)
