import unittest

import numpy as np

from fd_portal_entropy_hot_start_finite_crossing import locate, observable


class HotStartCrossing(unittest.TestCase):
    def test_same_plane_coordinates_remain_in_the_observable(self):
        Q=np.array([[2.,1.],[1.,3.]])
        def form(t):
            return np.diag([1.+.1*t,1.-.2*t])
        f=observable(lambda _:Q,form,0.)
        self.assertAlmostEqual(f(.1),np.sqrt(1.01)-1,places=14)

    def test_bracket_is_signed_and_directional(self):
        self.assertAlmostEqual(locate(lambda x:.3-x,.1,.6),.3)
        with self.assertRaises(ArithmeticError):
            locate(lambda x:x-.3,.1,.6)
