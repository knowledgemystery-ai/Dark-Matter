import unittest

import numpy as np

from fd_portal_entropy_eos_interpolation_check import min_piece


class SplineBounds(unittest.TestCase):
    def test_interior_minimum_not_lost(self):
        self.assertAlmostEqual(min_piece(np.array([0.,1.,-3.,4.]),2.),1.75)
        self.assertAlmostEqual(min_piece(np.array([0.,1.,-3.,4.]),1.),2.)

    def test_cubic_sign_change_inside_positive_endpoints(self):
        self.assertLess(min_piece(np.array([1.,0.,-3.,1.]),2.),0.)
        self.assertEqual(min_piece(np.array([0.,0.,0.,2.]),1.),2.)
        self.assertEqual(min_piece(np.array([0.,0.,0.,0.]),1.),0.)


if __name__=='__main__':
    unittest.main()
