import unittest

import numpy as np

from fd_portal_heat_capacity_slope_threshold import threshold


class HeatCapacityThreshold(unittest.TestCase):
    def test_rank_one_critical_slope_and_reoptimized_eigenvector(self):
        k=3.5
        S0=np.diag([-1.,-3.])
        R=np.ones((2,2))
        G=np.diag([2.,1.])
        r=threshold(G,S0+(k-3)*R,R,k)
        self.assertAlmostEqual(r['critical_slope_above_three'],3.75,places=12)
        self.assertAlmostEqual(r['independent_bisection_slope'],3.75,places=12)
        self.assertLess(max(r['rate_at_archived_slope']),0)
        np.testing.assert_allclose(r['rate_at_cubic_slope'],[-3.,-.5],rtol=2e-15)

    def test_nonnegative_reference_does_not_receive_spurious_threshold(self):
        V=np.eye(2)
        r=threshold(V,np.diag([1.,-1.]),np.diag([1.,0.]),3.)
        self.assertIsNone(r['critical_slope_above_three'])
