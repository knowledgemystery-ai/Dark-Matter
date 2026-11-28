import unittest

import numpy as np

from fd_portal_cubic_rank_factor_attribution import cross_thresholds


class CubicRankFactorAttribution(unittest.TestCase):
    def test_proportional_rank_changes_distance_from_three_not_whole_slope(self):
        v=cross_thresholds(4.,5.,.8)
        np.testing.assert_allclose(v,[[4.,4.25],[4.6,5.]],rtol=2e-15)

    def test_rank_factor_must_be_finite_and_positive(self):
        for r in (0.,-1.,float('inf')):
            with self.assertRaises(ArithmeticError):
                cross_thresholds(4.,5.,r)
