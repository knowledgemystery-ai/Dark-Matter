import unittest
import numpy as np
from fd_portal_entropy_finite_rank_threshold import threshold

class FiniteRankThreshold(unittest.TestCase):
    def test_reference_loss_first_reaches_zero_at_generalized_threshold(self):
        B=-np.eye(2)
        R=np.diag([2.,1.])
        t=threshold(B,R)
        self.assertAlmostEqual(t,.5)
        self.assertAlmostEqual(np.linalg.eigvalsh(B+t*R)[-1],0.)
        self.assertLess(np.linalg.eigvalsh(B+(t-.01)*R)[-1],0)
    def test_negative_increments_do_not_get_a_monotone_threshold(self):
        with self.assertRaises(ArithmeticError):
            threshold(-np.eye(2),np.diag([2.,-1.]))
