import unittest

import numpy as np

from fd_portal_cubic_remainder_attribution import swap_game


class CubicRemainderAttribution(unittest.TestCase):
    def test_all_orders_reproduce_an_additive_threshold_shift(self):
        a=np.array([np.diag([-.5,-.5]),np.diag([-.5,-.5])])
        b=a+np.array([np.diag([0.,-.25]),np.diag([-.2,-.05])])
        root,terms,single=swap_game(a,b,np.diag([0.,1.]))
        np.testing.assert_allclose(root,[4.,4.3])
        np.testing.assert_allclose(terms,[.25,.05])
        np.testing.assert_allclose(single,terms)

    def test_unsafe_mixed_form_is_rejected(self):
        a=[-np.eye(2),-np.eye(2)]
        b=[4*np.eye(2),-6*np.eye(2)]
        with self.assertRaises(ArithmeticError):
            swap_game(a,b,np.diag([0.,1.]))
