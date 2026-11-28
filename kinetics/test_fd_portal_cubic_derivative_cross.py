import unittest

import numpy as np

from fd_portal_cubic_derivative_cross import two_by_two


class CubicDerivativeCross(unittest.TestCase):
    def test_cross_without_replacing_base_form(self):
        G=np.eye(2)
        R=np.diag([0.,1.])
        p=(G,np.diag([-2.,1.]),R,5.)
        c=(G,np.diag([-2.,-.5]),R,4.25)
        a,b=two_by_two((p,c))
        self.assertAlmostEqual(a['critical_heat_capacity_slope'],4.)
        self.assertAlmostEqual(b['critical_heat_capacity_slope'],4.75)
        np.testing.assert_allclose(a['rates_with_pchip_and_cubic_slope'],[1.,.25])
        np.testing.assert_allclose(b['rates_with_pchip_and_cubic_slope'],[.25,-.5])
