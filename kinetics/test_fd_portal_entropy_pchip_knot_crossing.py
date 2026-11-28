import unittest
from types import SimpleNamespace

import numpy as np
from scipy.interpolate import CubicSpline, PchipInterpolator

from fd_portal_entropy_pchip_knot_crossing import slopes


class TabulatedHeatCapacitySlopes(unittest.TestCase):
    def test_cubic_is_continuous_but_pchip_need_not_be(self):
        x=np.array([0.,.3,.7,1.])
        g=np.array([8.,9.,12.,16.])
        for interp,expected in ((CubicSpline,True),(PchipInterpolator,False)):
            s=interp(x,g)
            e=SimpleNamespace(gs=s,slope=s.derivative())
            left,right=slopes(e,.3)
            self.assertEqual(abs(left-right)<1e-12,expected)
            def reference(t):
                return 3+(3*s(t,1)+s(t,2))/(3*s(t)+s(t,1))
            self.assertAlmostEqual(left,reference(.3-1e-9),places=6)
            self.assertAlmostEqual(right,reference(.3+1e-9),places=6)

    def test_rejects_not_a_knot(self):
        s=PchipInterpolator([0.,.5,1.],[1.,2.,4.])
        e=SimpleNamespace(gs=s,slope=s.derivative())
        with self.assertRaises(ArithmeticError):
            slopes(e,.4)


if __name__=='__main__':
    unittest.main()
