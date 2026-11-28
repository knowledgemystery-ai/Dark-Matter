import unittest

import numpy as np

from early_pole_endpoint_charge_cold_count_gap_ladder import residual, split_terms


class FixedScalarGap(unittest.TestCase):
    def test_old_shape_rescaling_and_warm_ledger(self):
        q=np.array([.1,.2,1.1])
        p=np.array([1.,2.,4.])
        r=np.array([1.3,2.6,3.5])
        m=.4
        a,b,c,d=split_terms(q,p,r,m,np.array([True,True,False]))
        self.assertAlmostEqual(b,0.)
        self.assertAlmostEqual(c,-.5/np.hypot(1.1,m))
        self.assertAlmostEqual(a+b+c,d)

    def test_same_cold_count_different_shape(self):
        q=np.array([.1,.2,1.1])
        p=np.array([1.,2.,4.])
        r=np.array([2.,1.,4.])
        a,b,c,d=split_terms(q,p,r,.4,np.array([True,True,False]))
        self.assertEqual(a,0.)
        self.assertEqual(c,0.)
        self.assertGreater(b,0.)
        self.assertAlmostEqual(b,d)

    def test_lower_monotone_derivative(self):
        m,u=.7,.8
        q=np.array([.1,.3,1.5])
        p=np.array([1.,.5,.1])
        v=(residual(.03,q,p,u)-residual(.03,q,p,m))/(u-m)
        self.assertGreater(v,2*m)


if __name__=='__main__':
    unittest.main()
