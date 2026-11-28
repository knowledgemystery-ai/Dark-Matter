import unittest

import numpy as np

from early_pole_endpoint_scalar_gauss2_energy_reconstruction import scaled_energy
from early_pole_endpoint_warm_streaming_pseudopressure_constraint import cell,feature,kernel


class PseudoPressureStreaming(unittest.TestCase):
    def test_regularized_fourth_feature_matches_rational_identity(self):
        E,h=.38,.16
        r=h/E
        t=np.array([0.,.11,.34,.77,1.])
        z=1+r*t
        rhs=(6/z-5+3*r*t-1/z**3)/r**3
        np.testing.assert_allclose(feature(t,E,h)[3],rhs,atol=2e-14,rtol=3e-14)

    def test_exact_null_direction_preserves_both_pressure_moments(self):
        a,b,m,A,p=.011,.022,.00087,100.,38.
        q=np.array([.012,.0139,.016,.0198,.021])
        w=np.array([.16,.16,.22,.27,.19])
        E=np.hypot(a,m)
        H=(b-a)*(b+a)/(E+np.hypot(b,m))
        Z=feature(scaled_energy(q,a,b,m),E,H)
        _,_,V=np.linalg.svd(Z)
        e=np.hypot(q,m)
        outer,inner,ref,error,atoms=cell(a,b,m,A,p,q,w,128)
        self.assertLess(error,2e-13)
        self.assertTrue(all(len(v['q'])<=4 for v in atoms))
        self.assertLessEqual(inner[0],ref)
        self.assertGreaterEqual(inner[1],ref)
        for h in (w-.08*V[-1],w+.08*V[-1]):
            self.assertGreater(np.min(h),0.)
            np.testing.assert_allclose(np.array([np.ones_like(q),e,q*q/e,q**4/e**3])@h,
                np.array([np.ones_like(q),e,q*q/e,q**4/e**3])@w,rtol=2e-14,atol=2e-16)
            v=float(h@kernel(q,m,A,p))
            self.assertLessEqual(outer[0],v)
            self.assertGreaterEqual(outer[1],v)


if __name__=='__main__':
    unittest.main()
