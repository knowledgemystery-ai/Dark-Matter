import unittest

import numpy as np

from early_pole_endpoint_scalar_gauss2_energy_reconstruction import scaled_energy
from early_pole_endpoint_warm_streaming_pressure_constraint import cell,feature,kernel


class PressureConstrainedStreaming(unittest.TestCase):
    def test_stable_feature_is_exact_pressure_constraint(self):
        a,b,m=.011,.023,.00087
        E=np.hypot(a,m)
        H=(b-a)*(b+a)/(E+np.hypot(b,m))
        t=np.array([0.,.11,.34,.89,1.])
        L=1./(E+H*t)-((1.-t)/E+t/(E+H))
        expect=-L*E*E*(E+H)/H**2
        np.testing.assert_allclose(feature(t,E,H)[2],expect,rtol=2e-14,atol=4e-16)

    def test_positive_moment_preserving_perturbations_inside_dual_enclosure(self):
        a,b,m,A,p=.016,.031,.00087,100.,24.
        q=np.array([.017,.019,.022,.03])
        w=np.array([.2,.1,.4,.3])
        E,H=np.hypot(a,m),np.hypot(b,m)-np.hypot(a,m)
        Z=feature(scaled_energy(q,a,b,m),E,H)
        _,_,V=np.linalg.svd(Z)
        step=.08*V[-1]
        moved=[w+step,w-step]
        self.assertGreater(np.min(moved),0.)
        for h in moved:
            np.testing.assert_allclose(Z@h,Z@w,atol=4e-16,rtol=1e-14)
            e=np.hypot(q,m)
            self.assertAlmostEqual(float(h@(q*q/e)),float(w@(q*q/e)),places=15)
        outer,inner,ref,defect,atoms=cell(a,b,m,A,p,q,w,96)
        self.assertLess(defect,2e-13)
        self.assertTrue(all(len(v["q"])<=3 for v in atoms))
        self.assertLessEqual(inner[0],ref)
        self.assertGreaterEqual(inner[1],ref)
        for witness in (h@kernel(q,m,A,p) for h in moved):
            self.assertLessEqual(outer[0],witness)
            self.assertGreaterEqual(outer[1],witness)


if __name__=='__main__':
    unittest.main()
