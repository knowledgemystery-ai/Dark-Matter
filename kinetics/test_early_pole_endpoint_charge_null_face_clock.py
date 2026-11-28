import unittest
import numpy as np
from early_pole_endpoint_charge_nested_moment_remap import local_q
from early_pole_endpoint_charge_null_face_clock import adjust


class ChargeNullFaceControl(unittest.TestCase):
    def test_motion_preserves_both_moments_and_all_caps(self):
        q=np.array([.001,.02,.02,.04])
        w=np.array([.15,.23,.31,.31])
        f=np.array([.39,.31,.31,.17])
        m=.04
        e=local_q(q,q[0],q[-1],m)
        for mu in (-.05,.05):
            z,h=adjust(q,w,f,m,q[[0,-1]],mu)
            self.assertGreaterEqual(min(z),0.)
            self.assertLessEqual(max(z),1.)
            np.testing.assert_allclose([w@z,(w*e)@z],[w@f,(w*e)@f],rtol=2e-15)
            self.assertLess(h,5e-15)
            self.assertGreater(mu*(z[1]-f[1]),0.)
            self.assertAlmostEqual(z[1],z[2],places=14)


if __name__=='__main__':
    unittest.main()
