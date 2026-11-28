import unittest
import numpy as np
from early_pole_endpoint_charge_nested_moment_remap import fit,local_q,transfer


class NestedChargeMomentRemap(unittest.TestCase):
    def test_empty_and_infeasible_capacity(self):
        w=np.ones(2)
        q=np.array([0.,1.])
        np.testing.assert_array_equal(fit(w,q,0.,0.),[0.,0.])
        with self.assertRaises(ArithmeticError):
            fit(w,q,1.8,.1)

    def test_stable_centered_energy_and_conservative_round_trip(self):
        c=(np.array([.0001,.004]),np.array([.0001,.004]),np.array([.4,.6]))
        f=(np.array([.0001,.002,.004]),np.array([.0001,.002,.002,.004]),np.array([.09,.31,.23,.37]))
        z=np.array([.41,.18])
        p,e=transfer(c,f,z,.04)
        q,h=transfer(f,c,p,.04)
        np.testing.assert_allclose(q,z,rtol=1e-13,atol=1e-14)
        self.assertLess(max(e.values()),3e-14)
        self.assertLess(max(h.values()),3e-14)
        self.assertEqual(local_q(f[1],c[0][0],c[0][-1],.04)[0],0.)
        self.assertEqual(local_q(f[1],c[0][0],c[0][-1],.04)[-1],1.)


if __name__=='__main__':
    unittest.main()
