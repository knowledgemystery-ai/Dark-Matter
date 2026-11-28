import unittest

import numpy as np

from early_pole_endpoint_charge_nested_moment_remap import local_q
from early_pole_endpoint_charge_symmetric_null_terminal import balanced_nodes


class TrulySymmetricChargeNull(unittest.TestCase):
    def test_opposite_single_vector_preserves_both_moments_and_caps(self):
        p=np.array([.001,.011,.011,.040])
        w=np.array([.16,.21,.39,.24])
        f=np.array([.42,.31,.31,.16])
        d=balanced_nodes(p,w,f,.04,p[[0,-1]])
        e=local_q(p,p[0],p[-1],.04)
        self.assertGreater(d[1],0.)
        self.assertEqual(d[1],d[2])
        np.testing.assert_allclose([w@d,(w*e)@d],0.,rtol=0.,atol=7e-18)
        np.testing.assert_allclose(((f+d)+(f-d))*.5,f,rtol=0.,atol=6e-17)
        self.assertTrue(np.all((f-abs(d)>=0)&(f+abs(d)<=1)))

    def test_zero_cell_has_no_direction(self):
        p=np.array([.001,.011,.011,.04])
        np.testing.assert_array_equal(balanced_nodes(p,np.ones(4),np.zeros(4),.04,p[[0,-1]]),np.zeros(4))


if __name__=='__main__':
    unittest.main()
