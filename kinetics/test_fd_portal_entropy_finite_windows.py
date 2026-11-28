import unittest

import numpy as np

from fd_portal_entropy_finite_windows import finite_columns, gain_on_plane, window_indices


class FiniteWindowTest(unittest.TestCase):
    def test_reparameterization_keeps_physical_gain(self):
        p=np.array([[1.,0.],[0.,1.],[1.,-2.]])
        G=np.diag([1.,3.,2.])
        A=np.array([[3.,2.],[1.,4.]])
        np.testing.assert_allclose(gain_on_plane(p@A,2*p@A,G,G),[2,2],atol=2e-14)
        np.testing.assert_allclose(gain_on_plane(p,.5*p,G,4*G),[1,1],atol=2e-14)

    def test_rotating_growing_plane_has_singular_gain(self):
        p=np.eye(3)[:,:2]
        q=np.array([[0.,.8],[0.,0.],[1.2,0.]])
        np.testing.assert_allclose(gain_on_plane(p,q,np.eye(3),np.eye(3)),[.8,1.2],atol=2e-14)
        with self.assertRaises(ArithmeticError):
            gain_on_plane(p*0,q,np.eye(3),np.eye(3))

    def test_finite_columns_keep_total_free_after_initial_state(self):
        t=np.array([.99,1.])
        def model(z):
            d=z-np.ones(3)
            return lambda x: np.ones((3,len(x)))+np.outer(d,[1.,2.])+np.outer([0,0,d[0]-d[1]],[2.,3.])
        out=finite_columns(model,model(np.ones(3)),t,.001)
        np.testing.assert_allclose(out[1],[[2,0],[0,2],[3,-3]],atol=2e-13)
        self.assertEqual(window_indices()[-1],(0,3))
