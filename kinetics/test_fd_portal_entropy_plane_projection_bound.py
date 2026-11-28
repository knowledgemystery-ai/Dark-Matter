import unittest

import numpy as np

from fd_portal_entropy_plane_projection_bound import comparison, projection_rate_bound_from_fractional_sine


class PlaneRateBound(unittest.TestCase):
    def test_identical_planes_with_different_basis(self):
        E=np.eye(3)[:,:2]
        T=np.array([[2.,1.],[-3.,1.]])
        z=comparison(np.diag([1.,2.,4.]),np.diag([0.,-2.,3.]),E,E@T)
        self.assertLess(z['largest_principal_sine'],2e-15)
        np.testing.assert_allclose(z['rate_first'],z['rate_second'],rtol=0.,atol=2e-15)

    def test_off_diagonal_metric_invariance(self):
        e=np.eye(3)
        G=np.array([[2.,.1,.2],[.1,1.5,-.1],[.2,-.1,.9]])
        S=np.array([[-3.,.5,1.],[.5,.3,.2],[1.,.2,1.6]])
        U,V=e[:,:2],e[:,1:]
        a=comparison(G,S,U,V)
        R=np.array([[1.,.2,-.3],[-1.,1.4,.2],[.1,-.4,.8]])
        b=comparison(R.T@G@R,R.T@S@R,np.linalg.solve(R,U),np.linalg.solve(R,V))
        for k in a:
            np.testing.assert_allclose(a[k],b[k],rtol=1e-13,atol=3e-15)

    def test_raw_angle_bound_is_safe_after_metric_rescaling(self):
        e=np.eye(3)
        U=e[:,:2]
        V=np.column_stack((e[:,0],e[:,1]+.2*e[:,2]))
        G=np.diag([.5,2.,10.])
        S=np.diag([-5.,2.,1.])
        d=comparison(e,S,U,V)['largest_principal_sine']
        a,b=projection_rate_bound_from_fractional_sine(G,S,d)
        exact=comparison(G,S,U,V)
        self.assertGreaterEqual(a,exact['largest_principal_sine'])
        self.assertGreaterEqual(b,exact['rayleigh_rate_change_bound'])


if __name__=='__main__':
    unittest.main()
