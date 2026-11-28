import unittest

import numpy as np

from fd_portal_direct_transient_entropy_witness import centered_rhs, pair


class CenteredPair(unittest.TestCase):
    def test_matches_independent_affine_trajectories(self):
        a=np.array([-.2,-2.,-16.])
        b=np.array([.05,.3,2.])
        d=np.array([.8,-.6,0.])
        dt=.42
        m,s,g=pair(lambda _,y:a*y+b,.55,.55+dt,d,1e-4,.008)
        fac=np.exp(a*dt)
        np.testing.assert_allclose(m,fac+(b/a)*(fac-1.),rtol=2e-13,atol=5e-14)
        np.testing.assert_allclose(s,fac*d,rtol=2e-12,atol=5e-13)
        np.testing.assert_allclose(g,a*fac*d,rtol=2e-11,atol=5e-12)

    def test_nonlinear_even_and_odd_part_are_preserved(self):
        m=np.array([.9,1.1,.25])
        z=np.array([.8,-.6,.2])
        e=2e-4
        a,b=centered_rhs(lambda _,x:x*x,.9,m,z,e)
        np.testing.assert_allclose(a,m*m+e**2*z*z,rtol=3e-15)
        np.testing.assert_allclose(b,2*m*z,rtol=1e-11)


if __name__=='__main__':
    unittest.main()
