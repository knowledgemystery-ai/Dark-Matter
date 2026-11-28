import unittest

import numpy as np

from fd_portal_composed_memory_entropy_metric import gain


class EntropyMetric(unittest.TestCase):
    def test_linear_coordinate_invariance_of_same_metrics(self):
        J=np.array([[.6,-1.],[.4,.5],[.12,.9]])
        L=np.array([[1.,.3],[.3,2.]])
        H=np.array([[2.,.2,.4],[.2,1.5,.1],[.4,.1,1.2]])
        P=np.array([[.4,2.],[-.6,1.2]])
        Q=np.array([[1.,.1,1.2],[.6,1.1,.4],[-.3,.4,2.]])
        z=gain(J,L,H)
        w=gain(np.linalg.solve(Q,J@P),P.T@L@P,Q.T@H@Q)
        np.testing.assert_allclose(z,w,rtol=3e-15,atol=1e-15)

    def test_metric_sphere_is_euclidean_when_identity(self):
        J=np.array([[.6,-.5],[.5,.5],[.9,.1]])
        np.testing.assert_allclose(gain(J,np.eye(2),np.eye(3)),np.linalg.svd(J,compute_uv=False),rtol=1e-14)


if __name__=='__main__':
    unittest.main()
