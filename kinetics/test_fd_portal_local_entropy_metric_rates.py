import unittest
from unittest.mock import patch

import numpy as np

from fd_portal_local_entropy_metric_rates import measure


class LocalEntropyRates(unittest.TestCase):
    def measure(self,f,J,G,B):
        p=np.array([.9,1.,0.,0.])
        motion=np.array([1.,0.,0.,0.])
        def full(_,__):
            return lambda v: np.stack((f+J@(v-1.),np.zeros(3)))
        def metric(_,z):
            return G+(z[0]-p[0])*B,0.
        target='fd_portal_local_entropy_metric_rates.'
        with patch(target+'full',new=full), patch(target+'matrix',new=metric),                 patch(target+'motion',return_value=(motion,0.)):
            return measure(None,p,2e-4,1e-5)

    def test_stationary_metric_local_fractional_frame(self):
        g=np.diag([2.,3.,5.])
        f=np.array([1.,2.,3.])
        J=np.diag([1.5,1.75,4.])
        z=self.measure(f,J,g,np.zeros((3,3)))
        # Moving fraction generator is J minus diag(f); rates for squared norm
        # double those for a norm in the steady metric.
        np.testing.assert_allclose(z['signed_squared_norm_rate'],[-.5,1.,2.],rtol=0.,atol=2e-10)
        np.testing.assert_allclose(z['instantaneous_fixed_total_plane_rate'],[-.5,1.],rtol=0.,atol=2e-10)

    def test_changing_metric_can_reverse_local_growth_witness(self):
        J=np.diag([-.5,-1.,-1.])
        z=self.measure(np.zeros(3),J,np.eye(3),np.diag([2.,0.,0.]))
        np.testing.assert_allclose(z['signed_squared_norm_rate'],[-2.,-2.,1.],rtol=0.,atol=2e-10)
        np.testing.assert_allclose(z['instantaneous_fixed_total_plane_rate'],[-2.,1.],rtol=0.,atol=2e-10)


if __name__=='__main__':
    unittest.main()
