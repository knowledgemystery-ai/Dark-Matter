import unittest

import numpy as np
from scipy.integrate import quad

from early_pole_symmetric_charge_pressure import bounds, cell_average, jensen_loss, kernel


class PressureTest(unittest.TestCase):
    def test_jensen_bound_replays_exact_two_edge_mixture(self):
        a,b,node=map(np.array,([1.,2.],[4.,8.],[2.,4.]))
        m=3.
        A,B,E=[np.hypot(q,m) for q in (a,b,node)]
        t=(E-A)/(B-A)
        chord=(1-t)*kernel(a,m)+t*kernel(b,m)
        np.testing.assert_allclose(kernel(node,m)-jensen_loss(a,b,node,m),chord,rtol=2e-14)

    def test_histogram_quadrature_has_independent_quadrature(self):
        a,b=map(np.array,([1e-9,.001,.05,5.],[2e-9,.006,1.,25.]))
        m=.003
        want=[quad(lambda x:x*x*kernel(x,m),x,y,epsabs=1e-100,epsrel=2e-12)[0]*
              3/((y-x)*(y*y+y*x+x*x)) for x,y in zip(a,b,strict=True)]
        np.testing.assert_allclose(cell_average(a,b,m,64),want,rtol=2e-12)

    def test_counts_do_not_certify_pressure_without_shape(self):
        a,b,node=map(np.array,([1.],[2.],[1.4]))
        v=np.array([[1.],[.95],[.9],[1.05],[1.1]])
        r=bounds(a,b,node,.1,v)[0]
        self.assertLess(r['shared_shape'][1],0.)
        self.assertLess(r['independent_fixed_mean_energies'][1],0.)
        self.assertLess(r['independent_counts_only'][0],0.)
        self.assertGreater(r['independent_counts_only'][1],0.)
