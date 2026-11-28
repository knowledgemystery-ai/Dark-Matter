import unittest
from pathlib import Path

import numpy as np

from early_pole_native_cold_bose_reweight import address, reweight, run, strip


class BoseReweighting(unittest.TestCase):
    def test_finite_factor_not_large_occupation_approximation(self):
        P=np.array([[3.,1.,2.],[2.,.5,1.5]])
        f=np.array([.1,10.])
        W=reweight(P,f)
        np.testing.assert_allclose(W[:,2],[3.2,17.])
        np.testing.assert_allclose(strip(W,f),P)
        with self.assertRaises(ArithmeticError):
            address([0.,1.,2.],np.array([[.5,1.5]]))

    def test_complete_full_and_upper_panels_respond_differently(self):
        a,b,c=run(Path(__file__).resolve().parent/'results')['rows']
        self.assertLess(a['weighted_native_d4_over_d2'][0],a['unweighted_d4_over_d2'][0])
        self.assertLess(a['all_orders_raw_chart_and_explicit_occupation_factors'][1][0],1.)
        self.assertGreater(c['all_orders_raw_chart_and_explicit_occupation_factors'][1][0],1.1)
        self.assertTrue(np.all(np.asarray(a['conditional_net_at_ladder_held_occupation'])[0]<0.))
        self.assertTrue(np.all(np.asarray(a['conditional_net_at_ladder_held_occupation'])[1]>0.))
