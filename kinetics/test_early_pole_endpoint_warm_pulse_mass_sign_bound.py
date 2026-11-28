import math
import unittest

import numpy as np

from early_pole_endpoint_warm_pulse_mass_sign_bound import bounds, cover


class NodalPressureMassBounds(unittest.TestCase):
    def test_both_signs_are_in_the_same_mixed_weight_enclosure(self):
        q=np.array([.1,1.])
        d=np.array([1.,-.05])
        lo,hi=bounds(q,d,.01,10.)
        values=[math.fsum(d*q*q/np.hypot(q,m)) for m in np.geomspace(.01,10.,31)]
        self.assertLess(lo,min(values))
        self.assertGreater(hi,max(values))
        z=cover(q,d,.01,10.)
        self.assertEqual(z['boundary_signs'],[1,-1])
        self.assertFalse(z['complete'])

    def test_uniform_positive_or_negative_occupation_delta(self):
        q=np.array([.02,.08,.55])
        for sign in (-1,1):
            d=sign*np.array([.7,.14,.1])
            z=cover(q,d,.2,.8)
            self.assertTrue(z['complete'])
            self.assertEqual(z['boundary_signs'],[sign,sign])
            self.assertEqual(z['mass_panel_count'],1)
            self.assertTrue(all(sign*p>0. for p in z['difference_enclosure']))

    def test_invalid_mass(self):
        with self.assertRaises(ArithmeticError):
            bounds(np.array([1.]),np.array([.2]),0.,1.)


if __name__=='__main__':
    unittest.main()
