import unittest

import numpy as np

from early_pole_symmetric_charge_cutoff_audit import contrasts, interior_enclosures


class CutoffAuditTest(unittest.TestCase):
    def test_same_cell_mass_can_hide_resolved_edge_contrast(self):
        e=[1.,2.,3.]
        a=[1.,2.,1.]
        V=np.array([[1.,1.,.01],[.9,1.,.01],[1.1,1.,.01],[.8,1.,.01],[1.2,1.,.01]])
        r=contrasts(e,V,a)
        i=interior_enclosures(e,V,a,r)
        self.assertEqual(r[1]['signs'],[-1,1,-1,1])
        self.assertEqual(i[1]['certified_signs'],[0,0,0,0])
        self.assertEqual(i[2]['certified_signs'],[-1,1,-1,1])
        np.testing.assert_allclose(i[2]['lower'],[-.11,.09,-.21,.19],atol=2e-16)
        np.testing.assert_allclose(i[2]['upper'],[-.09,.11,-.19,.21],atol=2e-16)
        self.assertEqual(i[2]['q_interval'],[2.,3.])

    def test_rejects_non_native_or_nonphysical_input(self):
        with self.assertRaises(ArithmeticError):
            contrasts([1,1],np.ones((5,2)),[1,1])
        with self.assertRaises(ArithmeticError):
            contrasts([1,2],-np.ones((5,2)),[1,1])
