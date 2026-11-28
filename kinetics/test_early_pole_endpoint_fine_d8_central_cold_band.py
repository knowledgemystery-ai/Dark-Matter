import copy
import unittest

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import (
    COLD_SPAN, D8_CENTRAL_COLD_BAND, FIRST_SIDE_BAND, special_band_tags)
from early_pole_endpoint_fine_d8_central_cold_receipt import decisions, reduce_one


class D8CentralBand(unittest.TestCase):
    def test_special_band_scopes_and_fingerprint(self):
        a=dict(grid='fine',x_end=.05,steps=2048,subcells_per_archived_cell=8)
        b=special_band_tags(a,D8_CENTRAL_COLD_BAND)
        self.assertEqual(b['d8_central_cold_q_span'],list(D8_CENTRAL_COLD_BAND))
        self.assertEqual(len(b['d8_central_cold_state_sha256']),64)
        self.assertNotEqual(b,special_band_tags(dict(a,reason='test'),D8_CENTRAL_COLD_BAND))
        self.assertEqual(special_band_tags(a,COLD_SPAN),{})
        for c in (dict(a,subcells_per_archived_cell=4),dict(a,scalar_div=2),dict(a,grid='third'),dict(a,steps=1024)):
            with self.assertRaises(ArithmeticError):
                special_band_tags(c,D8_CENTRAL_COLD_BAND)
        self.assertIn('sideband_state_sha256',special_band_tags(dict(a,subcells_per_archived_cell=2),FIRST_SIDE_BAND))

    def test_reducer_requires_the_complete_identical_native_panel_ledger(self):
        t=dict(scalar_mass=2.,edges=[1.,2.,3.],total=2)
        r=[dict(scalar_mass=2.,total=2,panel=i,bounds=[i+1.,i+2.],
                coeff=[[3.,1.,2.],[9.,3.,6.],[3.,1.,2.]]) for i in range(2)]
        self.assertEqual(reduce_one(r,t)[0].tolist(),[6.,2.,4.])
        with self.assertRaises(ArithmeticError): reduce_one(r[:1],t)
        z=copy.deepcopy(r)
        z[0]['coeff'][0][0]=2.5
        with self.assertRaises(ArithmeticError): reduce_one(z,t)

    def test_decision_includes_both_meshes_and_exact_d4_ranges(self):
        a=np.array([[3.,1.,2.]])
        b=np.array([[6.,2.,4.]])
        d4=dict(forward=[4.,8.],inverse=[4.,9.],conditional_net=[.01,.03])
        out=decisions(a,b,d4,.6)
        self.assertEqual(out['observed_two_mesh_d8_over_d4_forward'],[3/8,6/4])
        self.assertEqual(out['observed_two_mesh_d8_over_d4_inverse'],[1/9,2/4])
        self.assertEqual(out['observed_two_mesh_gross_order_vs_d4'],[0,-1])
        self.assertTrue(all(x>0 for x in out['d8_conditional_net_both_meshes']))
        self.assertEqual(out['d8_conditional_net_sign_from_observed_coefficients'],0)
        self.assertEqual(out['d4_conditional_net_sign'],1)
        self.assertEqual(decisions(a,b,d4,.1)['d8_conditional_net_sign_from_observed_coefficients'],-1)

    def test_unordered_or_negative_reference_rejected(self):
        a=np.array([[3.,1.,2.]])
        for d4 in (dict(forward=[4,3],inverse=[1,2]),dict(forward=[1,3],inverse=[0,2])):
            with self.assertRaises(ArithmeticError):
                decisions(a,a,d4,.1)


if __name__=='__main__':
    unittest.main()
