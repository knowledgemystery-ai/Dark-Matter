import copy
import unittest

from early_pole_endpoint_continuum_cold_band_pilot import (
    COLD_SPAN, D8_CENTRAL_COLD_BAND, FIRST_SIDE_BAND, special_band_tags)
from early_pole_endpoint_fine_d8_central_cold_receipt import reduce_one


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


if __name__=='__main__':
    unittest.main()
