import copy
import unittest

import numpy as np

from early_pole_native_cold_bose_reweight import reweight
from early_pole_native_d8_raw_band_receipt import reduce_rows


class D8RawReceipt(unittest.TestCase):
    def setUp(self):
        self.top=dict(edges=[1.,2.,3.,4.],scalar_mass=2.,total=3,mesh=[8,10,16,24])
        self.raw=np.array([[3.,1.,2.],[9.,3.,6.],[3.,1.,2.]])
        self.rows=[dict(panel=i,bounds=self.top['edges'][i:i+2],scalar_mass=2.,total=3,
                        mesh=self.top['mesh'],coeff=reweight(self.raw,.5+i).tolist()) for i in (1,2)]

    def test_exact_complete_moments(self):
        a=reduce_rows(self.rows,self.top,[1,2],[1.,2.,3.,4.],[.5,1.5,2.5])
        np.testing.assert_allclose(a,2*self.raw)

    def test_rejects_partial_extra_mixed_nonstock_and_stale(self):
        changed=copy.deepcopy(self.rows)
        changed[0]['mesh']=[4,6,10,18]
        foreign=[dict(row,d8_central_cold_q_span=[2.,4.]) for row in self.rows]
        for bad in (self.rows[:1],self.rows+[self.rows[-1]],changed,foreign,self.rows[::-1]):
            with self.assertRaises(ArithmeticError):
                reduce_rows(bad,self.top,[1,2],[1.,2.,3.,4.],[.5,1.5,2.5])
