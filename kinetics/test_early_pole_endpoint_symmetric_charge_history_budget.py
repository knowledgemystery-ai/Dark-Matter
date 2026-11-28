import unittest
from early_pole_endpoint_symmetric_charge_history_budget import classify


def now(values,x=(1.,2.,3.)):
    return [dict(x=t,scalar_yield=y,continuity_residual=0.,pair_number_error=0.) for t,y in zip(x,values,strict=True)]


class SavedSymmetricHistory(unittest.TestCase):
    def test_exact_sign_brackets(self):
        a,b,c=now([3.,2.,4.]),now([3.,3.,3.]),now([2.,4.,1.])
        d=classify(a,b,c)
        self.assertEqual(d['brackets_with_opposite_saved_odd_signs'],[[1.,2.],[2.,3.]])
        self.assertEqual([x['symmetric_even_yield'] for x in d['rows']],[-.5,0.,-.5])

    def test_detects_unmatched_checkpoint_clocks(self):
        with self.assertRaises(ArithmeticError):
            classify(now([1.,1.,1.]),now([1.,1.,1.]),now([1.,1.,1.],x=(1.,2.1,3.)))


if __name__=='__main__':
    unittest.main()
