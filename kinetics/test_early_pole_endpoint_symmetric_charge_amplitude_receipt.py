import unittest

from early_pole_endpoint_symmetric_charge_amplitude_receipt import compare


class HalfAmplitudeContrast(unittest.TestCase):
    def test_exact_linear_odd_quadratic_even(self):
        z=compare([-2.,.5],[-4.,2.])
        self.assertEqual(z['half_minus_linear_quadratic_target'],[0.,0.])
        self.assertEqual(z['ratio_to_linear_quadratic_target'],[1.,1.])

    def test_small_baseline_does_not_create_a_spurious_ratio(self):
        a=compare([1.,2.],[0.,1e-20])
        self.assertEqual(a['ratio_to_linear_quadratic_target'],[None,None])


if __name__=='__main__':
    unittest.main()
