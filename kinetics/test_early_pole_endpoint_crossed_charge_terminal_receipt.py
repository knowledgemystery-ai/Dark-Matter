import unittest
from early_pole_endpoint_crossed_charge_terminal_receipt import split


class CrossedTerminalContrast(unittest.TestCase):
    def test_interaction_and_path_averaged_change(self):
        z=split(1.,3.,5.,10.)
        self.assertEqual(z['average_charge_replacement'],5.5)
        self.assertEqual(z['average_scalar_replacement'],3.5)
        self.assertEqual(z['factorial_interaction'],3.)
        self.assertEqual(z['net_diagonal_change'],9.)

    def test_zero_diagonal_need_not_mean_no_replacement(self):
        z=split(2.,5.,1.,2.)
        self.assertEqual(z['net_diagonal_change'],0.)
        self.assertEqual(z['average_scalar_replacement'],-z['average_charge_replacement'])
        self.assertGreater(z['average_scalar_replacement'],0.)


if __name__=='__main__':
    unittest.main()
