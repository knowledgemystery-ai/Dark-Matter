import unittest

from early_pole_endpoint_fine_d4_cold_subband_enclosure import box, ratio, sign


class PositivePanelBounds(unittest.TestCase):
    def test_entire_known_panel_and_unknown_partial_are_distinct(self):
        panels=[[.0,1.],[1.,2.],[2.,3.]]
        masses=[5.,10.,20.]
        self.assertEqual(box(panels,masses,1.,2.),[10.,10.])
        self.assertEqual(box(panels,masses,.9,2.1),[10.,35.])
        self.assertEqual(box(panels,masses,.9,1.1),[0.,15.])

    def test_admissible_positivity_and_coverage(self):
        with self.assertRaises(ArithmeticError):
            box([[0.,1.],[2.,3.]],[1.,2.],0.,3.)
        with self.assertRaises(ArithmeticError):
            box([[0.,1.]],[-1.],0.,1.)
        self.assertIsNone(ratio([1.,2.],[0.,2.]))
        self.assertEqual(ratio([1.,2.],[2.,4.]),[.25,1.])
        self.assertEqual([sign(z) for z in ([0.,1.],[-2.,-1.],[1.,2.])],[0,-1,1])


if __name__=='__main__':
    unittest.main()
