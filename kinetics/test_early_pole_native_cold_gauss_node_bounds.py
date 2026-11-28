import unittest
from pathlib import Path

import numpy as np
from scipy.optimize import linprog

from early_pole_native_cold_gauss_node_bounds import DATE, fraction, run


class GaussNodeBounds(unittest.TestCase):
    def test_sharp_against_independent_linear_program(self):
        nodes = np.array([.15, .2, .45, .65, .9, 1.2])
        means = [.15, .21, .43, .78, 1.2]
        pad = .0331
        for a, b in [(0, 0), (0, 6), (2, 6), (0, 4), (1, 4)]:
            mask = (np.arange(len(nodes)) >= a) & (np.arange(len(nodes)) < b)
            c = mask.astype(float)
            x = fraction(nodes, mask, means, pad)
            for m, found in zip(means, x, strict=True):
                vals = []
                for sign in (1, -1):
                    sol = linprog(sign * c, A_ub=[nodes, -nodes], b_ub=[m+pad, pad-m],
                                  A_eq=[np.ones_like(nodes)], b_eq=[1], bounds=(0., None))
                    self.assertTrue(sol.success)
                    vals.append(sign * sol.fun)
                np.testing.assert_allclose(found, vals, atol=2e-14, rtol=0.)

    def test_no_impossible_mean(self):
        with self.assertRaises(ArithmeticError):
            fraction(np.array([0., 1.]), np.array([True, False]), [1.1], 1e-4)

    def test_exact_archived_resolution_is_not_continuum_resolution(self):
        root = Path(__file__).parent / 'results'
        ref = root / f'early-pole-native-cold-centered-energy-corridors-{DATE}.json'
        z = run(root, ref)
        self.assertEqual([[c['held_sign'] for c in t['charts']] for t in z['rows']],
                         [[-1, 1], [0, 0], [-1, 1], [-1, 1], [-1, 0], [-1, -1]])
        d4 = z['rows'][3]['charts'][1]
        self.assertEqual(d4['continuous_support_held_sign'], 0)
        self.assertTrue(all(s['neutral_interval'][1] < z['diagnostic_occupation']
                            and s['partial_node_panels'] == 2 for s in d4['meshes']))
        self.assertGreater(d4['held_net_hull'][0], 7e-21)
