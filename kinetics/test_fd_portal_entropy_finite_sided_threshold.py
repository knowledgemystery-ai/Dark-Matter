import json
from pathlib import Path
import unittest

import numpy as np

from fd_portal_entropy_finite_sided_threshold import DATE, boundary


def determinant_crossings(H, fixed, moved, f):
    a, b = H + f * fixed, moved
    c2 = np.linalg.det(b)
    c1 = a[0, 0]*b[1, 1] + a[1, 1]*b[0, 0] - 2*a[0, 1]*b[0, 1]
    c0 = np.linalg.det(a)
    return np.sort(np.roots([c2, c1, c0]))


class FiniteSided(unittest.TestCase):
    def test_threshold_is_monotone_until_fixed_side_alone_expands(self):
        H = -np.diag([1., 3.])
        A = np.diag([2., .4])
        B = np.diag([.5, 1.])
        self.assertAlmostEqual(boundary(H, A, B, 0.), 2.)
        self.assertAlmostEqual(boundary(H, A, B, .25), 1.)
        self.assertEqual(boundary(H, A, B, .5), 0.)
        self.assertIsNone(boundary(H, A, B, .75))

    def test_congruence_preserves_boundary(self):
        H = -np.array([[1., .5], [.5, 2.]])
        A = np.array([[.7, -.1], [-.1, .2]])
        B = np.array([[.3, .1], [.1, .84]])
        T = np.array([[3., 2.], [5., 4.]])
        expected = boundary(H, A, B, .9)
        actual = boundary(*[T.T @ V @ T for V in (H, A, B)], .9)
        self.assertAlmostEqual(actual, expected, places=12)

    def test_report_matches_both_numerical_settings(self):
        p = Path(__file__).parent/'results'/f'fd-portal-entropy-finite-sided-threshold-{DATE}.json'
        Z = json.loads(p.read_bytes())
        ref = p.with_name(f'fd-portal-entropy-finite-rank-threshold-{DATE}.json')
        old = json.loads(ref.read_bytes())
        for a, prev in zip(Z['rows'], old['rows'], strict=True):
            self.assertLess(max(a['consistency_with_uniform_receipt']), 3e-12)
            R = np.array(a['hot_then_cool_integrated_excess_gram'])
            H = np.array(prev['integrated_reference_matrix'])
            for blocks, summary in zip((R, R[::-1]), a['fixed_hot_then_fixed_cool'], strict=True):
                for f, c in zip(Z['fixed_multipliers'], summary['critical_other_multiplier'], strict=True):
                    roots = determinant_crossings(H, *blocks, f)
                    if c is None:
                        self.assertLess(roots[0], 0.)
                    else:
                        self.assertAlmostEqual(roots[0], c, places=10)
        vals = [[v['critical_other_multiplier'][-1] for v in a['fixed_hot_then_fixed_cool']] for a in Z['rows']]
        self.assertEqual([v is None for v in vals[0]], [v is None for v in vals[1]])
        for a, b in zip(*vals, strict=True):
            if a is not None:
                self.assertAlmostEqual(a, b, places=6)
