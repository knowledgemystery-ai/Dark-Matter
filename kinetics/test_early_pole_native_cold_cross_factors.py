import json
import unittest
from itertools import combinations
from pathlib import Path

import numpy as np

from early_pole_native_cold_cross_factors import DATE, fractions


class NativeColdCross(unittest.TestCase):
    def test_all_order_factor_closure_for_an_interaction(self):
        w={S:[np.exp(len(S)+(1 in S and 2 in S))] for k in range(4) for S in combinations(range(3),k)}
        t,v=fractions(w)
        np.testing.assert_allclose(t.ravel(),[1.,1.5,1.5])
        np.testing.assert_allclose(v,[4.])

    def test_exact_complete_ledgers_replay_independent_d4_anchor(self):
        p=Path(__file__).resolve().parent/'results'
        r=json.loads((p/f'early-pole-native-d4-transverse-real-pauli-{DATE}.json').read_bytes())
        for q in ('41','82'):
            c=json.loads((p/f'early-pole-native-d2d4-cold-cross-q{q}-{DATE}.json').read_bytes())
            self.assertEqual(len(c['rows']),8)
            self.assertLess(max(z['angular_relative_error'] for z in c['rows']),1e-9)
            t=c['changes']['all_order_log']
            self.assertLess(t[0][0],-2.)
            self.assertLess(abs(t[1][0]),.002)
            if q=='41':
                row=next(x for x in c['rows'] if x['mask']==[0,1,2])
                self.assertEqual(row['rows'][0]['raw_affine_coefficients'],r['complex_and_real_affine_results'][1])
