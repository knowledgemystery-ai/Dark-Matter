import json
from pathlib import Path
import unittest

import numpy as np

from early_pole_native_cold_exact_banded_mean import DATE,run,scalar_sum


class BandedGrossMeans(unittest.TestCase):
    def test_negative_probability_differences_use_reversed_endpoints(self):
        x=scalar_sum(np.array([-.2,.2]),np.array([[1.,2.],[3.,5.]]))
        np.testing.assert_allclose(x,[.2,.8],atol=1e-15,rtol=0.)

    def test_exact_replay_distinguishes_lowest_weighted_mean(self):
        p=Path(__file__).resolve().parent/'results'
        a=run(p)
        self.assertEqual(a,json.loads((p/f'early-pole-native-cold-exact-banded-mean-{DATE}.json').read_bytes()))
        for z in a['by_mesh']:
            for key in ('forward','inverse'):
                w=z['weighted'][key]
                raw=z['raw'][key]
                self.assertLess(w['conditional_D4_minus_D2'][0][0],0.)
                self.assertGreater(w['conditional_D4_minus_D2'][0][1],0.)
                self.assertTrue(all(v[0]>0. for v in raw['conditional_D4_minus_D2']))
                for i in (w,raw):
                    self.assertGreater(i['between_band_weight_redistribution'][0],i['within_band_mean_change'][1])
