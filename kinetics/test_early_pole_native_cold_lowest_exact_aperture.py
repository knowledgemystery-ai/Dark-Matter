import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from early_pole_native_cold_lowest_exact_aperture import DATE,TOPOLOGY,contract,finish,positive_affine


class LowestExactAperture(unittest.TestCase):
    def test_cut_only_two_original_children_on_both_meshes(self):
        r=Path(__file__).resolve().parent/'results'
        with TemporaryDirectory() as d:
            p=Path(d)
            for k in (2,4):
                for m in (1,2):
                    _,_,_,t,_=contract(r,p,k,m)
                    j,cut=TOPOLOGY[k]
                    self.assertEqual(t['recalculated_panel_indices'],[j,j+1])
                    self.assertEqual(t['edges'][j:j+3],cut)
            self.assertEqual(list(p.iterdir()),[])

    def test_three_moments_have_independent_closures(self):
        m=.1
        base=np.array([3e-5,4e-10,3e-5-4e-10])
        center=base*1e-9
        a=np.vstack([base,m*base+center,center])
        positive_affine(a,m)
        b=a.copy(); b[1]*=1.00001
        with self.assertRaisesRegex(ArithmeticError,'centered energy'):
            positive_affine(b,m)

    def test_replay_complete_low_band(self):
        r=Path(__file__).resolve().parent/'results'
        a=r/f'cold-lowest-exact-aperture-patches-{DATE}'
        B=json.loads((r/f'early-pole-native-cold-lowest-exact-aperture-{DATE}.json').read_bytes())
        self.assertEqual(finish(r,a),B)
        for z in B['charts']:
            self.assertLess(z['maximum_relative_paired_difference'],6e-5)
            p,q=z['net_box']
            self.assertEqual(z['sign'],int(p>0.)-int(q<0.))
