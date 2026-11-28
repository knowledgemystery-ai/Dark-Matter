import unittest
from pathlib import Path

import numpy as np

from early_pole_native_d8_raw_band_plan import run, select


class D8RawBandPlan(unittest.TestCase):
    def test_enumeration_minimizes_native_panel_count_not_intracell_shape(self):
        e=[0.,1.,2.,3.]
        Q=np.c_[e[:-1],e[1:]]
        T=[(Q,np.array([.1,.2,.7]))]
        self.assertEqual(select(e,T,.75)['panel_indices'],[1,2])
        T.append((Q,np.array([.7,.2,.1])))
        self.assertEqual(select(e,T,.75)['panel_indices'],[0,1,2])

    def test_exact_topology_yields_a_plan_not_a_new_d8_source(self):
        p=Path(__file__).resolve().parent/'results'
        ts=run(p,p/'early-pole-native-cold-weight-target-2026-11-28.json')
        self.assertEqual(ts['full_d8_panel_count'],195)
        self.assertEqual(ts['planned']['panel_indices'],list(range(157,195)))
        self.assertGreater(min(x[0] for x in ts['planned']['coverage']),.8)
