import unittest
from pathlib import Path

from early_pole_native_d8_retained_aperture import run


class D8RetainedAperture(unittest.TestCase):
    def test_exact_state_matches_upper_plan_and_has_nonzero_clipped_inverse(self):
        p=Path(__file__).resolve().parent/'results'
        x=run(p,p/'early-pole-native-d8-raw-band-plan-2026-11-28.json')
        f,b=x['rows']
        self.assertEqual(f['within_selected_raw_plan']['retained_cosine_fraction_range'],[1.,1.])
        self.assertTrue(5e-6<b['over_full']['first_clipped_q']<x['planned_span'][0])
        l,h=b['within_selected_raw_plan']['retained_cosine_fraction_range']
        self.assertTrue(.79<l<h<.91)
