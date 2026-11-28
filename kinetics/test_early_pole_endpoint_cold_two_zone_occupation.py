import unittest
from pathlib import Path

from early_pole_endpoint_cold_two_zone_occupation import (
    infimum_and_witness, run, summarize)


class TwoZoneOccupation(unittest.TestCase):
    def test_exact_common_cut_rejects_unknown_split_panel(self):
        data=[dict(bounds=[0.,2e-5],coeff=[[3.,1.,2.]])]
        with self.assertRaises(ArithmeticError):
            summarize(data)

    def test_open_limit_has_known_two_zone_optimum(self):
        D2=[[[1,1],[4,1]]]
        D4=[[[1,1],[3,1]]]
        r=infimum_and_witness(D2,D4,.01)
        self.assertEqual(r['two_zone_limit'],[2,3])
        self.assertEqual(r['two_zone_infimum_of_upper_minus_lower'],1)
        self.assertGreater(r['witness_range'],1)
        self.assertTrue(all(x<0 for x in r['witness_d2_net_full_and_upper_by_mesh'][0]))
        self.assertTrue(all(x>0 for x in r['witness_d4_net_full_and_upper_by_mesh'][0]))

    def test_two_zone_construction_is_close_to_any_shape_lower_bound(self):
        p=Path(__file__).resolve().parent/'results'
        r=run(p,p/'early-pole-cold-neutral-occupation-corridors-2026-11-28.json')
        self.assertTrue(1<r['two_zone_infimum_over_arbitrary_profile_lower_bound']<1.000005)
