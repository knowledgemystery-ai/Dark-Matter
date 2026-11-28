import unittest
from pathlib import Path

from early_pole_endpoint_cold_affine_occupation import Q0, Q1, moments, run


class AffineOccupation(unittest.TestCase):
    def test_unknown_within_panel_slopes_bound_normalized_momentum(self):
        m=(Q0+Q1)/2
        p=[dict(bounds=[Q0,m],coeff=[[3.,1.,2.]]),
           dict(bounds=[m,Q1],coeff=[[5.,1.,4.]])]
        t=moments(p)
        self.assertEqual(t[0],1/3)
        self.assertAlmostEqual(t[1],1/3)
        self.assertAlmostEqual(t[2],5/6)

    def test_replays_robust_smooth_witness_on_all_accepted_panels(self):
        p=Path(__file__).resolve().parent/'results'
        s='early-pole-cold'
        r=run(p,p/f'{s}-neutral-occupation-corridors-2026-11-28.json',p/f'{s}-two-zone-occupation-2026-11-28.json')
        d2,d4=r['witness_net_intervals_by_chart_mesh_and_band']
        self.assertTrue(all(i[1]<0 for mesh in d2 for i in mesh))
        self.assertTrue(all(i[0]>0 for mesh in d4 for i in mesh))
        self.assertTrue(.006<r['witness_change_relative_to_low_endpoint']<.007)
