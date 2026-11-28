import unittest
from pathlib import Path
import numpy as np
from early_pole_native_cold_subband_from_centered_energy import enclosed_cdf,run

class MomentSubband(unittest.TestCase):
    def test_midpoint_energy_leads_to_expected_fraction_bound(self):
        q=np.array([1.,1.5,2.]);m=2.
        a,c,b=q*q/(np.hypot(q,m)+m)
        z=enclosed_cdf(1.,2.,1.5,np.array([(a+c)/2]),m)[0]
        np.testing.assert_allclose(z,[.5,1.],atol=1e-11)
        z=enclosed_cdf(1.,2.,1.5,np.array([(b+c)/2]),m)[0]
        np.testing.assert_allclose(z,[0.,.5],atol=1e-11)

    def test_exact_both_meshes_resolve_previous_unresolved_central_sign(self):
        p=Path(__file__).resolve().parent/'results'
        r=run(p,p/'early-pole-native-cold-mean-from-centered-energy-2026-11-28.json')
        self.assertEqual([[m['sign'] for m in c['meshes']] for c in r['rows']],[[-1,-1],[1,1]])
