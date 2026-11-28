import unittest
from pathlib import Path

from early_pole_native_cold_centered_energy_corridors import run,roots

class EnergyCorridors(unittest.TestCase):
    def test_neutral_uses_crossed_interval_endpoints(self):
        self.assertEqual(roots([[7,9],[2,3],[4,5]]),[.4,.75])
    def test_exact_admitted_moments_add_only_resolved_signs(self):
        p=Path(__file__).resolve().parent/'results'
        z=run(p,p/'early-pole-native-cold-subband-from-centered-energy-2026-11-28.json')
        self.assertEqual([[c['held_sign'] for c in b['charts']] for b in z['rows']],
                         [[-1,1],[0,0],[-1,1],[-1,0],[-1,0],[-1,-1]])
        a,b=z['full_central_and_next_simultaneous_sufficient']
        self.assertLess(a,2.595e-9)
        self.assertGreater(b,2.595e-9)
