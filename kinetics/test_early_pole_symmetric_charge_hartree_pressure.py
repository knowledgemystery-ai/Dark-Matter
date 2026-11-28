import unittest
from pathlib import Path

from early_pole_symmetric_charge_hartree_pressure import run


class HartreePressure(unittest.TestCase):
    def test_exact_trace_and_pressure_sign_reversal(self):
        p=Path(__file__).resolve().parent/'results'
        s='early-pole-symmetric-charge'
        q=run(p,p/f'{s}-cutoff-audit-2026-11-28.json.gz',p/f'{s}-terminal-energy-2026-11-28.json')
        for r in q['rows']:
            self.assertLess(r['differential_trace_identity_abs_defect'],1e-17)
        self.assertTrue(all(r['corrected_scalar_pressure']<0 for r in q['rows'][:-1]))
        last=q['rows'][-1]
        self.assertLess(last['fixed_nominal_mass_scalar_kinetic_pressure'],0)
        self.assertGreater(last['corrected_scalar_pressure'],0)
        self.assertGreater(last['corrected_total_model_pressure'],0)
        self.assertLess(last['corrected_total_energy'],0)
