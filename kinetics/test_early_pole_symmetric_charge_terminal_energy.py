import unittest
from pathlib import Path

from early_pole_symmetric_charge_terminal_energy import run


class AdmittedEnergyLedger(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        root=Path(__file__).resolve().parent/'results'
        cls.rows=run(root,root/'early-pole-symmetric-charge-cutoff-audit-2026-11-28.json.gz')['rows']

    def test_pressure_drop_does_not_imply_scalar_energy_drop(self):
        row=self.rows[-1]
        self.assertLess(row['scalar_fixed_mass_three_pressure'],0.)
        self.assertGreater(row['scalar_fixed_mass_energy'],0.)
        self.assertGreater(row['scalar_on_own_pole_energy'],0.)
        self.assertLess(row['total_corrected_energy'],0.)
        self.assertLess(row['parent_energy'],0.)
        self.assertGreater(row['daughter_energy'],0.)

    def test_all_paired_expansion_trace_ledgers_close(self):
        for r in self.rows:
            self.assertEqual(len(r['cumulative_saved_trace_differences']),8)
            self.assertLess(abs(r['final_total_minus_initial_minus_trace']),3e-15)
            self.assertLess(abs(r['initial_total_corrected_energy']),1e-16)
            self.assertLess(abs(r['stable_state_total_absolute_defect']),2e-18)
