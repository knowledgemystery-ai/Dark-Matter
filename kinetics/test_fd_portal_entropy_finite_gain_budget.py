import json
import unittest
from math import fsum
from pathlib import Path

from fd_portal_entropy_finite_gain_budget import DATE,TERMS

class FiniteEntropyBudget(unittest.TestCase):
    def test_archived_signed_terms_close_in_both_integrations(self):
        f=Path(__file__).resolve().parent/'results'/f'fd-portal-entropy-finite-gain-budget-{DATE}.json'
        r=json.loads(f.read_bytes())
        self.assertEqual(len(r['rows']),2)
        for v in r['rows']:
            self.assertEqual(list(v['terms']),list(TERMS))
            self.assertAlmostEqual(fsum(v['terms'].values()),v['sum_of_integrals'],places=16)
            self.assertLess(abs(v['sum_of_integrals']-v['endpoint_squared_gain_minus_one']),3e-10)
            self.assertGreater(v['endpoint_gain_minus_one'],0)
