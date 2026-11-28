import hashlib
import json
from pathlib import Path
from shutil import copytree
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from early_pole_native_d8_upper_fine_canary import CellChart, DATE, SPECS, load, read, run


class D8UpperFineCanary(unittest.TestCase):
    def setUp(self):
        self.results = Path(__file__).resolve().parent / 'results'
        self.work = self.results / f'd8-upper-fine-canary-checkpoints-{DATE}'
        self.z = json.loads((self.results / f'early-pole-native-d8-upper-fine-canary-{DATE}.json').read_bytes())

    def test_exact_raw_files_reconstruct_all_nine_moments(self):
        stock = (self.work / 'stock.jsonl').read_bytes()
        self.assertEqual(json.loads(stock),self.z['raw_stock'])
        key = self.z['contract']
        self.assertEqual(key['exact_stock_row_sha256'],hashlib.sha256(stock).hexdigest())
        self.assertEqual(key['panel'],183)
        self.assertEqual(key['mesh'],[8,10,16,24])
        nodes = [json.loads(p.read_bytes()) for p in sorted(self.work.glob('real-node-*.json'))]
        self.assertEqual([r['contract']['node'] for r in nodes],list(range(8)))
        card,h = read(self.results / f'early-pole-cell-endpoint-fine-d8-x05-state-{DATE}.json.gz')
        self.assertEqual(key['exact_native_state_sha256'],h)
        d = CellChart(load('fine')[0],SPECS['fine'],8,card['reference_x'])
        m = self.z['raw_stock']['scalar_mass']
        q = np.array([r['contract']['q'] for r in nodes])
        g = np.array([r['contract']['outer_weight'] for r in nodes])
        self.assertTrue(key['bounds'][0]<q.min()<q.max()<key['bounds'][1])
        self.assertAlmostEqual(g.sum()/(key['bounds'][1]-key['bounds'][0]),1.)
        f = np.asarray(card['endpoint_occupation'])[np.searchsorted(d.ds.ed,q,side='right')-1]
        G,J,S = np.asarray([r['raw_angular_values'] for r in nodes]).T
        np.testing.assert_allclose(S,G-J,rtol=5e-13,atol=0.)
        e = np.hypot(q,m)
        A = np.array([np.ones_like(q),e,q*q/(e+m)])
        B = np.array([(1+f)*G,f*J,f*S+G])
        C = (A*(q*q*g/e)) @ B.T
        np.testing.assert_allclose(C,self.z['raw_stock']['coeff'],rtol=2e-12,atol=0.)
        np.testing.assert_allclose(C,self.z['exact_real_pauli_mirror_coefficients'],rtol=2e-15,atol=0.)

    def test_full_checkpoint_replay_and_tamper_rejection(self):
        with TemporaryDirectory() as d:
            p = Path(d)/'cache'
            copytree(self.work,p)
            self.assertEqual(run(self.results,p),self.z)
            a = p/'real-node-00.json'
            z = json.loads(a.read_bytes())
            z['contract']['q'] *= 1.00001
            a.write_text(json.dumps(z,allow_nan=False))
            with self.assertRaisesRegex(ArithmeticError,'Mixing|mixing'):
                run(self.results,p)
