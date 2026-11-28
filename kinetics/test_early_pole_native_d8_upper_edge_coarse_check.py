import json
from pathlib import Path
from tempfile import TemporaryDirectory
from shutil import copytree
import unittest

import numpy as np

from early_pole_native_d8_upper_edge_coarse_check import DATE, run


class UpperD8EdgeCoarse(unittest.TestCase):
    def test_replay_complete_archived_coarse_panel(self):
        root = Path(__file__).resolve().parent/'results'
        cache = root/f'd8-upper-edge-coarse-checkpoints-{DATE}'
        R = json.loads((root/f'early-pole-native-d8-upper-edge-coarse-check-{DATE}.json').read_bytes())
        self.assertEqual(R['binding']['panel'],157)
        self.assertEqual(R['binding']['mesh'],[4,6,10,18])
        self.assertLess(R['componentwise_max_relative_difference'],2e-12)
        raw = [json.loads(p.read_bytes()) for p in sorted(cache.glob('edge-real-*.json'))]
        self.assertEqual([q['binding']['real_node'] for q in raw],list(range(4)))
        np.testing.assert_allclose(np.asarray(R['stock']['coeff']),np.asarray(R['real']),rtol=2e-12,atol=0.)
        with TemporaryDirectory() as d:
            W = Path(d)/'checkpoints'
            copytree(cache,W)
            self.assertEqual(run(root,W),R)
            p = W/'edge-real-0.json'
            q = json.loads(p.read_bytes())
            q['binding']['plan_sha256'] = '0' * 64
            p.write_text(json.dumps(q))
            with self.assertRaisesRegex(ArithmeticError,'mix edge checkpoints'):
                run(root,W)
