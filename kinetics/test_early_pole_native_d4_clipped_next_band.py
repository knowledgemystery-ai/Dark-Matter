import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

import numpy as np

from early_pole_native_d4_clipped_next_band import IDS, SPLITS, prepare, verify_children


class D4Clipping(unittest.TestCase):
    def test_original_topologies_and_patch_ids(self):
        root = Path(__file__).resolve().parent / 'results'
        with TemporaryDirectory() as d:
            dest = Path(d)
            for m in (1, 2):
                plan = prepare(root, m, dest)[-1]
                self.assertEqual(plan['total'], 110)
                self.assertEqual(plan['recalculated_panel_indices'], IDS)
                self.assertEqual([plan['edges'][i:i+2] for i in IDS],
                                 [[lo, q] if j == 0 else [q, hi] for lo, q, hi in SPLITS for j in (0, 1)])
            self.assertEqual(list(dest.iterdir()), [])

    def test_parent_reclosure_checks_small_inverse_column(self):
        key = np.array([[2e-11, 4e-20, 2e-11-4e-20]]) * np.array([[1.], [1e-3], [1e-8]])
        with TemporaryDirectory() as d:
            p = Path(d) / 'base.jsonl'
            p.write_text(''.join(json.dumps(dict(bounds=[lo, hi], coeff=key.tolist()))+'\n'
                                 for lo, _, hi in SPLITS))
            rec = [dict(panel=i, coeff=(key / 2).tolist()) for i in IDS]
            self.assertEqual(verify_children(p, rec), [0., 0.])
            with self.assertRaisesRegex(ArithmeticError, 'all four'):
                verify_children(p, rec[:-1])
            rec[1]['coeff'][0][1] *= 1.001
            with self.assertRaisesRegex(ArithmeticError, 'local closure'):
                verify_children(p, rec)

    def test_completed_archived_patch_replays_without_additional_quadrature(self):
        from shutil import copy2
        from early_pole_native_d4_clipped_next_band import DATE, finish

        root = Path(__file__).resolve().parent / 'results'
        source = root / f'd4-clipped-next-band-patches-{DATE}'
        reported = json.loads((root / f'early-pole-native-d4-clipped-next-band-{DATE}.json').read_bytes())
        with TemporaryDirectory() as d:
            p = Path(d)
            for m in (1, 2):
                copy2(source / f'd4-cold-split-{m}.jsonl',p)
            self.assertEqual(finish(root,p),reported)
        for m in (1, 2):
            rows = [json.loads(z) for z in (source / f'd4-cold-split-{m}.jsonl').read_bytes().splitlines()]
            self.assertEqual([z['panel'] for z in rows], IDS)
            for z in rows:
                X = np.array(z['coeff'])
                F, B, S = X.T
                np.testing.assert_allclose(F,B+S,rtol=3e-13,atol=0.)
                np.testing.assert_allclose(X[1],z['scalar_mass']*X[0]+X[2],rtol=3e-13,atol=0.)
