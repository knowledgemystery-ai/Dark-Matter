import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from early_pole_native_d8_planned_upper_real import context, load_row, reduce, save, work


class D8PlannedUpper(unittest.TestCase):
    def test_checksummed_contract_rejects_changes(self):
        with TemporaryDirectory() as d:
            p = Path(d)/'node.json'
            meta = dict(panel=157,request='disjoint')
            save(p,dict(contract=meta,angle=[2.,1.,1.]))
            self.assertEqual(load_row(p,meta),[2.,1.,1.])
            with self.assertRaisesRegex(ArithmeticError,'Mixed or corrupted'):
                load_row(p,dict(meta,panel=158))
            q=json.loads(p.read_bytes())
            q['record']['angle'][0]+=1
            p.write_text(json.dumps(q))
            with self.assertRaisesRegex(ArithmeticError,'Mixed or corrupted'):
                load_row(p,meta)

    def test_both_seed_checkpoints_and_fail_closed_on_partial(self):
        root=Path(__file__).resolve().parent/'results'
        with TemporaryDirectory() as d:
            p=Path(d)
            for k,j in ((1,157),(2,183)):
                W=p/f'mesh{k}'
                meta,edges,*_=context(root,k)
                self.assertEqual(len(edges),196)
                work(root,W,k,True)
                files=sorted(W.glob(f'upper-{j}-*.json'))
                self.assertEqual(len(files),meta['mesh'][0])
                self.assertEqual(len(list(W.glob('upper-*.json'))),meta['mesh'][0])
            with self.assertRaisesRegex(ArithmeticError,'partial'):
                reduce(root,p)
