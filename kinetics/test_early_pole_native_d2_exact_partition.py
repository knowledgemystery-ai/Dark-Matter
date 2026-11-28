import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from early_pole_native_d2_exact_partition import DATE,EDGES,IDS,SPLITS,finish,prepare


class D2ExactFiveBands(unittest.TestCase):
    def test_exact_planned_later_split_preserves_low_child_positions(self):
        root=Path(__file__).resolve().parent/'results'
        with TemporaryDirectory() as d:
            for m in (1,2):
                _,_,low,_,new,_=prepare(root,Path(d),m)
                self.assertEqual(low['edges'][:11],new['edges'][:11])
                self.assertEqual(new['recalculated_panel_indices'],IDS)
                self.assertEqual([new['edges'][j:j+3] for j in IDS[::2]],SPLITS)
                self.assertEqual(new['total'],52)
            self.assertFalse(list(Path(d).iterdir()))

    def test_admitted_partition_replay(self):
        r=Path(__file__).resolve().parent/'results'
        x=r/f'd2-exact-partition-patches-{DATE}'
        B=json.loads((r/f'early-pole-native-d2-exact-partition-{DATE}.json').read_bytes())
        self.assertEqual(finish(r,x),B)
        self.assertEqual(B['edges'],list(EDGES))
        for a in B['mesh_runs']:
            self.assertEqual(sum(a['counts']),52)
