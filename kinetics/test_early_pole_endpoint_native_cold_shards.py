import copy
import json
from pathlib import Path
from tempfile import TemporaryDirectory
import unittest

from early_pole_endpoint_native_cold_shards import join,plan,verify


class DisjointNativeShards(unittest.TestCase):
    def test_complete_shards_restore_exact_bytes_and_reject_mutations(self):
        t=dict(grid='fine',mesh=[4,6,10,18],total=3,edges=[1.,2.,3.,4.])
        m=plan(t,2)
        self.assertEqual(m['shards'],[[0,2],[1]])
        p=copy.deepcopy(m)
        p['shards'][0]=[0,1]
        with self.assertRaises(ArithmeticError):
            verify(p)
        z=[(json.dumps(dict(grid='fine',mesh=t['mesh'],total=3,panel=i,bounds=t['edges'][i:i+2],coeff=[[0.,0.,0.]]*3))+'\n').encode() for i in range(3)]
        with TemporaryDirectory() as tmp:
            files=[Path(tmp)/str(i) for i in range(2)]
            files[0].write_bytes(z[0]+z[2])
            files[1].write_bytes(z[1])
            self.assertEqual(join(m,files),b''.join(z))
            with self.assertRaises(ArithmeticError):
                join(m,files[::-1])
            files[1].write_bytes(z[2])
            with self.assertRaises(ArithmeticError):
                join(m,files)

    def test_rejects_invalid_native_geometry(self):
        with self.assertRaises(ArithmeticError):
            plan(dict(mesh=[4,6,10,18],total=2,edges=[1.,2.,2.]),2)


if __name__=='__main__':
    unittest.main()
