"""Rebuild the two complete D4 cold ledgers without recomputing their integrals."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path
from tempfile import TemporaryDirectory

from early_pole_endpoint_continuum_cold_band_pilot import run as pilot
from early_pole_endpoint_native_cold_shards import join,plan,work
from early_pole_endpoint_scalar_subshell_chain_receipt import read

D='2026-11-28'


def run(root,ladder):
    d=json.loads(ladder.read_bytes())
    native=root/f'early-pole-cell-endpoint-fine-d4-x05-state-{D}.json.gz'
    a,h=read(native)
    if h!=d['cards'][2]['state_sha256_uncompressed']:
        raise ArithmeticError('Must replay the exact archived native D4 state')
    rows=[]
    with TemporaryDirectory() as v:
        tmp=Path(v)
        for i in (1,2):
            t=pilot(a,i,None,topology_only=True)
            m=plan(t,4)
            data=gzip.decompress((root/f'fine-d4-endpoint-continuum-cold-{i}.jsonl.gz').read_bytes())
            sha=hashlib.sha256(data).hexdigest()
            if sha!=d['raw_uncompressed_sha256'][2][i-1]:
                raise ArithmeticError('Original exact complete D4 ledger required')
            raw=data.splitlines(keepends=True)
            names=[]
            for j,indices in enumerate(m['shards']):
                p=tmp/f'mesh{i}-part{j}.jsonl'
                p.write_bytes(b''.join(raw[k] for k in indices))
                names.append(p)
            n=names[1].read_bytes()
            work(m,native,1,names[1])
            if n!=names[1].read_bytes() or join(m,names)!=data:
                raise ArithmeticError('Native ledgers must be exactly preserved')
            rows.append(dict(mesh=t['mesh'],full_panel_count=t['total'],
                         counts_per_shard=[len(x) for x in m['shards']],topology_sha256=m['topology_sha256'],
                         manifest_sha256=m['payload_sha256'],
                         archive_uncompressed_and_join_sha256=sha,byte_exact=True,resumed_full_shard_noop=1))
    return dict(rows=rows,scope='Independent resumption and arithmetic-free exact assembly on complete D4 ledgers; no intervention in running D8 calculations or new physical result.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('ladder',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
