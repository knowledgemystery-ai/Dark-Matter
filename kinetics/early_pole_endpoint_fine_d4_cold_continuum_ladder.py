"""Compare complete cold-band continuum integrals for the matched fine D1/D2/D4 histories."""
import argparse
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_band_panel_contrast import CUTS, PROBS, pair

DATE='2026-11-28'


def unpacked_hash(p):
    b=p.read_bytes()
    return hashlib.sha256(gzip.decompress(b) if p.suffix=='.gz' else b).hexdigest()


def run(root,d4):
    archive=root/f'finite-mass-qed-fixed-pole-cold-scalar-band-{DATE}.json'
    paths=[root/f'early-pole-cell-endpoint-fine{level}-x05-state-{DATE}.json.gz' for level in ('','-d2','-d4')]
    meshes=[tuple(root/f'fine{label}-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2))
            for label in ('','-d2')]+[d4]
    cards=[pair(a,b,archive,c) for a,b,c in zip(meshes,paths,(1,2,4),strict=True)]
    if [v['panels'] for v in cards]!=[27,49,108]:
        raise ArithmeticError('The archived support and all independently determined panel edges are required')
    z=np.array([v['number_energy_and_excess'] for v in cards])
    n=np.array([v['number_neutral'] for v in cards])
    if np.any(z<=0) or np.any(n<=0):
        raise ArithmeticError('Need positive gross terms and occupation roots')
    samples=[]
    for i,j in ((0,1),(1,2),(0,2)):
        x=float(np.sqrt(n[i]*n[j]))
        a,b=z[[i,j],0,1:].T
        net=(x*b-a)/a
        if n[i]!=n[j] and np.sign(net[0])==np.sign(net[1]):
            raise ArithmeticError('The root midpoint must have the predicted two diagnostic net signs')
        samples.append(dict(levels=[cards[k]['charge_subcells'] for k in (i,j)],
                            common_diagnostic_occupation=x,net_relative_to_inverse=net.tolist()))
    return dict(x=.05,k=.0015,q_band=[1e-6,1.5e-5],probabilities=PROBS,tail_cut_q=CUTS,
                cards=cards,d2_over_d1=(z[1]/z[0]).tolist(),d4_over_d2=(z[2]/z[1]).tolist(),
                d4_over_d1=(z[2]/z[0]).tolist(),diagnostic_root_midpoints=samples,
                raw_uncompressed_sha256=[[unpacked_hash(p) for p in v] for v in meshes],
                original_archive_sha256=unpacked_hash(archive),
                scope='Three evolved charge-grid states with fixed scalar-shell geometry and self-matched poles; numerical comparison on a retained band only, not spatial convergence or a physical source.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('d4',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
