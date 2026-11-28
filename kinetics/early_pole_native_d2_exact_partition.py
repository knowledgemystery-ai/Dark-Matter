"""Complete the matched D2 five-band partition with two exact native cuts."""
import argparse
import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_endpoint_cold_targeted_patch_join import merge
from early_pole_endpoint_continuum_cold_band_pilot import run as stock
from early_pole_endpoint_continuum_cold_band_receipt import read,reduce
from early_pole_endpoint_fine_d4_cold_continuum_ladder import DATE,unpacked_hash
from early_pole_native_cold_lowest_exact_aperture import TOPOLOGY,contract,positive_affine
from early_pole_native_d4_exact_partition import EDGES,totals

SPLITS = [[4.851225875113678e-6,5e-6,5.400633696589143e-6],
          [7.990475682151432e-6,8e-6,8.489575496766237e-6]]
IDS = [17,18,27,28]
TOL = 6e-5


def prepare(root,work,m,compute=False):
    low=root/f'cold-lowest-exact-aperture-patches-{DATE}'
    s0,original,card,C,n=contract(root,low,2,m)
    S=json.loads((root/f'early-pole-native-cold-lowest-exact-aperture-{DATE}.json').read_bytes())
    if n!=S['held_occupation'] or unpacked_hash(s0)!=S['charts'][0]['rows'][m-1]['patch_sha256']:
        raise ArithmeticError('The exact D2 low-cut basis must be preserved')
    p=work/f'd2-middle-exact-patch-{m}.jsonl'
    q=stock(card,m,p,topology_only=True,splits=SPLITS,prior_splits=[TOPOLOGY[2][1]],patch_only=True)
    if q['recalculated_panel_indices']!=IDS or q['total']!=52:
        raise ArithmeticError('Need exact new geometry and only four changed radial children')
    if compute:
        stock(card,m,p,splits=SPLITS,prior_splits=[TOPOLOGY[2][1]],patch_only=True)
    return s0,original,C,p,q,n


def finish(root,work):
    arrays=[]; overview=[]
    with TemporaryDirectory() as d:
        for m in (1,2):
            s0,original,C,p,q,n=prepare(root,work,m)
            A=Path(d)/f'start-{m}.jsonl'
            F=Path(d)/f'final-{m}.jsonl'
            merge(original,s0,C,A)
            V=read(p)
            if [r['panel'] for r in V]!=IDS:
                raise ArithmeticError('Four complete D2 children required')
            parents={tuple(v['bounds']):np.array(v['coeff']) for v in read(A)}
            diff=[]
            for j,(l,mid,h) in enumerate(SPLITS):
                x,y=[np.array(z['coeff']) for z in V[2*j:2*j+2]]
                positive_affine(x,q['scalar_mass'])
                positive_affine(y,q['scalar_mass'])
                diff.append(float(np.max(abs((x+y)/parents[(l,h)]-1.))))
            if max(diff)>TOL:
                raise ArithmeticError('Individual D2 parent components do not close')
            merge(A,p,q,F)
            W,b,_=reduce(F)
            W0,_,_=reduce(original)
            closure=float(np.max(abs(W/W0-1.)))
            if closure>TOL:
                raise ArithmeticError('Full D2 combined patch does not close')
            all_rows=read(F)
            groups=[[v for v in all_rows if lo<=v['bounds'][0] and v['bounds'][1]<=hi]
                    for lo,hi in zip(EDGES[:-1],EDGES[1:],strict=True)]
            if [v for g in groups for v in g]!=all_rows or any(
                (g[0]['bounds'][0],g[-1]['bounds'][1])!=(lo,hi)
                for g,lo,hi in zip(groups,EDGES[:-1],EDGES[1:],strict=True)):
                raise ArithmeticError('No D2 original support gaps or overlaps permitted')
            Z=[np.array([r['coeff'] for r in g]) for g in groups]
            rows=[n*t[0,2]-t[0,1] for t in map(totals,Z)]
            if abs(fsum(rows)-(n*W[0,2]-W[0,1]))>8e-29:
                raise ArithmeticError('Held combined D2 partition must reclose')
            arrays.append(Z)
            overview.append(dict(mesh=q['mesh'],counts=list(map(len,Z)),parent_componentwise_reclosure=diff,
                                 full_patched_moment_change_relative=closure,patch_sha256=unpacked_hash(p),
                                 signed_net_by_band=rows))
    bands=[]
    for a,b in zip(*arrays,strict=True):
        err=float(np.max(abs(a/b-1.)))
        if err>TOL:
            raise ArithmeticError('D2 exact subpartition paired mesh exceeds the committed gate')
        lo,hi=np.minimum(a,b)[:,0,:],np.maximum(a,b)[:,0,:]
        net=[fsum(n*lo[:,2]-hi[:,1]),fsum(n*hi[:,2]-lo[:,1])]
        neutral=[fsum(lo[:,1])/fsum(hi[:,2]),fsum(hi[:,1])/fsum(lo[:,2])]
        bands.append(dict(net_box=net,neutral_box=neutral,sign=int(net[0]>0)-int(net[1]<0),
                          max_relative_paired_component_difference=err))
    return dict(edges=list(EDGES),held_occupation=n,bands=bands,mesh_runs=overview,
                scope='Identical five-band geometry as held D4. Neither a physical photon population nor a rigorous quadrature error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('scratch',type=Path)
    p.add_argument('mode',choices=('plan','mesh1','mesh2','finish'))
    a=p.parse_args();a.scratch.mkdir(parents=True,exist_ok=True)
    if a.mode=='finish':
        print(json.dumps(finish(a.root,a.scratch),indent=2,allow_nan=False))
    else:
        for m in ((1,2) if a.mode=='plan' else (int(a.mode[-1]),)):
            v=prepare(a.root,a.scratch,m,a.mode!='plan')[-2]
            print(json.dumps(dict(mesh=v['mesh'],ids=v['recalculated_panel_indices'],total=v['total'])),flush=True)
