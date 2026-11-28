"""Patch only the unresolved native D2/D4 cut at q=3e-6."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_endpoint_cold_targeted_patch_join import merge
from early_pole_endpoint_continuum_cold_band_pilot import run as stock
from early_pole_endpoint_continuum_cold_band_receipt import read, reduce
from early_pole_endpoint_fine_d4_cold_continuum_ladder import DATE, unpacked_hash

Q = [1e-6,3e-6]
TOPOLOGY = {2:(7,[2.459599833865599e-6,3e-6,3.1356591257655634e-6]),
            4:(18,[2.97457405133568e-6,3e-6,3.1916537495062577e-6])}
TOL = 6e-5


def contract(root,w,k,m):
    p = root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json'
    C = json.loads(p.read_bytes())
    i = {2:1,4:2}[k]
    original = root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz'
    state = root/f'early-pole-cell-endpoint-fine-d{k}-x05-state-{DATE}.json.gz'
    if (unpacked_hash(original)!=C['raw_uncompressed_sha256'][i][m-1]
            or unpacked_hash(state)!=C['cards'][i]['state_sha256_uncompressed']):
        raise ArithmeticError('Exact original D2 or D4 archive required')
    a = json.loads(gzip.decompress(state.read_bytes()))
    path = w/f'lowest-d{k}-mesh{m}.jsonl'
    _, edges, meta = reduce(original)
    t = stock(a,m,path,topology_only=True)
    if t['edges']!=[edges[0][0],*[b for _,b in edges]] or any(t[key]!=v for key,v in meta.items()):
        raise ArithmeticError('Original exact stock geometry must replay')
    j,cut = TOPOLOGY[k]
    if t['edges'][j:j+2]!=[cut[0],cut[2]]:
        raise ArithmeticError('Original host of the single cut changed')
    g = stock(a,m,path,topology_only=True,splits=[cut],patch_only=True)
    if g['recalculated_panel_indices']!=[j,j+1] or g['edges'][j:j+3]!=cut:
        raise ArithmeticError('Only the two exact low-cut children are eligible')
    n = next(x['common_diagnostic_occupation'] for x in C['diagnostic_root_midpoints'] if x['levels']==[2,4])
    return path,original,a,g,n


def run_one(root,w,k,m):
    w.mkdir(parents=True,exist_ok=True)
    path,_,a,g,_ = contract(root,w,k,m)
    print(json.dumps(dict(level=k,mesh=g['mesh'],children=g['recalculated_panel_indices'])),flush=True)
    stock(a,m,path,splits=[TOPOLOGY[k][1]],patch_only=True)


def positive_affine(C,mass):
    F,B,D = C.T
    if not np.isfinite(C).all() or np.min(C)<=0. or max(abs((B+D)/F-1.))>5e-13:
        raise ArithmeticError('Independent affine entries must close')
    if max(abs((mass*C[0]+C[2])/C[1]-1.))>3e-13:
        raise ArithmeticError('Independent centered energy must close')


def finish(root,w):
    old = json.loads((root/f'early-pole-native-cold-patch-adjacent-corridors-{DATE}.json').read_bytes())
    prior = next(x for x in old['rows'] if x['q']==Q)
    result = []
    with TemporaryDirectory() as d:
        for k in (2,4):
            rows=[]; M=[]; support=None
            for m in (1,2):
                path,original,a,g,n = contract(root,w,k,m)
                j,cut = TOPOLOGY[k]
                received = read(path)
                if [r['panel'] for r in received]!=[j,j+1]:
                    raise ArithmeticError('Both complete target children required')
                C = np.array([r['coeff'] for r in received])
                for x in C:
                    positive_affine(x,g['scalar_mass'])
                oldparent = np.array(read(original)[j]['coeff'])
                rc = float(np.max(abs(C.sum(axis=0)/oldparent-1.)))
                if rc>TOL:
                    raise ArithmeticError('Local parent correction exceeds prospective componentwise tolerance')
                out = Path(d)/f'd{k}m{m}.jsonl'
                p = merge(original,path,g,out)
                keep = [r for r in read(out) if Q[0]<=r['bounds'][0] and r['bounds'][1]<=Q[1]]
                b = [r['bounds'] for r in keep]
                if b[0][0]!=Q[0] or b[-1][1]!=Q[1] or (support is not None and b!=support):
                    raise ArithmeticError('Exact common q cut and matching mesh topology required')
                support=b
                X = np.array([r['coeff'] for r in keep])
                M.append(X)
                rows.append(dict(mesh=g['mesh'],patch_sha256=unpacked_hash(path),selected_panels=len(keep),
                                 parent_componentwise_reclosure=rc,total_joined_panels=p['result_panels']))
            r0,r1 = M
            delta = float(np.max(abs(r0/r1-1.)))
            if delta>TOL:
                raise ArithmeticError('Paired lowest-cut result is underresolved')
            lo,hi=np.minimum(r0,r1)[:,0,:],np.maximum(r0,r1)[:,0,:]
            net = [fsum(n*lo[:,2]-hi[:,1]),fsum(n*hi[:,2]-lo[:,1])]
            low,high = prior['held_nets_D2_and_D4'][{2:0,4:1}[k]]
            ramp = 5e-12*(high-low)
            if net[0]<low-ramp or net[1]>high+ramp:
                raise ArithmeticError('Observed exact clip must agree with prior loose bounds')
            neutral = [fsum(lo[:,1])/fsum(hi[:,2]),fsum(hi[:,1])/fsum(lo[:,2])]
            result.append(dict(level=k,rows=rows,net_box=net,neutral_box=neutral,
                total_forward_and_inverse_by_mesh=[[fsum(x[:,0,i]) for i in (0,1)] for x in M],
                sign=int(net[0]>0.)-int(net[1]<0.),maximum_relative_paired_difference=delta))
    return dict(q=Q,held_occupation=n,charts=result,
                scope='Fixed-state held diagnostic on the exact previously unresolved cut. Observed two-mesh box, not continuum error or physical photon occupation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('work',type=Path)
    p.add_argument('mode',choices=('plan','d2m1','d2m2','d4m1','d4m2','finish'))
    a=p.parse_args()
    if a.mode=='finish':
        print(json.dumps(finish(a.root,a.work),indent=2,allow_nan=False))
    elif a.mode=='plan':
        for k in (2,4):
            for m in (1,2):
                _,_,_,c,_=contract(a.root,a.work,k,m)
                print(k,m,c['recalculated_panel_indices'])
    else:
        run_one(a.root,a.work,int(a.mode[1]),int(a.mode[3]))
