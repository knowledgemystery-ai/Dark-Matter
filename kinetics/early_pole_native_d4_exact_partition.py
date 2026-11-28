"""Rejoin all frozen cuts into one exact D4 partition and close its net budget."""
import argparse
import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_endpoint_cold_targeted_patch_join import merge
from early_pole_endpoint_continuum_cold_band_pilot import run as stock
from early_pole_endpoint_continuum_cold_band_receipt import read,reduce
from early_pole_native_cold_lowest_exact_aperture import DATE,TOPOLOGY,contract,finish as lowest
from early_pole_native_d4_clipped_next_band import IDS,SPLITS,prepare,unpacked_hash

EDGES = (1e-6,3e-6,5e-6,8e-6,1e-5,1.5e-5)
COUNTS = (19,15,24,17,36)


def totals(a):
    return np.array([[fsum(a[:,i,j]) for j in range(3)] for i in range(3)])


def selected(root,m,work,n,low_source):
    d4lo,orig,card,pa,n1=contract(root,low_source,4,m)
    midbase=root/f'd4-clipped-next-band-patches-{DATE}'
    _,orig2,d4mi,pb=prepare(root,m,midbase)
    if orig!=orig2 or n!=n1:
        raise ArithmeticError('Both recipes must share one exact original chart and held occupation')
    prob=stock(card,m,None,topology_only=True,splits=SPLITS,prior_splits=[TOPOLOGY[4][1]],patch_only=True)
    ids=prob['recalculated_panel_indices']
    expected=json.loads((root/f"early-pole-native-d4-clipped-next-band-{DATE}.json").read_bytes())
    if unpacked_hash(d4mi)!=expected["rows"][m-1]["stock_split_patch_sha256"]:
        raise ArithmeticError("Prior exact D4 cut must match the accepted archive bytes")
    old=read(d4mi)
    if [v['panel'] for v in old]!=IDS or any(v['bounds']!=pb['edges'][k:k+2] for v,k in zip(old,IDS,strict=True)):
        raise ArithmeticError('Preserve all original mid child geometries')
    if ids!=[k+1 for k in IDS] or any(v['bounds']!=prob['edges'][k:k+2] for v,k in zip(old,ids,strict=True)):
        raise ArithmeticError('Only explicit panel index translation is permitted')
    A,B,F=[work/f'{k}-{m}.jsonl' for k in ('first','second-patch','full')]
    merge(orig,d4lo,pa,A)
    B.write_text(''.join(json.dumps(dict(v,panel=k,total=prob['total']),allow_nan=False)+'\n'
                         for v,k in zip(old,ids,strict=True)))
    merge(A,B,prob,F)
    Q,bounds,_=reduce(F)
    original,_,_=reduce(orig)
    err=float(np.max(abs(Q/original-1.)))
    if err>6e-5 or len(bounds)!=sum(COUNTS):
        raise ArithmeticError('Exact patched total must reclose against original')
    rows=read(F)
    masks=[[v for v in rows if l<=v['bounds'][0] and v['bounds'][1]<=r] for l,r in zip(EDGES[:-1],EDGES[1:],strict=True)]
    if (tuple(map(len,masks))!=COUNTS or [v for b in masks for v in b]!=rows
           or any((b[0]['bounds'][0],b[-1]['bounds'][1])!=(l,r) for b,l,r in zip(masks,EDGES[:-1],EDGES[1:],strict=True))):
        raise ArithmeticError('All exact panels must be counted once')
    Z=[np.array([r['coeff'] for r in v]) for v in masks]
    net=np.array([n*A[0,2]-A[0,1] for A in map(totals,Z)])
    closure=abs(fsum(net)-(n*Q[0,2]-Q[0,1]))
    if closure>3e-29:
        raise ArithmeticError('D4 partition does not close the complete held net')
    return Z,dict(full_patched_moment_change_relative=err,held_net_partition_closure=closure,
                  signed_net_by_band=net.tolist(),central_share_of_positive_total=net[1]/fsum(net[net>0]),
                  opposing_upper_fraction_of_positive_total=-fsum(net[net<0])/fsum(net[net>0]))


def run(root):
    p=root/f'early-pole-native-cold-lowest-exact-aperture-{DATE}.json'
    inp=root/f'cold-lowest-exact-aperture-patches-{DATE}'
    L=json.loads(p.read_bytes())
    if lowest(root,inp)!=L:
        raise ArithmeticError('Lowest exact two-chart receipt must reproduce')
    n=L['held_occupation']
    with TemporaryDirectory() as d:
        pair=[selected(root,m,Path(d),n,inp) for m in (1,2)]
    bands=[]
    for a,b in zip(pair[0][0],pair[1][0],strict=True):
        err=float(np.max(abs(a/b-1.)))
        if err>6e-5:
            raise ArithmeticError('Unresolved final exact support')
        lo,hi=np.minimum(a,b)[:,0,:],np.maximum(a,b)[:,0,:]
        z=[fsum(n*lo[:,2]-hi[:,1]),fsum(n*hi[:,2]-lo[:,1])]
        B=[fsum(lo[:,1])/fsum(hi[:,2]),fsum(hi[:,1])/fsum(lo[:,2])]
        bands.append(dict(net_box=z,neutral_box=B,sign=int(z[0]>0)-int(z[1]<0),
                          maximum_relative_paired_difference=err))
    if ([z['sign'] for z in bands]!=[1,1,1,-1,-1]
          or any(a['neutral_box'][1]>=b['neutral_box'][0] for a,b in zip(bands[:-1],bands[1:],strict=True))):
        raise ArithmeticError('Need distinct ordered held band neutrals on exact partition')
    return dict(edges=list(EDGES),counts=list(COUNTS),held_occupation=n,bands=bands,by_mesh=[v[1] for v in pair],
                scope='Complete exact five-band partition of one held D4 numerical diagnostic; no pointwise ordering, spatial convergence, rigorous continuum error, or physical occupation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
