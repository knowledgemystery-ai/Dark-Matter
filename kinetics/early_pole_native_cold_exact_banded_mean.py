"""Separate between-band redistribution and within-band mean-q movement."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_endpoint_cold_targeted_patch_join import merge
from early_pole_endpoint_continuum_cold_band_receipt import read
from early_pole_native_cold_bose_reweight import address,strip
from early_pole_native_cold_cross_factors import CellChart, SPECS, load, read as read_card
from early_pole_native_cold_exact_paired_bands import run as verify
from early_pole_native_cold_mean_from_centered_energy import conditional_q
from early_pole_native_d2_exact_partition import DATE,prepare
from early_pole_native_d4_exact_partition import EDGES,selected


def scalar_sum(weights,boxes):
    return np.array([fsum(weights*np.where(weights>=0.,boxes[:,i],boxes[:,1-i])) for i in (0,1)])


def moments(root,path,chart,raw,column):
    rows=read(path)
    q=np.array([v['bounds'] for v in rows])
    a=np.array([v['coeff'] for v in rows])
    m=rows[0]['scalar_mass']
    if raw:
        z,_=read_card(root/f'early-pole-cell-endpoint-fine-d{chart}-x05-state-{DATE}.json.gz')
        d=CellChart(load('fine')[0],SPECS['fine'],chart,z['reference_x'])
        f=np.asarray(z['endpoint_occupation'][:d.ns])[address(d.ds.ed,q)]
        a=np.stack([strip(a[:,i],f) for i in range(3)],axis=1)
    if np.max(abs((m*a[:,0]+a[:,2])/a[:,1]-1.))>2e-13:
        raise ArithmeticError('Independent energy and centered moment must match')
    w=a[:,0,column]
    u=conditional_q(q,w,a[:,2,column],m)
    masked=[(q[:,0]>=l)&(q[:,1]<=h) for l,h in zip(EDGES[:-1],EDGES[1:],strict=True)]
    if not np.all(np.sum(masked,axis=0)==1):
        raise ArithmeticError('All panels appear in exactly one shared aperture')
    p=[fsum(w[b]) for b in masked]
    g=np.array([scalar_sum(w[b],u[b])/v for b,v in zip(masked,p,strict=True)])
    return g,np.array(p)/fsum(p)


def partition(m2,m4,p2,p4):
    shift=m4-m2[:,::-1]
    between=scalar_sum(p4-p2,(m2+m4)/2)
    within=scalar_sum((p2+p4)/2,shift)
    direct=scalar_sum(p4,m4)-scalar_sum(p2,m2)[::-1]
    a,b=m2.mean(axis=1),m4.mean(axis=1)
    check=fsum(p4*b-p2*a)-fsum((p4-p2)*(a+b)/2+(p4+p2)*(b-a)/2)
    if abs(check)>2e-20:
        raise ArithmeticError('Symmetric attribution failed algebraic closure')
    if not between[0]>within[1]>0. or within[0]<=0. or direct[0]<=0.:
        raise ArithmeticError('Between-band and within-band ordering must resolve')
    return dict(D2_then_D4_band_probability=[p2.tolist(),p4.tolist()],
                D2_then_D4_band_conditional_mean=[m2.tolist(),m4.tolist()],conditional_D4_minus_D2=shift.tolist(),
                whole_D4_minus_D2_direct=direct.tolist(),
                between_band_weight_redistribution=between.tolist(),within_band_mean_change=within.tolist(),
                algebraic_center_closure=abs(check))


def run(root):
    f=root/f'early-pole-native-cold-exact-paired-bands-{DATE}.json'
    data=f.read_bytes()
    C=json.loads(data)
    if verify(root)!=C:
        raise ArithmeticError('Need exactly matching complete accepted D2 and D4 bands')
    out=[]
    for m in (1,2):
        with TemporaryDirectory() as dd:
            w=Path(dd)
            L,O,Q,B,H,n=prepare(root,root/f'd2-exact-partition-patches-{DATE}',m)
            A,F=w/'d2-low',w/'d2-complete'
            merge(O,L,Q,A)
            merge(A,B,H,F)
            selected(root,m,w,n,root/f'cold-lowest-exact-aperture-patches-{DATE}')
            files=(F,w/f'full-{m}.jsonl')
            rec={}
            for raw in (False,True):
                cards={}
                for col,label in enumerate(('forward','inverse')):
                    (m2,p2),(m4,p4)=[moments(root,path,k,raw,col) for path,k in zip(files,(2,4),strict=True)]
                    z=partition(m2,m4,p2,p4)
                    signs=[int(a>0)-int(b<0) for a,b in z['conditional_D4_minus_D2']]
                    if signs!=([1]*5 if raw else [0,1,1,1,1]):
                        raise ArithmeticError('Interesting raw/weighted difference must resolve on both meshes')
                    cards[label]=z
                rec['raw' if raw else 'weighted']=cards
            out.append(rec)
    return dict(band_edges=list(EDGES),by_mesh=out,exact_paired_receipt_sha256=hashlib.sha256(data).hexdigest(),
                scope='Positive gross measures only, not signed net, physical momentum history or a D8 result. Between/within attribution is an exact symmetric decomposition, not a causal factorization.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
