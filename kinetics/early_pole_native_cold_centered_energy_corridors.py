"""Combine archived centered energy with the five held neutral-corridor checks."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_neutral_occupation_corridors import corridor, hull, overlap, run as support_only
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_subband_from_centered_energy import interval,ledger,run as verify
from early_pole_native_cold_cross_factors import DATE


def roots(x):
    B,D=np.asarray(x)[1:]
    if min(B[0],D[0])<=0:
        raise ArithmeticError('Positive resolved inverse and slope intervals required')
    return [B[0]/D[1],B[1]/D[0]]


def run(root,prior):
    raw=prior.read_bytes()
    Z=json.loads(raw)
    if verify(root,root/f'early-pole-native-cold-mean-from-centered-energy-{DATE}.json')!=Z:
        raise ArithmeticError('Exact previous central moment receipts required')
    a_raw=(root/f'early-pole-cold-neutral-occupation-corridors-{DATE}.json').read_bytes()
    a=json.loads(a_raw)
    if support_only(root)!=a:
        raise ArithmeticError('Exact archived support-only reference required')
    J=json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    datasets=[]
    for level,i in ((2,1),(4,2)):
        p=[root/f'fine-d{level}-endpoint-continuum-cold-{m}.jsonl.gz' for m in (1,2)]
        if [unpacked_hash(s) for s in p]!=J['raw_uncompressed_sha256'][i]:
            raise ArithmeticError('Only the complete two accepted ledgers are eligible')
        datasets.append([ledger(s) for s in p])
    n=Z['diagnostic_occupation']
    out=[]
    for old in a['rows']:
        q=old['q']
        data=[[interval(r,*q) for r in v] for v in datasets]
        rr=[hull([roots(z) for z in x]) for x in data]
        u=[]
        for i,V in enumerate(data):
            y=hull([(n*t[2]-t[1,::-1]).tolist() for t in V])
            x=rr[i]
            old_r=old[f'neutral_D{(2,4)[i]}_union_hull']
            if x[0]<old_r[0]*(1-2e-12) or x[1]>old_r[1]*(1+2e-12):
                raise ArithmeticError('New moment bounds contradict support-only bounds')
            u.append(dict(level=(2,4)[i],neutral_by_mesh=[roots(z) for z in V],
                          refined_neutral_hull=x,previous_support_only_hull=old_r,
                          held_net_hull=y,held_sign=int(y[0]>0)-int(y[1]<0)))
        out.append(dict(q=q,charts=u,**corridor(*rr)))
    q=[3e-6,5e-6]
    C=next(v for v in out if v['q']==q)
    if [v['held_sign'] for v in C['charts']]!=[v['meshes'][1]['sign'] for v in Z['rows']]:
        raise ArithmeticError('New all-band treatment must preserve the certified central signs')
    flag='d2_loss_d4_gain_sufficient'
    next_band=next(v for v in out if v['q']==[5e-6,8e-6])
    together=overlap(overlap(C[flag],next_band[flag]),out[0][flag])
    return dict(rows=out,diagnostic_occupation=n,full_central_and_next_simultaneous_sufficient=together,
                original_support_reference_sha256=hashlib.sha256(a_raw).hexdigest(),original_centered_subband_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Neutral corridors with positive panel centered moments, both admitted D2/D4 meshes; no D8 rate, physical occupation or modified D8 protocol.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','prior'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
