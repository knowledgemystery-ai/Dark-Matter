"""Use held centered energy to improve rigorous partial cold-band bounds."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import ledger
from early_pole_native_cold_mean_from_centered_energy import conditional_q, run as verify
from early_pole_native_cold_cross_factors import DATE


def enclosed_cdf(l,r,cut,mu,m):
    if cut<=l:
        return np.zeros((len(mu),2))
    if cut>=r:
        return np.ones((len(mu),2))
    q=np.array([l,cut,r])
    a,c,b=q*q/(np.hypot(q,m)+m)
    pad=3e-13*b
    lo=(c-(mu+pad))/(c-a)
    hi=(b-(mu-pad))/(b-c)
    return np.clip(np.c_[lo,hi],0.,1.)


def interval(rows,a,b):
    table=[]
    for v in rows:
        l,r=v['bounds']
        A,E,U=np.array(v['coeff'])
        m=v['scalar_mass']
        if np.max(abs((m*A+U)/E-1.))>2e-13 or min(A)<=0.:
            raise ArithmeticError('Independent energy and centered moment must close')
        mu=U/A
        conditional_q(np.array([[l,r]]*len(A)),A,U,m)
        lo,hi=[enclosed_cdf(l,r,x,mu,m) for x in (a,b)]
        table.append(A[:,None]*np.clip(hi-lo[:,::-1],0.,1.))
    x=np.array(table)
    return np.array([[fsum(x[:,j,i]) for i in (0,1)] for j in range(3)])


def run(root,prior):
    data=prior.read_bytes()
    S=root/f'early-pole-native-cold-shape-attribution-{DATE}.json'
    if verify(root,S)!=json.loads(data):
        raise ArithmeticError('Exact moment calibration result required')
    L=json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    q=[3e-6,5e-6]
    n=next(x['common_diagnostic_occupation'] for x in L['diagnostic_root_midpoints'] if x['levels']==[2,4])
    out=[]
    for level,i in ((2,1),(4,2)):
        V=[]
        for m in (1,2):
            p=root/f'fine-d{level}-endpoint-continuum-cold-{m}.jsonl.gz'
            if unpacked_hash(p)!=L['raw_uncompressed_sha256'][i][m-1]:
                raise ArithmeticError('Exact complete D2/D4 quadrature required')
            R=ledger(p)
            x=interval(R,*q)
            z=n*x[2]-x[1,::-1]
            V.append(dict(gross_forward_inverse_slope=x.tolist(),conditional_net=z.tolist(),
                          sign=int(z[0]>0.)-int(z[1]<0.)))
        if V[0]['sign']!=V[1]['sign']:
            raise ArithmeticError('Observed sign is not stable on both admitted meshes')
        out.append(dict(level=level,meshes=V))
    return dict(q=q,diagnostic_occupation=n,rows=out,mean_receipt_sha256=hashlib.sha256(data).hexdigest(),
                scope='Positive fractions from actual held centered and number panel moments on selected support. Separate from the prespecified D8 comparison, not its result or a physical sign.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','prior'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
