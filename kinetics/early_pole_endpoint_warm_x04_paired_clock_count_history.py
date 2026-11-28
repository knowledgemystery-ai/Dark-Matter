"""Replay count cancellation at all common outputs of three actual G9/mixed300 clocks."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read


def run(g9_report,mixed_report,g9,mixed,comparison):
    A,B=[read(x) for x in (g9_report,mixed_report)]
    raw=[p.read_bytes() for p in (g9,mixed,comparison)]
    G,M,C=map(json.loads,raw)
    h=[hashlib.sha256(v).hexdigest() for v in raw]
    X=[]
    for (r,rh),x,q in zip((A,B),(G,M),C['rows'],strict=True):
        if (r['coarse_steps']!=1024 or x['steps']!=[1024,2048,4096] or q['sha256'][1]!=rh
            or x['input_hashes'][-1]!=h[2]):
            raise ArithmeticError('Need exact report bytes and three-clock outputs')
        X.append(np.array([[[z['x'],z['scalar_yield']] for z in cc] for cc in (*r['checkpoint_trajectories'],x['checks'])]))
    a,b=X
    if (a.shape!=(3,8,2) or b.shape!=a.shape or np.any(a[:,:,1]<=0.) or np.any(b[:,:,1]<=0.) or
         any(not np.array_equal(u,a[0,:,0]) for u in np.r_[a[:,:,0],b[:,:,0]])):
        raise ArithmeticError('All count samples must share the same clock and be positive')
    U,V=a[:,:,1],b[:,:,1]
    pair=V/U-1.
    per=np.stack([U[:-1]/U[1:]-1.,V[:-1]/V[1:]-1.],axis=1)
    dp=np.diff(pair,axis=0)
    rows=[dict(x=x,matched_contrasts=pair[:,i].tolist(),within_chart_shifts=per[:,:,i].tolist(),
               change_in_contrast=dp[:,i].tolist(),
               absolute_suppression_factor=(abs(dp[:,i])/np.max(abs(per[:,:,i]),axis=1)).tolist())
          for i,x in enumerate(a[0,:,0])]
    expected=np.array([G['previous_over_next_minus_one'],M['previous_over_next_minus_one']])[:,:,0].T
    if np.max(abs(per[:,:,-1]-expected))>2e-14:
        raise ArithmeticError('Last output must match accepted endpoint totals')
    return dict(steps=[1024,2048,4096],fields=['G9','mixed300'],rows=rows,input_sha256=[A[1],B[1],*h],
                scope='Three actual matched clocks on each fixed-start history, empirical local common-mode cancellation only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('g9_report','mixed_report','g9','mixed','comparison'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
