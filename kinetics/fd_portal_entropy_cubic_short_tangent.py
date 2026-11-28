"""Resolve the four cubic finite windows using its own archived carried plane."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_entropy_short_tangent_windows import gain_difference, propagate
from fd_portal_entropy_spline_reintegration import integrate
from fd_portal_heat_capacity_slope_threshold import DATE


def run(eos,root):
    names=('entropy-spline-reintegration','spline-entropy-endpoint-gain','entropy-finite-windows')
    raw=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    B,A,F=map(json.loads,raw)
    dig=[hashlib.sha256(x).hexdigest() for x in raw]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (any(z['table_input_sha256']!=h for z in (B,A,F)) or A['input_sha256'][0]!=dig[0] or F['input_sha256'][2]!=dig[0]):
        raise ArithmeticError('Exact already re-evolved smooth controls required')
    P=CombinedExpansion(SplineEOS(eos),PositivePortal((17,44,32)),168)
    r,s=[next(z['alternate_reintegrated_state'] for z in B['rows'] if z['x']==x) for x in (.96,.99)]
    old=F['rows']['cubic']['windows']
    t0,t1=old[0]['interval']
    t2=old[2]['interval'][0]; t3=old[2]['interval'][1]
    times=[.96,.97,.98,t0,t1,t2,t3]
    C=np.array(next(x for x in A['rows'] if x['x']==.96)['alternate_end_fractional_columns'])
    Q,R=np.linalg.qr(C,mode='reduced')
    if min(abs(R.diagonal()))<=1e-6:
        raise ArithmeticError('Need independent full-rank plane columns from own background')
    out=[]
    for mesh,df,tol,step in ((4,4e-5,8e-12,.005),(6,2e-5,2e-12,.0025)):
        x=np.unique(np.concatenate([np.linspace(a,b,mesh) for a,b in zip(times[:-1],times[1:],strict=True)]))
        states,_=integrate(P,r,x,tol,step)
        p=dict(zip(x,states,strict=True))
        if max(abs(p[t0]-s))>1e-9:
            raise ArithmeticError('Cubic short initial state must replay prior alternate integration')
        out.append(propagate(P,p.__getitem__,times,Q,mesh,df))
    indices=((3,4),(4,5),(5,6),(3,6))
    rows=[]
    for (i,j),base in zip(indices,old,strict=True):
        if base['interval']!=[times[i],times[j]]:
            raise ArithmeticError('Window endpoints changed')
        low,high=[gain_difference(z[i],z[j]) for z in out]
        err=float(np.max(abs(low-high)))
        if err>2e-9 or abs(high[-1])<8*err:
            raise ArithmeticError('Cubic finite short-window sign did not resolve')
        rows.append(dict(x=base['interval'],ordered_gain_minus_one=high.tolist(),two_setting_defect=err,
                         prior_direct_long_history_ode_difference=base['ode_gain_difference']))
    return dict(windows=rows,start_x=.96,initial_rebased_fractional_columns=Q.tolist(),table_input_sha256=h,input_sha256=dig,
                scope='Cubic continuation on its own held alternate trajectory with an archived transported plane, not a full new evolution or physical contraction theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
