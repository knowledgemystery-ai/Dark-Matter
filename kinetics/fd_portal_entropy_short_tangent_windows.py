"""Transport the archived fixed-initial plane across a short PCHIP knot window."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff,full
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_heat_capacity_slope_threshold import DATE
from sm_dark_entropy_exchange_reference import StandardModel


def gain_difference(X,Y):
    l=eigh((Y-X+(Y-X).T)/2,X,eigvals_only=True)
    if min(l)<=-1.:
        raise ArithmeticError('Endpoint form must remain positive definite')
    return l/(np.sqrt(1+l)+1.)


def propagate(P,path,times,columns,mesh,ds,skip=()):
    forms=[]
    C=[columns]
    for a,b in zip(times[:-1],times[1:],strict=True):
        ts=np.linspace(a,b,mesh)
        L=[]
        for x in ts:
            p=path(x)
            f=full(P,p[1:4])
            L.append(diff(f,ds).sum(axis=0)-np.diag(f(np.ones(3)).sum(axis=0)))
        spline=CubicSpline(ts,L)
        def rhs(x,y):
            return (spline(x)@y.reshape(3,2)).ravel()
        s=solve_ivp(rhs,(a,b),C[-1].ravel(),method='DOP853',rtol=3e-13,atol=2e-15,max_step=(b-a)/3)
        if not s.success:
            raise ArithmeticError('Direct local tangent propagation failed')
        C.append(s.y[:,-1].reshape(3,2))
    for x,c in zip(times,C,strict=True):
        if x in skip:
            forms.append(None)
            continue
        G=matrix(P,path(x))[0]
        R=c.T@G@c
        forms.append((R+R.T)/2)
    return forms


def run(eos,root):
    names=('direct-expansion-096-120','transported-entropy-plane','entropy-pchip-knot-crossing','entropy-finite-windows')
    raws=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    B,T,K,F=map(json.loads,raws)
    sh=[hashlib.sha256(r).hexdigest() for r in raws]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (any(z['table_input_sha256']!=h for z in (B,T,K,F)) or T['input_sha256'][1]!=sh[0]
        or K['input_sha256'][0]!=sh[0] or F['input_sha256'][1]!=sh[0] or F['input_sha256'][-1]!=sh[2]):
        raise ArithmeticError('Exact archive, carried plane and completed negative control required')
    hist=next(r['points'] for r in B['rows'] if r['operator']=='combined')
    seed=next(a for a in hist if a[0]==.99)
    cols=np.array(next(z for z in T['rows'] if z['x']==.99)['endpoint_fractional_plane_columns'])
    Q,R=np.linalg.qr(cols,mode='reduced')
    if min(abs(R.diagonal()))<=1e-7:
        raise ArithmeticError('Exact archived initial-total plane must have rank two')
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    p=local_path(P,seed)
    x=K['crossing_x']
    if abs(p(x)[1]-K['visible_temperature_knot'])>3e-14:
        raise ArithmeticError('Stored knot must replay on actual held path')
    times=[.99,x-2.5e-5,x,x+2.5e-5,x+8e-5]
    a,b=[propagate(P,p,times,Q,n,d,skip=(x,)) for n,d in ((4,4e-5),(6,2e-5))]
    inds=((0,1),(1,3),(3,4),(0,4))
    old=F['rows']['pchip']['windows']
    out=[]
    for (i,j),base in zip(inds,old,strict=True):
        if base['interval']!=[times[i],times[j]]:
            raise ArithmeticError('Must use the same four finite windows')
        u,v=[gain_difference(s[i],s[j]) for s in (a,b)]
        err=float(np.max(abs(u-v)))
        if err>1e-8 or abs(v[-1])<8*err:
            raise ArithmeticError('Local tangent gain sign not resolved by both tested meshes')
        out.append(dict(x=[times[i],times[j]],ordered_gain_minus_one=v.tolist(),mesh_after_max_absolute_difference=err,
                        long_history_estimate_minus_one=base['largest_and_smallest'][0]-1.,
                        direct_long_history_ode_difference=base['ode_gain_difference']))
    return dict(windows=out,split_at_knot=x,initial_fractional_columns=cols.tolist(),rebased_fractional_columns=Q.tolist(),
                archived_column_condition_number=float(np.linalg.cond(cols)),table_input_sha256=h,
                input_sha256=sh,scope='Short original PCHIP tangent initialized by the archived numerical transported plane, not independent nonlinear evolution or a cubic estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
