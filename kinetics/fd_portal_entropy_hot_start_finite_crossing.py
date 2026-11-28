"""Locate a held finite-gain crossing from the same hot-side window origin."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff,full
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_entropy_short_tangent_windows import DATE,gain_difference
from sm_dark_entropy_exchange_reference import StandardModel


def dense(P,path,times,C,mesh,eps):
    sol=[]
    for a,b in zip(times[:-1],times[1:],strict=True):
        x=np.linspace(a,b,mesh)
        z=[]
        for t in x:
            p=path(t)
            f=full(P,p[1:4])
            z.append(diff(f,eps).sum(axis=0)-np.diag(f(np.ones(3)).sum(axis=0)))
        it=CubicSpline(x,z)
        def rhs(t,y):
            return (it(t)@y.reshape(3,2)).ravel()
        s=solve_ivp(rhs,(a,b),C.ravel(),method='DOP853',rtol=3e-13,atol=2e-15,max_step=(b-a)/3,dense_output=True)
        if not s.success:
            raise ArithmeticError('Local tangent must solve each prescribed short segment')
        C=s.y[:,-1].reshape(3,2)
        sol.append(s.sol)
    def frame(t):
        if not times[0]<=t<=times[-1]:
            raise ArithmeticError('Outside held clock window')
        i=min(np.searchsorted(times,t,side='right')-1,len(sol)-1)
        return sol[i](t).reshape(3,2)
    return frame


def observable(frame,metric,origin):
    def form(t):
        C=frame(t)
        M=C.T@metric(t)@C
        return (M+M.T)/2
    X=form(origin)
    return lambda t: float(gain_difference(X,form(t))[-1])


def locate(g,a,b):
    if not g(a)>0.>g(b):
        raise ArithmeticError('A positive-to-negative finite-gain bracket is required')
    return brentq(g,a,b,xtol=2e-14)


def run(eos,root):
    names=('direct-expansion-096-120','entropy-short-tangent-windows','entropy-short-nonlinear-witness')
    raw=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    L,T,N=map(json.loads,raw)
    ids=[hashlib.sha256(s).hexdigest() for s in raw]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (T['input_sha256'][0]!=ids[0] or N['input_sha256']!=ids[:2]
        or any(s['table_input_sha256']!=h for s in (L,T,N))):
        raise ArithmeticError('Exact short tangent and independent nonlinear witness required')
    ini,hot=T['windows'][0]['x']
    cool,last=T['windows'][2]['x']
    knot=T['split_at_knot']
    seed=next(p for p in next(r['points'] for r in L['rows'] if r['operator']=='combined') if p[0]==ini)
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    path=local_path(P,seed)
    C=np.asarray(T['rebased_fractional_columns'])
    rows=[]
    for n,eps in ((4,4e-5),(6,2e-5)):
        frame=dense(P,path,[ini,hot,knot,cool,last],C,n,eps)
        metric=lambda t:matrix(P,path(t))[0]
        hotg=observable(frame,metric,hot)
        inig=observable(frame,metric,ini)
        if abs(hotg(cool)-T['windows'][1]['ordered_gain_minus_one'][-1])>2e-11:
            raise ArithmeticError('Replay the original hot-to-cool direct tangent ratio')
        z=locate(hotg,cool,last)
        a,b=z-1e-8,z+1e-8
        if not cool<a<z<b<last or not hotg(a)>0.>hotg(b) or inig(z)<=0.:
            raise ArithmeticError('Need a local resolving bracket, distinct from the .99 origin')
        rows.append(dict(x_crossing=z,T_at_crossing=path(z)[1],sign_bracket_both_sides=[a,b],
                         hot_origin_gain_minus_one_on_sign_bracket=[hotg(a),hotg(b)],
                         hot_origin_gain_minus_one_at_original_last_clock=hotg(last),
                         x99_origin_gain_minus_one_at_same_crossing=inig(z),
                         hot_origin_gain_minus_one_at_prior_cool_clock=hotg(cool)))
    if abs(rows[0]['x_crossing']-rows[1]['x_crossing'])>1e-9:
        raise ArithmeticError('Independent input meshes did not reproduce the finite clock crossing')
    return dict(hot_origin=hot,earlier_origin=ini,cover=[cool,last],rows=rows,table_input_sha256=h,input_sha256=ids,
                scope='At least one numerical finite gain crossing in a sampled held bracket; no uniqueness, physical threshold or continuous stability claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
