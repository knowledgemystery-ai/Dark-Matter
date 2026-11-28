"""Resolve one-sided local entropy-rate signs at the nearby PCHIP EOS knot."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq, root

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import full
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_local_entropy_metric_rates import measure
from sm_dark_entropy_exchange_reference import StandardModel


def local_path(stage,seed):
    x,T,eta,lt=seed[:4]
    ref=stage.bath([eta,lt])
    scale=np.array([ref['n'],ref['rho'],ref['rho']+stage.eos.state(T)[0]])
    f=full(stage,[T,eta,lt])
    def rhs(s,y):
        a=np.exp(3*(s-x))
        return a*f(y/a).sum(axis=0)
    sol=solve_ivp(rhs,(x,x+.0011),np.ones(3),rtol=3e-13,atol=8e-15,max_step=.00016,
                  method='DOP853',dense_output=True)
    if not sol.success:
        raise ArithmeticError('Original short continuation through knot did not resolve')
    def state(s):
        u,y,z=scale*sol.sol(s)/np.exp(3*(s-x))
        t=brentq(lambda q:stage.eos.state(q)[0]/(z-y)-1.,stage.lo,stage.hi,xtol=5e-15)
        def err(q):
            d=stage.bath(q)
            return np.array([d['n']/u-1.,d['rho']/y-1.])
        d=root(err,[eta,lt],tol=2e-12).x
        if max(abs(err(d)))>2e-12:
            raise ArithmeticError('Archived-state continuation did not invert')
        return [s,t,*d]
    if max(abs(np.array(state(x))-seed[:4]))>1e-13:
        raise ArithmeticError('Demand the exact held original x=.99 state')
    return state


def slopes(E,k):
    idx=int(np.argmin(abs(E.gs.x-k)))
    if idx==0 or idx==len(E.gs.x)-1 or abs(E.gs.x[idx]-k)>1e-13:
        raise ArithmeticError('Need interior exact entropy knot')
    h=k-E.gs.x[idx-1]
    # One-sided limiting values of second derivative from the cubic polynomials.
    g2=[6*E.gs.c[0,idx-1]*h+2*E.gs.c[1,idx-1],2*E.gs.c[1,idx]]
    gs,g1=E.gs(k),E.slope(k)
    return [float(3+(3*g1+d)/(3*gs+g1)) for d in g2]


def run(eos,late,interp):
    raw=[r.read_bytes() for r in (late,interp)]
    B,I=map(json.loads,raw)
    H=[hashlib.sha256(z).hexdigest() for z in raw]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if B['table_input_sha256']!=h or I['table_input_sha256']!=h or I['input_sha256'][1]!=H[0]:
        raise ArithmeticError('Need exact late original branch and spline contrast')
    seed=next(p for p in next(z['points'] for z in B['rows'] if z['operator']=='combined') if p[0]==.99)
    eos0,eos1=StandardModel(eos),SplineEOS(eos)
    A,C=[CombinedExpansion(E,PositivePortal((17,44,32)),168) for E in (eos0,eos1)]
    k=np.log(.3595)
    path=local_path(A,seed)
    x=brentq(lambda z:np.log(path(z)[1])-k,.99,.9911,xtol=5e-15)
    a,b=[slopes(E,k) for E in (eos0,eos1)]
    if x-.99<1e-4 or abs(a[0]-a[1])<.04 or abs(b[0]-b[1])>2e-13:
        raise ArithmeticError('The targeted near .99 knot did not separate the interpolation regularities')
    rows=[]
    for dt in (8e-5,2.5e-5):
        before,after=[path(x+d) for d in (-dt,dt)]
        if not before[1]>.3595>after[1]:
            raise ArithmeticError('Require two controlled sides of the exact cooling knot')
        data=[]
        for p in (before,after):
            vv={}
            for tag,E in zip(('pchip','cubic'),(A,C),strict=True):
                hi=measure(E,p,2e-5,1e-6)
                lo=measure(E,p,4e-5,2e-6)
                err=max(abs(np.asarray(hi['signed_squared_norm_rate'])-lo['signed_squared_norm_rate']))
                if err>1e-5:
                    raise ArithmeticError(('One-sided rate still unresolved',p[0],err))
                vv[tag]=dict(max_rate=hi['signed_squared_norm_rate'][-1],defect=err)
            data.append(dict(x=p[0],visible_T=p[1],eos_results=vv))
        rows.append(dict(delta_x=dt,before_and_after=data))
    return dict(original_sample_x=.99,crossing_x=x,visible_temperature_knot=.3595,
                dlogcv_dlogT_cooler_and_hotter=dict(pchip=a,cubic=b),rows=rows,input_sha256=H,
                table_input_sha256=h,scope='Original evolved trajectory continued across one fixed normalized EOS knot; cubic comparison at exactly matched held states. No differentiability or physical sign claim at the knot.')
    

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','late','interp'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
