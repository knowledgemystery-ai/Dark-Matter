"""Measure finite-window entropy gain on the carried fixed-initial-total plane."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq, root

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_propagator import solve
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from sm_dark_entropy_exchange_reference import StandardModel


def gain_on_plane(P0, P1, G0, G1):
    """The same initial perturbation coordinates must parameterize both times."""
    X,Y=(P.T@G@P for P,G in ((P0,G0),(P1,G1)))
    if min(np.linalg.eigvalsh(X))<=0 or min(np.linalg.eigvalsh(Y))<=0:
        raise ArithmeticError('Nondegenerate positive restricted forms required')
    return np.sqrt(eigh(Y,X,eigvals_only=True))


def window_indices():
    return ((0,1),(1,2),(2,3),(0,3))


def state(P, start, x, q):
    n,r,u=q*np.exp(-3*(x-start[0]))*np.array(start[4:7])
    if min(n,r,u-r)<=0:
        raise ArithmeticError('Positive baseline coordinates required')
    T=brentq(lambda v:P.eos.state(v)[0]/(u-r)-1,P.lo,P.hi,xtol=5e-15)
    def residual(z):
        b=P.bath(z)
        return [b['n']/n-1,b['rho']/r-1]
    z=root(residual,start[2:4],tol=2e-12).x
    if max(abs(np.array(residual(z))))>4e-13:
        raise ArithmeticError('Inversion failed')
    return [x,T,*z]


def finite_columns(nonlinear,control,times,eps):
    cols=[]
    for d in np.eye(3)[:2]:
        p=nonlinear(np.ones(3)+eps*d)(times)
        m=nonlinear(np.ones(3)-eps*d)(times)
        cols.append(((p-m)/(2*eps)/control(times)).T)
    return np.stack(cols,axis=2)


def evaluate(P,path,query,tol,step):
    _,nonlinear,base,background=solve(P,path,tol,step)
    dark=P.bath(path[0][2:4])
    begin=[*path[0][:4],dark['n'],dark['rho'],dark['rho']+P.eos.state(path[0][1])[0]]
    points=[state(P,begin,x,q) for x,q in zip(query,base(query).T,strict=True)]
    metrics,checks=zip(*(matrix(P,p) for p in points))
    wide,narrow=[finite_columns(nonlinear,base,query,e) for e in (.001,.0005)]
    richardson=(4*narrow-wide)/3
    ans=[]
    for i,j in window_indices():
        values=[gain_on_plane(M[i],M[j],metrics[i],metrics[j]) for M in (wide,narrow,richardson)]
        ans.append(dict(interval=[query[i],query[j]],largest_and_smallest=list(values[-1][::-1]),
                        amplitude_largest=[z[-1] for z in values]))
    return dict(selected_states=points,windows=ans,archived_control_max_defect=background,
                max_metric_check=float(max(checks)))


def run(eos,early,late,spline,knot):
    paths=(early,late,spline,knot)
    raw=[p.read_bytes() for p in paths]
    S=[hashlib.sha256(b).hexdigest() for b in raw]
    A,B,C,K=map(json.loads,raw)
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if any(z['table_input_sha256']!=h for z in (A,B,C,K)) or C['input_sha256'][:2]!=S[:2] or K['input_sha256'][0]!=S[1]:
        raise ArithmeticError('Require coherent exact restricted history and knot evidence')
    p,q=[next(z['points'] for z in a['rows'] if z['operator']=='combined') for a in (A,B)]
    u=[z['alternate_reintegrated_state'] for z in C['rows']]
    if p[-1]!=q[0] or [s[0] for s in p+q[1:]]!=[r[0] for r in u]:
        raise ArithmeticError('Expected exact late bridge and paired clocks')
    x=K['crossing_x']
    t=np.array([.99,x-2.5e-5,x+2.5e-5,x+8e-5])
    rows={}
    for label,fit,history in (('pchip',StandardModel,p+q[1:]),('cubic',SplineEOS,u)):
        P=CombinedExpansion(fit(eos),PositivePortal((17,44,32)),168)
        lo,hi=[evaluate(P,history,t,*args) for args in ((8e-12,.012),(2e-12,.006))]
        for a,b in zip(lo['windows'],hi['windows'],strict=True):
            v=np.array(b['largest_and_smallest'])
            b['ode_gain_difference']=float(np.max(abs(v-a['largest_and_smallest'])))
            b['amplitude_gain_difference']=float(np.ptp(b['amplitude_largest']))
        rows[label]=hi
    return dict(rows=rows,crossing_x=x,input_sha256=S,table_input_sha256=h,
                initial_signed_fractional_perturbations=[.001,.0005],time_refinement=[[8e-12,.012],[2e-12,.006]],
                scope='Four finite time windows in two restricted numerical FD histories. Both endpoints use the same two initial-energy-zero perturbation coordinates; no continuous stability or physical claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','early','late','spline','knot'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
