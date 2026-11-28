"""Resolve the held-EOS threshold change outside the rank-one amplitude."""
import argparse
import hashlib
import json
from itertools import combinations
from math import comb, fsum
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_cubic_derivative_cross import slope_and_rank
from fd_portal_cubic_rank_factor_attribution import density
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff, full
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_entropy_knot_jump_audit import analytic_entropy_metric, forms
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_heat_capacity_slope_threshold import DATE
from fd_portal_local_entropy_metric_rates import motion
from fd_portal_visible_curvature_rate_attribution import visible
from sm_dark_entropy_exchange_reference import StandardModel

NAMES=('dark_metric_motion','visible_metric_motion_excluding_k','expansion_flow','collision_flow')


def swap_game(base,new,rank):
    N=len(base)
    D=np.asarray(new)-np.asarray(base)
    B=np.sum(base,axis=0)
    def value(selected):
        X=B+sum((D[i] for i in selected),start=np.zeros_like(B))
        if max(eigh(X,eigvals_only=True))>=0.:
            raise ArithmeticError('All held hybrids must be dissipative at reference slope 3')
        return 3.+1./max(eigh(rank,-X,eigvals_only=True))
    whole=tuple(range(N))
    allv={S:value(S) for k in range(N+1) for S in combinations(whole,k)}
    ans=[]
    for j in whole:
        ans.append(fsum((allv[tuple(sorted((*S,j)))]-v)/(N*comb(N-1,len(S))) for S,v in allv.items() if j not in S))
    total=allv[whole]-allv[()]
    if abs(fsum(ans)-total)>2e-15:
        raise ArithmeticError('Symmetric replacement allocations must close')
    return [allv[()],allv[whole]],ans,[allv[(j,)]-allv[()] for j in whole]


def components(P,p):
    k,R=slope_and_rank(P,p)
    f=full(P,p[1:4])
    a=f(np.ones(3))
    L=diff(f,2e-5)
    if len(a)!=2 or len(L)!=2:
        raise ArithmeticError('Need expansion and collision blocks')
    v,_=motion(P,p,a.sum(axis=0))
    p0=np.array(p[:4])
    G,S=forms(P,p)
    def derivative(function):
        h=1e-6
        C=[function(P,p0+j*h*v) for j in (-2,-1,1,2)]
        return (C[0]-8*C[1]+8*C[2]-C[3])/(12*h)
    V=derivative(visible)
    F=derivative(analytic_entropy_metric)
    D=F-V
    T=[D,V+(3.-k)*R]
    for u,z in zip(L,a,strict=True):
        M=u-np.diag(z)
        T.append(M.T@G+G@M)
    T=[(A+A.T)*.5 for A in T]
    if np.linalg.norm(np.sum(T,axis=0)-(S+(3.-k)*R),2)>2e-12*np.linalg.norm(S,2):
        raise ArithmeticError('Individual terms must reconstruct the exact held B form')
    return T,R


def run(eos,root):
    names=('direct-expansion-096-120','cubic-derivative-cross','cubic-rank-factor-attribution')
    raws=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    L,C,A=map(json.loads,raws)
    ids=[hashlib.sha256(s).hexdigest() for s in raws]
    he=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (C['input_sha256'][1]!=ids[0] or A['input_sha256'][0]!=ids[0] or A['input_sha256'][-1]!=ids[1]
        or any(x['table_input_sha256']!=he for x in (L,C,A))):
        raise ArithmeticError('Exact prior two-factor result required')
    P,Q=[CombinedExpansion(f(eos),PositivePortal((17,44,32)),168) for f in (StandardModel,SplineEOS)]
    hist=next(r['points'] for r in L['rows'] if r['operator']=='combined')
    walk=local_path(P,next(p for p in hist if p[0]==.99))
    out=[]
    for z in A['rows']:
        p=walk(z['x'])
        if abs(p[1]-z['T'])>3e-14:
            raise ArithmeticError('Exact held knot states required')
        B,R=components(P,p)
        C,_=components(Q,p)
        D=np.diag(density(P,p)/density(Q,p))
        C=[D.T@T@D for T in C]
        thresholds,shapley,single=swap_game(B,C,R)
        if max(abs(np.asarray(thresholds)-np.asarray(z['critical_slope_by_base_rows_and_rank_columns'])[:,0]))>5e-9:
            raise ArithmeticError('Independent split must preserve both previous thresholds')
        out.append(dict(x=p[0],T=p[1],thresholds_at_pchip_rank=thresholds,symmetric_all_orders_allocations=dict(zip(NAMES,shapley,strict=True)),
                        single_swaps_at_pchip_base=dict(zip(NAMES,single,strict=True))))
    return dict(rows=out,replace_order_scope='The allocations average increments over every order of the four additive form terms.',
                input_sha256=ids,table_input_sha256=he,
                scope='Algebraic threshold comparison at common density perturbations and fixed PCHIP rank amplitude; not a causal EOS or stability result.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
