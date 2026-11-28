"""Track the truly transported initial fixed-total-energy plane."""
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
from fd_portal_direct_comoving_variational_jacobian import diff, full
from fd_portal_direct_comoving_variational_propagator import solve
from fd_portal_local_entropy_metric_rates import motion
from sm_dark_entropy_exchange_reference import StandardModel

DATE='2026-11-28'


def flow(times, table, query, atol):
    f=CubicSpline(times,table)
    E=np.eye(3)[:,:2]
    def rhs(t,v):
        return (f(t)@v.reshape(3,2)).ravel()
    r=solve_ivp(rhs,(times[0],times[-1]),E.ravel(),t_eval=query,atol=atol,rtol=2e-11,
                method='DOP853',max_step=.006)
    if not r.success or not np.array_equal(r.t,query):
        raise ArithmeticError('Transported plane propagation incomplete')
    return r.y.T.reshape(-1,3,2)


def forms(stage,p):
    g,k=matrix(stage,p)
    f=full(stage,p[1:4])
    q=f(np.ones(3)).sum(axis=0)
    j=diff(f,2e-4).sum(axis=0)-np.diag(q)
    v,_=motion(stage,p,q)
    p0=np.asarray(p[:4])
    h=1e-5
    z=[matrix(stage,p0+k*v)[0] for k in (-2*h,-h,h,2*h)]
    s=(z[0]-8*z[1]+8*z[2]-z[3])/(12*h)+j.T@g+g@j
    return g,(s+s.T)*.5


def run(eos,root,local):
    names=('direct-expansion-055-096','direct-expansion-096-120','comoving-response-055-096',
           'direct-comoving-variational-propagator','comoving-matched-response-055-120')
    raw=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]+[local.read_bytes()]
    A,B,E,L,C,W=map(json.loads,raw)
    sha=[hashlib.sha256(x).hexdigest() for x in raw]
    se=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (W['input_sha256'][:2]!=sha[:2] or C['input_sha256']!=sha[2:4]
        or not A['table_input_sha256']==B['table_input_sha256']==W['table_input_sha256']==se):
        raise ArithmeticError('Archived FD histories, local forms and ordered witnesses must match')
    p,q=[next(r['points'] for r in x['rows'] if r['operator']=='combined') for x in (A,B)]
    if p[-1]!=q[0] or (p[0][0],q[-1][0])!=(.55,1.2):
        raise ArithmeticError('Require exact archived bridge state')
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    rates,_,base,baseline_error=solve(P,p+q[1:],8e-12,.012)
    # Compute the tangent in fixed initial units; convert to endpoint fractions
    # before projecting the local metric. Differentiation scale is local baseline.
    def fixed(x):
        y=base(x)
        return diff(lambda v:rates(x,v*y).sum(axis=0),2e-4)/y[None,:]
    query=np.array([.61,.96,.99,1.2])
    # A cusp at a normalized EOS knot can affect inverse thermodynamic
    # derivatives. Do not assume a formal spline order across all knots.
    t=np.linspace(.55,1.2,129)
    tab=np.array([fixed(s) for s in t])
    spans=[flow(t[::k],tab[::k],query,2e-13)/base(query).T[:,:,None] for k in (2,1)]
    err=max(np.max(abs(spans[1][i]-np.asarray(C['rows'][j]['fixed_total_initial_to_all_fractions'])))
            for i,j in ((1,0),(3,1)))
    if err>2e-5:
        raise ArithmeticError(('Transported columns fail endpoint match',err))
    result=[]
    for x,i in ((.61,0),(.99,2)):
        point=next(z for z in p+q[1:] if z[0]==x)
        G,S=forms(P,point)
        exact=next(row for row in W['rows'] if row['x']==x)
        replay=eigh(S,G,eigvals_only=True)
        if max(abs(replay-exact['signed_squared_norm_rate']))>2e-10:
            raise ArithmeticError('Local full space form must replay')
        Us=[np.linalg.qr(z[i],mode='reduced')[0] for z in spans]
        vs=[eigh(U.T@S@U,U.T@G@U,eigvals_only=True) for U in Us]
        proj=np.linalg.norm(Us[0]@Us[0].T-Us[1]@Us[1].T,2)
        defect=float(max(abs(vs[0]-vs[1])))
        result.append(dict(x=x,full_space_maximum_squared_norm_rate=replay[-1],
                           transport_fixed_initial_total_plane_squared_norm_rate=vs[1].tolist(),
                           same_time_fixed_total_plane_squared_norm_rate=exact['instantaneous_fixed_total_plane_rate'],
                           max_rate_discrepancy_between_meshes=defect,max_projector_discrepancy_two_norm=proj,
                           endpoint_fractional_plane_columns=spans[1][i].tolist()))
    return dict(rows=result,root_input_names=[*names,'local'],input_sha256=sha,table_input_sha256=se,
                background_defect=baseline_error,known_endpoint_column_defect=err,sampled_jacobian_mesh_sizes=[65,129],
                scope='Projected local rates for the historically fixed-initial-total-energy linear tangent plane at two archived x; no continuum bracket or nonlinear contraction claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root','local'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
