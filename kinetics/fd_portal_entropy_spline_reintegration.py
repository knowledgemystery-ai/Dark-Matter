"""Re-evolve the restricted background under the alternate positive spline EOS."""
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


def integrate(stage, seed, times, tol, step):
    f=full(stage,seed[1:4])
    x0=seed[0]
    def rhs(x,y):
        a=np.exp(3*(x-x0))
        return a*f(y/a).sum(axis=0)
    sol=solve_ivp(rhs,(times[0],times[-1]),np.ones(3),method='DOP853',atol=.03*tol,rtol=tol,
                  t_eval=times,max_step=step)
    if not sol.success or not np.array_equal(sol.t,times):
        raise ArithmeticError('Alternate EOS background solve incomplete')
    d=stage.bath(seed[2:4])
    scale=np.array([d['n'],d['rho'],d['rho']+stage.eos.state(seed[1])[0]])
    points=[]
    for x,Q in zip(sol.t,sol.y.T,strict=True):
        n,r,u=scale*Q/np.exp(3*(x-x0))
        if min(n,r,u-r)<=0.:
            raise ArithmeticError('Physical alternate snapshots required')
        T=brentq(lambda T:stage.eos.state(T)[0]/(u-r)-1.,stage.lo,stage.hi,xtol=5e-15)
        def invert(z):
            v=stage.bath(z)
            return np.array([v['n']/n-1.,v['rho']/r-1.])
        z=root(invert,points[-1][2:4] if points else seed[2:4],tol=2e-12).x
        if max(abs(invert(z)))>3e-12:
            raise ArithmeticError('Alternate dark-state inversion unresolved')
        points.append([x,T,*z])
    return np.asarray(points),scale


def run(eos,early,late,held):
    data=[r.read_bytes() for r in (early,late,held)]
    ids=[hashlib.sha256(v).hexdigest() for v in data]
    A,B,D=map(json.loads,data)
    eos_id=hashlib.sha256(eos.read_bytes()).hexdigest()
    if D['input_sha256'][:2]!=ids[:2] or not A['table_input_sha256']==B['table_input_sha256']==D['table_input_sha256']==eos_id:
        raise ArithmeticError('Reintegration requires the exact held sensitivity study and original samples')
    u,v=[next(a['points'] for a in z['rows'] if a['operator']=='combined') for z in (A,B)]
    if u[-1]!=v[0] or (u[0][0],v[-1][0])!=(.55,1.2):
        raise ArithmeticError('Wrong original reference path')
    old=np.asarray([z[:4] for z in u+v[1:]])
    t=old[:,0]
    P=CombinedExpansion(SplineEOS(eos),PositivePortal((17,44,32)),168)
    baseline=StandardModel(eos)
    g,r=[integrate(P,old[0],t,*p) for p in ((4e-11,.016),(7e-12,.008))]
    de=max(abs((g[0]-r[0]).ravel()))
    if de>8e-8 or max_abs_state_difference(r[0][0],old[0])>2e-13:
        raise ArithmeticError(('Reintegrated spline snapshots failed to resolve',de))
    rows=[]
    for original,p,held_row in zip(old,r[0],D['rows'],strict=True):
        v0=measure(P,p,2e-4,1e-5)
        a=v0['signed_squared_norm_rate'][-1]
        rows.append(dict(x=p[0],original_state=original.tolist(),alternate_reintegrated_state=p.tolist(),
                         original_max=held_row['original_max_rate'],alternate_held_state_max=held_row['alternate_max_rate'],
                         alternate_reintegrated_max=a,alternate_reintegrated_fixed_total_max=v0['instantaneous_fixed_total_plane_rate'][-1]))
    return dict(rows=rows,background_ode_refinement_absolute_defect=de,input_sha256=ids,table_input_sha256=eos_id,
                alternate_minus_original_initial_total_energy_relative=r[1][2]/(baseline.state(old[0,1])[0]+P.bath(old[0,2:4])['rho'])-1.,
                initial_condition='Same visible T, dark eta/logT; each EOS determines its own visible density and Friedmann rate.',
                scope='Separate restricted background with one alternate interpolation, no physical EOS error interval or claim between archived x samples.')


def max_abs_state_difference(a,b):
    return max(abs(a-b))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','early','late','held'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
