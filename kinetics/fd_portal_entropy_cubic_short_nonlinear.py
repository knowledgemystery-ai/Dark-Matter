"""Independently check the smooth cubic finite signs with short nonlinear pairs."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import full
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_entropy_hot_start_finite_crossing import dense
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_entropy_short_nonlinear_witness import endpoint_columns
from fd_portal_entropy_short_tangent_windows import DATE,gain_difference
from fd_portal_entropy_spline_reintegration import integrate


def run(eos,root):
    names=('entropy-spline-reintegration','spline-entropy-endpoint-gain','entropy-cubic-short-tangent')
    raw=[(root/f'fd-portal-{n}-{DATE}.json').read_bytes() for n in names]
    A,B,T=map(json.loads,raw)
    ids=[hashlib.sha256(x).hexdigest() for x in raw]
    sha=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (T['input_sha256'][:2]!=ids[:2] or B['input_sha256'][0]!=ids[0]
        or any(s['table_input_sha256']!=sha for s in (A,B,T))):
        raise ArithmeticError('Require exact completed cubic tangent and alternate state')
    r,s=[next(z['alternate_reintegrated_state'] for z in A['rows'] if z['x']==x) for x in (.96,.99)]
    P=CombinedExpansion(SplineEOS(eos),PositivePortal((17,44,32)),168)
    walk=local_path(P,s)
    r0=T['windows']
    clocks=np.array([*r0[0]['x'],r0[1]['x'][1],r0[2]['x'][1]])
    g=[matrix(P,walk(x))[0] for x in clocks]
    C=np.asarray(T['initial_rebased_fractional_columns'])
    ref=np.array([z['ordered_gain_minus_one'] for z in r0])
    idx=((0,1),(1,2),(2,3),(0,3))
    f=full(P,s[1:4])
    def rhs(t,y):
        a=np.exp(3*(t-clocks[0]))
        return a*f(y/a).sum(axis=0)
    cases=[]
    rebase=[]
    for mesh,jstep,eps,step in ((4,4e-5,2e-5,4e-5),(6,2e-5,1e-5,2e-5)):
        way=[.96,.97,.98,.99]
        x=np.unique(np.concatenate([np.linspace(a,b,mesh) for a,b in zip(way[:-1],way[1:])]))
        Y,_=integrate(P,r,x,2e-12,.0025)
        p=dict(zip(x,Y,strict=True))
        if max(abs(p[.99]-s))>1e-9:
            raise ArithmeticError('Re-evolved smooth state at .99 did not replay')
        C99=dense(P,p.__getitem__,way,C,mesh,jstep)(.99)
        Q,R=np.linalg.qr(C99,mode='reduced')
        if min(abs(R.diagonal()))<1e-5:
            raise ArithmeticError('Local smooth plane is degenerate')
        def evolve(y):
            z=solve_ivp(rhs,(clocks[0],clocks[-1]),y,method='DOP853',t_eval=clocks,rtol=7e-14,atol=5e-16,max_step=step)
            if not z.success or not np.array_equal(z.t,clocks) or np.min(z.y)<=0.:
                raise ArithmeticError('Smooth nonlinear short pair incomplete')
            return z.y.T
        V=endpoint_columns(evolve,Q,eps,evolve(np.ones(3)))
        W=[v.T@m@v for v,m in zip(V,g,strict=True)]
        cases.append([gain_difference(W[i],W[j]).tolist() for i,j in idx])
        rebase.append(Q.tolist())
    x,y=map(np.asarray,cases)
    defect=float(np.max(abs(y-ref))); repeat=float(np.max(abs(y-x)))
    if max(defect,repeat)>2e-8 or max(y[:,-1])>=0.:
        raise ArithmeticError(('Short smooth nonlinear control did not confirm retained signs',defect,repeat))
    return dict(windows=[z['x'] for z in r0],tangent_and_two_direct_gain_minus_one=[ref.tolist(),*cases],
                max_nonlinear_vs_tangent_defect=defect,max_two_setting_defect=repeat,
                local_rebased_planes=rebase,input_sha256=ids,table_input_sha256=sha,
                scope='Independent local nonlinear pairs at the alternate .99 state on its own transported plane, not a new full history or physical stability.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
