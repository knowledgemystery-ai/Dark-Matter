"""Check the local tangent with independent short nonlinear perturbations."""
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
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_entropy_short_tangent_windows import DATE, gain_difference
from sm_dark_entropy_exchange_reference import StandardModel


def endpoint_columns(evolve,Q,eps,control):
    z=[]
    for w in Q.T:
        a,b=[evolve(1+sgn*eps*w) for sgn in (-1,1)]
        z.append((b-a)/(2*eps*control))
    return np.stack(z,axis=2)


def run(eos,root):
    names=('direct-expansion-096-120','entropy-short-tangent-windows')
    raw=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    B,W=map(json.loads,raw)
    ids=[hashlib.sha256(r).hexdigest() for r in raw]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if W['input_sha256'][0]!=ids[0] or not B['table_input_sha256']==W['table_input_sha256']==h:
        raise ArithmeticError('Need exact original local tangent and archived path')
    from_receipt=lambda i: W['windows'][i]['x']
    times=np.array([*from_receipt(0),from_receipt(1)[1],from_receipt(2)[1]])
    seed=next(z for z in next(a['points'] for a in B['rows'] if a['operator']=='combined') if z[0]==times[0])
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    f=full(P,seed[1:4]); p=local_path(P,seed)
    G=[matrix(P,p(t))[0] for t in times]
    Q=np.asarray(W['rebased_fractional_columns'])
    def rhs(t,y):
        s=np.exp(3*(t-times[0]))
        return s*f(y/s).sum(axis=0)
    inds=((0,1),(1,2),(2,3),(0,3))
    result=[]
    for eps,step in ((2e-5,4e-5),(1e-5,2e-5)):
        def evolve(initial):
            s=solve_ivp(rhs,(times[0],times[-1]),initial,method='DOP853',rtol=7e-14,atol=5e-16,t_eval=times,max_step=step)
            if not s.success or not np.array_equal(s.t,times) or min(s.y.ravel())<=0.:
                raise ArithmeticError('Independent short nonlinear trajectory incomplete')
            return s.y.T
        z=evolve(np.ones(3))
        r=endpoint_columns(evolve,Q,eps,z)
        S=[a.T@b@a for a,b in zip(r,G,strict=True)]
        result.append([gain_difference(S[i],S[j]).tolist() for i,j in inds])
    v=np.array(result); ref=np.array([w['ordered_gain_minus_one'] for w in W['windows']])
    err=float(np.max(abs(v[1]-ref)))
    a=float(np.max(abs(v[1]-v[0])))
    if err>2e-7 or a>2e-7 or not np.array_equal(np.sign(v[1,:,-1]),np.sign(ref[:,-1])):
        raise ArithmeticError(('Independent short nonlinear sign control did not resolve',err,a))
    return dict(ordered_and_finite_gain_minus_one=[ref.tolist(),*result],maximum_tangent_comparison_defect=err,
                maximum_two_finite_setting_defect=a,input_sha256=ids,table_input_sha256=h,
                scope='Independent short nonlinear continuations initialized on the archived tangent plane at .99, not full independent .55-to-.99 trajectories.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
