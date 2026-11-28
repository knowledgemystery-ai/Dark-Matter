"""Independent signed nonlinear pulse in the selected transported entropy direction."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import full
from fd_portal_local_entropy_metric_rates import motion
from fd_portal_transported_entropy_plane import forms
from sm_dark_entropy_exchange_reference import StandardModel

D='2026-11-28'


def centered_rhs(f,t,m,z,eps):
    up=f(t,m+eps*z)
    dn=f(t,m-eps*z)
    return (up+dn)*.5,(up-dn)/(2*eps)


def pair(f,start,end,direction,eps,step):
    # Exactly equivalent to two independently signed nonlinear states in exact
    # arithmetic. Conditioning the difference as a state avoids subtracting
    # two separately accumulated O(1) integration errors at the endpoint.
    def rhs(t,w):
        mid,z=w.reshape(2,3)
        a,b=centered_rhs(f,t,mid,z,eps)
        return np.r_[a,b]
    a=solve_ivp(rhs,(start,end),np.r_[np.ones(3),direction],rtol=4e-11,atol=8e-13,
                method='DOP853',max_step=step)
    if not a.success:
        raise ArithmeticError('Signed nonlinear midpoint/difference integration did not finish')
    m,z=a.y[:,-1].reshape(2,3)
    return m,z,centered_rhs(f,end,m,z,eps)[1]


def run(eos,r):
    names=('direct-expansion-055-096','direct-expansion-096-120','transported-entropy-plane')
    raws=[(r/f'fd-portal-{k}-{D}.json').read_bytes() for k in names]
    h=[hashlib.sha256(x).hexdigest() for x in raws]
    A,B,C=map(json.loads,raws)
    e=hashlib.sha256(eos.read_bytes()).hexdigest()
    if C['input_sha256'][:2]!=h[:2] or not A['table_input_sha256']==B['table_input_sha256']==C['table_input_sha256']==e:
        raise ArithmeticError('Exact trajectory and exact previously selected plane required')
    u,v=[next(z['points'] for z in p['rows'] if z['operator']=='combined') for p in (A,B)]
    if u[-1]!=v[0] or u[0][0]!=.55:
        raise ArithmeticError('Wrong background')
    p0=u[0]
    x,vp=.99,next(z for z in v if z[0]==.99)
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    J=np.array(next(z['endpoint_fractional_plane_columns'] for z in C['rows'] if z['x']==x))
    G,S=forms(P,vp)
    U,R=np.linalg.qr(J,mode='reduced')
    rates,vec=eigh(U.T@S@U,U.T@G@U)
    d=np.linalg.solve(R,vec[:,-1])
    d/=np.linalg.norm(d)
    if d[0]<0.:
        d*=-1.
    expected=J@d
    f=full(P,p0[1:4])
    def fixed(t,y):
        z=np.exp(3*(t-p0[0]))
        return z*f(y/z).sum(axis=0)
    ctrl=solve_ivp(fixed,(p0[0],x),np.ones(3),rtol=8e-12,atol=3e-13,method='DOP853',max_step=.012)
    if not ctrl.success:
        raise ArithmeticError('Independent unperturbed trajectory did not finish')
    y=ctrl.y[:,-1]
    dy= fixed(x,y)
    z=np.array(vp[:4])
    a=full(P,vp[1:4])(np.ones(3)).sum(axis=0)
    dp, _=motion(P,vp,a)
    if max(abs(dy/y-a))>1e-6:
        raise ArithmeticError('Unperturbed archived local motion did not replay')
    h0=1e-5
    m=[matrix(P,z+j*h0*dp)[0] for j in (-2,-1,1,2)]
    Gdot=(m[0]-8*m[1]+8*m[2]-m[3])/(12*h0)
    rows=[]
    for eps in (2e-4,1e-4,5e-5):
        b,z,dz=pair(fixed,p0[0],x,np.r_[d,0.],eps,.012)
        w=z/y
        wdot=dz/y-w*dy/y
        n=w@G@w
        if n<=0 or not np.all(np.isfinite(wdot)):
            raise ArithmeticError('Finite normalized perturbation invalid')
        g=(w@Gdot@w+2*wdot@G@w)/n
        rows.append(dict(epsilon=eps,squared_norm_rate=g,dark_and_total_fractional_response=w.tolist(),
                         finite_vs_tangent_absolute=max(abs(w-expected)),
                         pair_midpoint_vs_control_fractional_defect=max(abs(b/y-1.))))
    return dict(projected_tangent_maximum=rates[-1],x_start=p0[0],x_probe=x,
                fixed_total_initial_unit_fractional_direction=np.r_[d,0.].tolist(),
                expected_tangent_fractional_output=expected.tolist(),rows=rows,input_sha256=h,table_input_sha256=e,
                scope='Finite signed probes with one original EOS and small perturbations. No physical stability or EOS robustness claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos',type=Path)
    p.add_argument('r',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
