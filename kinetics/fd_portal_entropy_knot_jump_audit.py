"""Compare the measured moving metric jump with its predicted rank-one term."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff, full
from fd_portal_entropy_pchip_knot_crossing import local_path, slopes
from fd_portal_local_entropy_metric_rates import motion
from fd_portal_visible_curvature_rate_attribution import visible
from sm_dark_entropy_exchange_reference import StandardModel


def analytic_entropy_metric(A,p):
    x,T,eta,logt=p[:4]
    td=np.exp(logt)
    z=A.bath([eta,logt])
    n,r,B,M,C=[z[k] for k in ('n','rho','dn_deta','chemical_energy','cv_fixed_number')]
    J=np.array([[B,B*M/td],[B*M,td*C+B*M*M/td]])
    K=np.diag([1.,1./td])@np.linalg.inv(J)
    out=visible(A,p)
    q=np.array([n,r])
    out[:2,:2]+=np.exp(3*x)*(q[:,None]*q[None,:])*K
    return (out+out.T)*.5


def forms(P,p):
    G,e=matrix(P,p)
    f=full(P,p[1:4])
    a=f(np.ones(3)).sum(axis=0)
    L=diff(f,2e-5).sum(axis=0)-np.diag(a)
    v,_=motion(P,p,a)
    x=np.asarray(p[:4])
    d=1e-6
    c,b,u,w=[analytic_entropy_metric(P,x+k*d*v) for k in (-2,-1,1,2)]
    S=(c-8*b+8*u-w)/(12*d)+L.T@G+G@L
    if np.max(abs(G-analytic_entropy_metric(P,x))) > 3e-12*np.max(abs(G)):
        raise ArithmeticError('Analytic metric value does not match the independent endpoint metric')
    return G,(S+S.T)*.5


def run(eos,late,knot):
    a,b=[r.read_bytes() for r in (late,knot)]
    A,B=map(json.loads,(a,b))
    h=[hashlib.sha256(u).hexdigest() for u in (a,b)]
    t=hashlib.sha256(eos.read_bytes()).hexdigest()
    if B['input_sha256'][0]!=h[0] or A['table_input_sha256']!=t or B['table_input_sha256']!=t:
        raise ArithmeticError('Exact short knot run and original path required')
    p=next(u for u in next(r['points'] for r in A['rows'] if r['operator']=='combined') if u[0]==.99)
    E=StandardModel(eos)
    P=CombinedExpansion(E,PositivePortal((17,44,32)),168)
    x=B['crossing_x']
    walk=local_path(P,p)
    state=walk(x)
    tc=B['visible_temperature_knot']
    if abs(state[1]-tc)>3e-14:
        raise ArithmeticError('Exact thermal knot reconstruction failed')
    # Only the derivative of the visible heat-capacity coefficient jumps;
    # its value, the motion, dark metric drift, and fixed-point flows are continuous.
    slopes_cold_hot=slopes(E,np.log(tc))
    f=full(P,state[1:4])(np.ones(3)).sum(axis=0)
    vr,_,_,cv=E.state(state[1])
    rd=P.bath(state[2:4])['rho']
    Tdot=((vr+rd)*(f[2]-3.)-rd*(f[1]-3.))/cv
    # At the interpolation knot the metric slope can fail a high-order smooth
    # difference, but the temperature rate itself must match the actual path.
    h0=1e-7
    check=(walk(x+h0)[1]-walk(x-h0)[1])/(2*h0)
    if abs(Tdot-check)>1e-7:
        raise ArithmeticError('One-sided EOS heat slope cannot change the continuous physical temperature rate')
    k=-(Tdot/state[1])*(slopes_cold_hot[0]-slopes_cold_hot[1])
    vis=visible(P,state)
    jump=k*vis
    G0=analytic_entropy_metric(P,state)
    data=[]
    for z in B['rows']:
        h0=z['delta_x']
        pb,pa=[walk(x+i*h0) for i in (-1,1)]
        (Gb,Sb),(Ga,Sa)=[forms(P,w) for w in (pb,pa)]
        norm_rel=np.linalg.norm(Sa-Sb-jump,2)/np.linalg.norm(jump,2)
        estimate=eigh(Sb+jump,G0,eigvals_only=True)[-1]
        exact=eigh(Sa,Ga,eigvals_only=True)[-1]
        if norm_rel>1.5e-2 or abs(estimate-exact)>.003:
            raise ArithmeticError('Observed jump not captured by the rank-one metric term')
        data.append(dict(delta_x=h0,operator_relative_rank_one_residual=norm_rel,
                         reported_cooler_side_max=exact,rank_one_adjusted_hot_form_max=estimate))
    return dict(rows=data,expected_scalar_multiplier_of_visible_metric=k,input_sha256=h,table_input_sha256=t,
                scope='One-sided limits of a numerical EOS metric across one knot; compares finite probes, not a continuum precision certificate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for z in ('eos','late','knot'):
        p.add_argument(z,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
