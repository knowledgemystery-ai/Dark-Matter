"""Attribute the signed entropy-norm rate on the already selected local directions."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff, full
from fd_portal_local_entropy_metric_rates import motion
from sm_dark_entropy_exchange_reference import StandardModel


def visible(stage,p):
    x,T,eta,logt=p[:4]
    v,_,_,cv=stage.eos.state(T)
    z=stage.bath([eta,logt])
    q=np.array([0.,-z['rho'],v+z['rho']])
    return np.exp(3*x)*np.outer(q,q)/(T*T*cv)


def source(stage,p,row):
    F=full(stage,p[1:4])
    src=F(np.ones(3))
    J=diff(F,2e-4)
    L=[a-np.diag(b) for a,b in zip(J,src,strict=True)]
    directed,_=motion(stage,p,src.sum(axis=0))
    p=np.asarray(p[:4])
    G,_=matrix(stage,p)
    V=visible(stage,p)
    x,T=p[:2]
    k=np.log(T)
    E=stage.eos
    g,g1,g2=E.gs(k),E.slope(k),E.slope.derivative()(k)
    cv_log_slope=3.+(3*g1+g2)/(3*g+g1)
    shift=directed[1]/T
    visible_heat_capacity_rate_over_cv=-cv_log_slope*shift
    w=np.asarray(row['largest_direction_fractional_coordinates_G_normalized'])
    e=w@G@w
    h=1e-5
    def drift(f):
        return (f(p-2*h*directed)-8*f(p-h*directed)+8*f(p+h*directed)-f(p+2*h*directed))/(12*h)
    Md=drift(lambda z:matrix(stage,z)[0])
    Vd=drift(lambda z:visible(stage,z))
    heat=(w@V@w/e)*visible_heat_capacity_rate_over_cv
    # Exact radiation-limit reference d log(c_v)/d log(T)=3, while all
    # trajectory values and other terms remain at the archived EOS.
    heat_departure=(w@V@w/e)*(-(cv_log_slope-3.)*shift)
    a,b=[float(w@(m.T@G+G@m)@w/e) for m in L]
    c=float(w@(Md-Vd)@w/e)
    rest=float(w@Vd@w/e-heat)
    z=[a,b,c,rest,heat]
    total=sum(z)
    prior=row['signed_squared_norm_rate'][-1]
    if abs(total-prior)>3e-9 or abs(e-1.)>2e-12:
        raise ArithmeticError('Rate attribution did not close on the accepted direction')
    return dict(x=x,visible_temperature=T,dlog_visible_heat_capacity_dlog_temperature=cv_log_slope,
                relative_visible_temperature_rate=shift,squared_entropy_norm_rate=total,
                terms=dict(fractional_flow_expansion=a,fractional_flow_collision=b,dark_metric_motion=c,
                           other_visible_metric_motion=rest,visible_heat_capacity_motion=heat),
                visible_heat_capacity_departure_from_T_cubed=heat_departure,
                diagnostic_rate_only_replacing_visible_heat_capacity_log_slope_by_3=total-heat_departure,
                original_direction= w.tolist())


def run(eos,early,late,local):
    raw=[r.read_bytes() for r in (early,late,local)]
    ha=[hashlib.sha256(r).hexdigest() for r in raw]
    A,B,W=map(json.loads,raw)
    he=hashlib.sha256(eos.read_bytes()).hexdigest()
    if W['input_sha256'][:2]!=ha[:2] or not A['table_input_sha256']==B['table_input_sha256']==W['table_input_sha256']==he:
        raise ArithmeticError('Exact published trajectories and local directions required')
    u,v=[next(a['points'] for a in z['rows'] if a['operator']=='combined') for z in (A,B)]
    p=u+v[1:]
    if u[-1]!=v[0] or [s[0] for s in p]!=[r['x'] for r in W['rows']]:
        raise ArithmeticError('Archive and local direction grid must match exactly')
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    return dict(rows=[source(P,x,r) for x,r in zip(p,W['rows'],strict=True)],
                input_sha256=ha,table_input_sha256=he,
                scope='Algebraic metric/flow attribution on the exact trajectory. Changing only a derivative is a diagnostic, not an alternative EOS or self-consistent evolution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','early','late','local'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
