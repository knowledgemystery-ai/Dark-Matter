"""Check instantaneous entropy-metric norm rates along the restricted FD trajectory."""
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
from sm_dark_entropy_exchange_reference import StandardModel


def motion(stage, point, f):
    x,T,e,l=point[:4]
    dark=stage.bath([e,l])
    v,_,_,cv=stage.eos.state(T)
    n,r,B,M,C=[dark[k] for k in ('n','rho','dn_deta','chemical_energy','cv_fixed_number')]
    td=np.exp(l)
    J=np.array([[B,B*M/td],[B*M,td*C+B*M*M/td]])
    q=np.array([n,r,r+v])
    # RHS rate is for abundance multiplied by a^3. The state inversion uses densities.
    dq=q*(f-3.)
    dz=np.r_[1.,(dq[2]-dq[1])/cv,np.linalg.solve(J,dq[:2])]
    # Differentiate the independent EOS/bath map directly to audit the motion conversion.
    def densities(p):
        v=stage.eos.state(p[1])[0]
        b=stage.bath(p[2:4])
        return np.array([b['n'],b['rho'],b['rho']+v])
    step=1e-5
    p=np.array(point[:4])
    g=(densities(p+step*dz)-densities(p-step*dz))/(2*step)
    err=max(abs((g-dq)/q))
    if err>3e-7:
        raise ArithmeticError(('Flow motion conversion unresolved',x,err))
    return dz,err


def measure(P, p, h, k):
    r=full(P,p[1:4])
    q_f=r(np.ones(3)).sum(axis=0)
    J=diff(r,h).sum(axis=0)-np.diag(q_f)
    v,flow_err=motion(P,p,q_f)
    p0=np.asarray(p[:4])
    G,err=matrix(P,p0)
    def mat(t):
        return matrix(P,p0+t*v)[0]
    def eval_at(dx):
        g=(mat(-2*dx)-8*mat(-dx)+8*mat(dx)-mat(2*dx))/(12*dx)
        S=g+J.T@G+G@J
        return (S+S.T)*.5
    S=eval_at(k)
    S0=eval_at(2*k)
    a,U=eigh(S,G)
    a0=eigh(S0,G,eigvals_only=True)
    b=eigh(S[:2,:2],G[:2,:2],eigvals_only=True)
    w=U[:,-1]
    # Independent symmetric secant using the linearized flow and moving endpoint
    # metrics. This is a consistency check, not a finite nonlinear propagation.
    h_t=7e-6
    xp=w+h_t*(J@w)
    xm=w-h_t*(J@w)
    check=(xp@mat(h_t)@xp-xm@mat(-h_t)@xm)/(2*h_t)/(w@G@w)
    if abs(check-a[-1])>2e-6 or max(abs(a0-a))>2e-6:
        raise ArithmeticError(('Metric rate does not resolve',p[0],check,a[-1],a0-a))
    return dict(x=p[0],signed_squared_norm_rate=a.tolist(),instantaneous_fixed_total_plane_rate=b.tolist(),
                largest_direction_fractional_coordinates_G_normalized=w.tolist(),
                entropy_metric_fd_relative_defect=err,density_motion_defect=flow_err,
                metric_derivative_mesh_defect=max(abs(a-a0)),independent_direction_rate_defect=abs(check-a[-1]))


def run(eos,early,late,receipt):
    raw=[f.read_bytes() for f in (early,late,receipt)]
    A,B,M=map(json.loads,raw)
    sh=[hashlib.sha256(s).hexdigest() for s in raw]
    table_hash=hashlib.sha256(eos.read_bytes()).hexdigest()
    if ([A['table_input_sha256'],B['table_input_sha256'],M['table_input_sha256']]!=[table_hash]*3
            or M['input_sha256'][:2]!=sh[:2] or B['previous_segment_sha256']!=sh[0]):
        raise ArithmeticError('Exact shared trajectories and entropy norm receipt required')
    u,v=[next(z['points'] for z in k['rows'] if z['operator']=='combined') for k in (A,B)]
    if u[-1]!=v[0] or (u[0][0],v[0][0],v[-1][0])!=(.55,.96,1.2):
        raise ArithmeticError('Physical background segments must share the exact same bridge')
    eos_obj=StandardModel(eos)
    stages=[CombinedExpansion(eos_obj,PositivePortal(z),n) for z,n in (((13,36,26),112),((17,44,32),168))]
    result=[]
    for point in u+v[1:]:
        low,high=[measure(z,point,h,k) for z,h,k in zip(stages,(4e-4,2e-4),(2e-5,1e-5),strict=True)]
        # Store separate absolute derivative/mesh deltas rather than silently accepting a sign.
        delta={key:float(max(abs(np.asarray(high[key])-low[key]))) for key in
               ('signed_squared_norm_rate','instantaneous_fixed_total_plane_rate')}
        if max(delta.values())>8e-5:
            raise ArithmeticError(('Quadrature / perturbation refinement rate did not resolve',point[0],delta))
        result.append(dict(**high,combined_mesh_rate_defects=delta))
    return dict(x_first=.55,x_last=1.2,burst='combined',norm='negative comoving entropy Hessian in moving fractional coordinates',
                input_sha256=sh,table_input_sha256=table_hash,rows=result,
                scope='Instantaneous squared-norm growth exponents at archived samples; locally fixed total plane is not the transported initial plane. Not an interval monotonicity theorem or physical stability result.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','early','late','receipt'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
