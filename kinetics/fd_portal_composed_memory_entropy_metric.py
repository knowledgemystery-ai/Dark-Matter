"""Check the physical entropy-curvature norm of the same fixed-total FD memory."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq,root

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from sm_dark_entropy_exchange_reference import StandardModel


def matrix(A,p):
    x,T,e,l=p[:4]
    td=np.exp(l)
    v,_,_,cv=A.eos.state(T)
    z=A.bath([e,l])
    n,r,B,M,C=[z[k] for k in ('n','rho','dn_deta','chemical_energy','cv_fixed_number')]
    J=np.array([[B,B*M/td],[B*M,td*C+B*M*M/td]])
    q=np.array([n,r,v+r])
    t=np.array([0.,-1.,1.])
    H=np.outer(t,t)/(T*T*cv)
    H[:2,:2]+=np.diag([1.,1./td])@np.linalg.inv(J)
    if np.max(abs(H-H.T))>2e-12*np.max(abs(H)):
        raise ArithmeticError('Physical entropy Hessian must be symmetric')
    cov=(q[:,None]*H*q[None,:]) * np.exp(3*x)

    def negative_gradient(f):
        N,R,V=q*f
        if V-R<=0:
            raise ArithmeticError('Visible density must stay positive')
        tv=brentq(lambda u:A.eos.state(u)[0]/(V-R)-1.,A.lo,A.hi,xtol=5e-15)
        def F(u):
            a=A.bath(u)
            return np.array([a['n']/N-1.,a['rho']/R-1.])
        a=root(F,[e,l],tol=2e-12).x
        if max(abs(F(a)))>4e-13:
            raise ArithmeticError('Independent perturbed dark state did not invert')
        u,vv=1/np.exp(a[1]),1/tv
        return np.array([a[0],vv-u,-vv])*q*np.exp(3*x)
    eps=2e-5
    z0=np.ones(3)
    finite=np.column_stack([(negative_gradient(z0-2*s)-8*negative_gradient(z0-s)
                             +8*negative_gradient(z0+s)-negative_gradient(z0+2*s))/(12*eps)
                            for s in eps*np.eye(3)])
    err=np.max(abs(cov-finite))/np.max(abs(cov))
    if err>3e-8 or min(np.linalg.eigvalsh(cov))<=0.:
        raise ArithmeticError('Independent positive thermodynamic metric required')
    return cov,err


def gain(J,start,end):
    # L L^T = metric. Stating the two norms via Cholesky removes arbitrary
    # intermediate coordinate scalings while retaining the physical endpoint metric.
    S=np.linalg.cholesky(start).T
    F=np.linalg.cholesky(end).T
    return np.linalg.svd(F@np.linalg.solve(S.T,J.T).T,compute_uv=False)


def run(eos_path,early,late,emap,lmap,composed):
    raw=[p.read_bytes() for p in (early,late,emap,lmap,composed)]
    A,B,E,L,C=map(json.loads,raw)
    h=[hashlib.sha256(z).hexdigest() for z in raw]
    eos=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if (A['table_input_sha256']!=eos or B['table_input_sha256']!=eos
        or E['input_sha256'][1]!=h[0] or L['input_sha256']!=[h[1],h[0]] or C['input_sha256']!=h[2:4]):
        raise ArithmeticError('Only the exact coherent prior histories and maps are eligible')
    u=next(r['points'] for r in A['rows'] if r['operator']=='combined')
    v=next(r['points'] for r in B['rows'] if r['operator']=='combined')
    P=CombinedExpansion(StandardModel(eos_path),PositivePortal((17,44,32)),168)
    start,e0=matrix(P,u[0])
    rows=[]
    for p,z in zip((u[-1],v[-1]),C['rows'],strict=True):
        at,e=matrix(P,p)
        w=np.array(z['fixed_total_initial_to_all_fractions'])
        direct=np.array(z['independent_finite_perturbed_composition'])
        full= gain(w,start[:2,:2],at)
        dark= gain(w[:2],start[:2,:2],at[:2,:2])
        check=gain(direct,start[:2,:2],at)
        if max(abs(full-check))>9e-6:
            raise ArithmeticError('Independently perturbed map must resolve gain in this norm')
        rows.append(dict(x_end=p[0],output_full_entropy_metric_gain=full.tolist(),output_dark_projected_gain=dark.tolist(),
                         independent_full_entropy_metric_gain=check.tolist(),end_metric=at.tolist(),
                         endpoint_gradient_fd_relative_defect=e))
    return dict(x_start=u[0][0],fields=C['fields'],start_metric=start.tolist(),start_gradient_fd_relative_defect=e0,rows=rows,
                table_input_sha256=eos,input_sha256=h,
                scope='Endpoint negative entropy Hessian weighted norm on one numerical FD history. Initial total energy fixed; final total may vary. Not a Lyapunov or physical stability theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos_path','early','late','emap','lmap','composed'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
