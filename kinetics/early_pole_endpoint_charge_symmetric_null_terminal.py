"""Use a truly symmetric cell-moment-null change, not equal fractions of asymmetric ranges."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_nested_moment_remap import X,local_q,prepare
from early_pole_endpoint_crossed_live_clock import case


def balanced_nodes(q,w,F,m,edges):
    # Return a node-wise signed change for 5% of the smaller of the two
    # feasible absolute center-mass displacements in every common cell.
    delta=np.zeros_like(F)
    for i,(l,r) in enumerate(zip(edges[:-1],edges[1:],strict=True)):
        b=slice(4*i,4*i+4)
        e=local_q(q[b],l,r,m)
        if e[0]!=0. or e[-1]!=1. or not 0.<e[1]==e[2]<1.:
            raise ArithmeticError('Exact nested two-endpoint subcell geometry required')
        p=e[1]; w0=w[b]; f=F[b]
        a=w0*f
        N,K=fsum(a),fsum(a*e)
        h=a[1]+a[2]
        c=w0[1]+w0[2]
        low=max(0.,(N-K-w0[0])/(1-p),(K-w0[3])/p)
        high=min(c,(N-K)/(1-p),K/p)
        rad=min(h-low,high-h)
        if rad < -2e-14*max(N,1e-300):
            raise ArithmeticError('Nominal occupation not in its pair budget')
        z=.05*max(0.,rad)
        v=np.array([-(1-p)/w0[0],1./c,1./c,-p/w0[3]])*z
        err=max(abs(fsum(v*w0)),abs(fsum(v*w0*e)))
        if err>1e-14*max(N,1e-300) or min(f-v)<0. or min(f+v)<0. or max(f+v)>1. or max(f-v)>1.:
            raise ArithmeticError('Symmetric move fails exact local moments or Pauli caps')
        delta[b]=v
    return delta


def setup(root):
    d,y,p=prepare(root,2)
    b=np.zeros_like(y)
    for name,sector,mass in (('d',slice(d.ns,d.iP),X*d.ratio),('p',slice(d.iP,d.n),X)):
        q,w,e=[getattr(d,u) for u in ('q'+name,'w'+name,name+'_edges')]
        b[sector]=balanced_nodes(q,w,y[sector],mass,e[::2])
    for s in (-1,1):
        a=d.state(X,y+s*b)
        v=d.state(X,y)
        if a[0]!=v[0] or abs(a[2]-v[2])>6e-15:
            raise ArithmeticError('Independently require the same scalar root and total energy at setup')
    return d,y,b,p


def run(root,sign,checkpoint,setup_only=False):
    d,y,b,data=setup(root)
    signed=y+sign*b
    args=dict(direction=sign,x_start=X,x_end=.05,steps=2048,initial_log_rate=case(d,signed,X),
              nominal_initial_log_rate=case(d,y,X),setup=data,
              declared_max_abs_occupation_shift=float(max(abs(b))))
    if setup_only:
        return args
    key=hashlib.sha256((data['provenance']+hashlib.sha256(Path(__file__).read_bytes()).hexdigest()+str(sign)).encode()).hexdigest()
    prev=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    checks,a=drive(d,.05,2048,x_start=X,initial=signed,stepper=step,restart=prev,
                   checkpoint=lambda x:write_drive_checkpoint(checkpoint,x),checkpoint_stride=128,provenance=key)
    t=json.loads(checkpoint.read_bytes())
    if t['completed']!=2048 or t['checks']!=checks or t['payload_sha256']!=drive_checkpoint_digest(t) or max(abs(z['continuity_residual']) for z in checks)>3e-10:
        raise ArithmeticError('Only complete admitted actual symmetric terminal trajectories')
    p=d.ds.w*a[:d.ns]
    return dict(**args,provenance=key,checks=checks,payload_sha256=t['payload_sha256'],
                final_cold_count=fsum(p[d.ds.ed[1:]<=1e-5]),final_whole_scalar_count=fsum(p),final_scalar_pole=d.state(.05,a)[0],
                scope='One prespecified pair of exactly symmetric feasible perturbations on a held remap, not a physical bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('sign',type=int,choices=(-1,1))
    p.add_argument('checkpoint',type=Path)
    p.add_argument('--setup-only',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
