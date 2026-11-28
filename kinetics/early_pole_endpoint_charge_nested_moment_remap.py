"""Two-moment Pauli-respecting transfer on the nested native D2/D4 charged grids."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import expit,logit

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,load

D='2026-11-28'
X=.04


def capacity_extreme(w,q,n):
    n0=n
    v=[]
    for a,x in zip(w,q,strict=True):
        z=min(n,a)
        v.append(z*x)
        n-=z
    if n>3e-14*n0:
        raise ArithmeticError('Requested cell exceeds Pauli capacity')
    return fsum(v)/n0


def fit(w,q,n,mean):
    if (len(w)!=len(q) or min(w)<=0 or not np.all(np.diff(q)>=0.)
        or not 0.<=n<fsum(w) or not 0.<=mean<=1.):
        raise ArithmeticError('Invalid positive normalized cell constraints')
    if n==0.:
        return np.zeros_like(w)
    if len(w)==2 and np.array_equal(q,[0.,1.]):
        z=n*np.array([1.-mean,mean])/w
        if max(z)>1. or min(z)<0.:
            raise ArithmeticError('No feasible coarse two-moment endpoint projection')
        return z
    lo,hi=capacity_extreme(w,q,n),capacity_extreme(w[::-1],q[::-1],n)
    if not lo<mean<hi:
        raise ArithmeticError(('Need strictly interior entropy-fit moments',lo,mean,hi))
    target=logit(n/fsum(w))
    def at(b):
        a=brentq(lambda a:fsum(w*expit(a+b*q))/n-1.,target-max(0.,b)-8.,target-min(0.,b)+8.,xtol=2e-13)
        z=expit(a+b*q)
        return z,fsum(w*z*q)/n
    B=4.
    while B<=1024. and not at(-B)[1]<mean<at(B)[1]:
        B*=2.
    if B>1024.:
        raise ArithmeticError('Interior moment fit unresolved at finite tilt')
    b=brentq(lambda p:at(p)[1]-mean,-B,B,xtol=5e-13)
    return at(b)[0]


def local_q(p,l,r,m):
    E0,E1=np.hypot([l,r],m)
    e=np.hypot(p,m)
    x=((p-l)*(p+l)/(r-l)/(r+l))*(E1+E0)/(e+E0)
    return np.where(p==l,0.,np.where(p==r,1.,x))


def transfer(S,T,F,m):
    es,ps,ws=S
    et,pt,wt=T
    panels=min(len(es),len(et))-1
    a,b=(len(es)-1)//panels,(len(et)-1)//panels
    if (max(a,b)!=2 or min(a,b)!=1 or len(F)!=len(ps)
        or not np.allclose(es[::a],et[::b],rtol=0.,atol=2e-17)):
        raise ArithmeticError('Need a strictly nested pair of native charged grids')
    E=es[::a]
    Y=np.zeros_like(pt)
    err=[]
    for i,(l,r) in enumerate(zip(E[:-1],E[1:],strict=True)):
        A=slice(2*a*i,2*a*(i+1))
        B=slice(2*b*i,2*b*(i+1))
        p,t=[local_q(z,l,r,m) for z in (ps[A],pt[B])]
        z=ws[A]*F[A]
        n=fsum(z)
        mean=fsum(z*p)/n if n else 0.
        v=fit(wt[B],t,n,mean)
        out=wt[B]*v
        Y[B]=v
        err.append(max(abs(fsum(out)-n),abs(fsum(out*t)-n*mean))/max(n,1e-300))
    if max(err)>1e-11 or min(Y)<0. or max(Y)>1.:
        raise ArithmeticError('Per-coarse-cell count-energy and Pauli closure failed')
    old=ws*F; new=wt*Y
    return Y,dict(max_cell_moment_relative_error=max(err),
                  normalized_total_number_error=abs(fsum(new)/fsum(old)-1.),
                  normalized_total_true_energy_error=abs(fsum(new*np.hypot(pt,m))/fsum(old*np.hypot(ps,m))-1.))


def prepare(root,source):
    dest=6-source
    base=load('fine')[0]
    Ds=[CellChart(base,SPECS['fine'],n,.005) for n in (source,dest)]
    rows=[read(root/f'early-pole-cell-endpoint-fine-d{n}-x04-state-{D}.json.gz') for n in (source,dest)]
    (A,h),(B,g)=rows
    N,k=read(root/f'early-pole-cell-endpoint-fine-d2-x04-state-{D}.json.gz')
    for R,n,J in zip((A,B),(source,dest),Ds,strict=True):
        if (R['grid'],R['subcells_per_archived_cell'],R['counts'],R['x_end'],R['steps'],R['reference_x'])!=('fine',n,[J.ns,J.nf,J.np],X,2048,.005):
            raise ArithmeticError('Exact .04 native terminal-start cards required')
    if k!=(h if source==2 else g) or not np.array_equal(Ds[0].ds.ed,Ds[1].ds.ed):
        raise ArithmeticError('Both branches need the same native D2 scalar state')
    V=np.asarray(A['endpoint_occupation'])
    scalar=np.asarray(N['endpoint_occupation'])[:Ds[0].ns]
    src=np.r_[scalar,V[Ds[0].ns:]]
    ch=[]
    for name,field in (('daughter','d'),('parent','p')):
        C=[(getattr(J,f'{field}_edges'),getattr(J,f'q{field}'),getattr(J,f'w{field}')) for J in Ds]
        a=V[Ds[0].ns:Ds[0].iP] if field=='d' else V[Ds[0].iP:]
        u,z=transfer(*C,a,X*Ds[0].ratio if field=='d' else X)
        ch.append((u,z))
    tgt=np.r_[scalar,ch[0][0],ch[1][0]]
    old,new=[J.state(X,y) for J,y in zip(Ds,(src,tgt),strict=True)]
    scale=max(abs(old[2]),1e-50)
    if old[0]!=new[0] or abs(old[2]-new[2])/scale>3e-14:
        raise ArithmeticError('Total energy and exact shared scalar root must replay at remap')
    code=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    P=dict(source_charge_grid=source,destination_charge_grid=dest,held_scalar_donor=2,x_start=X,x_end=.05,steps=2048,
           scalar_pole=new[0],total_energy_before_after=[old[2],new[2]],remap_errors=[z[1] for z in ch],
           input_state_hashes=[h,g,k],source_sha256=code)
    key=hashlib.sha256(json.dumps(P,separators=(',',':')).encode()).hexdigest()
    return Ds[1],tgt,dict(**P,provenance=key)


def run(root,source,checkpoint,setup_only=False):
    chart,y,doc=prepare(root,source)
    if setup_only:
        return doc
    P=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    checks,u=drive(chart,.05,2048,x_start=X,initial=y,stepper=step,restart=P,checkpoint=lambda x:write_drive_checkpoint(checkpoint,x),
                   checkpoint_stride=128,provenance=doc['provenance'])
    R=json.loads(checkpoint.read_bytes())
    if R['completed']!=2048 or R['payload_sha256']!=drive_checkpoint_digest(R) or R['checks']!=checks or max(abs(q['continuity_residual']) for q in checks)>3e-10:
        raise ArithmeticError('Only a completed conservative alternate charge-grid history eligible')
    m=chart.ds.ed[1:]<=1e-5
    z=chart.ds.w*u[:chart.ns]
    return dict(**doc,cold_count=fsum(z[m]),scalar_count=fsum(z),scalar_pole_final=chart.state(.05,u)[0],
                final_payload_sha256=R['payload_sha256'],checks=checks,
                scope='One held two-moment remap, not a unique distribution or physical convergence bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('source',type=int,choices=(2,4))
    p.add_argument('checkpoint',type=Path)
    p.add_argument('--setup-only',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
