"""Isolate the scalar cold count term in the charge-ladder Hartree root shift."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS, load

D='2026-11-28'
X=.05
COLD=1e-5


def residual(lam,q,p,m):
    return m*m-(1e-9*X)**2-(lam/2.)*fsum(p/np.hypot(q,m))


def split_terms(q,p,z,m,mask):
    E=np.hypot(q,m)
    stable=q*q/(m*E*(E+m))
    a,b=[fsum(u[mask]) for u in (p,z)]
    c,d=[fsum((u*stable)[mask]) for u in (p,z)]
    total=fsum((z-p)/E)
    count=(b-a)*(1/m-d/b)
    shape=a*(c/a-d/b)
    warm=fsum(((z-p)/E)[~mask])
    if abs(total-count-shape-warm)>2e-15*max(fsum(p/E),fsum(z/E)):
        raise ArithmeticError('Exact cold-count, shape and outside attribution failed')
    return count,shape,warm,total


def run(root):
    src=load('fine')[0]
    entries=[]
    for label,n in zip(('','-d2','-d4','-d8'),(1,2,4,8),strict=True):
        A,h=read(root/f'early-pole-cell-endpoint-fine{label}-x05-state-{D}.json.gz')
        if (A['grid'],A['x_end'],A['subcells_per_archived_cell'],A['steps'])!=('fine',X,n,2048):
            raise ArithmeticError('Exact native fine charge ladder required')
        d=CellChart(src,SPECS['fine'],n,A['reference_x'])
        y=np.asarray(A['endpoint_occupation'])
        if A['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1.:
            raise ArithmeticError('Invalid endpoint occupations')
        entries.append((A,h,n,d,y))
    S=entries[0][3].ds
    q=np.asarray(S.p)
    mask=S.ed[1:]<=COLD
    if (np.sum(mask)==0 or np.sum(~mask)==0 or S.ed[np.sum(mask)]!=COLD
        or any(not np.array_equal(k.ds.ed,S.ed) or not np.array_equal(k.ds.w,S.w) or k.ds.lam!=S.lam for _,_,_,k,_ in entries)):
        raise ArithmeticError('Only charge resolution may change; the scalar chart must be exactly common')
    seen=[]
    for A,h,n,d,y in entries:
        m=d.state(X,y)[0]
        p=S.w*y[:d.ns]
        eps=residual(S.lam,q,p,m)
        if abs(eps)>2e-11*m*m:
            raise ArithmeticError('The independently evaluated gap equation did not close')
        seen.append((n,d,y,p,m,h))
    result=[]
    for (i,d,y,a,m,h),(j,e,z,b,u,k) in zip(seen[:-1],seen[1:],strict=True):
        count,shape,tail,net=split_terms(q,a,b,m,mask)
        sc=fsum(b[mask])/fsum(a[mask])
        new=y[:d.ns].copy()
        new[mask]*=sc
        p=S.w*new
        m_count=d.ds.gap(X,new)
        independent=brentq(lambda v:residual(S.lam,q,p,v),.7*m,1.3*m,xtol=2e-17)
        if abs(independent-m_count)>2e-14*m:
            raise ArithmeticError('Independent count-only scalar root failed')
        f0=residual(S.lam,q,p,u)
        bound=abs(f0)/(2*min(m_count,u))
        if abs(m_count-u)>bound+2e-16*m:
            raise ArithmeticError('Minimum monotone derivative should bound replacement error')
        denom=fsum(a/np.hypot(q,m))
        result.append(dict(levels=[i,j],pole_old=m,pole_full_new=u,pole_cold_count_substitution=m_count,
                           cold_count_ratio=sc,fractional_observed_pole_change=(u-m)/m,
                           relative_error_of_pole_shift=abs(m_count-u)/abs(u-m),
                           derivative_upper_bound_for_shift_error=bound/abs(u-m),
                           relative_fixed_mass_gap_terms=dict(cold_count=count/denom,cold_shape=shape/denom,
                                                            outside_cold=tail/denom,total=net/denom)))
    return dict(x=X,cold_q_upper=COLD,states_sha256_uncompressed=[v[-1] for v in seen],rows=result,
                states=[dict(charge_level=n,pole=m,cold_count=fsum(p[mask]),whole_scalar_count=fsum(p))
                        for n,d,y,p,m,h in seen],
                scope='Fixed-snapshot Hartree root attribution on the exact shared scalar chart; not an attribution of time evolution, a source-freeze error, or grid convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
