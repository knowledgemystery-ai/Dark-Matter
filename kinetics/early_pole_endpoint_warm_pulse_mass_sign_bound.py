"""Enclose the nodal pressure sign continuously between the archived x04 and x05 masses."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_warm_propagated_closure_contrast import run as replay

DATE='2026-11-28'


def bounds(q,d,a,b):
    """Use monotone positive weights, handling negative occupation differences separately."""
    if not 0.<a<=b or q.shape!=d.shape or min(q)<0 or not np.all(np.isfinite(q+d)):
        raise ArithmeticError('Require a positive mass interval and finite nodal inputs')
    u=q*q/np.hypot(q,a)
    v=q*q/np.hypot(q,b)
    lo=d*np.where(d>=0,v,u)
    hi=d*np.where(d>=0,u,v)
    # Conservative pad for this normal-range float evaluation, not an uncertainty
    # budget for the occupations, their archived differences or the integrator.
    scale=fsum(abs(d*u))
    pad=64*np.finfo(float).eps*scale
    return float(np.nextafter(fsum(lo)-pad,-np.inf)),float(np.nextafter(fsum(hi)+pad,np.inf))


def cover(q,d,a,b):
    # Measure the actual contrast at sampled points even if a coarse enclosure
    # initially overlaps zero. Do not turn lack of certification into a zero.
    a0,b0=bounds(q,d,a,a)
    c0,d0=bounds(q,d,b,b)
    first=int(a0>0.)-int(b0<0.)
    last=int(c0>0.)-int(d0<0.)
    result=dict(boundary_signs=[first,last],complete=False)
    if first==0 or last==0 or first!=last:
        return result
    stack=[(a,b,0)]
    out=[]
    while stack:
        aa,bb,level=stack.pop()
        lo,hi=bounds(q,d,aa,bb)
        if (first==1 and lo>0.) or (first==-1 and hi<0.):
            out.append((lo,hi,level))
        elif level>=14:
            return dict(result,unresolved_mass_interval=[aa,bb],unresolved_difference_interval=[lo,hi])
        else:
            mid=(aa+bb)/2
            stack.extend(((mid,bb,level+1),(aa,mid,level+1)))
    return dict(result,complete=True,mass_panel_count=len(out),max_depth=max(o[2] for o in out),
                difference_enclosure=[min(y[0] for y in out),max(y[1] for y in out)])


def run(root):
    bs=[(root/f'early-pole-warm-checkpoint-clock-comparison-x{k}-{DATE}.json').read_bytes() for k in ('04','05')]
    A,B=map(json.loads,bs)
    if (A['x'],B['x'])!=(.04,.05) or A['fields']!=B['fields']:
        raise ArithmeticError('Matched archived pressure mass range required')
    a,b=A['mass'],B['mass']
    rows=[]
    def record(x,t,q,d,normal):
        z=cover(q,d,a,b)
        # Normalize only to the previously measured positive total pressure.
        if z['complete']:
            z['difference_over_fixed_reference_pressure']=[v/normal for v in z.pop('difference_enclosure')]
        rows.append(dict(x_insertion=x,x_measurement=t,**z))
    W=replay(root,observer=record)
    if any(x['common_pressure_mass']!=b for x in W['rows']) or len(rows)!=6:
        raise ArithmeticError('Expected three exact held arms and both measurement times')
    return dict(archived_mass_interval=[a,b],mass_sha256=[hashlib.sha256(x).hexdigest() for x in bs],
                input_sha256_uncompressed=W['input_uncompressed_sha256'],rows=rows,
                scope='Continuous mass-weight check of the alternate nodal pressure for fixed archived occupations. No rerun at modified mass, no assertion about reconstructed pressure at other masses, and no physical uncertainty bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
