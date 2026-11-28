"""Bound how many warm cells would need exact shape to preserve a four-moment phase sign."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_streaming_pseudopressure_constraint import cell
from early_pole_endpoint_warm_streaming_pin_budget import solve_budget
from early_pole_endpoint_warm_three_chart_stress_hull import conditioned


def run(native,pressure,states):
    (u,uh)=read(native)
    raw=pressure.read_bytes()
    D=json.loads(raw)
    V=[read(f) for f in states]
    if D['input_sha256'][0]!=uh or [h for _,h in V]!=D['input_sha256'][-3:-1]:
        raise ArithmeticError('Need exact G8/G9 baseline used to build the pressure-conditioned receipt')
    d,y=unpack(u,False,.05)
    m=d.state(.05,y)[0]
    if D['formal_mass']!=m:
        raise ArithmeticError('Reference mass changed')
    tables=[conditioned(d,q,m,.05,None)[1] for q,_ in V]
    total=[fsum(z[2] for z in q) for q in tables]
    out=[]
    for R in D['rows']:
        A,p=R['A'],R['phase']
        nodal=[k['nodal'] for k in R['charts'][:2]]
        signed=nodal[1]-nodal[0]
        origin,target=(0,1) if signed>0. else (1,0)
        margin=abs(signed)
        slots=[]; L=[]; U=[]
        for i,pair in enumerate(zip(*tables,strict=True)):
            if pair[0][:2]!=pair[1][:2] or pair[0][-1]!=pair[1][-1]:
                raise ArithmeticError('Warm coordinate split differs')
            if not pair[0][-1]:
                continue
            values=[]
            for q,N in zip(pair,total,strict=True):
                a,b,z,_,x,w,_=q
                if z:
                    O,I,P,*_=cell(a,b,m,A,p,x,w/z,1536)
                    values.append((O*z/N,I*z/N,P*z/N))
                else:
                    values.append((np.zeros(2),np.zeros(2),0.))
            a,b=[values[i] for i in (origin,target)]
            low=max(a[1][1],a[2])-a[2]+b[2]-min(b[1][0],b[2])
            high=max(a[0][1],a[2])-a[2]+b[2]-min(b[0][0],b[2])
            if low<0. or high<low:
                raise ArithmeticError('Invalid additive uncertainty loss')
            L.append(low); U.append(high)
            slots.append(dict(cell=i,q_bounds=list(pair[0][:2]),constructive_opposite_loss=low,bounded_opposite_loss=high))
        r,t,Lreq,_=solve_budget(L,margin)
        s,u,_,Usuff=solve_budget(U,margin)
        if Lreq>Usuff:
            raise ArithmeticError('Incompatible cardinality bounds')
        out.append(dict(A=A,phase=p,nodal_G9_minus_G8=signed,n_warm_cells=len(slots),
                        constructive_required_at_least=Lreq,outer_sufficient_at_most=Usuff,
                        worst_margin_no_pins=[margin-t[0],margin-u[0]],
                        guaranteed_margin_at_sufficient_pins=margin-u[Usuff],
                        outer_prioritized_pins=[slots[i]['cell'] for i in s[:Usuff]],
                        constructive_loss_rank=[slots[i]['cell'] for i in r[:max(Lreq,Usuff)]],
                        cells=slots))
    return dict(x=.05,formal_mass=m,warm_grid_subintervals=1536,rows=out,input_sha256=[uh,hashlib.sha256(raw).hexdigest(),*[x for _,x in V]],
                scope='Counterfactual cardinality budget: fix both charts’ original shape in a selected common warm cell; cold already fixed. Chart-local pressure and pseudo-pressure moments fixed elsewhere. No calibration or physical transfer claim.')


if __name__=='__main__':
    c=argparse.ArgumentParser(description=__doc__)
    for a in ('native','pressure'):
        c.add_argument(a,type=Path)
    c.add_argument('states',type=Path,nargs=2)
    print(json.dumps(run(**vars(c.parse_args())),indent=2,allow_nan=False))
