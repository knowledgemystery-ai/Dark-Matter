"""Audit serialized pressure-preserving sign witnesses without an optimizer."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_scalar_gauss2_energy_reconstruction import scaled_energy,shell_rule
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_three_chart_stress_hull import conditioned
from early_pole_endpoint_warm_three_chart_streaming_hull import kernel
from early_pole_endpoint_warm_streaming_pressure_constraint import feature


def run(native,witness,states):
    u,uh=read(native); raw=witness.read_bytes(); v=json.loads(raw)
    C=[read(p) for p in states]
    if uh!=v['input_sha256'][0] or [h for _,h in C]!=v['input_sha256'][-3:-1]:
        raise ArithmeticError('Witness must bind exact G8/G9 x=.05 states')
    base,y=unpack(u,False,.05)
    m=base.state(.05,y)[0]
    if m!=v['formal_mass']:
        raise ArithmeticError('Witness and common reference mass disagree')
    R=next(k for k in v['rows'] if k['A']==10.)
    out=[]
    for (c,_),row in zip(C,R['charts'][:2],strict=True):
        tag,cells,_=conditioned(base,c,m,.05,None)
        by={x['cell']:x for x in row['positive_warm_atoms']}
        vals=[]; defect=0.; E=0.
        for j,(a,b,N,_,q,w,warm) in enumerate(cells):
            if not N:
                continue
            E+=N
            if not warm:
                X,_,L=shell_rule(a,b,m,256)
                W=np.exp(L)
                cold=fsum(W*kernel(X,m,10.,R['phase']))/fsum(W)
                vals.append([N*cold,N*cold])
                continue
            z=by.pop(j)
            if z['bounds']!=[a,b] or z['count']!=N:
                raise ArithmeticError('Witness cannot move probability between cells')
            V=[]
            E0=np.hypot(a,m); H=(b-a)*(b+a)/(np.hypot(b,m)+E0)
            target=feature(scaled_energy(q,a,b,m),E0,H)@(w/N)
            for h in (z['lower'],z['upper']):
                X,W=[np.asarray(h[k]) for k in ('q','weights')]
                if not 1<=len(X)<=3 or min(W)<0. or np.any(X<a*(1-2e-15)) or np.any(X>b*(1+2e-15)):
                    raise ArithmeticError('Infeasible positive support')
                t=scaled_energy(X,a,b,m)
                error=np.max(abs(feature(t,E0,H)@W-target))
                defect=max(defect,error)
                V.append(N*fsum(W*kernel(X,m,10.,R['phase'])))
            vals.append(V)
        result=np.array([fsum(v) for v in zip(*vals,strict=True)])/E
        if by or defect>2e-13 or np.max(abs(result-row['constructive']))>2e-14:
            raise ArithmeticError('Independent archived positive supports did not replay')
        out.append(dict(chart=tag,values=result.tolist(),maximum_held_moment_defect=defect))
    A,B=(x['values'] for x in out)
    delta=[B[0]-A[1],B[1]-A[0]]
    if not delta[0]<-1e-12<1e-12<delta[1]:
        raise ArithmeticError('Both signs are not independently resolved')
    return dict(A=10.,phase=R['phase'],rows=out,G9_minus_G8=delta,witness_sha256=hashlib.sha256(raw).hexdigest())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('native',type=Path)
    p.add_argument('witness',type=Path)
    p.add_argument('states',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
