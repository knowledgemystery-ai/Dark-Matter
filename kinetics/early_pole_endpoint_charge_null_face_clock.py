"""Vary the remaining charged-node allocation while preserving each cell's count and energy."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_charge_nested_moment_remap import X,local_q,prepare
from early_pole_endpoint_crossed_live_clock import case
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def adjust(q,w,F,M,span,fraction):
    out=np.zeros_like(F)
    max_def=0.
    for i,(l,r) in enumerate(zip(span[:-1],span[1:],strict=True)):
        p=slice(4*i,4*i+4)
        e=local_q(q[p],l,r,M)
        if e[0]!=0. or e[-1]!=1. or not 0.<e[1]==e[2]<1.:
            raise ArithmeticError('Expected two shared interior nodes on the exact fine subdivision')
        d=e[1]
        m=w[p]*F[p]
        N=fsum(m)
        K=fsum(e*m)
        c=w[p]
        C=c[1]+c[2]
        lo=max(0.,(N-K-c[0])/(1.-d),(K-c[3])/d)
        hi=min(C,(N-K)/(1.-d),K/d)
        old=m[1]+m[2]
        if lo>old+2e-14*N or hi<old-2e-14*N:
            raise ArithmeticError('Nominal split not inside the exact moment-preserving range')
        mid=old+abs(fraction)*((hi if fraction>0. else lo)-old)
        val=np.array([N-K-(1.-d)*mid,mid*c[1]/C,mid*c[2]/C,K-d*mid])
        y=val/c
        if min(y)<0. or max(y)>1.:
            raise ArithmeticError('Moment-preserving alternative left Pauli range')
        max_def=max(max_def,abs(fsum(val)-N)/max(N,1e-300),abs(fsum(val*e)-K)/max(N,1e-300))
        out[p]=y
    if max_def>3e-14:
        raise ArithmeticError('Count and centered energy should be preserved per shared cell')
    return out,max_def


def run(root,setup_card):
    d,y,P=prepare(root,2)
    H=setup_card.read_bytes()
    if P!=json.loads(H):
        raise ArithmeticError('Exact declared up-remap setup required')
    native,sha=read(root/'early-pole-cell-endpoint-fine-d4-x04-state-2026-11-28.json.gz')
    if sha!=P['input_state_hashes'][1]:
        raise ArithmeticError('Expected exact native D4 charges')
    O=np.asarray(native['endpoint_occupation'])
    Q=np.r_[y[:d.ns],O[d.ns:]]
    ref=case(d,Q,X)
    base=case(d,y,X)
    out=[]
    for f in (-.05,+.05):
        z=y.copy()
        marks=[]
        for which,k,m in (('d',slice(d.ns,d.iP),X*d.ratio),('p',slice(d.iP,d.n),X)):
            p,w,e=[getattr(d,n) for n in ('q'+which,'w'+which,which+'_edges')]
            z[k],err=adjust(p,w,y[k],m,e[::2],f)
            marks.append(err)
        a=case(d,z,X)
        out.append(dict(signed_fraction_of_feasible_move=f,max_local_moment_error=max(marks),
                        clock=a,net_change_from_nominal=a['clock_net_forward_inverse'][0]-base['clock_net_forward_inverse'][0],
                        whole_state_energy_difference=d.state(X,z)[2]-d.state(X,y)[2]))
    return dict(rows=out,held_scalar_donor=2,destination_charge_grid=4,source_charge_donor=2,x=X,
                nominal_two_moment_remap=base,matched_native_D4_charge_counterfactual=ref,
                setup_sha256=hashlib.sha256(H).hexdigest(),
                scope='Frozen chart/Pauli-preserving shape control for this chosen 5% moment-null move only; not a terminal evolution or unique higher-moment reconstruction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('setup_card',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
