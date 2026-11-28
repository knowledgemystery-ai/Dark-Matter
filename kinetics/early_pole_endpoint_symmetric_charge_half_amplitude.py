"""Test half amplitude at the same 4096-step clock and exact remap."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_symmetric_null_terminal import X,setup
from early_pole_endpoint_crossed_live_clock import case


def run(root,prior,sign,checkpoint,setup_only=False):
    prev=prior.read_bytes()
    P=json.loads(prev)
    if P['mesh_counts']!=[2048,4096] or sign not in (-1,1):
        raise ArithmeticError('Exact original 4096 amplitude study required')
    d,y,b,pro=setup(root)
    start=y+sign*b/2.
    x0,x1=d.state(X,y),d.state(X,start)
    if (not np.array_equal(start[:d.ns],y[:d.ns]) or x0[0]!=x1[0] or abs(x0[2]-x1[2])>6e-15
        or min(start)<0. or max(start[d.ns:])>1.):
        raise ArithmeticError('Half-amplitude restart must preserve scalar seed, energy and Pauli positivity')
    H=[hashlib.sha256(s).hexdigest() for s in (prev,Path(__file__).read_bytes(),Path(setup.__code__.co_filename).read_bytes())]
    key=hashlib.sha256((pro['provenance']+''.join(H)+str(sign)).encode()).hexdigest()
    meta=dict(sign=sign,amplitude_ratio=0.5,steps=4096,setup=pro,prior_and_source_hashes=H,provenance=key,
              initial_clock=case(d,start,X),initial_nominal_clock=case(d,y,X))
    if setup_only:
        return meta
    old=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    c,f=drive(d,.05,4096,x_start=X,initial=start,stepper=step,restart=old,
              checkpoint=lambda p:write_drive_checkpoint(checkpoint,p),checkpoint_stride=128,provenance=key)
    C=json.loads(checkpoint.read_bytes())
    if (C['completed']!=4096 or C['checks']!=c or C['payload_sha256']!=drive_checkpoint_digest(C)
        or len(c)!=8 or abs(C['x_reached']-.05)>2e-17 or max(abs(a['continuity_residual']) for a in c)>3e-10):
        raise ArithmeticError('Full same-clock advanced numerical state required')
    q=d.ds.w*f[:d.ns]
    return dict(**meta,checks=c,final_cold=fsum(q[d.ds.ed[1:]<=1e-5]),final_total=fsum(q),
                final_pole=d.state(.05,f)[0],payload_sha256=C['payload_sha256'],
                scope='Finite factor-of-two amplitude check for one moment-null shape direction, not a physical or asymptotic error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','prior'):
        p.add_argument(k,type=Path)
    p.add_argument('sign',type=int,choices=(-1,1))
    p.add_argument('checkpoint',type=Path)
    p.add_argument('--setup-only',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
