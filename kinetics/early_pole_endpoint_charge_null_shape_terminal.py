"""Evolve the two prespecified 5%-of-feasible two-moment-null charge-shape alternatives."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_nested_moment_remap import X,prepare
from early_pole_endpoint_charge_null_face_clock import adjust


def run(root,clock,sign,checkpoint):
    raw=clock.read_bytes()
    a=json.loads(raw)
    d,y,P=prepare(root,2)
    if a['x']!=X or a['destination_charge_grid']!=4 or sign not in (-1,1):
        raise ArithmeticError('Exact D2 charged donor to D4 grid at fixed D2 scalar seed required')
    old=d.state(X,y)
    f=.05*sign
    z0=y.copy()
    err=[]
    for which,k,m in (('d',slice(d.ns,d.iP),X*d.ratio),('p',slice(d.iP,d.n),X)):
        p,w,e=[getattr(d,n) for n in ('q'+which,'w'+which,which+'_edges')]
        z0[k],epsilon=adjust(p,w,y[k],m,e[::2],f)
        err.append(epsilon)
    e=d.state(X,z0)
    if e[0]!=old[0] or not np.array_equal(z0[:d.ns],y[:d.ns]) or abs(e[2]-old[2])>5e-15:
        raise ArithmeticError('Alternate remap must preserve scalar seed and moment-defined total energy')
    h=[hashlib.sha256(z).hexdigest() for z in (raw,Path(__file__).read_bytes(),Path(adjust.__code__.co_filename).read_bytes())]
    key=hashlib.sha256((P['provenance']+''.join(h)+str(sign)).encode()).hexdigest()
    prev=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    checks,z=drive(d,.05,2048,x_start=X,initial=z0,stepper=step,restart=prev,
                   checkpoint=lambda x:write_drive_checkpoint(checkpoint,x),checkpoint_stride=128,provenance=key)
    t=json.loads(checkpoint.read_bytes())
    if t['completed']!=2048 or t['payload_sha256']!=drive_checkpoint_digest(t) or t['checks']!=checks or max(abs(u['continuity_residual']) for u in checks)>3e-10:
        raise ArithmeticError('Complete actual alternative charge-shape propagation required')
    a=d.ds.w*z[:d.ns]
    w=d.ds.ed[1:]<=1e-5
    return dict(sign=sign,fraction_of_feasible_move=f,setup=P,
                maximum_relative_moment_defect=max(err),source_inputs=h,provenance=key,
                final_count=fsum(a),final_cold_count=fsum(a[w]),final_scalar_pole=d.state(.05,z)[0],
                payload_sha256=t['payload_sha256'],checks=checks,
                scope='Prespecified equal-charge-moment shape perturbation of one remap. No all-shape robustness, physical rate, or convergence conclusion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('clock',type=Path)
    p.add_argument('sign',type=int,choices=(-1,1))
    p.add_argument('checkpoint',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
