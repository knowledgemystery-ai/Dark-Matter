"""Reintegrate the exact nominal and equal-opposite charge shapes at twice the terminal resolution."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_symmetric_null_terminal import X,setup


def run(root,prior,sign,checkpoint):
    a=prior.read_bytes()
    T=json.loads(a)
    if sign not in (-1,0,1) or tuple(T['initial_clock'])!=('low','nominal','high','finite_odd','finite_even'):
        raise ArithmeticError('Need the accepted original exactly symmetric setup')
    d,y,b,S=setup(root)
    y0=y+sign*b
    if not np.array_equal(y0[:d.ns],y[:d.ns]):
        raise ArithmeticError('Scalar seed must be held bitwise')
    H=[hashlib.sha256(z).hexdigest() for z in (a,Path(__file__).read_bytes(),Path(setup.__code__.co_filename).read_bytes())]
    key=hashlib.sha256((S['provenance']+''.join(H)+str(sign)).encode()).hexdigest()
    old=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    obs,f=drive(d,.05,4096,x_start=X,initial=y0,stepper=step,restart=old,
                checkpoint=lambda p:write_drive_checkpoint(checkpoint,p),checkpoint_stride=128,provenance=key)
    d0=json.loads(checkpoint.read_bytes())
    if (d0['completed']!=4096 or d0['payload_sha256']!=drive_checkpoint_digest(d0) or d0['checks']!=obs
        or len(obs)!=8 or max(abs(p['continuity_residual']) for p in obs)>3e-10):
        raise ArithmeticError('Validated complete high-resolution continuation required')
    a=d.ds.w*f[:d.ns]
    mask=d.ds.ed[1:]<=1e-5
    return dict(sign=sign,x_window=[X,.05],steps=4096,checks=obs,final_cold=fsum(a[mask]),final_total=fsum(a),
                final_pole=d.state(.05,f)[0],setup=S,provenance=key,source_hashes=H,
                complete_payload_sha256=d0['payload_sha256'],
                scope='Exact 4096-step time refinement of this finite symmetric charge family, not evidence of spatial convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('prior',type=Path)
    p.add_argument('sign',choices=(-1,0,1),type=int)
    p.add_argument('checkpoint',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
