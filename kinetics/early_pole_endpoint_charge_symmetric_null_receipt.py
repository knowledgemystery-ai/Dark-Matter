"""Audit the finite symmetric even/odd cold response around the nominal remap."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_charge_symmetric_null_terminal import run as replay,setup


def run(root,remap,records,payloads):
    raw=remap.read_bytes()
    A=json.loads(raw)['results']
    d,_,_,_=setup(root)
    result=[]
    for sign,r,p in zip((-1,1),records,payloads,strict=True):
        a,b=r.read_bytes(),p.read_bytes()
        u,v=json.loads(a),json.loads(b)
        expected=replay(root,sign,p,setup_only=True)
        if (any(u[k]!=q for k,q in expected.items()) or v['completed']!=2048 or v['checks']!=u['checks']
            or abs(v['x_reached']-.05)>2e-17 or len(v['checks'])!=8
            or v['payload_sha256']!=drive_checkpoint_digest(v) or u['payload_sha256']!=v['payload_sha256']
            or max(abs(q['continuity_residual']) for q in v['checks'])>3e-10):
            raise ArithmeticError('Require exact complete symmetric endpoints and restart setup')
        y=np.asarray(v['occupation'])
        if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0. or max(y[d.ns:])>1.:
            raise ArithmeticError('The exact archived endpoint must respect Pauli and positivity')
        z=d.ds.w*y[:d.ns]
        p=dict(cold_count=fsum(z[d.ds.ed[1:]<=1e-5]),whole_scalar_count=fsum(z),pole=d.state(.05,y)[0])
        if any(p[k]!=u[f] for k,f in (('cold_count','final_cold_count'),('whole_scalar_count','final_whole_scalar_count'),('pole','final_scalar_pole'))):
            raise ArithmeticError('Independent endpoint observable replay failed')
        result.append((p,hashlib.sha256(a).hexdigest(),u['initial_log_rate']['clock_net_forward_inverse'][0],
                       u['nominal_initial_log_rate']['clock_net_forward_inverse'][0]))
    J={}
    for k,o in A.items():
        ref=o['grid4_donor2']
        a,b=[q[0][k] for q in result]
        J[k]=dict(endpoint_low=a,nominal=ref,endpoint_high=b,finite_odd=(b-a)/2.,finite_even=(b+a)/2.-ref,
                  largest_change_over_nominal_grid_replacement=max(abs(a-ref),abs(b-ref))/abs(o['average_grid_and_representation_replacement']))
    a,b=[q[2] for q in result]; ref=result[0][3]
    if ref!=result[1][3]:
        raise ArithmeticError('Common frozen reference required')
    return dict(results=J,initial_clock=dict(low=a,nominal=ref,high=b,finite_odd=(b-a)/2.,finite_even=(b+a)/2.-ref),
                raw_receipt_sha256=[q[1] for q in result],nominal_remap_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Finite even and odd contrasts of truly equal and opposite moment-null displacements; not an extrapolation to infinitesimal variation or other charge histories.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('remap',type=Path)
    for key in ('records','payloads'):
        p.add_argument('--'+key,nargs=2,type=Path,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
