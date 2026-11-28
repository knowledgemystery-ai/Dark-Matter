"""Hold native scalar shape while giving the cold group its crossed donor count."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_crossed_charge_terminal_history import prepare
from early_pole_endpoint_scalar_subshell_chain_receipt import read

D='2026-11-28'


def run(root,history,full_cross,charge,checkpoint,setup_only=False):
    scalar=6-charge
    d,v,mask,setup=prepare(root,history,charge,scalar)
    doc=full_cross.read_bytes()
    C=json.loads(doc)
    if any(C[k]!=a for k,a in setup.items()):
        raise ArithmeticError('Same complete full-scalar cross required')
    y,h=read(root/f'early-pole-cell-endpoint-fine-d{charge}-x04-state-{D}.json.gz')
    if h!=setup['canonical_input_states_sha256'][0]:
        raise ArithmeticError('Exact native scalar shape required')
    v0=np.asarray(y['endpoint_occupation'])
    query=v0.copy()
    target=fsum((d.ds.w*v[:d.ns])[mask])
    norm=fsum((d.ds.w*v0[:d.ns])[mask])
    ix=np.flatnonzero(mask)
    scale=target/norm
    query[ix]*=scale
    if (not np.array_equal(query[d.ns:],v[d.ns:]) or not np.array_equal(query[:d.ns][~mask],v0[:d.ns][~mask])
        or abs(fsum((d.ds.w*query[:d.ns])[mask])-target)>3e-17):
        raise ArithmeticError('Cold count, native shape and frozen charged occupation must match')
    edge=d.state(.04,query)
    kk=hashlib.sha256(doc).hexdigest()
    sk=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    pr=hashlib.sha256((setup['provenance']+kk+sk).encode()).hexdigest()
    report=dict(charge_donor=charge,cold_scalar_count_donor=scalar,x_start=.04,x_end=.05,steps=2048,
                native_initial_cold_count=norm,intervened_initial_cold_count=target,initial_cold_rescale=scale,
                initial_mass=edge[0],initial_energy=edge[2],initial_H=edge[3],
                native_initial_sha256=h,source_sha256=sk,full_cross_receipt_sha256=kk,provenance=pr)
    if setup_only:
        return report
    prev=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    checks,Y=drive(d,.05,2048,x_start=.04,initial=query,stepper=step,restart=prev,
                   checkpoint=lambda x:write_drive_checkpoint(checkpoint,x),checkpoint_stride=128,provenance=pr)
    x=json.loads(checkpoint.read_bytes())
    if x['completed']!=2048 or x['checks']!=checks or len(checks)!=8 or x['payload_sha256']!=drive_checkpoint_digest(x) or max(abs(z['continuity_residual']) for z in checks)>3e-10:
        raise ArithmeticError('Only complete conservative crossed count evolution accepted')
    c=d.ds.w*Y[:d.ns]
    return dict(**report,endpoint_cold_count=fsum(c[mask]),endpoint_whole_scalar_count=fsum(c),
                endpoint_pole=d.state(.05,Y)[0],completed_payload_sha256=x['payload_sha256'],checks=checks,
                scope='Both donors at .04 but only cold scalar count swapped; not a physical production or convergence statement.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('history',type=Path)
    p.add_argument('full_cross',type=Path)
    p.add_argument('charge',choices=(2,4),type=int)
    p.add_argument('checkpoint',type=Path)
    p.add_argument('--setup-only',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
