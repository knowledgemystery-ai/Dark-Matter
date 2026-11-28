"""Continue one exact frozen scalar/charged crossed restart on an unchanged terminal clock."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_startup_pilot import drive, drive_checkpoint_digest, write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step as fast
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS, load
from early_pole_resolved_parent_reference import S

D='2026-11-28'


def prepare(root,history,charge,scalar):
    a=history.read_bytes()
    H=json.loads(a)
    row=next(r for r in H['rows'] if r['levels']==[2,4])
    if charge==scalar or set((charge,scalar))!={2,4} or row['epochs']!=[.01,.02,.04,.05]:
        raise ArithmeticError('Only the two declared non-native terminal controls are eligible')
    d=CellChart(load('fine')[0],SPECS['fine'],charge,.005)
    raw=[]; ids=[]
    for p in (charge,scalar):
        for x in (4,5):
            v,h=read(root/f'early-pole-cell-endpoint-fine-d{p}-x0{x}-state-{D}.json.gz')
            if (v['grid'],v['subcells_per_archived_cell'],v['x_end'],v['steps'],v['reference_x'])!=('fine',p,x/100,2048,.005):
                raise ArithmeticError('Exact fine D2 and D4 states from matched clocks required')
            if h!=H['time_states_sha256'][{2:1,4:2}[p]][x-2]:
                raise ArithmeticError('The state must be the exact accepted time-budget input')
            raw.append(v)
            ids.append(h)
    d_other=CellChart(load('fine')[0],SPECS['fine'],scalar,.005)
    if not all(np.array_equal(getattr(d.ds,k),getattr(d_other.ds,k)) for k in ('ed','p','w')):
        raise ArithmeticError('Crossing scalar state must not change the scalar grid')
    y=np.asarray(raw[0]['endpoint_occupation'])
    other=np.asarray(raw[2]['endpoint_occupation'])
    v=np.r_[other[:d.ns],y[d.ns:]]
    if y.shape!=(d.n,) or other.shape!=(d_other.n,) or len(v)!=d.n or not np.all(np.isfinite(v)) or min(v)<0. or max(v[d.ns:])>1.:
        raise ArithmeticError('Crossed state must preserve exact Pauli and positivity gates')
    src=[Path(__file__),Path(drive.__code__.co_filename),Path(fast.__code__.co_filename)]
    src_ids=[hashlib.sha256(f.read_bytes()).hexdigest() for f in src]
    key=hashlib.sha256(json.dumps([*ids,hashlib.sha256(a).hexdigest(),*src_ids,int(charge),int(scalar)],
                                 separators=(',',':')).encode()).hexdigest()
    mask=d.ds.ed[1:]<=1e-5
    fi=np.asarray(raw[1]['endpoint_occupation'])
    xa=d.state(.04,y)
    xb=d.state(.04,v)
    det=dict(charge_donor=charge,scalar_donor=scalar,start=.04,end=.05,steps=2048,
             native_start_scalar_yield=fsum(d.ds.w*y[:d.ns])/S,
             crossed_start_scalar_yield=fsum(d.ds.w*v[:d.ns])/S,
             native_endpoint_scalar_yield=fsum(d.ds.w*fi[:d.ns])/S,
             fixed_native_start_energy_clock=[xa[0],xa[2],xa[3]],crossed_start_energy_clock=[xb[0],xb[2],xb[3]],
             canonical_input_states_sha256=ids,source_sha256=src_ids,
             terminal_history_budget_sha256=hashlib.sha256(a).hexdigest(),provenance=key)
    return d,v,mask,det


def run(root,history,charge,scalar,checkpoint,setup_only=False):
    d,v,mask,data=prepare(root,history,charge,scalar)
    if setup_only:
        return data
    old=json.loads(checkpoint.read_bytes()) if checkpoint.exists() else None
    checks,z=drive(d,.05,2048,x_start=.04,initial=v,stepper=fast,restart=old,
                   checkpoint=lambda x:write_drive_checkpoint(checkpoint,x),checkpoint_stride=128,
                   provenance=data['provenance'])
    p=json.loads(checkpoint.read_bytes())
    if (p['completed']!=2048 or p['payload_sha256']!=drive_checkpoint_digest(p)
        or abs(p['x_reached']-.05)>2e-17 or p['checks']!=checks or len(checks)!=8
        or max(abs(t['continuity_residual']) for t in checks)>3e-10):
        raise ArithmeticError('Exact completed continuation and checked conservation required')
    r=d.state(.05,z)
    return dict(**data,endpoint_cold_count=fsum((d.ds.w*z[:d.ns])[mask]),
                endpoint_scalar_count=fsum(d.ds.w*z[:d.ns]),endpoint_scalar_pole=r[0],
                complete_checkpoint_payload_sha256=p['payload_sha256'],checks=checks,
                scope='One actual crossed scalar/charged restart at x=.04 for this numerical model, not a unique decomposition of the whole earlier history.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('history',type=Path)
    p.add_argument('charge',type=int,choices=(2,4))
    p.add_argument('scalar',type=int,choices=(2,4))
    p.add_argument('checkpoint',type=Path)
    p.add_argument('--setup-only',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
