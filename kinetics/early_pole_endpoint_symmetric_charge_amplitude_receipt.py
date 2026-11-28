"""Audit exact equal-opposite half amplitudes against the prior same-clock controls."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_charge_symmetric_null_terminal import setup
from early_pole_endpoint_symmetric_charge_half_amplitude import run as setup_run
from early_pole_endpoint_symmetric_charge_history_budget import classify


def compare(half,full):
    x,y=np.asarray(half),np.asarray(full)
    if x.shape!=(2,) or y.shape!=(2,) or not np.all(np.isfinite(np.r_[x,y])):
        raise ArithmeticError('Two finite odd/even quantities required')
    a=np.array([.5,.25])*y
    return dict(half_odd_even=x.tolist(),full_odd_even=y.tolist(),linear_quadratic_target=a.tolist(),
                half_minus_linear_quadratic_target=(x-a).tolist(),
                ratio_to_linear_quadratic_target=[float(z/t) if abs(t)>1e-18 else None for z,t in zip(x,a,strict=True)])


def run(root,prior,nominal,records,checkpoints):
    raw=prior.read_bytes()
    P=json.loads(raw)
    h=hashlib.sha256(raw).hexdigest()
    H=nominal.read_bytes()
    N=json.loads(H)
    d,_,_,C=setup(root)
    if (P['mesh_counts']!=[2048,4096] or N['setup']!=C or N['steps']!=4096 or N['sign']!=0
        or hashlib.sha256(H).hexdigest()!=P['fine_receipts'][1]['receipt_sha256']
        or len(N['checks'])!=8 or N['final_cold']!=P['fine_terminal_cold'][1]):
        raise ArithmeticError('Need the corresponding exact archived 4096 nominal control')
    result=[]; states=[]; ids=[]
    for sign,record,cp in zip((-1,1),records,checkpoints,strict=True):
        r,b=record.read_bytes(),cp.read_bytes()
        D,B=json.loads(r),json.loads(b)
        X=setup_run(root,prior,sign,cp,setup_only=True)
        if (any(D[k]!=v for k,v in X.items()) or B['completed']!=4096 or B['checks']!=D['checks']
            or abs(B['x_reached']-.05)>2e-17 or B['payload_sha256']!=drive_checkpoint_digest(B)
            or B['payload_sha256']!=D['payload_sha256'] or len(B['checks'])!=8
            or max(abs(s['continuity_residual']) for s in B['checks'])>3e-10):
            raise ArithmeticError('Exact full 4096 half-amplitude checkpoint required')
        o=np.asarray(B['occupation'])
        if o.shape!=(d.n,) or not np.all(np.isfinite(o)) or min(o)<0. or max(o[d.ns:])>1.:
            raise ArithmeticError('Final finite half-amplitude state insufficient')
        f=d.ds.w*o[:d.ns]
        a=fsum(f[d.ds.ed[1:]<=1e-5])
        if (a!=D['final_cold'] or fsum(f)!=D['final_total'] or d.state(.05,o)[0]!=D['final_pole']):
            raise ArithmeticError('Recovered exact final scalar observables must agree')
        result.append(D['checks'])
        states.append(a)
        ids.append(dict(sign=sign,receipt_sha256=hashlib.sha256(r).hexdigest(),payload_sha256=B['payload_sha256']))
    G=classify(result[0],N['checks'],result[1])
    rows=[]
    for a,b in zip(G['rows'],P['rows'],strict=True):
        if a['x']!=b['x']:
            raise ArithmeticError('Shared saved clock mismatch')
        rows.append(dict(x=a['x'],**compare([a['symmetric_odd_yield'],a['symmetric_even_yield']],b['n4096_odd_even'])))
    u,v,w=states[0],N['final_cold'],states[1]
    return dict(rows=rows,half_cold_terminal=[u,v,w],terminal_cold=compare([(w-u)/2.,(w+u)/2.-v],P['cold_odd_even_at_4096']),
                half_odd_brackets=G['brackets_with_opposite_saved_odd_signs'],full_odd_brackets=P['fine_native_history_brackets'],
                prior_time_refinement_sha256=h,nominal_receipt_sha256=hashlib.sha256(H).hexdigest(),half_receipts=ids,
                scope='One finite factor-of-two amplitude change and exact 4096 time clock. No asymptotic, physical stability or convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','prior','nominal'):
        p.add_argument(k,type=Path)
    for k in ('records','checkpoints'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
