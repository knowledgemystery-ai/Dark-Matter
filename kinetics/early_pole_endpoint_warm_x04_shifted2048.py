"""Run the locked shifted log-time schedule at 2048 steps for either warm chart."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_warm_g9_prescribed_clock import clock,integrate


def run(previous,fine,uniform,shifted,checkpoint,control=None):
    (A,ha),(B,hb)=map(read,(previous,fine))
    raw=[p.read_bytes() for p in (uniform,shifted)]
    U,S=map(json.loads,raw)
    h=[hashlib.sha256(r).hexdigest() for r in raw]
    tag,d,_=chart_for(B,budget_control=control)
    c=hashlib.sha256(control.read_bytes()).hexdigest() if control is not None else None
    if (tag not in (9,'mixed300') or (A['x_end'],B['x_end'])!=(.02,.04) or B['actual_previous_raw_sha256']!=ha
        or B['steps']!=2048 or U['input_hashes'][:2]!=[ha,hb] or S['input_sha256']!=[ha,hb,h[0]]+([c] if c else [])
        or (U.get('control_sha256')!=c) or U['steps']!=[1024,2048,4096]
        or S['steps']!=1024 or S['exponent']!=1.2 or (U['x_start'],U['x_end'])!=(.02,.04)):
        raise ArithmeticError('Require exact chart start, accepted 4096 result and shift prescription')
    files=(Path(__file__),Path(clock.__code__.co_filename),Path(step.__code__.co_filename))
    keys=[hashlib.sha256(p.read_bytes()).hexdigest() for p in files]
    lock=hashlib.sha256(''.join([ha,hb,*h,*keys,*([c] if c else [])]).encode()).hexdigest()
    n=2048
    y,checks,sig=integrate(d,np.asarray(A['nodal_occupation']),clock(n),checkpoint,lock)
    q=json.loads(checkpoint.read_bytes())
    if (q['payload_sha256']!=drive_checkpoint_digest(q) or q['signature']!=sig or q['completed']!=n
        or checks!=q['checks'] or len(checks)!=8):
        raise ArithmeticError('Missing exact final shifted2048 checkpoint')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    def measure(card):
        return gauss(card,base,U['mass'],.04,expected_steps=n,budget_control=control)[0]
    if not np.array_equal(measure(B),U['terminal_moments'][1]):
        raise ArithmeticError('Uniform state and measurement reference differ')
    k=measure(dict(B,nodal_occupation=y.tolist()))
    r=np.asarray(U['terminal_moments'])
    out=np.array([r[0]/r[2]-1.,np.array(S['relative_to_uniform4096'][-1]),r[1]/r[2]-1.,k/r[2]-1.])
    return dict(chart=tag,steps=n,exponent=1.2,terminal_moments=k.tolist(),reference_moments=r[2].tolist(),
                source_sha256=[ha,hb,*h,*([c] if c else [])],signature_source_sha256=keys,
                rows=['uniform1024','shifted1024','uniform2048','shifted2048'],fields=U['fields'],
                relative_to_uniform4096=out.tolist(),shifted_over_uniform2048_discrepancy=(abs(out[-1])/abs(out[-2])).tolist(),
                checks=checks,checkpoint_payload_sha256=q['payload_sha256'],
                scope='Exact same-start half-width refinement on a prescribed shifted grid; uniform 4096 is a comparison only, not a continuum solution.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('previous','fine','uniform','shifted','checkpoint'):
        p.add_argument(k,type=Path)
    p.add_argument('--control',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
