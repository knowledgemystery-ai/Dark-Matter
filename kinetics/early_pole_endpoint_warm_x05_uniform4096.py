"""Independently extend the shorter x05 interval to a uniform 4096-step control."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive,drive_checkpoint_digest,write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(previous,fine,uniform,shifted,checkpoint,control=None):
    (a,ha),(b,hb)=map(read,(previous,fine))
    raw=[p.read_bytes() for p in (uniform,shifted)]
    U,S=map(json.loads,raw)
    ids=[hashlib.sha256(r).hexdigest() for r in raw]
    tag,d,_=chart_for(b,budget_control=control)
    name='G9' if tag==9 else tag
    c=hashlib.sha256(control.read_bytes()).hexdigest() if control is not None else None
    if (name not in ('G9','mixed300') or (a['x_end'],b['x_end'])!=(.04,.05) or b['actual_previous_raw_sha256']!=ha
        or b['steps']!=2048 or U['steps']!=[1024,2048] or U['x']!=.05 or S['inputs']!=[ha,hb,ids[0]]+([c] if c else [])
        or S['chart']!=name or (S['x_start'],S['x_end'],S['clock_exponent'])!=(.04,.05,1.2)):
        raise ArithmeticError('Exact fixed x04 restart and validated x05 contrasting clocks required')
    q=next(i for i in U['rows'] if i['chart']==name)
    if [ha,hb]!=[q['sha256'][0],q['sha256'][2]] or [S['moments'][0],S['moments'][-1]]!=[q['coarse'],q['fine']]:
        raise ArithmeticError('Previous comparison must replay exact measured baseline')
    keys=[hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),Path(drive.__code__.co_filename),Path(step.__code__.co_filename))]
    lock=hashlib.sha256(''.join([ha,hb,*ids,*keys,*([c] if c else [])]).encode()).hexdigest()
    n=4096
    old=json.loads(checkpoint.read_bytes()) if checkpoint.is_file() else None
    checks,y=drive(d,.05,n,x_start=.04,initial=np.asarray(a['nodal_occupation']),stepper=step,
                   checkpoint=lambda v:write_drive_checkpoint(checkpoint,v),restart=old,checkpoint_stride=128,provenance=lock)
    z=json.loads(checkpoint.read_bytes())
    if z['completed']!=n or z['payload_sha256']!=drive_checkpoint_digest(z) or z['checks']!=checks or len(checks)!=8 or any(abs(r['continuity_residual'])>3e-10 for r in checks):
        raise ArithmeticError('Final uniform4096 checkpoint must pass all prior gates')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    def measure(card,k):
        return gauss(card,base,U['mass'],.05,expected_steps=k,budget_control=control)[0]
    if not np.array_equal(measure(b,2048),q['fine']):
        raise ArithmeticError('Baseline no longer reproduces')
    f=measure(dict(b,steps=n,nodal_occupation=y.tolist()),n)
    R=np.array(S['moments'])/f-1.
    return dict(chart=name,x_start=.04,x_end=.05,fields=U['fields'],steps=n,terminal_moments=f.tolist(),
                rows=['uniform1024','shifted1024','uniform2048'],relative_to_uniform4096=R.tolist(),
                shifted_over_unshifted_1024_discrepancy=(abs(R[1])/abs(R[0])).tolist(),
                source_sha256=[ha,hb,*ids,*([c] if c else [])],signature_source_sha256=keys,
                checks=checks,checkpoint_payload_sha256=z['payload_sha256'],
                scope='Same fixed x04 restart and native mass. An independent short-window comparison against uniform4096, not a continuum or physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('previous','fine','uniform','shifted','checkpoint'):
        p.add_argument(k,type=Path)
    p.add_argument('--control',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
