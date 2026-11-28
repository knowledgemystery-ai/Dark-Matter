"""Same-start one-chart third-clock check of the actual mixed300 count and pressure."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive,write_drive_checkpoint,drive_checkpoint_digest
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(previous,fine,accepted,comparison,control,checkpoint):
    (a,ha),(b,hb)=[read(p) for p in (previous,fine)]
    t_raw,c_raw=[p.read_bytes() for p in (accepted,comparison)]
    t,c=json.loads(t_raw)['rows'][-1],json.loads(c_raw)
    tag,d,_=chart_for(b,budget_control=control)
    source=hashlib.sha256(control.read_bytes()).hexdigest()
    old=c['rows'][1]
    h=hashlib.sha256(t_raw).hexdigest()
    if (tag!='mixed300' or a['x_end']!=.02 or b['x_end']!=.04 or b['steps']!=2048 or t['x']!=.04
            or b['actual_previous_raw_sha256']!=ha or hb!=t['states_sha256'][2]
            or c['accepted_sha256']!=h or c['steps']!=[1024,2048] or old['chart']!='mixed300'
            or [old['sha256'][0],old['sha256'][2]]!=[ha,hb]):
        raise ArithmeticError('Need exact matched mixed300 x=.02 restart and published x=.04 comparison')
    paths=[Path(__file__),Path(step.__code__.co_filename),Path(drive.__code__.co_filename)]
    hashes=[hashlib.sha256(p.read_bytes()).hexdigest() for p in paths]
    lock=hashlib.sha256(''.join([ha,hb,h,hashlib.sha256(c_raw).hexdigest(),source,*hashes]).encode()).hexdigest()
    n=4096
    previous_checkpoint=json.loads(checkpoint.read_bytes()) if checkpoint.is_file() else None
    checks,y=drive(d,.04,n,x_start=.02,initial=np.asarray(a['nodal_occupation']),stepper=step,
                   checkpoint=lambda pkt:write_drive_checkpoint(checkpoint,pkt),restart=previous_checkpoint,
                   checkpoint_stride=128,provenance=lock)
    q=json.loads(checkpoint.read_bytes())
    if (q['payload_sha256']!=drive_checkpoint_digest(q) or q['completed']!=n or q['checks']!=checks
            or len(checks)!=8 or max(abs(r['continuity_residual']) for r in checks)>3e-10
            or max(r['pair_number_error'] for r in checks)>3e-13):
        raise ArithmeticError('Incomplete or ineligible third-clock record')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    mass=t['common_native_mass']
    fine_again=gauss(b,base,mass,.04,expected_steps=2048,budget_control=control)[0]
    if not np.array_equal(fine_again,old['fine']):
        raise ArithmeticError('Old archived moments did not reproduce in the third-clock measurement')
    hi=gauss(dict(b,steps=n,nodal_occupation=y.tolist()),base,mass,.04,expected_steps=n,budget_control=control)[0]
    c0,c1,c2=np.array([old['coarse'],fine_again,hi])
    delta=np.array([c0/c1-1.,c1/c2-1.])
    return dict(steps=[1024,2048,4096],x_start=.02,x_end=.04,fields=c['fields'],mass=mass,
                terminal_moments=[c0.tolist(),c1.tolist(),c2.tolist()],previous_over_next_minus_one=delta.tolist(),
                input_hashes=[ha,hb,h,hashlib.sha256(c_raw).hexdigest()],driver_source_hashes=hashes,
                control_sha256=source,checkpoint_payload_sha256=q['payload_sha256'],provenance=lock,checks=checks,
                scope='One chart with fixed exact restart only. No full-history or cross-chart time-error claim.')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for k in ('previous','fine','accepted','comparison','control','checkpoint'):
        cli.add_argument(k,type=Path)
    print(json.dumps(run(**vars(cli.parse_args())),indent=2,allow_nan=False))
