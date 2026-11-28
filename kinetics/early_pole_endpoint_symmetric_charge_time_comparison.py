"""Compare the actual 4096 histories against the exact 2048 finite symmetric run."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_charge_symmetric_null_terminal import setup
from early_pole_endpoint_symmetric_charge_history_budget import classify


def run(root,previous_summary,previous_history,records,checkpoints):
    a,b=[p.read_bytes() for p in (previous_summary,previous_history)]
    ha,hb=[hashlib.sha256(v).hexdigest() for v in (a,b)]
    A,B=map(json.loads,(a,b))
    if B['input_summary_sha256']!=ha:
        raise ArithmeticError('Need matched 2048 terminal and trajectory audits')
    d,y,z,S=setup(root)
    mass=[]
    checks=[]
    refs=[]
    for sign,g,p in zip((-1,0,1),records,checkpoints,strict=True):
        w=p.read_bytes()
        raw=g.read_bytes()
        Z,T=map(json.loads,(raw,w))
        if (Z['steps']!=4096 or T['completed']!=4096 or Z['sign']!=sign or Z['setup']!=S
            or Z['source_hashes'][0]!=ha or Z['complete_payload_sha256']!=T['payload_sha256']
            or T['payload_sha256']!=drive_checkpoint_digest(T) or Z['checks']!=T['checks']
            or len(T['checks'])!=8 or abs(T['x_reached']-.05)>2e-17
            or max(abs(u['continuity_residual']) for u in T['checks'])>3e-10):
            raise ArithmeticError('Only exact complete matched 4096 histories permitted')
        Y=np.asarray(T['occupation'])
        if Y.shape!=(d.n,) or not np.all(np.isfinite(Y)) or min(Y)<0 or max(Y[d.ns:])>1:
            raise ArithmeticError('Fine checkpoint inadmissible')
        g=d.ds.w*Y[:d.ns]
        v=fsum(g[d.ds.ed[1:]<=1e-5])
        if v!=Z['final_cold'] or fsum(g)!=Z['final_total'] or d.state(.05,Y)[0]!=Z['final_pole']:
            raise ArithmeticError('Exact final observable did not replay')
        mass.append(v)
        checks.append(Z['checks'])
        refs.append(dict(sign=sign,receipt_sha256=hashlib.sha256(raw).hexdigest(),payload_sha256=T['payload_sha256']))
    V=classify(*checks)
    pairs=[]
    for lo,hi in zip(B['rows'],V['rows'],strict=True):
        if abs(lo['x']-hi['x'])>1.5e-17:
            raise ArithmeticError('Compare only the same saved times')
        a=np.array([lo['symmetric_odd_yield'],lo['symmetric_even_yield']])
        b=np.array([hi['symmetric_odd_yield'],hi['symmetric_even_yield']])
        pairs.append(dict(x=hi['x'],n2048_odd_even=a.tolist(),n4096_odd_even=b.tolist(),n4096_minus_n2048=(b-a).tolist(),
                          preserved_sample_signs=(np.sign(a)==np.sign(b)).tolist()))
    old=A['results']['cold_count']
    v=np.array([old['endpoint_low'],old['nominal'],old['endpoint_high']])
    g=np.array(mass)
    metrics=lambda q: np.array([(q[2]-q[0])/2.,(q[2]+q[0])/2.-q[1]])
    return dict(rows=pairs,old_terminal_cold=v.tolist(),fine_terminal_cold=g.tolist(),
                cold_odd_even_at_2048=metrics(v).tolist(),cold_odd_even_at_4096=metrics(g).tolist(),
                fine_native_history_brackets=V['brackets_with_opposite_saved_odd_signs'],
                baseline_brackets=B['brackets_with_opposite_saved_odd_signs'],
                mesh_counts=[2048,4096],exact_2048_sha256=[ha,hb],fine_receipts=refs,
                scope='Actual same-start time refinement of one symmetric charge intervention; not a physical sign, spatial convergence, or higher-order temporal estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','previous_summary','previous_history'):
        p.add_argument(k,type=Path)
    for k in ('records','checkpoints'):
        p.add_argument('--'+k,type=Path,nargs=3,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
