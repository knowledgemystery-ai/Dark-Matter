"""Bind and measure actual uniform4096 interior states at frozen observable mass."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_warm_embedded_step_probe import point


def run(previous,reference,final,uniform,paired,control,samples):
    (p,ph),(r,rh),(F,fh)=map(read,(previous,reference,final))
    raw=[z.read_bytes() for z in (uniform,paired)]
    u,v=map(json.loads,raw)
    h=[hashlib.sha256(s).hexdigest() for s in raw]
    required=(1536,2048,2560,3072,3584)
    tag=u['chart']
    if ((p['x_end'],r['x_end'],u['x_start'],u['x_end'],v['x'])!=(.04,.05,.04,.05,.05)
        or u['steps']!=4096 or F['completed']!=4096 or F['target_steps']!=4096 or F['payload_sha256']!=drive_checkpoint_digest(F)
        or u['checkpoint_payload_sha256']!=F['payload_sha256'] or u['source_sha256'][:3]!=[ph,rh,h[1]]
        or r['actual_previous_raw_sha256']!=ph or tag not in ('G9','mixed300')):
        raise ArithmeticError('Require exact completed independent uniform4096 trajectory')
    c=hashlib.sha256(control.read_bytes()).hexdigest()
    if tag=='mixed300' and u['source_sha256'][-1]!=c:
        raise ArithmeticError('Mixed chart control did not replay')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    t=np.linspace(0.,np.log(1.25),4097)
    vs=[read(f) for f in samples]
    if any(z['completed']!=n for (z,_),n in zip(vs,required,strict=True)):
        raise ArithmeticError('Use predetermined common observation points only')
    for (q,_),n in zip(vs,required,strict=True):
        if (q['target_steps']!=4096 or q['signature']!=F['signature'] or q['payload_sha256']!=drive_checkpoint_digest(q)
            or q['trace']!=F['trace'][:n] or q['caps']!=F['caps'][:n] or q['checks']!=F['checks'][:n//512]
            or q['x_reached']!=float(.04*np.exp(t[n]))):
            raise ArithmeticError('Archived progress card is not an exact prefix of the final control')
    states=[(0,p['x_end'],p['nodal_occupation']),*[(q['completed'],q['x_reached'],q['occupation']) for q,_ in vs],
            (4096,r['x_end'],F['occupation'])]
    out=[]
    for k,x,y in states:
        A=dict(r,x_end=x,steps=4096,nodal_occupation=y)
        row=point(A,tag,control if tag=='mixed300' else None,v['mass'],base)
        out.append(dict(uniform4096_prefix_steps=k,**row))
    return dict(chart=tag,common_native_x05_mass=v['mass'],fields=v['fields'],rows=out,
                input_sha256=[ph,rh,fh,*h,c,*[j for _,j in vs]],
                scope='Actual interior states, local embedded one/two/four probes at one fixed native x05 pressure mass. Not a propagated global error budget.')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    for k in ('previous','reference','final','uniform','paired','control'):
        cli.add_argument(k,type=Path)
    cli.add_argument('samples',type=Path,nargs=5)
    print(json.dumps(run(**vars(cli.parse_args())),indent=2,allow_nan=False))
