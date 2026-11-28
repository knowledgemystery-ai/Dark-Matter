"""Partition the shifted-clock discrepancy into cold and warm cells."""
import argparse
import hashlib
import json
from pathlib import Path
from math import fsum

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(card,shifted,ultra,shifted_receipt,ultra_receipt,control=None):
    (F,hf),(A,ha),(B,hb)=map(read,(card,shifted,ultra))
    raw=[p.read_bytes() for p in (shifted_receipt,ultra_receipt)]
    S,U=map(json.loads,raw)
    tag,d,rep=chart_for(F,budget_control=control)
    if (S['input_sha256'][1]!=hf or U['input_hashes'][1]!=hf or U['x_end']!=.04 or S['x_end']!=.04
        or [A['completed'],B['completed']]!=[1024,4096]
        or [x['payload_sha256'] for x in (A,B)]!=[S['checkpoint_payload_sha256'],U['checkpoint_payload_sha256']]
        or any(x['payload_sha256']!=drive_checkpoint_digest(x) for x in (A,B))):
        raise ArithmeticError('Exact finished shifted and uniform-state receipts required')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    p,q=[gauss(dict(F,steps=n,nodal_occupation=D['occupation']),base,U['mass'],.04,
               expected_steps=n,budget_control=control) for D,n in ((A,1024),(B,4096))]
    v=(p[2]-q[2])/q[0][:2]
    z=np.array([fsum(a) for a in v.T])
    if not np.array_equal(q[0],U['terminal_moments'][-1]) or np.max(abs(z-S['relative_to_uniform4096'][-1][:2]))>2e-14:
        raise ArithmeticError('Terminal measure does not sum to exact clock receipt')
    R={}
    L=np.array([fsum(abs(a)) for a in v.T])
    for name,g in (('cold',rep==1),('warm',rep>1)):
        R[name]=dict(signed=[fsum(a) for a in v[g].T],absolute_share=[fsum(abs(a))/b for a,b in zip(v[g].T,L,strict=True)])
    ids=np.argmax(abs(v),axis=0)
    return dict(chart=tag,x=.04,fields=['number','pressure'],signed=z.tolist(),absolute=L.tolist(),groups=R,
                largest_cells=[dict(q=d.ds.ed[i:i+2].tolist(),signed=v[i,j]) for j,i in enumerate(ids)],
                input_sha256=[hf,ha,hb,*[hashlib.sha256(a).hexdigest() for a in raw]],
                scope='Actual shifted 1024 against fixed uniform 4096 endpoint, not asymptotic numerical error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('card','shifted','ultra','shifted_receipt','ultra_receipt'):
        p.add_argument(k,type=Path)
    p.add_argument('--control',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
