"""Decompose exact uniform1024 and shifted1024 discrepancies relative to uniform4096."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def run(card,coarse,shifted,ultra,paired,shifted_receipt,ultra_receipt,control=None):
    (F,hf),(A,ha),(B,hb),(C,hc)=map(read,(card,coarse,shifted,ultra))
    raw=[p.read_bytes() for p in (paired,shifted_receipt,ultra_receipt)]
    P,S,U=map(json.loads,raw)
    h=[hashlib.sha256(u).hexdigest() for u in raw]
    tag,d,rep=chart_for(F,budget_control=control)
    name='G9' if tag==9 else tag
    row=next(u for u in P['rows'] if u['chart']==name)
    if (P['x']!=.05 or F['x_end']!=.05 or U['chart']!=name or S['chart']!=name or U['source_sha256'][1:4]!=[hf,*h[:2]]
        or row['sha256'][3]!=ha or [x['completed'] for x in (A,B,C)]!=[1024,1024,4096]
        or B['payload_sha256']!=S['checkpoint_payload_sha256'] or C['payload_sha256']!=U['checkpoint_payload_sha256']
        or any(x['payload_sha256']!=drive_checkpoint_digest(x) for x in (A,B,C))):
        raise ArithmeticError('Exact completed reference and 1024 checkpoints required')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    def eval(z,n):
        return gauss(dict(F,steps=n,nodal_occupation=z['occupation']),base,P['mass'],.05,
                     expected_steps=n,budget_control=control)
    va,vb,vc=[eval(z,k) for z,k in ((A,1024),(B,1024),(C,4096))]
    if not np.array_equal([va[0],vb[0],vc[0]],[row['coarse'],S['moments'][1],U['terminal_moments']]):
        raise ArithmeticError('Each independently measured accepted receipt must replay')
    rows=[]
    for i,(z,k) in enumerate(((va,'uniform1024'),(vb,'shifted1024'))):
        q=(z[2]-vc[2])/vc[0][:2]
        L=[fsum(abs(t)) for t in q.T]
        total=[fsum(t) for t in q.T]
        if max(abs(np.array(total)-U['relative_to_uniform4096'][i][:2]))>2e-14:
            raise ArithmeticError('Cell totals do not reproduce endpoint comparison')
        groups={}
        for label,mask in (('cold',rep==1),('warm',rep>1)):
            groups[label]=dict(signed=[fsum(t) for t in q[mask].T],
                               share_of_absolute=[fsum(abs(t))/r for t,r in zip(q[mask].T,L,strict=True)])
        top=np.argsort(-abs(q),axis=0)[:4]
        rows.append(dict(clock=k,signed=total,absolute_l1=L,groups=groups,
                         largest_cells_by_field=[[dict(q=d.ds.ed[j:j+2].tolist(),signed=q[j,i])
                                                for j in top[:,i]] for i in (0,1)]))
    return dict(chart=name,x=.05,fields=['count','pressure'],rows=rows,input_sha256=[hf,ha,hb,hc,*h],
                scope='Exact per-cell endpoint discrepancy bookkeeping versus actual uniform4096; not an independent global error or continuum field estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('card','coarse','shifted','ultra','paired','shifted_receipt','ultra_receipt'):
        p.add_argument(k,type=Path)
    p.add_argument('--control',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
