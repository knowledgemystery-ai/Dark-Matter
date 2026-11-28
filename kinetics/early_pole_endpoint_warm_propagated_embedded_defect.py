"""Measure the terminal effect of just one differing local step, then identical evolution."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive,drive_checkpoint_digest,write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_endpoint_warm_embedded_step_probe import advance


def run(card,pair,local,checkpoints):
    a,ah=read(card)
    raw=[p.read_bytes() for p in (pair,local)]
    P,L=map(json.loads,raw)
    h=[hashlib.sha256(t).hexdigest() for t in raw]
    old=next(r for r in L['rows'] if r['x_start']==.04 and r['chart']=='G9')
    if a['x_end']!=.04 or a['steps']!=2048 or P['x']!=.05 or P['rows'][0]['sha256'][0]!=ah or ah not in L['input_sha256']:
        raise ArithmeticError('Require exact held start of both prior probes')
    tag,d,repeats=chart_for(a)
    if tag!=9:
        raise ArithmeticError('One actual G9 chart only')
    y=np.asarray(a['nodal_occupation'])
    h0=np.log(1.25)/1024
    xx=.04*np.exp(h0)
    edge=[advance(d,y,.04,h0,k)[0] for k in (1,2)]
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    def held(z,m,x):
        return gauss(dict(a,x_end=x,nodal_occupation=z.tolist()),base,m,x,expected_steps=2048)
    initial=[held(z,old['common_start_mass'],.04)[0] for z in edge]
    if not np.array_equal(initial,old['levels'][-1]['moments'][:2]):
        raise ArithmeticError('Exact embedded local pair is not the one in the previous audit')
    m=P['mass']
    keys=[hashlib.sha256(p.read_bytes()).hexdigest() for p in
          (Path(__file__),Path(advance.__code__.co_filename),Path(drive.__code__.co_filename),Path(step.__code__.co_filename))]
    r=[]
    records=[]
    for j,(z,C) in enumerate(zip(edge,checkpoints,strict=True)):
        token=hashlib.sha256(('|'.join([ah,*h,*keys,str(j)])).encode()).hexdigest()
        prev=json.loads(C.read_bytes()) if C.exists() else None
        checks,Y=drive(d,.05,512,x_start=float(xx),initial=z,stepper=step,restart=prev,
                       checkpoint=lambda q:write_drive_checkpoint(C,q),checkpoint_stride=128,provenance=token)
        q=json.loads(C.read_bytes())
        if (q['completed']!=512 or q['payload_sha256']!=drive_checkpoint_digest(q)
            or q['checks']!=checks or len(checks)!=8 or max(abs(z['continuity_residual']) for z in checks)>3e-10):
            raise ArithmeticError('Both propagation arms need exact complete local conservation')
        r.append(Y)
        records.append(dict(initial_substeps=j+1,first_interval=q['checks'][0]['x'],
                            checkpoint_payload_sha256=q['payload_sha256'],max_continuity_residual=max(abs(z['continuity_residual']) for z in checks),
                            max_pair_defect=max(z['pair_number_error'] for z in checks)))
    def difference(v,x):
        V,W=[held(q,m,x) for q in v]
        E=(W[2]-V[2])/V[0][:2]
        G=(W[0]-V[0])/V[0]
        if max(abs(np.array([fsum(k) for k in E.T])-G[:2]))>2e-14:
            raise ArithmeticError('Same pressure mass and global cell sum required')
        L1=[fsum(abs(z)) for z in E.T]
        out={}
        for label,mask in (('cold',repeats==1),('warm',repeats>1)):
            out[label]=dict(signed=[fsum(z) for z in E[mask].T],absolute_share=[fsum(abs(z))/n for z,n in zip(E[mask].T,L1,strict=True)])
        return dict(x=x,fine_minus_coarse=G.tolist(),groups=out,absolute_l1=L1)
    return dict(chart='G9',x_start=.04,local_log_interval=h0,propagation_x_start=float(xx),x_end=.05,
                common_native_x05_mass=m,fields=['count','pressure','warm_pressure'],before_after=[difference(edge,float(xx)),difference(r,.05)],
                records=records,input_sha256=[ah,*h],source_sha256=keys,
                scope='Controlled finite propagated first-step contrast only; not decomposition of whole-grid error or a physical endpoint.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('card','pair','local'):
        p.add_argument(k,type=Path)
    p.add_argument('checkpoints',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
