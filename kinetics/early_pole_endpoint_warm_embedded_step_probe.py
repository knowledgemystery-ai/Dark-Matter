"""Independent local one/two/four-step number and pressure defects at held starts."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load


def advance(d,y,x,width,k):
    init=y.copy()
    t=np.linspace(0.,width,k+1)
    for a,b in zip(t[:-1],t[1:],strict=True):
        xm=x*np.exp((a+b)/2)
        M,_,_,H=d.state(xm,y)
        mesh,_=d.chart(xm,M,3)
        z,_=step(mesh,d,y,H*(b-a)/2)
        M,_,_,H=d.state(xm,z)
        mesh,_=d.chart(xm,M,3)
        y,_=step(mesh,d,y,H*(b-a),rate_state=z)
    s,f,p=[fsum(a) for a in np.split(d.w*(y-init),(d.ns,d.iP))]
    err=max(abs(s-f),abs(s+p))
    if min(y)<0 or max(y[d.ns:])>1. or err>3e-13:
        raise ArithmeticError('Local probe violates the count or occupation gate')
    return y,err


def point(A,tag,control,m,base):
    name,d,repeat=chart_for(A,budget_control=control)
    if (tag,name) not in (('G9',9),('mixed300','mixed300')):
        raise ArithmeticError('Chart identity is not the one in the lock')
    x=A['x_end']
    y=np.asarray(A['nodal_occupation'])
    width=np.log(1.25)/1024
    n=A['steps']
    def f(z):
        return gauss(dict(A,nodal_occupation=z.tolist()),base,m,x,expected_steps=n,budget_control=control)
    ref=f(y)[0]
    levels=[]
    for h in (4*width,width):
        data=[advance(d,y,x,h,k) for k in (1,2,4)]
        cards=[f(z) for z,_ in data]
        q=np.array([a for a,_,_ in cards])
        diffs=np.diff(q,axis=0)/ref
        cells=(cards[1][2]-cards[0][2])/ref[:2]
        if max(abs(np.array([fsum(t) for t in cells.T])-diffs[0,:2]))>3e-14:
            raise ArithmeticError('Local cell accounting does not close')
        group={}
        L=np.array([fsum(abs(t)) for t in cells.T])
        for key,mask in (('cold',repeat==1),('warm',repeat>1)):
            group[key]=dict(signed=[fsum(t) for t in cells[mask].T],
                            absolute_share=[fsum(abs(t))/a for t,a in zip(cells[mask].T,L,strict=True)])
        ix=np.argmax(abs(cells),axis=0)
        if np.any(~np.isfinite(diffs)):
            raise ArithmeticError('Nonfinite local pressure or count difference')
        levels.append(dict(log_interval=h,runs=[1,2,4],moments=q.tolist(),fine_minus_coarse_relative_to_start=diffs.tolist(),
                           maximum_pair_number_error=max(r for _,r in data),
                           fine_minus_coarse_cell_fields=['count','pressure'],cell_l1=L.tolist(),groups=group,
                           largest_cells=[dict(q=d.ds.ed[i:i+2].tolist(),signed=cells[i,j]) for j,i in enumerate(ix)]))
    return dict(x_start=x,chart=tag,common_start_mass=m,held_start_moments=ref.tolist(),levels=levels)


def run(root,control,receipts):
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    xids=['02','04','05']
    out=[]; hashes=[]
    d='2026-11-28'
    for k,f in zip(xids,receipts,strict=True):
        u=f.read_bytes(); p=json.loads(u)
        x=int(k)/100
        if p['x']!=x or p['steps']!=[1024,2048] or [r['chart'] for r in p['rows']]!=['G9','mixed300']:
            raise ArithmeticError('Unexpected prior checkpoint-clocks document')
        hashes.append(hashlib.sha256(u).hexdigest())
        for name,tag,ref in zip(('gauss9','mixed300'),('G9','mixed300'),p['rows'],strict=True):
            A,h=read(root/f'early-pole-fine-charge2-scalar4-{name}warm-x{k}-state-{d}.json.gz')
            if h!=ref['sha256'][2] or A['x_end']!=x or A['steps']!=2048:
                raise ArithmeticError('Exact directly evolved archived start required')
            hashes.append(h)
            ans=point(A,tag,control if name=='mixed300' else None,p['mass'],base)
            if ans['held_start_moments']!=ref['fine']:
                raise ArithmeticError('Reference observable did not replay')
            out.append(ans)
    return dict(fields=['count','pressure','warm_pressure'],rows=out,input_sha256=hashes,
                control_sha256=hashlib.sha256(control.read_bytes()).hexdigest(),
                scope='Local same-start embedded discrepancies with frozen measurement mass; neither attribution of global error nor a physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('control',type=Path)
    p.add_argument('receipts',type=Path,nargs=3)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
