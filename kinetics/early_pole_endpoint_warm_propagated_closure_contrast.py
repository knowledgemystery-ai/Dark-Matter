"""Compare reconstructed and native nodal pressure along the three exact G9 pulse pairs."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_gauss_warm_comparison import chart_for, gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart, SPECS, load
from early_pole_endpoint_warm_embedded_step_probe import advance
from early_pole_endpoint_warm_pressure_cells import partition

DATE='2026-11-28'
W='early-pole-warm-gauss9-'
P='early-pole-fine-charge2-scalar4-gauss9warm-'


def run(root, observer=None):
    root=Path(root)
    card, acol=read(root/f'{P}x04-state-{DATE}.json.gz')
    endcard, bcol=read(root/f'{P}x05-state-{DATE}.json.gz')
    h={'start_card':acol,'end_card':bcol}
    order, chart, repeat=chart_for(card)
    if order!=9 or chart_for(endcard)[0]!=9 or card['x_end']!=.04 or endcard['x_end']!=.05:
        raise ArithmeticError('Exact G9 reference geometry required')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    rows=[]
    for label,prefix,count in (('',None,512),('late-',2560,192),('near-end-',3584,64)):
        r=(root/f'{W}{label}propagated-step-defect-{DATE}.json').read_bytes()
        h[label+'receipt']=hashlib.sha256(r).hexdigest()
        D=json.loads(r)
        if D['chart']!='G9' or D['x_end']!=.05:
            raise ArithmeticError('Published G9 propagation receipt required')
        m=D['common_native_x05_mass']
        x=D['x_start']
        if prefix is None:
            width=D['local_log_interval']
            starts=D['before_after']
            old=[q['fine_minus_coarse'] for q in starts]
            checks=D['records']
            if acol!=D['input_sha256'][0]:
                raise ArithmeticError('Wrong actual first-pulse state')
            y=card['nodal_occupation']
        else:
            G, h[label+'prefix']=read(root/f'{W}x05-uniform4096-progress-{prefix}-{DATE}.json.gz')
            if (D['input_sha256'][:2]!=[bcol,h[label+'prefix']]
                or G['payload_sha256']!=drive_checkpoint_digest(G) or G['completed']!=prefix or G['x_reached']!=x):
                raise ArithmeticError('Wrong held interior state')
            y=G['occupation']
            width=D['log_pulse_width']
            starts=D['late_before_after']
            old=[q['relative_fine_minus_coarse'] for q in starts]
            checks=D['checks']
        local=[advance(chart,np.asarray(y),x,width,k)[0] for k in (1,2)]
        terminal=[]
        for arm,c in zip(('one','two'),checks,strict=True):
            B,h[label+arm]=read(root/f'{W}{label}propagated-step-{arm}-checkpoint-{DATE}.json.gz')
            if (B['payload_sha256']!=c['checkpoint_payload_sha256']
                or B['payload_sha256']!=drive_checkpoint_digest(B) or B['completed']!=count or abs(B['x_reached']-.05)>np.spacing(.05)):
                raise ArithmeticError('Wrong completed propagation branch')
            terminal.append(B['occupation'])
        q=chart.ds.p
        multiplier=q*q/np.hypot(q,m)
        data=[]
        for stamp,occ,reference in zip(starts,(local,terminal),old,strict=True):
            x0=stamp['x']
            (u,_,a),(v,_,b)=[gauss(dict(card,x_end=x0,nodal_occupation=np.asarray(Y).tolist()),base,m,x0,
                                            expected_steps=2048) for Y in occ]
            reconstructed=(v-u)/u
            if not np.allclose(reconstructed,reference,rtol=0.,atol=5e-15):
                raise ArithmeticError('Reconstructed pressure must reproduce the published pair')
            N=[chart.ds.w*np.asarray(Y)[:chart.ns] for Y in occ]
            nodal_total=fsum(N[0]*multiplier)
            nodal_relative=fsum((N[1]-N[0])*multiplier)/nodal_total
            if observer is not None:
                observer(x,x0,q,N[1]-N[0],nodal_total)
            M=partition(a,b)/u[1]
            groups={name:[fsum(z) for z in M[mask].T] for name,mask in
                    (('cold',repeat==1),('warm',repeat>1))}
            if abs(fsum(M[:,0])-reconstructed[1])>3e-15 or abs(fsum(M[:,1])-reconstructed[1]+fsum(M[:,2]))>3e-15:
                raise ArithmeticError('Pressure decomposition must close')
            data.append(dict(x=x0,reconstructed_relative_pressure=reconstructed[1],nodal_relative_pressure=nodal_relative,
                             total_count_relative=reconstructed[0],relative_components_by_cell=groups))
        rows.append(dict(x_insertion=x,common_pressure_mass=m,pulse_width=width,after_pulse_and_endpoint=data))
    return dict(rows=rows,component_fields=['total_pressure','cell_count','mean_pressure_within_cell'],
                input_uncompressed_sha256=h,scope='Three exact paired runs, fixed mass on each; nodal pressure is an alternate diagnostic, not a physical pressure estimate or a time-convergence certificate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
