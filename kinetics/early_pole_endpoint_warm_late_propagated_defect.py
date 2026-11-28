"""Compare a late embedded pulse on an actual uniform4096 progress state."""
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


def run(card,progress,final,receipt,original,checkpoints):
    (A,ha),(G,hg),(F,hf)=map(read,(card,progress,final))
    raw=[p.read_bytes() for p in (receipt,original)]
    U,R=map(json.loads,raw)
    h=[hashlib.sha256(p).hexdigest() for p in raw]
    n=2560
    if (A['x_end']!=.05 or A['steps']!=2048 or U['chart']!='G9' or U['source_sha256'][1]!=ha
        or G['completed']!=n or F['completed']!=4096 or G['target_steps']!=4096 or F['target_steps']!=4096
        or any(z['payload_sha256']!=drive_checkpoint_digest(z) for z in (G,F)) or G['signature']!=F['signature']
        or F['payload_sha256']!=U['checkpoint_payload_sha256'] or G['trace']!=F['trace'][:n]
        or G['caps']!=F['caps'][:n] or G['checks']!=F['checks'][:n//512]):
        raise ArithmeticError('Require the actual held interior prefix bound to the completed x05 uniform4096')
    tag,d,rep=chart_for(A)
    if tag!=9 or R['chart']!='G9' or R['x_end']!=.05:
        raise ArithmeticError('The prior first-step contrast must be the same G9 geometry and reference mass')
    m=R['common_native_x05_mass']
    width=float(np.log(1.25)/1024)
    x=float(G['x_reached'])
    if x!=float(.04*np.exp(np.linspace(0.,np.log(1.25),4097)[n])) or R['local_log_interval']!=width:
        raise ArithmeticError('Require preset insertion time and equal pulse clock')
    y=np.asarray(G['occupation'])
    x1=float(x*np.exp(width))
    start=[advance(d,y,x,width,k)[0] for k in (1,2)]
    N=192
    L=np.log(.05/x1)
    relative_step=(L/N)/(np.log(.05/R['propagation_x_start'])/512)
    if not .98<=relative_step<=1.02:
        raise ArithmeticError('Continuation resolution must remain comparable to the first insertion')
    keys=[hashlib.sha256(p.read_bytes()).hexdigest() for p in
          (Path(__file__),Path(advance.__code__.co_filename),Path(step.__code__.co_filename),Path(drive.__code__.co_filename))]
    data=[]; gates=[]
    for j,(y,q) in enumerate(zip(start,checkpoints,strict=True),1):
        key=hashlib.sha256(''.join([ha,hg,hf,*h,*keys,str(j)]).encode()).hexdigest()
        prior=json.loads(q.read_bytes()) if q.exists() else None
        checks,z=drive(d,.05,N,x_start=x1,initial=y,stepper=step,restart=prior,checkpoint_stride=64,provenance=key,
                       checkpoint=lambda v:write_drive_checkpoint(q,v))
        payload=json.loads(q.read_bytes())
        e=max(abs(a['continuity_residual']) for a in checks)
        if (payload['completed']!=N or payload['payload_sha256']!=drive_checkpoint_digest(payload)
            or payload['checks']!=checks or len(checks)!=8 or e>3e-10):
            raise ArithmeticError('Exact eligible continuation only')
        gates.append(dict(trajectory=j,checkpoint_payload_sha256=payload['payload_sha256'],max_continuity=e,
                          max_pair_defect=max(a['pair_number_error'] for a in checks)))
        data.append(z)
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    def measure(y,x):
        return gauss(dict(A,x_end=x,nodal_occupation=y.tolist()),base,m,x,expected_steps=2048)
    def diff(Y,x):
        a,b=[measure(v,x) for v in Y]
        z=(b[0]-a[0])/a[0]
        C=(b[2]-a[2])/a[0][:2]
        Z=np.array([fsum(u) for u in C.T])
        if np.max(abs(Z-z[:2]))>3e-14:
            raise ArithmeticError('Signed cell differences failed to close')
        L=np.array([fsum(abs(v)) for v in C.T])
        return dict(x=x,relative_fine_minus_coarse=z.tolist(),sum_abs_cell_changes=L.tolist(),
                    warm_share=[fsum(abs(v))/l for v,l in zip(C[rep>1].T,L,strict=True)])
    return dict(chart='G9',fields=R['fields'],common_native_x05_mass=m,x_start=x,x_end=.05,
                log_pulse_width=width,continuation_steps=N,continuation_step_ratio_to_first_insertion=relative_step,
                late_before_after=[diff(start,x1),diff(data,.05)],first_before_after=R['before_after'],checks=gates,
                input_sha256=[ha,hg,hf,*h],signature_source_sha256=keys,
                scope='One fixed-width perturbation inserted at a different actual point, subsequent updates identical within each pair. Not a global discretization error or physical inference.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for a in ('card','progress','final','receipt','original'):
        p.add_argument(a,type=Path)
    p.add_argument('checkpoints',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
