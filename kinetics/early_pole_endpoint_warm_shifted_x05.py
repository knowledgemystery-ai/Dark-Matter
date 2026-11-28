"""Apply the prespecified shifted clock on an independent restart interval."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import checkpoint_signature,drive_checkpoint_digest,write_drive_checkpoint
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_gauss_warm_comparison import chart_for,gauss
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart,SPECS,load
from early_pole_live_bounded_face_extents_pilot import hold_trace
from early_pole_resolved_parent_reference import S


def clock(n):
    # Fixed before seeing this candidate's result; shift resolution toward the late segment.
    b=1.2
    L=np.log(.05/.04)
    z=np.log1p(np.expm1(b)*np.linspace(0.,1.,n+1))*(L/b)
    z[0]=0.; z[-1]=L
    if np.min(np.diff(z))<=0.:
        raise ArithmeticError('Prescribed time grid is not strictly positive')
    return z


def integrate(d,h,c,file,key):
    x0,x1=.04,.05
    n=len(c)-1
    if h.shape!=(d.n,) or np.any(~np.isfinite(h)) or min(h)<0. or max(h[d.ns:])>1.:
        raise ArithmeticError('The exact initial state is not eligible')
    signature=checkpoint_signature(d,h,c,x0,x1,n,step,key)
    start=0
    y=h.copy()
    trace,caps,checks=[],[],[]
    if file.is_file():
        z=json.loads(file.read_bytes())
        start=z['completed']
        y=np.asarray(z['occupation'])
        trace,caps,checks=[z[k] for k in ('trace','caps','checks')]
        if (not 1<=start<=n or z['signature']!=signature or z['target_steps']!=n
            or z['payload_sha256']!=drive_checkpoint_digest(z) or z['x_reached']!=x0*np.exp(c[start])
            or len(trace)!=start or len(caps)!=start or len(checks)!=start//(n//8)
            or y.shape!=h.shape or np.any(~np.isfinite(np.r_[y,trace,caps])) or min(y)<0. or max(y[d.ns:])>1.):
            raise ArithmeticError('Invalid shifted-clock checkpoint')
    for i,(a,b) in enumerate(zip(c[start:-1],c[start+1:],strict=True),start+1):
        mid=.5*(a+b)
        x=x0*np.exp(mid)
        M,_,_,H=d.state(x,y)
        mesh,_=d.chart(x,M,3)
        pred,l1=step(mesh,d,y,H*(b-a)/2.)
        M,_,_,H=d.state(x,pred)
        mesh,_=d.chart(x,M,3)
        next_y,l2=step(mesh,d,y,H*(b-a),rate_state=pred)
        trace.append(hold_trace(d,x0,a,mid,y)+hold_trace(d,x0,mid,b,next_y))
        caps.append(max(l1,l2))
        y=next_y
        if i%(n//8)==0:
            u,v,w=[fsum(k) for k in np.split(d.w*(y-h),(d.ns,d.iP))]
            pair=max(abs(u-v),abs(u+w))
            x2=x0*np.exp(b)
            balance=d.energy_difference(x0,h,x2,y)-fsum(trace)
            if pair>3e-13 or abs(balance)>3e-10 or min(y)<0. or max(y[d.ns:])>1.:
                raise ArithmeticError('Shifted time grid violates the existing number, continuity or Pauli gates')
            checks.append(dict(x=float(x2),scalar_yield=fsum(d.w[:d.ns]*y[:d.ns])/S,
                               pair_number_error=pair,continuity_residual=balance,max_limited_extent=float(max(caps))))
        if i%128==0:
            write_drive_checkpoint(file,dict(signature=signature,completed=i,target_steps=n,x_reached=float(x0*np.exp(b)),
                                             occupation=y.tolist(),trace=trace,caps=caps,checks=checks))
    return y,checks,signature


def run(previous,fine,paired,checkpoint,control=None):
    (a,ha),(b,hb)=map(read,(previous,fine))
    raw=paired.read_bytes(); p=json.loads(raw)
    m=p['mass']
    tag,d,_=chart_for(b,budget_control=control)
    name='G9' if tag==9 else tag
    if (name not in ('G9','mixed300') or (a['x_end'],b['x_end'])!=(.04,.05)
        or b['steps']!=2048 or b['actual_previous_raw_sha256']!=ha or p['x']!=.05 or p['steps']!=[1024,2048]):
        raise ArithmeticError('Require exact accepted x04 restart and independent x05 uniform reference')
    q=next(i for i in p['rows'] if i['chart']==name)
    if [ha,hb]!=[q['sha256'][0],q['sha256'][2]]:
        raise ArithmeticError('Archived branch does not match the x05 uniform receipt')
    h=hashlib.sha256(raw).hexdigest()
    c=hashlib.sha256(control.read_bytes()).hexdigest() if control else None
    keys=[hashlib.sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),Path(step.__code__.co_filename),Path(hold_trace.__code__.co_filename))]
    lock=hashlib.sha256(''.join([ha,hb,h,*([c] if c else []),*keys]).encode()).hexdigest()
    n=1024
    t=clock(n)
    y,checks,signature=integrate(d,np.asarray(a['nodal_occupation']),t,checkpoint,lock)
    z=json.loads(checkpoint.read_bytes())
    if z['completed']!=n or z['checks']!=checks or z['signature']!=signature or z['payload_sha256']!=drive_checkpoint_digest(z) or len(checks)!=8:
        raise ArithmeticError('Complete eligible replay required')
    base=ScalarSubshellChart(load('fine')[0],SPECS['fine'],2,4,1e-5)
    def measure(d,n):
        return gauss(d,base,m,.05,expected_steps=n,budget_control=control)[0]
    fine_held=measure(b,2048)
    if not np.array_equal(fine_held,q['fine']):
        raise ArithmeticError('Uniform reference does not replay')
    ans=measure(dict(b,steps=n,nodal_occupation=y.tolist()),n)
    r=np.array([q['coarse'],ans])/fine_held-1.
    return dict(x_start=.04,x_end=.05,chart=name,fields=p['fields'],clock_exponent=1.2,rows=['uniform1024','shifted1024'],
                relative_to_uniform2048=r.tolist(),shifted_over_uniform1024_abs_discrepancy=(abs(r[1])/abs(r[0])).tolist(),
                moments=[q['coarse'],ans.tolist(),fine_held.tolist()],checks=checks,
                inputs=[ha,hb,h,*([c] if c else [])],checkpoint_payload_sha256=z['payload_sha256'],provenance=lock,
                scope='Independent short-interval portability check at unchanged shift exponent; uniform 2048 is a reference only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for a in ('previous','fine','paired','checkpoint'):
        p.add_argument(a,type=Path)
    p.add_argument('--control',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
