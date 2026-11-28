"""Audit fixed-mass scalar pressure against the exact full endpoint energy ledger."""
import argparse
import gzip
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import checkpoint_signature
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_charge_symmetric_null_terminal import setup
from early_pole_symmetric_charge_cutoff_audit import DATE, run as cutoff

AMPLITUDE=(0.,-1.,1.,-.5,.5)


def decompose(chart,x,y,z):
    """Stable same-clock difference with the explicit Hartree subtraction."""
    m,E,e0,_=chart.state(x,y)
    u,F,e1,_=chart.state(x,z)
    M=(u-m)*(u+m)
    n=chart.ns
    N=chart.w
    D=N*(z-y)
    s,d,p=[fsum(v) for v in np.split(D*E,(n,chart.iP))]
    qp=chart.ds.p
    p3=fsum(D[:n]*qp*qp/E[:n])
    sc_trace=fsum(D[:n]*(m*m/E[:n]))
    dm=fsum(N[:n]*z[:n]*M/(E[:n]+F[:n]))
    mf=-M*(u*u+m*m-2*(1e-9*x)**2)/(2*chart.ds.lam)
    v=fsum([s,d,p,dm,mf])
    chk=chart.energy_difference(x,y,x,z)
    if abs(fsum([p3,sc_trace,-s]))>2e-19:
        raise ArithmeticError('Fixed-mass kinetic and mass-weighted moments do not close')
    if abs(v-chk)>2e-18 or abs(v-(e1-e0))>9e-16:
        raise ArithmeticError('Independent stable and direct endpoint energy differences must close')
    return dict(scalar_number_shift=fsum(D[:n]),scalar_fixed_mass_energy=s,
                scalar_fixed_mass_three_pressure=p3,scalar_fixed_mass_energy_minus_three_pressure=sc_trace,
                scalar_energy_change_from_pole_only=dm,mean_field_term_shift=mf,
                scalar_on_own_pole_energy=fsum([s,dm]),daughter_energy=d,parent_energy=p,
                charged_pair_energy_sum=fsum([d,p]),total_corrected_energy=v,
                stable_state_total_absolute_defect=abs(v-chk),direct_subtraction_absolute_defect=abs(v-(e1-e0)))


def run(root,cutoff_receipt):
    raw=cutoff_receipt.read_bytes()
    U=cutoff(root)
    if U!=json.loads(gzip.decompress(raw)):
        raise ArithmeticError('Require exact original five-arm cutoff receipt')
    d,a,delta,_=setup(root)
    x0,x1,steps=.04,.05,4096
    t=np.linspace(0.,np.log(x1/x0),steps+1)
    arms=[]
    for v,b in zip(U['inputs'],AMPLITUDE,strict=True):
        stem=root/f"early-pole-{v['name']}-x05"
        receipt=Path(f'{stem}-receipt-{DATE}.json').read_bytes()
        content=Path(f'{stem}-checkpoint-{DATE}.json.gz').read_bytes()
        if (hashlib.sha256(receipt).hexdigest()!=v['receipt_sha256']
            or hashlib.sha256(content).hexdigest()!=v['compressed_checkpoint_sha256']):
            raise ArithmeticError('Input hashes must match verified cutoff receipt')
        A,B=json.loads(receipt),json.loads(gzip.decompress(content))
        start=a+b*delta
        if B['signature']!=checkpoint_signature(d,start,t,x0,x1,steps,step,A['provenance']):
            raise ArithmeticError('Saved ledger must be bound to exact symmetric initial condition')
        tr=np.asarray(B['trace'])
        q=np.asarray(B['occupation'])
        if tr.shape!=(steps,) or not np.isfinite(tr).all():
            raise ArithmeticError('Complete finite trace required')
        r=B['checks'][-1]['continuity_residual']
        if abs(d.energy_difference(x0,start,x1,q)-fsum(tr)-r)>1e-18:
            raise ArithmeticError('Each saved continuity-ledger residual must independently replay')
        arms.append((start,q,tr,r))
    out=[]
    for card,b,Q in zip(U['inputs'][1:],AMPLITUDE[1:],arms[1:],strict=True):
        a0,a1,V0,R0=arms[0]
        b0,b1,V1,R1=Q
        D=decompose(d,x1,a1,b1)
        initial=d.energy_difference(x0,a0,x0,b0)
        tr=V1-V0
        work=fsum(tr)
        eps=D['total_corrected_energy']-initial-work
        if abs(eps-(R1-R0))>2e-18 or abs(initial)>1e-16:
            raise ArithmeticError('Initial-energy and between-arm trace closure failed')
        s0=d.w[:d.ns]*a1[:d.ns]
        m,F,*_=d.state(x1,a1)
        normalized=D['scalar_fixed_mass_three_pressure']/fsum(s0*d.ds.p**2/F[:d.ns])
        out.append(dict(arm=card['name'],signed_amplitude=b,**D,initial_total_corrected_energy=initial,
                        final_fixed_mass_nodal_pressure_relative=normalized,
                        cumulative_saved_trace_differences=[fsum(tr[:n]) for n in range(steps//8,steps+1,steps//8)],
                        final_total_minus_initial_minus_trace=eps,
                        difference_of_archived_continuity_residuals=R1-R0))
    return dict(rows=out,clock=[x0,x1,steps],nominal_corrected_energy=d.state(x1,arms[0][1])[2],
                cutoff_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Same native terminal quadrature and archived trace; stable differences of a closed numerical model, not a separate physical energy transfer or stress-energy prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('cutoff_receipt',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.cutoff_receipt),indent=2,allow_nan=False))
