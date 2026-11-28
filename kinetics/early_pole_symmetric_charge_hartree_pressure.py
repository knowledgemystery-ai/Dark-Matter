"""Distinguish the fixed-mass pressure moment from the model Hartree pressure."""
import argparse
import gzip
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import X
from early_pole_endpoint_charge_symmetric_null_terminal import setup
from early_pole_symmetric_charge_cutoff_audit import DATE
from early_pole_symmetric_charge_terminal_energy import run as energy_audit


def moments(d,x,u):
    N=d.w*u
    m,E,rho,_=d.state(x,u)
    q=np.r_[d.ds.p,d.qd,d.qp]
    bare2=(1e-9*x)**2
    U=(m*m-bare2)**2/(2*d.ds.lam)
    P=[fsum(v) for v in np.split(N*q*q/(3*E),(d.ns,d.iP))]
    T=d.trace(np.log(x/X),u)
    gap= m*m-bare2-.5*d.ds.lam*fsum(N[:d.ns]/E[:d.ns])
    if abs(gap)>5e-18 or abs(fsum([rho,-3*fsum([*P,U]),-T]))>1e-15:
        raise ArithmeticError('Hartree gap and model trace identity must both hold')
    h=2e-4
    def rate(h):
        return d.energy_difference(x*np.exp(-h),u,x*np.exp(h),u)/(2*h)
    check=abs((4*rate(h/2)-rate(h))/3-T)
    if check>1e-12:
        raise ArithmeticError('Independent fixed-occupation expansion derivative failed')
    return dict(m=m,E=E,q=q,N=N,rho=rho,P=P,U=U,T=T,trace_derivative_defect=check,gap=gap)


def run(root,cutoff,energy):
    raw=energy.read_bytes()
    witness=energy_audit(root,cutoff)
    if witness!=json.loads(raw):
        raise ArithmeticError('Must replay exact terminal energy and continuity audit')
    d,*_=setup(root)
    x=.05
    old=json.loads(gzip.decompress((root/f'early-pole-charge-symmetric-time-nominal-x05-checkpoint-{DATE}.json.gz').read_bytes()))
    f0=np.asarray(old['occupation'])
    A=moments(d,x,f0)
    ref=fsum([A['P'][0],A['U']])
    out=[]
    for R in witness['rows']:
        f=np.asarray(json.loads(gzip.decompress((root/f"early-pole-{R['arm']}-x05-checkpoint-{DATE}.json.gz").read_bytes()))['occupation'])
        B=moments(d,x,f)
        sl=slice(None,d.ns)
        m,u=A['m'],B['m']
        a,b=A['E'][sl],B['E'][sl]
        mass_term=fsum(B['N'][sl]*A['q'][sl]**2/3*((m-u)*(m+u))/(a*b*(a+b)))
        comp=np.array([fsum(v) for v in np.split(d.w*(f-f0)*A['q']**2/(3*A['E']),(d.ns,d.iP))])
        mf=-R['mean_field_term_shift']
        scalar=fsum([comp[0],mass_term,mf])
        total=fsum([*comp,mass_term,mf])
        z=B['T']-A['T']
        trace_defect=abs(fsum([R['total_corrected_energy'],-3*total,-z]))
        if (trace_defect>1e-15 or abs(comp[0]*3-R['scalar_fixed_mass_three_pressure'])>5e-21
            or abs(total-fsum([*(np.array(B['P'])-A['P']),B['U']-A['U']]))>2e-16):
            raise ArithmeticError('Pressure decomposition disagrees with native trace or direct arithmetic')
        out.append(dict(arm=R['arm'],count_shift=R['scalar_number_shift'],
                        fixed_nominal_mass_scalar_kinetic_pressure=comp[0],
                        change_from_scalar_kinetic_pole_only=mass_term,mean_field_pressure=mf,
                        corrected_scalar_pressure=scalar,corrected_scalar_pressure_relative=scalar/ref,
                        daughter_pressure=comp[1],parent_pressure=comp[2],corrected_total_model_pressure=total,
                        corrected_total_energy=R['total_corrected_energy'],bare_trace_change=z,
                        differential_trace_identity_abs_defect=trace_defect,
                        expansion_derivative_against_trace_abs_defect=B['trace_derivative_defect'],
                        hartree_gap_equation_abs_defect=abs(B['gap'])))
    return dict(rows=out,x=x,nominal_corrected_scalar_pressure=ref,
                nominal_corrected_total_pressure=fsum([*A['P'],A['U']]),
                nominal_expansion_derivative_against_trace_abs_defect=A['trace_derivative_defect'],
                energy_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Differentiable, vacuum-subtracted Hartree pressure for the exact held numerical model, not a physical plasma equation of state or convergence result.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','cutoff','energy'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
