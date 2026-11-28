"""Transform exact fixed-total-energy memory to eta/log-T coordinates at both ends."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from sm_dark_entropy_exchange_reference import StandardModel


def pullback(stage,p):
    e,l=p[2:4]
    b=stage.bath([e,l])
    N,R,B,M,C=[b[k] for k in ('n','rho','dn_deta','chemical_energy','cv_fixed_number')]
    T=np.exp(l)
    A=np.array([[B,B*M/T],[B*M,T*C+B*M*M/T]])/np.array([N,R])[:,None]
    E=2e-5
    def flux(x):
        u=stage.bath(x)
        return np.array([u['n']/N,u['rho']/R])
    h=np.array([e,l])
    D=np.column_stack([(flux(h-2*d)-8*flux(h-d)+8*flux(h+d)-flux(h+2*d))/(12*E)
                       for d in E*np.eye(2)])
    err=max(np.max(abs(A-D)),0.)
    if err>3e-9 or np.linalg.cond(A)>1000.:
        raise ArithmeticError('Thermodynamic chain rule did not independently resolve')
    return A,err


def run(eos_path,earlier,late,emap,lmap,comp):
    raw=[a.read_bytes() for a in (earlier,late,emap,lmap,comp)]
    A,B,C,D,E=map(json.loads,raw)
    h=[hashlib.sha256(a).hexdigest() for a in raw]
    eos=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if (C['input_sha256'][1]!=h[0] or D['input_sha256']!=[h[1],h[0]] or E['input_sha256']!=h[2:4]
        or A['table_input_sha256']!=eos or B['table_input_sha256']!=eos):
        raise ArithmeticError('Need exact shared FD EOS and direct archived path')
    left=next(r['points'] for r in A['rows'] if r['operator']=='combined')
    right=next(r['points'] for r in B['rows'] if r['operator']=='combined')
    Q=CombinedExpansion(StandardModel(eos_path),PositivePortal((17,44,32)),168)
    J0,e0=pullback(Q,left[0])
    rows=[]
    for t,z in zip((left[-1],right[-1]),E['rows'],strict=True):
        J,e=pullback(Q,t)
        v=np.asarray(z['fixed_total_initial_to_all_fractions'])[:2]
        q=np.asarray(z['independent_finite_perturbed_composition'])[:2]
        W,H=[np.linalg.solve(J,x@J0) for x in (v,q)]
        delta=np.max(abs(W-H))
        if delta>9e-6:
            raise ArithmeticError('Independent perturbed columns do not support intensive-coordinate audit')
        rows.append(dict(x_end=t[0],fractional_gain=np.linalg.svd(v,compute_uv=False).tolist(),
                         eta_logT_gain=np.linalg.svd(W,compute_uv=False).tolist(),
                         independent_eta_logT_gain=np.linalg.svd(H,compute_uv=False).tolist(),
                         eta_logT_jacobian=W.tolist(),maximum_composed_column_defect=delta,
                         endpoint_thermodynamic_derivative_defect=e))
    return dict(x_start=left[0][0],initial_total_energy_held=True,fields=['eta_dark','log_T_dark'],rows=rows,
                start_thermodynamic_derivative_defect=e0,input_sha256=h,table_input_sha256=eos,
                scope='Coordinate/norm sensitivity of the same fixed-start numerical tangent, not a change in physical dynamics or stability criterion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos_path','earlier','late','emap','lmap','comp'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
