"""Separate the rank-one amplitude change from the remaining EOS rate form."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_cubic_derivative_cross import slope_and_rank
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_entropy_knot_jump_audit import forms
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_heat_capacity_slope_threshold import DATE, threshold
from sm_dark_entropy_exchange_reference import StandardModel


def cross_thresholds(kp,kc,r):
    if min(kp,kc)<=3. or not np.isfinite(r) or r<=0.:
        raise ArithmeticError('Expected positive rank amplitudes and finite threshold increments')
    return [[kp,3.+(kp-3.)/r],[3.+(kc-3.)*r,kc]]


def density(P,p):
    z=P.bath(p[2:4])
    return np.array([z['n'],z['rho'],P.eos.state(p[1])[0]+z['rho']])


def run(eos,root):
    names=('direct-expansion-096-120','entropy-pchip-knot-crossing','cubic-derivative-cross')
    raws=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    L,K,A=map(json.loads,raws)
    ids=[hashlib.sha256(s).hexdigest() for s in raws]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (K['input_sha256'][0]!=ids[0] or A['input_sha256'][1]!=ids[0] or A['input_sha256'][5]!=ids[1]
        or any(x['table_input_sha256']!=h for x in (L,K,A))):
        raise ArithmeticError('The exact previous crossed-derivative diagnostic is required')
    P,Q=[CombinedExpansion(f(eos),PositivePortal((17,44,32)),168) for f in (StandardModel,SplineEOS)]
    hist=next(r['points'] for r in L['rows'] if r['operator']=='combined')
    walk=local_path(P,next(p for p in hist if p[0]==.99))
    out=[]
    for row in A['rows'][2:]:
        p=walk(row['x'])
        if abs(p[1]-row['T'])>3e-14:
            raise ArithmeticError('Held temperature must match')
        D=np.diag(density(P,p)/density(Q,p))
        G0,S0=forms(P,p)
        G1,S1=forms(Q,p)
        k0,R0=slope_and_rank(P,p)
        k1,R1=slope_and_rank(Q,p)
        G1,S1,R1=[D.T@X@D for X in (G1,S1,R1)]
        ratio=float(np.sum(R0*R1)/np.sum(R0*R0))
        rank_error=np.linalg.norm(R1-ratio*R0,2)/np.linalg.norm(R1,2)
        if rank_error>2e-13:
            raise ArithmeticError('Rank direction must agree in the common density coordinates')
        cv_ratio=Q.eos.state(p[1])[3]/P.eos.state(p[1])[3]
        kp,kc=[row[f'base_{x}']['critical_heat_capacity_slope'] for x in ('pchip','cubic')]
        M=cross_thresholds(kp,kc,ratio)
        checks=[]
        for G,S,k,R in ((G0,S0,k0,R0),(G1,S1,k1,R1)):
            B=S+(3.-k)*R
            checks.append([threshold(G,B,Rj,3.)['critical_slope_above_three'] for Rj in (R0,R1)])
        if max(abs(np.array(checks)-M).ravel())>6e-10:
            raise ArithmeticError('Proportional-rank formulas must match direct crossed eigensolves')
        out.append(dict(x=p[0],T=p[1],cubic_over_pchip_rank_coefficient=ratio,
                        cubic_over_pchip_heat_capacity=cv_ratio,
                        cubic_over_pchip_cooling_magnitude=ratio*cv_ratio,
                        common_rank_direction_relative_defect=rank_error,
                        critical_slope_by_base_rows_and_rank_columns=M,
                        only_rank_amplitude_shift_at_pchip_base=M[0][1]-M[0][0],
                        only_base_shift_at_pchip_rank=M[1][0]-M[0][0],
                        native_full_critical_slope_shift=kc-kp))
    return dict(rows=out,input_sha256=ids,table_input_sha256=h,
                scope='Local two-factor change in common density coordinates. Splits the rank-one amplitude from the rest of the held rate form, not a causal separation of all EOS or physical effects.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
