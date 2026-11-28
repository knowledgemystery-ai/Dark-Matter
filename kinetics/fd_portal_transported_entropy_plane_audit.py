"""Audit the measured mesh contrast via entropy-metric principal angles."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_entropy_plane_projection_bound import projection_rate_bound_from_fractional_sine
from fd_portal_transported_entropy_plane import forms
from sm_dark_entropy_exchange_reference import StandardModel


def run(eos,early,late,witness):
    inputs=[p.read_bytes() for p in (early,late,witness)]
    H=[hashlib.sha256(p).hexdigest() for p in inputs]
    A,B,W=map(json.loads,inputs)
    s=hashlib.sha256(eos.read_bytes()).hexdigest()
    if W['input_sha256'][:2]!=H[:2] or W['table_input_sha256']!=s:
        raise ArithmeticError('Require exact transported-plane witness')
    p,q=[next(z['points'] for z in R['rows'] if z['operator']=='combined') for R in (A,B)]
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    rows=[]
    for z in W['rows']:
        x=z['x']
        S=next(t for t in p+q[1:] if t[0]==x)
        G,C=forms(P,S)
        U,R=np.linalg.qr(np.asarray(z['endpoint_fractional_plane_columns']),mode='reduced')
        rates=eigh(U.T@C@U,U.T@G@U,eigvals_only=True)
        if min(abs(R.diagonal()))<=1e-10 or max(abs(rates-z['transport_fixed_initial_total_plane_squared_norm_rate']))>3e-10:
            raise ArithmeticError('Need a resolved two-dimensional plane that replays the reported result')
        sine,delta=projection_rate_bound_from_fractional_sine(G,C,z['max_projector_discrepancy_two_norm'])
        if delta<z['max_rate_discrepancy_between_meshes'] or abs(rates[-1])<=delta:
            raise ArithmeticError('Finite numerical planes do not certify the reported same-sign contrast')
        rows.append(dict(x=x,replayed_fine_rates=rates.tolist(),entropy_sine_bound_for_computed_planes=sine,
                         rate_change_bound_for_computed_planes=delta,verified_two_mesh_maximum_sign=int(np.sign(rates[-1]))))
    return dict(rows=rows,input_sha256=H,table_input_sha256=s,
                scope='Bounds the effect of the observed finite mesh plane difference, not the unknown remaining discretization error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','early','late','witness'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
