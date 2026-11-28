"""Refine the weak full-window FD response in quadrature, integrator and pulse size."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_propagator import solve
from sm_dark_entropy_exchange_reference import StandardModel

D='2026-11-28'


def run(table,root):
    dep=('direct-expansion-055-096','direct-expansion-096-120','direct-full-window-memory-witness')
    blobs=[(root/f'fd-portal-{s}-{D}.json').read_bytes() for s in dep]
    hashes=[hashlib.sha256(b).hexdigest() for b in blobs]
    E,L,W=map(json.loads,blobs)
    eos=hashlib.sha256(table.read_bytes()).hexdigest()
    if W['input_sha256'][:2]!=hashes[:2] or not E['table_input_sha256']==L['table_input_sha256']==W['table_input_sha256']==eos:
        raise ArithmeticError('Exact existing trajectories, weak direction, and EOS required')
    p,q=[next(t['points'] for t in r['rows'] if t['operator']=='combined') for r in (E,L)]
    if p[-1]!=q[0] or (p[0][0],q[0][0],q[-1][0])!=(.55,.96,1.2):
        raise ArithmeticError('Exact shared .96 state required')
    d=W['rows'][1]
    f=np.array(d['initial_fractional_direction'])
    if f[2]!=0 or abs(f@f-1.)>2e-15:
        raise ArithmeticError('Initial total energy fixed and direction normalized')
    stage=CombinedExpansion(StandardModel(table),PositivePortal((20,52,38)),192)
    _,flows,base,bgdef=solve(stage,p+q[1:],3e-12,.008)
    u=base(1.2)
    rows=[]
    for e in (2e-4,1e-4):
        left,right=[flows(np.ones(3)+np.sign(k)*e*f)(1.2) for k in (-1,+1)]
        z=(right-left)/(2*e*u)
        rows.append(dict(epsilon=e,central_fractional_column=z.tolist(),
                         norm_dark=np.linalg.norm(z[:2]),early_dark_energy_fixed_total_direction=f.tolist(),
                         endpoint_pair_midpoint_fractional_displacement=((right+left)/(2*u)-1.).tolist()))
    a,b=[np.array(r['central_fractional_column']) for r in rows]
    nominal=np.array(d['finite_half_step_fractional_output'])
    ordered=np.array(d['ordered_fractional_output'])
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ArithmeticError('Refined independent perturbations must be finite')
    return dict(rows=rows,table_input_sha256=eos,input_sha256=dict(zip(dep,hashes,strict=True)),
                baseline_path_defect=bgdef,quadrature_settings=[20,52,38,192],rtol=3e-12,max_step=.008,
                same_pulse_new_minus_nominal=(a-nominal).tolist(),halved_pulse_change=(b-a).tolist(),
                finest_minus_ordered=(b-ordered).tolist(),scope='Independent weak-direction finite secants on one restricted FD history. Combined settings contrast is not an isolated quadrature or temporal error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('table',type=Path)
    p.add_argument('root',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
