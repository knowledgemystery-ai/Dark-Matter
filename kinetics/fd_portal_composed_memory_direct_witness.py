"""Directly perturb the whole .55→1.20 branch in its two selected memory directions."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_propagator import solve
from sm_dark_entropy_exchange_reference import StandardModel


def run(eos_path,earlier,later,early_map,late_map,composed):
    raw=[f.read_bytes() for f in (earlier,later,early_map,late_map,composed)]
    E,L,A,B,C=map(json.loads,raw)
    h=[hashlib.sha256(z).hexdigest() for z in raw]
    eos=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if (A['input_sha256'][1]!=h[0] or B['input_sha256']!=[h[1],h[0]] or C['input_sha256']!=h[2:4]
        or E['table_input_sha256']!=eos or L['table_input_sha256']!=eos
        or [A['operator'],B['operator'],C['operator']]!=['combined']*3):
        raise ArithmeticError('Exact five bridge/ordered inputs and published EOS required')
    u=next(r['points'] for r in E['rows'] if r['operator']=='combined')
    v=next(r['points'] for r in L['rows'] if r['operator']=='combined')
    if u[-1]!=v[0] or (u[0][0],u[-1][0],v[-1][0])!=(.55,.96,1.2):
        raise ArithmeticError('Both direct source paths must have the same exact .96 state')
    points=u+v[1:]
    X=np.asarray(C['rows'][-1]['fixed_total_initial_to_all_fractions'])
    _,sv,R=np.linalg.svd(X[:2])
    for i in range(2):
        j=np.argmax(abs(R[i]))
        if R[i,j]<0.:
            R[i]*=-1.
    stage=CombinedExpansion(StandardModel(eos_path),PositivePortal((17,44,32)),168)
    _,nl,base,defect=solve(stage,points,8e-12,.012)
    rows=[]
    for i,f in enumerate(R):
        direction=np.r_[f,0.]
        expected=X@f
        samples=[]
        for e in (4e-4,2e-4):
            u,v=[nl(np.ones(3)+d*e*direction)(1.2) for d in (1.,-1.)]
            z=(u-v)/(2*e*base(1.2))
            samples.append(z)
        delta=np.max(abs(samples[-1]-expected))
        mesh=np.max(abs(samples[0]-samples[1]))
        if delta>9e-6 or mesh>8e-6:
            raise ArithmeticError(('Whole-window finite witness did not replay ordered result',i,delta,mesh))
        rows.append(dict(initial_fractional_direction=direction.tolist(),ordered_fractional_output=expected.tolist(),
                         expected_dark_two_norm=sv[i],finite_half_step_fractional_output=samples[-1].tolist(),
                         finite_dark_two_norm=np.linalg.norm(samples[-1][:2]),absolute_column_defect=delta,
                         finite_perturbation_size_defect=mesh))
    return dict(x_start=.55,x_end=1.2,operator='combined',fields=C['fields'],same_start_reference_path_defect=defect,
                perturbation_sizes=[4e-4,2e-4],rows=rows,input_sha256=h,table_input_sha256=eos,
                scope='Direct finite-epsilon perturbation of one restricted Fermi–Dirac numerical path, not a physical stability or freeze-out statement.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos_path','earlier','later','early_map','late_map','composed'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
