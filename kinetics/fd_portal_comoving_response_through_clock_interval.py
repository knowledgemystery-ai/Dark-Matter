"""Follow physical initial fixed-total-energy perturbations through the direct fast-rate interval."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_propagator import solve,ordered
from sm_dark_entropy_exchange_reference import StandardModel


def run(eos_path,early,through):
    raw=[p.read_bytes() for p in (early,through)]
    h=[hashlib.sha256(r).hexdigest() for r in raw]
    A,B=map(json.loads,raw)
    table=hashlib.sha256(eos_path.read_bytes()).hexdigest()
    if table!=A['table_input_sha256'] or table!=B['table_input_sha256'] or B['restart_and_clock_sha256'][0]!=h[0]:
        raise ArithmeticError('Use exact direct histories and same accepted EOS')
    x=next(r['points'][-1] for r in A['rows'] if r['operator']=='combined')
    data=next(r['points'] for r in B['rows'] if r['operator']=='combined')
    if (data[0][0],data[-1][0])!=(.55,.96) or max(abs(np.asarray(data[0])-x))>2e-15:
        raise ArithmeticError('Original .55 endpoint must exactly join the later branch')
    E=StandardModel(eos_path)
    P,Q=[CombinedExpansion(E,PositivePortal(m),n) for m,n in (((13,36,26),112),((17,44,32),168))]
    g,_,c0,b0=solve(P,data,4e-11,.016)
    f,nl,c1,b1=solve(Q,data,8e-12,.012)
    def integral(field,ctrl,k,epsilon):
        return ordered(field,ctrl,np.linspace(.55,.96,k),epsilon)
    a,_,_=integral(g,c0,33,6e-4)
    b,_,_=integral(f,c1,65,3e-4)
    C,error,trace=integral(f,c1,129,3e-4)
    err=max(np.max(abs(a-b)),np.max(abs(b-C)))
    if err>7e-5:
        raise ArithmeticError(('Insufficient earlier tangent sampling',err))
    initial=np.ones(3)
    q=2e-4*np.eye(3)[:2]
    D=np.column_stack([(nl(initial+d)(.96)-nl(initial-d)(.96))/(4e-4) for d in q])
    delta=np.max(abs(C[:,:2]-D))
    if delta>9e-6:
        raise ArithmeticError(('Independent fixed-total-energy perturbations disagree',delta))
    endpoint=c1(.96)
    frac=C/endpoint[:,None]
    sig=np.linalg.svd(frac[:,:2],compute_uv=False)
    s=np.linalg.svd(frac[:2,:2],compute_uv=False)
    return dict(operator='combined',interval=[.55,.96],fields=['comoving_number','comoving_dark_energy','comoving_total_energy'],
                initial_fixed_total_energy_columns=2,ordered_response=C.tolist(),independent_perturbed_columns=D.tolist(),
                max_trajectory_archive_defects=[b0,b1],mesh_defects=[np.max(abs(a-b)),np.max(abs(b-C))],
                liouville_abs_defect=error,integrated_trace=trace,independent_perturbed_max_defect=delta,
                final_basis=endpoint.tolist(),fixed_total_initial_to_all_fractional_outputs=frac[:,:2].tolist(),
                singular_values_fixed_total_to_all_fractions=sig.tolist(),
                singular_values_fixed_total_to_dark_fractions=s.tolist(),input_sha256=h,table_input_sha256=table,
                scope='Numerical finite-time memory test of the restricted direct combined Fermi–Dirac history; not a physical stability or decoupling criterion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos_path','early','through'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
