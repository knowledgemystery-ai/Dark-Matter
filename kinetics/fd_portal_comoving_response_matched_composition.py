"""Compose two actual FD tangent intervals in the shared fractional coordinates."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(early,late):
    raw=[p.read_bytes() for p in (early,late)]
    A,B=map(json.loads,raw)
    if (A['operator']!=B['operator'] or (A['interval'],B['interval'])!=([.55,.96],[.96,1.2])
        or A['fields']!=B['fields'] or A['input_sha256'][1]!=B['input_sha256'][1]
        or A['table_input_sha256']!=B['table_input_sha256'] or A['initial_fixed_total_energy_columns']!=2):
        raise ArithmeticError('Both exact charts must share the .96 bridge and FD EOS')
    X=np.array(A['ordered_response'])/np.array(A['final_basis'])[:,None]
    Y=np.array(B['ordered_response'])/np.array(B['final_baseline_in_start_units'])[:,None]
    Q=np.array(A['independent_perturbed_columns'])/np.array(A['final_basis'])[:,None]
    R=np.array(B['finite_perturbation_response'])/np.array(B['final_baseline_in_start_units'])[:,None]
    C=Y@X[:,:2]
    D=R@Q
    e=np.max(abs(C-D))
    if e>=5e-6 or not np.array_equal(X[:,:2],A['fixed_total_initial_to_all_fractional_outputs']):
        raise ArithmeticError('Both independently perturbed factorizations should reproduce the match')
    out=[]
    for t,a,b in ((.96,X[:,:2],Q),(1.2,C,D)):
        v=[np.linalg.svd(z,compute_uv=False).tolist() for z in (a,a[:2],b,b[:2])]
        out.append(dict(x_end=t,fixed_total_initial_to_all_fractions=a.tolist(),independent_finite_perturbed_composition=b.tolist(),
                         output_gain_all=v[0],output_gain_dark_only=v[1],independent_gain_all=v[2],independent_gain_dark_only=v[3]))
    return dict(x_start=.55,x_bridge=.96,operator=A['operator'],fields=A['fields'],rows=out,
                maximum_composed_column_defect=e,input_sha256=[hashlib.sha256(u).hexdigest() for u in raw],
                scope='Two ordered numerical tangent maps with shared fractional geometry; independent composition uses separately integrated local finite perturbed columns, not a global nonlinear sensitivity experiment or physical stability assertion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('early','late'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
