"""Evaluate direct full-window FD secants in the archived physical entropy metrics."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_composed_memory_entropy_metric import gain

D='2026-11-28'
P='fd-portal-'


def run(root):
    deps=('direct-expansion-055-096','direct-expansion-096-120','comoving-response-055-096',
          'direct-comoving-variational-propagator','comoving-matched-response-055-120',
          'comoving-response-entropy-metric','direct-full-window-memory-witness')
    contents=[(root/f'{P}{n}-{D}.json').read_bytes() for n in deps]
    h=[hashlib.sha256(a).hexdigest() for a in contents]
    E,L,A,B,C,M,W=map(json.loads,contents)
    if (W['input_sha256']!=h[:5] or M['input_sha256']!=h[:5] or C['input_sha256']!=h[2:4]
        or not E['table_input_sha256']==L['table_input_sha256']==M['table_input_sha256']==W['table_input_sha256']
        or not M['x_start']==W['x_start']==.55 or not M['rows'][-1]['x_end']==W['x_end']==1.2
        or W['operator']!='combined' or W['fields']!=M['fields'] or W['perturbation_sizes']!=[4e-4,2e-4]):
        raise ArithmeticError('Exact matching entropy, witness and source receipts required')
    rows=W['rows']
    R=np.array([r['initial_fractional_direction'] for r in rows])
    S=np.array([r['finite_half_step_fractional_output'] for r in rows]).T
    O=np.array([r['ordered_fractional_output'] for r in rows]).T
    H0=np.array(M['start_metric'])[:2,:2]
    H1=np.array(M['rows'][-1]['end_metric'])
    J=np.array(C['rows'][-1]['fixed_total_initial_to_all_fractions'])
    R2=R[:,:2]
    if not np.array_equal(R[:,2],np.zeros(2)) or not np.allclose(R2@R2.T,np.eye(2),rtol=0.,atol=5e-15):
        raise ArithmeticError('Fixed-total initial plane and exact orthogonal directions required')
    if not np.allclose(J@R2.T,O,rtol=0.,atol=5e-15):
        raise ArithmeticError('Ordered archived columns disagree with directions')
    if not np.allclose(np.max(abs(S-O),axis=0),[a['absolute_column_defect'] for a in rows],rtol=0.,atol=5e-15):
        raise ArithmeticError('Finite witness column defects do not reproduce')
    J_secant=S@R2
    g= gain(J,H0,H1)
    f= gain(J_secant,H0,H1)
    if not np.allclose(g,M['rows'][-1]['output_full_entropy_metric_gain'],rtol=0.,atol=4e-14):
        raise ArithmeticError('Ordered entropy gain must replay')
    err=float(gain(J_secant-J,H0,H1)[0])
    if max(abs(g-f)) > err+4e-15:
        raise ArithmeticError('Gain comparison violates the spectral norm defect bound')
    sf=np.linalg.svd(J_secant[:2],compute_uv=False)
    return dict(x_start=.55,x_end=1.2,fields=W['fields'],input_sha256=dict(zip(deps,h,strict=True)),
                endpoint_entropy_gains=dict(ordered=g.tolist(),direct_sampled_secant=f.tolist(),
                                            operator_norm_of_difference=err),
                direct_sampled_secant_dark_fractional_gains=sf.tolist(),
                dark_fractional_response_norm_on_ordered_weak_direction=rows[1]['finite_dark_two_norm'],
                ordered_weak_dark_fractional_gain=rows[1]['expected_dark_two_norm'],
                scope='Two independently sampled finite secants represented in the same endpoint norms; not a global nonlinear contraction or continuum convergence bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
