"""Measure comparison-clock changes across the two actual shifted refinements."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(g9,mixed,previous,uniform):
    raw=[p.read_bytes() for p in (g9,mixed,previous,uniform)]
    A,B,P,U=map(json.loads,raw)
    h=[hashlib.sha256(r).hexdigest() for r in raw]
    if [r['chart'] for r in (A,B)]!=[9,'mixed300'] or any(q['fields']!=U['fields'] or q['steps']!=2048 or q['exponent']!=1.2 for q in (A,B)):
        raise ArithmeticError('Require both exact shifted2048 controls')
    for j,q in enumerate((A,B)):
        if (q['source_sha256'][2]!=U['input_sha256'][j] or q['source_sha256'][3]!=P['input_sha256'][j]
            or U['input_sha256'][j]!=P['input_sha256'][2+j]):
            raise ArithmeticError('New records do not bind exact previous schedules and references')
    X=np.array(A['relative_to_uniform4096'])
    Y=np.array(B['relative_to_uniform4096'])
    z=(1+np.asarray(U['mixed300_over_G9_minus_one'][-1]))*(Y-X)/(1+X)
    if np.max(abs(z[:3]-np.asarray(P['relative_paired_contrast_minus_uniform4096'])[[0,2,1]]))>3e-16:
        raise ArithmeticError('Previous matched comparison must replay')
    a=np.stack([X,Y],axis=1)
    return dict(fields=U['fields'],charts=U['charts'],rows=A['rows'],within_chart_relative_to_uniform4096=a.tolist(),
                paired_minus_uniform4096=z.tolist(),
                shifted_over_uniform4096_reference_discrepancy_at_2048=(abs(a[-1])/abs(a[-2])).tolist(),
                shifted_over_uniform_paired_discrepancy_at_2048=(abs(z[-1])/abs(z[-2])).tolist(),
                shifted_1024_to_2048_absolute_discrepancy_ratio=(abs(a[-1])/abs(a[1])).tolist(),
                input_sha256=h,
                scope='Only same-start actual time-clock diagnostics versus fixed uniform 4096, no inferred asymptotic order or continuum bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('g9','mixed','previous','uniform'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
