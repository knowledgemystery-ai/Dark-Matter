"""Compare transported clock prescription against the exact shorter-start reference."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(g9,mixed,paired):
    raw=[p.read_bytes() for p in (g9,mixed,paired)]
    A,B,P=map(json.loads,raw)
    h=[hashlib.sha256(z).hexdigest() for z in raw]
    for r,p in zip((A,B),P['rows'],strict=True):
        if (r['x_start'],r['x_end'],r['clock_exponent'],r['chart'])!=(.04,.05,1.2,p['chart']) or r['inputs'][2]!=h[2]:
            raise ArithmeticError('Expected exact independent shorter-start attempt on the prescribed same clock')
        if r['moments'][0]!=p['coarse'] or r['moments'][-1]!=p['fine'] or r['fields']!=P['fields']:
            raise ArithmeticError('Do not change the baseline pressure measurement or uniform reference')
    x,y=[np.array(z['relative_to_uniform2048']) for z in (A,B)]
    old=np.asarray(P['mixed300_over_G9_by_clock'])
    diff=(1.+old[-1])*(y-x)/(1.+x)
    if np.max(abs(diff[0]-(old[0]-old[1])))>3e-16:
        raise ArithmeticError('Original shorter-segment matched clock contrast did not replay')
    return dict(x_start=.04,x_end=.05,fields=P['fields'],charts=['G9','mixed300'],
                rows=['uniform1024','shifted1024'],
                relative_to_uniform2048=np.stack([x,y],axis=1).tolist(),
                shifted_discrepancy_over_uniform1024=(abs(np.stack([x[1],y[1]]))/abs(np.stack([x[0],y[0]]))).tolist(),
                paired_contrast_minus_uniform2048=diff.tolist(),
                paired_shifted_over_uniform1024=(abs(diff[1])/abs(diff[0])).tolist(),input_sha256=h,
                scope='Out-of-interval numerical portability with a single fixed exponent and exact starts; not a global error, continuum estimate, or physical transfer.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('g9','mixed','paired'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
