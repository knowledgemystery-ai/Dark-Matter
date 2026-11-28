"""Compare exact uniform time-clock discrepancies and paired contrasts independently."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(g9,mixed,paired):
    raw=[p.read_bytes() for p in (g9,mixed,paired)]
    G,M,P=map(json.loads,raw)
    h=[hashlib.sha256(x).hexdigest() for x in raw]
    for i,t in enumerate((G,M)):
        if (t['input_hashes'][-1]!=h[2] or t['steps']!=[1024,2048,4096] or t['mass']!=P['mass']
            or t['fields']!=P['fields'] or t['x_start']!=.02 or t['x_end']!=.04
            or t['terminal_moments'][:2]!=[P['rows'][i][j] for j in ('coarse','fine')]):
            raise ArithmeticError('Not two exact same geometry family, start clock, pressure mass and lock inputs')
    A=np.array(G['terminal_moments'])
    B=np.array(M['terminal_moments'])
    ratio=B/A-1.
    if not np.array_equal(ratio[:2],P['mixed300_over_G9_by_clock']):
        raise ArithmeticError('Original two-clock contrast must replay exactly')
    ds=np.stack([A[:-1]/A[1:]-1., B[:-1]/B[1:]-1.],axis=1)
    contrasts=np.diff(ratio,axis=0)
    relative=np.abs(contrasts)/np.max(np.abs(ds),axis=1)
    return dict(x_start=.02,x_end=.04,steps=[1024,2048,4096],fields=P['fields'],
                charts=['G9','mixed300'],mixed300_over_G9_minus_one=ratio.tolist(),
                fine_minus_coarse_contrast=contrasts.tolist(),within_chart_coarse_over_fine_minus_one=ds.tolist(),
                abs_contrast_clock_change_over_larger_within_chart_change=relative.tolist(),
                input_sha256=h,scope='Two genuinely evolved matched-clock segments from each chart’s own exact .02 restart; not full-history convergence or physical error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('g9','mixed','paired'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
