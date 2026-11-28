"""Check whether the shared shifted schedule preserves G9/mixed300 contrast accuracy."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(g9,mixed,ug,um,pair):
    raw=[p.read_bytes() for p in (g9,mixed,ug,um,pair)]
    hashes=[hashlib.sha256(b).hexdigest() for b in raw]
    G,M,A,B,R=map(json.loads,raw)
    for i,(p,q) in enumerate(((G,A),(M,B))):
        if (p['x_start']!=.02 or p['x_end']!=.04 or p['steps']!=1024 or p['exponent']!=1.2
            or p['input_sha256'][:2]!=q['input_hashes'][:2] or p['input_sha256'][2]!=hashes[2+i]
            or R['input_sha256'][i]!=hashes[2+i] or p['fields']!=R['fields']):
            raise ArithmeticError('Requires the same prescribed clock and locked uniform returns')
        stored=np.array(q['terminal_moments'])
        if np.max(abs(np.asarray(p['relative_to_uniform4096'][:2])-(stored[:2]/stored[2]-1.)))>2e-16:
            raise ArithmeticError('Different 4096 time reference or reference measurement')
    if G['relative_min_max_step']!=M['relative_min_max_step']:
        raise ArithmeticError('The two actual shifted schedules need to have equal widths')
    a,b=(np.array(p['relative_to_uniform4096']) for p in (G,M))
    baseline=np.array(R['mixed300_over_G9_minus_one'][-1])
    # Subtract the paired high-clock ratio without cancellation from repeatedly subtracting 1.
    changes=(1+baseline)*(b-a)/(1+a)
    if np.max(abs(changes[:2]-(np.array(R['mixed300_over_G9_minus_one'][:2])-baseline)))>3e-16:
        raise ArithmeticError('Pair relation does not recover either exact older clock')
    return dict(x_start=.02,x_end=.04,fields=R['fields'],rows=G['rows'],
                uniform4096_contrast=baseline.tolist(),
                relative_paired_contrast_minus_uniform4096=changes.tolist(),
                within_chart_relative_changes_from_uniform4096=np.stack([a,b],axis=1).tolist(),
                shifted_paired_discrepancy_over_uniform1024=(abs(changes[-1])/abs(changes[0])).tolist(),
                input_sha256=hashes,
                scope='Matched fixed geometry and fixed starts only; unsigned gain uses the uniform 4096 reference, not a continuum or physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('g9','mixed','ug','um','pair'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
