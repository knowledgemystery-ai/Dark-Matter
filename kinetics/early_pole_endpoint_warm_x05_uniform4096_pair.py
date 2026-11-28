"""Measure shifted x05 count and pressure effects using both genuine uniform4096 references."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def run(g9,mixed,base,prior):
    raw=[p.read_bytes() for p in (g9,mixed,base,prior)]
    A,B,R,P=map(json.loads,raw)
    h=[hashlib.sha256(s).hexdigest() for s in raw]
    for j,v in enumerate((A,B)):
        if (v['source_sha256'][2]!=h[2] or v['source_sha256'][3]!=P['input_sha256'][j]
            or v['steps']!=4096 or (v['x_start'],v['x_end'],v['chart'])!=(.04,.05,P['charts'][j])
            or v['fields']!=P['fields']):
            raise ArithmeticError('Must use exact previous x05 charts and independently evolved uniform controls')
    if P['input_sha256'][2]!=h[2]:
        raise ArithmeticError('Earlier portability comparison does not bind this baseline')
    c=np.asarray(B['terminal_moments'])/A['terminal_moments']-1.
    X,Y=[np.asarray(v['relative_to_uniform4096']) for v in (A,B)]
    D=(1+c)*(Y-X)/(1+X)
    # The same signed contrast change against the uniform2048 reference must replay.
    target=np.asarray(P['paired_contrast_minus_uniform2048'])
    if np.max(abs(D[:2]-D[2]-target))>3e-16 or np.max(abs(c+D[[0,2]]-R['mixed300_over_G9_by_clock']))>3e-16:
        raise ArithmeticError('Earlier signed response does not replay')
    z=np.stack([X,Y],axis=1)
    return dict(x_start=.04,x_end=.05,rows=A['rows'],fields=A['fields'],charts=P['charts'],uniform4096_matched_contrast=c.tolist(),
                paired_contrast_minus_uniform4096=D.tolist(),relative_to_uniform4096=z.tolist(),
                shifted_over_uniform1024_within_chart=(abs(z[1])/abs(z[0])).tolist(),
                shifted_over_uniform1024_paired=(abs(D[1])/abs(D[0])).tolist(),
                input_sha256=h,
                scope='Exact same-start ordered refinements against finite reference. No physical prediction or asymptotic continuum conclusion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('g9','mixed','base','prior'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
