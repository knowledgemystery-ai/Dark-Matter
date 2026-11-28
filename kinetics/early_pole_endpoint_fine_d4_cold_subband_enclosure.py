"""Compare absolute and conditional signed sources without guessing subpanel shapes."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_receipt import read
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash


def box(bounds, positive, a, b):
    l,r=np.asarray(bounds).T
    y=np.asarray(positive)
    if (y.shape!=l.shape or not np.all(np.isfinite(y)) or np.any(y<0)
        or not np.all(l<r) or not np.array_equal(l[1:],r[:-1])
        or not l[0]<=a<b<=r[-1]):
        raise ArithmeticError('Need contiguous panels and positive measure')
    return [fsum(y[(l>=a)&(r<=b)]),fsum(y[(r>a)&(l<b)])]


def ratio(a,b):
    if b[0]<=0.: return None
    return [a[0]/b[1],a[1]/b[0]]


def sign(interval):
    return int(interval[0]>0.)-int(interval[1]<0.)


def run(ladder,raw):
    data=ladder.read_bytes()
    L=json.loads(data)
    ids=[unpacked_hash(p) for p in raw]
    if ids!=[r[1] for r in L['raw_uncompressed_sha256'][1:]]:
        raise ArithmeticError('Use the exact accepted fine meshes from D2 and D4')
    parts=[read(f) for f in raw]
    n=next(z['common_diagnostic_occupation'] for z in L['diagnostic_root_midpoints'] if z['levels']==[2,4])
    boxes=[np.array([u['bounds'] for u in p]) for p in parts]
    coeffs=[np.array([u['coeff'][0] for u in p]) for p in parts]
    cuts=[1e-6,3e-6,5e-6,8e-6,1e-5,1.5e-5]
    rows=[]
    for a,b in zip(cuts[:-1],cuts[1:],strict=True):
        z=[]
        for x,Y in zip(boxes,coeffs,strict=True):
            A,B,D=[box(x,f,a,b) for f in Y.T]
            # A partial panel has no known signed-shape moment. Bound the
            # two gross contributions separately rather than clipping sign.
            net=[n*D[0]-B[1],n*D[1]-B[0]]
            z.append(dict(forward=A,inverse=B,conditional_net=net,net_sign=sign(net)))
        rows.append(dict(q=[a,b],D2=z[0],D4=z[1],
                         d4_over_d2_forward=ratio(z[1]['forward'],z[0]['forward']),
                         d4_over_d2_inverse=ratio(z[1]['inverse'],z[0]['inverse'])))
    return dict(rows=rows,retained_span=[cuts[0],cuts[-1]],conditional_photon_occupation=n,
                ladder_sha256=hashlib.sha256(data).hexdigest(),raw_uncompressed_sha256=ids,
                scope='Positive-panel containment bounds; no intrapanel interpolation, physical occupation determination, or spatial convergence assertion.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('ladder',type=Path)
    p.add_argument('raw',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
