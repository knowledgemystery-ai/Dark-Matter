"""An independent raw point control in the disjoint upper D8 aperture."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_native_cold_cross_factors import DATE, K, MESHES, angle, native, transverse_real_pauli
from early_pole_native_cold_retained_aperture import aperture
from early_pole_native_d8_retained_aperture import run as held


def select(target):
    a,b=[r['raw_median_possible'] for r in target['fine_mesh']['charts']]
    l,h=max(a[0],b[0]),min(a[1],b[1])
    if not l<h:
        raise ArithmeticError('A nonempty common median enclosure is required')
    return (l+h)/2,[l,h]


def run(root):
    names=('native-cold-weight-target','native-d8-raw-band-plan','native-d8-retained-aperture')
    p=[root/f'early-pole-{k}-{DATE}.json' for k in names]
    raw=[z.read_bytes() for z in p]
    T,P,H=map(json.loads,raw)
    if held(root,p[1])!=H or P['previous_target_sha256']!=hashlib.sha256(raw[0]).hexdigest():
        raise ArithmeticError('Replay exact independent geometry and median selection')
    q,window=select(T)
    if not P['planned']['span'][0]<q<P['planned']['span'][1] or q<=5e-6:
        raise ArithmeticError('The check must be disjoint from the live central band')
    rows=[]
    for level in (2,4,8):
        c=native(root,level)
        local=[]
        E=np.hypot(c['scalar_mass'],q)
        args=(c['charges'],K,c['omega'],c['residue'],E,'T')
        for m in MESHES:
            s,t=[angle(*args,m[0],m[1:],q,occupation_callbacks=c['callbacks'],**kw)
                 for kw in ({},dict(mirror_azimuth=True,transverse_vertex_kernel=transverse_real_pauli))]
            v=np.asarray(s[0])
            err=float(np.max(abs(t[0]/v-1.)))
            if min(v)<=0. or not np.isfinite(v).all() or err>2e-12 or not np.array_equal(s[1],t[1]) or s[2:]!=t[2:]:
                raise ArithmeticError('Independent upper angular evaluations did not resolve')
            local.append(dict(mesh=list(m),stock=v.tolist(),real_pauli_mirror=t[0].tolist(),
                              maximum_algorithm_relative_difference=err,angular_support_counts=s[1].tolist(),
                              angular_panels=s[2]))
        change=float(np.max(abs(np.asarray(local[0]['stock'])/local[1]['stock']-1.)))
        if change>1e-5:
            raise ArithmeticError('Held point must resolve on both angular meshes')
        fraction=aperture(q,c,float(c['charges'][1].E[-1]),1)[0]
        rows.append(dict(native_level=level,native_state_sha256_uncompressed=c['card_hash'],
                         inverse_geometric_cosine_fraction=fraction,
                         meshes=local,maximum_observed_mesh_relative_difference=change))
    a,b,c=[np.asarray(z['meshes'][-1]['stock']) for z in rows]
    return dict(q=q,intersection_of_raw_median_enclosures=window,results=rows,
                d4_over_d2=(b/a).tolist(),d8_over_d4=(c/b).tolist(),
                input_sha256=[hashlib.sha256(z).hexdigest() for z in raw],
                scope='One selected fixed native upper q and k=.0015; direct raw angular coefficients and exact geometry, not a panel or full-band result or spatial convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
