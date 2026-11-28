"""Reduce the exact independently evolved D8 central-band paired quadratures."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import D8_CENTRAL_COLD_BAND, run as topology
from early_pole_endpoint_continuum_cold_band_receipt import read
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_endpoint_scalar_subshell_chain_receipt import read as read_state


def reduce_one(rows, fresh):
    edges=fresh['edges']
    if (len(rows)!=fresh['total'] or [r['panel'] for r in rows]!=list(range(fresh['total']))
        or [r['bounds'] for r in rows]!=[edges[i:i+2] for i in range(len(edges)-1)]
        or any(any(r[k]!=v for k,v in fresh.items() if k!='edges') for r in rows)):
        raise ArithmeticError('Complete ordered exact native D8 panels required')
    c=np.array([r['coeff'] for r in rows])
    if c.shape!=(len(rows),3,3) or not np.all(np.isfinite(c)) or np.min(c)<0:
        raise ArithmeticError('Nonnegative finite coefficient ledgers required')
    led=np.max(abs(c[:,:,0]-c[:,:,1]-c[:,:,2])/np.maximum(np.max(c,axis=2),1e-50))
    z=np.array([[fsum(c[:,i,j]) for j in range(3)] for i in range(3)])
    if min(z.ravel())<=0 or led>4e-13:
        raise ArithmeticError('Integrated totals must be positive and their affine ledger close')
    if max(abs((fresh['scalar_mass']*z[0]+z[2])/z[1]-1.))>2e-14:
        raise ArithmeticError('Independent energy-centering identity failed')
    return z


D4_REFERENCE_SHA256='6a5b61404351598840849795d1c9df36a5a8f14b6ec9264d7bc6c4fe223855c1'


def sign_interval(x):
    return int(x[0]>0.)-int(x[1]<0.)


def decisions(a,b,d4,n):
    v=np.array([a[0],b[0]])
    if not np.isfinite(n) or n<0. or not np.all(np.isfinite(v)) or min(v.ravel())<=0.:
        raise ArithmeticError('Positive finite number ledgers and finite nonnegative occupation required')
    intervals=np.array([d4['forward'],d4['inverse']])
    if (intervals.shape!=(2,2) or not np.all(np.isfinite(intervals))
        or min(intervals.ravel())<=0. or np.any(intervals[:,0]>intervals[:,1])):
        raise ArithmeticError('Ordered positive finite D4 intervals required')
    gross=np.array([v[:,:2].min(axis=0)/intervals[:,1],v[:,:2].max(axis=0)/intervals[:,0]]).T
    neutral=v[:,1]/v[:,2]
    # First evaluate both complete rows separately. Then also use the
    # coefficient box so accidental correlation of the pair is not required
    # for the reported sign relative to the observed numerical variation.
    paired_net=[fsum((n*row[2],-row[1])) for row in v]
    box=[fsum((n*min(v[:,2]),-max(v[:,1]))),
         fsum((n*max(v[:,2]),-min(v[:,1])))]
    if not np.all(np.isfinite(np.r_[gross.ravel(),neutral,box])):
        raise ArithmeticError('Comparison overflow')
    return dict(observed_two_mesh_d8_over_d4_forward=gross[0].tolist(),
                observed_two_mesh_d8_over_d4_inverse=gross[1].tolist(),
                observed_two_mesh_gross_order_vs_d4=[sign_interval(z-1.) for z in gross],
                d8_number_neutral_both_meshes=neutral.tolist(),
                d8_conditional_net_both_meshes=paired_net,
                d8_conditional_net_box_from_observed_coefficients=box,
                d8_conditional_net_sign_from_observed_coefficients=sign_interval(box),
                d4_conditional_net_interval=d4['conditional_net'],
                d4_conditional_net_sign=sign_interval(d4['conditional_net']))


def run(card,subband,raw):
    s,sh=read_state(card)
    if (s['grid'],s['x_end'],s['steps'],s['subcells_per_archived_cell'])!=('fine',.05,2048,8):
        raise ArithmeticError('Independently evolved native fine D8 endpoint required')
    span=D8_CENTRAL_COLD_BAND
    info=[topology(s,m,None,topology_only=True,span=span) for m in (1,2)]
    if info[0]['edges']!=info[1]['edges'] or (info[0]['edges'][0],info[0]['edges'][-1])!=span:
        raise ArithmeticError('Common exact D8 support required')
    a,b=[reduce_one(read(z),g) for z,g in zip(raw,info,strict=True)]
    d=max(abs(a/b-1.).ravel())
    if d>6e-5:
        raise ArithmeticError(('Held central D8 quadrature does not resolve fixed-state source',d))
    h=subband.read_bytes()
    subband_sha256=hashlib.sha256(h).hexdigest()
    if subband_sha256!=D4_REFERENCE_SHA256:
        raise ArithmeticError('Require the exact archived D4 subband enclosure')
    U=json.loads(h)
    p=next(v for v in U['rows'] if tuple(v['q'])==span)
    n=float(U['conditional_photon_occupation'])
    q=p['D4']
    return dict(q=span,checked_mesh=[t['mesh'] for t in info],panels=info[1]['total'],
                coefficient_fields=['forward','inverse','slope'],row_fields=['number','energy','centered_energy'],
                d8_ledger=b.tolist(),d8_max_mesh_relative_difference=d,d8_number_neutral=b[0,1]/b[0,2],
                d8_over_d4_forward=[b[0,0]/q['forward'][1],b[0,0]/q['forward'][0]],
                d8_over_d4_inverse=[b[0,1]/q['inverse'][1],b[0,1]/q['inverse'][0]],
                held_photon_occupation=n,d8_conditional_net_relative_to_inverse=n*b[0,2]/b[0,1]-1.,
                **decisions(a,b,q,n),
                state_digest_uncompressed=sh,subband_sha256=subband_sha256,
                raw_uncompressed_sha256=[unpacked_hash(x) for x in raw],
                scope='Held native D8 band [3e-6,5e-6] and D4 panel enclosure, paired-mesh variation and coefficient boxes describe observed numerical variation, not an error bound on the continuum or complete D8 area.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('card',type=Path)
    p.add_argument('subband',type=Path)
    p.add_argument('raw',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
