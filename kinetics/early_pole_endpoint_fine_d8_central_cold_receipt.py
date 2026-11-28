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
    U=json.loads(h)
    p=next(v for v in U['rows'] if tuple(v['q'])==span)
    n=float(U['conditional_photon_occupation'])
    q=p['D4']
    if min(q['forward'][0],q['inverse'][0])<=0.:
        raise ArithmeticError('Independent positive D4 central-band interval required')
    return dict(q=span,checked_mesh=[t['mesh'] for t in info],panels=info[1]['total'],
                coefficient_fields=['forward','inverse','slope'],row_fields=['number','energy','centered_energy'],
                d8_ledger=b.tolist(),d8_max_mesh_relative_difference=d,d8_number_neutral=b[0,1]/b[0,2],
                d8_over_d4_forward=[b[0,0]/q['forward'][1],b[0,0]/q['forward'][0]],
                d8_over_d4_inverse=[b[0,1]/q['inverse'][1],b[0,1]/q['inverse'][0]],
                held_photon_occupation=n,d8_conditional_net_relative_to_inverse=n*b[0,2]/b[0,1]-1.,
                state_digest_uncompressed=sh,subband_sha256=hashlib.sha256(h).hexdigest(),
                raw_uncompressed_sha256=[unpacked_hash(x) for x in raw],
                scope='Held native D8 band [3e-6,5e-6] and D4 panel enclosure, not complete D8 area or a physical prediction.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('card',type=Path)
    p.add_argument('subband',type=Path)
    p.add_argument('raw',type=Path,nargs=2)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
