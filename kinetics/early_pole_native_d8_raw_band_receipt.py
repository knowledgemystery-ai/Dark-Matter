"""Validate complete paired ledgers for a future exact D8 raw-band result."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import run as native_topology
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import address, ledger, strip
from early_pole_native_cold_cross_factors import CellChart, DATE, SPECS, load, read
from early_pole_native_d8_raw_band_plan import run as plan


def reduce_rows(rows,topology,indices,scalar_edges,scalar_occupations):
    e=topology['edges']
    metadata={k:v for k,v in topology.items() if k!='edges'}
    if ([z.get('panel') for z in rows]!=indices
        or [z.get('bounds') for z in rows]!=[e[i:i+2] for i in indices]
        or any({k:v for k,v in z.items() if k not in ('panel','bounds','coeff')}!=metadata for z in rows)):
        raise ArithmeticError('Need all exact ordered future stock D8 panels, with no central-only or foreign rows')
    bounds=np.array([r['bounds'] for r in rows])
    weights=np.array([r['coeff'] for r in rows])
    if (weights.shape!=(len(indices),3,3) or not np.isfinite(weights).all() or np.any(weights<=0.)
        or np.max(abs(weights[:,:,0]-weights[:,:,1]-weights[:,:,2])/weights[:,:,0])>5e-13):
        raise ArithmeticError('Positive finite three-by-three weights with independent affine closure required')
    bin_ids=address(scalar_edges,bounds)
    f=np.asarray(scalar_occupations)[bin_ids]
    raw=np.stack([strip(weights[:,j],f) for j in range(3)],axis=1)
    total=np.array([[fsum(raw[:,i,j]) for j in range(3)] for i in range(3)])
    if np.max(abs((metadata['scalar_mass']*total[0]+total[2])/total[1]-1.))>2e-13:
        raise ArithmeticError('Unweighted centered energy must reproduce the separate energy moment')
    return total


def run(root,candidate,raw):
    record=candidate.read_bytes()
    A=json.loads(record)
    T=root/f'early-pole-native-cold-weight-target-{DATE}.json'
    if plan(root,T)!=A:
        raise ArithmeticError('Require an exact byte-tracked future raw band plan')
    a,h=read(root/f'early-pole-cell-endpoint-fine-d8-x05-state-{DATE}.json.gz')
    if h!=A['d8_state_sha256_uncompressed']:
        raise ArithmeticError('The candidate must refer to the same exact native state')
    top=[native_topology(a,m,None,topology_only=True) for m in (1,2)]
    d=CellChart(load('fine')[0],SPECS['fine'],8,a['reference_x'])
    n=np.asarray(a['endpoint_occupation'][:d.ns])
    i=A['planned']['panel_indices']
    t,u=[reduce_rows(ledger(f),s,i,d.ds.ed,n) for f,s in zip(raw,top,strict=True)]
    delta=float(np.max(abs(t/u-1.)))
    if delta>6e-5:
        raise ArithmeticError('The two complete raw selected meshes are not resolved')
    return dict(q=A['planned']['span'],native_stock_indices=i,raw_coefficients_by_mesh=[t.tolist(),u.tolist()],
                raw_moments=['number','energy','centered_energy'],affine_columns=['forward','inverse','slope'],
                maximum_observed_raw_relative_difference=delta,plan_sha256=hashlib.sha256(record).hexdigest(),
                raw_ledger_uncompressed_sha256=[unpacked_hash(f) for f in raw],
                scope='Exact selected upper-band raw D8 only. Not central weighted gate, full-band, physical occupation, or spatial error estimate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('candidate',type=Path)
    p.add_argument('raw',nargs=2,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
