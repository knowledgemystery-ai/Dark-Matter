"""Choose a future raw-focused D8 span without starting another quadrature."""
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
from early_pole_native_cold_weight_target import mass_interval, run as prior_target


def select(edges,tables,fraction):
    if not 0.<fraction<1. or not 2<=len(edges) or np.any(np.diff(edges)<=0.):
        raise ArithmeticError('Ordered D8 edges and a strict coverage target required')
    p,q=np.triu_indices(len(edges),1)
    good=np.full(len(p),True)
    fences=[]
    for S,F in tables:
        if (len(S)!=len(F) or np.any(F<0.) or fsum(F)<=0.
            or not np.array_equal(S[1:,0],S[:-1,1]) or list(S[[0,-1],[0,1]])!=list(np.asarray(edges)[[0,-1]])):
            raise ArithmeticError('Contiguous positive complete raw panel source required')
        r=np.searchsorted(S[:,0],np.asarray(edges)[p],side='left')
        s=np.searchsorted(S[:,1],np.asarray(edges)[q],side='right')
        mass=np.r_[0.,np.cumsum(F/fsum(F))]
        gap=mass[s]-mass[r]
        fences.append(np.min(abs(gap-fraction)))
        good &= gap>=fraction
    if min(fences)<2e-13 or not np.any(good) or not np.isfinite(fences).all():
        raise ArithmeticError('Enumerated coverage ordering unresolved near the numerical floor')
    i,j=min(zip(p[good],q[good],strict=True),key=lambda x:(x[1]-x[0],edges[x[1]]-edges[x[0]],x[0]))
    band=[edges[i],edges[j]]
    result=[mass_interval(b,m,*band) for b,m in tables]
    if min(v[0] for v in result)<fraction:
        raise ArithmeticError('Direct positive-sum replay failed selected coverage')
    return dict(span=band,panel_indices=list(range(int(i),int(j))),candidate_interval_count=len(p),
                intervals_meeting_coverage=int(np.count_nonzero(good)),coverage=result,
                minimum_candidate_gap_to_target_observed=min(fences))


def run(root,previous):
    raw=previous.read_bytes()
    T=json.loads(raw)
    if prior_target(root,root/f'early-pole-native-cold-bose-reweight-{DATE}.json')!=T:
        raise ArithmeticError('Expected exact previous raw/weighted gate')
    J=json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    a,g=read(root/f'early-pole-cell-endpoint-fine-d8-x05-state-{DATE}.json.gz')
    if (a['grid'],a['subcells_per_archived_cell'],a['steps'],a['x_end'])!=('fine',8,2048,.05):
        raise ArithmeticError('Exact independent D8 endpoint required')
    u,v=[native_topology(a,m,None,topology_only=True) for m in (1,2)]
    edges=u['edges']
    if edges!=v['edges'] or [edges[0],edges[-1]]!=J['q_band'] or len(edges)-1!=u['total'] or u['total']!=v['total']:
        raise ArithmeticError('Both planned D8 topologies must share the full exact retained span')
    tables=[]
    labels=[]
    for k,i in ((2,1),(4,2)):
        s,h=read(root/f'early-pole-cell-endpoint-fine-d{k}-x05-state-{DATE}.json.gz')
        if h!=J['cards'][i]['state_sha256_uncompressed']:
            raise ArithmeticError('Exact scalar occupations from accepted D2/D4 required')
        S=CellChart(load('fine')[0],SPECS['fine'],k,s['reference_x'])
        f=np.asarray(s['endpoint_occupation'][:S.ns])
        for m in (1,2):
            p=root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz'
            if unpacked_hash(p)!=J['raw_uncompressed_sha256'][i][m-1]:
                raise ArithmeticError('Exact complete archived occupations and ledgers required')
            D=ledger(p)
            B=np.array([r['bounds'] for r in D])
            Y=strip(np.array([r['coeff'][0] for r in D]),f[address(S.ds.ed,B)])
            for chan,j in (('forward',0),('inverse',1)):
                tables.append((B,Y[:,j]))
                labels.append([k,m,chan])
    z=select(edges,tables,.8)
    d8=CellChart(load('fine')[0],SPECS['fine'],8,a['reference_x'])
    e=np.asarray(edges)
    bins=address(d8.ds.ed,np.c_[e[:-1],e[1:]])
    return dict(planned=z,coverage_labels=labels,full_d8_panel_count=u['total'],
                planned_scalar_bins=sorted({int(bins[i]) for i in z['panel_indices']}),
                d8_state_sha256_uncompressed=g,topology_sha256=hashlib.sha256(json.dumps(edges,separators=(',',':')).encode()).hexdigest(),
                previous_target_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Topology-only future raw-span plan based on complete D2/D4 coverage; no new D8 kernels, no predicted D8 coverage, and no modification of the separate live central gate.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','previous'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
