"""Compare raw and Bose-weighted momentum targets on exact accepted ledgers."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import compare
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import address, ledger, reweight, run as verify, strip
from early_pole_native_cold_cross_factors import CellChart, DATE, SPECS, load, read


def mass_interval(bounds,mass,a,b):
    l,r=bounds.T
    t=fsum(mass)
    return [fsum(mass[(l>=a)&(r<=b)])/t,fsum(mass[(r>a)&(l<b)])/t]


def median_interval(bounds,mass):
    i=int(np.searchsorted(np.cumsum(mass/fsum(mass)),.5))
    return bounds[i].tolist()


def summary(rows):
    z=[]
    for bounds,raw,weighted in rows:
        U,V=[u/fsum(u) for u in (raw,weighted)]
        z.append(dict(raw_median_possible=median_interval(bounds,raw),
                      weighted_median_possible=median_interval(bounds,weighted),
                      exact_raw_vs_weighted_TV=fsum(abs(U-V))/2.,
                      central_raw_fraction=mass_interval(bounds,raw,3e-6,5e-6),
                      central_weighted_fraction=mass_interval(bounds,weighted,3e-6,5e-6),
                      upper_raw_fraction=mass_interval(bounds,raw,1e-5,1.5e-5),
                      upper_weighted_fraction=mass_interval(bounds,weighted,1e-5,1.5e-5)))
    e=[np.r_[b[:,0],b[-1,1]] for b,*_ in rows]
    metrics={}
    for t,i in (('raw',1),('weighted',2)):
        r=compare(e[0],e[1],rows[0][i],rows[1][i])
        metrics[t]=dict(cdf_cut_with_largest_guaranteed_difference=r['strongest_q'],
                        cdf_lower_bound=r['normalized_TV_lower_bound_via_cdf'],
                        minimum_compatible_TV=r['sharp_compatible_TV']['minimum_TV_bounds'])
    return dict(charts=z,chart_distance=metrics)


def run(root,audit):
    raw=audit.read_bytes()
    if verify(root)!=json.loads(raw):
        raise ArithmeticError('Need exact verified Bose reweight receipt')
    path=root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json'
    J=json.loads(path.read_bytes())
    charts=[]
    for i,k in enumerate((2,4),1):
        S,h=read(root/f'early-pole-cell-endpoint-fine-d{k}-x05-state-{DATE}.json.gz')
        if J['cards'][i]['state_sha256_uncompressed']!=h:
            raise ArithmeticError('Exact accepted scalar histogram required')
        C=CellChart(load('fine')[0],SPECS['fine'],k,S['reference_x'])
        f=np.asarray(S['endpoint_occupation'][:C.ns])
        b=[]
        for m in (1,2):
            P=root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz'
            if unpacked_hash(P)!=J['raw_uncompressed_sha256'][i][m-1]:
                raise ArithmeticError('Need both complete accepted chart ledgers')
            r=ledger(P)
            bounds=np.array([x['bounds'] for x in r])
            if not np.array_equal(bounds[1:,0],bounds[:-1,1]):
                raise ArithmeticError('Need contiguous exact native panels')
            n=f[address(C.ds.ed,bounds)]
            W=np.array([x['coeff'][0] for x in r])
            A=strip(W,n)
            b.append((bounds,A[:,0],reweight(A,n)[:,0]))
        charts.append(b)
    out=[summary([z[m] for z in charts]) for m in (0,1)]
    for u,v in zip(out[0]['charts'],out[1]['charts'],strict=True):
        if max(np.max(abs(np.asarray(u[k])-v[k])) for k in u)>2e-10:
            raise ArithmeticError('Observed chart shapes disagree between accepted meshes')
    for k in ('raw','weighted'):
        if max(abs(np.asarray(out[0]['chart_distance'][k]['minimum_compatible_TV'])-
                   out[1]['chart_distance'][k]['minimum_compatible_TV']))>2e-9:
            raise ArithmeticError('Observed inter-chart shape contrasts are not mesh-stable')
    return dict(fine_mesh=out[1],coarse_mesh=out[0],central_q=[3e-6,5e-6],upper_q=[1e-5,1.5e-5],
                finite_bose_audit_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Support-only contrasts on finite archived raw and Bose-weighted forward panels. No physical population, converged spatial shape, or extrapolation to D8.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','audit'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
