"""Algebraically remove and cross the finite scalar occupation on exact panels."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_receipt import read as ledger
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_cross_factors import CellChart, DATE, SPECS, load, read


def reweight(raw, f):
    A,B,D=np.asarray(raw).T
    return np.array([(1+f)*A,f*B,A+f*D]).T


def strip(weighted,f):
    if min(f)<=0.:
        raise ArithmeticError('Retained scalar occupation must be positive')
    A,B,D=np.asarray(weighted).T
    A=A/(1+f)
    z=np.array([A,B/f,(D-A)/f]).T
    if min(z.ravel())<=0. or max(abs(reweight(z,f)/weighted-1.).ravel())>2e-15:
        raise ArithmeticError('Original finite factors must reconstruct exactly on the same panel')
    return z


def address(edges,bounds):
    edges,bounds=map(np.asarray,(edges,bounds))
    if (edges.ndim!=1 or len(edges)<2 or np.any(np.diff(edges)<=0.)
        or bounds.ndim!=2 or len(bounds)==0 or bounds.shape[1]!=2
        or np.any(bounds[:,0]>=bounds[:,1]) or not np.isfinite(np.r_[edges,bounds.ravel()]).all()):
        raise ArithmeticError('Finite ordered scalar grid and panel bounds required')
    i=np.searchsorted(edges,bounds[:,0],side='right')-1
    if min(i)<0 or max(i)>=len(edges)-1 or np.any(bounds[:,1]>edges[i+1]):
        raise ArithmeticError('Every complete panel must lie within one exact occupied scalar bin')
    return i


def summarize(inputs,endpoint,n,a,b):
    vals=[]
    for bounds,raw,cell in inputs:
        l,r=bounds.T
        keep=(l>=a)&(r<=b)
        if not any(keep) or np.any(((r>a)&(l<b))!=keep):
            raise ArithmeticError('Use only exact complete retained spans, no intrapanel redistribution')
        def tot(z):
            return np.array([fsum(t) for t in z[keep].T])
        rec=tot(raw)
        pair=[tot(reweight(raw,z[cell])) for z in endpoint]
        vals.append((rec,pair))
    r,p=map(np.asarray,zip(*vals,strict=True))
    A,B,C,D=p[0,0],p[0,1],p[1,0],p[1,1]
    dis=np.log(D/A)
    occ=(np.log(B/A)+np.log(D/C))/2
    kernel=(np.log(C/A)+np.log(D/B))/2
    if max(abs(occ+kernel-dis))>8e-15:
        raise ArithmeticError('Crossed exact occupations and cold panel kernels must close')
    return dict(q=[a,b],unweighted_coefficients=r.tolist(),
                weighted_coefficients_by_raw_chart_then_scalar_donor=p.tolist(),
                unweighted_d4_over_d2=(r[1]/r[0]).tolist(),weighted_native_d4_over_d2=(D/A).tolist(),
                raw_chart_change_at_either_fixed_scalar_histogram=[(C/A).tolist(),(D/B).tolist()],
                all_orders_raw_chart_and_explicit_occupation_factors=np.exp([kernel,occ]).tolist(),
                conditional_net_at_ladder_held_occupation=(n*p[:,:,2]-p[:,:,1]).tolist())


def run(root):
    L=root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json'
    data=L.read_bytes()
    J=json.loads(data)
    state=[]; grid=None
    rows=[]
    for ix,lev in enumerate((2,4),1):
        S,h=read(root/f'early-pole-cell-endpoint-fine-d{lev}-x05-state-{DATE}.json.gz')
        if h!=J['cards'][ix]['state_sha256_uncompressed'] or S['subcells_per_archived_cell']!=lev:
            raise ArithmeticError('Need the exact native state for the accepted cold quadrature')
        d=CellChart(load('fine')[0],SPECS['fine'],lev,S['reference_x'])
        if grid is not None and not np.array_equal(grid,d.ds.ed):
            raise ArithmeticError('Scalar occupations must share the exact native scalar grid')
        grid=d.ds.ed
        f=np.array(S['endpoint_occupation'][:d.ns])
        state.append(f)
        mesh=[]
        for i in (1,2):
            p=root/f'fine-d{lev}-endpoint-continuum-cold-{i}.jsonl.gz'
            if unpacked_hash(p)!=J['raw_uncompressed_sha256'][ix][i-1]:
                raise ArithmeticError('Need both exact complete admitted cold ledgers')
            R=ledger(p)
            b=np.array([r['bounds'] for r in R])
            j=address(grid,b)
            W=np.array([r['coeff'][0] for r in R])
            mesh.append((b,strip(W,f[j]),j))
        rows.append(mesh)
    n=next(z['common_diagnostic_occupation'] for z in J['diagnostic_root_midpoints'] if z['levels']==[2,4])
    lo,hi=J['q_band']
    spans=((lo,hi),(lo,1e-5),(1e-5,hi))
    r=[]
    for a,b in spans:
        x,y=[summarize([R[m] for R in rows],state,n,a,b) for m in (0,1)]
        err=max(abs(np.asarray(x[k])/y[k]-1.).ravel().max() for k in ('unweighted_coefficients','weighted_coefficients_by_raw_chart_then_scalar_donor'))
        if err>3e-10:
            raise ArithmeticError('Independent meshes do not resolve reweighting')
        r.append(dict(**y,maximum_observed_two_mesh_relative_difference=err))
    return dict(rows=r,held_photon_occupation=n,ladder_sha256=hashlib.sha256(data).hexdigest(),
                scope='Finite Bose factors removed and crossed exactly on panel supports. The unweighted kernel still includes the native pole, charges, and q²/E measure; no causal trajectory decomposition or physical source inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
