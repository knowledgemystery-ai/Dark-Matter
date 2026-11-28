"""Bound the mean q using archived centered-energy moments, not interpolation."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import address, ledger, reweight, strip
from early_pole_native_cold_cross_factors import CellChart, DATE, SPECS, load, read
from early_pole_native_cold_shape_attribution import moment_interval, run as prior_shape


def conditional_q(bounds,n,center,m):
    u=bounds*bounds/(np.hypot(bounds,m)+m)
    t=center/n
    scale=np.max(u,axis=1)*3e-13
    if (m<=0. or min(n)<=0. or not np.isfinite(t).all() or np.any(np.diff(u,axis=1)<=0.)
        or np.any(t<u[:,0]-scale) or np.any(t>u[:,1]+scale)):
        raise ArithmeticError('Centered energy must lie inside each exact panel')
    # q(u)=sqrt(u(2m+u)) is strictly concave. Use the chord for a
    # conservative lower bound, Jensen for an upper one, and widen for roundoff.
    tl=np.maximum(t-scale,u[:,0])
    th=np.minimum(t+scale,u[:,1])
    a=bounds[:,0]+(bounds[:,1]-bounds[:,0])*(tl-u[:,0])/(u[:,1]-u[:,0])
    b=np.sqrt(th*(2*m+th))
    return np.c_[np.maximum(bounds[:,0],a-2e-13*bounds[:,1]),
                 np.minimum(bounds[:,1],b+2e-13*bounds[:,1])]


def compare(x):
    result={}
    for label,j in (('forward',0),('inverse',1)):
        raw=[]; weighted=[]; shape=[]; kernel=[]; direct=[]
        for i,(q,Y,nh) in enumerate(x):
            m=Y[0]
            A,B=Y[1]
            u=conditional_q(q,A[:,j],B[:,j],m)
            w=A[:,j]/fsum(A[:,j])
            C,D=[reweight(A,h)[:,j] for h in nh]
            C,D=[z/fsum(z) for z in (C,D)]
            sign=(-1.)**(1-i)
            raw.append(moment_interval(u,w))
            weighted.append(moment_interval(u,(C,D)[i]))
            shape.append(moment_interval(u,(D-C)/2))
            kernel.append(moment_interval(u,sign*(C+D)/2))
            direct.append(moment_interval(u,sign*(C,D)[i]))
        add=lambda z: np.sum(z,axis=0).tolist()
        result[label]=dict(raw_mean_q_by_chart=np.asarray(raw).tolist(),native_weighted_mean_q_by_chart=np.asarray(weighted).tolist(),
                           raw_native_mean_difference=(raw[1]-raw[0][::-1]).tolist(),
                           weighted_native_mean_difference=add(direct),
                           occupation_only_mean_difference=add(shape),kernel_only_mean_difference=add(kernel))
    return result


def run(root,shape):
    raw=shape.read_bytes()
    T=root/f'early-pole-native-cold-bose-reweight-{DATE}.json'
    if prior_shape(root,T)!=json.loads(raw):
        raise ArithmeticError('Require exact stored normalized shape comparison')
    J=json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    parts=[]; hist=[]
    for k,i in ((2,1),(4,2)):
        a,h=read(root/f'early-pole-cell-endpoint-fine-d{k}-x05-state-{DATE}.json.gz')
        if h!=J['cards'][i]['state_sha256_uncompressed']:
            raise ArithmeticError('Exact previously admitted state required')
        d=CellChart(load('fine')[0],SPECS['fine'],k,a['reference_x'])
        hist.append(np.asarray(a['endpoint_occupation'][:d.ns]))
        out=[]
        for m in (1,2):
            L=root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz'
            if unpacked_hash(L)!=J['raw_uncompressed_sha256'][i][m-1]:
                raise ArithmeticError('Both complete exact meshes required')
            R=ledger(L); q=np.array([r['bounds'] for r in R])
            w=np.array([r['coeff'] for r in R]); n=R[0]['scalar_mass']
            if np.max(abs((n*w[:,0,:]+w[:,2,:])/w[:,1,:]-1.))>2e-13:
                raise ArithmeticError('Centered energy and independent energy must close')
            idx=address(d.ds.ed,q)
            out.append((q,n,w,idx))
        parts.append(out)
    result=[]
    for m in (0,1):
        X=[]
        for i,r in enumerate(parts):
            q,n,w,idx=r[m]
            weights=[strip(w[:,j],hist[i][idx]) for j in (0,2)]
            X.append((q,(n,weights),[y[idx] for y in hist]))
        result.append(compare(X))
    err=max(np.max(abs(np.asarray(result[0][k][j])-result[1][k][j])) for k in result[0] for j in result[0][k])
    if err>2e-13:
        raise ArithmeticError('Independent meshes must resolve the moment enclosures')
    return dict(by_mesh=result,maximum_mesh_difference=err,exact_shape_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Concavity bounds from complete held panel masses and centered energy, not an assumed intrapanel shape or physical convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','shape'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
