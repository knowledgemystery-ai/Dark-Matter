"""Find the minimum local rotation from the instantaneous fixed-total plane."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_transported_entropy_plane import forms
from sm_dark_entropy_exchange_reference import StandardModel

DATE='2026-11-28'


def cone_distance(K,n,grid=262144):
    """Bound the closest entropy-whitened nonnegative-rate direction to a plane."""
    v,Q=eigh(K)
    if not v[0]<=v[1]<0.<v[2] or not 8<=grid:
        raise ArithmeticError('Exactly one expanding direction and two strictly negative directions required')
    n=n/np.linalg.norm(n)
    P=np.linalg.svd(n[None,:],full_matrices=True)[2][1:].T
    if not eigh(P.T@K@P,eigvals_only=True)[-1]<0.:
        raise ArithmeticError('Reference plane already intersects the nonnegative cone')
    radius=np.sqrt(v[-1]/(-v[:2]))
    t=np.arange(grid)*(2*np.pi/grid)
    a,b=radius[0]*np.cos(t),radius[1]*np.sin(t)
    vals=abs(((n@Q)@np.array([a,b,np.ones_like(t)]))/np.sqrt(1+a*a+b*b))
    g=float(min(vals))
    eps=float(max(radius)*np.pi/grid+3e-14)
    encl=np.rad2deg(np.arcsin([max(g-eps,0.),min(g+3e-14,1.)]))
    return dict(spectral_rates=v.tolist(),degrees_interval=encl.tolist(),sample_minimal_abs_normal_dot=g,
                sampling_slack_in_dot=eps,angular_samples=grid)


def plane_angle(P,normal):
    U=np.linalg.qr(P,mode='reduced')[0]
    return float(np.rad2deg(np.arcsin(min(np.linalg.norm(U.T@normal),1.))))


def run(eos,root):
    names=('direct-expansion-055-096','direct-expansion-096-120','local-entropy-metric-rates','transported-entropy-plane')
    D=[(root/f'fd-portal-{n}-{DATE}.json').read_bytes() for n in names]
    A,B,L,T=map(json.loads,D)
    ids=[hashlib.sha256(t).hexdigest() for t in D]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if any(v['table_input_sha256']!=h for v in (A,B,L,T)) or T['input_sha256'][:2]!=ids[:2] or T['input_sha256'][-1]!=ids[2]:
        raise ArithmeticError('Require same EOS, exact modal and transported-plane witnesses')
    u,v=[next(d['points'] for d in z['rows'] if d['operator']=='combined') for z in (A,B)]
    if u[-1]!=v[0]:
        raise ArithmeticError('Exact held bridge required')
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    xcs=[]
    for x in (.61,.99):
        p=next(r for r in u+v[1:] if r[0]==x)
        l,t=[next(r for r in z['rows'] if r['x']==x) for z in (L,T)]
        G,S=forms(P,p)
        vl,W=eigh(G)
        W_half=(W*np.sqrt(vl))@W.T
        W_inv=(W/np.sqrt(vl))@W.T
        K=W_inv@S@W_inv
        N=W_inv[:,2]/np.linalg.norm(W_inv[:,2])
        cols=np.asarray(t['endpoint_fractional_plane_columns'])
        U=np.linalg.qr(W_half@cols,mode='reduced')[0]
        c=cone_distance(K,N)
        rho=eigh(U.T@K@U,eigvals_only=True)
        if (max(abs(np.array(c['spectral_rates'])-l['signed_squared_norm_rate']))>5e-10
            or max(abs(rho-t['transport_fixed_initial_total_plane_squared_norm_rate']))>5e-10):
            raise ArithmeticError('Whitened ambient and transported spectra must independently close')
        xcs.append(dict(x=x,threshold_cone_distance=c,
                        carried_vs_current_plane_entropy_metric_angle_degrees=plane_angle(W_half@cols,N),
                        carried_vs_current_plane_coordinate_angle_degrees=plane_angle(cols,np.eye(3)[:,2]),
                        actual_carried_plane_rates=rho.tolist()))
    return dict(rows=xcs,input_sha256=ids,table_input_sha256=h,
                scope='Local geometry for two held exact numerical FD rate forms; no time interpolation, physical stability, or coordinate-plane equivalence inferred.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
