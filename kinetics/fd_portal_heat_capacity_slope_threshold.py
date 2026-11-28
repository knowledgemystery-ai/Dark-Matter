"""Reoptimize FD local entropy rates under a single visible metric derivative change."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh
from scipy.optimize import brentq

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import full
from fd_portal_entropy_knot_jump_audit import forms as knot_forms
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_local_entropy_metric_rates import motion
from fd_portal_transported_entropy_plane import forms
from fd_portal_visible_curvature_rate_attribution import visible
from sm_dark_entropy_exchange_reference import StandardModel

DATE='2026-11-28'


def threshold(G,S,A,kref):
    """Threshold for largest rate zero when replacing only the derivative k."""
    B=S+(3.-kref)*A
    unit=eigh(B,G,eigvals_only=True)
    # If the reference form is already nonnegative, a positive rank-one
    # increment cannot test the intended restoration of strict dissipation.
    if max(unit)>=0.:
        return dict(rate_at_cubic_slope=unit.tolist(),critical_slope_above_three=None,
                    rate_at_archived_slope=eigh(S,G,eigvals_only=True).tolist())
    eps=eigh(A,-B,eigvals_only=True)
    if max(eps)<=1e-14 or min(eps)<-1e-12*max(eps):
        raise ArithmeticError('Expected a resolved positive rank-one perturbation')
    k=3.+1./max(eps)
    def r(v):
        return eigh(B+(v-3.)*A,G,eigvals_only=True)[-1]
    checked=brentq(r,3.,3.+2*(k-3.),xtol=3e-13)
    if abs(checked-k)>5e-10*max(1.,abs(k)) or not r(k-.01)<0.<r(k+.01):
        raise ArithmeticError('Determinant-lemma root must resolve an independent sign crossing')
    return dict(rate_at_cubic_slope=unit.tolist(),critical_slope_above_three=k,
                independent_bisection_slope=checked,rate_at_archive_minus_threshold=kref-k,
                rate_at_archived_slope=eigh(S,G,eigvals_only=True).tolist())


def near_knot(P,K,start):
    walk=local_path(P,start)
    rows=[]
    for probe in K['rows']:
        for r in probe['before_and_after']:
            p=walk(r['x'])
            if abs(p[1]-r['visible_T'])>3e-14:
                raise ArithmeticError('Exact prior short continuation must replay')
            G,S=knot_forms(P,p)
            rates=full(P,p[1:4])(np.ones(3)).sum(axis=0)
            direction,_=motion(P,p,rates)
            u=np.log(p[1])
            e=P.eos
            g,g1,g2=e.gs(u),e.slope(u),e.slope.derivative()(u)
            k=3.+(3.*g1+g2)/(3.*g+g1)
            z=threshold(G,S,(-direction[1]/p[1])*visible(P,p),k)
            if abs(z['rate_at_archived_slope'][-1]-r['eos_results']['pchip']['max_rate'])>3e-8:
                raise ArithmeticError('Spectra must independently reproduce both sides of the knot')
            rows.append(dict(x=r['x'],T=p[1],actual_heat_capacity_log_slope=k,**z))
    return rows


def run(eos,root):
    names=('direct-expansion-055-096','direct-expansion-096-120','local-entropy-metric-rates',
           'visible-curvature-rate-attribution','transported-entropy-plane','entropy-pchip-knot-crossing')
    data=[(root/f'fd-portal-{f}-{DATE}.json').read_bytes() for f in names]
    A,B,L,C,T,K=map(json.loads,data)
    ids=[hashlib.sha256(a).hexdigest() for a in data]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (C['input_sha256']!=ids[:3] or T['input_sha256'][:2]!=ids[:2]
        or T['input_sha256'][-1]!=ids[2] or K['input_sha256'][0]!=ids[1]
        or any(z['table_input_sha256']!=h for z in (A,B,L,C,T,K))):
        raise ArithmeticError('Expected the exact EOS, trajectory, transported plane and local forms')
    u,v=[next(r['points'] for r in z['rows'] if r['operator']=='combined') for z in (A,B)]
    if u[-1]!=v[0]:
        raise ArithmeticError('Expected exact archived bridge')
    E=StandardModel(eos)
    P=CombinedExpansion(E,PositivePortal((17,44,32)),168)
    Q=CombinedExpansion(E,PositivePortal((13,36,26)),112)
    rows=[]
    for x in (.61,.99):
        p=next(r for r in u+v[1:] if r[0]==x)
        l,c,t=[next(r for r in d['rows'] if r['x']==x) for d in (L,C,T)]
        G,S=forms(P,p)
        H=-c['relative_visible_temperature_rate']*visible(P,p)
        g,s=forms(Q,p)
        t0,_=motion(Q,p,full(Q,p[1:4])(np.ones(3)).sum(axis=0))
        h0=(-t0[1]/p[1])*visible(Q,p)
        k=c['dlog_visible_heat_capacity_dlog_temperature']
        U=np.linalg.qr(np.array(t['endpoint_fractional_plane_columns']),mode='reduced')[0]
        d={}
        for label,space,reference in [('full',np.eye(3),l['signed_squared_norm_rate']),
                                     ('fixed_total_now',np.eye(3)[:,:2],l['instantaneous_fixed_total_plane_rate']),
                                     ('carried_fixed_total_initial',U,t['transport_fixed_initial_total_plane_squared_norm_rate'])]:
            z=threshold(*(space.T@r@space for r in (G,S,H)),kref=k)
            if np.max(abs(np.array(z['rate_at_archived_slope'])-reference))>3e-10:
                raise ArithmeticError('Must replay prior rate spectra')
            low=threshold(*(space.T@r@space for r in (g,s,h0)),kref=k)
            delta=abs(z['critical_slope_above_three']-low['critical_slope_above_three'])
            if delta>2e-6:
                raise ArithmeticError('Independent collision quadrature cannot resolve the threshold')
            d[label]=dict(**z,coarser_quadrature_critical_slope=low['critical_slope_above_three'],
                          quadrature_threshold_abs_difference=delta)
        rows.append(dict(x=x,actual_heat_capacity_log_slope=k,relative_visible_temperature_rate=c['relative_visible_temperature_rate'],spaces=d))
    probes=near_knot(P,K,next(w for w in v if w[0]==.99))
    return dict(rows=rows,knot_full_space_probes=probes,input_sha256=ids,table_input_sha256=h,
                scope='Algebraic rank-one derivative substitution at held states and held tangent planes; not an EOS alternative, physical stability criterion, or unobserved time guarantee.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
