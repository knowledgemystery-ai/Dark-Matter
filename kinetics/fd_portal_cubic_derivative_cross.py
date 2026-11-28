"""Cross EOS heat-capacity log slopes while retaining each base local form."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import full
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_entropy_knot_jump_audit import forms as knot_forms
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_heat_capacity_slope_threshold import DATE, threshold
from fd_portal_local_entropy_metric_rates import motion
from fd_portal_transported_entropy_plane import forms
from fd_portal_visible_curvature_rate_attribution import visible
from sm_dark_entropy_exchange_reference import StandardModel


def slope_and_rank(stage,p):
    g=stage.eos
    t=np.log(p[1])
    a,b,c=g.gs(t),g.slope(t),g.slope.derivative()(t)
    k=3.+(3.*b+c)/(3.*a+b)
    rates=full(stage,p[1:4])(np.ones(3)).sum(axis=0)
    v,_=motion(stage,p,rates)
    return k,-v[1]/p[1]*visible(stage,p)


def two_by_two(forms_and_slopes):
    out=[]
    for G,S,R,k in forms_and_slopes:
        t=threshold(G,S,R,k)
        rates=[eigh(S+(u-k)*R,G,eigvals_only=True)[-1] for _,_,_,u in forms_and_slopes]
        out.append(dict(critical_heat_capacity_slope=t['critical_slope_above_three'],
                        native_heat_capacity_slope=k,rates_with_pchip_and_cubic_slope=rates))
    return out


def run(eos,root):
    names=('direct-expansion-055-096','direct-expansion-096-120','local-entropy-metric-rates',
           'visible-curvature-rate-attribution','entropy-eos-interpolation-check',
           'entropy-pchip-knot-crossing','heat-capacity-slope-threshold')
    D=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    E,L,W,V,C,K,T=map(json.loads,D)
    ids=[hashlib.sha256(s).hexdigest() for s in D]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if (C['input_sha256']!=ids[:3] or V['input_sha256']!=ids[:3]
        or K['input_sha256'][0]!=ids[1] or T['input_sha256'][:4]!=ids[:4]
        or T['input_sha256'][-1]!=ids[5] or any(d['table_input_sha256']!=h for d in (E,L,W,V,C,K,T))):
        raise ArithmeticError('Exact archived forms, spline, knot and threshold receipts required')
    u,v=[next(z['points'] for z in a['rows'] if z['operator']=='combined') for a in (E,L)]
    if u[-1]!=v[0]:
        raise ArithmeticError('Expected exact archived bridge')
    P,Q=[CombinedExpansion(f(eos),PositivePortal((17,44,32)),168) for f in (StandardModel,SplineEOS)]
    walk=local_path(P,next(p for p in v if p[0]==.99))
    rows=[]
    def evaluate(p,method,prior,critical):
        packets=[]
        for stage in (P,Q):
            G,S=method(stage,p)
            k,R=slope_and_rank(stage,p)
            packets.append((G,S,R,k))
        z=two_by_two(packets)
        actual=[a['rates_with_pchip_and_cubic_slope'][i] for i,a in enumerate(z)]
        if (max(abs(np.array(actual)-prior))>5e-8 or abs(z[0]['critical_heat_capacity_slope']-critical)>5e-9):
            raise ArithmeticError('Both native columns must reproduce independent receipts')
        pred=z[0]['rates_with_pchip_and_cubic_slope'][1]
        rows.append(dict(x=p[0],T=p[1],base_pchip=z[0],base_cubic=z[1],
                         full_cubic_minus_pchip_with_only_cubic_derivative=actual[1]-pred,
                         opposite_sign_from_swapping_only_derivative=bool(pred*actual[1]<0.)))
    for x in (.61,.99):
        p=next(s for s in u+v[1:] if s[0]==x)
        old=next(r for r in C['rows'] if r['x']==x)
        new=next(r for r in T['rows'] if r['x']==x)
        evaluate(p,forms,[old['original_max_rate'],old['alternate_max_rate']],new['spaces']['full']['critical_slope_above_three'])
    for probe in K['rows']:
        for old in probe['before_and_after']:
            p=walk(old['x'])
            r=next(t for t in T['knot_full_space_probes'] if t['x']==old['x'])
            evaluate(p,knot_forms,[old['eos_results'][k]['max_rate'] for k in ('pchip','cubic')],r['critical_slope_above_three'])
    return dict(rows=rows,input_sha256=ids,table_input_sha256=h,
                scope='At held background T, eta and dark T; swap derivatives only inside each EOS local form. No counterfactual EOS trajectory, physical attribution or stability theorem.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
