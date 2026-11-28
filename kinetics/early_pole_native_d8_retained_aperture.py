"""A topology-only finite aperture check for the native D8 upper plan."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_native_cold_cross_factors import DATE, K, native
from early_pole_native_cold_retained_aperture import M, aperture, p, summarize
from early_pole_native_d8_raw_band_plan import run as plan


def monotone_interval(card,cap,channel,span):
    a,b=span
    A,B=[aperture(q,card,cap,channel) for q in (a,b)]
    x=A[2]+K-a
    y=B[2]+K-b
    j=2 if channel==0 else 1
    m=card['daughter_mass'] if channel==0 else M
    Emin=min(A[3][j],B[3][j])
    vprime=(Emin/p(Emin,m))*b/np.hypot(card['scalar_mass'],b)
    # For N=K^2+q^2-v^2, v decreases; N rises and hence N(a) is a lower bound.
    cprime_guard=b*(2*b+2*x*vprime)-((K-x)*(K+x)+a*a)
    if not x>=y>0. or cprime_guard>=0. or x>=K+a:
        raise ArithmeticError('Cannot certify positive aperture and monotone cone bounds by endpoints')
    return dict(retained_cosine_fraction_range=[B[0],A[0]],unclipped_cosine_limit_endpoints=[A[1],B[1]],
                upper_bound_for_q_times_numerator_derivative_minus_numerator=cprime_guard)


def run(root,planning):
    raw=planning.read_bytes()
    Q=plan(root,root/f'early-pole-native-cold-weight-target-{DATE}.json')
    if Q!=json.loads(raw):
        raise ArithmeticError('Need identical independently replayed D8 raw plan')
    card=native(root,8)
    if card['card_hash']!=Q['d8_state_sha256_uncompressed']:
        raise ArithmeticError('Plan and retained cap must refer to same exact native D8 state')
    rows=[]
    for i,s in enumerate(card['charges']):
        if s.f[-1]<=0.:
            raise ArithmeticError('Tail cell must be occupied within the finite chart')
        cap=float(s.E[-1])
        rows.append(dict(channel=('forward_parent','inverse_daughter')[i],over_full=summarize(card,cap,i),
                         within_selected_raw_plan=monotone_interval(card,cap,i,Q['planned']['span'])))
    cut=rows[1]['over_full']['first_clipped_q']
    if cut is None or cut>=Q['planned']['span'][0]:
        raise ArithmeticError('Entire planned inverse aperture must be after the verified first cut')
    return dict(rows=rows,planned_span=Q['planned']['span'],x=.05,k=K,plan_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Exact finite numerical supports and kinematic cosine enclosures; no collision weights, new D8 kernels, physical exclusion, or change to live central run.')


if __name__=='__main__':
    P=argparse.ArgumentParser(description=__doc__)
    for k in ('root','planning'):
        P.add_argument(k,type=Path)
    print(json.dumps(run(**vars(P.parse_args())),indent=2,allow_nan=False))
