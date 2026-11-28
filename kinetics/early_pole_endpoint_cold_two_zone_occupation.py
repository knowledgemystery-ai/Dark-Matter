"""Can one nonconstant held spectrum reverse signs on two nested bands?"""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

from early_pole_endpoint_cold_neutral_occupation_corridors import DATE, run as corridor
from early_pole_endpoint_continuum_cold_band_receipt import read

Q_CUT=1e-5


def summarize(rows):
    lo=[r for r in rows if r['bounds'][1]<=Q_CUT]
    hi=[r for r in rows if r['bounds'][0]>=Q_CUT]
    if not lo or not hi or lo+hi!=rows or lo[-1]['bounds'][1]!=hi[0]['bounds'][0]:
        raise ArithmeticError('Exact common native edge required, without an unknown split panel')
    return [[fsum(r['coeff'][0][k] for r in zone) for k in (1,2)] for zone in (lo,hi)]


def allowed_low_band(pair, u):
    return [(fsum([B0,B1,-D1*u]))/D0 for ((B0,D0),(B1,D1)) in pair]


def net(p,lo,hi):
    (B0,D0),(B1,D1)=p
    return [fsum([D0*lo,-B0,D1*hi,-B1]),fsum([D1*hi,-B1])]


def infimum_and_witness(d2,d4,slack=1e-16):
    if slack<=0. or any(v<=0 for P in (d2,d4) for r in P for W in r for v in W):
        raise ArithmeticError('Positive diagnostic slopes, inverse weights and margin required')
    u0=max(B/D for _,(B,D) in d4)
    lmax=min(allowed_low_band(d2,u0))
    lmin=max(allowed_low_band(d4,u0),default=0.)
    ucap=min(B/D for _,(B,D) in d2)
    # For u>=u0, the D2 upper constraint on l decreases strictly with u.
    # If the other constraints have positive slack at the limiting boundary,
    # the resulting lower bound on u-l is the exact two-zone infimum.
    if not (0.<lmax<u0<ucap and max(lmin,0.)<lmax):
        raise ArithmeticError('These ledgers do not support the prespecified open-limit construction')
    hi=u0+slack
    lo=min(allowed_low_band(d2,hi))-slack
    a,b=[[net(r,lo,hi) for r in pair] for pair in (d2,d4)]
    if not (0.<lo<hi and max(v for r in a for v in r)<0.<min(v for r in b for v in r)):
        raise ArithmeticError('Strict four-way sign witness failed')
    return dict(two_zone_infimum_of_upper_minus_lower=u0-lmax,two_zone_limit=[lmax,u0],
                strict_constructive_witness=[lo,hi],witness_range=hi-lo,
                witness_d2_net_full_and_upper_by_mesh=a,witness_d4_net_full_and_upper_by_mesh=b,
                other_constraints_limit_room=[lmax-max(lmin,0.),ucap-u0])


def run(root,receipt):
    R=corridor(root)
    raw=receipt.read_bytes()
    if R!=json.loads(raw):
        raise ArithmeticError('Exact original neutral corridors must replay')
    V=[]
    for k in (2,4):
        r=[read(root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz') for m in (1,2)]
        V.append([summarize(v) for v in r])
    a=infimum_and_witness(*V)
    # Any common shape must reach above every accepted D4 upper neutral
    # threshold and below every accepted D2 full neutral threshold.
    fn=R['rows'][0]['neutral_D2_union_hull'][0]
    un=R['rows'][-1]['neutral_D4_union_hull'][1]
    bound=un-fn
    if not 0.<bound<=a['two_zone_infimum_of_upper_minus_lower']:
        raise ArithmeticError('Construction must be consistent with the arbitrary-profile lower bound')
    return dict(q_cut=Q_CUT,arbitrary_common_profile_range_strict_lower_bound=bound,
                two_zone_infimum_over_arbitrary_profile_lower_bound=a['two_zone_infimum_of_upper_minus_lower']/bound,
                **a,retained_fields=['inverse','slope'],zones_by_mesh=dict(D2=V[0],D4=V[1]),
                neutral_corridor_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Common piecewise held photon diagnostic on two complete fixed-state meshes per chart. Neither a live transport solution nor a spatial convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('receipt',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.receipt),indent=2,allow_nan=False))
