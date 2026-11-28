"""A smooth conditional occupation contrast using only native slope supports."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

from early_pole_endpoint_cold_neutral_occupation_corridors import DATE
from early_pole_endpoint_cold_two_zone_occupation import Q_CUT, run as two_zone
from early_pole_endpoint_continuum_cold_band_receipt import read

Q0,Q1=1e-6,1.5e-5


def moments(rows):
    if (not rows or any(not Q0<=r['bounds'][0]<r['bounds'][1]<=Q1 for r in rows)
        or any(r['coeff'][0][2]<=0 for r in rows)
        or any(a['bounds'][1]!=b['bounds'][0] for a,b in zip(rows,rows[1:]))):
        raise ArithmeticError('Positive ordered native slope measure required')
    S=fsum(r['coeff'][0][2] for r in rows)
    B=fsum(r['coeff'][0][1] for r in rows)
    limits=[fsum(r['coeff'][0][2]*((r['bounds'][i]-Q0)/(Q1-Q0)) for r in rows)/S for i in (0,1)]
    return [B/S,*limits,S]


def decisions(d2,d4,margin=1e-15):
    # b>=0; enforce D2 loss using the right endpoint means, D4 gain using left.
    caps=[(r,hi) for row in d2 for r,_,hi,_ in row]
    floors=[(r,lo) for row in d4 for r,lo,_,_ in row]
    # Active impossibility and lower bound come from pairing D2 full and D4 upper.
    crossing=[(u[0]-f[0],u[1]-f[2]) for f,_ in d2 for _,u in d4]
    if any(num<=0. or den<=0. for num,den in crossing) or margin<=0:
        raise ArithmeticError('Required positive neutral and support gaps')
    b0=max(num/den for num,den in crossing)
    limlo=max(r-b0*lo for r,lo in floors)
    limhi=min(r-b0*hi for r,hi in caps)
    # Reject another constraint removing the limiting boundary; both known
    # certificates coincide there up to floating-point evaluation.
    if abs(limlo-limhi)>5e-24 or limlo<=0.:
        raise ArithmeticError('Pairwise lower bound is not the actual affine infimum')
    b=b0+margin
    l=max(r-b*lo for r,lo in floors)
    h=min(r-b*hi for r,hi in caps)
    a=(l+h)/2
    if not 0.<l<a<h:
        raise ArithmeticError('Strict common affine support witness did not resolve')
    z=[[[[S*(a+b*lo-r),S*(a+b*hi-r)] for r,lo,hi,S in pair] for pair in k] for k in (d2,d4)]
    if not (max(i[1] for p in z[0] for i in p)<0.<min(i[0] for p in z[1] for i in p)):
        raise ArithmeticError('Every native support enclosure must resolve all signs')
    return dict(affine_robust_range_infimum=b0,affine_range_witness=b,
                witness_occupation_on_lower_and_upper_endpoint=[a,a+b],
                witness_change_relative_to_low_endpoint=b/a,
                witness_net_intervals_by_chart_mesh_and_band=z)


def run(root,corridor,zone):
    raw=zone.read_bytes()
    Z=two_zone(root,corridor)
    if Z!=json.loads(raw):
        raise ArithmeticError('Exact two-zone evidence must replay first')
    T=[]
    for k in (2,4):
        t=[]
        for m in (1,2):
            a=read(root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz')
            b=[r for r in a if r['bounds'][0]>=Q_CUT]
            if a[0]['bounds'][0]!=Q0 or a[-1]['bounds'][1]!=Q1 or b[0]['bounds'][0]!=Q_CUT:
                raise ArithmeticError('Exact same full and exact upper panel edges required')
            t.append([moments(v) for v in (a,b)])
        T.append(t)
    ans=decisions(*T)
    return dict(**ans,slope_weighted_normalized_momentum_enclosures=T,
                retained_fields=['neutral','weighted_q_unit_lower','weighted_q_unit_upper','slope'],
                q_span=[Q0,Q1],upper_span=[Q_CUT,Q1],
                range_ratio_to_optimal_two_zone=ans['affine_robust_range_infimum']/Z['two_zone_infimum_of_upper_minus_lower'],
                two_zone_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Affine held occupation sufficient for all within-panel positive slope shapes consistent with retained masses. No physical occupation, interpolation, or spatial convergence asserted.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','corridor','zone'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
