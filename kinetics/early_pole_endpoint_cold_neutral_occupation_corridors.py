"""Held-occupation sign corridors from accepted D2/D4 cold panel ledgers."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_receipt import read
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_endpoint_fine_d4_cold_subband_enclosure import box, run as subband

DATE='2026-11-28'


def strictly(a,b):
    return [float(a),float(b)] if a<b else None


def root_bounds(bounds, coeff, lo, hi):
    inverse,slope=[box(bounds,c,lo,hi) for c in coeff.T[1:]]
    if min(inverse[0],slope[0])<=0.:
        raise ArithmeticError('Positive enclosed inverse and slope are required')
    return [inverse[0]/slope[1],inverse[1]/slope[0]]


def corridor(d2,d4):
    """Open values n that make D2 lose and D4 gain, or the converse."""
    return dict(d2_loss_d4_gain_sufficient=strictly(d4[1],d2[0]),
                d2_loss_d4_gain_possible_within_enclosures=strictly(d4[0],d2[1]),
                d2_gain_d4_loss_sufficient=strictly(d2[1],d4[0]),
                d2_gain_d4_loss_possible_within_enclosures=strictly(d2[0],d4[1]))


def hull(rows):
    return [min(z[0] for z in rows),max(z[1] for z in rows)]


def overlap(a,b):
    return strictly(max(a[0],b[0]),min(a[1],b[1])) if a is not None and b is not None else None


def sign(n,roots):
    return int(n>roots[1])-int(n<roots[0])


def run(root):
    L=root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json'
    B=root/f'early-pole-fine-d4-cold-subband-enclosure-{DATE}.json'
    external=json.loads(B.read_bytes())
    ladder=json.loads(L.read_bytes())
    cards=[[root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz' for m in (1,2)] for k in (2,4)]
    if ([[unpacked_hash(z) for z in meshes] for meshes in cards]!=ladder['raw_uncompressed_sha256'][1:]
        or subband(L,[z[1] for z in cards])!=external):
        raise ArithmeticError('Both meshes and the D4 subband receipt must match the accepted ladder')
    tables=[]
    for meshes in cards:
        pair=[]
        for f in meshes:
            v=read(f)
            g=np.array([r['coeff'][0] for r in v])
            if not np.isfinite(g).all() or np.any(g<=0.):
                raise ArithmeticError('Need strictly positive and complete gross panel weights')
            pair.append((np.array([u['bounds'] for u in v]),g))
        tables.append(pair)
    intervals=[external['retained_span'],*[u['q'] for u in external['rows']]]
    n=external['conditional_photon_occupation']
    out=[]
    for lo,hi in intervals:
        v=[[root_bounds(x,y,lo,hi) for x,y in pair] for pair in tables]
        R,S=map(hull,v)
        out.append(dict(q=[lo,hi],neutral_D2_each_mesh=v[0],neutral_D4_each_mesh=v[1],
                        neutral_D2_union_hull=R,neutral_D4_union_hull=S,
                        original_diagnostic_held_occupation_signs=[sign(n,R),sign(n,S)],
                        **corridor(R,S)))
    F=out[0]
    U=out[-1]
    center=next(z for z in out if z['q']==[3e-6,5e-6])
    suf='d2_loss_d4_gain_sufficient'
    poss='d2_loss_d4_gain_possible_within_enclosures'
    upper_low=min(U['neutral_D2_union_hull'][0],U['neutral_D4_union_hull'][0])
    full_hi=F[poss][1]
    return dict(rows=out,old_diagnostic_occupation=n,
                simultaneous_full_and_central_sufficient=overlap(F[suf],center[suf]),
                simultaneous_full_and_upper_possible=overlap(F[poss],U[poss]),
                both_upper_negative_throughout_full_opposite_sign_corridor=full_hi<=upper_low,
                gap_from_full_corridor_ceiling_to_first_upper_neutral=upper_low-full_hi,
                comparison_meshes=[1,2],ladder_sha256=hashlib.sha256(L.read_bytes()).hexdigest(),
                subband_sha256=hashlib.sha256(B.read_bytes()).hexdigest(),
                scope='Bounds only from exact admitted positive panel weights and two sampled quadratures. A constant held occupation is a diagnostic, not a live photon solution or spatial limit.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
