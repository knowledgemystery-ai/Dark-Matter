"""Use the exact D4 cuts to shrink neighboring held-occupation corridors."""
import argparse
import hashlib
import json
from pathlib import Path
from shutil import copy2
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_endpoint_cold_neutral_occupation_corridors import hull, sign, strictly
from early_pole_endpoint_fine_d4_cold_continuum_ladder import DATE, unpacked_hash
from early_pole_native_cold_centered_energy_corridors import roots
from early_pole_native_cold_subband_from_centered_energy import interval, ledger
from early_pole_native_d4_clipped_next_band import finish

EXPECTED_PATCH = (True, False, False, True, True, True)
EXPECTED_SIGNS = ([-1,1],[0,0],[-1,1],[-1,1],[-1,-1],[-1,-1])


def run(root):
    p = root / f'early-pole-native-d4-clipped-next-band-{DATE}.json'
    q = root / f'early-pole-native-cold-centered-energy-corridors-{DATE}.json'
    C,L = [json.loads(f.read_bytes()) for f in (p,q)]
    X = json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    n = C['diagnostic_occupation']
    if n != L['diagnostic_occupation'] or len(L['rows']) != len(EXPECTED_PATCH):
        raise ArithmeticError('Preserve the admitted held diagnostic and partition')
    D2 = []
    for i in (1,2):
        f = root/f'fine-d2-endpoint-continuum-cold-{i}.jsonl.gz'
        if unpacked_hash(f) != X['raw_uncompressed_sha256'][1][i-1]:
            raise ArithmeticError('Exact accepted native D2 full ledgers required')
        D2.append(ledger(f))
    with TemporaryDirectory() as d:
        W = Path(d)
        for i in (1,2):
            f = root/f'd4-clipped-next-band-patches-{DATE}'/f'd4-cold-split-{i}.jsonl'
            copy2(f,W)
        if finish(root,W) != C:
            raise ArithmeticError('D4 exact-cut acceptance must replay from the original and raw patch')
        D4 = [ledger(W/f'd4-cold-clipped-joined-{i}.jsonl') for i in (1,2)]
    result = []
    for old, want, signs in zip(L['rows'],EXPECTED_PATCH,EXPECTED_SIGNS,strict=True):
        l,r = old['q']
        U = [interval(v,l,r) for v in D2]
        if hull([roots(e) for e in U]) != old['charts'][0]['refined_neutral_hull']:
            raise ArithmeticError('Original D2 centered-moment comparator must replay')
        V = [interval(v,l,r) for v in D4]
        full = all(l in [z['bounds'][0] for z in v] and r in [z['bounds'][1] for z in v] for v in D4)
        if full != want:
            raise ArithmeticError('Exact D4 face geometry does not match the split')
        # Each component keeps its own observed extrema. This does not
        # extrapolate an angular error beyond the admitted two meshes.
        a,b = [np.stack([f[:,:,0].min(axis=0),f[:,:,1].max(axis=0)],axis=-1) for f in map(np.stack,(U,V))]
        RR = [roots(f) for f in (a,b)]
        found = [sign(n,f) for f in RR]
        if found!=signs:
            raise ArithmeticError('Expected same held diagnostic and resolved sign pattern')
        old4 = old['charts'][1]['refined_neutral_hull']
        if RR[1][0] < old4[0]*(1-4e-12) or RR[1][1] > old4[1]*(1+4e-12):
            raise ArithmeticError('Clipped information should not widen the existing D4 root bounds')
        result.append(dict(q=old['q'],neutral_D2_and_D4=RR,D4_both_edges_exact_after_patch=full,
                           D4_prior_neutral_hull=old4,held_signs_D2_D4=found,
                           held_nets_D2_and_D4=[(n*v[2]-v[1,::-1]).tolist() for v in (a,b)]))
    W = [(j,k,s) for j,r in enumerate(result) for k,s in enumerate(r['held_signs_D2_D4']) if s]
    low=max((result[j]['neutral_D2_and_D4'][k][1],j,k) for j,k,s in W if s>0)
    high=min((result[j]['neutral_D2_and_D4'][k][0],j,k) for j,k,s in W if s<0)
    common=strictly(low[0],high[0])
    if common is None or not common[0] < n < common[1]:
        raise ArithmeticError('Resolved sign vector lacks a shared held corridor')
    sequence = [row['neutral_D2_and_D4'][1] for row in result[2:]]
    gaps = [b[0]-a[1] for a,b in zip(sequence[:-1],sequence[1:],strict=True)]
    if min(gaps)<=0.:
        raise ArithmeticError('D4 neutral root ordering above 3e-6 must be separated')
    return dict(rows=result,original_diagnostic_occupation=n,common_strict_occupation_corridor=common,
                adjoining_D4_neutral_gap_lower_above_q3=gaps,
                relative_occupation_margin_below_and_above=[(n-low[0])/n,(high[0]-n)/n],
                limiting_band_q=[result[j]['q'] for _,j,k in (low,high)],
                limiting_levels=[(2,4)[k] for _,j,k in (low,high)],
                compared_sha256=[hashlib.sha256(f.read_bytes()).hexdigest() for f in (p,q)],
                scope='Reuse an accepted fixed D4 exact-split ledger; all statements are held affine numerical signs, not live occupation or continuous integration error control.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
