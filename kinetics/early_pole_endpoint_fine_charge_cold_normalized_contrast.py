"""Separate normalized cold source shapes from the complete charge-ladder magnitude change."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_non_nested_cdf_contrast import compare
from early_pole_endpoint_cold_panel_transport_upper import most
from early_pole_endpoint_cold_profile_interval_transport import fixed, transport
from early_pole_endpoint_cold_normalized_profile_envelope import geometry
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash

D='2026-11-28'


def one(e,f,a,b):
    u,v=(z/fsum(z) for z in (a,b))
    low=transport(fixed(dict(base=u)),fixed(dict(base=v)),e,f,'fixed',independent_bracket=True)
    us=most(e,f,u,v)
    shape=compare(e,f,a,b)
    w=low['independent_cost_lower_upper_q']
    x=us['sharp_upper_q']
    if w[1]>x:
        raise ArithmeticError('Independent minimum transport enclosure exceeds compatible maximum')
    return dict(integral_ratio=shape['upper_over_lower_band_response'],
                strongest_cdf_cut=shape['strongest_q'],
                guaranteed_cdf_gap=shape['normalized_TV_lower_bound_via_cdf'],
                minimum_compatible_TV_bounds=shape['sharp_compatible_TV']['minimum_TV_bounds'],
                minimum_compatible_W1_q_bounds=w,maximum_compatible_W1_q=x,
                exact_shared_inner_cuts=np.intersect1d(e[1:-1],f[1:-1]).tolist())


def run(root,ladder):
    payload=ladder.read_bytes()
    L=json.loads(payload)
    paths=[[root/f'fine{tag}-endpoint-continuum-cold-{i}.jsonl.gz' for i in (1,2)] for tag in ('','-d2','-d4')]
    hashes=[[unpacked_hash(f) for f in z] for z in paths]
    if hashes!=L['raw_uncompressed_sha256'] or L['q_band']!=[1e-6,1.5e-5]:
        raise ArithmeticError('Only the exact complete accepted D1/D2/D4 ledgers are eligible')
    plans=[[geometry(f) for f in p] for p in paths]
    output=[]
    for l,r in ((0,1),(1,2)):
        stage=[]
        for (e,A),(f,B) in zip(plans[l],plans[r],strict=True):
            results={}
            for key,i,j in (('forward_number',0,0),('inverse_number',0,1),('forward_excess',2,0)):
                results[key]=one(e,f,A[:,i,j],B[:,i,j])
            stage.append(results)
        err={key:0. for key in ('minimum_compatible_TV_bounds','minimum_compatible_W1_q_bounds','maximum_compatible_W1_q')}
        for key in stage[1]:
            for attr in err:
                err[attr]=max(err[attr],float(np.max(abs(np.asarray(stage[1][key][attr])-stage[0][key][attr]))))
        if (err['minimum_compatible_TV_bounds']>2e-9 or err['minimum_compatible_W1_q_bounds']>2e-12
            or err['maximum_compatible_W1_q']>2e-12):
            raise ArithmeticError('Shape bounds are not stable across both independent integrals')
        output.append(dict(charge_levels=[2**l,2**r],finer_mesh=stage[1],maximum_two_mesh_defects=err))
    return dict(q_band=L['q_band'],pairs=output,input_raw_uncompressed_sha256=hashes,
                input_ladder_sha256=hashlib.sha256(payload).hexdigest(),
                scope='All normalized placements compatible with held positive native panels, not a physical shape error or magnitude-convergence statement.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('ladder',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
