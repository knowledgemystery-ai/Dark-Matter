"""Probe what the cold-count-only scalar root predicts for fixed point apertures."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field, continuum_moments, model_callback, piecewise_logit, root
from early_pole_endpoint_cold_aperture_attribution import K, coefficients
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS, load

D='2026-11-28'


def run(directory,gap,ladder):
    raw=gap.read_bytes()
    z=ladder.read_bytes()
    G,L=map(json.loads,(raw,z))
    if G['states_sha256_uncompressed'][:3]!=[card['state_sha256_uncompressed'] for card in L['cards']]:
        raise ArithmeticError('The same native evolved states must be used throughout')
    base=load('fine')[0]
    rows=[]
    for j,i in ((4,2),(8,3)):
        p=directory/f'early-pole-cell-endpoint-fine-d{j}-x05-state-{D}.json.gz'
        A,sha=read(p)
        if sha!=G['states_sha256_uncompressed'][i] or A['subcells_per_archived_cell']!=j or A['x_end']!=.05:
            raise ArithmeticError('Need exact new charge states from the gap study')
        at=G['rows'][i-1]
        if at['levels'][-1]!=j:
            raise ArithmeticError('Count substitution pair mismatch')
        d=CellChart(base,SPECS['fine'],j,A['reference_x'])
        y=np.asarray(A['endpoint_occupation'])
        m=d.state(.05,y)[0]
        if m!=at['pole_full_new']:
            raise ArithmeticError('Held new scalar mass does not replay')
        fits=[]; ch=[]; cb=[]
        for mass,e,n in ((.05,d.p_edges,d.wp*y[d.iP:]),(.05*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
            p,F,_,_,a,b=continuum_moments(e,mass,n,36)
            fits.append((p,F,mass))
            cb.append(model_callback(mass,e,a,b))
            ch.append(piecewise_logit(mass,e,a,b,n,1))
        _,w,R=root(lambda s:sum(charge_field(*u,s,36) for u in fits),K)
        if j==4 and L['cards'][2]['pole']['omega']!=w:
            raise ArithmeticError('D4 native photon state changed')
        for q in (3.5e-6,8.2e-6):
            # Same exact late charged state, pole, and q in both arms.
            old,old_err=coefficients(ch,m,w,R,q,cb)
            new,new_err=coefficients(ch,at['pole_cold_count_substitution'],w,R,q,cb)
            rows.append(dict(charge_subcells=j,q=q,true_scalar_pole=m,substituted_scalar_pole=at['pole_cold_count_substitution'],
                             photon_pole=w,k=K,raw_true=old.tolist(),raw_substituted=new.tolist(),
                             relative_difference=(new/old-1.).tolist(),coarse_fine_mesh_defects=[old_err,new_err]))
    return dict(rows=rows,fields=['forward','inverse','slope'],gap_sha256=hashlib.sha256(raw).hexdigest(),
                ladder_sha256=hashlib.sha256(z).hexdigest(),
                scope='Held cold count scalar mass substitute, two selected q tests, not a whole-band D8 result or evolutionary cause claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for v in ('directory','gap','ladder'):
        p.add_argument(v,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
