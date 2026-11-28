"""Reduce the two exact shape-null terminal trials against the nominal remap."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_charge_nested_moment_remap import prepare


def run(root,clock,remap,records,payloads):
    c,b=clock.read_bytes(),remap.read_bytes()
    hc,hb=[hashlib.sha256(z).hexdigest() for z in (c,b)]
    R=json.loads(b)['results']
    d,_,P=prepare(root,2)
    out=[]
    for s,g,h in zip((-1,1),records,payloads,strict=True):
        z=g.read_bytes()
        a,t=json.loads(z),json.loads(h.read_bytes())
        if (a['sign']!=s or a['fraction_of_feasible_move']!=s*.05 or a['setup']!=P or a['source_inputs'][0]!=hc
            or t['completed']!=2048 or t['payload_sha256']!=drive_checkpoint_digest(t) or t['payload_sha256']!=a['payload_sha256']
            or t['checks']!=a['checks'] or len(t['checks'])!=8 or abs(t['x_reached']-.05)>2e-17):
            raise ArithmeticError('Require both exact complete alternate shape histories')
        u=np.asarray(t['occupation'])
        if u.shape!=(d.n,) or not np.all(np.isfinite(u)) or min(u)<0. or max(u[d.ns:])>1.:
            raise ArithmeticError('Final column violates occupation admissibility')
        p=d.ds.w*u[:d.ns]
        V=dict(cold_count=fsum(p[d.ds.ed[1:]<=1e-5]),whole_scalar_count=fsum(p),pole=d.state(.05,u)[0])
        if any(V[key]!=a[k] for key,k in (('cold_count','final_cold_count'),('whole_scalar_count','final_count'),('pole','final_scalar_pole'))):
            raise ArithmeticError('Independent nonlinear terminal moments failed replay')
        diff={k:q-R[k]['grid4_donor2'] for k,q in V.items()}
        out.append(dict(fraction=s*.05,endpoint=V,difference_from_nominal_up_remap=diff,
                        absolute_difference_fraction_of_grid_replacement={k:abs(v/R[k]['average_grid_and_representation_replacement'])
                                                                          for k,v in diff.items()},
                        raw_receipt_sha256=hashlib.sha256(z).hexdigest(),checkpoint_payload_sha256=t['payload_sha256']))
    return dict(x_window=[.04,.05],rows=out,static_clock_sha256=hc,nominal_terminal_remap_sha256=hb,
                scope='Two specific statically prespecified moment-preserving charged-shape moves. Not an exhaustive higher-moment remap bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','clock','remap'):
        p.add_argument(k,type=Path)
    for k in ('records','payloads'):
        p.add_argument('--'+k,type=Path,nargs=2,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
