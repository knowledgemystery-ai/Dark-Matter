"""Separate the count part of the inherited scalar donor effect after actual restart."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_count_only_cross_terminal import run as setup_run
from early_pole_endpoint_crossed_charge_terminal_history import prepare


def run(root,history,cross,full_cross,records,payloads):
    crossed_bytes=cross.read_bytes()
    B=json.loads(crossed_bytes)
    K=B['quantities']
    collect=[]
    for c,f,r,p in zip((2,4),full_cross,records,payloads,strict=True):
        setup=setup_run(root,history,f,c,p,setup_only=True)
        d,_,mask,_=prepare(root,history,c,6-c)
        blob=r.read_bytes()
        doc=json.loads(blob)
        q=json.loads(p.read_bytes())
        if (B['crossed_case_receipt_sha256'][c//2-1]!=hashlib.sha256(f.read_bytes()).hexdigest()
            or any(doc[k]!=v for k,v in setup.items()) or q['completed']!=2048
            or q['payload_sha256']!=drive_checkpoint_digest(q) or doc['completed_payload_sha256']!=q['payload_sha256']
            or doc['checks']!=q['checks'] or len(doc['checks'])!=8 or abs(q['x_reached']-.05)>2e-17):
            raise ArithmeticError('Exact count-only and full crossed receipts required')
        y=np.asarray(q['occupation'])
        if y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0. or max(y[d.ns:])>1.:
            raise ArithmeticError('Complete physically feasible numeric checkpoint required')
        z=d.ds.w*y[:d.ns]
        V=dict(cold_count=fsum(z[mask]),whole_scalar_count=fsum(z),pole=d.state(.05,y)[0])
        if any(doc[q]!=V[k] for q,k in (('endpoint_cold_count','cold_count'),('endpoint_whole_scalar_count','whole_scalar_count'),('endpoint_pole','pole'))):
            raise ArithmeticError('Exact count-only endpoint did not replay')
        collect.append((V,hashlib.sha256(blob).hexdigest(),doc['completed_payload_sha256']))
    out={}
    for key,y in K.items():
        a,b=[z[0][key] for z in collect]
        shifted=fsum((a-y['native_D2'],y['native_D4']-b))*.5
        full=y['average_scalar_replacement']
        out[key]=dict(count_only_endpoint_native_charge2=a,count_only_endpoint_native_charge4=b,
                      path_averaged_count_only_replacement=shifted,full_scalar_replacement=full,
                      count_share_of_full_scalar_replacement=shifted/full,
                      path_averaged_residual_shape_and_warm_state=full-shifted,
                      count_only_minus_full_scalar_endpoints=[a-y['charge_D2_scalar_D4'],b-y['charge_D4_scalar_D2']])
    return dict(quantities=out,x_window=[.04,.05],complete_cross_sha256=hashlib.sha256(crossed_bytes).hexdigest(),
                count_only_source_receipt_sha256=[v[1] for v in collect],count_only_source_payload_sha256=[v[2] for v in collect],
                scope='Actual held cold-number intervention; the residual combines full warm and scalar shape replacement under one terminal restart only.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','history','cross'):
        p.add_argument(k,type=Path)
    for k in ('full_cross','records','payloads'):
        p.add_argument('--'+k.replace('_','-'),type=Path,nargs=2,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
