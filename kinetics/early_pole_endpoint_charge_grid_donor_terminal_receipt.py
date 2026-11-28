"""Compare charge-grid and charged-state replacements with the D2 scalar start fixed."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_charge_nested_moment_remap import prepare
from early_pole_endpoint_crossed_charge_terminal_receipt import split


def run(root,cross,records,payloads):
    D=cross.read_bytes()
    old=json.loads(D)['quantities']
    row={}
    refs=[]
    for j,r,p in zip((2,4),records,payloads,strict=True):
        chart,_,setup=prepare(root,j)
        data,blob=r.read_bytes(),p.read_bytes()
        a,b=map(json.loads,(data,blob))
        if (any(a[k]!=v for k,v in setup.items()) or b['completed']!=2048
            or b['payload_sha256']!=drive_checkpoint_digest(b) or a['final_payload_sha256']!=b['payload_sha256']
            or b['checks']!=a['checks'] or len(b['checks'])!=8
            or max(abs(z['continuity_residual']) for z in b['checks'])>3e-10):
            raise ArithmeticError('The exact two complete terminal remap arcs are required')
        occ=np.asarray(b['occupation'])
        if occ.shape!=(chart.n,) or not np.all(np.isfinite(occ)) or min(occ)<0. or max(occ[chart.ns:])>1.:
            raise ArithmeticError('Final Pauli/positivity validation failed')
        z=chart.ds.w*occ[:chart.ns]
        t=chart.ds.ed[1:]<=1e-5
        row[j]=dict(cold_count=fsum(z[t]),whole_scalar_count=fsum(z),pole=chart.state(.05,occ)[0])
        if any(a[i]!=row[j][k] for i,k in (('cold_count','cold_count'),('scalar_count','whole_scalar_count'),('scalar_pole_final','pole'))):
            raise ArithmeticError('Recomputed occupation moments failed exact replay')
        refs.append(dict(record_sha256=hashlib.sha256(data).hexdigest(),payload_sha256=b['payload_sha256']))
    results={}
    for k,v in old.items():
        Q=split(v['native_D2'],row[4][k],row[2][k],v['charge_D4_scalar_D2'])
        results[k]=dict(grid2_donor2=Q['native_D2'],grid2_donor4=Q['charge_D2_scalar_D4'],
                        grid4_donor2=Q['charge_D4_scalar_D2'],grid4_donor4=Q['native_D4'],
                        diagonal_charge_difference=Q['net_diagonal_change'],
                        average_grid_and_representation_replacement=Q['average_charge_replacement'],
                        average_inherited_charged_state_replacement=Q['average_scalar_replacement'],
                        interaction=Q['factorial_interaction'])
    return dict(results=results,held_scalar_donor=2,x_window=[.04,.05],remapped_references=refs,
                input_completed_cross_sha256=hashlib.sha256(D).hexdigest(),
                scope='Separates destination grid with its selected moment-preserving representation from inherited charge state under one controlled D2 scalar start. No unique higher-moment remap or physical convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('cross',type=Path)
    p.add_argument('--records',nargs=2,type=Path,required=True)
    p.add_argument('--payloads',nargs=2,type=Path,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
