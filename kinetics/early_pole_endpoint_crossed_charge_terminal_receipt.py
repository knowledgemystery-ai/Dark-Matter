"""Reduce only two complete matched terminal cross trajectories."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_crossed_charge_terminal_history import prepare


def split(a,b,c,d):
    # a=D2+D2, b=D2 charges + D4 scalar, c=D4 charges + D2 scalar, d=D4+D4.
    charge=((c-a)+(d-b))*.5
    scalar=((b-a)+(d-c))*.5
    if abs(fsum((charge,scalar))-(d-a))>5e-16*max(abs(a),abs(b),abs(c),abs(d),1e-60):
        raise ArithmeticError('Path-averaged contrast did not telescope')
    return dict(native_D2=a,charge_D2_scalar_D4=b,charge_D4_scalar_D2=c,native_D4=d,
                net_diagonal_change=d-a,average_charge_replacement=charge,average_scalar_replacement=scalar,
                factorial_interaction=fsum((d,a,-c,-b)))


def run(root,history,gap,records,final_checkpoints):
    g=gap.read_bytes()
    G=json.loads(g)
    pairs=[]
    for ch,sc,rec,cp in zip((2,4),(4,2),records,final_checkpoints,strict=True):
        d,_,mask,setup=prepare(root,history,ch,sc)
        f,h=rec.read_bytes(),cp.read_bytes()
        a,b=map(json.loads,(f,h))
        if any(a[k]!=v for k,v in setup.items()):
            raise ArithmeticError('Crossed record does not use the exact predetermined setting')
        if (b['payload_sha256']!=drive_checkpoint_digest(b) or a['complete_checkpoint_payload_sha256']!=b['payload_sha256']
            or b['target_steps']!=2048 or b['completed']!=2048 or abs(b['x_reached']-.05)>2e-17
            or b['checks']!=a['checks'] or len(b['checks'])!=8 or max(abs(u['continuity_residual']) for u in b['checks'])>3e-10):
            raise ArithmeticError('Two complete final checkpoints required')
        occ=np.asarray(b['occupation'])
        if occ.shape!=(d.n,) or not np.all(np.isfinite(occ)) or min(occ)<0 or max(occ[d.ns:])>1.:
            raise ArithmeticError('Invalid final state')
        pop=d.ds.w*occ[:d.ns]
        x=dict(cold_count=fsum(pop[mask]),whole_scalar_count=fsum(pop),pole=d.state(.05,occ)[0])
        for k,v in (('cold_count','endpoint_cold_count'),('whole_scalar_count','endpoint_scalar_count'),('pole','endpoint_scalar_pole')):
            if x[k]!=a[v]:
                raise ArithmeticError('Final optical/continuity work does not reproduce the receipt')
        pairs.append((x,hashlib.sha256(f).hexdigest(),b['payload_sha256'],a['canonical_input_states_sha256']))
    states={s['charge_level']:s for s in G['states']}
    if G['states_sha256_uncompressed'][1:3]!=[p[3][1] for p in pairs]:
        raise ArithmeticError('The exact native diagonal endpoints must be the matched gap states')
    return dict(x=[.04,.05],source_gap_sha256=hashlib.sha256(g).hexdigest(),
                crossed_case_receipt_sha256=[x[1] for x in pairs],crossed_state_payload_sha256=[x[2] for x in pairs],
                quantities={k:split(states[2][k],pairs[0][0][k],pairs[1][0][k],states[4][k]) for k in pairs[0][0]},
                scope='Conditional two-factor crossed terminal calculation; does not decompose earlier physical production or constrain convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','history','gap'):
        p.add_argument(k,type=Path)
    p.add_argument('--records',type=Path,nargs=2,required=True)
    p.add_argument('--final-checkpoints',type=Path,nargs=2,required=True)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
