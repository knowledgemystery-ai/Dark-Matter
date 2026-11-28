"""Replay the two native diagonal endpoints with the same driver as the crossed arms."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_startup_pilot import drive
from early_pole_endpoint_capacity_fast_step import step
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,load


def run(r,n):
    (a,ha),(b,hb)=[read(r/f'early-pole-cell-endpoint-fine-d{n}-x{x}-state-2026-11-28.json.gz') for x in ('04','05')]
    if not (a['x_end'],b['x_end'])==(.04,.05) or not a['steps']==b['steps']==2048:
        raise ArithmeticError('Exact archived .04 to .05 clock required')
    d=CellChart(load('fine')[0],SPECS['fine'],n,.005)
    u=np.asarray(a['endpoint_occupation']); v=np.asarray(b['endpoint_occupation'])
    C,z=drive(d,.05,2048,x_start=.04,initial=u,stepper=step)
    return dict(charge=n,native_exact_input_hashes=[ha,hb],endpoint_bitwise_equal=bool(np.array_equal(z,v)),
                maximum_occupation_absolute_defect=float(max(abs(z-v))),
                integrated_count_absolute_defect=abs(fsum(d.ds.w*(z-v)[:d.ns])),
                maximum_continuity_residual=max(abs(q['continuity_residual']) for q in C))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('r',type=Path)
    p.add_argument('n',type=int,choices=(2,4))
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
