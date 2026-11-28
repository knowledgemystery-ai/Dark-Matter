"""Locate the cold-count charge-resolution deficit in the archived time and scalar bins."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS, load

D='2026-11-28'


def run(root,gap):
    g=gap.read_bytes()
    W=json.loads(g)
    fixed=load('fine')[0]
    epochs=(.01,.02,.04,.05)
    levels=(1,2,4,8)
    d=[]; state=[]; ph=[]
    for n,tag in zip(levels,('','-d2','-d4','-d8'),strict=True):
        d.append(CellChart(fixed,SPECS['fine'],n,.005))
        u=[]; h=[]
        for x,steps in zip(epochs,(1024,2048,2048,2048),strict=True):
            s,k=read(root/f'early-pole-cell-endpoint-fine{tag}-x{round(100*x):02d}-state-{D}.json.gz')
            if (s['grid'],s['subcells_per_archived_cell'],s['counts'],s['x_end'],s['steps'],s['reference_x'])!=(
                 'fine',n,[d[-1].ns,d[-1].nf,d[-1].np],x,steps,.005):
                raise ArithmeticError('Only matched completed charge clocks eligible')
            y=np.asarray(s['endpoint_occupation'])
            if y.shape!=(d[-1].n,) or not np.all(np.isfinite(y)) or min(y)<0. or max(y[d[-1].ns:])>1.:
                raise ArithmeticError('Endpoint occupation not physical on the discrete chart')
            u.append(y)
            h.append(k)
        state.append(u)
        ph.append(h)
    if [v[-1] for v in ph]!=W['states_sha256_uncompressed'] or W['cold_q_upper']!=1e-5:
        raise ArithmeticError('Exact terminal gap study inputs required')
    edges=d[0].ds.ed
    cold=edges[1:]<=1e-5
    eq=np.flatnonzero(cold)
    w=d[0].ds.w
    if any(not np.array_equal(c.ds.ed,edges) or not np.array_equal(c.ds.w,w) for c in d):
        raise ArithmeticError('Scalar bin geometry must be identical through the entire ladder')
    comparisons=[]
    for k,(aa,bb) in enumerate(zip(state[:-1],state[1:],strict=True)):
        L,R=d[k:k+2]
        E=np.asarray([w*(b[:R.ns]-a[:L.ns]) for a,b in zip(aa,bb,strict=True)])
        total=np.array([fsum(v) for v in E[:,cold]])
        if total[-1]>=0 or not np.isclose(total[-1],W['states'][k+1]['cold_count']-W['states'][k]['cold_count'],
                                          rtol=0.,atol=2e-16):
            raise ArithmeticError('Cold count must replay current gap data')
        increments=np.diff(E,axis=0)
        accu=np.r_[total[0],[fsum(v) for v in increments[:,cold]]]
        if abs(fsum(accu)-total[-1])>4e-16:
            raise ArithmeticError('Common-bin charge deficit timeline does not telescope')
        abs_last=fsum(abs(x) for x in E[-1,cold])
        largest=np.argsort(-abs(E[-1,eq]),kind='stable')[:5]
        tail=E[-1,~cold]
        compare=[]
        for i,(x,z) in enumerate(zip(epochs[:-1],epochs[1:],strict=True)):
            def sector(A,b,c):
                return [fsum(z) for z in np.split(A.w*(c-b),(A.ns,A.iP))]
            a,b=np.asarray(sector(L,aa[i],aa[i+1])),np.asarray(sector(R,bb[i],bb[i+1]))
            J=b-a
            err=max(abs(a[0]-a[1]),abs(a[0]+a[2]),abs(b[0]-b[1]),abs(b[0]+b[2]))
            if err>3e-13 or abs(J[0]-fsum(increments[i]))>3e-13:
                raise ArithmeticError('Independent charged count replay failed')
            compare.append(dict(x=[x,z],scalar_daughter_parent_difference_of_net_changes=J.tolist(),
                                absolute_pair_count_residual=err))
        comparisons.append(dict(levels=list(levels[k:k+2]),epochs=list(epochs),cold_count_differences=total.tolist(),
                                timeline_increment_fraction_of_terminal_cold_difference=(accu/total[-1]).tolist(),
                                cold_final_total_absolute_bin_difference=abs_last,
                                terminal_signed_to_absolute_cold_bin_ratio=total[-1]/abs_last,
                                warm_final_count_difference=fsum(tail),
                                largest_terminal_cold_differences=[dict(q=edges[j:j+2].tolist(),delta=E[-1,j],
                                                                        fraction_signed=E[-1,j]/total[-1])
                                                                   for j in eq[largest]],conserved_sector_checks=compare))
    return dict(rows=comparisons,initial_cold_count_group_end=1e-5,time_states_sha256=ph,
                gap_sha256=hashlib.sha256(g).hexdigest(),
                scope='Archived independently evolved native trajectories; endpoint differences telescope in time and shared scalar bins but are not errors or causal collision derivatives.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('gap',type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
