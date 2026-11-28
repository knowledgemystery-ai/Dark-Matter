"""Cross normalized Bose-weighted shapes and bound their signed mean shifts."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import address, ledger, reweight, run as verify, strip
from early_pole_native_cold_cross_factors import CellChart, DATE, SPECS, load, read


def moment_interval(bounds,signed):
    x=bounds*signed[:,None]
    return np.array([fsum(x.min(axis=1)),fsum(x.max(axis=1))])


def summarize(charts,grid):
    bs=sorted(set(np.concatenate([r[1] for r in charts]).tolist()))
    H=[]
    for _,ids,N in charts:
        H.append(np.array([[fsum(w[ids==i]) for i in bs] for w in N]))
    (A,B),(C,D)=H
    K=(C+D-A-B)/2
    F=(B+D-A-C)/2
    V=D-A
    tv=lambda v: fsum(abs(v))/2
    if max(abs(K+F-V))>2e-15:
        raise ArithmeticError('Signed normalized crossed contributions did not close')
    p=[]
    for j,r in enumerate(charts):
        b,_,N=r
        p.append([moment_interval(b,z) for z in (N[1]-N[0],((-1.)**(1-j))*(N[0]+N[1]),((-1.)**(1-j))*N[j])])
    P=np.sum(p,axis=0)
    j=int(np.argmin(V))
    return dict(shared_bin_spans=[[max(float(grid[i]),1e-6),min(float(grid[i+1]),1.5e-5)] for i in bs],
                native_bin_difference=V.tolist(),kernel_component_by_bin=K.tolist(),occupation_component_by_bin=F.tolist(),
                native_TV=tv(V),kernel_component_TV=tv(K),occupation_component_TV=tv(F),
                smallest_absolute_signed_bin_components=[min(abs(K)),min(abs(F))],
                all_bin_components_same_sign=bool(np.array_equal(np.sign(K),np.sign(F))),
                biggest_loss_bin=int(j),biggest_fraction_of_all_loss=-V[j]/tv(V),
                occupation_mean_shift_interval=(P[0]/2.).tolist(),
                kernel_mean_shift_interval=(P[1]/2.).tolist(),
                native_mean_shift_interval=P[2].tolist())


def run(root,receipt):
    raw=receipt.read_bytes()
    if verify(root)!=json.loads(raw):
        raise ArithmeticError('Require exact previous complete D2/D4 reweight receipt')
    lad=json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    grids=[]; hist=[]; path=[]
    for ix,k in enumerate((2,4),1):
        A,H=read(root/f'early-pole-cell-endpoint-fine-d{k}-x05-state-{DATE}.json.gz')
        if H!=lad['cards'][ix]['state_sha256_uncompressed']:
            raise ArithmeticError('Exact accepted chart required')
        d=CellChart(load('fine')[0],SPECS['fine'],k,A['reference_x'])
        hist.append(np.asarray(A['endpoint_occupation'][:d.ns]))
        grids.append(d.ds.ed)
        L=[root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz' for m in (1,2)]
        if [unpacked_hash(z) for z in L]!=lad['raw_uncompressed_sha256'][ix]:
            raise ArithmeticError('Both exact complete paired ledgers required')
        path.append(L)
    if not np.array_equal(*grids):
        raise ArithmeticError('Scalar bins must be common')
    report=[]
    for m in (0,1):
        out={}
        for label,j,k in (('forward_number',0,0),('inverse_number',0,1),('forward_centered_energy',2,0)):
            charts=[]
            for ix in (0,1):
                R=ledger(path[ix][m])
                b=np.array([r['bounds'] for r in R])
                ids=address(grids[0],b)
                u=strip(np.array([r['coeff'][j] for r in R]),hist[ix][ids])
                N=[reweight(u,n[ids])[:,k] for n in hist]
                charts.append((b,ids,[z/fsum(z) for z in N]))
            out[label]=summarize(charts,grids[0])
        report.append(out)
    for k in report[0]:
        a,b=[r[k] for r in report]
        if abs(a['native_TV']-b['native_TV'])>2e-10 or max(abs(np.array(a['occupation_mean_shift_interval'])-b['occupation_mean_shift_interval']))>2e-15:
            raise ArithmeticError('Paired normalized profiles do not resolve')
    return dict(by_mesh=report,reweight_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Normalized exact shared-bin effects and support-only mean bounds on two held numerical kernels; no unique physical allocation or spatial convergence.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','receipt'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
