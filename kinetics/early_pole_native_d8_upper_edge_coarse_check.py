"""Check reflected real-Pauli arithmetic on the lower edge of the D8 raw plan."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_continuum_cold_band_pilot import run as pilot
from early_pole_native_cold_cross_factors import (
    DATE, K, CellChart, SPECS, angle, load, native, read, transverse_real_pauli)

PLAN_HASH = '9339f2abe2245acaa58cca77c4b37d6c99de76a4b3a96ea1df3c11dd91d830f9'
PANEL = 157


def run(root, scratch):
    p = (root/f'early-pole-native-d8-raw-band-plan-{DATE}.json').read_bytes()
    if hashlib.sha256(p).hexdigest()!=PLAN_HASH:
        raise ArithmeticError('Exact preregistered upper span required')
    R = json.loads(p)
    a,h = read(root/f'early-pole-cell-endpoint-fine-d8-x05-state-{DATE}.json.gz')
    M = pilot(a, 1, None, topology_only=True)
    topo = hashlib.sha256(json.dumps(M['edges'],separators=(',',':')).encode()).hexdigest()
    lo,hi = M['edges'][PANEL:PANEL+2]
    if (R['planned']['panel_indices'][0]!=PANEL or lo!=R['planned']['span'][0]
            or h!=R['d8_state_sha256_uncompressed'] or topo!=R['topology_sha256']
            or M['mesh']!=[4,6,10,18] or M['total']!=195):
        raise ArithmeticError('Need the held lower-edge upper panel on mesh one')
    c = native(root,8)
    mass,w,res = [c[k] for k in ('scalar_mass','omega','residue')]
    if (mass,w,res)!=(M['scalar_mass'],M['omega'],M['residue']):
        raise ArithmeticError('Independent reconstruction disagrees on fixed pole')
    d = CellChart(load('fine')[0],SPECS['fine'],8,a['reference_x'])
    occ = np.asarray(a['endpoint_occupation'][:d.ns])
    scratch.mkdir(parents=True,exist_ok=True)
    stock = scratch/'stock-linear.jsonl'
    pilot(a,1,stock,panel_indices=[PANEL])
    row = json.loads(stock.read_bytes())
    if ({k:v for k,v in row.items() if k not in ('panel','bounds','coeff')}
            != {k:v for k,v in M.items() if k!='edges'} or row['panel']!=PANEL or row['bounds']!=[lo,hi]):
        raise ArithmeticError('Wrong stock reference metadata')
    method = hashlib.sha256(Path(transverse_real_pauli.__code__.co_filename).read_bytes()).hexdigest()
    sig = dict(plan_sha256=PLAN_HASH,state_sha256=h,topology_sha256=topo,method_sha256=method,
               source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               mesh=M['mesh'],panel=PANEL,bounds=[lo,hi],stock_sha256=hashlib.sha256(stock.read_bytes()).hexdigest())
    V = np.zeros((3,3))
    x,y = roots_legendre(4)
    for i,(q,z) in enumerate(zip((lo+hi)/2+(hi-lo)*x/2,(hi-lo)*y/2,strict=True)):
        binding = dict(**sig,real_node=i,q=q,weight=z)
        s = scratch/f'edge-real-{i}.json'
        E = np.hypot(q,mass)
        if s.exists():
            raw=json.loads(s.read_bytes())
            if raw['binding']!=binding:
                raise ArithmeticError('Cannot mix edge checkpoints')
            A,B,D=raw['angle']
        else:
            A,B,D=angle(c['charges'],K,w,res,E,'T',6,(10,18),q,occupation_callbacks=c['callbacks'],
                        mirror_azimuth=True,transverse_vertex_kernel=transverse_real_pauli)[0]
            t=s.with_suffix('.tmp')
            t.write_text(json.dumps(dict(binding=binding,angle=[A,B,D]),allow_nan=False)+'\n')
            t.replace(s)
        if min(A,B,D)<=0. or abs((B+D)/A-1.)>6e-13:
            raise ArithmeticError('Independent affine node does not close')
        f=occ[np.searchsorted(d.ds.ed,q,side='right')-1]
        V += z*q*q/E * np.outer([1.,E,q*q/(E+mass)],[(1+f)*A,f*B,A+f*D])
    S=np.array(row['coeff'])
    diff=float(np.max(abs(V/S-1.)))
    if np.any(S<=0.0) or not np.isfinite(diff) or diff>2e-12:
        raise ArithmeticError('Coarse lower-edge canary failed arithmetic comparison')
    return dict(binding=sig,stock=row,real=V.tolist(),componentwise_max_relative_difference=diff,
                scope='One distinct first upper-plan panel at coarse mesh; no full-band or spatial claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','scratch'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
