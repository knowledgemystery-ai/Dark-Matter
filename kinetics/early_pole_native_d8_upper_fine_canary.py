"""Compare stock and exact real Pauli on one disjoint native upper D8 panel."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_continuum_cold_band_pilot import run as pilot
from early_pole_native_cold_cross_factors import (
    CellChart, DATE, K, SPECS, angle, load, native, read, transverse_real_pauli)


def run(root,scratch):
    names=('native-d8-raw-band-plan','native-d8-upper-point')
    raw=[(root/f'early-pole-{s}-{DATE}.json').read_bytes() for s in names]
    P,Q=map(json.loads,raw)
    ids=[hashlib.sha256(x).hexdigest() for x in raw]
    card,h=read(root/f'early-pole-cell-endpoint-fine-d8-x05-state-{DATE}.json.gz')
    if h!=P['d8_state_sha256_uncompressed'] or ids[0]!=Q['input_sha256'][1]:
        raise ArithmeticError('Exact admitted native state and upper point required')
    M=pilot(card,2,None,topology_only=True)
    H=hashlib.sha256(json.dumps(M['edges'],separators=(',',':')).encode()).hexdigest()
    if H!=P['topology_sha256']:
        raise ArithmeticError('Exact full native D8 topology must replay')
    q=Q['q']
    i=int(np.searchsorted(M['edges'],q,side='right')-1)
    a,b=M['edges'][i:i+2]
    if i not in P['planned']['panel_indices'] or not 5e-6<a<q<b:
        raise ArithmeticError('Disjoint upper canary panel required')
    scratch.mkdir(parents=True,exist_ok=True)
    stock=scratch/'stock.jsonl'
    pilot(card,2,stock,panel_indices=[i])
    y=stock.read_bytes()
    if len(y.splitlines())!=1 or not y.endswith(b'\n'):
        raise ArithmeticError('Exactly one stock row must finish')

    c=native(root,8)
    m,w,R=[c[k] for k in ('scalar_mass','omega','residue')]
    if (m,w,R)!=(M['scalar_mass'],M['omega'],M['residue']):
        raise ArithmeticError('Independent state preparation disagrees with stock pole')
    d=CellChart(load('fine')[0],SPECS['fine'],8,card['reference_x'])
    f=np.asarray(card['endpoint_occupation'][:d.ns])
    n,outer,inner,az=M['mesh']
    x,g=roots_legendre(n)
    qn=(a+b)/2+(b-a)*x/2
    gn=(b-a)*g/2
    method=hashlib.sha256(Path(transverse_real_pauli.__code__.co_filename).read_bytes()).hexdigest()
    sig=dict(source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             exact_native_state_sha256=h,full_topology_sha256=H,real_pauli_method_sha256=method,
             exact_stock_row_sha256=hashlib.sha256(y).hexdigest(),panel=i,bounds=[a,b],mesh=M['mesh'])
    W=np.zeros((3,3))
    for j,(q,z) in enumerate(zip(qn,gn,strict=True)):
        f0=f[np.searchsorted(d.ds.ed,q,side='right')-1]
        path=scratch/f'real-node-{j:02d}.json'
        contract=dict(**sig,node=j,q=q,outer_weight=z)
        if path.exists():
            s=json.loads(path.read_bytes())
            if s['contract']!=contract:
                raise ArithmeticError('No mixing canary checkpoint versions')
            A,B,D=s['raw_angular_values']
        else:
            A,B,D=angle(c['charges'],K,w,R,np.hypot(q,m),'T',outer,(inner,az),q,
                          occupation_callbacks=c['callbacks'],mirror_azimuth=True,transverse_vertex_kernel=transverse_real_pauli)[0]
            tmp=path.with_suffix('.tmp')
            tmp.write_text(json.dumps(dict(contract=contract,raw_angular_values=[A,B,D]),allow_nan=False)+'\n')
            tmp.replace(path)
        E=np.hypot(q,m)
        u=z*q*q/E*np.array([(1+f0)*A,f0*B,f0*D+A])
        W[0]+=u; W[1]+=E*u; W[2]+=q*q/(E+m)*u

    v=json.loads(y)
    if (v['panel']!=i or v['bounds']!=[a,b] or
        {k:x for k,x in v.items() if k not in ('panel','bounds','coeff')}!={k:x for k,x in M.items() if k!='edges'}):
        raise ArithmeticError('Stock metadata does not match selected upper canary')
    V=np.asarray(v['coeff'])
    defect=float(np.max(abs(W/V-1.)))
    if min(V.ravel())<=0. or min(W.ravel())<=0. or not np.isfinite(defect) or defect>2e-12:
        raise ArithmeticError(('Exact upper panel arithmetic is not yet corroborated',defect))
    return dict(contract=sig,raw_stock=v,exact_real_pauli_mirror_coefficients=W.tolist(),
                maximum_algorithm_relative_difference=defect,selected_reference_point=Q['q'],input_sha256=ids,
                scope='One complete disjoint exact weighted upper panel, no upper-band rate or change to active central calculations.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','scratch'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
