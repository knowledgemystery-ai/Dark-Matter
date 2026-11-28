"""Run the already planned upper D8 support on separate real-Pauli meshes."""
import argparse
import fcntl
import hashlib
import json
from math import fsum
from pathlib import Path
from tempfile import TemporaryDirectory
from shutil import copytree

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_continuum_cold_band_pilot import run as stock
from early_pole_native_cold_bose_reweight import address, strip
from early_pole_native_cold_cross_factors import CellChart, DATE, K, SPECS, angle, load, native, read, transverse_real_pauli
from early_pole_native_d8_upper_edge_coarse_check import PLAN_HASH, run as coarse_canary
from early_pole_native_d8_upper_fine_canary import run as fine_canary
from finite_mass_qed_timelike_pole_angular_reference import radial

INDICES = tuple(range(157,195))
COMPLETE_TOL = 6e-5


def digest(b):
    return hashlib.sha256(b).hexdigest()


def wire(v):
    return json.dumps(v,sort_keys=True,separators=(',',':'),allow_nan=False).encode()


def save(p, z):
    payload = dict(record=z,sha256=digest(wire(z)))
    tmp = p.with_suffix('.tmp')
    tmp.write_bytes(wire(payload)+b'\n')
    tmp.replace(p)


def load_row(p, contract):
    z = json.loads(p.read_bytes())
    if z['sha256']!=digest(wire(z['record'])) or z['record']['contract']!=contract:
        raise ArithmeticError('Mixed or corrupted upper-band checkpoint')
    return z['record']['angle']


def context(root, mesh):
    p = (root/f'early-pole-native-d8-raw-band-plan-{DATE}.json').read_bytes()
    h0 = digest(p)
    if h0 != PLAN_HASH:
        raise ArithmeticError('Exact predeclared raw band plan required')
    P = json.loads(p)
    card,h = read(root/f'early-pole-cell-endpoint-fine-d8-x05-state-{DATE}.json.gz')
    M = stock(card,mesh,None,topology_only=True)
    edges = M['edges']
    hh = digest(json.dumps(edges,separators=(',',':')).encode())
    if (h!=P['d8_state_sha256_uncompressed'] or hh!=P['topology_sha256'] or M['total']!=195
        or P['planned']['panel_indices']!=list(INDICES) or P['planned']['span']!=[edges[157],edges[195]]):
        raise ArithmeticError('D8 upper native state, topology or span changed')
    C = native(root,8)
    if (C['card_hash']!=h or any(C[k]!=M[k] for k in ('scalar_mass','omega','residue'))):
        raise ArithmeticError('Independent reconstruction must match stock pole')
    d = CellChart(load('fine')[0],SPECS['fine'],8,card['reference_x'])
    f = np.asarray(card['endpoint_occupation'][:d.ns])
    src = [Path(__file__),*[Path(x.__code__.co_filename) for x in (stock,native,angle,radial,transverse_real_pauli)]]
    method = dict(**{k:v for k,v in M.items() if k!='edges'},plan_sha256=h0,state_sha256=h,topology_sha256=hh,
                  mirror_azimuth=True,source_sha256_by_module={p.name:digest(p.read_bytes()) for p in src})
    return method,edges,d.ds.ed,f,C


def panel(root, mesh, work, j, c, compute, seeded=None):
    meta,edges,scalar_edges,f,charges = c
    l,r = edges[j:j+2]
    mass,w,res = [meta[k] for k in ('scalar_mass','omega','residue')]
    n,outer,inner,az = meta['mesh']
    x,v = roots_legendre(n)
    bin_ = address(scalar_edges,[[l,r]])[0]
    f0 = f[bin_]
    out = np.zeros((3,3))
    for i,(q,z) in enumerate(zip((l+r)/2+(r-l)*x/2,(r-l)*v/2,strict=True)):
        path = work/f'upper-{j:03d}-{i:02d}.json'
        bind = dict(**meta,bounds=[l,r],panel=j,node=i,q=q,weight=z)
        E = np.hypot(q,mass)
        if path.exists():
            a = load_row(path,bind)
        elif seeded is not None:
            a = seeded[i]
            save(path,dict(contract=bind,angle=a))
        elif compute:
            a = angle(charges['charges'],K,w,res,E,'T',outer,(inner,az),q,
                      occupation_callbacks=charges['callbacks'],mirror_azimuth=True,
                      transverse_vertex_kernel=transverse_real_pauli)[0].tolist()
            save(path,dict(contract=bind,angle=a))
        else:
            raise FileNotFoundError(path)
        F,B,S = a
        if min(a)<=0. or not np.isfinite(a).all() or abs((B+S)/F-1.)>7e-13:
            raise ArithmeticError('Invalid positive affine upper node')
        out += z*q*q/E * np.outer([1.,E,q*q/(E+mass)],[(1.+f0)*F,f0*B,F+f0*S])
    return out,float(f0)


def seeds(root, mesh, work, c):
    prefixes = {1:('upper-edge-coarse','edge-real-*.json','angle','stock-linear.jsonl',157),
                2:('upper-fine-canary','real-node-*.json','raw_angular_values','stock.jsonl',183)}
    stem,pattern,key,stockfile,j = prefixes[mesh]
    if mesh==1:
        name=f'd8-{stem}-checkpoints-{DATE}'
        run=coarse_canary
        receipt=f'early-pole-native-d8-{stem}-check-{DATE}.json'
        comparison='stock'
    else:
        name=f'd8-{stem}-checkpoints-{DATE}'
        run=fine_canary
        receipt=f'early-pole-native-d8-{stem}-{DATE}.json'
        comparison='raw_stock'
    R=json.loads((root/receipt).read_bytes())
    with TemporaryDirectory() as td:
        src=root/name
        paths=sorted(src.glob(pattern))
        if len(paths)!=c[0]['mesh'][0] or len((src/stockfile).read_bytes().splitlines())!=1:
            raise ArithmeticError('Complete accepted reference checkpoints required')
        tmp=Path(td)/'canary'
        copytree(src,tmp)
        if run(root,tmp)!=R:
            raise ArithmeticError('Canary receipt must reproduce exactly before any upper run')
        data=[json.loads(p.read_bytes())[key] for p in paths]
    out,_=panel(root,mesh,work,j,c,False,data)
    if np.max(abs(out/np.asarray(R[comparison]['coeff'])-1.))>2e-12:
        raise ArithmeticError('Converted exact real-node array did not match the accepted stock panel')
    return j


def work(root, scratch, mesh, only_seed=False):
    scratch.mkdir(parents=True,exist_ok=True)
    C=context(root,mesh)
    with (scratch/'worker.lock').open('a+b') as f:
        fcntl.flock(f.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        j0=seeds(root,mesh,scratch,C)
        print(json.dumps(dict(mesh=mesh,seed=j0,contract=C[0])),flush=True)
        if not only_seed:
            # Alternate edges to expose both ends early without spawning more workers.
            order=[j for pair in zip(INDICES[:19],INDICES[19:][::-1],strict=True) for j in pair]
            for j in order:
                if j!=j0:
                    panel(root,mesh,scratch,j,C,True)
                print(json.dumps(dict(mesh=mesh,complete_panel=j)),flush=True)


def reduce(root, scratch):
    W=[]; V=[]; meshes=[]; hashes=[]
    for m in (1,2):
        d=scratch/f'mesh{m}'
        C=context(root,m)
        if len(list(d.glob('upper-*.json')))!=len(INDICES)*C[0]['mesh'][0]:
            raise ArithmeticError('Do not reduce a partial real-Pauli upper support')
        rec=[panel(root,m,d,j,C,False) for j in INDICES]
        A=np.array([v[0] for v in rec])
        B=np.array([v[1] for v in rec])
        raw=np.stack([strip(A[:,i],B) for i in range(3)],axis=1)
        if (np.max(abs((C[0]['scalar_mass']*raw[:,0]+raw[:,2])/raw[:,1]-1.))>3e-13
                or np.max(abs((raw[:,:,1]+raw[:,:,2])/raw[:,:,0]-1.))>6e-13):
            raise ArithmeticError('Raw moments violate affine or centered-energy closure')
        W.append(A); V.append(raw); meshes.append(C[0])
        hashes.append({p.name:digest(p.read_bytes()) for p in sorted(d.glob('upper-*.json'))})
    collapse=lambda z: [[fsum(z[:,i,j]) for j in range(3)] for i in range(3)]
    raw=np.array([collapse(z) for z in V])
    delta=float(np.max(abs(raw[0]/raw[1]-1.)))
    panel_delta=float(np.max(abs(V[0]/V[1]-1.)))
    if max(delta,panel_delta)>COMPLETE_TOL:
        raise ArithmeticError(('Underresolved preregistered two-mesh upper quadrature',delta,panel_delta))
    return dict(panel_ids=list(INDICES),mesh_contracts=meshes,by_mesh_raw_totals=raw.tolist(),
                by_mesh_weighted_totals=[collapse(z) for z in W],validated_node_file_sha256=hashes,
                raw_total_max_relative_difference=delta,raw_per_panel_max_relative_difference=panel_delta,
                scope='Only the complete selected upper raw band on fixed D8 state; neither central-weighted nor full-native total, no D8 coverage fraction, physical occupation or spatial error.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('scratch',type=Path)
    p.add_argument('mode',choices=('mesh1','mesh2','seed1','seed2','reduce'))
    a=p.parse_args()
    if a.mode=='reduce':
        print(json.dumps(reduce(a.root,a.scratch),indent=2,allow_nan=False))
    else:
        m=int(a.mode[-1])
        work(a.root,a.scratch/f'mesh{m}',m,a.mode.startswith('seed'))
