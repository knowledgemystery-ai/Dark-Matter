"""Cross scalar pole, photon pole and charged shape on exact D2 and D4 states."""
import argparse
import hashlib
import json
from itertools import combinations, product
from math import comb, fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import (
    charge_field, continuum_moments, model_callback, piecewise_logit, root as pole)
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS, load
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_timelike_pole_batched_absorption_reference import spacelike_parent_minimum
from finite_mass_qed_transverse_real_pauli_reference import transverse_real_pauli

DATE='2026-11-28'
K=.0015
MESHES=((4,10,18),(8,16,24))
FACTORS=('scalar_pole','photon_energy_and_residue','charged_distributions')


def native(root,n):
    path=root/f'early-pole-cell-endpoint-fine-d{n}-x05-state-{DATE}.json.gz'
    A,h=read(path)
    if (A['grid'],A['subcells_per_archived_cell'],A['steps'],A['x_end'])!=('fine',n,2048,.05):
        raise ArithmeticError('Exact independently evolved native fine state required')
    d=CellChart(load('fine')[0],SPECS['fine'],n,A['reference_x'])
    v=np.asarray(A['endpoint_occupation'])
    mu=d.state(.05,v)[0]
    ch,cb,fit=[],[],[]
    for m,e,w in ((.05,d.p_edges,d.wp*v[d.iP:]),(.05*d.ratio,d.d_edges,d.wd*v[d.ns:d.iP])):
        q,f,_,_,a,b=continuum_moments(e,m,w,36)
        fit.append((q,f,m))
        ch.append(piecewise_logit(m,e,a,b,w,1))
        cb.append(model_callback(m,e,a,b))
    _,w,R=pole(lambda z:sum(charge_field(*row,z,36) for row in fit),K)
    return dict(card_hash=h,scalar_mass=mu,omega=w,residue=R,daughter_mass=.05*d.ratio,charges=ch,callbacks=cb)


def fractions(values):
    w={S:np.log(np.asarray(z)) for S,z in values.items()}
    order=tuple(range(3))
    if set(w)!={S for k in range(4) for S in combinations(order,k)} or not np.isfinite(list(w.values())).all():
        raise ArithmeticError('Complete positive three-factor outcomes required')
    a=[]
    for j in order:
        a.append(np.array([fsum((w[tuple(sorted((*S,j)))][i]-z[i])/(3*comb(2,len(S)))
                                for S,z in w.items() if j not in S) for i in range(len(w[()]))]))
    full=w[order]-w[()]
    if np.max(abs(np.sum(a,axis=0)-full))>1e-13:
        raise ArithmeticError('All-order crossed allocations fail product closure')
    return np.asarray(a),full


def evaluate(x,mask,mesh,q):
    m,w,c=[x[int(i in mask)] for i in range(3)]
    E=np.hypot(m['scalar_mass'],q)
    row=angle(c['charges'],K,w['omega'],w['residue'],E,'T',mesh[0],mesh[1:],q,
              occupation_callbacks=c['callbacks'],mirror_azimuth=True,
              transverse_vertex_kernel=transverse_real_pauli)[0]
    if not np.isfinite(row).all() or min(row)<=0:
        raise ArithmeticError('Crossed coefficient row is not positive finite')
    t=E-w['omega']; V=np.hypot(K,q)
    cut=spacelike_parent_minimum(.05,c['daughter_mass'],t,V,(t-V)*(t+V))
    return dict(raw_affine_coefficients=row.tolist(),parent_threshold_at_cosine_zero=cut)


def run(root,q,scratch):
    cards=[native(root,n) for n in (2,4)]
    identity=[{k:v for k,v in z.items() if k not in ('charges','callbacks')} for z in cards]
    a=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    key=dict(source_sha256=a,raw_native_states=identity,q=q,k=K,angular_meshes=[list(m) for m in MESHES])
    scratch.mkdir(parents=True,exist_ok=True)
    data=[{},{}]
    rec=[]
    for maskbits in product((0,1),repeat=3):
        mask=tuple(j for j,b in enumerate(maskbits) if b)
        p=''.join(map(str,maskbits)); rows=[]
        for k,mesh in enumerate(MESHES):
            f=scratch/f'factor-{p}-mesh{k+1}.json'
            card=dict(**key,mask=list(mask),mesh=list(mesh))
            if f.exists():
                saved=json.loads(f.read_bytes())
                if saved['contract']!=card:
                    raise ArithmeticError('Independent point checkpoint refuses mismatched sources')
                z=saved['result']
            else:
                z=evaluate(cards,mask,mesh,q)
                temp=f.with_suffix('.tmp')
                temp.write_text(json.dumps(dict(contract=card,result=z),allow_nan=False)+'\n')
                temp.replace(f)
            rows.append(z)
            data[k][mask]=z['raw_affine_coefficients']
        x,y=[np.asarray(z['raw_affine_coefficients']) for z in rows]
        delta=float(max(abs(x/y-1.)))
        rec.append(dict(mask=list(mask),rows=rows,angular_relative_error=delta))
    if max(v['angular_relative_error'] for v in rec)>1e-5:
        raise ArithmeticError('Full crossed angular quadrature is not resolved')
    A,F=fractions(data[1])
    B,_=fractions(data[0])
    return dict(contract=key,factors=list(FACTORS),rows=rec,
                changes=dict(all_order_log=A.tolist(),all_order_multipliers=np.exp(A).tolist(),
                             total_log=F.tolist(),total_multipliers=np.exp(F).tolist(),
                             max_two_mesh_log_allocation_defect=float(np.max(abs(A-B)))),
                scope='Frozen D2/D4 point comparisons with mathematical crossed inputs; not trajectory attribution, whole-band convergence, or physical production.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('q',type=float,choices=(4.1e-6,8.2e-6))
    p.add_argument('scratch',type=Path)
    a=p.parse_args()
    print(json.dumps(run(a.root,a.q,a.scratch),indent=2,allow_nan=False))
