"""Compare real algebra against the verified complex factorization at native q settings."""
import argparse
import json
import time
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field, continuum_moments, model_callback, piecewise_logit, root as pole_root
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_timelike_pole_angular_reference import angle
from finite_mass_qed_transverse_factorized_spinor_reference import transverse_pair
from finite_mass_qed_transverse_real_pauli_reference import transverse_real_pauli


def run(root,level,timing=False):
    q,mesh={4:(4.1e-6,(4,10,18)),8:(4.35e-6,(8,16,24))}[level]
    A,h=read(root/f'early-pole-cell-endpoint-fine-d{level}-x05-state-2026-11-28.json.gz')
    if (A['grid'],A['subcells_per_archived_cell'],A['steps'],A['x_end'])!=('fine',level,2048,.05):
        raise ArithmeticError('Use the exact independently evolved fine native endpoint')
    d=CellChart(load('fine')[0],SPECS['fine'],level,A['reference_x'])
    v=np.asarray(A['endpoint_occupation'])
    mass=d.state(.05,v)[0]
    ch,cb,fit=[],[],[]
    for m,e,n in ((.05,d.p_edges,d.wp*v[d.iP:]),(.05*d.ratio,d.d_edges,d.wd*v[d.ns:d.iP])):
        p,w,_,_,a,b=continuum_moments(e,m,n,36)
        fit.append((p,w,m))
        ch.append(piecewise_logit(m,e,a,b,n,1))
        cb.append(model_callback(m,e,a,b))
    k=.0015
    _,w,R=pole_root(lambda x:sum(charge_field(*y,x,36) for y in fit),k)
    outcomes=[]; times=[]
    for f in (transverse_pair,transverse_real_pauli):
        start=time.perf_counter()
        outcomes.append(angle(ch,k,w,R,np.hypot(mass,q),'T',mesh[0],mesh[1:],q,occupation_callbacks=cb,
                              mirror_azimuth=True,transverse_vertex_kernel=f))
        times.append(time.perf_counter()-start)
    a,b=outcomes
    err=float(np.max(abs(b[0]/a[0]-1.)))
    if min(a[0])<=0. or not np.array_equal(a[1],b[1]) or a[2:]!=b[2:] or err>2e-12:
        raise ArithmeticError(('Real Pauli algebra did not replay the previous factorization',level,err))
    return dict(charge_level=level,q=q,k=k,mesh=mesh,mirror_azimuth=True,physical_native_state_sha256_uncompressed=h,
                photonic_pole=w,scalar_mass=mass,complex_and_real_affine_results=[a[0].tolist(),b[0].tolist()],
                unchanged_support=a[1].tolist(),relative_max_difference=err,
                **({'wall_seconds_complex_real':times} if timing else {}),
                scope='Same nodes and azimuth mirror, two algebraic spin sums; not a quadrature, cold band, or physics claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('level',type=int,choices=(4,8))
    p.add_argument('--timing',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
