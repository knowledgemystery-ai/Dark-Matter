"""Confirm pointwise transverse factorization on exact native D4 and D8 continuum states."""
import argparse
import hashlib
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

D='2026-11-28'


def run(root,timing=False,mirror_factorized=False):
    rows=[]
    cases=(((8,4.2e-6,(4,10,18)),) if mirror_factorized else
           ((4,3.5e-6,(4,10,18)),(8,4e-6,(4,10,18)),(8,4.8e-6,(8,16,24))))
    for c,q,mesh in cases:
        A,h=read(root/f'early-pole-cell-endpoint-fine-d{c}-x05-state-{D}.json.gz')
        if (A['grid'],A['subcells_per_archived_cell'],A['steps'],A['x_end'])!=('fine',c,2048,.05):
            raise ArithmeticError('Exact native fine states for this narrow benchmark required')
        d=CellChart(load('fine')[0],SPECS['fine'],c,A['reference_x'])
        y=np.asarray(A['endpoint_occupation'])
        M=d.state(.05,y)[0]
        ch=[]; cb=[]; fits=[]
        for m,b,z in ((.05,d.p_edges,d.wp*y[d.iP:]),(.05*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
            p,w,_,_,a,beta=continuum_moments(b,m,z,36)
            fits.append((p,w,m))
            ch.append(piecewise_logit(m,b,a,beta,z,1))
            cb.append(model_callback(m,b,a,beta))
        k=.0015
        _,w,R=pole_root(lambda v:sum(charge_field(*j,v,36) for j in fits),k)
        time_log=[]; ret=[]
        for f in (None,transverse_pair):
            now=time.perf_counter()
            ret.append(angle(ch,k,w,R,np.hypot(M,q),'T',mesh[0],mesh[1:],q,occupation_callbacks=cb,transverse_vertex_kernel=f,mirror_azimuth=mirror_factorized and f is not None))
            time_log.append(time.perf_counter()-now)
        a,b=ret
        error=float(np.max(abs(a[0]/b[0]-1.)))
        if (not np.array_equal(a[1],b[1]) or a[2:]!=b[2:] or np.min(a[0])<=0. or error>1e-12):
            raise ArithmeticError(('Factorized transverse spin sum did not replay exact support and coefficients',c,q,error))
        entry=dict(charge_level=c,q=q,mesh=mesh,state_sha256_uncompressed=h,k=k,photon_pole=w,scalar_mass=M,
                   before_after_coefficients=[a[0].tolist(),b[0].tolist()],relative_defect=error,
                   identical_radial_node_support=a[1].tolist(),angular_panels=a[2])
        if timing:
            entry['wall_seconds_reference_factorized']=time_log
        rows.append(entry)
    return dict(rows=rows,**({'factorized_uses_mirror_azimuth':True} if mirror_factorized else {}),scope='Algebraic factorization of the same transverse discrete angular quadratures, not a new quadrature validation or D8 band result.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('--timing',action='store_true')
    p.add_argument('--mirror-factorized',action='store_true')
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
