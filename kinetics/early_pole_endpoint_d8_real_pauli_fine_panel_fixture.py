"""Verify one complete native D8 fine panel with independent optional shadow arithmetic."""
import argparse
import hashlib
import json
import time
from pathlib import Path
from tempfile import TemporaryDirectory

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import D8_CENTRAL_COLD_BAND,run
from early_pole_endpoint_scalar_subshell_chain_receipt import read


def fixture(state,stock,alternative,compute=False):
    A,h=read(state)
    raw,alt=[p.read_bytes() for p in (stock,alternative)]
    if any(len(q.splitlines())!=1 or not q.endswith(b'\n') for q in (raw,alt)):
        raise ArithmeticError('One exact completed raw line per method required')
    u,v=map(json.loads,(raw,alt))
    if u['panel']!=4 or v['panel']!=4 or tuple(u['mesh'])!=(8,10,16,24):
        raise ArithmeticError('Designated native stock panel unchanged')
    t,x=[run(A,2,None,topology_only=True,span=D8_CENTRAL_COLD_BAND,real_pauli_transverse=f,mirror_azimuth=f) for f in (False,True)]
    for row,top in ((u,t),(v,x)):
        if ({k:z for k,z in row.items() if k not in ('panel','bounds','coeff')}!={k:z for k,z in top.items() if k!='edges'}
            or row['bounds']!=top['edges'][4:6]):
            raise ArithmeticError('Source state, integrator, radial bins, or coefficient flags differ')
    a,b=map(np.asarray,[u['coeff'],v['coeff']])
    err=float(np.max(abs(b/a-1.)))
    if min(a.ravel())<=0. or min(b.ravel())<=0. or not np.isfinite(err) or err>2e-12:
        raise ArithmeticError('Same discrete quadrature did not reproduce exact panel coefficients')
    d=dict(exact_native_state_sha256=h,stock_raw_row_sha256=hashlib.sha256(raw).hexdigest(),
           real_pauli_row_sha256=hashlib.sha256(alt).hexdigest(),kept_native_mesh=u['mesh'],
           panel=4,bounds=u['bounds'],stock_coefficients=u['coeff'],real_pauli_mirror_coefficients=v['coeff'],max_relative_difference=err,
           method_fingerprint=x['real_pauli_transverse_sha256'],
           scope='One already completed stock native panel independently compared to the optional kernel. No D8 total, no ongoing worker change.')
    if not compute:
        return d
    with TemporaryDirectory() as tmp:
        out=Path(tmp)/'shadow.jsonl'
        start=time.perf_counter()
        run(A,2,out,span=D8_CENTRAL_COLD_BAND,panel_indices=[4],real_pauli_transverse=True,mirror_azimuth=True)
        if out.read_bytes()!=alt:
            raise ArithmeticError('Fresh optional full-panel replay is not byte-exact')
        return dict(**d,fresh_real_pauli_wall_seconds=time.perf_counter()-start,exact_replayed_raw_bytes=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('state','stock','alternative'):
        p.add_argument(name,type=Path)
    p.add_argument('--compute',action='store_true')
    print(json.dumps(fixture(**vars(p.parse_args())),indent=2,allow_nan=False))
