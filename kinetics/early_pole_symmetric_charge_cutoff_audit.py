"""Test exact terminal cumulative count contrasts at every native scalar edge."""
import argparse
import gzip
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_cell_endpoint_startup_pilot import drive_checkpoint_digest
from early_pole_endpoint_charge_symmetric_null_terminal import setup

DATE='2026-11-28'
NAMES=('charge-symmetric-time-nominal','charge-symmetric-time-lower',
       'charge-symmetric-time-upper','symmetric-half-lower','symmetric-half-upper')


def contrasts(edges, data, cell_scale):
    e=np.asarray(edges)
    V=np.asarray(data)
    a=np.asarray(cell_scale)
    if (e.ndim!=1 or len(e)==0 or np.any(np.diff(e)<=0) or np.any(e<=0.)
        or a.shape!=e.shape or V.shape!=(5,len(e)) or not np.isfinite(np.r_[e,a,V.ravel()]).all()
        or np.any(a<=0.) or np.any(V<0.)):
        raise ArithmeticError('Expected strictly ordered native cells and five nonnegative occupations')
    N=V[0]*a
    delta=(V[1:]-V[0])*a
    out=[]
    for i,x in enumerate(e,1):
        d=np.array([fsum(y[:i]) for y in delta])
        scale=max(fsum(N[:i]),max(fsum(y[:i]*a[:i]) for y in V[1:]),np.finfo(float).tiny)
        # This floor distinguishes resolvable saved-data sign from near cancellation.
        # It is not a trajectory, cell discretization, or independent physical error bound.
        floor=2e-14*scale
        code=lambda y: int(y>floor)-int(y<-floor)
        odd=[(d[j+1]-d[j])/2 for j in (0,2)]
        even=[(d[j+1]+d[j])/2 for j in (0,2)]
        out.append(dict(q=float(x),nominal=fsum(N[:i]),delta=d.tolist(),
                        full_and_half_odd=odd,full_and_half_even=even,
                        signs=[code(v) for v in d],sign_floor=floor))
    return out


def interior_enclosures(edges, data, cell_scale, rows):
    """Allow the five nonnegative distributions to split each cell independently."""
    masses=np.asarray(data)*np.asarray(cell_scale)
    out=[]
    for i,x in enumerate(edges):
        prev=np.asarray(rows[i-1]['delta']) if i else np.zeros(4)
        l,u=prev-masses[0,i],prev+masses[1:,i]
        floor=rows[i]['sign_floor']
        out.append(dict(q_interval=[float(edges[i-1]) if i else 0.,float(x)],
                        lower=l.tolist(),upper=u.tolist(),
                        certified_signs=[int(a>floor)-int(b<-floor) for a,b in zip(l,u,strict=True)]))
    return out


def run(root):
    d,_,_,D=setup(root)
    ids=[]; occ=[]; receipts=[]
    for name,sign in zip(NAMES,(0,-1,1,-1,1),strict=True):
        stem=root/f'early-pole-{name}-x05'
        rec=Path(f'{stem}-receipt-{DATE}.json').read_bytes()
        gz=Path(f'{stem}-checkpoint-{DATE}.json.gz').read_bytes()
        R,B=json.loads(rec),json.loads(gzip.decompress(gz))
        z=np.asarray(B['occupation'])
        h=B['payload_sha256']
        if (R['setup']!=D or R['sign']!=sign or R['steps']!=4096 or B['completed']!=4096
            or B['target_steps']!=4096 or abs(B['x_reached']-.05)>2e-17 or len(B['checks'])!=8
            or B['checks']!=R['checks'] or h!=drive_checkpoint_digest(B)
            or h!=R.get('complete_payload_sha256',R.get('payload_sha256'))
            or z.shape!=(d.n,) or not np.isfinite(z).all() or min(z)<0. or max(z[d.ns:])>1.):
            raise ArithmeticError(('Need the exact completed symmetric runs',name))
        f=d.ds.w*z[:d.ns]
        if (fsum(f)!=R['final_total'] or fsum(f[d.ds.ed[1:]<=1e-5])!=R['final_cold']
            or d.state(.05,z)[0]!=R['final_pole']):
            raise ArithmeticError(('Archived terminal state did not independently replay',name))
        ids.append(dict(name=name,receipt_sha256=hashlib.sha256(rec).hexdigest(),
                        compressed_checkpoint_sha256=hashlib.sha256(gz).hexdigest(),payload_sha256=h))
        occ.append(z[:d.ns])
        receipts.append(R)
    rows=contrasts(d.ds.ed[1:],occ,d.ds.w)
    k=max(i for i,z in enumerate(rows) if z['q']<=1e-5)
    L=rows[k]
    a=receipts[0]['final_cold']
    if max(abs(u-(v['final_cold']-a)) for u,v in zip(L['delta'],receipts[1:],strict=True))>L['sign_floor']:
        raise ArithmeticError('Cold convention must independently replay on the last whole native cell')
    return dict(rows=rows,within_cell_enclosures=interior_enclosures(d.ds.ed[1:],occ,d.ds.w,rows),
                reference_index=k,reference_convention_q=1e-5,
                scalar_fields=list(NAMES),inputs=ids,
                scope='Exact saved native-cell totals and worst-case nonnegative redistribution within one cell. No new trajectory, observation-level response, or independent spatial discretization error bound.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
