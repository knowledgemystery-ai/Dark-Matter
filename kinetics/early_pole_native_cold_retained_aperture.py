"""Finite archived parent/daughter supports and the cold angular aperture."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from early_pole_endpoint_cold_native_aperture_rectangle import threshold
from early_pole_native_cold_cross_factors import DATE, K, native

M=.05
SPAN=(1e-6,1.5e-5)


def p(E,m):
    return np.sqrt((E-m)*(E+m))


def sep(t,M,m,ep,ed):
    if not (t<0. and ep>M>m and ed>m):
        raise ArithmeticError('Expected the past-directed heavy-to-light angular branch')
    return (M*M-m*m-t*(ep+ed))/(p(ep,M)+p(ed,m))


def aperture(q,s,pole_extent,channel):
    mu,w,m=s['scalar_mass'],s['omega'],s['daughter_mass']
    t=np.hypot(mu,q)-w
    ep,ed=(pole_extent,pole_extent-t) if channel==0 else (pole_extent+t,pole_extent)
    v=sep(t,M,m,ep,ed)
    if not p(ep,M)+p(ed,m)>K+q:
        raise ArithmeticError('The momentum-sum branch must not constrain retained cap')
    c=((K-v)*(K+v)+q*q)/(2*K*q)
    return float((np.clip(c,-1.,1.)+1.)/2.),float(c),float(v-(K-q)),(t,ep,ed)


def summarize(s,cap,channel):
    lo,hi=SPAN
    # The support shrinks in cosine if the gap sep(q)-(K-q) becomes positive.
    # Bound its derivative below independently of a q sample: the moving
    # species has momentum derivative E/p, decreasing in its energy.
    a=aperture(lo,s,cap,channel)
    b=aperture(hi,s,cap,channel)
    j=2 if channel==0 else 1
    moving_mass=s['daughter_mass'] if channel==0 else M
    E=min(a[3][j],b[3][j])
    derivative=1.-(E/p(E,moving_mass))*hi/np.hypot(s['scalar_mass'],hi)
    if derivative<=0. or a[2]>=0.:
        raise ArithmeticError('Independently need a monotone crossing and open lower band')
    cut=brentq(lambda q: aperture(q,s,cap,channel)[2],lo,hi,xtol=2e-19) if b[2]>0. else None
    keys=[lo,4.1e-6,8.2e-6,hi]
    out=[]
    for q in keys:
        f,c,_,(t,ep,ed)=aperture(q,s,cap,channel)
        if -1.<c<1.:
            V=np.sqrt(K*K+q*q-2*K*q*c)
            exact=threshold(M,s['daughter_mass'],t,V)[0]
            if abs(exact/ep-1.)>2e-13:
                raise ArithmeticError('Independent collinear threshold does not close')
        out.append(dict(q=q,unclipped_cosine_limit=c,geometric_fraction_of_cosine_interval=f))
    return dict(first_clipped_q=cut,cut_gap_derivative_lower_bound=derivative,finite_energy_cap=cap,samples=out)


def run(root,neutral):
    A=neutral.read_bytes()
    v=json.loads(A)
    s=[native(root,k) for k in (2,4)]
    if v['contract']['q']!=4.1e-6 or any(v['contract']['raw_native_states'][j]['card_hash']!=u['card_hash'] for j,u in enumerate(s)):
        raise ArithmeticError('Exact completed D2/D4 crossed point for identical held states required')
    rows=[]
    for i,card in enumerate(s):
        for j,sp in enumerate(card['charges']):
            cap=float(sp.E[-1])
            if sp.f[-1]<=0.:
                raise ArithmeticError('Need nonzero occupation inside the last retained cell')
            rows.append(dict(native_chart=(2,4)[i],channel=('forward_parent','inverse_daughter')[j],
                             **summarize(card,cap,j)))
    return dict(rows=rows,q_span=list(SPAN),photon_momentum=K,native_state_hashes=[u['card_hash'] for u in s],
                crossed_point_sha256=hashlib.sha256(A).hexdigest(),
                scope='Geometric fraction of the angular domain on finite numerical charts, not fraction of rate or physical absence of occupation beyond the retained caps.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    for k in ('root','neutral'):
        parser.add_argument(k,type=Path)
    print(json.dumps(run(**vars(parser.parse_args())),indent=2,allow_nan=False))
