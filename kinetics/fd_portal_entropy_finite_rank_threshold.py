"""Threshold for all fixed carried directions under an integrated rate-term swap."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_entropy_finite_gain_budget import terms
from fd_portal_entropy_hot_start_finite_crossing import DATE,dense
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_entropy_short_tangent_windows import gain_difference
from sm_dark_entropy_exchange_reference import StandardModel


def sym(x):
    return (x+x.T)/2.


def threshold(B,R):
    v=eigh(R,-B,eigvals_only=True)
    if v[-1]<=0 or v[0]<-1e-12*v[-1]:
        raise ArithmeticError('Requires negative reference and a positive integrated metric term')
    return 1/v[-1]


def run(eos,root):
    names=('direct-expansion-096-120','entropy-short-tangent-windows','entropy-finite-gain-budget')
    raw=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    A,T,B=map(json.loads,raw)
    ids=[hashlib.sha256(s).hexdigest() for s in raw]
    he=hashlib.sha256(eos.read_bytes()).hexdigest()
    if T['input_sha256'][0]!=ids[0] or B['input_sha256'][:2]!=ids[:2] or any(z['table_input_sha256']!=he for z in (A,T,B)):
        raise ArithmeticError('Exact completed budget and original tangent required')
    a,b=B['window']
    K=[a,T['windows'][0]['x'][1],T['split_at_knot'],T['windows'][2]['x'][0],b]
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    init=next(p for q in A['rows'] if q['operator']=='combined' for p in q['points'] if p[0]==a)
    path=local_path(P,init)
    Q=np.asarray(T['rebased_fractional_columns'])
    rows=[]
    for n,e,k in ((4,4e-5,3),(6,2e-5,6)):
        C=dense(P,path,K,Q,n,e)
        def G(t):
            W=C(t)
            return sym(W.T@matrix(P,path(t))[0]@W)
        X,Y=G(a),G(b)
        H=np.zeros((2,2)); R=np.zeros((2,2))
        x,w=leggauss(k)
        for s,t in zip(K[:-1],K[1:],strict=True):
            h=(t-s)/2
            for y,z in zip((s+t)/2+h*x,h*w,strict=True):
                M=C(y)
                V=terms(P,path(y),e)
                R+=z*sym(M.T@V[2]@M)
                H+=z*sym(M.T@sum(V[j] for j in (0,1,3,4))@M)
        closure=float(np.max(abs(eigh(sym(H+R-Y+X),X,eigvals_only=True))))
        u=threshold(H,R)
        a0,b0=[gain_difference(X,X+H+c*R) for c in (0.,1.)]
        if (closure>3e-10 or not 0.<u<1. or a0[-1]>=0. or b0[-1]<=0.
            or abs(gain_difference(X,X+H+u*R)[-1])>5e-13):
            raise ArithmeticError('Integrated all-direction relative Gram form failed closure')
        rows.append(dict(critical_multiplier_for_held_integrated_log_slope_term=u,
                         largest_gain_minus_one_at_no_term=a0[-1],largest_gain_minus_one_at_observed_term=b0[-1],
                         direct_metric_endpoint_closure=closure,
                         integrated_reference_matrix=H.tolist(),integrated_slope_matrix=R.tolist()))
    return dict(window=B['window'],rows=rows,table_input_sha256=he,input_sha256=ids,
                scope='Counterfactual integrated rate form on one fixed numerical tangent and background, not an independently evolved entropy metric or EOS.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
