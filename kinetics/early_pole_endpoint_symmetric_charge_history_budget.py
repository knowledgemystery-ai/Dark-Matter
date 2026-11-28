"""Audit true charge-pair odd and even responses at every common saved time."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def classify(lo,ref,hi):
    if not (len(lo)==len(ref)==len(hi)>1):
        raise ArithmeticError('Need all matched histories')
    x=np.array([[r['x'] for r in q] for q in (lo,ref,hi)])
    Y=np.array([[v['scalar_yield'] for v in q] for q in (lo,ref,hi)])
    if (not np.array_equal(x[0],x[1]) or not np.array_equal(x[2],x[1]) or np.any(np.diff(x[0])<=0.)
        or np.any(Y<=0.) or not np.all(np.isfinite(Y))):
        raise ArithmeticError('Same exact increasing checkpoint times and positive yields required')
    a,b=Y[0]-Y[1],Y[2]-Y[1]
    odd,even=(b-a)*.5,(a+b)*.5
    good_gap=np.abs(odd)>1e-12
    flips=[x[1,i:i+2].tolist() for i in range(len(odd)-1) if good_gap[i] and good_gap[i+1] and odd[i]*odd[i+1]<0.]
    return dict(rows=[dict(x=float(t),low_minus_nominal=float(l),high_minus_nominal=float(h),
                           symmetric_odd_yield=float(u),symmetric_even_yield=float(v))
                      for t,l,h,u,v in zip(x[1],a,b,odd,even,strict=True)],
                brackets_with_opposite_saved_odd_signs=flips,
                max_continuity_residual=max(abs(q['continuity_residual']) for row in (lo,ref,hi) for q in row),
                max_absolute_pair_count_residual=max(abs(q['pair_number_error']) for row in (lo,ref,hi) for q in row))


def run(summary,lower,nominal,upper):
    rows=[p.read_bytes() for p in (lower,nominal,upper)]
    Q=[hashlib.sha256(x).hexdigest() for x in rows]
    p=summary.read_bytes()
    M=json.loads(p)
    A,B,C=map(json.loads,rows)
    if ([Q[0],Q[2]]!=M['raw_receipt_sha256'] or A['direction']!=-1 or C['direction']!=1
        or A['setup']!=C['setup'] or any(B[k]!=v for k,v in A['setup'].items())
        or len(B['checks'])!=8):
        raise ArithmeticError('Must use both exact symmetric and nominal held-remap histories')
    T=classify(A['checks'],B['checks'],C['checks'])
    scale=A['final_whole_scalar_count']/A['checks'][-1]['scalar_yield']
    E=M['results']['whole_scalar_count']
    for k,m in (('finite_odd','symmetric_odd_yield'),('finite_even','symmetric_even_yield')):
        if abs(T['rows'][-1][m]*scale-E[k])>5e-17:
            raise ArithmeticError('Independent terminal whole-count and yield ledgers must agree')
    return dict(**T,input_receipt_sha256=Q,input_summary_sha256=hashlib.sha256(p).hexdigest(),yield_to_count=scale,
                scope='Signs only at saved matched times; this does not count all zeros, resolve time-step effects, or establish physical nonlinearity.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('summary','lower','nominal','upper'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
