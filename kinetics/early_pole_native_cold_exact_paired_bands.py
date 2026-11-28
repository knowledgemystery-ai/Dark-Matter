"""Compare held D2/D4 signs and neutral intervals on exactly matched bands."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

from early_pole_native_d2_exact_partition import DATE,finish
from early_pole_native_d4_exact_partition import run as d4


def run(root):
    d2_path=root/f'early-pole-native-d2-exact-partition-{DATE}.json'
    d4_path=root/f'early-pole-native-d4-exact-partition-{DATE}.json'
    data=[f.read_bytes() for f in (d2_path,d4_path)]
    X,Y=map(json.loads,data)
    q=X['edges']
    n=X['held_occupation']
    if (q!=Y['edges'] or n!=Y['held_occupation'] or
        finish(root,root/f'd2-exact-partition-patches-{DATE}')!=X or d4(root)!=Y):
        raise ArithmeticError('Both exact whole partitions and held state must replay')
    e=[list(v) for v in zip(q[:-1],q[1:],strict=True)]
    rows=[]
    limits=[]; lo=[]
    for i,(x,y) in enumerate(zip(X['bands'],Y['bands'],strict=True)):
        q0,q1=[a['neutral_box'] for a in (x,y)]
        ss=[a['sign'] for a in (x,y)]
        if ss!=[-1,1 if i<3 else -1] or not q1[1]<q0[0]:
            raise ArithmeticError('Distinct held sign and neutral ordering needed')
        for chart,b,sign in zip((2,4),(q0,q1),ss,strict=True):
            if sign>0:
                lo.append((b[1],i,chart))
            else:
                limits.append((b[0],i,chart))
        rows.append(dict(q=e[i],neutral_D2=q0,neutral_D4=q1,
                         D2_minus_D4_neutral_lower=q0[0]-q1[1],
                         net_box_D2=x['net_box'],net_box_D4=y['net_box'],held_signs_D2_D4=ss))
    a,b=max(lo),min(limits)
    if not a[0]<n<b[0]:
        raise ArithmeticError('Need a nonempty common held corridor')
    N=[v['signed_net_by_band'] for v in X['mesh_runs']]
    return dict(rows=rows,held_occupation=n,common_open_corridor=[a[0],b[0]],
                limiting_charts=[a[2],b[2]],limiting_spans=[e[a[1]],e[b[1]]],
                relative_margin_below_above=[(n-a[0])/n,(b[0]-n)/n],
                D2_central_share_of_loss=[u[1]/fsum(u) for u in N],
                input_sha256=[hashlib.sha256(s).hexdigest() for s in data],
                scope='Two independently evolved held charts with exactly matching retained momentum bands and observed two-mesh boxes; not a spatial-limit or photon-history calculation.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    print(json.dumps(run(p.parse_args().root),indent=2,allow_nan=False))
