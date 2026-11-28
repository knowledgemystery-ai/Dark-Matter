"""Recompute the local diagnostic on held states using a second thermodynamically consistent EOS interpolation."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_local_entropy_metric_rates import measure
from sm_dark_entropy_exchange_reference import C, StandardModel


def min_piece(c,h):
    p=np.trim_zeros(c,'f')
    if not len(p):
        return 0.
    roots=np.roots(np.polyder(p)) if len(p)>2 else []
    cand=[0.,h,*[z.real for z in roots if abs(z.imag)<2e-13 and 0.<z.real<h]]
    return float(np.min(np.polyval(p,cand)))


class SplineEOS(StandardModel):
    def __init__(self, table):
        super().__init__(table)
        u=self.gs.x
        z=self.gs(u)
        self.gs=CubicSpline(u,z,extrapolate=False)
        self.slope=self.gs.derivative()
        ent,dens=[],[]
        for i,h in enumerate(np.diff(u)):
            c=self.gs.c[:,i]
            ent.append(min_piece(c,h))
            g=3*c+np.r_[0.,np.polyder(c)]
            dens.append(min_piece(g,h))
        self.positivity=(min(ent),min(dens))
        if min(self.positivity)<=0:
            raise ArithmeticError('Independent spline must retain positive entropy and heat capacity on all intervals')
        a=self.press.sol(u[0])[0]
        self.press=solve_ivp(lambda x,_: [4*C/3*self.gs(x)*np.exp(4*x)],(u[0],u[-1]),[a],
                            rtol=3e-11,atol=min(1e-19,1e-10*abs(a)),max_step=.006,dense_output=True)
        if not self.press.success:
            raise ArithmeticError('Need thermodynamically integrated alternate EOS pressure')


def run(eos,early,late,local):
    js=[p.read_bytes() for p in (early,late,local)]
    h=[hashlib.sha256(p).hexdigest() for p in js]
    A,B,W=map(json.loads,js)
    he=hashlib.sha256(eos.read_bytes()).hexdigest()
    if W['input_sha256'][:2]!=h[:2] or not A['table_input_sha256']==B['table_input_sha256']==W['table_input_sha256']==he:
        raise ArithmeticError('Exact held trajectory states and base EOS evidence required')
    u,v=[next(z['points'] for z in y['rows'] if z['operator']=='combined') for y in (A,B)]
    pts=u+v[1:]
    if u[-1]!=v[0] or [p[0] for p in pts]!=[r['x'] for r in W['rows']]:
        raise ArithmeticError('Exact archived diagnostic samples required')
    base=StandardModel(eos)
    new=SplineEOS(eos)
    P=CombinedExpansion(new,PositivePortal((17,44,32)),168)
    results=[]
    for p,row in zip(pts,W['rows'],strict=True):
        old=np.array(base.state(p[1]))
        now=np.array(new.state(p[1]))
        m=measure(P,p,2e-4,1e-5)
        a=row['signed_squared_norm_rate'][-1]
        b=m['signed_squared_norm_rate'][-1]
        results.append(dict(x=p[0],visible_T=p[1],original_max_rate=a,alternate_max_rate=b,
                            original_local_fixed_total_max=row['instantaneous_fixed_total_plane_rate'][-1],
                            alternate_local_fixed_total_max=m['instantaneous_fixed_total_plane_rate'][-1],
                            visible_thermodynamic_relative_change=(now/old-1.).tolist(),
                            alternate_motion_defect=m['density_motion_defect'],
                            alternate_rate_replay_defect=m['independent_direction_rate_defect']))
    return dict(rows=results,positive_entropy_and_capacity_polynomial_lower_values=new.positivity,
                visible_fields=['rho','pressure','entropy_density','heat_capacity'],
                input_sha256=h,table_input_sha256=he,
                scope='Local EOS sensitivity at the same held T and dark physical state of the original trajectory, not a re-evolution in the alternate EOS or a physical uncertainty interval.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','early','late','local'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
