"""Check whether finite-horizon entropy gains are robust to the alternate smooth EOS."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import gain, matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff
from fd_portal_direct_comoving_variational_propagator import solve
from fd_portal_entropy_eos_interpolation_check import SplineEOS
from fd_portal_transported_entropy_plane import flow


def run(eos,spline,metric):
    raw=[p.read_bytes() for p in (spline,metric)]
    B,M=map(json.loads,raw)
    h=[hashlib.sha256(b).hexdigest() for b in raw]
    te=hashlib.sha256(eos.read_bytes()).hexdigest()
    if B['input_sha256'][:2]!=M['input_sha256'][:2] or not B['table_input_sha256']==M['table_input_sha256']==te:
        raise ArithmeticError('The re-evolved alternate background and nominal metric must match')
    path=[r['alternate_reintegrated_state'] for r in B['rows']]
    if (path[0][0],path[-1][0])!=(.55,1.2):
        raise ArithmeticError('Need the full held alternate window')
    P=CombinedExpansion(SplineEOS(eos),PositivePortal((17,44,32)),168)
    rates,_,base,err=solve(P,path,8e-12,.012)
    def jac(x):
        q=base(x)
        return diff(lambda u:rates(x,q*u).sum(axis=0),2e-4)/q[None,:]
    times=np.linspace(.55,1.2,129)
    grid=np.array([jac(x) for x in times])
    qry=np.array([.96,1.2])
    U,V=[flow(times[::k],grid[::k],qry,2e-13)/base(qry).T[:,:,None] for k in (2,1)]
    begin,a=matrix(P,path[0])
    ans=[]
    for i,(x,row) in enumerate(zip(qry,M['rows'],strict=True)):
        G,m=matrix(P,next(p for p in path if p[0]==x))
        coarse,alt=[gain(z[i],begin[:2,:2],G) for z in (U,V)]
        refr=np.array(row['output_full_entropy_metric_gain'])
        d=np.max(abs(coarse-alt))
        if d>1e-5:
            raise ArithmeticError(('Alternate numerical entropy gain is still underresolved',x,d))
        ans.append(dict(x=x,original_endpoint_entropy_gains=refr.tolist(),
                        alternate_endpoint_entropy_gains=alt.tolist(),alternate_minus_original=(alt-refr).tolist(),
                        alternate_mesh_gain_difference=d,alternate_endpoint_metric_fd_relative_error=m,
                        alternate_end_fractional_columns=V[i].tolist()))
    return dict(rows=ans,input_sha256=h,table_input_sha256=te,replayed_alternate_background_defect=err,
                initial_metric_fd_relative_error=a,sampled_jacobian_mesh_sizes=[65,129],
                scope='Two different EOS interpolants, each evaluated on its own re-evolved restricted history and endpoint entropy metrics; not EOS uncertainty or physical stability bounds.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','spline','metric'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
