"""Bound one-dimensional transport without interpolating stored native panels."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import address, ledger, strip
from early_pole_native_cold_cross_factors import CellChart, DATE, SPECS, load, read
from early_pole_native_cold_mean_from_centered_energy import conditional_q, run as verify


def cdf(bounds,mass,center,m,s):
    conditional_q(bounds,mass,center,m)
    e=np.r_[bounds[:,0],bounds[-1,1]]
    if not np.array_equal(bounds[1:,0],bounds[:-1,1]) or min(s)<e[0] or max(s)>e[-1]:
        raise ArithmeticError('Need exact contiguous admitted support')
    j=np.clip(np.searchsorted(e,s,side='right')-1,0,len(mass)-1)
    low,high=bounds[j].T
    interior=(s>low)&(s<high)
    yl=np.where(s>=high,1.,0.)
    yr=yl.copy()
    def u(q):
        return q*q/(np.hypot(q,m)+m)
    a,t,b=map(u,(low[interior],s[interior],high[interior]))
    mu=center[j[interior]]/mass[j[interior]]
    pad=3e-13*b
    yl[interior]=(t-mu-pad)/(t-a)
    yr[interior]=(b-mu+pad)/(b-t)
    P=np.r_[0.,np.cumsum(mass)]
    x=(P[j]+mass[j]*np.clip(np.array([yl,yr]),0.,1.))/fsum(mass)
    return np.clip(x+np.array([[-3e-14],[3e-14]]),0,1)


def integral(s,a,b):
    # Each true CDF is monotone: on one whole cut interval the value lies
    # between its lower enclosure at the left and upper enclosure at the right.
    lo=np.maximum(0,np.maximum(a[0,:-1]-b[1,1:],b[0,:-1]-a[1,1:]))
    hi=np.maximum(a[1,1:]-b[0,:-1],b[1,1:]-a[0,:-1])
    q=np.diff(s)
    return [fsum(q*lo),fsum(q*hi)]


def run(root,mean):
    raw=mean.read_bytes()
    W=json.loads(raw)
    shape=root/f'early-pole-native-cold-shape-attribution-{DATE}.json'
    if verify(root,shape)!=W:
        raise ArithmeticError('Exact independent centered bounds on raw and weighted mean required')
    L=json.loads((root/f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    B=json.loads((root/f'early-pole-fine-charge-cold-normalized-contrast-{DATE}.json').read_bytes())
    out=[]
    for k,i in ((2,1),(4,2)):
        C,h=read(root/f'early-pole-cell-endpoint-fine-d{k}-x05-state-{DATE}.json.gz')
        if h!=L['cards'][i]['state_sha256_uncompressed']:
            raise ArithmeticError('Exact accepted native state required')
        d=CellChart(load('fine')[0],SPECS['fine'],k,C['reference_x'])
        f=np.asarray(C['endpoint_occupation'][:d.ns])
        J=[]
        for m in (1,2):
            path=root/f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz'
            if unpacked_hash(path)!=L['raw_uncompressed_sha256'][i][m-1]:
                raise ArithmeticError('Both complete paired ledgers must match')
            a=ledger(path); q=np.array([r['bounds'] for r in a])
            z=np.array([r['coeff'] for r in a])[:,(0,2)]
            Z=np.stack([strip(z[:,j],f[address(d.ds.ed,q)]) for j in (0,1)],axis=1)
            J.append((q,[z,Z],a[0]['scalar_mass']))
        out.append(J)
    report=[]
    for k in (0,1):
        r={}
        cuts=[np.unique(np.concatenate([np.linspace(*L['q_band'],n),*[t[k][0].ravel() for t in out]])) for n in (1025,4097)]
        for j,label in enumerate(('forward','inverse')):
            s={}
            for i,kind in enumerate(('weighted','raw')):
                vals=[]
                for q in cuts:
                    p=[cdf(t[k][0],t[k][1][i][:,0,j],t[k][1][i][:,1,j],t[k][2],q) for t in out]
                    vals.append(integral(q,*p))
                l,a=vals
                if a[0]<l[0]-2e-16 or a[1]>l[1]+2e-16:
                    raise ArithmeticError('Nested cut refinement must not loosen monotone bounds')
                mean_bound=W['by_mesh'][k][label][f'{kind}_native_mean_difference']
                lower=max(0.,mean_bound[0],-mean_bound[1],a[0])
                if lower>a[1]:
                    raise ArithmeticError('Independent mean and CDF enclosures disagree')
                s[kind]=dict(cdf_distance_bounds=a,coarser_cut_distance_bounds=l,distance_with_independent_mean_lower=[lower,a[1]])
                if kind=='weighted':
                    old=B['pairs'][-1]['finer_mesh'][f'{label}_number']
                    s[kind]['previous_support_only_distance_bounds']=[old['minimum_compatible_W1_q_bounds'][0],old['maximum_compatible_W1_q']]
            r[label]=s
        report.append(r)
    err=max(abs(np.array(report[0][j][k]['distance_with_independent_mean_lower'])-
                report[1][j][k]['distance_with_independent_mean_lower']).max() for j in report[0] for k in report[0][j])
    if err>3e-13:
        raise ArithmeticError('Independent retained mesh transport bounds must agree')
    return dict(by_mesh=report,additional_uniform_cut_counts=[1025,4097],max_two_mesh_distance_bound_defect=err,
                centered_mean_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Conservative W1 bounds for complete held native raw and weighted forward/inverse q measures, not a physical profile or exact optimizer.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('root','mean'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
