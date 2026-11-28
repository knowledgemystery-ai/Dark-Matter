"""Pressure contrasts on the five matched terminal scalar states at held mass."""
import argparse
import gzip
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss

from early_pole_endpoint_charge_symmetric_null_terminal import setup
from early_pole_symmetric_charge_cutoff_audit import DATE, run as audit


def kernel(q,m):
    return q*q/np.hypot(q,m)


def cell_average(a,b,m,order):
    t,w=leggauss(order)
    q=(a[:,None]+b[:,None])/2+(b-a)[:,None]*t/2
    return (kernel(q,m)*q*q)@w*(1.5/(a*a+a*b+b*b))


def jensen_loss(a,b,node,m):
    """Difference between a concave value at node energy and its edge chord."""
    A,B,E=[np.hypot(q,m) for q in (a,b,node)]
    # Rationalize differences even when both nodes are near the held mass.
    left=(node-a)*(node+a)/(E+A)
    right=(b-node)*(b+node)/(B+E)
    return m*m*left*right/(A*B*E)


def bounds(a,b,node,m,mass):
    V=np.asarray(mass)
    if (V.shape!=(5,len(a)) or not np.isfinite(V).all() or min(V.ravel())<0.
        or not np.isfinite(m) or m<=0. or np.any(a<=0.)
        or np.any(a>=node) or np.any(node>=b)):
        raise ArithmeticError('Positive mass, ordered cell edges, interior nodes and nonnegative states required')
    K=[kernel(q,m) for q in (a,node,b)]
    D=V[1:]-V[0]
    lower=np.minimum(D*K[0],D*K[2])
    upper=np.maximum(D*K[0],D*K[2])
    loss=jensen_loss(a,b,node,m)
    out=[]
    for i,d in enumerate(D,1):
        nodal=fsum(d*K[1])
        out.append(dict(shared_shape=[fsum(lower[i-1]),fsum(upper[i-1])],
                        independent_counts_only=[fsum(d*K[0]-V[0]*(K[2]-K[0])),
                                                 fsum(d*K[0]+V[i]*(K[2]-K[0]))],
                        independent_fixed_mean_energies=[nodal-fsum(V[i]*loss),nodal+fsum(V[0]*loss)],
                        nodal=nodal))
    return out


def run(root,prior):
    raw=prior.read_bytes()
    checked=audit(root)
    if checked!=json.loads(gzip.decompress(raw)):
        raise ArithmeticError('Require exact replay of the original five cutoff-audit fixtures')
    d,*_=setup(root)
    edges=d.ds.ed
    a,b=edges[:-1],edges[1:]
    if np.any(b<=a) or len(a)!=len(checked['rows']):
        raise ArithmeticError('Same complete native chart required')
    node=np.sqrt(a*b)
    V=[]; M=[]
    for z in checked['inputs']:
        stem=root/f"early-pole-{z['name']}-x05"
        rec=Path(f'{stem}-receipt-{DATE}.json').read_bytes()
        gz=Path(f'{stem}-checkpoint-{DATE}.json.gz').read_bytes()
        if hashlib.sha256(rec).hexdigest()!=z['receipt_sha256'] or hashlib.sha256(gz).hexdigest()!=z['compressed_checkpoint_sha256']:
            raise ArithmeticError('Terminal fixture has changed since cutoff replay')
        V.append(np.asarray(json.loads(gzip.decompress(gz))['occupation'][:d.ns])*d.ds.w)
        M.append(json.loads(rec)['final_pole'])
    m=M[0]
    V=np.asarray(V)
    Q,S=[cell_average(a,b,m,n) for n in (16,32)]
    E=kernel(node,m)
    ref=fsum(V[0]*S)
    rnode=fsum(V[0]*E)
    diff=(V[1:]-V[0])*S
    rows=[]
    for j,(z,diff_bin) in enumerate(zip(bounds(a,b,node,m,V),diff,strict=True),1):
        dN=fsum(V[j]-V[0])/fsum(V[0])
        dP=fsum(diff_bin)/ref
        mm=M[j]; u,v=[np.hypot(node,R) for R in (m,mm)]
        dk=node*node*((m-mm)*(m+mm))/(u*v*(u+v))
        i=int(np.argmin(diff_bin))
        rows.append(dict(arm=checked['scalar_fields'][j],count_relative=dN,pressure_relative=dP,
                         pressure_per_particle_relative=(dP-dN)/(1+dN),
                         quadrature_relative_discrepancy=abs(fsum((V[j]-V[0])*Q)/fsum(V[0]*Q)-dP),
                         common_shape_bounds=[t/ref for t in z['shared_shape']],
                         independent_count_only_bounds=[t/rnode for t in z['independent_counts_only']],
                         independent_count_and_node_mean_energy_bounds=[t/rnode for t in z['independent_fixed_mean_energies']],
                         nodal_pressure_relative=z['nodal']/rnode,
                         nodal_same_occupation_mass_only_effect=fsum(V[j]*dk)/rnode,
                         cold_q_at_most_1e_minus_5_number_weighted_pressure_contrast=fsum(diff_bin[b<=1e-5])/ref,
                         most_negative_cell=dict(q=[a[i],b[i]],pressure_relative_contribution=diff_bin[i]/ref)))
    if max(r['quadrature_relative_discrepancy'] for r in rows)>1e-12:
        raise ArithmeticError('Independent quadrature settings fail to reproduce pressure contrasts')
    f=np.array([r['pressure_relative'] for r in rows]).reshape(2,2)
    odd=(f[:,1]-f[:,0])/2
    even=f.mean(axis=1)
    return dict(rows=rows,held_measurement_mass=m,nominal_histogram_pressure_without_common_third=ref,
                nominal_nodal_pressure_without_common_third=rnode,
                full_and_half_odd_pressure_relative=odd.tolist(),full_and_half_even_pressure_relative=even.tolist(),
                half_odd_over_linear_full=odd[1]/(.5*odd[0]),half_even_over_quadratic_full=even[1]/(.25*even[0]),
                quadrature_orders=[16,32],cutoff_receipt_sha256=hashlib.sha256(raw).hexdigest(),
                scope='Fixed-mass momentum-weighted diagnostic at the exact saved terminal states. Arbitrary independent shapes need held per-cell node mean energies to certify pressure signs; no physical pressure or spatial convergence claim.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('root',type=Path)
    p.add_argument('prior',type=Path)
    args=p.parse_args()
    print(json.dumps(run(args.root,args.prior),indent=2,allow_nan=False))
