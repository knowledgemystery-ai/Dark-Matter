"""Complete cosine integration at fixed timelike pole momentum."""
import argparse
import json

import numpy as np
from scipy.special import roots_legendre

from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q
from finite_mass_qed_timelike_pole_affine_source_reference import setup
from finite_mass_qed_timelike_pole_batched_absorption_reference import load, pole, radial, spacelike_parent_minimum


def angular_knots(ch,k,w,es,qscalar=Q,*,allow_sum_boundaries=False,parent_energy_edges=None):
    if qscalar<0:
        raise ValueError('A scalar momentum magnitude cannot be negative')
    if qscalar==0:
        return np.array([-1.,1.]),0.
    parent,daughter=ch
    m0,m1=parent.m,daughter.m
    t=es-w
    extra=[] if parent_energy_edges is None else parent_energy_edges
    e=np.unique(np.r_[parent.E,daughter.E+t,extra])
    e=e[(e>m0)&(e-t>m1)]
    p=np.sqrt((e-m0)*(e+m0))
    b=np.sqrt((e-t-m1)*(e-t+m1))
    sep=(m0*m0-m1*m1-2*e*t+t*t)/(b+p)
    c=((k-sep)*(k+sep)+qscalar*qscalar)/(2*k*qscalar)
    upper=((k-p-b)*(k+p+b)+qscalar*qscalar)/(2*k*qscalar)
    sum_arm=(upper>-1)&(upper<1)
    if np.any(sum_arm) and not allow_sum_boundaries:
        raise ArithmeticError('Unexpected sum-of-momenta angular boundary in this soft branch')
    ce=c[(c>-1)&(c<1)]
    ee=e[(c>-1)&(c<1)]
    if allow_sum_boundaries:
        ce=np.r_[ce,upper[sum_arm]]
        ee=np.r_[ee,e[sum_arm]]
    err=[]
    for u,E in zip(ce,ee,strict=True):
        if qscalar==Q:
            _,_,_,edge,_=setup(m0,m1,k,w,es,u)
        else:
            V=np.hypot(qscalar*np.sqrt((1-u)*(1+u)),qscalar*u-k)
            s=(t-V)*(t+V)
            A=m0*m0+s-m1*m1
            edge=(spacelike_parent_minimum(m0,m1,t,V,s) if t>=0 else
                  (-A*t+V*np.sqrt(A*A-4*m0*m0*s))/(-2*s))
        err.append(abs(edge/E-1))
    defect=max(err,default=0.)
    if defect>3e-12:
        raise ArithmeticError(('Exact angular energy support knot',k,w,defect))
    past=((k+t)*(k-t)+qscalar*qscalar)/(2*k*qscalar)
    return np.unique(np.clip(np.r_[-1.,1.,ce,past],-1.,1.)),defect


def angle(ch,k,w,R,es,name,order,inner,qscalar=Q,*,allow_future_spacelike=False,occupation_callbacks=None,mirror_azimuth=False,parent_energy_edges=None,transverse_vertex_kernel=None):
    bins,defect=angular_knots(ch,k,w,es,qscalar,allow_sum_boundaries=allow_future_spacelike,
                             parent_energy_edges=parent_energy_edges)
    x,v=roots_legendre(order)
    counts=np.zeros(3,dtype=int)
    A=np.zeros(3) if parent_energy_edges is None else np.zeros((len(parent_energy_edges)-1,3))
    for a,b in zip(bins[:-1],bins[1:],strict=True):
        for u,z in zip((a+b)/2+(b-a)*x/2,(b-a)*v/4,strict=True):
            kwargs={} if occupation_callbacks is None else {'occupation_callbacks':occupation_callbacks}
            if mirror_azimuth:
                kwargs['mirror_azimuth']=True
            if parent_energy_edges is not None:
                kwargs['parent_energy_edges']=parent_energy_edges
            if transverse_vertex_kernel is not None:
                kwargs['transverse_vertex_kernel']=transverse_vertex_kernel
            C,_,_=radial(ch,k,w,R,es,u,name,*inner,qscalar=qscalar,
                         allow_future_spacelike=allow_future_spacelike,**kwargs)
            tot=C if parent_energy_edges is None else C.sum(axis=0)
            counts[0]+=tot[0]>0
            counts[1]+=tot[1]>0
            counts[2]+=(tot[0]>0)&(tot[1]==0)
            A+=z*C
    return A,counts,len(bins)-1,defect


def tabulate(g):
    d,ch,es=load(g)
    n=d['endpoint_scalar_occupations'][np.searchsorted(d['scalar_momentum_edges'],Q)-1]
    output=[]
    near={'fine':(.001395,.001405),'third':(.001415,.001425)}
    targets=[('T',.00135),*[('T',x) for x in near[g]],('T',.0015),('T',.00256),('L',.0012),('L',.00256)]
    for name,k in targets:
        w,R=pole(ch,name,k)
        A,C,N,J=angle(ch,k,w,R,es,name,3,(3,8))
        B,D,_,K=angle(ch,k,w,R,es,name,6,(6,12))
        if np.any(np.sign(A)!=np.sign(B)):
            raise ArithmeticError(('Unsupported angular directional claim',g,name,k,A,B))
        err=float(max(abs(A/B-1.))) if np.all(B!=0.) else float(max(abs(A-B))/max(1e-50,max(abs(B))))
        if err>3e-7:
            raise ArithmeticError(('Joint angular estimate did not converge',g,name,k,err))
        center,_,_=radial(ch,k,w,R,es,0.,name,6,12)
        output.append(dict(grid=g,mode=name,k=k,omega=w,accepted_inner_mesh=[6,12],outer_order=6,
                           complete_cosine_average=B.tolist(),cosine_zero=center.tolist(),
                           max_joint_mesh_relative_defect=err,angular_panels=N,accepted_gain_nodes=int(D[0]),
                           accepted_inverse_nodes=int(D[1]),accepted_forward_only_nodes=int(D[2]),
                           max_analytic_angle_knot_defect=max(J,K),
                           conditional_balancing_population=B[1]/(B[2]+B[0]/n)
                            if (B[1]>0. and B[2]+B[0]/n>0.) else None))
    if not (output[0]['accepted_gain_nodes']==0 and output[1]['accepted_forward_only_nodes']>0
            and all(q['accepted_inverse_nodes']==0 for q in output[1:3])):
        raise ArithmeticError('Past-null or one-sided support fixture failed')
    return output


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('grid',choices=('fine','third'))
    print(json.dumps(tabulate(p.parse_args().grid),indent=2,allow_nan=False))
