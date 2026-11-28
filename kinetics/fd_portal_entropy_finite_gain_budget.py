"""Close the finite entropy gain along one maximizer with a split rate integral."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

from fd_portal_combined_entropy_split_expansion import CombinedExpansion
from fd_portal_composed_memory_entropy_metric import matrix
from fd_portal_differential_reciprocity_reference import PositivePortal
from fd_portal_direct_comoving_variational_jacobian import diff,full
from fd_portal_entropy_hot_start_finite_crossing import DATE,dense
from fd_portal_entropy_pchip_knot_crossing import local_path
from fd_portal_visible_curvature_rate_attribution import visible
from sm_dark_entropy_exchange_reference import StandardModel

TERMS=('dark_metric_motion','visible_metric_motion_at_cubic_slope','visible_heat_capacity_slope_above_three',
       'expansion_flow','collision_flow')


def dark_form(P,p):
    x,_,eta,lt=p[:4]
    T=np.exp(lt)
    z=P.bath([eta,lt])
    n,r,B,M,C=[z[k] for k in ('n','rho','dn_deta','chemical_energy','cv_fixed_number')]
    J=np.array([[B,B*M/T],[B*M,T*C+B*M*M/T]])
    H=np.diag([1.,1./T])@np.linalg.inv(J)
    out=np.zeros((3,3))
    q=np.array([n,r])
    out[:2,:2]=np.exp(3*x)*(q[:,None]*H*q[None,:])
    return (out+out.T)/2


def terms(P,point,eps):
    p=np.asarray(point[:4])
    x,T,e,_=p[:4]
    state=P.bath(p[2:4])
    r=state['rho']
    v,_,_,cv=P.eos.state(T)
    phi=full(P,p[1:4])
    con=phi(np.ones(3))
    L=diff(phi,eps)-np.array([np.diag(z) for z in con])
    z=con.sum(axis=0)
    n,r,B,D,C=[state[k] for k in ('n','rho','dn_deta','chemical_energy','cv_fixed_number')]
    td=np.exp(p[3])
    J=np.array([[B,B*D/td],[B*D,td*C+B*D*D/td]])
    q0=np.array([n,r,r+v])
    slopes=q0*(z-3.)
    vel=np.r_[1.,(slopes[2]-slopes[1])/cv,np.linalg.solve(J,slopes[:2])]
    # Use an independent within-side density probe. The fixed wider stencil
    # would straddle the numerical EOS knot at the closest Gauss nodes.
    dt=min(1e-6,.2*abs(T-.3595)/abs(vel[1]))
    def rho(w):
        a=P.bath(w[2:4])
        return np.array([a['n'],a['rho'],a['rho']+P.eos.state(w[1])[0]])
    if dt<1e-10 or np.max(abs(((rho(p+dt*vel)-rho(p-dt*vel))/(2*dt)-slopes)/q0))>3e-7:
        raise ArithmeticError('Same-side density direction check failed')
    U=-vel[1]/T
    a=np.array([0.,-r,r+v])
    b=a*(z-3.)
    V=visible(P,p)
    G=dark_form(P,p)+V
    h=1e-6
    j=[dark_form(P,p+k*h*vel) for k in (-2,-1,1,2)]
    dark=(j[0]-8*j[1]+8*j[2]-j[3])/(12*h)
    base=(3.+5.*U)*V+np.exp(3*x)*(np.outer(a,b)+np.outer(b,a))/(T*T*cv)
    t=np.log(T)
    g,g1,g2=P.eos.gs(t),P.eos.slope(t),P.eos.slope.derivative()(t)
    rank=((3*g1+g2)/(3*g+g1))*U*V
    return [dark,base,rank,*[m.T@G+G@m for m in L]]


def integrate(P,path,frame,breaks,c,order,eps):
    t,w=leggauss(order)
    total=[]
    for a,b in zip(breaks[:-1],breaks[1:],strict=True):
        hh=(b-a)/2
        for x,z in zip((a+b)/2+hh*t,hh*w,strict=True):
            v=frame(x)@c
            total.append([z*(v@r@v) for r in terms(P,path(x),eps)])
    return np.array([fsum(v) for v in np.asarray(total).T])


def run(eos,root):
    names=('direct-expansion-096-120','entropy-short-tangent-windows','entropy-short-nonlinear-witness')
    raw=[(root/f'fd-portal-{k}-{DATE}.json').read_bytes() for k in names]
    A,B,C=map(json.loads,raw)
    ds=[hashlib.sha256(z).hexdigest() for z in raw]
    h=hashlib.sha256(eos.read_bytes()).hexdigest()
    if B['input_sha256'][0]!=ds[0] or C['input_sha256']!=ds[:2] or any(z['table_input_sha256']!=h for z in (A,B,C)):
        raise ArithmeticError('Exact original finite tangent and nonlinear evidence required')
    x,t=B['windows'][-1]['x']
    hot=B['windows'][0]['x'][1]; cool=B['windows'][2]['x'][0]
    bk=[x,hot,B['split_at_knot'],cool,t]
    P=CombinedExpansion(StandardModel(eos),PositivePortal((17,44,32)),168)
    seed=next(z for z in next(r['points'] for r in A['rows'] if r['operator']=='combined') if z[0]==x)
    path=local_path(P,seed)
    initial=np.asarray(B['rebased_fractional_columns'])
    data=[]
    for n,eps,k in ((4,4e-5,3),(6,2e-5,6)):
        frame=dense(P,path,bk,initial,n,eps)
        def G(z):
            u=frame(z)
            W=u.T@matrix(P,path(z))[0]@u
            return (W+W.T)/2
        X,Y=G(x),G(t)
        lam,Q=eigh(Y-X,X)
        v=Q[:,-1]
        V=integrate(P,path,frame,bk,v,k,eps)
        gain=lam[-1]/(np.sqrt(1.+lam[-1])+1.)
        delta=abs(fsum(V)-lam[-1])
        if abs(gain-B['windows'][-1]['ordered_gain_minus_one'][-1])>1e-11 or delta>3e-10:
            raise ArithmeticError(('Integrated finite gain does not close',delta))
        data.append(dict(terms=dict(zip(TERMS,V.tolist(),strict=True)),sum_of_integrals=fsum(V),
                         endpoint_squared_gain_minus_one=lam[-1],endpoint_gain_minus_one=gain,closure_defect=delta))
    return dict(window=[x,t],rows=data,table_input_sha256=h,input_sha256=ds,
                scope='Budget integrated along the exact final-gain maximizing carried original tangent, not an allocation of physical causes or global dynamics.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('eos','root'):
        p.add_argument(k,type=Path)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
