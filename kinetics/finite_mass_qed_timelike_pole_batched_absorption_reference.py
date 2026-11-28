"""Batched fixed-aperture spinors with exact two-arm local histogram support."""
import gzip
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq
from scipy.special import roots_legendre

from early_pole_hard_yukawa_photoproduction_reference import I, dot, norm_spinors, slash
from finite_mass_qed_frozen_transverse_plasmon_gate_reference import Q, Species
from finite_mass_qed_longitudinal_static_spectral_reference import DensityResponse
from finite_mass_qed_timelike_pole_affine_source_reference import run_one


def pole(ch,name,k):
    if name=='L':
        return DensityResponse(ch).pole(k)
    def eq(w):
        return (w-k)*(w+k)-sum(s.field(k/w)[0] for s in ch)
    w=brentq(eq,np.nextafter(k,np.inf),k+.006,xtol=5e-17)
    return w,1/(2*w+k/w**2*sum(s.field(k/w)[1] for s in ch))


def reduce_vertex(p,b,K,m0,m1,e):
    # All momenta are contravariant; spinors keep the full angular layout.
    w,k=K[0],K[3]
    mu2=(w-k)*(w+k)
    v=[]
    for o,m in ((p,m0),(b,m1)):
        spatial=o[...,1:]
        z=spatial[...,2]
        q=np.linalg.norm(spatial,axis=-1)
        turn=q-z
        mask=z>0
        turn[mask]=np.sum(spatial[mask,:2]**2,axis=-1)/(q[mask]+z[mask])
        v.append((m*m*w*w+mu2*q*q)/(np.hypot(m,q)*w+q*k)+k*turn)
    left=(2*dot(p,e)[...,None,None]*I+slash(K)@slash(e))/(mu2+2*v[0])[...,None,None]
    right=(2*dot(b,e)[...,None,None]*I-slash(e)@slash(K))/(mu2-2*v[1])[...,None,None]
    return left+right


def full_vertex(p,b,K,m0,m1,e):
    s=dot(K,K)
    return ((slash(p+K)+m0*I)@slash(e)/(s+2*dot(p,K))[...,None,None]
             +slash(e)@(slash(b-K)+m1*I)/(s-2*dot(b,K))[...,None,None])


def spacelike_parent_minimum(M,m,t,v,s):
    # The upper quadratic root is the positive-energy threshold. For
    # positive A*t, rationalize the subtractive numerator near s=0.
    A=M*M+s-m*m
    z=v*np.sqrt(A*A-4*M*M*s)
    if A*t>0:
        return (A*A+4*M*M*v*v)/(2*(z+A*t))
    return (-A*t+z)/(-2*s)


def radial(ch,k,w,R,es,c,name,n,az,audit=False,extra_energy_breaks=None,qscalar=Q,
           *,allow_future_spacelike=False,occupation_callbacks=None,mirror_azimuth=False,parent_energy_edges=None,transverse_vertex_kernel=None):
    P,F=ch
    M,m=P.m,F.m
    T=np.array([es-w,qscalar*np.sqrt((1-c)*(1+c)),0.,qscalar*c-k])
    v=np.linalg.norm(T[1:])
    t=T[0]
    hi=max(P.E[-1],F.E[-1]+t)
    edges=None if parent_energy_edges is None else np.asarray(parent_energy_edges,dtype=float)
    if edges is not None and (edges.ndim!=1 or len(edges)<2 or not np.all(np.isfinite(edges)) or
            edges[0]>M or edges[-1]<hi or np.any(np.diff(edges)<=0)):
        raise ValueError('Parent energy bands must partition the entire retained support')
    empty=lambda: np.zeros(3) if edges is None else np.zeros((len(edges)-1,3))
    s=(t-v)*(t+v)
    if t>=0 and not allow_future_spacelike:
        raise ValueError('This branch requires a past-directed transfer energy')
    if allow_future_spacelike:
        # K heavier than the external scalar implies that a timelike
        # q-K is past-directed, and therefore cannot equal P-b for M>m.
        if (w-k)*(w+k)<=es*es-qscalar*qscalar or not M>m:
            raise ValueError('Timelike-transfer exclusion does not hold here')
    if s>=0:
        return empty(),False,0.
    A=M*M+s-m*m
    lo=(spacelike_parent_minimum(M,m,t,v,s) if t>=0 else
        (-A*t+v*np.sqrt(A*A-4*M*M*s))/(-2*s))
    if not lo<hi:
        return empty(),False,0.
    extras=[] if extra_energy_breaks is None else extra_energy_breaks(M,m,T,s,A,v,lo,hi)
    if edges is not None:
        extras=np.r_[extras,edges]
    ee=np.unique(np.clip(np.r_[lo,hi,P.E,F.E+t,extras],lo,hi))
    a,b=ee[:-1],ee[1:]
    c0=(a+b)/2
    ff=np.r_[0.,P.f,0.][np.searchsorted(P.E,c0)]
    gg=np.r_[0.,F.f,0.][np.searchsorted(F.E,c0-t)]
    keep=(ff>0)|(gg>0)
    if not any(keep):
        return empty(),False,0.
    a,b,ff,gg,c0=(x[keep] for x in (a,b,ff,gg,c0))
    x,h=roots_legendre(n)
    e=(c0[:,None]+(b-a)[:,None]*x/2).ravel()
    z=(e*t-A/2)/v
    p2=(e-M)*(e+M)
    radicand=p2-z*z
    if np.min(radicand/(1+p2))<-3e-14:
        raise ArithmeticError('No valid radial angle')
    p=np.sqrt(np.maximum(radicand,0.))
    axis=T[1:]/v
    cross=np.array([axis[2],0.,-axis[0]])
    if mirror_azimuth and (az<2 or az%2):
        raise ValueError('Mirror reduction requires a positive even full azimuth order')
    # With k and q in the x-z plane, reflection in y is complex conjugation
    # of spinors and slashed real vectors in this basis. For e_y there is
    # one additional overall sign. Every squared amplitude is thus even.
    aphi=(np.arange(az//2 if mirror_azimuth else az)+.5)*2*np.pi/az
    transverse=np.cos(aphi)[:,None]*cross+np.sin(aphi)[:,None]*np.array([0.,1.,0.])
    mom=z[:,None,None]*axis+p[:,None,None]*transverse[None,:,:]
    P0=np.concatenate([np.broadcast_to(e[:,None,None],(*mom.shape[:2],1)),mom],axis=-1)
    B0=P0-T
    K=np.array([w,0.,0.,k])
    pol=np.eye(4)[1:3] if name=='T' else np.eye(4)[:1]
    integral=0.
    max_def=0.
    if transverse_vertex_kernel is not None:
        if name!='T' or audit:
            raise ValueError('The optional transverse kernel is only valid for the unaudited transverse pair')
        values=transverse_vertex_kernel(P0,B0,K,M,m)
        if values.shape!=P0.shape[:-1] or not np.all(np.isfinite(values)) or np.any(values<0.):
            raise ArithmeticError('Invalid optional transverse spin sum')
        integral=values.mean(axis=-1)
    else:
        for u in pol:
            A=norm_spinors(P0,B0,M,m,reduce_vertex(P0,B0,K,M,m,u))
            if audit:
                B=norm_spinors(P0,B0,M,m,full_vertex(P0,B0,K,M,m,u))
                max_def=max(max_def,float(np.sum(abs(A-B))/np.sum(B)))
            integral=integral+A.mean(axis=-1)
    if occupation_callbacks is None:
        e_int=((b-a)/2)*((integral.reshape(-1,n))@h)*R/v
        rows=np.asarray([ff*(1-gg),gg*(1-ff),ff-gg])
        w0=e_int
        node_bin=c0
    else:
        # A number/energy-conserving continuum reconstruction carries
        # occupation variation *inside* each of the same support panels.
        fp,fb=[fn(z) for fn,z in zip(occupation_callbacks,(e,e-t),strict=True)]
        if fp.shape!=e.shape or fb.shape!=e.shape or np.min([fp.min(),fb.min()])<0 or np.max([fp.max(),fb.max()])>1:
            raise ArithmeticError('Continuum charged occupation left its Pauli domain')
        w0=(((b-a)/2)[:,None]*h).ravel()*integral*R/v
        rows=np.asarray([fp*(1-fb),fb*(1-fp),fp-fb])
        node_bin=e
    if edges is None:
        C=rows@w0
    else:
        kbin=np.searchsorted(edges,node_bin,side='right')-1
        if min(kbin)<0 or max(kbin)>=len(edges)-1:
            raise ArithmeticError('An energy sample left its declared retained band')
        C=np.column_stack([np.bincount(kbin,v*w0,minlength=len(edges)-1) for v in rows])
    D=np.atleast_2d(C)
    if np.min(D[:,:2])<0 or np.any(abs(D[:,0]-D[:,1]-D[:,2])>5e-13*np.max(D[:,:2],axis=1)):
        raise ArithmeticError('Impossible signed population support')
    return C,True,max_def


def load(g):
    folder=Path(__file__).parent/'results'
    with gzip.open(folder/f'early-pole-spectrum-{g}-m06-2026-11-26.json.gz') as s:d=json.load(s)
    M=d['x_end']
    ch=[Species(m,d[f'endpoint_{j}_momentum_edges'],d[f'endpoint_{j}_occupations'])
        for m,j in ((M,'parent'),(M*d['mass_ratio'],'fermion'))]
    es=np.hypot(Q,d['checkpoints'][-1]['scalar_pole_over_T'])
    return d,ch,es


def fixture():
    out=[]
    for g in ('fine','third'):
        d,ch,es=load(g)
        for name,k in (('T',.0015),('T',.00256),('L',.0012),('L',.00256)):
            w,R=pole(ch,name,k)
            for c in (0.,-.77,.82):
                A,ok,trace=radial(ch,k,w,R,es,c,name,3,6,True)
                B=run_one(d,ch,g,name,k,c,3,6)['raw_kernel']
                err=float(max(abs(A/B-1.)))
                if not ok or max(err,trace)>2e-7:
                    raise ArithmeticError((g,name,k,c,err,trace))
                out.append(dict(grid=g,pol=name,k=k,c=c,batch_vs_fixed_rel=err,vertex_trace_rel=trace))
    return out


if __name__=='__main__':
    print(json.dumps(fixture(),indent=2,allow_nan=False))
