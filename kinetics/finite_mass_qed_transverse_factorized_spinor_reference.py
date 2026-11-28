"""Factor constant transverse vertices into two-by-two spinor contractions."""
import numpy as np

from early_pole_hard_yukawa_photoproduction_reference import G, slash, spinors


def denominator(o,m,w,k):
    spatial=o[...,1:]
    z=spatial[...,2]
    q=np.linalg.norm(spatial,axis=-1)
    turn=q-z
    mask=z>0.
    turn[mask]=np.sum(spatial[mask,:2]**2,axis=-1)/(q[mask]+z[mask])
    mu2=(w-k)*(w+k)
    return (m*m*w*w+mu2*q*q)/(np.hypot(m,q)*w+q*k)+k*turn


def transverse_pair(p,b,K,M,m):
    if K.shape!=(4,) or not np.array_equal(K[1:3],[0.,0.]) or not K[0]>K[3]>0.:
        raise ValueError('Require a forward timelike pole in the exact z-direction')
    w,k=K[0],K[3]
    mu2=(w-k)*(w+k)
    if min(M,m)<=0. or np.sqrt(mu2)>=2.*m:
        raise ValueError("Require the nonsingular finite-mass transverse branch")
    A=1./(mu2+2.*denominator(p,M,w,k))
    B=1./(mu2-2.*denominator(b,m,w,k))
    U=spinors(p,M)
    V=np.swapaxes(spinors(b,m).conj(),-1,-2)@G[0]
    S=V@U
    T=A+B
    Kslash=slash(K)
    out=np.zeros_like(A)
    for j in (1,2):
        # {Kslash,epsslash}=0 for the x and y unit transverse vectors.
        # Complete the amplitude *before* taking its absolute square.
        epsilon=np.eye(4)[j]
        C=Kslash@slash(epsilon)
        a=-2*(p[...,j]*A+b[...,j]*B)
        value=a[...,None,None]*S+T[...,None,None]*(V@(C@U))
        out+=np.sum(abs(value)**2,axis=(-1,-2))
    return out
