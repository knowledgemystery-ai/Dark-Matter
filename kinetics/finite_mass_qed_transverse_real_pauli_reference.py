"""Real Pauli-coordinate spin sum on the nonsingular finite-mass transverse branch."""
import numpy as np

from finite_mass_qed_transverse_factorized_spinor_reference import denominator


def inner(u,v):
    return np.sum(u*v,axis=-1)


def transverse_real_pauli(p,b,K,M,m):
    if K.shape!=(4,) or not np.array_equal(K[1:3],[0.,0.]) or not K[0]>K[3]>0.:
        raise ValueError('Require a forward timelike pole in the exact z-direction')
    w,k=K[0],K[3]
    mu2=(w-k)*(w+k)
    if min(M,m)<=0. or np.sqrt(mu2)>=2.*m:
        raise ValueError('Require the nonsingular finite-mass transverse branch')
    A=1./(mu2+2.*denominator(p,M,w,k))
    B=1./(mu2-2.*denominator(b,m,w,k))
    t=A+B
    P=np.sqrt(p[...,0]+M)
    F=np.sqrt(b[...,0]+m)
    v=p[...,1:]/P[...,None]
    u=b[...,1:]/F[...,None]
    uv=inner(u,v)
    WF=P*F
    iS=WF-uv
    sS=-np.cross(u,v)
    out=np.zeros_like(A)
    for j in (0,1):
        e=np.eye(3)[j]
        h=np.cross(np.array([0.,0.,1.]),e)
        ia=-2.*(A*p[...,j+1]+B*b[...,j+1])
        # For C=Kslash epsslash, write ubar_b C u_p = cI + i σ·g.
        c=w*(P*u[...,j]-F*v[...,j])-k*inner(np.cross(u,h),v)
        g=(-k*(WF+uv)[...,None]*h-w*F[...,None]*np.cross(e,v)+w*P[...,None]*np.cross(u,e)
           +k*(inner(u,h)[...,None]*v+inner(h,v)[...,None]*u))
        alpha=ia*iS+t*c
        beta=ia[...,None]*sS+t[...,None]*g
        out+=2.*(alpha*alpha+inner(beta,beta))
    return out
