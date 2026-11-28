"""Bound the change in a compressed symmetric rate between two equal-dimensional planes."""
import numpy as np


def comparison(G,S,U,V):
    """Reported angle in the entropy Euclidean frame, not the raw fraction frame.

    The Rayleigh change bound uses the full-form spectral spread times the sine
    of the largest principal angle. It bounds differences between these two
    finite subspaces, not their distance from an unknown exact plane.
    """
    R=np.linalg.cholesky(G).T
    rhs=np.linalg.solve(R.T,S)
    X=np.linalg.solve(R.T,rhs.T).T
    X=(X+X.T)*.5
    A,B=[np.linalg.qr(R@z,mode='reduced')[0] for z in (U,V)]
    P,Q=A@A.T,B@B.T
    d=np.linalg.norm(P-Q,ord=2)
    L=np.linalg.eigvalsh(X)
    E,F=[np.linalg.eigvalsh(z.T@X@z) for z in (A,B)]
    bound=(L[-1]-L[0])*d
    if max(abs(E-F))>bound+1e-10*max(1.,max(abs(L))):
        raise ArithmeticError('Compressed rates violate their principal-angle change bound')
    return dict(rate_first=E.tolist(),rate_second=F.tolist(),largest_principal_sine=float(d),
                spectral_spread=float(L[-1]-L[0]),rayleigh_rate_change_bound=float(bound))


def projection_rate_bound_from_fractional_sine(G,S,d):
    """Conservative entropy-rate sensitivity bound from a raw fractional plane angle."""
    if not np.isfinite(d) or not 0.<=d<=1.+2e-14:
        raise ValueError('Expected a principal-angle sine')
    R=np.linalg.cholesky(G).T
    rhs=np.linalg.solve(R.T,S)
    Z=np.linalg.solve(R.T,rhs.T).T
    lam=np.linalg.eigvalsh((Z+Z.T)*.5)
    ent_sine_bound=min(1.,np.linalg.cond(R,2)*d)
    return float(ent_sine_bound),float(ent_sine_bound*(lam[-1]-lam[0]))
