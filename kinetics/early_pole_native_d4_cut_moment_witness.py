"""Exhibit smooth nonnegative D4 moment realizations of both exact-cut signs."""
import argparse
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_continuum_cold_band_receipt import read
from early_pole_endpoint_fine_d4_cold_continuum_ladder import DATE
from early_pole_native_d4_clipped_next_band import CUTS, SPLITS, prepare
from early_pole_native_cold_subband_from_centered_energy import interval

WIDTH = 1e-10


def center(q, mass):
    return q*q/(np.hypot(q,mass)+mass)


def bump_means(q, scalar_mass, order):
    l, c, r = q
    w = WIDTH
    if min(c-l, r-c) <= 2*w:
        raise ArithmeticError('Need four disjoint bump interiors')
    x, v = roots_legendre(order)
    t = (x+1.)/2.
    g = 15.*v*(t*(1.-t))**2
    if abs(fsum(g)-1.) > 1e-14:
        raise ArithmeticError('Normalized quartic bump failed')
    return np.array([fsum(g*center(a+w*t, scalar_mass)) for a in (l, c-w, c, r-w)])


def mixture(u, target, pair):
    j, k = pair
    t = (target-u[j])/(u[k]-u[j])
    p = np.zeros(4)
    p[[j,k]] = [1.-t,t]
    if not 0.<t<1. or abs(fsum(p*u)/target-1.) > 4e-15:
        raise ArithmeticError('Positive smooth moments must match exactly to rounding')
    return p


def nearby_aperture(rows, inside, n):
    """Exact branch zero of the old conservative lower bound, not a kernel sign."""
    p = rows[32]
    mass = p['scalar_mass']
    a,b = p['bounds']
    lo,hi = center(np.array([a,b]),mass)
    A,_,U = np.array(p['coeff'])
    mu = U/A
    # Remain on the branch with zero admissible right-side slope and
    # nonzero, unsaturated right-side inverse.
    pad = 3e-13*hi
    if np.max(mu[1:]+pad)>=center(CUTS[0],mass) or min(mu[1:])<=lo:
        raise ArithmeticError('Local upper-tail moment-bound branch changed')
    upper = interval([rows[55]],*CUTS)
    floor = inside + n*upper[2,0]-upper[1,1]
    if floor<=0.:
        raise ArithmeticError('Positive remaining interior required for this branch')
    cut_energy = lo + A[1]*(mu[1]+pad-lo)/floor
    q = float(np.sqrt(cut_energy*(cut_energy+2*mass)))
    x = [CUTS[0],q,5.001e-6]
    if not CUTS[0]<q<x[-1]<b:
        raise ArithmeticError('Aperture crossing left its expected local bracket')
    y=[]
    for a0 in x:
        W=interval(rows,a0,CUTS[1])
        v=n*W[2,0]-W[1,1]
        check=floor-A[1]*(mu[1]+pad-lo)/(center(a0,mass)-lo)
        if abs(v-check)>1e-33:
            raise ArithmeticError('Analytic branch and general conservative ledger disagree')
        y.append(float(v))
    if not y[0]<0.<y[-1]:
        raise ArithmeticError('Changed-aperture sensitivity must bracket zero')
    return dict(lower_cut_at_zero=q,lower_cut_and_conservative_net_lower=list(map(list,zip(x,y,strict=True))))



def run(root, scratch):
    results = []
    for m in (1,2):
        L, source, _, _ = prepare(root, m, scratch)
        n = next(z['common_diagnostic_occupation'] for z in L['diagnostic_root_midpoints'] if z['levels']==[2,4])
        rows = read(source)
        wholly = []
        for r in rows:
            a, b = r['bounds']
            S = np.array(r['coeff'])
            M = r['scalar_mass']
            if (not np.all(S[:,1:] > 0.) or
                    np.max(abs((M*S[0]+S[2])/S[1]-1.)) > 2e-13 or
                    np.max(abs((S[:,1]+S[:,2])/S[:,0]-1.)) > 2e-13):
                raise ArithmeticError('Independent slope/inverse measures must reproduce the stored gross moments')
            p = S[2,1:]/S[0,1:]
            if not np.all((p > center(a,M)) & (p < center(b,M))):
                raise ArithmeticError('Every uncut mean must be realizable by a smooth interior mixture')
            if CUTS[0]<=a<b<=CUTS[1]:
                wholly.append(n*S[0,2]-S[0,1])
        witnesses = []
        for alt in (0,1):
            parts = []
            for j, (id0, cuts) in enumerate(zip((32,55), SPLITS, strict=True)):
                r = rows[id0]
                if r['panel']!=id0 or r['bounds']!=[cuts[0],cuts[2]]:
                    raise ArithmeticError('Wrong exact cut panel')
                _, c, _ = cuts
                i = np.array([False,False,True,True]) if c == CUTS[0] else np.array([True,True,False,False])
                S = np.array(r['coeff'])[:,[2,1]]
                means = [bump_means(cuts, r['scalar_mass'], order) for order in (16,32)]
                if np.max(abs(means[0]/means[1]-1.))>1e-13:
                    raise ArithmeticError('Independent bump energy evaluations must agree')
                g = ((0,1),(0,2)) if (j+alt)%2==0 else ((0,2),(0,1))
                mix = [mixture(means[1], C/A, ix) for A,C,ix in zip(S[0], S[2], g, strict=True)]
                z = S[0]*np.array([p[i].sum() for p in mix])
                parts.append(dict(panel=id0, slope_then_inverse_mixture=[p.tolist() for p in mix],
                                  slope_then_inverse_mass_in_band=z.tolist(), net=n*z[0]-z[1]))
            y=fsum([*wholly,*[v['net'] for v in parts]])
            witnesses.append(dict(net=y, cut_panels=parts))
        if not witnesses[0]['net'] < -1e-22 < 0. < witnesses[1]['net']:
            raise ArithmeticError('Smooth same-moment witness signs must separate zero')
        results.append(dict(mesh=m, exact_band_interior_net=fsum(wholly), witnesses=witnesses,
                            nearby_aperture_sensitivity=nearby_aperture(rows,fsum(wholly),n)))
    return dict(q=list(CUTS), held_occupation=n, standardized_bump_width=WIDTH, rows=results,
                scope='Feasible nonnegative C1 densities for the archived slope and inverse panel moments, not evaluated collision kernels. The exact-cut new quadrature is a separate test.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('scratch', type=Path)
    print(json.dumps(run(**vars(p.parse_args())), indent=2, allow_nan=False))
