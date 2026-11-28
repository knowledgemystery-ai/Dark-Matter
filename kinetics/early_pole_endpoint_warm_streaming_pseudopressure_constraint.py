"""Add actual native pseudo-pressure to the held streaming shell moments."""
import argparse
import hashlib
import json
from pathlib import Path
from math import fsum

import numpy as np
from scipy.optimize import linprog

from early_pole_endpoint_scalar_gauss2_energy_reconstruction import scaled_energy,shell_rule
from early_pole_endpoint_scalar_subshell_chain_receipt import read
from early_pole_endpoint_scalar_warm_evolution_audit import unpack
from early_pole_endpoint_warm_three_chart_stress_hull import conditioned
from early_pole_endpoint_warm_three_chart_streaming_hull import curvature,kernel


def feature(t,E,h):
    r=h/E
    return np.array([np.ones_like(t),t,t*(1.-t)/(1.+r*t),t**3*(4.+3.*r*t)/(1.+r*t)**3])


def cell(a,b,m,A,phase,q,w,n):
    E=np.hypot(a,m)
    H=(b-a)*(b+a)/(np.hypot(b,m)+E)
    t=np.linspace(0.,1.,n+1)
    Q=np.sqrt(a*a+t*H*(2*E+t*H))
    qbar=scaled_energy(q,a,b,m)
    T=np.r_[t,qbar]
    U=feature(T,E,H)
    raw=feature(qbar,E,H)
    target=np.array([fsum(s*w) for s in raw])
    Qall=np.r_[Q,q]
    f=kernel(Qall,m,A,phase)
    center=kernel(q,m,A,phase)
    witness=fsum(w*center)
    # Subtract terms already determined by the constraints, then scale the
    # remaining objective so the LP's absolute objective tolerance is harmless.
    known=np.linalg.lstsq(U.T,f,rcond=None)[0]
    residue=f-known@U
    fact=max(np.max(abs(residue)),2e-15)
    P=residue/fact
    extrema=[]; outer=[]; defects=[]; atoms=[]
    for sense in (1.,-1.):
        z=linprog(sense*P,A_eq=U,b_eq=target,bounds=(0.,None),method='highs',
                  options={'dual_feasibility_tolerance':1e-9,'primal_feasibility_tolerance':1e-9})
        if not z.success:
            raise ArithmeticError(('Positive archive-containing mesh must be feasible',z.message))
        idx=np.flatnonzero(z.x>1e-12)
        if not 1<=len(idx)<=4:
            raise ArithmeticError('Expected a basic four-constraint witness')
        weights=np.linalg.lstsq(U[:,idx],target,rcond=None)[0]
        if min(weights)<-2e-13:
            raise ArithmeticError('Renormalized atomic witness is not positive')
        weights=np.maximum(weights,0.)
        err=np.max(abs(U[:,idx]@weights-target))
        if err>2e-13:
            raise ArithmeticError('Positive witness does not preserve all four constraints')
        y=z.eqlin.marginals*fact+sense*known
        rest=sense*f-y@U
        tau=max(0.,float(-np.min(rest)))
        # In normalized excitation energy, J''=-2 E0^2 (E0+H)/E(t)^3.
        B=curvature(a,b,m,A,phase)+abs(y[2])*2*(E+H)/E+abs(y[3])*12*(2.+H/E)
        outer.append(sense*(float(y@target)-tau-B/(8*n*n)))
        extrema.append(float(f[idx]@weights))
        defects.append(err)
        atoms.append(dict(q=Qall[idx].tolist(),weights=weights.tolist()))
    if not outer[0]-3e-14<=min(extrema[0],witness)<=max(extrema[1],witness)<=outer[1]+3e-14:
        raise ArithmeticError('Positive/dual cell interval audit did not close')
    return np.array(outer),np.array(extrema),witness,max(defects),atoms


def average(c,m,A,phase,n,record=False):
    out=[]; inner=[]; node=[]; errors=[]; mass=[]; supports=[]
    for j,(a,b,z,_,q,p,warm) in enumerate(c):
        if not z:
            continue
        if warm:
            O,I,w,e,atoms=cell(a,b,m,A,phase,q,p/z,n)
            if record:
                supports.append(dict(cell=j,bounds=[a,b],count=z,lower=atoms[0],upper=atoms[1]))
        else:
            Q,_,L=shell_rule(a,b,m,n)
            P=np.exp(L); w=fsum(P*kernel(Q,m,A,phase))/fsum(P)
            O=I=np.array([w,w]); e=0.
        out.append(z*O); inner.append(z*I); node.append(z*w); errors.append(e); mass.append(z)
    N=fsum(mass)
    return [np.array([fsum(v) for v in zip(*x,strict=True)])/N for x in (out,inner)]+[fsum(node)/N,max(errors),supports]


def run(native,accepted,control,stream,prior,states):
    raw=[p.read_bytes() for p in (accepted,stream,prior)]
    A,B,C=map(json.loads,raw)
    h=[hashlib.sha256(a).hexdigest() for a in raw]
    u,uh=read(native); cr=hashlib.sha256(control.read_bytes()).hexdigest()
    V=[read(p) for p in states]
    R=A['rows'][-1]
    if (R['x']!=.05 or A['charts']!=['G8','G9','mixed300'] or R['states_sha256']!=[q for _,q in V]
        or C['input_sha256'][3]!=h[1] or B['input_sha256']!=[uh,h[0],cr,*R['states_sha256']]):
        raise ArithmeticError('Require exact prior mature charts and their identical phase convention')
    d,y=unpack(u,False,.05)
    m=d.state(.05,y)[0]
    if m!=C['formal_mass'] or m!=R['common_native_mass']:
        raise ArithmeticError('Common measurement mass changed')
    charts=[conditioned(d,x,m,.05,c) for (x,_),c in zip(V,(None,None,control),strict=True)]
    if [k for k,_,_ in charts]!=[8,9,'mixed300']:
        raise ArithmeticError('Chart identities changed')
    rows=[]
    for O in C['rows']:
        a,b=O['A'],O['phase']
        coarse,fine=[[average(c,m,a,b,n,record=(a==100. and n==768 and i<2))
                      for i,(_,c,_) in enumerate(charts)] for n in (192,768)]
        cards=[]
        for i,(g,l,w,e,atoms) in enumerate(fine):
            if max(abs(coarse[i][2]-w),abs(O['charts'][i]['nodal']-w))>2e-14:
                raise ArithmeticError('Independent and already accepted positive atom reference disagree')
            if not coarse[i][0][0]-2e-14<=l[0]<=l[1]<=coarse[i][0][1]+2e-14:
                raise ArithmeticError('Fine feasible positive witness violates coarse valid envelope')
            p0,p1=O['charts'][i]['outer']
            g=np.array([max(g[0],p0),min(g[1],p1)])
            if not g[0]-2e-14<=l[0]<=l[1]<=g[1]+2e-14:
                raise ArithmeticError('Four-moment witnesses violate the intersected three-moment envelope')
            cards.append(dict(outer=g.tolist(),constructive=l.tolist(),nodal=w,maximum_witness_moment_defect=e))
            if atoms:
                cards[-1]['positive_warm_atoms']=atoms
        pairs=[]
        for j,k in ((0,1),(0,2),(1,2)):
            q,v=cards[j],cards[k]
            def gap(field):
                l,h=q[field]; a,b=v[field]
                return [a-h,b-l]
            R,I=gap('outer'),gap('constructive')
            sign=('positive' if R[0]>0 else 'negative' if R[1]<0 else 'both_witnessed' if I[0]<0<I[1] else 'inconclusive')
            pairs.append(dict(charts=[A['charts'][k],A['charts'][j]],outer=R,constructive=I,sign=sign))
        rows.append(dict(A=a,phase=b,charts=cards,pairs=pairs))
    return dict(x=.05,formal_mass=m,rows=rows,input_sha256=[uh,cr,*h,*A['rows'][-1]['states_sha256']],
                scope='Native pressure and pseudo-pressure fixed independently per warm cell; finite constructive witnesses plus ordinary-floating-point dual/interpolation bounds. Not cosmological transfer.')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('native','accepted','control','stream','prior'):
        p.add_argument(k,type=Path)
    p.add_argument('states',type=Path,nargs=3)
    print(json.dumps(run(**vars(p.parse_args())),indent=2,allow_nan=False))
