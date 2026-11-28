"""Integrated fixed-pole cold response on an opt-in continuum endpoint reconstruction."""
import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_cell_endpoint_collision_pilot import CellChart
from early_pole_cell_endpoint_transverse_audit import charge_field,continuum_moments,model_callback,piecewise_logit,root
from early_pole_endpoint_stencil_remap_pilot import SPECS,load
from finite_mass_qed_fixed_pole_scalar_q_band_reference import q_knots
from finite_mass_qed_timelike_pole_angular_reference import angle


def refine_panels(edges,splits):
    """Add exactly one strictly interior knot in each requested existing panel."""
    split=np.asarray(splits,dtype=float)
    if (split.ndim!=2 or split.shape[1]!=3 or not np.all(np.isfinite(split))
            or np.any(split[:,0]>=split[:,1]) or np.any(split[:,1]>=split[:,2])):
        raise ArithmeticError('Need finite strictly ordered split intervals')
    j=np.searchsorted(edges,split[:,0])
    if (np.any(j>=len(edges)-1) or len(j)!=len(np.unique(j))
            or np.any(edges[j]!=split[:,0]) or np.any(edges[j+1]!=split[:,2])):
        raise ArithmeticError('Requested split does not match an original panel exactly')
    out=np.union1d(edges,split[:,1])
    v=np.searchsorted(out,split[:,1])
    if np.any(out[v-1]!=split[:,0]) or np.any(out[v+1]!=split[:,2]):
        raise ArithmeticError('Split panel geometry did not close')
    return out,set(np.r_[v-1,v].tolist())


COLD_SPAN=(1e-6,1.5e-5)
FIRST_SIDE_BAND=(1.5e-5,1.8e-5)
D8_CENTRAL_COLD_BAND=(3e-6,5e-6)


def special_band_tags(card,span):
    if span==COLD_SPAN:
        return {}
    if (span==FIRST_SIDE_BAND and 'scalar_div' not in card and card.get('subcells_per_archived_cell')==2
            and card.get('steps')==2048 and card.get('x_end')==.05):
        stem='sideband'
    elif (span==D8_CENTRAL_COLD_BAND and card.get('grid')=='fine' and 'scalar_div' not in card
            and card.get('subcells_per_archived_cell')==8 and card.get('steps')==2048 and card.get('x_end')==.05):
        stem='d8_central_cold'
    else:
        raise ArithmeticError('Requested special band does not match its independently evolved native chart')
    serial=json.dumps(card,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
    return {stem+'_q_span':list(span),stem+'_state_sha256':hashlib.sha256(serial).hexdigest()}


def run(card,mesh,output,topology_only=False,mirror_azimuth=False,splits=None,patch_only=False,prior_splits=None,
        *,span=COLD_SPAN,panel_indices=None,factorized_transverse=False):
    g=card['grid']
    bundle=load(g)[0]
    tags={}
    if 'scalar_div' in card:
        from early_pole_endpoint_scalar_subshell_pilot import ScalarSubshellChart
        d=ScalarSubshellChart(bundle,SPECS[g],card['charge_div'],card['scalar_div'],card['cut'])
        bundle=dict(bundle,scalar_momentum_edges=d.ds.ed)
        serial=json.dumps(card,sort_keys=True,separators=(',',':'),allow_nan=False).encode()
        tags=dict(scalar_div=card['scalar_div'],charge_div=card['charge_div'],scalar_cut=card['cut'],
                  state_sha256=hashlib.sha256(serial).hexdigest())
    else:
        d=CellChart(bundle,SPECS[g],card['subcells_per_archived_cell'],reference_x=card['reference_x'])
    y=np.asarray(card['endpoint_occupation'])
    if card['counts']!=[d.ns,d.nf,d.np] or y.shape!=(d.n,) or not np.all(np.isfinite(y)) or min(y)<0 or max(y[d.ns:])>1:
        raise ArithmeticError('Invalid endpoint occupation card')
    x=card['x_end']
    m=d.state(x,y)[0]
    C=[]; fit=[]; callbacks=[]
    for mass,edges,n in ((x,d.p_edges,d.wp*y[d.iP:]),(x*d.ratio,d.d_edges,d.wd*y[d.ns:d.iP])):
        p,w,err,s,alpha,beta=continuum_moments(edges,mass,n,36)
        fit.append((p,w,mass))
        callbacks.append(model_callback(mass,edges,alpha,beta))
        C.append(piecewise_logit(mass,edges,alpha,beta,n,1))
    k=.0015
    _,omega,res=root(lambda b:sum(charge_field(*a,b,36) for a in fit),k)
    span=tuple(span)
    tags.update(special_band_tags(card,span))
    kernel=None
    if factorized_transverse:
        if span!=D8_CENTRAL_COLD_BAND:
            raise ArithmeticError('Optional factorized runner is initially restricted to the verified native D8 central band')
        from finite_mass_qed_transverse_factorized_spinor_reference import transverse_pair as kernel
        tags['factorized_transverse_sha256']=hashlib.sha256(Path(kernel.__code__.co_filename).read_bytes()).hexdigest()
    if d.ds.ed[0]>span[0] or d.ds.ed[-1]<span[1]:
        raise ArithmeticError('Missing retained scalar occupation support')
    u,v=[q_knots(bundle,C,k,omega,m,span,s) for s in (49,89)]
    if u.shape!=v.shape or max(abs(u-v))>9e-16:
        raise ArithmeticError('Shifted cold topology not resolved')
    if prior_splits:
        v,_=refine_panels(v,prior_splits)
    if splits is not None:
        v,changed=refine_panels(v,splits)
    elif patch_only:
        raise ArithmeticError('Patch-only mode requires a split plan')
    meshes=((2,3,6,12),(4,6,10,18),(8,10,16,24))
    order,angular,inn,az=meshes[mesh]
    xs,ws=roots_legendre(order)
    expected={'grid':g,'time_x':x,'scalar_mass':m,'omega':omega,'residue':res,'mesh':list(meshes[mesh]),'total':len(v)-1}
    expected.update(tags)
    if mirror_azimuth:
        expected['mirror_azimuth']=True
    if topology_only:
        return dict(expected,edges=v.tolist(),**({} if splits is None else
                    dict(recalculated_panel_indices=sorted(changed))))
    permitted=changed if patch_only else set(range(len(v)-1))
    if panel_indices is not None:
        if (patch_only or splits is not None or prior_splits is not None
            or not panel_indices or any(type(i) is not int or i not in permitted for i in panel_indices)
            or sorted(set(panel_indices))!=list(panel_indices)):
            raise ArithmeticError('Native shard indices must be ordered, unique and inside the original topology')
        permitted=set(panel_indices)
    done=set()
    if output.exists():
        for line in output.read_text().splitlines():
            r=json.loads(line)
            j=r['panel']
            if j not in permitted or j in done or r['bounds']!=v[j:j+2].tolist() or r.get('mirror_azimuth',False)!=mirror_azimuth or r.get('factorized_transverse_sha256')!=expected.get('factorized_transverse_sha256') or any(r[k]!=a for k,a in expected.items()):
                raise ArithmeticError('Changed endpoint or mesh in resume')
            done.add(j)
    with output.open('a',encoding='utf-8',buffering=1) as f:
        for j,(a,b) in enumerate(zip(v[:-1],v[1:],strict=True)):
            if j in done or j not in permitted:
                continue
            coeff=np.zeros((3,3))
            for q,z in zip((a+b)/2+(b-a)*xs/2,(b-a)*ws/2,strict=True):
                E=np.hypot(q,m)
                n=y[np.searchsorted(d.ds.ed,q,side='right')-1]
                kw={} if not mirror_azimuth else {'mirror_azimuth':True}
                if kernel is not None:
                    kw['transverse_vertex_kernel']=kernel
                A,B,D=angle(C,k,omega,res,E,'T',angular,(inn,az),q,occupation_callbacks=callbacks,**kw)[0]
                data=z*q*q/E*np.array([(1+n)*A,n*B,n*D+A])
                coeff[0]+=data
                coeff[1]+=E*data
                coeff[2]+=(q*q/(E+m))*data
            print(json.dumps(dict(expected,panel=j,bounds=[a,b],coeff=coeff.tolist()),allow_nan=False),
                  file=f,flush=True)
            os.fsync(f.fileno())


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('state',type=Path)
    p.add_argument('mesh',type=int,choices=(0,1,2))
    p.add_argument('output',type=Path,nargs='?')
    p.add_argument('--topology-only',action='store_true')
    p.add_argument('--mirror-azimuth',action='store_true')
    p.add_argument('--split-plan',type=Path)
    p.add_argument('--patch-only',action='store_true')
    p.add_argument('--prior-split-plan',type=Path)
    band=p.add_mutually_exclusive_group()
    band.add_argument('--first-sideband',action='store_true')
    band.add_argument('--d8-central-cold-band',action='store_true')
    p.add_argument('--factorized-transverse',action='store_true')
    o=p.parse_args()
    if o.output is None and not o.topology_only:
        p.error('output is required except in topology-only mode')
    splits=None if o.split_plan is None else json.loads(o.split_plan.read_text())['targeted_parent_panel_splits']
    prior=None if o.prior_split_plan is None else json.loads(o.prior_split_plan.read_text())['targeted_parent_panel_splits']
    result=run(json.loads(o.state.read_text()),o.mesh,o.output,o.topology_only,o.mirror_azimuth,splits,o.patch_only,prior_splits=prior,
               span=FIRST_SIDE_BAND if o.first_sideband else D8_CENTRAL_COLD_BAND if o.d8_central_cold_band else COLD_SPAN,
               factorized_transverse=o.factorized_transverse)
    if result is not None:
        print(json.dumps(result,indent=2,allow_nan=False))
