"""Plan, run, and exactly reduce disjoint cold-panel workers on one approved topology."""
import argparse
import gzip
import fcntl
import hashlib
import json
from pathlib import Path

import numpy as np

from early_pole_endpoint_continuum_cold_band_pilot import COLD_SPAN, run as pilot
from early_pole_endpoint_scalar_subshell_chain_receipt import read

MESHES=((2,3,6,12),(4,6,10,18),(8,10,16,24))


def digest(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,separators=(',',':'),allow_nan=False).encode()).hexdigest()


def plan(topology,number):
    edge=topology['edges']
    n=topology['total']
    if (not isinstance(number,int) or not 1<=number<=n or len(edge)!=n+1
        or not np.all(np.isfinite(edge)) or any(a>=b for a,b in zip(edge[:-1],edge[1:],strict=True))
        or tuple(topology['mesh']) not in MESHES):
        raise ArithmeticError('Valid exact topology and nonempty shard count required')
    p=dict(topology=topology,topology_sha256=digest(topology),shards=[list(range(i,n,number)) for i in range(number)])
    return dict(**p,payload_sha256=digest(p))


def verify(m):
    want=plan(m['topology'],len(m['shards']))
    if want!=m:
        raise ArithmeticError('A shard manifest must bind the exact panel partition and topology')
    return m['topology']


def work(manifest,state,idx,output):
    t=verify(manifest)
    if not 0<=idx<len(manifest['shards']):
        raise ArithmeticError('Select one valid shard index only')
    a,_=read(state)
    mesh=MESHES.index(tuple(t['mesh']))
    x=tuple(t.get('d8_central_cold_q_span',t.get('sideband_q_span',COLD_SPAN)))
    kw=dict(span=x,mirror_azimuth=t.get('mirror_azimuth',False))
    if 'factorized_transverse_sha256' in t:
        kw['factorized_transverse']=True
    fresh=pilot(a,mesh,None,topology_only=True,**kw)
    if fresh!=t:
        raise ArithmeticError('Original state, pole, geometry or mesh do not match the manifest')
    # Prevent two processes on this cloud workspace writing the same shard.
    with output.with_suffix(output.suffix+'.lock').open('a+b') as guard:
        fcntl.flock(guard.fileno(),fcntl.LOCK_EX|fcntl.LOCK_NB)
        pilot(a,mesh,output,panel_indices=manifest['shards'][idx],**kw)


def lines(f):
    op=gzip.open if f.suffix=='.gz' else open
    with op(f,'rb') as q:
        return q.readlines()


def join(manifest,src):
    t=verify(manifest)
    observed={}
    if len(src)!=len(manifest['shards']):
        raise ArithmeticError('Exactly one completed file per declared shard is required')
    for allowed,f in zip(manifest['shards'],src,strict=True):
        rows=lines(f)
        if len(rows)!=len(allowed):
            raise ArithmeticError('Incomplete or oversized native shard')
        for j,u in zip(allowed,rows,strict=True):
            v=json.loads(u)
            if (not u.endswith(b'\n') or v['panel']!=j or j in observed or v['bounds']!=t['edges'][j:j+2]
                or any(v.get(k)!=x for k,x in t.items() if k!='edges')
                or v.get('mirror_azimuth',False)!=t.get('mirror_azimuth',False)
                or set(v.keys())!=((set(t)-{'edges'})|{'panel','bounds','coeff'})
                or np.shape(v['coeff'])!=(3,3) or not np.all(np.isfinite(v['coeff']))):
                raise ArithmeticError('Wrong native shard panel, repeat, state or support')
            observed[j]=u
    if set(observed)!=set(range(t['total'])):
        raise ArithmeticError('Do not reduce an incomplete native band')
    return b''.join(observed[i] for i in range(t['total']))


def main(args):
    if args.action=='plan':
        data=plan(json.loads(args.topology.read_bytes()),args.count)
        print(json.dumps(data,indent=2,allow_nan=False))
    elif args.action=='run':
        work(json.loads(args.manifest.read_bytes()),args.state,args.index,args.output)
    else:
        out=join(json.loads(args.manifest.read_bytes()),args.shards)
        # A complete check precedes the first output write.
        tmp=args.output.with_suffix(args.output.suffix+'.tmp')
        tmp.write_bytes(out)
        tmp.replace(args.output)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='action',required=True)
    a=sub.add_parser('plan'); a.add_argument('topology',type=Path); a.add_argument('count',type=int)
    b=sub.add_parser('run'); b.add_argument('manifest',type=Path); b.add_argument('state',type=Path); b.add_argument('index',type=int); b.add_argument('output',type=Path)
    c=sub.add_parser('join'); c.add_argument('manifest',type=Path); c.add_argument('output',type=Path); c.add_argument('shards',type=Path,nargs='+')
    main(p.parse_args())
