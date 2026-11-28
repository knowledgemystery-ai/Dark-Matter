"""Patch only the two native D4 panels crossing q=5e-6 or 8e-6."""
import argparse
import gzip
import json
from math import fsum
from pathlib import Path

import numpy as np

from early_pole_endpoint_cold_targeted_patch_join import merge
from early_pole_endpoint_continuum_cold_band_pilot import run as pilot
from early_pole_endpoint_continuum_cold_band_receipt import read, reduce
from early_pole_endpoint_fine_d4_cold_continuum_ladder import DATE, unpacked_hash

CUTS = (5e-6, 8e-6)
EXPECTED_MESHES = ([4, 6, 10, 18], [8, 10, 16, 24])
SPLITS = [[4.778274331905574e-6, CUTS[0], 5.083582387749736e-6],
          [7.863047969834538e-6, CUTS[1], 8.04609200708498e-6]]
IDS = [32, 33, 56, 57]
TOL = 6e-5


def prepare(root, mesh, work, compute=False):
    L = json.loads((root / f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    state_file = root / f'early-pole-cell-endpoint-fine-d4-x05-state-{DATE}.json.gz'
    source = root / f'fine-d4-endpoint-continuum-cold-{mesh}.jsonl.gz'
    expected = L['raw_uncompressed_sha256'][2][mesh - 1]
    if unpacked_hash(source) != expected or unpacked_hash(state_file) != L['cards'][2]['state_sha256_uncompressed']:
        raise ArithmeticError('Expected exact admitted D4 endpoint and native ledger')
    _, bounds, meta = reduce(source)
    state = json.loads(gzip.decompress(state_file.read_bytes()))
    scratch = work / f'd4-cold-split-{mesh}.jsonl'
    stock = pilot(state, mesh, scratch, topology_only=True)
    if stock['edges'] != [bounds[0][0], *[b for _, b in bounds]] or any(stock[k] != v for k, v in meta.items()):
        raise ArithmeticError('Changed original stock quadrature topology or metadata')
    plan = pilot(state, mesh, scratch, topology_only=True, splits=SPLITS, patch_only=True)
    if plan['recalculated_panel_indices'] != IDS or plan['mesh'] != EXPECTED_MESHES[mesh - 1]:
        raise ArithmeticError('Expected the four preregistered children')
    if compute:
        pilot(state, mesh, scratch, splits=SPLITS, patch_only=True)
    return L, source, scratch, plan


def verify_children(source, children):
    original = {tuple(r['bounds']): np.array(r['coeff']) for r in read(source)}
    if [r['panel'] for r in children] != IDS:
        raise ArithmeticError('Expected all four children in exact order')
    vals = []
    for (lo, _, hi), r, s in zip(SPLITS, children[::2], children[1::2], strict=True):
        a = np.array([r['coeff'], s['coeff']])
        v = original[(lo, hi)]
        if np.any(v <= 0) or np.any(a <= 0) or not np.isfinite(a).all():
            raise ArithmeticError('Need positive finite gross and affine moments')
        rel = abs(a.sum(axis=0) / v - 1.)
        if np.max(rel) > TOL:
            raise ArithmeticError('Original number, energy or centered component fails local closure')
        vals.append(float(rel.max()))
    return vals


def finish(root, work):
    R, rows, local = [], [], []
    for mesh in (1, 2):
        L, base, path, plan = prepare(root, mesh, work)
        dat = read(path)
        local.append(verify_children(base, dat))
        joined = work / f'd4-cold-clipped-joined-{mesh}.jsonl'
        if joined.exists():
            raise FileExistsError(f'Choose a new output directory for {joined}')
        info = merge(base, path, plan, joined)
        eligible = [r for r in read(joined) if r['bounds'][0] >= CUTS[0] and r['bounds'][1] <= CUTS[1]]
        if not eligible or eligible[0]['bounds'][0] != CUTS[0] or eligible[-1]['bounds'][1] != CUTS[1]:
            raise ArithmeticError('Selected band must end at the exact installed edges')
        x = np.array([r['coeff'] for r in eligible])
        totals = np.array([[fsum(x[:, i, j]) for j in range(3)] for i in range(3)])
        R.append(x)
        rows.append(dict(mesh=plan['mesh'], selected_panels=len(x), moments=totals.tolist(),
                         stock_split_patch_sha256=unpacked_hash(path), join=info))
    if R[0].shape != R[1].shape or np.any(np.asarray(R) <= 0):
        raise ArithmeticError('The two exact clipped support arrays must match and be positive')
    dx = float(np.max(abs(R[0] / R[1] - 1.)))
    if dx > TOL:
        raise ArithmeticError('Panel componentwise two-mesh difference is underresolved')
    n = next(z['common_diagnostic_occupation'] for z in L['diagnostic_root_midpoints'] if z['levels'] == [2, 4])
    X = np.asarray(R)[:, :, 0, :]
    lo, hi = X.min(axis=0), X.max(axis=0)
    b = [fsum(n * lo[:, 2] - hi[:, 1]), fsum(n * hi[:, 2] - lo[:, 1])]
    return dict(q=list(CUTS), old_panels=SPLITS, children=IDS, rows=rows,
                maximum_relative_componentwise_mesh_difference=dx,
                relative_parent_reclosure_by_mesh_and_panel=local,
                diagnostic_occupation=n, two_mesh_coefficient_box_net=b,
                held_sign=int(b[0] > 0) - int(b[1] < 0),
                scope='Exact clipping and an observed radial/angular two-mesh box on the held D4 state. This is not a proved continuum integration error bound or a D8 result.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('work', type=Path)
    p.add_argument('mode', choices=('plan', 'mesh1', 'mesh2', 'finish'))
    a = p.parse_args()
    if a.mode == 'finish':
        print(json.dumps(finish(a.root, a.work), indent=2, allow_nan=False))
    else:
        for i in ((1, 2) if a.mode == 'plan' else (int(a.mode[-1]),)):
            print(json.dumps(prepare(a.root, i, a.work, a.mode != 'plan')[-1], allow_nan=False), flush=True)
