"""Sharp partial-mass bounds for the archived discrete radial quadratures."""
import argparse
import hashlib
import json
from math import fsum
from pathlib import Path

import numpy as np
from scipy.special import roots_legendre

from early_pole_endpoint_continuum_cold_band_receipt import reduce
from early_pole_endpoint_fine_d4_cold_continuum_ladder import unpacked_hash
from early_pole_native_cold_bose_reweight import ledger
from early_pole_native_cold_cross_factors import DATE


def fraction(nodes, mask, mu, pad):
    """Bounds on in-band mass using an interval for the mean centered energy."""
    mu = np.asarray(mu)
    if (min(mu) < nodes[0] - pad or max(mu) > nodes[-1] + pad
            or not np.isfinite(mu).all()):
        raise ArithmeticError('Mean is incompatible with the discrete support')
    lo = np.zeros_like(mu)
    hi = np.ones_like(mu)
    if not mask.any():
        return np.c_[lo, lo]
    if mask.all():
        return np.c_[hi, hi]
    i, o = nodes[mask], nodes[~mask]
    # Feasible centers are between the minima and maxima of the conditional
    # supports. Each extremum is attained by mixing the corresponding nodes.
    d = i[0] - o[0]
    if d > 0:
        hi = np.minimum(hi, (mu + pad - o[0]) / d)
    elif d < 0:
        lo = np.maximum(lo, (mu + pad - o[0]) / d)
    d = i[-1] - o[-1]
    if d > 0:
        lo = np.maximum(lo, (mu - pad - o[-1]) / d)
    elif d < 0:
        hi = np.minimum(hi, (mu - pad - o[-1]) / d)
    if np.any(lo > hi):
        raise ArithmeticError('No feasible restricted mass fraction')
    return np.c_[lo, hi]


def panel(v, a, b):
    l, r = v['bounds']
    A, E, U = np.array(v['coeff'])
    m = v['scalar_mass']
    if min(A[:2]) <= 0 or np.max(abs((m * A + U) / E - 1.)) > 2e-13:
        raise ArithmeticError('Gross moments must be positive with energy closure')
    x = roots_legendre(v['mesh'][0])[0]
    q = (l + r) / 2 + (r - l) * x / 2
    center = q * q / (np.hypot(q, m) + m)
    mask = (q >= a) & (q < b)
    if min(np.min(abs(q-a)), np.min(abs(q-b))) <= 32 * np.finfo(float).eps * max(abs(q)):
        raise ArithmeticError('Band cut is numerically coincident with an underlying node')
    box = A[:2, None] * fraction(center, mask, U[:2] / A[:2], 3e-13 * center[-1])
    return box, bool(mask.any() and not mask.all())


def run(root, prior):
    raw = prior.read_bytes()
    h0 = hashlib.sha256(raw).hexdigest()
    if h0 != "50e3a0675bfe1d77a8b71d7f64d5bc1b20c2cca7829886358f6dff9cf5758f27":
        raise ArithmeticError("Expected the exact continuous-support comparison")
    Z = json.loads(raw)
    L = json.loads((root / f'early-pole-fine-d4-cold-continuum-ladder-{DATE}.json').read_bytes())
    n = next(z['common_diagnostic_occupation'] for z in L['diagnostic_root_midpoints'] if z['levels'] == [2, 4])
    if n != Z['diagnostic_occupation']:
        raise ArithmeticError('Held occupation changed')
    data, hashes = [], []
    for k, ix in ((2, 1), (4, 2)):
        p = [root / f'fine-d{k}-endpoint-continuum-cold-{m}.jsonl.gz' for m in (1, 2)]
        h = [unpacked_hash(f) for f in p]
        if h != L['raw_uncompressed_sha256'][ix]:
            raise ArithmeticError('Need exact admitted full-band quadrature meshes')
        for f in p:
            reduce(f)
        data.append([ledger(f) for f in p])
        hashes.append(h)
    out = []
    for v in Z['rows']:
        charts = []
        for k, V, c in zip((2, 4), data, v['charts'], strict=True):
            values = []
            for R in V:
                p = [panel(w, *v['q']) for w in R]
                X = np.array([m[0] for m in p])
                F, B = [[fsum(X[:, j, t]) for t in (0, 1)] for j in (0, 1)]
                net = [n * F[0] - (1 + n) * B[1], n * F[1] - (1 + n) * B[0]]
                if F[0] <= B[1]:
                    raise ArithmeticError("Need a resolved positive slope for the neutral bound")
                neutral = [B[0] / (F[1] - B[0]), B[1] / (F[0] - B[1])]
                values.append(dict(forward=F, inverse=B, affine_net=net, neutral_interval=neutral,
                                   partial_node_panels=sum(m[1] for m in p)))
            hull = [min(s['affine_net'][0] for s in values), max(s['affine_net'][1] for s in values)]
            sign = int(hull[0] > 0) - int(hull[1] < 0)
            if c['held_sign'] and sign != c['held_sign']:
                raise ArithmeticError('Node-aware signs contradict admitted continuous-support signs')
            charts.append(dict(level=k, meshes=values, held_net_hull=hull, held_sign=sign,
                               continuous_support_held_net_hull=c['held_net_hull'],
                               continuous_support_held_sign=c['held_sign']))
        out.append(dict(q=v['q'], charts=charts))
    return dict(rows=out, diagnostic_occupation=n,
                raw_uncompressed_sha256=hashes,
                continuous_support_receipt_sha256=h0,
                scope='Restricted mass of archived Gauss-Legendre atomic quadratures. A band edge filters existing nodes; it is not a new integration on the clipped panel, an error bound for the continuous q integral, or a D8 result.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('root', type=Path)
    p.add_argument('prior', type=Path)
    print(json.dumps(run(**vars(p.parse_args())), indent=2, allow_nan=False))
