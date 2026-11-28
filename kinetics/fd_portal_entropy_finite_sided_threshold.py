"""Separate pre-knot and post-knot control of the held finite Gram gain."""
import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from numpy.polynomial.legendre import leggauss
from scipy.linalg import eigh

from fd_portal_entropy_finite_rank_threshold import (CombinedExpansion, DATE, PositivePortal, StandardModel,
                                                     dense, gain_difference, local_path, matrix, sym, terms, threshold)


def boundary(H, fixed, moved, fixed_multiplier):
    """Boundary in the nonnegative adjustable multiplier, or no feasible value."""
    base = H + fixed_multiplier * fixed
    top = np.linalg.eigvalsh(sym(base))[-1]
    if top > 0:
        return None
    if top == 0:
        return 0.
    return threshold(base, moved)


def run(eos, root):
    fragments = ('direct-expansion-096-120', 'entropy-short-tangent-windows',
                 'entropy-finite-gain-budget', 'entropy-finite-rank-threshold')
    raw = [(root / f'fd-portal-{k}-{DATE}.json').read_bytes() for k in fragments]
    A, T, B, R0 = map(json.loads, raw)
    ids = [hashlib.sha256(s).hexdigest() for s in raw]
    he = hashlib.sha256(eos.read_bytes()).hexdigest()
    if (T['input_sha256'][0] != ids[0] or B['input_sha256'][:2] != ids[:2]
            or R0['input_sha256'] != ids[:3] or any(z['table_input_sha256'] != he for z in (A, T, B, R0))):
        raise ArithmeticError('Need the exact trajectory, reference budget and uniform matrix receipt')
    a, b = B['window']
    K = [a, T['windows'][0]['x'][1], T['split_at_knot'], T['windows'][2]['x'][0], b]
    P = CombinedExpansion(StandardModel(eos), PositivePortal((17, 44, 32)), 168)
    init = next(p for q in A['rows'] if q['operator'] == 'combined' for p in q['points'] if p[0] == a)
    path = local_path(P, init)
    Q = np.asarray(T['rebased_fractional_columns'])
    out = []
    for (n, e, k), check in zip(((4, 4e-5, 3), (6, 2e-5, 6)), R0['rows'], strict=True):
        C = dense(P, path, K, Q, n, e)
        M = C(a)
        X = sym(M.T @ matrix(P, path(a))[0] @ M)
        H = np.zeros((2, 2))
        R = np.zeros((2, 2, 2))
        x, w = leggauss(k)
        for i, (s, t) in enumerate(zip(K[:-1], K[1:], strict=True)):
            for y, z in zip((s + t) / 2 + (t - s) * x / 2, (t - s) * w / 2, strict=True):
                M = C(y)
                V = terms(P, path(y), e)
                R[i // 2] += z * sym(M.T @ V[2] @ M)
                H += z * sym(M.T @ sum(V[j] for j in (0, 1, 3, 4)) @ M)
        defects = [float(max(abs(eigh(sym(u - v), X, eigvals_only=True))))
                   for u, v in ((H, check['integrated_reference_matrix']), (R.sum(axis=0), check['integrated_slope_matrix']))]
        if max(defects) > 3e-12 or np.linalg.eigvalsh(H)[-1] >= 0:
            raise ArithmeticError('Grouped integration disagrees with the accepted full window')
        moments = []
        for fixed, moved in (R, R[::-1]):
            if np.linalg.eigvalsh(fixed)[0] <= 0:
                raise ArithmeticError('Each half must independently have positive integrated excess form')
            # Fixed-side fraction is a prescribed grid, whether or not it
            # permits an adjustable-side contraction in [0, 1].
            caps = [boundary(H, fixed, moved, f) for f in (0., .25, .5, .75, 1.)]
            if caps[-1] is not None and abs(gain_difference(X, X + H + fixed + caps[-1]*moved)[-1]) > 7e-13:
                raise ArithmeticError('Side-specific threshold did not reproduce endpoint neutrality')
            moments.append(dict(critical_other_multiplier=caps,
                                largest_gain_minus_one_if_other_removed=gain_difference(X, X + H + fixed)[-1]))
        out.append(dict(hot_then_cool_duration=[K[2]-a, b-K[2]],
                        hot_then_cool_integrated_excess_gram=R.tolist(),
                        consistency_with_uniform_receipt=defects,
                        fixed_hot_then_fixed_cool=moments))
    return dict(window=[a, b], pchip_knot=K[2], fixed_multipliers=[0., .25, .5, .75, 1.], rows=out,
                input_sha256=ids, table_input_sha256=he,
                scope='Independent multiplier on each side of the EOS knot, on one unchanged tangent and temperature path. No change to the EOS, new evolution, or physical stability claim.')


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('eos', type=Path)
    p.add_argument('root', type=Path)
    print(json.dumps(run(**vars(p.parse_args())), indent=2, allow_nan=False))
