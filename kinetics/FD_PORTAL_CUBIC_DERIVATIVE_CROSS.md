# The smooth-EOS sign near the knot is not explained by its heat-capacity slope alone

The [rank-one slope treatment](FD_PORTAL_HEAT_CAPACITY_SLOPE_THRESHOLD.md) accounts for the original PCHIP sign split at its knot. A different question is whether it also explains the behavior of the *complete* smooth cubic EOS at those same physical states. Hold visible temperature and both dark variables fixed. Cross just d log(c_v)/d log(T) onto either held base rate form and leave that form's metric, motion, and collision Jacobian unchanged. Both uncrossed cases independently replay the archived local rates.

| Held state | Maximum rate, PCHIP form with only the cubic slope | Maximum rate, full cubic form |
|---|---:|---:|
| x=.61 | −.00695019 | −.00705701 |
| x=.99 | −.00062643 | −.00226625 |
| Wider hot-side probe | +.00017424 | −.00149115 |
| Wider cool-side probe | +.00036247 | −.00131335 |
| Closer hot-side probe | +.00041781 | −.00125539 |
| Closer cool-side probe | +.00047664 | −.00119982 |

The slope-only substitution changes the original maximum to negative at the two archived clock points, but predicts the **wrong sign in all four** probes bracketing the knot. The native cubic slope at those probes is just above the PCHIP-form threshold near 3.37253, but below the cubic-form threshold near 3.37455. Thus one must retain the rest of the form, not merely smooth the one derivative, to recover the smooth-EOS sign in this narrow layer. The roughly .00202 change in the critical slope reflects the remaining numerical EOS changes jointly; this test does not assign it to one term.

All entries are local rates on the original held trajectory, not finite-window amplification, a re-evolved EOS history, or a physical thermal stability result. No rate is assigned at the PCHIP knot itself.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_cubic_derivative_cross
python -W error fd_portal_cubic_derivative_cross.py "$EOS" results |
  cmp - results/fd-portal-cubic-derivative-cross-2026-11-28.json
```
