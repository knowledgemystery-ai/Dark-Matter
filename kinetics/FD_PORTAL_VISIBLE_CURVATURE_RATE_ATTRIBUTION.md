# Visible heat capacity contributes to the two local entropy-metric growth signs

On each archived largest-rate direction, separate the moving-metric term into dark, visible heat-capacity, and remaining visible terms, and the flow term into expansion and collision Jacobians in the same fractional basis. This is an algebraic split evaluated at the original FD state and local direction.

| x | Full squared-norm rate | Expansion flow | Collision flow | Visible heat-capacity departure from cubic scaling | Rate after substituting only that log slope |
|---:|---:|---:|---:|---:|---:|
| .61 | +0.01235354 | −0.04547850 | +0.00025934 | +0.24117303 | −0.22881949 |
| .99 | +0.04367067 | −0.07517074 | +0.00029766 | +0.35431658 | −0.31064592 |

The visible term in the entropy Hessian is proportional to 1/(T² c_v). Its dependence on the heat capacity contributes about +2.83 and +2.85 per unit x, largely cancelled by other contributions. On these two directions the parts due to departures of d log(c_v)/d log(T) from 3 are already larger than the small observed positive net rates. Substituting only the log slope by 3 switches each **same direction's** diagnostic derivative to negative; this does not change the underlying state, metric value, rate equations, or direction and is not a calculation with another EOS. It does not establish that the maximum rate after this substitution is negative in every direction. The contribution assigned to the direct collision-flow Jacobian is small on these two specified directions; this is not a claim about all tangent modes or a causal decomposition of the past history.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_visible_curvature_rate_attribution.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-local-entropy-metric-rates-$d.json"
```

Use the same exact normalized central EOS as the local rate receipt.
