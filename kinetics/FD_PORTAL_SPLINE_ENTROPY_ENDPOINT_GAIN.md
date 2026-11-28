# The finite endpoint entropy gain is insensitive to the tested interpolation change

The previous alternate cubic EOS reversed two sampled local entropy-rate signs after complete re-evolution. That result concerns an instantaneous derivative. Here independently time-order the alternate restricted FD tangent over the full x=.55 to 1.20 background, then measure endpoint entropy gains for fixed **initial** total energy. Use each history's own positive entropy-curvature metric at both endpoints.

| Final x | Original PCHIP largest gain | Alternate smooth cubic largest gain | Difference |
|---:|---:|---:|---:|
| .96 | .83610232023 | .83610218576 | −1.345e−7 |
| 1.20 | .80932701252 | .80932587405 | −1.138e−6 |

The two alternate ordered Jacobian samplings (65 and 129 points) change the two endpoint gains by at most 8.19e−9 in absolute terms. The alternate unperturbed history is independently replayed within 8.07e−11 in its normalized comoving coordinates. The smaller entropy gain differs by about 4.10e−7 at x=1.20, so its fractional sensitivity is much larger than that of the largest gain. Local derivative signs and finite window gains should not be used interchangeably.

The maximum endpoint gain stays well below 1 for these two numerical backgrounds. This is not a bound on the physical EOS family, a claim of monotone contraction at every time, or a convergence result for cosmological production.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_spline_entropy_endpoint_gain.py "$EOS" \
  "$p-entropy-spline-reintegration-$d.json" "$p-comoving-response-entropy-metric-$d.json"
```

Use the same normalized central EOS and exact independently re-evolved spline receipt as the earlier comparison.
