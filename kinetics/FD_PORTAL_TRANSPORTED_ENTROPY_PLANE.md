# The initially fixed-total-energy perturbation plane can grow locally

On the original PCHIP numerical history, propagate the two independent fixed-total-initial-energy tangent directions continuously from x=.55. Then evaluate the local entropy-metric rate on the resulting plane at the two archived positive full-space points. The current-total-energy-zero restriction is not the same plane: perturbations of comoving total energy subsequently acquire different pressure work.

| x | Full-space max rate | Same-time fixed-total max | Transported initial-plane max |
|---:|---:|---:|---:|
| .61 | +0.01235354 | −1.34177493 | −1.34177301 |
| .99 | +0.04367067 | −0.30642663 | +0.02768503 |

The rates are for the squared entropy norm per unit x. In fixed initial units the tangent uses 65 and 129 sampled chronological Jacobians along a directly reintegrated background. At the two known endpoints, the finer output columns match the independent published ordered composition within 1.66e−8 absolute. At x=.99 the change of the two computed projected rates between meshes is at most 1.97e−5. From the projector angle and full spectral spread we independently bound that difference by 8.51e−4; the conservative bound at x=.61 is 8.91e−3. These bounds compare the two **calculated** planes, not an exact continuum plane.

In the original interpolation there is therefore a numerical counterexample at x=.99 to assigning a monotonically decreasing entropy norm to every perturbation from the fixed-initial-energy plane, despite its net endpoint decay. The separate same-state EOS test has already shown that its full-space local positive sign is interpolation-sensitive. This is not a physical instability result or an unsampled interval claim.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_plane_projection_bound
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_transported_entropy_plane.py "$EOS" results "$p-local-entropy-metric-rates-$d.json"
python -W error fd_portal_transported_entropy_plane_audit.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-transported-entropy-plane-$d.json"
```
