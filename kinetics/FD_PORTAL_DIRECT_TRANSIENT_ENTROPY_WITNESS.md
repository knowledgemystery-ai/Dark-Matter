# The positive PCHIP local transient also appears in finite signed trajectories

Begin at x=.55 with fixed initial total energy. Select the initial unit fractional direction whose transported plane at x=.99 maximizes the previously published instantaneous PCHIP entropy rate. Evolve actual ±epsilon nonlinear trajectories from .55 to .99. Use their midpoint and centered difference as ODE state variables: this is an exact transformation of two signed trajectories in real arithmetic and avoids subtracting accumulated order-one endpoint integration errors.

| Initial epsilon | Instantaneous squared entropy-norm rate at x=.99 |
|---:|---:|
| 2e−4 | .02768385393 |
| 1e−4 | .02768387363 |
| 5e−5 | .02768397336 |
| Ordered tangent | .02768502789 |

Across the three finite probes the maximum response-column defect against the selected tangent is 1.32e−9 absolute. The rate spread is below 1.20e−7 and the finite rates differ from the ordered result by less than 1.18e−6. The checked quantity is the derivative of the norm of the signed secant, using the original moving entropy metric and nominal background. It does not determine the maximum over all finite input states.

This finite witness supports the original numerical PCHIP result. It does not remove the EOS interpolation dependence established by the smooth-spline reintegration and the nearby one-sided knot test.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -W error -m unittest test_fd_portal_direct_transient_entropy_witness
python -W error fd_portal_direct_transient_entropy_witness.py "$EOS" results
```

Use the same exact normalized central EOS table and transported-plane receipt as the previous study.
