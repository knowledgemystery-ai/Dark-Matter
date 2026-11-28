# The same direct memory in an entropy-curvature norm

Use the previously archived exact combined FD path and its two independent fixed-initial-total-energy response columns. The same input to the same output can have different endpoint gains if the norm changes. Here the metric is the local Hessian of minus the total comoving entropy. It is evaluated at both exact endpoints, not fitted to change the sign of an outcome.

Take positive physical variables (n, rho_dark, rho_total). The dark thermodynamics supply the gradient (eta, −1/T_dark) of minus the dark entropy. The visible component depends on rho_vis=rho_total−rho_dark, introducing the positive rank-one term proportional to (0,−1,1) outer itself, divided by T_vis² times visible heat capacity. Transform to fractional densities by multiplication by the diagonal densities on both sides, and include the same comoving factor a³. Four-point finite differences are obtained independently by inverting the physical bath densities and differentiating this negative-entropy gradient (relative agreement below 3.9e-11). All three metrics are positive definite.

At the input, take the two-dimensional fixed-total-energy subspace. At the output, keep all three variables, including the total energy response from work. The singular values of the response weighted by the metric square roots are:

| Endpoint | Largest in fractional dark coordinates (prior) | Largest in full entropy norm | Smaller in full entropy norm |
|---:|---:|---:|---:|
| .96 | 1.2976283 | 0.8361023 | 1.28846e-4 |
| 1.20 | 1.2387783 | 0.8093270 | 1.26718e-4 |

The separately perturbed interval maps give largest gains 0.8361022 and 0.8093268 in the same metric. If total energy output is dropped instead, the entropy norm is slightly smaller (0.8355811 and 0.8085027). In particular the earlier numerical amplification in fractional coordinates does **not** imply amplification in this entropy curvature metric. This is endpoint contraction at two sampled intervals of one restricted model, not a monotonicity theorem, Lyapunov criterion or freeze-out result.

```bash
cd kinetics
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28; p=results/fd-portal; e=eos-range.dat
python -W error -m unittest -v test_fd_portal_composed_memory_entropy_metric
python -W error fd_portal_composed_memory_entropy_metric.py "$e" \
 "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
 "$p-comoving-response-055-096-$d.json" "$p-direct-comoving-variational-propagator-$d.json" \
 "$p-comoving-matched-response-055-120-$d.json"
```
