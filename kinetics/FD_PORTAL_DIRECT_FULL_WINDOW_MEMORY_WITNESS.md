# Direct nonlinear check of the expanding FD memory from x=.55 to 1.20

The earlier linear result composes two variational integrations at their exact shared x=.96 state. This check instead starts at x=.55 and actually evolves independent +ε and −ε perturbations all the way to 1.20, without reinitializing them at .96. The two initial directions are the orthonormal input directions of largest and smallest ordered number/energy gain, with initial total energy fixed. The unperturbed control reproduces the archived background on the two segments to 4.41e-9 in start-normalized comoving coordinates.

| Initial direction | Ordered dark fractional output norm | Direct dark fractional output norm (ε=2e-4) | Maximum absolute column difference vs ordered | Maximum change from ε=4e-4 to half |
|:---|---:|---:|---:|---:|
| Strong | 1.23877833 | 1.23877845 | 1.16e-6 | 1.20e-6 |
| Weak | 8.49468e-6 | 8.12977e-6 | 5.45e-7 | 6.00e-7 |

The strong/weak contrast is reproduced by the independent trajectories; the weak value is resolved much more coarsely in relative terms. The weak dark-column norm is not a bound on the surviving perturbation in total energy.

The entropy audit expresses the ordered map and the two sampled finite secants in the same positive comoving entropy-curvature metrics. The largest gain is 0.80932701 for the ordered result and 0.80932762 for the sampled secant map, with entropy-metric norm 3.89e-6 for their difference. In contrast, the largest sampled gain in the fractional dark number/energy norm is 1.23877845, because the norms measure different quantities. This is a restricted expanding FD numerical history; it is not a global nonlinear contraction theorem, a freeze-out calculation, or a physical continuum result.

## Reproduce in the kinetics directory

Use the same normalized central EOS file as in the existing direct expansion and entropy-metric receipts (SHA-256 906d2b43bfd925c65aea3ff25ce2840cf5da7a49fcf31d8190199ae1004168ad).

```bash
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_composed_memory_direct_witness.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-comoving-response-055-096-$d.json" "$p-direct-comoving-variational-propagator-$d.json" \
  "$p-comoving-matched-response-055-120-$d.json"
python -W error fd_portal_full_window_witness_entropy_audit.py results
```
