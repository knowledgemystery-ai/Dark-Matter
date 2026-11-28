# Short nonlinear pairs on the separately evolved cubic plane

The [cubic tangent](FD_PORTAL_ENTROPY_CUBIC_SHORT_TANGENT.md) uses its own independently continued background. Carry its archived x=.96 tangent plane to the alternate x=.99 state in two resolutions, check that state against the earlier smooth re-evolution, and QR-rebase there. Starting at x=.99, run signed nonlinear initial-value pairs for both local plane directions, and divide their differences by an independently integrated unperturbed control. This procedure never borrows the original PCHIP plane.

The ratios are estimates from small finite pairs on the same four clock windows, not arbitrary finite-amplitude nonlinear stability bounds or a new evolution starting at .55. The smooth cubic fit is one alternative numerical EOS, not a measurement of the physical heat capacity.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_cubic_short_nonlinear
python -W error fd_portal_entropy_cubic_short_nonlinear.py "$EOS" results |
  cmp - results/fd-portal-entropy-cubic-short-nonlinear-2026-11-28.json
```

## Completed nonlinear check

The finer nonlinear estimate of maximum gain minus one is −2.239443e−6, −4.590945e−7, −5.063732e−7 and −3.204909e−6 across the four preserved windows. Both finite settings find all four signs negative. The finer estimate matches the tangent to within 1.87e−11 across both gains; the two finite settings differ by at most 1.90e−11. This independently confirms the numerical contrast with [the PCHIP nonlinear check](FD_PORTAL_ENTROPY_SHORT_NONLINEAR_WITNESS.md) on the matched clocks.
