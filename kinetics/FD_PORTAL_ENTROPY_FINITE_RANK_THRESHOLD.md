# A fixed-trajectory threshold for all carried finite-window directions

The [integrated rate budget](FD_PORTAL_ENTROPY_FINITE_GAIN_BUDGET.md) evaluated one maximizing direction. Instead integrate the entire two-dimensional Gram rate on the exact same moving carried plane. Begin with the original PCHIP path, not a fresh EOS replacement. Let H be the integrated form without the visible heat-capacity excess above slope 3, and R the integrated excess form. The candidate squared endpoint change is H+αR. Compare H+R independently with the original endpoint form difference.

When H is negative definite and R positive semidefinite, the first nonnegative direction occurs at α=1/λ_max(R,−H). This is a monotone scaling of a **single integrated term on a fixed trajectory**; it does not change the EOS, evolve a different dynamical state, or infer physical stability.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_finite_rank_threshold
python -W error fd_portal_entropy_finite_rank_threshold.py "$EOS" results |
  cmp - results/fd-portal-entropy-finite-rank-threshold-2026-11-28.json
```

## Completed comparison

Both settings give a critical multiplier of **0.95295863** (observed difference 1.12e−9). Removing the term entirely gives a largest counterfactual gain minus one of −4.86169e−5; retaining it gives +2.83050e−6. In this fixed-rate comparison, a reduction of just over 4.704% of the integrated excess eliminates every positive carried direction over the measured window. The endpoint form closure error is below 3.94e−12. This threshold describes scaling the excess term on an unchanged tangent; it is not a tolerance on the physical EOS.
