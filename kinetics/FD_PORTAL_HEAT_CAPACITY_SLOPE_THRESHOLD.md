# Which heat-capacity slopes permit local entropy-norm growth?

The [earlier curvature split](FD_PORTAL_VISIBLE_CURVATURE_RATE_ATTRIBUTION.md) replaced the visible heat-capacity log slope by 3 on one previously selected direction at a time. Here reoptimize the direction over the full perturbation space and each of two physically different planes. Hold the thermal state, all other metric and flow terms, and the [transported initial plane](FD_PORTAL_TRANSPORTED_ENTROPY_PLANE.md) fixed.

At either held state, replacing k=d log(c_v)/d log(T) changes the symmetric rate form by k times a positive rank-one coefficient, up to an additive constant. The reference form S(3) is negative definite in all three spaces. The critical increase above 3 is the inverse largest generalized eigenvalue of that coefficient against −S(3). Independent bisection of the **reoptimized** rate verifies the zero and both sides of the sign change. A separate collision quadrature confirms the thresholds at the two held states.

| Held x | Archived k | Full-space critical k | Critical k on transported initial plane | Largest full / transported rates at k=3 |
|---:|---:|---:|---:|---:|
| .61 | 3.27925 | 3.26495 | 91.80449 | −.22794 / −1.34582 |
| .99 | 3.42516 | 3.37256 | 3.39286 | −.28457 / −.28536 |

The corresponding thresholds on the **instantaneously** fixed-total plane are 94.15475 and 17.87426. The large extrapolated values quantify weak sensitivity of those held directions to this specific rank-one term; they are not plausible thermodynamic slopes or evaluations of other EOS tables. At x=.99 the original k exceeds the critical slopes of both the full and transported spaces, but not the instantaneous plane.

Across the four recorded offsets from the [interpolation knot](FD_PORTAL_ENTROPY_PCHIP_KNOT_CROSSING.md), the reoptimized full-space thresholds occupy a narrow band, 3.37252376 to 3.37253560. The actual hot-side slopes are above every one of those thresholds; the cool-side slopes are below every one. This resolves the earlier sign split as the slope crossing a nearly unchanged numerical threshold. No rate is assigned at the knot itself, and no interval or physical monotonicity is inferred.

This is a derivative-only diagnostic on held numerical states, not an independent EOS evolution, a statement of global physical stability, or a substitute for a finite-time perturbation test.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_heat_capacity_slope_threshold
python -W error fd_portal_heat_capacity_slope_threshold.py "$EOS" results |
  cmp - results/fd-portal-heat-capacity-slope-threshold-2026-11-28.json
```

Use the same exact normalized central EOS table as the earlier derivative and transport receipts.

**Cubic comparison.** The threshold is nearly unchanged across the PCHIP sides, not under an entire EOS change. [Crossing only the cubic derivative](FD_PORTAL_CUBIC_DERIVATIVE_CROSS.md) gives the wrong sign in all four near-knot probes. [The rest of the rate form](FD_PORTAL_CUBIC_REMAINDER_ATTRIBUTION.md) moves the threshold mainly through visible metric motion, partly offset by expansion.
