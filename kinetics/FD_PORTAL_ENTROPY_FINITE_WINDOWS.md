# Does the local interpolation-dependent entropy growth persist over a finite interval?

The earlier [transported initial plane](FD_PORTAL_TRANSPORTED_ENTROPY_PLANE.md) grew locally at x=.99 under the original numerical interpolation; the [knot audit](FD_PORTAL_ENTROPY_PCHIP_KNOT_CROSSING.md) showed a discontinuous local rate just later. Those results do not determine how much a transported perturbation gains over a finite interval. Here both the original and the [smooth re-evolved background](FD_PORTAL_ENTROPY_SPLINE_REINTEGRATION.md) use common time endpoints set by the **original** knot crossing. The alternate EOS does not use its own crossing to place a different interval.

Only the initial total comoving energy perturbation is zero. Each signed perturbation is evolved from the original starting clock x=.55 in its own numerical background, without resetting or projecting its current total energy at either endpoint. Two independent initial perturbations of fractional comoving dark number and dark energy span the plane. For each EOS, centered differences at two amplitudes estimate the three-by-two endpoint map, then Richardson extrapolation cancels the leading finite-amplitude correction. The maximum gain is the square root of the largest generalized eigenvalue of the two induced entropy curvature forms in the **same initial two coordinates**, not the product of unrelated full-space gains.

A gain above 1 means some initial-plane direction increases its entropy norm over that finite interval; below 1 means all directions on this numerical linearized plane decrease over the interval. Compare departures from 1 directly with the numerical differences measured between ODE settings and amplitudes. Those differences are diagnostics, not rigorous discretization bounds. They do not establish a global or physical stability property and do not replace the separate long-window audit.

```bash
cd kinetics
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
python -W error -m unittest test_fd_portal_entropy_finite_windows
d=2026-11-28
p=results/fd-portal
python -W error fd_portal_entropy_finite_windows.py "$EOS" \
  "$p-direct-expansion-055-096-$d.json" "$p-direct-expansion-096-120-$d.json" \
  "$p-entropy-spline-reintegration-$d.json" "$p-entropy-pchip-knot-crossing-$d.json"
```

## Completed result: the direct long-history differences do not resolve the signs

The [complete finite-window receipt](results/fd-portal-entropy-finite-windows-2026-11-28.json) is retained as a negative numerical control. In the four PCHIP windows the fine extrapolated largest gain differs from one by 2.46e−6 to 1.33e−5, whereas the two ODE settings disagree by 5.54e−5 to 3.05e−4. The two finite amplitudes also differ enough to reverse the estimated sign in every PCHIP window. The original full-history perturbation subtraction is therefore ill-conditioned for resolving these short windows. The smoother cubic has smaller differences, but no short-window sign is claimed from this calculation: its ODE and amplitude changes are comparable to the inferred gain deficits. The saved one-sided **instantaneous** results and longer-interval tangent endpoint controls do not settle these four finite ratios.

**Follow-up:** [Rebased short transported tangents](FD_PORTAL_ENTROPY_SHORT_TANGENT_WINDOWS.md) resolve the four finite windows in both interpolation choices without re-subtracting the full trajectories.
